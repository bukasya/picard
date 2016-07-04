/*
 * The MIT License
 *
 * Copyright (c) 2016 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package picard.sam.markduplicates;

import htsjdk.samtools.DuplicateSet;
import htsjdk.samtools.DuplicateSetIterator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import picard.PicardException;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

/**
 * Created by fleharty on 5/23/16.
 */
public class UmiAwareDuplicateSetIterator implements CloseableIterator<DuplicateSet> {
    private DuplicateSetIterator wrappedIterator;
    private Iterator<DuplicateSet> nextSetsIterator;
    private List<List<Integer>> adjacencyList;
    private List<Integer> groups;
    int editDistanceToJoin;
    boolean addInferredUmi;

    public UmiAwareDuplicateSetIterator(final DuplicateSetIterator wrappedIterator, final int editDistanceToJoin, boolean addInferredUmi) {
        this.wrappedIterator = wrappedIterator;
        this.editDistanceToJoin = editDistanceToJoin;
        this.addInferredUmi = addInferredUmi;
        nextSetsIterator = Collections.emptyIterator();
    }

    @Override
    public void close() {
        wrappedIterator.close();
    }

    @Override
    public boolean hasNext() {
        return nextSetsIterator.hasNext() || wrappedIterator.hasNext();
    }

    @Override
    public DuplicateSet next() {
        if(!nextSetsIterator.hasNext())
            process(wrappedIterator.next());

        return nextSetsIterator.next();
    }

    // Takes a duplicate set and breaks it up into possible smaller sets according to the UMI,
    // and updates nextSetsIterator to be an iterator on that set of DuplicateSets.
    private void process(final DuplicateSet set) {

        List<SAMRecord> records = set.getRecords();

        groups = new ArrayList<>();
        adjacencyList = new ArrayList<>();

        // If any records don't have the RX attribute... don't do anything special
        for(int i = 0;i < records.size();i++) {
            if(records.get(i).getAttribute("RX") == null) {
                nextSetsIterator = Collections.singleton(set).iterator();
                return;
            } else {
            }
        }

        // Sort records by RX tag
        Collections.sort(records, new Comparator<SAMRecord>() {
            @Override
            public int compare(final SAMRecord lhs, final SAMRecord rhs) {
                if(lhs == null || rhs == null) {
                    return 0;
                } else if(lhs.getAttribute("RX") == null || rhs.getAttribute("RX") == null) {
                    return 0;
                }
                else {
                    return ((String) lhs.getAttribute("RX")).compareTo((String) rhs.getAttribute("RX"));
                }
            }
        });

        int n = records.size();

        // Locate records that have identical UMI sequences
        List<String> uniqueObservedUMIs = new ArrayList<String>();

        uniqueObservedUMIs.add((String) records.get(0).getAttribute("RX"));
        for(int i = 1;i < n;i++) {
            // If the records differ (after sorting) we have a new duplicate set.
            if(!records.get(i).getAttribute("RX").equals(records.get(i-1).getAttribute("RX"))) {
                uniqueObservedUMIs.add((String) records.get(i).getAttribute("RX"));
            }
        }

        // Construct Adjacency List of UMIs that are close
        for(int i = 0;i < uniqueObservedUMIs.size();i++) {
            adjacencyList.add(i, new ArrayList<Integer>());
            groups.add(i, 0);
            for(int j = 0;j < uniqueObservedUMIs.size();j++) {
                if( getEditDistance(uniqueObservedUMIs.get(i), uniqueObservedUMIs.get(j)) <= editDistanceToJoin) {
                    adjacencyList.get(i).add(j);
                }
            }
        }

        // Join Groups
        int nGroups = 0;
        for(int i = 0;i < adjacencyList.size();i++) {
            // Have I visited this yet?
            if(groups.get(i) == 0) {
                // No, I haven't yet seen this
                nGroups++; // We've now seen a new group

                // Depth first search on adjacencyList, setting all the values to a group
                mergeGroups(i, nGroups);
            }
        }

        // Construct DuplicateSetList
        List<DuplicateSet> duplicateSetList= new ArrayList<>();
        for(int i = 0;i < nGroups;i++) {
            DuplicateSet e = new DuplicateSet();
            duplicateSetList.add(e);
        }

        // Assign each record to a duplicate set
        for(int i = 0;i < records.size();i++) {
            String umi = (String) records.get(i).getAttribute("RX");

            // Figure out which group this belongs to
            int recordGroup = groups.get(uniqueObservedUMIs.indexOf((String) umi));
            duplicateSetList.get(recordGroup-1).add(records.get(i));

        }

        // Optionally add the inferred Umi
        if(addInferredUmi) {
            for(int i = 0;i < nGroups;i++) {
                // For each duplicate set identify the most common UMI
                List<SAMRecord> recordsFromDuplicateSet = duplicateSetList.get(i).getRecords();

                // Count the number of times each UMI appears
                Map<String, Long> umiCounts = new HashMap<>();
                for(SAMRecord record : recordsFromDuplicateSet) {
                    String currentUmi = (String) record.getAttribute("RX");
                    Long currentCount = umiCounts.get(currentUmi);
                    if(currentCount == null) {
                        umiCounts.put(currentUmi, new Long(1));
                    } else {
                        umiCounts.put(currentUmi, currentCount + 1);
                    }
                }

                // Find most common UMI
                Map.Entry<String, Long> maxEntry = null;
                for(Map.Entry<String, Long> entry : umiCounts.entrySet()) {
                    if(maxEntry == null || entry.getValue().compareTo(maxEntry.getValue()) > 0) {
                        maxEntry = entry;
                    }
                }

                // Assign the inffered UMI to all reads in the current group
                for(SAMRecord record : recordsFromDuplicateSet) {
                    record.setAttribute("RI", maxEntry.getKey());
                }
            }
        }

        nextSetsIterator = duplicateSetList.iterator();
    }

    private void mergeGroups(final int index, final int group) {
        if(groups.get(index) == 0) {
            groups.set(index, group);
            for (int i = 0; i < adjacencyList.get(index).size(); i++) {
                mergeGroups(adjacencyList.get(index).get(i), group);
            }
        }
    }

    private int getEditDistance(String s1, String s2) {
        if(s1 == null && s2 == null) {
            return 0;
        }
        if(s1.length() != s2.length()) {
            throw new PicardException("Barcode " + s1 + " and " + s2 + " do not have matching lengths.");
        }
        int count = 0;
        for(int i = 0;i < s1.length();i++) {
            if(s1.charAt(i) != s2.charAt(i)) {
                count++;
            }
        }
        return count;
    }
}

