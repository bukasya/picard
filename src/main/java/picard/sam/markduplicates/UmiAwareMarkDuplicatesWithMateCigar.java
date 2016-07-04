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
import htsjdk.samtools.SAMRecordDuplicateComparator;
import htsjdk.samtools.util.*;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.programgroups.Alpha;

/**
 * This is a simple tool to mark duplicates using the DuplicateSetIterator, DuplicateSet, and SAMRecordDuplicateComparator.
 *
 * Users should continue to use MarkDuplicates in general.  The main motivation of this tool was the fact that 
 * MarkDuplicates has many, many, many useful test cases, but few unit tests for validating individual duplicate sets. To
 * test the DuplicateSetIterator, DuplicateSet, and SAMRecordDuplicateComparator, the most expedient method was to write
 * this tool and make sure it behaves similarly to MarkDuplicates.  Not the best, I know, but good enough.  NH 06/25/2015.
 *
 *
 * See MarkDuplicates for more details.
 *
 * @author fleharty
 */
@CommandLineProgramProperties(
        usage = "Examines aligned records in the supplied SAM or BAM file to locate duplicate molecules. " +
                "All records are then written to the output file with the duplicate records flagged.",
        usageShort = "Examines aligned records in the supplied SAM or BAM file to locate duplicate molecules.",
        programGroup = Alpha.class
)

public class UmiAwareMarkDuplicatesWithMateCigar extends SimpleMarkDuplicatesWithMateCigar {

    @Option(shortName = "EDIT_DISTANCE_TO_JOIN",
            doc = "This option specifies the edit distance of UMIs to join", optional = true)
    public int EDIT_DISTANCE_TO_JOIN = 1;


    @Option(shortName = "ADD_INFERRED_UMI",
            doc = "This option adds the inferred UMI to the bam in the RI tag", optional = true)
    public boolean ADD_INFERRED = false;


    private final Log log = Log.getInstance(UmiAwareMarkDuplicatesWithMateCigar.class);

    /** Stock main method. */
    public static void main(final String[] args) {
        new UmiAwareMarkDuplicatesWithMateCigar().instanceMainWithExit(args);
    }

    @Override
    protected CloseableIterator<DuplicateSet> getDuplicateSetIterator(SamHeaderAndIterator headerAndIterator, SAMRecordDuplicateComparator comparator) {
        return new UmiAwareDuplicateSetIterator(
                    new DuplicateSetIterator(headerAndIterator.iterator,
                    headerAndIterator.header,
                    false,
                    comparator), EDIT_DISTANCE_TO_JOIN, ADD_INFERRED);
    }
}