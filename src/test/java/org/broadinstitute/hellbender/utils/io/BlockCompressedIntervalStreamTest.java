package org.broadinstitute.hellbender.utils.io;

import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Iterator;

import org.broadinstitute.hellbender.utils.io.BlockCompressedIntervalStream.Reader;
import org.broadinstitute.hellbender.tools.walkers.sv.PairedEndAndSplitReadEvidenceCollection.LocusDepth;

public class BlockCompressedIntervalStreamTest extends GATKBaseTest {
    final static String inputFilename = "/home/tsharpe/data/baf/ac.bci";

    @Test
    public void testReader() {
        final Reader<LocusDepth> reader = new Reader<>(new GATKPath(inputFilename), dis -> {
                    try {
                        final int contigId = dis.readInt();
                        final int position = dis.readInt();
                        final int refIdx = dis.read();
                        final int altIdx = dis.read();
                        final int totalDepth = dis.readInt();
                        final int altDepth = dis.readInt();
                        return new LocusDepth(contigId, position, refIdx, altIdx, totalDepth, altDepth);
                    } catch ( final IOException ioe ) {
                        throw new UserException("unable to read from " + inputFilename);
                    }
                });
        final SVInterval svInterval = new SVInterval(0, 21500, 21900);
        final Iterator<LocusDepth> ldItr = reader.findOverlappers(svInterval, (ld, interval) -> {
            final SVInterval ldInterval =
                    new SVInterval(ld.getContigId(), ld.getPosition(), ld.getPosition() + 1);
            return ldInterval.overlaps(interval);
        });
        while ( ldItr.hasNext() ) {
            System.out.println(ldItr.next());
        }
    }
}
