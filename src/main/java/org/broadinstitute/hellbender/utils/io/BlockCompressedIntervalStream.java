package org.broadinstitute.hellbender.utils.io;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.seekablestream.SeekableBufferedStream;
import htsjdk.samtools.seekablestream.SeekablePathStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.BlockCompressedStreamConstants;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;

import java.io.*;
import java.nio.file.Path;
import java.util.*;
import java.util.function.BiFunction;
import java.util.function.Function;

public class BlockCompressedIntervalStream {

    // our final empty block adds to the extra information
    // it has a virtual file pointer to the start of the index
    public static final byte[] EMPTY_GZIP_BLOCK_WITH_INDEX_POINTER = {
            BlockCompressedStreamConstants.GZIP_ID1,
            (byte)BlockCompressedStreamConstants.GZIP_ID2,
            BlockCompressedStreamConstants.GZIP_CM_DEFLATE,
            BlockCompressedStreamConstants.GZIP_FLG,
            0, 0, 0, 0, // modification time
            BlockCompressedStreamConstants.GZIP_XFL,
            (byte)BlockCompressedStreamConstants.GZIP_OS_UNKNOWN,
            BlockCompressedStreamConstants.GZIP_XLEN + 12, 0,
            BlockCompressedStreamConstants.BGZF_ID1,
            BlockCompressedStreamConstants.BGZF_ID2,
            BlockCompressedStreamConstants.BGZF_LEN, 0,
            39, 0, // Total block size - 1 as little-endian short
            (byte)'I', (byte)'P', 8, 0, // index pointer extra data

            // 8-byte little-endian long representing file-pointer to beginning of index
            // (this data gets overwritten each time a stream is closed)
            1, 2, 3, 4, 5, 6, 7, 8,

            3, 0, // empty payload
            0, 0, 0, 0, // crc
            0, 0, 0, 0, // uncompressedSize
    };
    public static final int FILE_POINTER_OFFSET = 22;

    // each compressed block of data will have (at least) one of these as a part of the index
    // for each contig that appears in a compressed block the SVInterval tracks the smallest starting
    //   coordinate and largest end coordinate of any object in the block
    // the filePosition member is the virtual file offset of the first object in the block (or, the
    //   first object for a new contig, if there are multiple contigs represented within the block)
    public final static class IndexEntry {
        final SVInterval interval;
        long filePosition;

        public IndexEntry( final SVInterval interval, final long filePosition ) {
            this.interval = interval;
            this.filePosition = filePosition;
        }

        public IndexEntry( final DataInputStream dis ) throws IOException {
            final int contig = dis.readInt();
            final int start = dis.readInt();
            final int end = dis.readInt();
            this.interval = new SVInterval(contig, start, end);
            this.filePosition = dis.readLong();
        }

        public SVInterval getInterval() { return interval; }
        public long getFilePosition() { return filePosition; }

        public void write( final DataOutputStream dos ) throws IOException {
            dos.writeInt(interval.getContig());
            dos.writeInt(interval.getStart());
            dos.writeInt(interval.getEnd());
            dos.writeLong(filePosition);
        }
    }

    // a class for writing arbitrary objects to a block compressed stream with a self-contained index
    // the only restriction is that you must supply a lambda that writes enough state to a DataOutputStream
    //   to allow you to reconstitute the object when you read it back in later AND you have to
    //   return an SVInterval so that we can do indexing.
    public static class Writer <T> {
        final GATKPath path;
        final SAMSequenceDictionary dict;
        final BiFunction<T, DataOutputStream, SVInterval> writeFunc;
        final OutputStream os;
        final BlockCompressedOutputStream bcos;
        final DataOutputStream dos;
        SVInterval lastInterval;
        final List<IndexEntry> indexEntries;
        long blockFilePosition;
        int blockContig;
        int blockStart;
        int blockEnd;
        boolean firstBlockMember;

        public final static int DEFAULT_COMPRESSION_LEVEL = 6;

        public Writer( final GATKPath path,
                       final SAMSequenceDictionary dict,
                       final BiFunction<T, DataOutputStream, SVInterval> writeFunc ) {
            this(path, dict, writeFunc, DEFAULT_COMPRESSION_LEVEL);
        }

        public Writer( final GATKPath path,
                       final SAMSequenceDictionary dict,
                       final BiFunction<T, DataOutputStream, SVInterval> writeFunc,
                       final int compressionLevel ) {
            this.path = path;
            this.dict = dict;
            this.writeFunc = writeFunc;
            this.os = path.getOutputStream();
            this.bcos = new BlockCompressedOutputStream(os, (Path)null, compressionLevel);
            this.dos = new DataOutputStream(bcos);
            this.lastInterval = null;
            this.indexEntries = new ArrayList<>();
            this.firstBlockMember = true;
            writeDictionary();
        }

        private void writeDictionary() {
            try {
                dos.writeInt(dict.size());
                for ( final SAMSequenceRecord rec : dict.getSequences() ) {
                    dos.writeInt(rec.getSequenceLength());
                    dos.writeUTF(rec.getSequenceName());
                }
                dos.flush();
            } catch ( final IOException ioe ) {
                throw new UserException("unable to write dictionary to " + path, ioe);
            }
        }

        public void write( final T tee ) {
            final long prevFilePosition = bcos.getPosition();
            // write the object, and make sure the order is OK (coordinate-sorted intervals)
            final SVInterval interval = writeFunc.apply(tee, dos);
            if ( interval.compareTo(lastInterval) < 0 ) {
                throw new UserException("intervals are not coordinate sorted");
            }

            // if this is the first interval we've seen in a block, just capture the block-start data
            if ( firstBlockMember ) {
                startBlock(prevFilePosition, interval);
                return;
            }

            // if the contig changes emit a new index entry (for the previous contig) and
            //   restart tracking of the block
            if ( interval.getContig() != lastInterval.getContig() ) {
                addIndexEntry();
                startBlock(prevFilePosition, interval);
                return;
            }

            // extend the tracked interval, as necessary
            blockEnd = Math.max(blockEnd, interval.getEnd());
            lastInterval = interval;

            // if writing this element caused a new block to be compressed and added to the file
            if ( isNewBlock(prevFilePosition, bcos.getPosition()) ) {
                addIndexEntry();
                firstBlockMember = true;
            }
        }

        public void close() {
            // take care of any pending index entry, if necessary
            if ( !firstBlockMember ) {
                addIndexEntry();
            }

            try {
                dos.flush(); // complete the data block

                long indexPosition = bcos.getPosition(); // current position is the start of the index

                // write the index entries
                dos.writeInt(indexEntries.size());
                for ( final IndexEntry indexEntry : indexEntries ) {
                    indexEntry.write(dos);
                }
                dos.flush(); // and complete the block

                // write a 0-length terminator block at the end that captures the index position
                final byte[] emptyBlockWithIndexPointer =
                        Arrays.copyOf(EMPTY_GZIP_BLOCK_WITH_INDEX_POINTER,
                                EMPTY_GZIP_BLOCK_WITH_INDEX_POINTER.length);
                for ( int idx = FILE_POINTER_OFFSET; idx != FILE_POINTER_OFFSET + 8; ++idx ) {
                    emptyBlockWithIndexPointer[idx] = (byte)indexPosition;
                    indexPosition >>>= 8;
                }
                os.write(emptyBlockWithIndexPointer);

                bcos.close(false); // we've already handled the terminator block
            } catch ( final IOException ioe ) {
                throw new UserException("unable to add index and close " + path, ioe);
            }
        }

        private void startBlock( final long filePosition, final SVInterval interval ) {
            blockFilePosition = filePosition;
            lastInterval = interval;
            blockContig = interval.getContig();
            blockStart = interval.getStart();
            blockEnd = interval.getEnd();
            firstBlockMember = false;
        }

        private void addIndexEntry() {
            final SVInterval blockInterval = new SVInterval(blockContig, blockStart, blockEnd);
            indexEntries.add(new IndexEntry(blockInterval, blockFilePosition));
        }
    }

    // a class for reading arbitrary objects from a block compressed stream with a self-contained index
    // the only restriction is that you must supply a lambda that reads from a DataInputStream
    //   to reconstitute the object.
    public static final class Reader <T> implements Iterator<T> {
        final GATKPath path;
        final Function<DataInputStream, T> readFunc;
        final long indexFilePointer;
        final SeekableStream ss;
        final BlockCompressedInputStream bcis;
        final DataInputStream dis;
        final SAMSequenceDictionary dict;
        final long dataFilePointer;
        SVIntervalTree<Long> index;

        public Reader( final GATKPath path, final Function<DataInputStream, T> readFunc ) {
            this.path = path;
            this.readFunc = readFunc;
            try {
                this.ss = new SeekablePathStream(path.toPath());
                this.indexFilePointer = findIndexFilePointer(ss);
                this.bcis = new BlockCompressedInputStream(new SeekableBufferedStream(ss));
                this.dis = new DataInputStream(bcis);
            } catch ( final IOException ioe ) {
                throw new UserException("unable to read index location from " + path);
            }

            this.dict = readDictionary();
            this.dataFilePointer = bcis.getPosition(); // having read dictionary, we're pointing at the data
        }

        private SAMSequenceDictionary readDictionary() {
            try {
                final int nRecs = dis.readInt();
                final List<SAMSequenceRecord> seqRecs = new ArrayList<>(nRecs);
                for ( int idx = 0; idx != nRecs; ++idx ) {
                    final int contigSize = dis.readInt();
                    final String contigName = dis.readUTF();
                    seqRecs.add(new SAMSequenceRecord(contigName, contigSize));
                }
                return new SAMSequenceDictionary(seqRecs);
            } catch ( final IOException ioe ) {
                throw new UserException("unable to read dictionary from " + path, ioe);
            }
        }

        public SAMSequenceDictionary getDictionary() { return dict; }

        // restart the iteration from the beginning
        public void reset() {
            try {
                bcis.seek(dataFilePointer);
            } catch ( final IOException ioe ) {
                throw new UserException("failed to seek " + path, ioe);
            }
        }

        // iterate over all objects in the stream
        @Override public boolean hasNext() { return bcis.getPosition() < indexFilePointer; }
        @Override public T next() {
            if ( !hasNext() ) {
                throw new NoSuchElementException();
            }
            return readFunc.apply(dis);
        }

        // find all the objects in the stream inflating just those blocks that might have relevant objects
        // you supply a lambda to filter out objects that don't matter (for example you might reject
        //   objects that don't overlap the interval of interest, or those that aren't completely
        //   contained in the interval, or whatever).
        private class OverlapIterator implements Iterator<T> {
            final SVInterval svInterval;
            final BiFunction<T, SVInterval, Boolean> approveFunc;
            final Iterator<SVIntervalTree.Entry<Long>> indexEntryIterator;
            long blockStartPosition;
            T nextT;

            public OverlapIterator( final SVInterval svInterval,
                                    final BiFunction<T, SVInterval, Boolean> approveFunc ) {
                this.svInterval = svInterval;
                this.approveFunc = approveFunc;
                this.indexEntryIterator = Reader.this.index.overlappers(svInterval);
                advance();
            }

            @Override public boolean hasNext() { return nextT != null; }
            @Override public T next() { final T result = nextT; advance(); return result; }

            private void advance() {
                do {
                    if ( isNewBlock(blockStartPosition, Reader.this.bcis.getPosition()) ) {
                        if ( !indexEntryIterator.hasNext() ) {
                            nextT = null;
                            return;
                        }
                        blockStartPosition = indexEntryIterator.next().getValue();
                        try {
                            Reader.this.bcis.seek(blockStartPosition);
                        } catch ( final IOException ioe ) {
                            throw new UserException("can't find next block for " + path, ioe);
                        }
                    }
                    nextT = Reader.this.next();
                } while ( !approveFunc.apply(nextT, svInterval) );
            }
        }

        public Iterator<T> findOverlappers( final SVInterval svInterval,
                                            final BiFunction<T, SVInterval, Boolean> approveFunc ) {
            if ( index == null ) {
                loadIndex();
            }
            return new OverlapIterator(svInterval, approveFunc);
        }

        private long findIndexFilePointer( final SeekableStream ss ) {
            final int finalBlockLen = EMPTY_GZIP_BLOCK_WITH_INDEX_POINTER.length;
            final byte[] finalBlock = new byte[finalBlockLen];
            try {
                ss.seek(ss.length() - finalBlockLen);
                ss.readFully(finalBlock);
                ss.seek(0);
            } catch ( final IOException ioe ) {
                throw new UserException("unable to read final bgzip block from " + path, ioe);
            }
            for ( int idx = 0; idx != FILE_POINTER_OFFSET; ++idx ) {
                if ( EMPTY_GZIP_BLOCK_WITH_INDEX_POINTER[idx] != finalBlock[idx] ) {
                    throw new UserException(
                            "unable to recover index pointer from final block of " + path);
                }
            }
            for ( int idx = FILE_POINTER_OFFSET + 8; idx != finalBlockLen; ++idx ) {
                if ( EMPTY_GZIP_BLOCK_WITH_INDEX_POINTER[idx] != finalBlock[idx] ) {
                    throw new UserException(
                            "unable to recover index pointer from final block of " + path);
                }
            }
            long indexFilePointer = 0;
            int idx = FILE_POINTER_OFFSET + 8;
            while ( --idx >= FILE_POINTER_OFFSET ) {
                indexFilePointer <<= 8;
                indexFilePointer |= finalBlock[idx];
            }
            return indexFilePointer;
        }

        private void loadIndex() {
            final SVIntervalTree<Long> intervalTree = new SVIntervalTree<>();
            try {
                bcis.seek(indexFilePointer);
                int nEntries = dis.readInt();
                while ( nEntries-- > 0 ) {
                    final IndexEntry entry = new IndexEntry(dis);
                    intervalTree.put(entry.getInterval(), entry.getFilePosition());
                }
                bcis.seek(dataFilePointer);
            } catch ( final IOException ioe ) {
                throw new UserException("unable to read index from " + path, ioe);
            }
            index = intervalTree;
        }
    }

    public static boolean isNewBlock( final long filePosition1, final long filePosition2 ) {
        // upper 48 bits contain the block offset
        // check to see if there are any bit differences in those upper 48 bits
        return ((filePosition1 ^ filePosition2) & ~0xffffL) != 0;
    }
}
