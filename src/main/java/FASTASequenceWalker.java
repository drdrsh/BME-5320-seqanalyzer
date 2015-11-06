import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.io.FastaReader;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.biojava.nbio.core.sequence.io.GenericFastaHeaderParser;
import org.biojava.nbio.core.sequence.io.ProteinSequenceCreator;

import javax.sound.midi.Sequence;
import java.io.File;
import java.io.FileInputStream;
import java.util.*;

/**
 * Created by Mostafa on 11/5/2015.
 */
public class FASTASequenceWalker implements Iterable<FASTASequenceWalker.SequenceEntry> {

    public static class SequenceEntry {
        enum SeqType {TYPE_FULL, TYPE_5UTR, TYPE_3UTR, TYPE_EXON, TYPE_CDS, TYPE_INTRON, TYPE_UNKNOWN};
        String geneId;
        String seqId;
        String seqSubId;
        long length;
        SeqType seqType;
        DNASequence sequence;
        String typeAnnotation;
    }

    private ArrayList<SequenceEntry> mSequenceList;
    HashMap <String, SequenceEntry.SeqType> seqTypeMap;
    private String mFilePath = null;


    private FASTASequenceWalker(){}
    public FASTASequenceWalker(String filepath) {
        mFilePath = filepath;

        mSequenceList = new ArrayList<SequenceEntry>();
        seqTypeMap = new HashMap<String, SequenceEntry.SeqType>();
        seqTypeMap.put("exon", SequenceEntry.SeqType.TYPE_EXON);
        seqTypeMap.put("cds", SequenceEntry.SeqType.TYPE_CDS);
        seqTypeMap.put("utr5", SequenceEntry.SeqType.TYPE_5UTR);
        seqTypeMap.put("utr3", SequenceEntry.SeqType.TYPE_3UTR);

        readFile();
    }

    private SequenceEntry.SeqType translateSeqType(String type){

        if(seqTypeMap.containsKey(type)){
            return seqTypeMap.get(type);
        }

        return SequenceEntry.SeqType.TYPE_UNKNOWN;

    }

    private void readFile(){
        LinkedHashMap<String, DNASequence> a;
        try {
            File file = new File(mFilePath);
            a = FastaReaderHelper.readFastaDNASequence(file);
            for (  Map.Entry<String,DNASequence> entry : a.entrySet() ) {
                SequenceEntry s = createEntryFromHeader(entry.getValue().getOriginalHeader());
                s.sequence = entry.getValue();
                mSequenceList.add(s);
            }
        } catch (Exception e){
            e.printStackTrace();
        }
    }

    private SequenceEntry createEntryFromHeader(String header) {

        SequenceEntry s = new SequenceEntry();

        String [] parts = header.split(" ");
        int typeHeaderIndex = 1;
        if(parts.length == 3){
            typeHeaderIndex = 2;
        }

        String [] firstParts = parts[0].split(":");
        if(firstParts.length == 1) {
            s.seqId = header;
            s.geneId = null;
            s.seqType = SequenceEntry.SeqType.TYPE_FULL;
            return s;
        }

        s.geneId = firstParts[0];
        s.seqId = firstParts[1];
        s.seqSubId = parts[typeHeaderIndex];

        String [] endParts = parts[parts.length-1].split(":");
        s.seqType = translateSeqType(endParts[0]);
        s.typeAnnotation = endParts[1];
        return s;

    }

    public Iterator<SequenceEntry> iterator() {
        final Iterator<SequenceEntry> internalIterator = mSequenceList.iterator();

        return new Iterator<SequenceEntry>() {

            public boolean hasNext() {
                return internalIterator.hasNext();
            }

            public SequenceEntry next() {
                return internalIterator.next();
            }

            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }
}





