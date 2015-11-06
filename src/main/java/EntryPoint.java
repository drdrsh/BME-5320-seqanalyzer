import org.biojava.nbio.alignment.*;
import org.biojava.nbio.alignment.template.PairwiseSequenceAligner;
import org.biojava.nbio.alignment.template.PairwiseSequenceScorer;
import org.biojava.nbio.alignment.template.SequencePair;
import org.biojava.nbio.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.template.Sequence;

import java.util.Iterator;
import java.util.Random;

/**
 * Created by Mostafa on 11/5/2015.
 */
public class EntryPoint {

    public static final long RANDOM_SEED = 500;

    public static DNASequence ALIGNMENT_REF;

    private static DNASequence generateRandomSequence(long length) {
        String bases = "ATGC";
        Random rand = new Random(RANDOM_SEED);
        StringBuilder res = new StringBuilder();
        for (long i = 0; i < length; i++) {
            int randIndex=rand.nextInt(bases.length());
            res.append(bases.charAt(randIndex));
        }

        try {
            return new DNASequence(res.toString());
        } catch(Exception e){
            return null;
        }
    }

    public static int nwAlign(DNASequence a, DNASequence b){
        SubstitutionMatrix<NucleotideCompound> matrix = SubstitutionMatrixHelper.getNuc4_4();
        SimpleGapPenalty gapP = new SimpleGapPenalty();
        gapP.setOpenPenalty((short)5);
        gapP.setExtensionPenalty((short)2);


        PairwiseSequenceAligner<DNASequence, NucleotideCompound> pair =
                Alignments.getPairwiseAligner(
                        a,
                        b,
                        Alignments.PairwiseSequenceAlignerType.GLOBAL,
                        gapP,
                        matrix
                );

        System.out.println(pair.getScore());
        return (int)Math.round(pair.getScore());
    }
    public static void main(String [] args){
        ALIGNMENT_REF = generateRandomSequence(2500);

        DirectoryReader x = new DirectoryReader("./data/high_error/");
        x.open();
        for(DirectoryReader.GeneFileEntry o : x){

            FASTASequenceWalker fastaReader = new FASTASequenceWalker(o.seqFilename);
            Iterator<FASTASequenceWalker.SequenceEntry> it = fastaReader.iterator();

            while(it.hasNext()) {
                FASTASequenceWalker.SequenceEntry seq = it.next();
                nwAlign(seq.sequence, ALIGNMENT_REF);


            }

        }
    }

}
