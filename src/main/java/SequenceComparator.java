import org.apache.lucene.analysis.Analyzer;
import org.apache.lucene.analysis.Tokenizer;
import org.apache.lucene.analysis.ngram.EdgeNGramTokenizer;
import org.apache.lucene.analysis.ngram.NGramTokenizer;
import org.biojava.nbio.core.sequence.DNASequence;
import org.codelibs.minhash.MinHash;

import java.io.Reader;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Created by Mostafa on 11/5/2015.
 */
public class SequenceComparator {

    private DNASequence mFirstSequence;
    private DNASequence mSecondSequence;
    private int mMaxLevels = 0;
    private WordCombinationGenerator mWordGenerator;

    private SequenceComparator(){}
    public SequenceComparator(DNASequence a, DNASequence b, int maxLevels) {
        mFirstSequence = a;
        mSecondSequence = b;
        maxLevels = mMaxLevels;
        //mWordGenerator = new WordCombinationGenerator(maxLevels);
    }


    private HashSet<Integer> convertToHashSet(String seq, int wordLength) {
        HashSet<Integer> h = new HashSet<Integer>();
        for(int x=0;x<seq.length();x+=wordLength){
            String s = seq.substring(x, wordLength);
            h.add(mWordGenerator.lookup(s));
        }
        return h;
    }

    public double getSimilarity(){

        int minN = 6;
        int maxN = 12;

        int hashBit = 2;

        // A base seed for hash functions.
        int seed = 0;

        // The number of hash functions.
        int num = 256;

        // Analyzer for 1-bit 128 hash.
        Tokenizer aTokenizer = new NGramTokenizer(new StringReader(mFirstSequence.getSequenceAsString() ), minN, maxN);
        Tokenizer bTokenizer = new NGramTokenizer(new StringReader(mSecondSequence.getSequenceAsString()), minN, maxN);

        Analyzer aAnalyzer = MinHash.createAnalyzer(aTokenizer, hashBit, seed, num);
        Analyzer bAnalyzer = MinHash.createAnalyzer(bTokenizer, hashBit, seed, num);

        // Calculate a minhash value. The size is hashBit*num.
        try {
            byte[] aminhash = MinHash.calculate(aAnalyzer, mFirstSequence.getSequenceAsString());
            byte[] bminhash = MinHash.calculate(bAnalyzer, mSecondSequence.getSequenceAsString());
            return MinHash.compare(aminhash, bminhash);
        } catch(Exception e){
            return 0;
        }
    }

}





