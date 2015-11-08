import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;

import java.io.File;
import java.util.*;

/**
 * Created by Mostafa on 11/5/2015.
 */
public class WordCombinationGenerator {

    private HashMap<Integer, ArrayList<String>> mWordList = new HashMap<Integer, ArrayList<String>>();
    private HashMap<String, Integer> mLookupTables = new HashMap<String, Integer>();

    private static final char [] mLetters = new char[] {'A', 'G', 'C', 'T'};

    private int mCurrentLevel = 0;
    private int mRequiredLevels = 0;

    private WordCombinationGenerator(){}

    public WordCombinationGenerator(int k) {
        mRequiredLevels = k;
        recurse();
        buildLookupTables();
    }

    private void buildLookupTables() {
        for(int i=0;i<mRequiredLevels;i++) {
            ArrayList<String> x = mWordList.get(i);
            for(int j=0;j<x.size();j++){
                mLookupTables.put(x.get(j), j);
            }
        }
    }

    public int lookup(String s){
        return mLookupTables.get(s);
    }

    public String [] getCombinations(int level) {
        ArrayList<String> e = mWordList.get(level);
        String[] strArr = new String[e.size()];
        strArr = e.toArray(strArr);
        return strArr;
    }

    void recurse() {
        int n = mLetters.length;

        if(mCurrentLevel > mRequiredLevels){
            return;
        }

        ArrayList<String> c = new ArrayList<String>();

        if (mCurrentLevel == 0) {
            mWordList.put(0, c);
            mCurrentLevel++;
            recurse();
        } else if(mCurrentLevel == 1) {
            mWordList.put(1, c);
            for (int i = 0; i < n; ++i) {
                c.add("" + mLetters[i]);
            }
            mCurrentLevel++;
            recurse();
        } else {
            int previousLevel = mCurrentLevel - 1;
            ArrayList<String> previousLevelList = mWordList.get(previousLevel);
            mWordList.put(mCurrentLevel, c);
            for (int i = 0; i < n; ++i) {
                for(int x=0;x<previousLevelList.size();x++) {
                    c.add( ("" + mLetters[i]) + previousLevelList.get(x));
                }
            }
            mCurrentLevel++;
            recurse();
        }
    }


    // The main recursive method to print all possible strings of length k
     void printAllKLengthRec(int k) {

    }

}





