import org.biojava.nbio.core.sequence.DNASequence;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;

import java.util.*;

public class TerminalHelper {

    public static class TerminalSegment {
        private String currentValue;
        private int id;
        public void update(String v) {
            currentValue = v;
            TerminalHelper.update();
        }
        public void remove() {
            mProgressParts.remove(id);
            TerminalHelper.update();
        }
    }

    private static int mCounter = 0;
    private static TreeMap<Integer, TerminalSegment> mProgressParts = new TreeMap<Integer, TerminalSegment>();
    private static boolean mIsProgressMode = false;
    public static void line(String s) {
        if(!mIsProgressMode) {
            System.out.println(s);
        }
    }

    public static void error(String s) {
        System.out.println(s);
    }


    public static TerminalSegment addProgressSegment(){
        TerminalSegment t = new TerminalSegment();
        t.id = mCounter++;
        mProgressParts.put(t.id, t);
        update();
        return t;
    }

    public static void update(){
        Iterator<Integer> keys = mProgressParts.keySet().iterator();
        StringBuilder s = new StringBuilder("\r");
        while(keys.hasNext()){
            Integer i = keys.next();
            s.append(mProgressParts.get(i).currentValue);
        }
        System.out.print(s);
    }


}

