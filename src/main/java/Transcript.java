import org.biojava.nbio.core.sequence.DNASequence;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;

import java.util.ArrayList;
import java.util.HashMap;

public class Transcript {
    String id;
    ArrayList<Exon> exons;
    DNASequence sequence;
    DNASequence fiveUTR;
    DNASequence threeUTR;
    double threeUTRComplexity;
    double threeUTRsimilarity;
    double fiveUTRcomplexity;
    double fiveUTRsimilarity;
    double overallComplexity;
    double overallSimilarity;




    public static Transcript fromJSON(JSONObject o, HashMap<String, Exon> exonMap) {
        Transcript t = new Transcript();
        t.id = JSONHelper.getString(o, "id");
        t.fiveUTR = JSONHelper.getSequence(o, "5utr");
        t.threeUTR = JSONHelper.getSequence(o, "3utr");
        t.exons = new ArrayList<Exon>();

        JSONArray eArray = (JSONArray)o.get("exons");
        for(int i=0;i<eArray.size();i++) {
            String x = eArray.get(i).toString();
            t.exons.add(exonMap.get(x));
        }

        return t;
    }

}

