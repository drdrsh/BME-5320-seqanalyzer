import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.TranscriptSequence;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import sun.reflect.annotation.ExceptionProxy;

import java.io.File;
import java.io.FileReader;
import java.lang.reflect.Array;
import java.util.*;

/**
 * Created by Mostafa on 11/5/2015.
 */
public class JSONGeneReader {

    public Gene readGene(String filepath){
        Gene g = new Gene();
        g.exons = new HashMap<String, Exon>();
        g.transcripts = new ArrayList<Transcript>();

        LinkedHashMap<String, DNASequence> a;
        try {
            File file = new File(filepath);
            JSONParser j = new JSONParser();
            JSONObject o = (JSONObject) j.parse(new FileReader(file));
            g.start = JSONHelper.getLong(o, "start");
            g.end = JSONHelper.getLong(o, "end");
            g.id = JSONHelper.getString(o, "id");
            g.sequence = JSONHelper.getSequence(o, "sequence");

            JSONObject exons = JSONHelper.getObject(o, "exons");
            Iterator<String> eit = exons.keySet().iterator();
            while(eit.hasNext()){
                String k = eit.next();
                g.exons.put(k, Exon.fromJSON((JSONObject)exons.get(k)));
            }

            JSONObject transcripts = JSONHelper.getObject(o, "transcripts");
            Iterator<String> tit = transcripts.keySet().iterator();
            while(tit.hasNext()){
                String k = tit.next();
                g.transcripts.add(Transcript.fromJSON( (JSONObject)transcripts.get(k), g.exons));
            }

            return g;

        } catch (Exception e) {
            e.printStackTrace();
            return null;
        }
    }


}





