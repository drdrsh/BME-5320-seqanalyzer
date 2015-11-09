import org.biojava.nbio.core.sequence.DNASequence;
import org.json.simple.JSONObject;

public class Exon {
    String id;
    boolean isConstitutive;
    DNASequence sequence;
    long start;
    long end;

    public static Exon fromJSON(JSONObject o) {
        Exon e = new Exon();
        e.id = JSONHelper.getString(o, "id");
        e.start = JSONHelper.getLong(o, "start");
        e.end = JSONHelper.getLong(o, "end");
        try {
            e.isConstitutive = JSONHelper.getBoolean(o, "constitutive");
        } catch (Exception ex) {
            e.isConstitutive = false;
        }
        e.sequence = JSONHelper.getSequence(o, "sequence");
        return e;
    }

}