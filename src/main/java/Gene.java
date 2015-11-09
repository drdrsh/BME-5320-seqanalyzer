import org.biojava.nbio.core.sequence.DNASequence;

import java.util.ArrayList;
import java.util.HashMap;

public class Gene {
    String id;
    long start;
    long end;
    DNASequence sequence;
    HashMap<String, Exon> exons;
    ArrayList<Transcript> transcripts;
}
