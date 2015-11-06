import java.io.File;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 * Created by Mostafa on 11/5/2015.
 */
public class  DirectoryReader implements Iterable<DirectoryReader.GeneFileEntry> {

    public Iterator<GeneFileEntry> iterator() {
        final Iterator<Map.Entry<String, GeneFileEntry>> internalIterator = mGeneFileList.entrySet().iterator();

        return new Iterator<GeneFileEntry>() {

            public boolean hasNext() {
                return internalIterator.hasNext();
            }

            public GeneFileEntry next() {
                Map.Entry<String, GeneFileEntry> e = internalIterator.next();
                return e.getValue();
            }

            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }

    public static class GeneFileEntry {
        String geneId;
        String seqFilename;
        String goFilename;
    }

    private static class FilenameParts {
        String geneId;
        String dataType;
    }

    private HashMap<String, GeneFileEntry> mGeneFileList;

    private String mDirectoryPath = "";

    private DirectoryReader(){}

    private FilenameParts getNameParts(String filename){
        FilenameParts p = new FilenameParts();
        String [] parts = filename.split("\\.");
        p.geneId = parts[0];
        p.dataType = parts[1];
        return p;
    }

    public void open(){
        File dir = new File(mDirectoryPath);
        File[] directoryListing = dir.listFiles();
        if (directoryListing != null) {
            for (File child : directoryListing) {

                FilenameParts parts =  getNameParts(child.getName());
                GeneFileEntry current;

                if(!mGeneFileList.containsKey(parts.geneId)){
                    current = new GeneFileEntry();
                    current.geneId = parts.geneId;
                    mGeneFileList.put(parts.geneId, current);
                } else {
                    current = mGeneFileList.get(parts.geneId);
                }

                if(parts.dataType.equals("seq")){
                    current.seqFilename = child.getAbsolutePath();
                }

                if(parts.dataType.equals("go")){
                    current.goFilename = child.getAbsolutePath();
                }

            }
        }
    }

    public DirectoryReader(String directory) {
        mDirectoryPath = directory;
        mGeneFileList = new HashMap<String, GeneFileEntry>();
    }

}
