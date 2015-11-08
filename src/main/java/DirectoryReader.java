import java.io.File;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 * Created by Mostafa on 11/5/2015.
 */
public class  DirectoryReader implements Iterable<DirectoryReader.GeneFileEntry> {

    private boolean isReady = false;
    public Iterator<GeneFileEntry> iterator() {
        return new Iterator<GeneFileEntry>() {
            private Iterator<Map.Entry<String, GeneFileEntry>> mInternalIt = null;

            public boolean hasNext() {
                if(!isReady){
                    readFiles();
                    mInternalIt = mGeneFileList.entrySet().iterator();
                }
                return mInternalIt.hasNext();
            }

            public GeneFileEntry next() {
                if(!isReady){
                    readFiles();
                    mInternalIt = mGeneFileList.entrySet().iterator();
                }
                Map.Entry<String, GeneFileEntry> e =  mInternalIt.next();
                return e.getValue();
            }

            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }

    public static class GeneFileEntry {
        String classId;
        String geneId;
        String seqFilename;
        String goFilename;
    }

    private static class FilenameParts {
        String geneId;
        String dataType;
    }

    private HashMap<String, GeneFileEntry> mGeneFileList;
    private HashMap<String, String> mDirectoryList;


    private FilenameParts getNameParts(String filename){
        FilenameParts p = new FilenameParts();
        String [] parts = filename.split("\\.");
        p.geneId = parts[0];
        p.dataType = parts[1];
        return p;
    }

    private void readFiles() {

        Iterator<Map.Entry<String, String>> it = mDirectoryList.entrySet().iterator();
        while(it.hasNext()) {
            Map.Entry<String, String> e = it.next();
            String directory = e.getValue();
            String classId = e.getKey();
            File dir = new File(directory);
            File[] directoryListing = dir.listFiles();
            if (directoryListing != null) {
                for (File child : directoryListing) {

                    FilenameParts parts = getNameParts(child.getName());
                    GeneFileEntry current;
                    if (!mGeneFileList.containsKey(parts.geneId)) {
                        current = new GeneFileEntry();
                        current.geneId = parts.geneId;
                        mGeneFileList.put(parts.geneId, current);
                    } else {
                        current = mGeneFileList.get(parts.geneId);
                    }
                    current.classId = classId;

                    if (parts.dataType.equals("seq")) {
                        current.seqFilename = child.getAbsolutePath();
                    }

                    if (parts.dataType.equals("go")) {
                        current.goFilename = child.getAbsolutePath();
                    }

                }
            }
        }
        isReady = true;
    }

    public DirectoryReader() {
        mGeneFileList = new HashMap<String, GeneFileEntry>();
        mDirectoryList = new HashMap<String, String>();
    }

    public void addDirectory(String key, String path){
        mDirectoryList.put(key, path);
    }

}
