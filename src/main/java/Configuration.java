import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;

import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

/**
 * Created by Mostafa on 12/6/2015.
 */
public class Configuration {

    private String suggestedSuffix;

    public String getDefaultSuffix(){
        return suggestedSuffix;
    }

    public static class DirectoryEntry {
        String classification;
        String path;
    }

    public static class Fileset {
        String filenamePrefix;
        DirectoryEntry [] directories;
    }

    private ArrayList<Fileset> filesets = new ArrayList<Fileset>();
    private HashMap<String, Object> settings = new HashMap<String, Object>();

    public Fileset [] getFileset(){
        Fileset [] ret = new Fileset[filesets.size()];
        filesets.toArray(ret);
        return ret;
    }

    public Object get(String key) {
        if(!settings.containsKey(key)){
            throw new IndexOutOfBoundsException();
        }
        return settings.get(key);
    }

    public boolean getBool(String key) {
        return (Boolean) get(key);
        /*
        String b = getString(key).toLowerCase();
        if(b.equals("true")){
            return true;
        } else {
            return false;
        }
        */
    }

    public String getString(String key) {
        return (String) get(key);
    }


    public int getInt(String key) {
        long l = (Long) get(key);
        return (int) l;
    }

    public double getDouble(String key) {
        //return Double.parseDouble(getString(key));
        try {
            return (Double) get(key);
        } catch (ClassCastException e) {
            return (double) ((Long)get(key));
        }
    }

    public static Configuration fromJSONFile(String filename) {

        StringBuilder sb = new StringBuilder();
        Configuration ret = new Configuration();
        try {
            File file = new File(filename);
            JSONParser j = new JSONParser();
            JSONObject o = (JSONObject) j.parse(new FileReader(file));

            JSONObject settings = (JSONObject)o.get("settings");
            Iterator<String> keyset = settings.keySet().iterator();
            while(keyset.hasNext()){
                String k = keyset.next();
                if(!k.equals("id")) {
                    sb.append("_" + settings.get(k));
                }
                ret.settings.put(k, settings.get(k));
            }

            JSONArray files = (JSONArray)o.get("files");
            for(int i=0;i<files.size(); i++) {
                JSONObject obj = (JSONObject) files.get(i);
                Fileset fs = new Fileset();
                fs.filenamePrefix = (String)obj.get("filename_prefix");

                JSONObject dir = (JSONObject) obj.get("directories");
                Iterator<String> ks = dir.keySet().iterator();
                ArrayList<DirectoryEntry> arrDe = new ArrayList<DirectoryEntry>();
                while(ks.hasNext()){
                    String k = ks.next();
                    DirectoryEntry de = new DirectoryEntry();
                    de.path = k;
                    de.classification = (String) dir.get(k);
                    arrDe.add(de);
                }
                fs.directories = new DirectoryEntry[arrDe.size()];
                arrDe.toArray(fs.directories);
                ret.filesets.add(fs);
            }


        } catch(Exception e) {
            return null;
        }
        ret.suggestedSuffix = sb.toString();
        return ret;
    }



}
