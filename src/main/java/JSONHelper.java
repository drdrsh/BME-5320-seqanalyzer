import com.sun.javaws.exceptions.InvalidArgumentException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.json.simple.JSONObject;

import java.io.File;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 * Created by Mostafa on 11/5/2015.
 */
public class JSONHelper {

    public static boolean getBoolean(JSONObject o, String key) throws InvalidArgumentException {
        String res = getString(o, key);

        if(res.equals("true")){
            return true;
        }

        if(res.equals("false")) {
            return false;
        }
        try {
            int x = Integer.parseInt(res);
            if(x == 0){
                return false;
            }
            return true;
        } catch(Exception e) {
            throw new InvalidArgumentException(new String [] {"Invalid format"});
        }
    }

    public static String getString(JSONObject o, String key) {
        return o.get(key).toString();
    }

    public static DNASequence getSequence(JSONObject o, String key) {
        try {
            return new DNASequence(getString(o, key));
        } catch (Exception e) {
            return null;
        }
    }

    public static JSONObject getObject(JSONObject o, String key) {
        return ((JSONObject)o.get(key));
    }

    public static long getLong(JSONObject o, String key) {
        return Long.parseLong(getString(o, key));
    }


}
