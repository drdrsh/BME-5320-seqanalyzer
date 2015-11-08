import org.apache.commons.lang3.StringUtils;
import org.apache.lucene.util.ArrayUtil;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.*;

/**
 * Created by Mostafa on 11/5/2015.
 */
public class Dataset {

    private ArrayList<HashMap<String, String>> mData;

    private int mCurrentRow = -1;
    private String [] mFields;
    private Dataset(){}

    public Dataset(String [] fieldIds){
        mFields = new String[fieldIds.length];
        System.arraycopy(fieldIds, 0, mFields, 0, fieldIds.length);
        mData = new ArrayList<HashMap<String, String>>();
    }

    private HashMap<String, String> createEmptyRow() {
        HashMap<String, String> row = new HashMap<String, String>();
        for(int i=0;i<mFields.length;i++) {
            row.put(mFields[i], "");
        }
        return row;
    }

    public HashMap<String, String> getHashMap(int id) {

        HashMap<String, String> map1 = mData.get(id);
        HashMap<String, String> map2 = new HashMap<String, String>();

        Set<Map.Entry<String, String>> set1 = map2.entrySet();
        for (Map.Entry<String, String> e : set1){
            map2.put(e.getKey(), e.getValue());
        }
        return map2;

    }

    public String [] getArray(int id){
        HashMap<String, String> row = mData.get(id);
        String [] arr = new String[mFields.length];
        for(int i=0;i<mFields.length;i++) {
            arr[i] = row.get(mFields[i]);
        }
        return arr;
    }

    public void newRow(){
        mCurrentRow++;
        mData.add(createEmptyRow());
    }

    public void newRow(Object [] arr){
        newRow();
        for(int i=0;i<mFields.length;i++) {
            setValue(mFields[i], arr[i].toString());
        }
    }

    public void setValue(String key, double value){
        setValue(key, Double.toString(value));
    }

    public void setValue(String key, float value){
        setValue(key, Float.toString(value));
    }

    public void setValue(String key, long value){
        setValue(key, Long.toString(value));
    }

    public void setValue(String key, int value){
       setValue(key, Integer.toString(value));
    }

    public void setValue(String key, String value){
        mData.get(mCurrentRow).put(key, value);
    }

    public String getString(String key) {
        return mData.get(mCurrentRow).get(key);
    }

    public double getDouble(String key) {
        return Double.parseDouble(mData.get(mCurrentRow).get(key));
    }

    public float getFloat(String key) {
        return Float.parseFloat(mData.get(mCurrentRow).get(key));
    }

    public int getInt(String key) {
        return Integer.parseInt(mData.get(mCurrentRow).get(key));
    }

    public void exportToCSV(String filename) throws FileNotFoundException, java.io.UnsupportedEncodingException {
        PrintWriter writer = new PrintWriter(filename, "UTF-8");
        writer.println(StringUtils.join(mFields, ","));
        for(int i=0;i<mData.size();i++){
            writer.println(StringUtils.join(getArray(i),","));
        }
        writer.close();
    }

}
