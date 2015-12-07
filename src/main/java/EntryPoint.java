import be.ac.ulg.montefiore.run.jahmm.*;
import be.ac.ulg.montefiore.run.jahmm.learn.BaumWelchLearner;
import be.ac.ulg.montefiore.run.jahmm.toolbox.MarkovGenerator;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.util.DoubleArray;
import org.apache.lucene.analysis.Analyzer;
import org.apache.lucene.analysis.Tokenizer;
import org.apache.lucene.analysis.ngram.NGramTokenizer;
import org.biojava.nbio.alignment.*;
import org.biojava.nbio.alignment.template.PairwiseSequenceAligner;
import org.biojava.nbio.alignment.template.PairwiseSequenceScorer;
import org.biojava.nbio.alignment.template.SequencePair;
import org.biojava.nbio.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.codelibs.minhash.MinHash;
import org.json.simple.JSONObject;

import javax.smartcardio.TerminalFactory;
import java.io.*;
import java.nio.charset.Charset;
import java.text.Format;
import java.text.NumberFormat;
import java.util.*;
import java.util.zip.GZIPOutputStream;

/**
 * Created by Mostafa on 11/5/2015.
 */
public class EntryPoint {

    public static Configuration activeConfiguration;
    public static Dataset activeDataset;

    public static DNASequence ALIGNMENT_REF;
    public static HashMap<String, Double> mComplexityMap = new HashMap<String, Double>();
    public static HashMap<String, Double> mSimilarityMap = new HashMap<String, Double>();

    public static double mSimilarityMin = Double.MAX_VALUE;
    public static double mSimilarityMax = Double.MIN_VALUE;
    public static double mComplexityMin = Double.MAX_VALUE;
    public static double mComplexityMax = Double.MIN_VALUE;

    /*
    public static final int RANDOM_SEED = 450;
    private static final String OUTPUT_FILE_NAME = RANDOM_SEED + "_-cs+b_k3a.csv";


    public static final boolean CALCULATE_COMPLEXITY = false;
    public static final boolean CALCULATE_SIMILARITY = false;



    public static final double COMPLEXITY_UPPER_LIMIT = 1000;
    public static final double COMPLEXITY_LOWER_LIMIT = 0;
    public static final double COMPLEXITY_CORRECTION  = 0;

    public static final double SIMILARITY_UPPER_LIMIT = 1000;
    public static final double SIMILARITY_LOWER_LIMIT = 0;
    public static final double SIMILARITY_CORRECTION  = 0;

    public static final int KMER_SCALE_FACTOR = 100000;
    public static final KMER_MODE KMODE = KMER_MODE.KMER_ABSOLUTE;

    public enum KMER_MODE {KMER_RELATIVE, KMER_ABSOLUTE};
    */
    /*
    private static class Obs extends ObservationDiscrete<Base> {

        private char c;
        private static Base charToBase(char c){
            c = Character.toLowerCase(c);
            switch(c){
                case 'a': return Base.NT_A;
                case 'g': return Base.NT_G;
                case 't': return Base.NT_T;
                case 'c': return Base.NT_C;
                default : return Base.NT_A;
            }
        }
        public Obs(char c){
            super(charToBase(c));
            this.c = c;
        }

        @Override
        public String toString(NumberFormat numberFormat) {
            return "" + c;
        }
    }
    public enum Base {
        NT_A, NT_G, NT_C, NT_T;

        public ObservationDiscrete<Base> observation() {
            return new ObservationDiscrete<Base>(this);
        }
    };
    static List<List<Obs>> generateSequences(DNASequence seq)
    {
        List<Obs> observations = new ArrayList<Obs>();
        String strSeq = seq.getSequenceAsString();
        for(int i = 0; i<strSeq.length(); i++) {
            observations.add(new Obs(strSeq.charAt(i)));
        }
        List<List<Obs>> sequences = new ArrayList<List<Obs>>();
        sequences.add(observations);
        return sequences;
    }


    public static double getHMMProbability(DNASequence seq) {


        Hmm<ObservationDiscrete<Base>> hmm =
                new Hmm<ObservationDiscrete<Base>>(4,
                        new OpdfDiscreteFactory<Base>(Base.class));

        hmm.setOpdf(0, new OpdfDiscrete(Base.class, new double[] {0.25,0.25,0.25,0.25}));
        hmm.setOpdf(1, new OpdfDiscrete(Base.class, new double[] {0.25,0.25,0.25,0.25}));
        hmm.setOpdf(2, new OpdfDiscrete(Base.class, new double[] {0.25,0.25,0.25,0.25}));
        hmm.setOpdf(3, new OpdfDiscrete(Base.class, new double[] {0.25,0.25,0.25,0.25}));

        hmm.setAij(0, 0, 0.25);
        hmm.setAij(0, 1, 0.25);
        hmm.setAij(0, 2, 0.25);
        hmm.setAij(0, 3, 0.25);

        hmm.setAij(1, 0, 0.25);
        hmm.setAij(1, 1, 0.25);
        hmm.setAij(1, 2, 0.25);
        hmm.setAij(1, 3, 0.25);

        hmm.setAij(2, 0, 0.25);
        hmm.setAij(2, 1, 0.25);
        hmm.setAij(2, 2, 0.25);
        hmm.setAij(2, 3, 0.25);

        hmm.setAij(3, 0, 0.25);
        hmm.setAij(3, 1, 0.25);
        hmm.setAij(3, 2, 0.25);
        hmm.setAij(3, 3, 0.25);


        BaumWelchLearner bwl = new BaumWelchLearner();
        List<List<Obs>> x = generateSequences(ALIGNMENT_REF);
        Hmm<ObservationDiscrete<Base>> learntHmm = bwl.learn(hmm, x);

        List<Obs> testSequence = new ArrayList<Obs>();
        double a = learntHmm.probability(generateSequences(seq).get(0));
        double b = learntHmm.probability(x.get(0));

        System.out.println("");
        System.out.println("prob " + a + "|| prob " + b);
        return 0.0;

    }
    */

    public static LinkedHashMap<String, Double> getKMers(DNASequence dnaSeq, int k) {

        boolean kmodeAbsolute = false;
        if(activeConfiguration.getString("kmers_mode").equals("absolute")){
            kmodeAbsolute = true;
        }


        WordCombinationGenerator m = new WordCombinationGenerator(k);
        String [] possibleKmers = m.getCombinations(k);

        LinkedHashMap<String, Double> map = new LinkedHashMap<String, Double>();
        for(int i=0;i<possibleKmers.length;i++) {
            map.put(possibleKmers[i], 0.0);
        }

        String seq = dnaSeq.getSequenceAsString();
        long seqLength = seq.length();

        for(int i=0;i<seq.length()-(k-1);i++){
            String kmer = seq.substring(i, i+k);
            if(map.containsKey(kmer)) {
                double previousValue = map.get(kmer);
                double addition = 1;
                if(!kmodeAbsolute) {
                    addition = ((double) activeConfiguration.getDouble("kmers_scale") / (double) seqLength);
                }
                map.put(kmer, previousValue + addition );
            } else {
                //To avoid getting kmers containing "N" characters.
                continue;
                //map.put(kmer, 0.0);
            }
        }
        return map;
    }

    public static double getCompressionRatio(String str) {
        if (str == null || str.length() == 0) {
            return 0;
        }

        BufferedWriter writer = null;
        ByteArrayOutputStream bos = null;
        try {
            bos = new ByteArrayOutputStream();
            GZIPOutputStream zip = new GZIPOutputStream(bos);
            writer = new BufferedWriter(new OutputStreamWriter(zip, "UTF-8"));
            writer.append(str);
        } catch (Exception e) {
            return -1;
        } finally {
            if (writer != null) {
                try {
                    writer.close();
                } catch (Exception e) {
                }
            }
        }
        return (((double) str.length() / (double) bos.size()));
    }

    private static double geScore(DNASequence x){
        return 0.0;
    }

    private static double getSimiliarity(DNASequence d) {

        if(!activeConfiguration.getBool("minhash_enabled")) {
            return 1;
        }

        //DNASequence ref = d.getLength()>ALIGNMENT_REF.getLength()?ALIGNMENT_REF:generateRandomSequence(d.getLength());
        //DNASequence ref = generateRandomSequence(d.getLength());
        //DNASequence ref = generateRandomSequence(d.getLength());
        //String refSeqString = ref.getSequenceAsString();

        if(d.getLength() == 0){
            return 0.0;
        }

        String refSeqString = StringUtils.repeat('A', d.getLength());

        int minN = activeConfiguration.getInt("minhash_min");
        int maxN = activeConfiguration.getInt("minhash_max");

        int hashBit = 2;

        int seed = activeConfiguration.getInt("seed");

        int num = activeConfiguration.getInt("minhash_hashfunctions");

        Tokenizer aTokenizer = new NGramTokenizer(new StringReader(d.getSequenceAsString() ), minN, maxN);
        Tokenizer bTokenizer = new NGramTokenizer(new StringReader(refSeqString), minN, maxN);

        Analyzer aAnalyzer = MinHash.createAnalyzer(aTokenizer, hashBit, seed, num);
        Analyzer bAnalyzer = MinHash.createAnalyzer(bTokenizer, hashBit, seed, num);

        try {
            byte[] aminhash = MinHash.calculate(aAnalyzer, d.getSequenceAsString());
            byte[] bminhash = MinHash.calculate(bAnalyzer, refSeqString);
            double val = MinHash.compare(aminhash, bminhash);
            return val;
        } catch(Exception e){
            return 0;
        }
    }


    private static double getComplexity(DNASequence d) {

        if(!activeConfiguration.getBool("complexity_enabled")) {
            return 1;
        }

        if(d.getLength() == 0){
            return 0;
        }
        String maxString = StringUtils.repeat('a', d.getLength());
        double max = getCompressionRatio(maxString);
        double res = getCompressionRatio(d.getSequenceAsString());
        //System.out.println("1.0 - " + res + " - " + max);
        double result = 1 - (res / max);
        if (Double.isInfinite(result)) {
            result = 0;
            //System.out.println(d.getSequenceAsString());
       }

        return result;

    }

    private static DNASequence generateRandomSequence(long length) {
        String bases = "ATGC";
        Random rand = new Random(activeConfiguration.getInt("seed"));
        StringBuilder res = new StringBuilder();
        for (long i = 0; i < length; i++) {
            int randIndex = rand.nextInt(bases.length());
            res.append(bases.charAt(randIndex));
        }

        try {
            return new DNASequence(res.toString());
        } catch (Exception e) {
            return null;
        }
    }

    public static int nwAlign(DNASequence a, DNASequence b) {
        SubstitutionMatrix<NucleotideCompound> matrix = SubstitutionMatrixHelper.getNuc4_4();
        SimpleGapPenalty gapP = new SimpleGapPenalty();
        gapP.setOpenPenalty((short) 5);
        gapP.setExtensionPenalty((short) 2);


        PairwiseSequenceAligner<DNASequence, NucleotideCompound> pair =
                Alignments.getPairwiseAligner(
                        a,
                        b,
                        Alignments.PairwiseSequenceAlignerType.GLOBAL,
                        gapP,
                        matrix
                );

        System.out.println(pair.getScore());
        return (int) Math.round(pair.getScore());
    }

    public static double transformValue(double value, double srcMin, double srcMax, double dstMin, double dstMax){
        double displacement = value - srcMin;
        double displacementPercent = displacement / (srcMax - srcMin);
        double destDisplacement = displacementPercent * (dstMax - dstMin);
        double result = dstMin + destDisplacement;
        return result;
    }

    public static void main(String[] args) throws InterruptedException {

        System.out.println("Gene Analyzer 0.0.1 by Mostafa Abdelraouf");


        if(args.length == 0) {
            System.out.println("No configuration file was specified, exiting");
            System.exit(1);
        }

        ArrayList<String> filesToProcess = new ArrayList<String>();
        filesToProcess.add("./config/simple.json");
        filesToProcess.add("./config/complexity.json");
        filesToProcess.add("./config/minhash.json");
        filesToProcess.add("./config/kmers.json");
        filesToProcess.add("./config/kmers_absolute.json");

        filesToProcess.add("./config/complexity_minhash.json");
        filesToProcess.add("./config/complexity_kmers.json");
        filesToProcess.add("./config/minhash_kmers.json");
        filesToProcess.add("./config/full.json");

        //filesToProcess.add("./config/quick_minhash.json");


        args = filesToProcess.toArray(args);

        long processStartTime = System.currentTimeMillis();
        TerminalHelper.TerminalSegment configConsoleSegment = TerminalHelper.addProgressSegment();
        for(int i =0;i<args.length;i++) {
            configConsoleSegment.update("Processing config file " + args[i]);
            runConfigurationFile(args[i]);
        }
        configConsoleSegment.remove();
        System.out.println("\nDone process in " + ((System.currentTimeMillis() - processStartTime) / 1000) + " seconds");


    }

    private static void runConfigurationFile(String configFile) {

        try {
            activeConfiguration = Configuration.fromJSONFile(configFile);
        } catch (Exception e) {
            System.out.println("Failed to read configuration file, " + e.getMessage());
            System.exit(2);
        }

        ALIGNMENT_REF = generateRandomSequence(2500);

        Configuration.Fileset[] fs = activeConfiguration.getFileset();
        TerminalHelper.TerminalSegment geneConsoleSegment = TerminalHelper.addProgressSegment();
        TerminalHelper.TerminalSegment featureConsoleSegment = TerminalHelper.addProgressSegment();
        for (int i = 0; i < fs.length; i++) {
            activeDataset = new Dataset();

            HashMap<String, ArrayList<Gene>> geneList = new HashMap<String, ArrayList<Gene>>();

            DirectoryReader directoryReader = new DirectoryReader();

            for (int j = 0; j < fs[i].directories.length; j++) {
                directoryReader.addDirectory(
                        fs[i].directories[j].classification,
                        fs[i].directories[j].path
                );
            }
            JSONGeneReader geneReader = new JSONGeneReader();

            int totalGeneCount = 0;
            int counter = 1;
            for (DirectoryReader.GeneFileEntry g : directoryReader) {

                if (!geneList.containsKey(g.classId)) {
                    geneList.put(g.classId, new ArrayList<Gene>());
                }

                Gene gene = geneReader.readGene(g.seqFilename);
                geneConsoleSegment.update(" - Preprocessing gene " + g.geneId + "(" + counter + "/" + directoryReader.getCount() + ")");
                preprocessGene(gene, featureConsoleSegment);

                geneList.get(g.classId).add(gene);
                counter++;
                totalGeneCount++;
            }

            Iterator<String> featureIds = mSimilarityMap.keySet().iterator();
            while (featureIds.hasNext()) {
                String featureId = featureIds.next();

                if (activeConfiguration.getBool("complexity_enabled")) {
                    double complexity = mComplexityMap.get(featureId);
                    complexity -= activeConfiguration.getDouble("complexity_scale_correction");
                    complexity = transformValue(complexity, mComplexityMin - activeConfiguration.getDouble("complexity_scale_correction"), mComplexityMax, activeConfiguration.getDouble("complexity_scale_min"), activeConfiguration.getDouble("complexity_scale_max"));
                    mComplexityMap.put(featureId, complexity);
                }

                if (activeConfiguration.getBool("minhash_enabled")) {
                    double similarity = mSimilarityMap.get(featureId);
                    similarity -= activeConfiguration.getDouble("minhash_scale_correction");
                    similarity = transformValue(similarity, mSimilarityMin - activeConfiguration.getDouble("minhash_scale_correction"), mSimilarityMax, activeConfiguration.getDouble("minhash_scale_min"), activeConfiguration.getDouble("minhash_scale_max"));
                    mSimilarityMap.put(featureId, similarity);
                }

            }

            counter = 1;
            Iterator<String> geneClasses = geneList.keySet().iterator();
            while (geneClasses.hasNext()) {
                String classId = geneClasses.next();
                ArrayList<Gene> list = geneList.get(classId);
                for (int l = 0; l < list.size(); l++) {
                    Gene g = list.get(l);
                    geneConsoleSegment.update(" - Processing gene " + g.id + "(" + counter + "/" + totalGeneCount + ")");
                    processGene(g, classId, featureConsoleSegment);
                }
            }


            try {
                //String filename = activeConfiguration.getString("id") + "_" + fs[i].filenamePrefix + activeConfiguration.getDefaultSuffix() + ".csv";
                String filename = activeConfiguration.getString("id") + "_" + fs[i].filenamePrefix + ".csv";
                activeDataset.exportToCSV("./output/" + filename);
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        geneConsoleSegment.remove();
        featureConsoleSegment.remove();
    }

    private static void addCommonStat(String paramName, DescriptiveStatistics obj)  {
        activeDataset.setValue("mean_" + paramName, obj.getMean());
        activeDataset.setValue("median_" + paramName, obj.getPercentile(50));
        activeDataset.setValue("sd_" + paramName, obj.getStandardDeviation());
    }

    private static void preprocessGene(Gene g, TerminalHelper.TerminalSegment featureConsoleSeg) {

        double complexity = 0.0;
        double similarity = 0.0;

        featureConsoleSeg.update(" - analyzing complexity");
        complexity = getComplexity(g.sequence);

        mComplexityMax = Math.max(mComplexityMax, complexity);
        mComplexityMin = Math.min(mComplexityMin, complexity);
        mComplexityMap.put(g.id, complexity);

        featureConsoleSeg.update(" - analyzing gene sequence");
        similarity = getSimiliarity(g.sequence);
        mSimilarityMax= Math.max(mSimilarityMax, similarity);
        mSimilarityMin = Math.min(mSimilarityMin, similarity);
        mSimilarityMap.put(g.id, similarity);

        int counter = 1;
        Iterator<Map.Entry<String, Exon>> it = g.exons.entrySet().iterator();
        while (it.hasNext()) {

            Map.Entry<String, Exon> e = it.next();
            Exon ex = e.getValue();
            featureConsoleSeg.update(String.format(" - exon %s (%d/%d)", ex.id, counter, g.exons.entrySet().size()));

            complexity = getComplexity(ex.sequence);
            mComplexityMax = Math.max(mComplexityMax, complexity);
            mComplexityMin = Math.min(mComplexityMin, complexity);
            mComplexityMap.put(ex.id, complexity);

            similarity = getSimiliarity(ex.sequence);
            mSimilarityMax= Math.max(mSimilarityMax, similarity);
            mSimilarityMin = Math.min(mSimilarityMin, similarity);
            mSimilarityMap.put(ex.id, similarity);

            counter++;
        }


        counter = 1;
        StringBuilder transSeq;
        Iterator<Transcript> tit = g.transcripts.iterator();
        while (tit.hasNext()) {

            Transcript tr = tit.next();
            transSeq = new StringBuilder();
            featureConsoleSeg.update(String.format(" - transcript %s (%d/%d)", tr.id, counter, g.transcripts.size()));

            for (int i = 0; i < tr.exons.size(); i++) {
                Exon currentExon = tr.exons.get(i);
                transSeq.append(currentExon.sequence.getSequenceAsString());
            }

            try {

                tr.sequence = new DNASequence(transSeq.toString());

                //Transcript Sequence Analysis
                complexity = getComplexity(tr.sequence);
                mComplexityMax = Math.max(mComplexityMax, complexity);
                mComplexityMin = Math.min(mComplexityMin, complexity);
                mComplexityMap.put(tr.id, complexity);

                similarity = getSimiliarity(tr.sequence);
                mSimilarityMax= Math.max(mSimilarityMax, similarity);
                mSimilarityMin = Math.min(mSimilarityMin, similarity);
                mSimilarityMap.put(tr.id, similarity);

                //5UTR sequence Analysis
                complexity = getComplexity(tr.fiveUTR);
                mComplexityMax = Math.max(mComplexityMax, complexity);
                mComplexityMin = Math.min(mComplexityMin, complexity);
                mComplexityMap.put(tr.id + "_fiveUTR", complexity);

                similarity = getSimiliarity(tr.fiveUTR);
                mSimilarityMax= Math.max(mSimilarityMax, similarity);
                mSimilarityMin = Math.min(mSimilarityMin, similarity);
                mSimilarityMap.put(tr.id + "_fiveUTR", similarity);

                //3UTR sequence Analysis
                complexity = getComplexity(tr.threeUTR);
                mComplexityMax = Math.max(mComplexityMax, complexity);
                mComplexityMin = Math.min(mComplexityMin, complexity);
                mComplexityMap.put(tr.id + "_threeUTR", complexity);

                similarity = getSimiliarity(tr.threeUTR);
                mSimilarityMax= Math.max(mSimilarityMax, similarity);
                mSimilarityMin = Math.min(mSimilarityMin, similarity);
                mSimilarityMap.put(tr.id + "_threeUTR", similarity);

            } catch (Exception e) {
                e.printStackTrace();
            }
            counter++;
        }
    }

    private static void processGene(Gene g, String classId, TerminalHelper.TerminalSegment featureConsoleSeg) {

        activeDataset.newRow();
        activeDataset.setValue("gene_id", g.id);

        if(activeConfiguration.getBool("kmers_enabled")) {
            LinkedHashMap<String, Double> kmerScores = getKMers(g.sequence, activeConfiguration.getInt("kmers_length"));
            Iterator<String> keysIterator = kmerScores.keySet().iterator();
            while (keysIterator.hasNext()) {
                String key = keysIterator.next();
                activeDataset.setValue("gene_kmer_" + key, kmerScores.get(key));
            }
        }

        int exon_count = 0;
        int constit_exons = 0;
        DescriptiveStatistics exonLength = new DescriptiveStatistics();
        DescriptiveStatistics exonScore = new DescriptiveStatistics();
        DescriptiveStatistics exonComplexity = new DescriptiveStatistics();
        DescriptiveStatistics exonGCCount = new DescriptiveStatistics();

        DescriptiveStatistics transLength = new DescriptiveStatistics();
        DescriptiveStatistics transConstitutiveFraction = new DescriptiveStatistics();
        DescriptiveStatistics transExonSpacing = new DescriptiveStatistics();
        DescriptiveStatistics transScore = new DescriptiveStatistics();
        DescriptiveStatistics transComplexity = new DescriptiveStatistics();
        DescriptiveStatistics transGCCount = new DescriptiveStatistics();
        DescriptiveStatistics transExonLength = new DescriptiveStatistics();
        DescriptiveStatistics transExonPerTrans = new DescriptiveStatistics();

        DescriptiveStatistics trans5UTRComplexity = new DescriptiveStatistics();
        DescriptiveStatistics trans5UTRScore = new DescriptiveStatistics();
        DescriptiveStatistics trans5UTRLength = new DescriptiveStatistics();
        DescriptiveStatistics trans5UTRGCCount = new DescriptiveStatistics();

        DescriptiveStatistics trans3UTRComplexity = new DescriptiveStatistics();
        DescriptiveStatistics trans3UTRScore = new DescriptiveStatistics();
        DescriptiveStatistics trans3UTRLength = new DescriptiveStatistics();
        DescriptiveStatistics trans3UTRGCCount = new DescriptiveStatistics();

        activeDataset.setValue("gene_length", g.sequence.getLength());
        activeDataset.setValue("gene_gccount", g.sequence.getGCCount());

        if(activeConfiguration.getBool("complexity_enabled")) {
            activeDataset.setValue("gene_complexity", mComplexityMap.get(g.id));
        }

        if(activeConfiguration.getBool("minhash_enabled")) {
            activeDataset.setValue("gene_score", mSimilarityMap.get(g.id));
        }
        int counter = 1;

        Iterator<Map.Entry<String, Exon>> it = g.exons.entrySet().iterator();
        while (it.hasNext()) {
            Map.Entry<String, Exon> e = it.next();
            Exon ex = e.getValue();
            featureConsoleSeg.update(String.format(" - exon %s (%d/%d)", ex.id, counter, g.exons.entrySet().size()));
            exon_count++;
            if(ex.isConstitutive){
                constit_exons++;
            }
            exonScore.addValue(mSimilarityMap.get(ex.id));
            exonComplexity.addValue(mComplexityMap.get(ex.id));
            exonLength.addValue(ex.sequence.getLength());
            exonGCCount.addValue(ex.sequence.getGCCount());

            counter++;

        }

        activeDataset.setValue("gene_exon_nonconstit_percent", ((double)(exon_count-constit_exons)/(double)exon_count) * 100 );
        activeDataset.setValue("exon_count", exon_count);

        addCommonStat("exon_length", exonLength);
        if(activeConfiguration.getBool("complexity_enabled")) {
            addCommonStat("exon_complexity", exonComplexity);
        }
        addCommonStat("exon_gccount", exonGCCount);
        if(activeConfiguration.getBool("minhash_enabled")) {
            addCommonStat("exon_score", exonScore);
        }

        counter = 1;
        StringBuilder transSeq;
        Iterator<Transcript> tit = g.transcripts.iterator();
        while (tit.hasNext()) {

            Transcript tr = tit.next();
            transSeq = new StringBuilder();
            featureConsoleSeg.update(String.format(" - transcript %s (%d/%d)", tr.id, counter, g.transcripts.size()));
            int transConstitExons = 0;

            for (int i = 0; i < tr.exons.size(); i++) {
                Exon currentExon = tr.exons.get(i);
                Exon previousExon = i == 0 ? null : tr.exons.get(i - 1);
                if(currentExon.isConstitutive) {
                    transConstitExons++;
                }

                transExonLength.addValue(currentExon.sequence.getLength());
                if (previousExon != null) {
                    double x = Math.abs(currentExon.start - previousExon.end);
                    transExonSpacing.addValue(x);
                    //transExonSpacing.addValue(currentExon.start - previousExon.end);
                }
                transSeq.append(currentExon.sequence.getSequenceAsString());
            }
            transConstitutiveFraction.addValue( ( (double)(tr.exons.size()-transConstitExons) / (double)tr.exons.size()) * 100 );
            transExonPerTrans.addValue(tr.exons.size());

            try {
                tr.sequence = new DNASequence(transSeq.toString());
                transGCCount.addValue(tr.sequence.getGCCount());
                transScore.addValue(mSimilarityMap.get(tr.id));
                transComplexity.addValue(mComplexityMap.get(tr.id));
                transLength.addValue(tr.sequence.getLength());

                trans5UTRGCCount.addValue(tr.fiveUTR.getGCCount());
                trans5UTRScore.addValue(mSimilarityMap.get(tr.id + "_fiveUTR"));
                trans5UTRComplexity.addValue(mComplexityMap.get(tr.id + "_fiveUTR"));
                trans5UTRLength.addValue(tr.fiveUTR.getLength());

                trans3UTRGCCount.addValue(tr.threeUTR.getGCCount());
                trans3UTRScore.addValue(mSimilarityMap.get(tr.id + "_threeUTR"));
                trans3UTRComplexity.addValue(mComplexityMap.get(tr.id + "_threeUTR"));
                trans3UTRLength.addValue(tr.threeUTR.getLength());

            } catch (Exception e) {
                e.printStackTrace();
            }
            counter++;
        }

        activeDataset.setValue("transcript_count", g.transcripts.size());

        addCommonStat("trans_exon_nonconstit_percent", transConstitutiveFraction);

        addCommonStat("distance_between_exons", transExonSpacing);

        addCommonStat("exons_per_trans", transExonPerTrans);

        addCommonStat("trans_length", transLength);
        addCommonStat("trans_gccount", transGCCount);

        if(activeConfiguration.getBool("complexity_enabled")) {
            addCommonStat("trans_complexity", transComplexity);
            addCommonStat("3utr_complexity", trans3UTRComplexity);
            addCommonStat("5utr_complexity", trans5UTRComplexity);
        }

        if(activeConfiguration.getBool("minhash_enabled")) {
            addCommonStat("trans_score", transScore);
            addCommonStat("5utr_score", trans5UTRScore);
            addCommonStat("3utr_score", trans3UTRScore);
        }

        addCommonStat("5utr_length", trans5UTRLength);
        addCommonStat("5utr_gccount", trans5UTRGCCount);

        addCommonStat("3utr_length", trans3UTRLength);
        addCommonStat("3utr_gccount", trans3UTRGCCount);

        activeDataset.setValue("error_level", classId);

        featureConsoleSeg.remove();
        featureConsoleSeg = null;
    }

}