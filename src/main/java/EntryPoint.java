import cc.mallet.fst.CRF;
import cc.mallet.types.Alphabet;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
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

import javax.smartcardio.TerminalFactory;
import java.io.*;
import java.nio.charset.Charset;
import java.text.Format;
import java.util.*;
import java.util.zip.GZIPOutputStream;

/**
 * Created by Mostafa on 11/5/2015.
 */
public class EntryPoint {

    public static Dataset data;
    public static final int RANDOM_SEED = 500;

    public static DNASequence ALIGNMENT_REF;

    public static HashMap<String, Double> mComplexityMap = new HashMap<String, Double>();
    public static HashMap<String, Double> mSimilarityMap = new HashMap<String, Double>();
    public static double mSimilarityMin = Double.MAX_VALUE;
    public static double mSimilarityMax = Double.MIN_VALUE;
    public static double mComplexityMin = Double.MAX_VALUE;
    public static double mComplexityMax = Double.MIN_VALUE;

    public static final double COMPLEXITY_UPPER_LIMIT = 1000;
    public static final double COMPLEXITY_LOWER_LIMIT = 0;
    public static final double COMPLEXITY_CORRECTION  = 0;

    public static final double SIMILARITY_UPPER_LIMIT = 1000;
    public static final double SIMILARITY_LOWER_LIMIT = 0;
    public static final double SIMILARITY_CORRECTION  = 0;


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

    private static double getSimiliarity(DNASequence d) {

        int minN = 6;
        int maxN = 12;

        int hashBit = 2;

        // A base seed for hash functions.
        int seed = EntryPoint.RANDOM_SEED;

        // The number of hash functions.
        int num = 256;

        // Analyzer for 1-bit 128 hash.
        Tokenizer aTokenizer = new NGramTokenizer(new StringReader(d.getSequenceAsString() ), minN, maxN);
        Tokenizer bTokenizer = new NGramTokenizer(new StringReader(ALIGNMENT_REF.getSequenceAsString()), minN, maxN);

        Analyzer aAnalyzer = MinHash.createAnalyzer(aTokenizer, hashBit, seed, num);
        Analyzer bAnalyzer = MinHash.createAnalyzer(bTokenizer, hashBit, seed, num);

        // Calculate a minhash value. The size is hashBit*num.
        try {
            byte[] aminhash = MinHash.calculate(aAnalyzer, d.getSequenceAsString());
            byte[] bminhash = MinHash.calculate(bAnalyzer, ALIGNMENT_REF.getSequenceAsString());
            return MinHash.compare(aminhash, bminhash);
        } catch(Exception e){
            return 0;
        }
    }


    private static double getComplexity(DNASequence d) {
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
        Random rand = new Random(RANDOM_SEED);
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
        return dstMin + destDisplacement;
    }

    public static void main(String[] args) throws InterruptedException {

        System.out.println("Gene Analyzer 0.0.1 by Mostafa Abdelraouf");

        ALIGNMENT_REF = generateRandomSequence(2500);

        data = new Dataset(
                new String[]{
                    "gene_id",
                    "gene_length", "gene_complexity", "gene_gccount", "gene_score",

                    "transcript_count",
                    "exon_count",

                    "mean_exon_length", "median_exon_length", "sd_exon_length",
                    "mean_exon_complexity", "median_exon_complexity", "sd_exon_complexity",
                    "mean_exon_score", "median_exon_score", "sd_exon_score",
                    "mean_exon_gccount", "median_exon_gccount", "sd_exon_gccount",

                    "mean_distance_between_exons", "median_distance_between_exons", "sd_distance_between_exons",
                    "mean_exons_per_trans", "median_exons_per_trans", "sd_exons_per_trans",

                    "mean_trans_length", "median_trans_length", "sd_trans_length",
                    "mean_trans_complexity", "median_trans_complexity", "sd_trans_complexity",
                    "mean_trans_score", "median_trans_score", "sd_trans_score",
                    "mean_trans_gccount", "median_trans_gccount", "sd_trans_gccount",

                    "mean_5utr_length", "median_5utr_length", "sd_5utr_length",
                    "mean_5utr_complexity", "median_5utr_complexity", "sd_5utr_complexity",
                    "mean_5utr_score", "median_5utr_score", "sd_5utr_score",
                    "mean_5utr_gccount", "median_5utr_gccount", "sd_5utr_gccount",

                    "mean_3utr_length", "median_3utr_length", "sd_3utr_length",
                    "mean_3utr_complexity", "median_3utr_complexity", "sd_3utr_complexity",
                    "mean_3utr_score", "median_3utr_score", "sd_3utr_score",
                    "mean_3utr_gccount", "median_3utr_gccount", "sd_3utr_gccount",

                    "error_level"
                }
        );

        HashMap<String, ArrayList<Gene>> geneList = new HashMap<String, ArrayList<Gene>>();

        DirectoryReader directoryReader = new DirectoryReader();

        directoryReader.addDirectory("high", "./data/high_error/");
        directoryReader.addDirectory("low", "./data/low_error/");
//      directoryReader.addDirectory("test", "./data/small_set/");

        long processStartTime = System.currentTimeMillis();
        JSONGeneReader geneReader = new JSONGeneReader();
        TerminalHelper.TerminalSegment geneConsoleSegment = TerminalHelper.addProgressSegment() ;

        int totalGeneCount = 0;
        int counter = 1;
        for (DirectoryReader.GeneFileEntry g : directoryReader) {

            if(!geneList.containsKey(g.classId)){
                geneList.put(g.classId, new ArrayList<Gene>());
            }

            Gene gene = geneReader.readGene(g.seqFilename);
            geneConsoleSegment.update("Preprocessing gene " + g.geneId + "(" + counter  + "/" + directoryReader.getCount() + ")");
            preprocessGene(gene);

            geneList.get(g.classId).add(gene);
            counter++;
            totalGeneCount++;
        }

        //TODO: Normalize
        Iterator<String> featureIds = mSimilarityMap.keySet().iterator();
        while(featureIds.hasNext()) {
            String featureId = featureIds.next();

            double complexity = mComplexityMap.get(featureId);
            System.out.println("Previous Complexity " + complexity);
            complexity -= COMPLEXITY_CORRECTION;
            complexity = transformValue(complexity, mComplexityMin - COMPLEXITY_CORRECTION, mComplexityMax, COMPLEXITY_LOWER_LIMIT, COMPLEXITY_UPPER_LIMIT);
            mComplexityMap.put(featureId, complexity);
            System.out.println("New Complexity " + complexity);

            double similarity = mSimilarityMap.get(featureId);
            System.out.println("Previous Similarity " + similarity);
            similarity -= SIMILARITY_CORRECTION;
            similarity = transformValue(similarity, mSimilarityMin - SIMILARITY_CORRECTION, mSimilarityMax, SIMILARITY_LOWER_LIMIT, SIMILARITY_UPPER_LIMIT);
            System.out.println("New Similarity " + similarity);
            mSimilarityMap.put(featureId, similarity);
        }

        counter = 1;
        Iterator<String> geneClasses = geneList.keySet().iterator();
        while(geneClasses.hasNext()) {
            String classId = geneClasses.next();
            ArrayList<Gene> list = geneList.get(classId);
            for(int i=0;i<list.size();i++) {
                Gene g = list.get(i);
                geneConsoleSegment.update("Processing gene " + g.id + "(" + counter + "/" + totalGeneCount + ")");
                processGene(g, classId);
            }
        }

        System.out.println("\nDone process in " + ((System.currentTimeMillis() - processStartTime) / 1000) + " seconds");

        try {
            data.exportToCSV("test.csv");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private static void addCommonStat(String paramName, DescriptiveStatistics obj)  {
        data.setValue("mean_" + paramName, obj.getMean());
        data.setValue("median_" + paramName, obj.getPercentile(50));
        data.setValue("sd_" + paramName, obj.getStandardDeviation());
    }

    private static void preprocessGene(Gene g) {

        double complexity = 0.0;
        double similarity = 0.0;
        TerminalHelper.TerminalSegment featureConsoleSeg = TerminalHelper.addProgressSegment() ;

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

    private static void processGene(Gene g, String classId) {

        data.newRow();
        data.setValue("gene_id", g.id);
        data.setValue("error_level", classId);

        TerminalHelper.TerminalSegment featureConsoleSeg = TerminalHelper.addProgressSegment() ;

        int exon_count = 0;
        DescriptiveStatistics exonLength = new DescriptiveStatistics();
        DescriptiveStatistics exonScore = new DescriptiveStatistics();
        DescriptiveStatistics exonComplexity = new DescriptiveStatistics();
        DescriptiveStatistics exonGCCount = new DescriptiveStatistics();

        DescriptiveStatistics transLength = new DescriptiveStatistics();
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



        data.setValue("gene_length", g.sequence.getLength());
        data.setValue("gene_gccount", g.sequence.getGCCount());
        data.setValue("gene_complexity", mComplexityMap.get(g.id));
        data.setValue("gene_score", mSimilarityMap.get(g.id));

        int counter = 1;

        Iterator<Map.Entry<String, Exon>> it = g.exons.entrySet().iterator();
        while (it.hasNext()) {
            Map.Entry<String, Exon> e = it.next();
            Exon ex = e.getValue();
            featureConsoleSeg.update(String.format(" - exon %s (%d/%d)", ex.id, counter, g.exons.entrySet().size()));
            exon_count++;
            exonScore.addValue(mSimilarityMap.get(ex.id));
            exonComplexity.addValue(mComplexityMap.get(ex.id));
            exonLength.addValue(ex.sequence.getLength());
            exonGCCount.addValue(ex.sequence.getGCCount());

            counter++;

        }

        data.setValue("exon_count", exon_count);

        addCommonStat("exon_length", exonLength);
        addCommonStat("exon_complexity", exonComplexity);
        addCommonStat("exon_gccount", exonGCCount);
        addCommonStat("exon_score", exonScore);


        counter = 1;
        StringBuilder transSeq;
        Iterator<Transcript> tit = g.transcripts.iterator();
        while (tit.hasNext()) {

            Transcript tr = tit.next();
            transSeq = new StringBuilder();
            featureConsoleSeg.update(String.format(" - transcript %s (%d/%d)", tr.id, counter, g.transcripts.size()));

            for (int i = 0; i < tr.exons.size(); i++) {
                Exon currentExon = tr.exons.get(i);
                Exon previousExon = i == 0 ? null : tr.exons.get(i - 1);

                transExonLength.addValue(currentExon.sequence.getLength());
                if (previousExon != null) {
                    transExonSpacing.addValue(currentExon.start - previousExon.end);
                }
                transSeq.append(currentExon.sequence.getSequenceAsString());
            }
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

        data.setValue("transcript_count", g.transcripts.size());

        addCommonStat("distance_between_exons", transExonSpacing);

        addCommonStat("exons_per_trans", transExonPerTrans);

        addCommonStat("trans_length", transLength);
        addCommonStat("trans_complexity", transComplexity);
        addCommonStat("trans_gccount", transGCCount);
        addCommonStat("trans_score", transScore);

        addCommonStat("5utr_length", trans5UTRLength);
        addCommonStat("5utr_complexity", trans5UTRComplexity);
        addCommonStat("5utr_gccount", trans5UTRGCCount);
        addCommonStat("5utr_score", trans5UTRScore);

        addCommonStat("3utr_length", trans3UTRLength);
        addCommonStat("3utr_complexity", trans3UTRComplexity);
        addCommonStat("3utr_gccount", trans3UTRGCCount);
        addCommonStat("3utr_score", trans3UTRScore);

        featureConsoleSeg.remove();
        featureConsoleSeg = null;
    }

}