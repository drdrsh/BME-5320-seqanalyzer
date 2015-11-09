import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
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

import java.io.*;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Map;
import java.util.Random;
import java.util.zip.GZIPOutputStream;

/**
 * Created by Mostafa on 11/5/2015.
 */
public class EntryPoint {

    public static Dataset data;
    public static final long RANDOM_SEED = 500;

    public static DNASequence ALIGNMENT_REF;
    public static double getCompressionRatio(String str) {
        if (str == null || str.length() == 0) {
            return 0;
        }

        BufferedWriter writer = null;
        ByteArrayOutputStream bos = null;
        try{
            bos = new ByteArrayOutputStream();
            GZIPOutputStream zip = new GZIPOutputStream(bos);
            writer = new BufferedWriter(new OutputStreamWriter(zip, "UTF-8"));
            writer.append(str);
        } catch (Exception e) {
            return -1;
        } finally {
            if(writer != null){
                try {
                    writer.close();
                } catch(Exception e) {
                }
            }
        }
        return ( ((double)str.length() / (double)bos.size()));
    }

    private static double getComplexity(DNASequence d) {
        String maxString = StringUtils.repeat('a', d.getLength());
        double max = getCompressionRatio(maxString);
        double res = getCompressionRatio(d.getSequenceAsString());
        //System.out.println("1.0 - " + res + " - " + max);
        double result = 1 - (res / max);
        if(Double.isInfinite(result)) {
            result = 0;
            System.out.println(d.getSequenceAsString());

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

    public static void main(String[] args) {

        ALIGNMENT_REF = generateRandomSequence(2500);

        data = new Dataset(
                new String[] {
                    "gene_id",
                    "gene_length", "gene_complexity", "gc_count", "gene_score",
                    "exon_count",
                    "mean_exon_length", "median_exon_length", "sd_exon_length",
                    "mean_exon_complexity", "median_exon_complexity", "sd_exon_complexity",
                    "mean_exon_score", "median_exon_score", "sd_exon_score",
                    "number_of_transcripts",
                    "error_level"
                }
        );

        DirectoryReader directoryReader = new DirectoryReader();
        directoryReader.addDirectory("high", "./data/high_error/");
        directoryReader.addDirectory("low",  "./data/low_error/");

        long processStartTime = System.currentTimeMillis();
        JSONGeneReader geneReader = new JSONGeneReader();
        for (DirectoryReader.GeneFileEntry g : directoryReader) {
            long geneStartTime = System.currentTimeMillis();
            processGene(geneReader.readGene(g.seqFilename), g.classId);
            System.out.println("Done gene in " + ((System.currentTimeMillis() - geneStartTime) / 1000) + " seconds");
        }
        System.out.println("Done process in " + ((System.currentTimeMillis() - processStartTime) / 1000) + " seconds");

        try {
            data.exportToCSV("test.csv");
        } catch(Exception e) {
            e.printStackTrace();
        }


    }

    private static void processGene(Gene g, String classId) {

        data.newRow();
        data.setValue("gene_id", g.id);
        data.setValue("error_level", classId);

        int exon_count = 0;
        DescriptiveStatistics exonLength = new DescriptiveStatistics();
        DescriptiveStatistics exonScore = new DescriptiveStatistics();
        DescriptiveStatistics exonComplexity = new DescriptiveStatistics();
        DescriptiveStatistics exonGCCount = new DescriptiveStatistics();


        System.out.println("\r Processing gene " + g.id);
        data.setValue("gene_length", g.sequence.getLength());
        data.setValue("gc_count", g.sequence.getGCCount());
        data.setValue("gene_complexity", getComplexity(g.sequence));
        data.setValue("gene_score", new SequenceComparator(g.sequence, ALIGNMENT_REF, 7).getSimilarity());

        int counter = 0;

        Iterator<Map.Entry<String, Exon>> it = g.exons.entrySet().iterator();
        while(it.hasNext()){
            Map.Entry<String, Exon> e = it.next();
            Exon ex =  e.getValue();
            System.out.println("\r Processing exon " + counter + " (" + ex.id + ")");

            exon_count++;
            exonScore.addValue(new SequenceComparator(ex.sequence, ALIGNMENT_REF, 7).getSimilarity());
            exonComplexity.addValue(getComplexity(ex.sequence));
            exonLength.addValue(ex.sequence.getLength());
            exonGCCount.addValue(ex.sequence.getGCCount());

            counter++;

        }

        data.setValue("exon_count", exon_count);

        data.setValue("mean_exon_length", exonLength.getMean());
        data.setValue("median_exon_length", exonLength.getPercentile(50));
        data.setValue("sd_exon_length", exonLength.getStandardDeviation());

        data.setValue("mean_exon_complexity", exonComplexity.getMean());
        data.setValue("median_exon_complexity", exonComplexity.getPercentile(50));
        data.setValue("sd_exon_complexity", exonComplexity.getStandardDeviation());

        data.setValue("mean_exon_score", exonScore.getMean());
        data.setValue("median_exon_score", exonScore.getPercentile(50));
        data.setValue("sd_exon_score", exonScore.getStandardDeviation());
        data.setValue("number_of_transcripts", g.transcripts.size());
    }

}