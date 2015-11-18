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

        DirectoryReader directoryReader = new DirectoryReader();

        directoryReader.addDirectory("high", "./data/high_error/");
        directoryReader.addDirectory("low", "./data/low_error/");

        long processStartTime = System.currentTimeMillis();
        JSONGeneReader geneReader = new JSONGeneReader();
        int counter = 1;
        TerminalHelper.TerminalSegment geneConsoleSegment = TerminalHelper.addProgressSegment() ;
        for (DirectoryReader.GeneFileEntry g : directoryReader) {
            long geneStartTime = System.currentTimeMillis();
            geneConsoleSegment.update("Processing gene " + g.geneId + "(" + counter  + "/" + directoryReader.getCount() + ")");
            processGene(geneReader.readGene(g.seqFilename), g.classId);
            //System.out.println("Done gene in " + ((System.currentTimeMillis() - geneStartTime) / 1000) + " seconds");
            counter++;
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
        featureConsoleSeg.update(" - analyzing complexity");
        data.setValue("gene_complexity", getComplexity(g.sequence));
        featureConsoleSeg.update(" - analyzing gene sequence");
        data.setValue("gene_score", new SequenceComparator(g.sequence, ALIGNMENT_REF, 7).getSimilarity());

        int counter = 1;

        Iterator<Map.Entry<String, Exon>> it = g.exons.entrySet().iterator();
        while (it.hasNext()) {
            Map.Entry<String, Exon> e = it.next();
            Exon ex = e.getValue();
            featureConsoleSeg.update(String.format(" - exon %s (%d/%d)", ex.id, counter, g.exons.entrySet().size()));
            exon_count++;
            exonScore.addValue(new SequenceComparator(ex.sequence, ALIGNMENT_REF, 7).getSimilarity());
            exonComplexity.addValue(getComplexity(ex.sequence));
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
                transComplexity.addValue(getComplexity(tr.sequence));
                transScore.addValue(new SequenceComparator(tr.sequence, ALIGNMENT_REF, 7).getSimilarity());
                transLength.addValue(tr.sequence.getLength());

                trans5UTRGCCount.addValue(tr.fiveUTR.getGCCount());
                trans5UTRComplexity.addValue(getComplexity(tr.fiveUTR));
                trans5UTRScore.addValue(new SequenceComparator(tr.fiveUTR, ALIGNMENT_REF, 7).getSimilarity());
                trans5UTRLength.addValue(tr.fiveUTR.getLength());

                trans3UTRGCCount.addValue(tr.threeUTR.getGCCount());
                trans3UTRComplexity.addValue(getComplexity(tr.threeUTR));
                trans3UTRScore.addValue(new SequenceComparator(tr.threeUTR, ALIGNMENT_REF, 7).getSimilarity());
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