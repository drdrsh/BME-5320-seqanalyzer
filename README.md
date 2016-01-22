# BME-5320-seqanalyzer

This BME-5320 course project for my group.

This piece of code will build a feature set for a set of genes that produced consistently high errors in the study by 
Fonseca et al. (http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0107026)

Here are some documentation of the tools used in this project

## search_genes.php

This script takes access the Ensembl database directly and returns genes that lie within a given range of values (with arbitrary tolerance) of gene length, transcript count and exon count. It doesn't have a command line interface. For convenience I will provide the SQL query directly should you ever need to look for a gene or a set of genes that satisfy specific criteria.

```sql
SELECT
	g.`gene_id` as gene_id,
	g.`stable_id` as gene_sid
FROM
	`gene` as g
WHERE
	g.`seq_region_end` - g.`seq_region_start` BETWEEN $gene_length_range[0] AND $gene_length_range[1]
AND
	(SELECT COUNT(*) FROM `transcript` as tx WHERE tx.`gene_id` = g.`gene_id`) 
		BETWEEN
	$transcript_range[0] AND $transcript_range[1]
AND
	(SELECT 
		COUNT(DISTINCT et2.exon_id) FROM `transcript` as tx2 
	RIGHT OUTER JOIN 
		`exon_transcript` as et2 ON et2.transcript_id = tx2.transcript_id 
	WHERE 
		tx2.`gene_id` = g.`gene_id`
	) 
	BETWEEN 
		$exon_range[0] and $exon_range[1]
AND
	g.`stable_id` NOT IN ($not_allowed_genes)
```

### Note
This query has been modified to allow searching for 3' and 5' UTR regions with length within a specific range, without this constrain a leaker was introduced into the dataset that led to 3' and 5' UTR length be used by classifiers to predict the state of the gene, controlling for this factor resulted in more realistic results.

## get_data.pl

This script takes Ensembl gene stable IDs as parameters, you can specify more than one ID. 

    perl get_data.pl ENSG00000109920 ENSG00000116898 ENSG00000126088

It will contact the Ensembl APIs, retrieve the selected Human genes and stores them in a json file under the current working directory, the file name will correspond to the gene ID, e.g. ENSG00000109920.seq.json.

To run this script you need to have perl, bioperl and the ensembl APIs installed and configured on your machine.  To do that you can follow the [Ensembl API installation instructions](http://www.ensembl.org/info/docs/api/api_installation.html)

    

    //JSON File Structure
    {
    	"id" : "(string) gene stable id",
    	"start": "(long) gene start offset",
    	"end" : "(long) gene end offset",
    	"nearest_gene": "(string) stable id of the nearest gene",
    	"sequence": "(string) the full sequence of the gene",
    	"exons": "(object) a hashmap of exons indexed by exon stable id",
    	"transcripts": "(object) a hashmap of transcripts indexed by transcript stable id"
    }

    //Exon Object 
	{
    	"id" : "(string) exon stable id",
    	"start" : "(long) exon start offset",
    	"end" : "(long) exon end offset",
    	"constitutive": "(int) whether or not the exon is constitutive (0 = false, 1 = true)",
    	"sequence": "full sequence of the exon"
    }
    
    //Transcript Object
    {
	    "id": "(string) transcript stable id",
	    "3utr": "sequence of the 3' untranslated region (empty string if none)",
	    "5utr": "sequence of the 5' untranslated region (empty string if none)",
	    "exons": "(array) an array of exon stable ids that constitute this transcript"
    }

 
## Data set builder

Dataset builder takes the raw gene information and compiles a series of CSV files that contain a variaty of features that will be fed into Weka for analysis. The program reads a JSON settings file that specifies features to include in the generated dataset, sample configuration files are stored under "config" directory.

### The full set of features are
* gene_length
* gene_complexity
* gene_minhash_score
* gene_***_kmer_count
* exon_count
* non_constitutive_exon_percent
* stat_exon_length
* stat_exon_gccount
* stat_exon_complexity
* stat_exon_minhash_score
* stat_transcript_count
* stat_transcript_length
* stat_transcript_gccount
* stat_exon_per_transcript
* stat_exon_to_exon_distance
* stat_non_constitutive_exon_count
* stat_3utr_length
* stat_3utr_gccount
* stat_5utr_length
* stat_5utr_gccount
* stat_5utr_complexity
* stat_transcript_complexity
* stat_3utr_complexity
* stat_5utr_complexity
* stat_transcript_minhash_score
* stat_3utr_minhash_score
* stat_5utr_minhash_score

### Notes 
===
- features with stat_ prefix will be further broken down into statistical summary parameters like (mean, median and standard deviation)
- Kmer feature is broken down into the counts of all possible k-mer 
  
