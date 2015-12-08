<?php

//define('OUTPUT_MODE', 'all_matches');
define('OUTPUT_MODE', 'one_best_match');


$dataset = array();

$gene_directory = './data/training_high_error/';
$files = scandir($gene_directory);
foreach($files as $f) {
    
    if(strtolower(pathinfo( $f,  PATHINFO_EXTENSION )) != 'json'){
        continue;
    }
    $json = json_decode(file_get_contents($gene_directory . $f), true);
    $dataset[$json['id']] = array();
    $dataset[$json['id']]['gene_length'] = abs(intval($json['end']) - intval($json['start']));
    $dataset[$json['id']]['transcript_count'] = count($json['transcripts']);
    $dataset[$json['id']]['exon_count'] = count($json['exons']);
    $utr5_total = 0;
    $utr3_total = 0;
    foreach($json['transcripts'] as $t){
        $utr5_total += strlen($t['5utr']);
        $utr3_total += strlen($t['3utr']);
    }
    $dataset[$json['id']]['mean_5utr'] =  floatval($utr5_total / $dataset[$json['id']]['transcript_count']);
    $dataset[$json['id']]['mean_3utr'] =  floatval($utr3_total / $dataset[$json['id']]['transcript_count']);
}

$trial_diffs = array(20, 25, 40, 80, 100, 250);
//Diff +/- 10%
/*
$gene_length = array(221162, 63496, 317190, 26452, 25863, 5170, 5715, 9236, 27781, 49322, 36523, 37336, 19861, 27384, 25927, 15063, 39427, 9707, 52484, 31805, 11071, 515541, 11392, 50927, 3429, 23904, 48128, 24169, 12903, 8126, 64902, 8720, 11169, 374224, 4397, 2483, 56312, 24976, 7381, 45988, 11603, 100893, 23914);

$transcript_count = array(15, 25, 30, 27, 9, 11, 9, 11, 7, 6, 9, 11, 5, 19, 19, 21, 19, 12, 6, 15, 7, 15, 8, 16, 20, 23, 19, 4, 29, 6, 6, 4, 3, 7, 6, 15, 16, 9, 10, 15, 8, 5, 15);

$exon_count = array(47, 73, 97, 64, 36, 43, 27, 36, 18, 34, 33, 55, 13, 62, 43, 58, 55, 40, 15, 56, 21, 96, 46, 61, 54, 48, 77, 15, 92, 16, 37, 16, 11, 45, 13, 47, 85, 31, 52, 52, 37, 47, 64);

$id_list = array(  'ENSG00000105971', 'ENSG00000101294', 'ENSG00000149294', 'ENSG00000125462', 'ENSG00000197375', 'ENSG00000159214', 'ENSG00000171863', 'ENSG00000157379', 'ENSG00000196329', 'ENSG00000203965', 'ENSG00000163541', 'ENSG00000168676', 'ENSG00000166295', 'ENSG00000089006', 'ENSG00000159176', 'ENSG00000177082', 'ENSG00000142920', 'ENSG00000197070', 'ENSG00000103187', 'ENSG00000106258', 'ENSG00000134548', 'ENSG00000133812', 'ENSG00000213213', 'ENSG00000109920', 'ENSG00000126088', 'ENSG00000237765', 'ENSG00000114857', 'ENSG00000127074', 'ENSG00000213930', 'ENSG00000111786', 'ENSG00000154438', 'ENSG00000116898', 'ENSG00000198912', 'ENSG00000082068', 'ENSG00000187066', 'ENSG00000204348', 'ENSG00000054690', 'ENSG00000117481', 'ENSG00000185189', 'ENSG00000011485', 'ENSG00000122566', 'ENSG00000187889', 'ENSG00000140553');
*/

$connection = mysqli_connect('ensembldb.ensembl.org', 'anonymous', '', 'homo_sapiens_core_82_38', '5306' );
//$connection = mysqli_connect('useastdb.ensembl.org', 'anonymous', '', 'homo_sapiens_core_82_38', '5306' );

if(!$connection){
    die("Could not connect to database");
}
/*
gene_id -> 10288909
utr5, utr3
219	1906
239	1136

*/

$gene_data = array();
$gene_utr_data = array();
  
function getUTRData($candidate_ids) {

    global $connection, $gene_data;

    for($i=0;$i<count($candidate_ids);$i++) {
        $candidate_ids[$i] = intval($candidate_ids[$i]);
    }
    $loaded_genes = array_keys($gene_data);
    $candidate_ids = array_diff($candidate_ids, $loaded_genes);
    
    if(count($candidate_ids) > 0){
        $transcript_ids = array();
        
        $gene_list = implode(',', $candidate_ids);
        $sql = "SELECT `translation`.*, `transcript`.`gene_id`, `transcript`.`transcript_id` FROM transcript LEFT JOIN translation ON translation.transcript_id = transcript.transcript_id WHERE `gene_id` IN ($gene_list)";
        $res = mysqli_query($connection, $sql);
        while ($row = mysqli_fetch_assoc($res)) {
            $gene_id = intval($row['gene_id']);
            $transcript_id = intval($row['transcript_id']);
            if(!isset($gene_data[$gene_id]['transcripts'])){
                $gene_data[$gene_id]['transcripts'] = array();
            }
            

            $transcript_ids[] = $transcript_id;
            $row['exons'] = array();
            $gene_data[$gene_id]['transcripts'][$transcript_id] = $row;
        }

        $exon_ids = array();
        $transcript_list = implode(',', $transcript_ids);
        $sql = "SELECT transcript.gene_id, exon_transcript.*, (exon.seq_region_end - exon.seq_region_start) as exon_length FROM `exon_transcript`  LEFT JOIN exon on exon.exon_id=exon_transcript.exon_id LEFT JOIN transcript ON transcript.transcript_id=exon_transcript.transcript_id WHERE exon_transcript.transcript_id IN ($transcript_list) ORDER BY exon_transcript.transcript_id, rank ASC";
        $res = mysqli_query($connection, $sql);
        while ($row = mysqli_fetch_assoc($res)) {
            $gene_id = intval($row['gene_id']);
            $transcript_id = intval($row['transcript_id']);
            $row['exon_length'] = intval($row['exon_length']);
            $row['rank'] = intval($row['rank']);
            $gene_data[$gene_id]['transcripts'][$transcript_id]['exons'][] = $row;
        }
        
    }
 
    $resultset = array();
    foreach($candidate_ids as $candidate_id) {
        if(isset($gene_utr_data[$candidate_id])){
            continue;
        }
        $gene = $gene_data[$candidate_id];
        $resultset[$candidate_id] = array(
            'avg_5utr' => 0,
            'avg_3utr' => 0
        );
        $cum_utr5_length = 0;
        $cum_utr3_length = 0;
        foreach($gene['transcripts'] as $trans) {
            $utr5_total_length = 0;
            $utr3_total_length = 0;
            if($trans['translation_id'] === NULL) {
                continue;
            }

            $start_exon_done = false;
            $end_exon_open = false;
            foreach($trans['exons'] as $exon){
                if(!$start_exon_done){
                    if($trans['start_exon_id'] == $exon['exon_id']){
                        $utr5_total_length += intval($trans['seq_start']) + 1;
                        $start_exon_done = true;
                    } else {
                        $utr5_total_length += intval($exon['exon_length']);
                    }
                }
                if($trans['end_exon_id'] == $exon['exon_id']){
                    $utr3_total_length += intval($exon['exon_length']) - intval($trans['seq_end']) + 1;
                    $end_exon_open = true;
                } else {
                    if($end_exon_open){
                        $utr3_total_length += intval($exon['exon_length']);
                    }
                }
            }
            $cum_utr5_length += $utr5_total_length;
            $cum_utr3_length += $utr3_total_length;
        }
        $resultset[$candidate_id]['avg_5utr'] = $cum_utr5_length / count($gene['transcripts']);
        $resultset[$candidate_id]['avg_3utr'] = $cum_utr3_length / count($gene['transcripts']);
    }
    return $resultset;
}


$id_list = array_keys($dataset);
foreach($id_list as &$id){
    $id = "'$id'";
}
$not_allowed_genes = implode(',', $id_list);

$result = array();

foreach($dataset as $gene_stable_id => $gene){
    foreach($trial_diffs as $trial_diff) {
        $len = $gene['gene_length'];
        $ts = $gene['transcript_count'];
        $ex = $gene['exon_count'];
        $utr5 = $gene['mean_5utr'];
        $utr3 = $gene['mean_3utr'];
        
        $len_margin = intval(ceil($len * ($trial_diff / 100)));
        $ts_margin  = intval(ceil($ts * ($trial_diff / 100)));
        $ex_margin  = intval(ceil($ex * ($trial_diff / 100)));
        $utr5_margin  = floatval(($utr5 * ($trial_diff / 100)));
        $utr3_margin  = floatval(($utr3 * ($trial_diff / 100)));

        $gene_length_range = array($len - $len_margin, $len + $len_margin);
        $transcript_range = array($ts - $ts_margin, $ts + $ts_margin);
        $exon_range = array($ex - $ex_margin, $ex + $ex_margin);
        $utr5_range = array($utr5 - $utr5_margin, $utr5 + $utr5_margin);
        $utr3_range = array($utr3 - $utr3_margin, $utr3 + $utr3_margin);
        $sql =
            <<<SQLQuery
            SELECT
        g.`gene_id` as gene_id,
        g.`stable_id` as gene_sid
    FROM
        `gene` as g
    WHERE
        g.`seq_region_end` - g.`seq_region_start` BETWEEN $gene_length_range[0] AND $gene_length_range[1]
            AND
        (SELECT COUNT(*) FROM `transcript` as tx WHERE tx.`gene_id` = g.`gene_id`) BETWEEN $transcript_range[0] AND $transcript_range[1]
            AND
        (SELECT COUNT(DISTINCT et2.exon_id) FROM `transcript` as tx2 RIGHT OUTER JOIN `exon_transcript` as et2 ON et2.transcript_id = tx2.transcript_id WHERE tx2.`gene_id` = g.`gene_id`) BETWEEN $exon_range[0] and $exon_range[1]
            AND
        g.`stable_id` NOT IN ($not_allowed_genes)
SQLQuery;
        $res = mysqli_query($connection, $sql);

        $candidates = array();
        $id_to_stable_map = array();
        
        $acceptedCandidates = array();
        
        while ($row = mysqli_fetch_assoc($res)) {
            $candidates[] = $row['gene_id'];
            $id_to_stable_map[$row['gene_id']] = $row['gene_sid'];
        }
        $utr_data = getUTRData($candidates);
        
        //echo "utr5 range $utr5_range[0] - $utr5_range[1]\n";
        //echo "utr3 range $utr3_range[0] - $utr3_range[1]\n";

        foreach($utr_data as $gene_id => $d) {
            $utr5_in_range = ($d['avg_5utr'] > $utr5_range[0] && $d['avg_5utr'] < $utr5_range[1]);
            $utr3_in_range = ($d['avg_3utr'] > $utr3_range[0] && $d['avg_3utr'] < $utr3_range[1]);
            //echo "Candidate has avg_5utr of $d[avg_5utr] and avg_3utr of $d[avg_3utr]\n";
            if($utr5_in_range && $utr3_in_range) {
                $acceptedCandidates[] = $id_to_stable_map[$gene_id];
            }
        }

        $result[$gene_stable_id] = $acceptedCandidates;
        if(count($result[$gene_stable_id]) != 0){
            echo "$gene_stable_id : Found matches at " . $trial_diff . "%\n";
            break;
        }
        echo "$gene_stable_id : Could not find matches at " . $trial_diff . "%\n";
    }
}

mysqli_close($connection);



foreach($result as $key => $value){
    if(OUTPUT_MODE == 'one_best_match') {
        echo $value[array_rand($value)] . " ";
    } else {
        $counter = 0;
        foreach($value as $v) {
            echo $v. " ";
        }
    }
}
//file_put_contents("C:/Users/Mostafa/IdeaProjects/sequence-analyzer/mapping.json", json_encode($result));
