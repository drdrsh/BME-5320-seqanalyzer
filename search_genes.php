<?php
$trial_diffs = array(0.5, 1, 5, 10, 15, 20, 25, 50);
//Diff +/- 10%
$gene_length = array(221162, 63496, 317190, 26452, 25863, 5170, 5715, 9236, 27781, 49322, 36523, 37336, 19861, 27384, 25927, 15063, 39427, 9707, 52484, 31805, 11071, 515541, 11392, 50927, 3429, 23904, 48128, 24169, 12903, 8126, 64902, 8720, 11169, 374224, 4397, 2483, 56312, 24976, 7381, 45988, 11603, 100893, 23914);

$transcript_count = array(15, 25, 30, 27, 9, 11, 9, 11, 7, 6, 9, 11, 5, 19, 19, 21, 19, 12, 6, 15, 7, 15, 8, 16, 20, 23, 19, 4, 29, 6, 6, 4, 3, 7, 6, 15, 16, 9, 10, 15, 8, 5, 15);

$exon_count = array(47, 73, 97, 64, 36, 43, 27, 36, 18, 34, 33, 55, 13, 62, 43, 58, 55, 40, 15, 56, 21, 96, 46, 61, 54, 48, 77, 15, 92, 16, 37, 16, 11, 45, 13, 47, 85, 31, 52, 52, 37, 47, 64);

$id_list = array(  'ENSG00000105971', 'ENSG00000101294', 'ENSG00000149294', 'ENSG00000125462', 'ENSG00000197375', 'ENSG00000159214', 'ENSG00000171863', 'ENSG00000157379', 'ENSG00000196329', 'ENSG00000203965', 'ENSG00000163541', 'ENSG00000168676', 'ENSG00000166295', 'ENSG00000089006', 'ENSG00000159176', 'ENSG00000177082', 'ENSG00000142920', 'ENSG00000197070', 'ENSG00000103187', 'ENSG00000106258', 'ENSG00000134548', 'ENSG00000133812', 'ENSG00000213213', 'ENSG00000109920', 'ENSG00000126088', 'ENSG00000237765', 'ENSG00000114857', 'ENSG00000127074', 'ENSG00000213930', 'ENSG00000111786', 'ENSG00000154438', 'ENSG00000116898', 'ENSG00000198912', 'ENSG00000082068', 'ENSG00000187066', 'ENSG00000204348', 'ENSG00000054690', 'ENSG00000117481', 'ENSG00000185189', 'ENSG00000011485', 'ENSG00000122566', 'ENSG00000187889', 'ENSG00000140553');

$connection = mysqli_connect('useastdb.ensembl.org', 'anonymous', '', 'homo_sapiens_core_82_38', '5306' );

if(!$connection){
    die("Could not connect to database");
}


foreach($id_list as &$id){
    $id = "'$id'";
}
$not_allowed_genes = implode(',', $id_list);


$result = array();

for($i=0;$i<count($gene_length);$i++) {
    foreach($trial_diffs as $trial_diff) {
        $len = $gene_length[$i];
        $ts = $transcript_count[$i];
        $ex = $exon_count[$i];
        $len_margin = intval(ceil($len * ($trial_diff / 100)));
        $ts_margin  = intval(ceil($ts * ($trial_diff / 100)));
        $ex_margin  = intval(ceil($ex * ($trial_diff / 100)));

        $gene_length_range = array($len - $len_margin, $len + $len_margin);
        $transcript_range = array($ts - $ts_margin, $ts + $ts_margin);
        $exon_range = array($ex - $ex_margin, $ex + $ex_margin);

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

        $stable_id = str_replace("'", '', $id_list[$i]);
        $result[$stable_id] = array();
        while ($row = mysqli_fetch_assoc($res)) {
            $result[$stable_id][] = $row['gene_sid'];
        }

        if(count($result[$stable_id]) != 0){
            echo "$stable_id : Found matches at " . $trial_diff . "%\n";
            break;
        }
        echo "$stable_id : Could not find matches at " . $trial_diff . "%\n";
    }
}

mysqli_close($connection);

foreach($result as $key => $value){
    echo $value[array_rand($value)] . " ";
}
//file_put_contents("C:/Users/Mostafa/IdeaProjects/sequence-analyzer/mapping.json", json_encode($result));
