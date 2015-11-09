#!/usr/bin/perl
use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use JSON;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org', # alternatively 'useastdb.ensembl.org'
    -user => 'anonymous'
);

my $gene_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Gene' );

my @gene_list = @ARGV;

foreach (@gene_list) {
    my $gene;
    my %gene_obj;
    my $nearest_gene;
    
    $gene = $gene_adaptor->fetch_by_stable_id($_);
    $nearest_gene = $gene->get_nearest_Gene();

    my $file = $gene->stable_id() . ".seq.json";
    
    $gene_obj{'id'} = $gene->stable_id();
    $gene_obj{'nearest_gene'} = $nearest_gene->stable_id();
    $gene_obj{'start'} = $gene->start();
    $gene_obj{'end'} = $gene->end();
    $gene_obj{'sequence'} = $gene->seq();
    
    my $exons = $gene->get_all_Exons();
    my %exon_map;
    while( my $exon = shift @{$exons} ) {
        my %exon_obj;
        $exon_obj{'constitutive'} = $exon->is_constitutive();
        $exon_obj{'id'} = $exon->stable_id();
        $exon_obj{'start'} = $exon->start();
        $exon_obj{'end'} = $exon->end();
        $exon_obj{'sequence'} = $exon->seq->seq;
        
        $exon_map{$exon->stable_id()} = \%exon_obj;
    }
    $gene_obj{'exons'} = \%exon_map;


    my $transcripts = $gene->get_all_Transcripts();
    my %transcript_map;
    while( my $transcript = shift @{$transcripts} ) {
        my %trans_obj;
        my $tmp;
        
        $trans_obj{'id'} = $transcript->stable_id();
        $trans_obj{'5utr'} = "";
        $trans_obj{'3utr'} = "";

        $tmp = $transcript->five_prime_utr();
        if(defined $tmp){
            $trans_obj{'5utr'} = $tmp->seq;
        }

        $tmp = $transcript->three_prime_utr();
        if(defined $tmp){
            $trans_obj{'3utr'} = $tmp->seq;
        }
        
        my @trans_exons = ();
        my $exons = $transcript->get_all_Exons();
        while( my $exon = shift @{$exons} ) {
            push @trans_exons, $exon->stable_id();
        }
        $trans_obj{'exons'} = \@trans_exons;
        $transcript_map{$transcript->stable_id()} = \%trans_obj;
    }
    $gene_obj{'transcripts'} = \%transcript_map;


    my $json = encode_json \%gene_obj;
    
    unless(open FILE, '>'.$file) {
        die "\nUnable to create $file\n";
    }
    print FILE $json;
    close FILE;

}