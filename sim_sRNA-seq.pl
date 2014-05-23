#!/usr/bin/perl -w
use strict;
use Getopt::Std;

my $ver_num = 0.1;

my $usage = "sim_sRNA-seq.pl version: $ver_num

Simulate plant small RNA-seq data

Depdency:

samtools

Usage:

sim_sRNA-seq.pl [options] -s soft-masked-genome.fa
  --- OR, LESS PREFERRED  ---
sim_sRNA-seq.pl [options] -d hard-masked-genome.fa -n not-masked-genome.fa

Options:

-s : path to soft-masked genome. Lower-case letters are assumed repeat-masked. 
     If -s specified, neither -d nor -n can be specified.
-d : path to hard-masked genome.
     If -d is specified, -n must also specified, and -s cannot be specified
-n : path to not-masked genome. 
     If -n is specified, -d must also be specified, and -s cannot be specified.
-v : print version number and quit
-h : print help message and quit
-r : desired total number of reads, in millions. Default: 5
-e : per-read probability of a single nt sequencing error (substitution). Default: 1E-4

";

# initialize defaults

our $opt_r = 5;
our $opt_e = 0.0001;
our $opt_h;
our $opt_v;
our $opt_s;
our $opt_d;
our $opt_n;


# get options and validate them
getopt('rehvsdn');

if ($opt_v) {
    print "sim_sRNA-seq.pl version: $ver_num\n";
    exit;
}
if ($opt_h) {
    print "$usage\n";
    exit;
}
unless ($opt_r =~ /^\d+$/) {
    print STDERR "Fatal: Option -r must be an integer.\n$usage\n";
    exit;
}
unless(($opt_e >=0 ) and ($opt_e < 1)) {
    print STDERR "Fatal: Invalid value for option -e.\n$usage\n";
    exit;
}


# convert option r to millions
my $total_reads = 1000000 * $opt_r;

# check for samtools installation
(open(SAMCHECK, "which samtools |")) || die "FATAL: samtools check failed on command \'which samtools\'\n$usage\n";
my $samtools_check = <SAMCHECK>;
close SAMCHECK;
unless($samtools_check) {
    print "FATAL: samtools check failed. samtools must be installed to use this script.\n$usage\n";
}

# verify genome files
if($opt_s) {
    # opt_d and opt_n not allowed
    if(($opt_d) or ($opt_n)) {
	print "FATAL: If option -s is specified, options -d and -n cannot be specified.\n$usage\n";
	exit;
    }
} elsif (($opt_d) and ($opt_n)) {
    # opt_s not allowed
    if($opt_s) {
	print "FATAL: If options -d and/or -n are specified, option -s cannot be specified.\n$usage\n";
	exit;
    }
} else {
    # no valid combination of genome(s) specified
    print "FATAL: Invalid genome options. Either specify -s or both -d and -n.\n$usage\n";
    exit;
}

# validate genomes. file paths are global
my $validated_genome = validate_genome();
unless($validated_genome) {
    die "\n$usage\n";
}

# build a hash of one-based cumulative positions for nts in the genome, based on order in the fai file
my %genome_positions = get_genome_positions();
my $gp_max = get_gp_max(\%genome_positions);

# test
#print "\ngp_max: $gp_max\n";
#exit;

# determine reads allowed for each class

my $mir_reads = int (0.3 * $total_reads);
my $tasi_reads = int (0.05 * $total_reads);
my $het_reads = $total_reads - $mir_reads - $tasi_reads;

# populate arrays with numbers for each locus
my $mir_num = 100;  ## NOMINALLY get 100 MIRNAs .. actual number could differ a little. This is really the number of divisions on my imaginary log10 plot.
my @mir_ns = get_nts($mir_reads,$mir_num);

my $tasi_num = 20;
my @tasi_ns = get_nts($tasi_reads,$tasi_num);

my $het_num = 10000;
my @het_ns = get_nts($het_reads,$het_num);

# Open output streams
my $out_summary = "simulated_sRNA-seq_reads_overview.txt";
my $out_fasta = "simulated_sRNA-seq_reads.fasta";


# no overwrites
if(-e $out_summary) {
    print "FATAL: Output file $out_summary already exists - no overwrites allowed\n";
    exit;
}
if(-e $out_fasta) {
    print "FATAL: Output file $out_fasta already exists - no overwrites allowed\n";
    exit;
}

(open(OUTS, ">$out_summary")) || die "Fatal: file open error\n";
(open(OUTF, ">$out_fasta")) || die "Fatal: file open error\n";
######################

# hash for occupied regions (those already selected)
# structure: keys: chromosome name, value: anonymous array of integers, where each pair is a coordinate pair. 
# this list is NOT sorted in any way
my %occupied = ();
my $simulated_read;
my $simulated_read_location;
my $simulated_read_errors;
my $fasta_header;

## miRNA simulation
my $mir_locus_n = 0;

foreach my $mir_count (@mir_ns) {
    my $mir_locus_size = 125;
    my $mir_genome_type = "masked";
    my ($mir_locus,$mir_locus_for_seq) = get_a_locus(\%occupied,\%genome_positions,\$gp_max, \$mir_locus_size, \$mir_genome_type);
    
    my $mir_strand = pick_a_strand();
    my($mature_mir,$star_mir) = pick_an_arm($mir_locus);
    
    # test
    #print "mir_locus: $mir_locus mir_count: $mir_count mir_strand: $mir_strand mature_mir: $mature_mir star_mir: $star_mir\n";
    
    ++$mir_locus_n;
    print OUTS "MIRNA_$mir_locus_n\t$mir_locus\t$mir_strand\t$mir_count\n";
	
    # simulate expression
    # get all possible perfect sub-sequences and their coordinates
    # keys: "miR0:0" or "star0:0" form. Values: tab-delimited coordinates, sequence
    my %mir_perfects = get_mir_perfects($mir_locus,$mir_locus_for_seq,$mir_strand,$mature_mir,$star_mir);
    my $mir_key;
    
    for(my $i = 1; $i <= $mir_count; ++$i) {
	$mir_key = sample_mir_keys();
	unless(exists($mir_perfects{$mir_key})) {
	    die "FATAL: key $mir_key not found in hash mir_perfects .. blame Axtell!\n";
	}
	my @mir_vals = split("\t", $mir_perfects{$mir_key});
	$simulated_read_location = $mir_vals[0];
	($simulated_read,$simulated_read_errors) = errorify_read($mir_vals[1]);
	
	# test
	#print "\tsimulated_read_location: $simulated_read_location simuated_read: $simulated_read simulated_read_errors: $simulated_read_errors\n";
	$fasta_header = ">MIRNA_$mir_locus_n" . "_" . "$i" . "_" . "$simulated_read_location" . "_" . "$mir_strand" . "_" . "$simulated_read_errors";
	print OUTF "$fasta_header\n$simulated_read\n";
    }
}

## tasiRNA simulation
my $tasi_locus_n = 0;

foreach my $tasi_count (@tasi_ns) {
    my $tasi_locus_size = 140;
    my $tasi_genome_type = "masked";
    my ($tasi_locus,$tasi_locus_for_seq) = get_a_locus(\%occupied,\%genome_positions,\$gp_max, \$tasi_locus_size, \$tasi_genome_type);
    
    ++$tasi_locus_n;
    print OUTS "TAS_$tasi_locus_n\t$tasi_locus\t\.\t$tasi_count\n";
    
    # simulate tasiRNA expression
    # get all possible tasiRNAs to be sampled from .. 
    # 6 phases of 21mers, with 3'end variants of 20mers, 22mers.
    # keys are in the form "top_[phase]_[size]", or "bot_[phase]_[size]", where phase is 0-5, and size is 20, 21, or 22
    # values are tab-delimited read coordinates, and sequence
    my %tasi_perfects = get_tasi_perfects($tasi_locus,$tasi_locus_for_seq);
    my $tasi_key;
    for(my $j = 1; $j <= $tasi_count; ++$j) {
	$tasi_key = sample_tasi_keys();
	unless(exists($tasi_perfects{$tasi_key})) {
	    die "FATAL: key $tasi_key not found in hash tasi_perfects .. blame Axtell!\n";
	}
	my @tasi_vals = split ("\t", $tasi_perfects{$tasi_key});
	$simulated_read_location = $tasi_vals[0];
	($simulated_read,$simulated_read_errors) = errorify_read($tasi_vals[1]);
	my $tasi_read_strand;
	if($tasi_key =~ /^top_/) {
	    $tasi_read_strand = "+";
	} else {
	    $tasi_read_strand = "-";
	}
	$fasta_header = ">TAS_$tasi_locus_n" . "_" . "$j" . "_" . "$simulated_read_location" . "_" . "$tasi_read_strand" . "_" . "$simulated_read_errors";
	print OUTF "$fasta_header\n$simulated_read\n";
    }
}

## heterochromatic siRNA simulation
my $het_locus_n = 0;
foreach my $het_count (@het_ns) {
    my $het_locus_size = 100;
    my $het_genome_type = "notmasked";
    my ($het_locus,$het_locus_for_seq) = get_a_locus(\%occupied,\%genome_positions,\$gp_max, \$het_locus_size, \$het_genome_type);
    
    ++$het_locus_n;
    print OUTS "HET_$het_locus_n\t$het_locus\t\.\t$het_count\n";
    
    # simulate het siRNA expresion
    # get all possible het siRNA to be sampled from .. 
    
    my %het_perfects = get_het_perfects($het_locus,$het_locus_for_seq);
    my $het_key;
    for(my $k = 1; $k <= $het_count; ++$k) {
	$het_key = sample_het_keys($het_locus_for_seq);
	unless(exists($het_perfects{$het_key})) {
	    die "FATAL: failed to find key $het_key in hash het_perfects .. blame Axtell!\n";
	}
	my @het_vals = split ("\t", $het_perfects{$het_key});
	$simulated_read_location = $het_vals[0];
	($simulated_read,$simulated_read_errors) = errorify_read($het_vals[1]);
	my $het_read_strand;
	if($het_key =~ /^top_/) {
	    $het_read_strand = "+";
	} else {
	    $het_read_strand = "-";
	}
	$fasta_header = ">HET_$het_locus_n" . "_" . "$k" . "_" . "$simulated_read_location" . "_" . "$het_read_strand" . "_" . "$simulated_read_errors";
	print OUTF "$fasta_header\n$simulated_read\n";
    }
}
close OUTF;
close OUTS;


########
sub sample_het_keys {
    my($seq) = @_;
    my $strand_pick = rand();
    my $strand;
    if($strand_pick < 0.5) {
	$strand = "top";
    } else {
	$strand = "bot";
    }
    my $pos_limit = (length $seq) - 25;
    my $pos = int(rand($pos_limit));
    my $size_pick = rand();
    my $size;
    if($size_pick < 0.9) {
	$size = 24;
    } elsif (($size_pick >= 0.9) and ($size_pick < 0.95)) {
	$size = 23;
    } elsif (($size_pick >= 0.95) and ($size_pick < 0.98)) {
	$size = 22;
    } else {
	$size = 21;
    }
    my $key = "$strand" . "_" . "$pos" . "_" . "$size";
    return $key;
}

sub get_het_perfects {
    my($locus,$for_seq) = @_;
    # parse locus
    my $chr;
    my $start;
    my $stop;
    if($locus =~ /^(\S+):(\d+)-(\d+)$/) {
	$chr = $1;
	$start = $2;
	$stop = $3;
    } else {
	die "FATAL in sub-routine get_het_perfects .. could not parse locus $locus\n";
    }
    my $read_start;
    my $read_stop;
    my $read_loc;
    my $read_seq;
    my $key;
    my $value;
    my %hash = ();
    
    for(my $i = 0; $i < ((length $for_seq) - 24); ++$i) {
	# 24
	$read_seq = substr($for_seq,$i,24);
	$read_start = $start + $i;
	$read_stop = $start + $i + 24 - 1;
	$read_loc = "$chr" . ":" . "$read_start" . "-" . "$read_stop";
	$key = "top_" . "$i" . "_" . "24";
	$value = "$read_loc\t$read_seq";
	$hash{$key} = $value;
	# now bot
	$read_seq = revcomp(substr($for_seq,$i,24));
	$key = "bot_" . "$i" . "_" . "24";
	$value = "$read_loc\t$read_seq";
	$hash{$key} = $value;
	
	# top strand, 23
	$read_seq = substr($for_seq,$i,23);
	$read_start = $start + $i;
	$read_stop = $start + $i + 23 - 1;
	$read_loc = "$chr" . ":" . "$read_start" . "-" . "$read_stop";
	$key = "top_" . "$i" . "_" . "23";
	$hash{$key} = $value;
	# now bot
	$read_seq = revcomp(substr($for_seq,$i,23));
	$key = "bot_" . "$i" . "_" . "23";
	$value = "$read_loc\t$read_seq";
	$hash{$key} = $value;
	
	# top strand, 22
	$read_seq = substr($for_seq,$i,22);
	$read_start = $start + $i;
	$read_stop = $start + $i + 22 - 1;
	$read_loc = "$chr" . ":" . "$read_start" . "-" . "$read_stop";
	$key = "top_" . "$i" . "_" . "22";
	$hash{$key} = $value;
	# now bot
	$read_seq = revcomp(substr($for_seq,$i,22));
	$key = "bot_" . "$i" . "_" . "22";
	$value = "$read_loc\t$read_seq";
	$hash{$key} = $value;

	# top strand, 21
	$read_seq = substr($for_seq,$i,21);
	$read_start = $start + $i;
	$read_stop = $start + $i + 21 - 1;
	$read_loc = "$chr" . ":" . "$read_start" . "-" . "$read_stop";
	$key = "top_" . "$i" . "_" . "21";
	$hash{$key} = $value;
	# now bot
	$read_seq = revcomp(substr($for_seq,$i,21));
	$key = "bot_" . "$i" . "_" . "21";
	$value = "$read_loc\t$read_seq";
	$hash{$key} = $value;
    }
    return %hash;
}
	

sub sample_tasi_keys {
    my $strand_pick = rand();
    my $strand;
    if($strand_pick < 0.5) {
	$strand = "top";
    } else {
	$strand = "bot";
    }
    my $phase = int(rand(5));
    my $size_pick = rand();
    my $size;
    if($size_pick < 0.8) {
	$size = 21;
    } elsif (($size_pick >= 0.8) and ($size_pick < 0.9)) {
	$size = 20;
    } else {
	$size = 22;
    }
    my $key = "$strand" . "_" . "$phase" . "_" . "$size";
    return $key;
}

sub get_tasi_perfects {
    my($locus,$locus_for_seq) = @_;
    my %hash = ();
    # parse locus info
    my $chr;
    my $start;
    my $stop;
    if($locus =~ /^(\S+):(\d+)-(\d+)$/) {
	$chr = $1;
	$start = $2;
	$stop = $3;
    } else {
	die "FATAL in sub-routine get_tasi_perfects .. failed to parse locus $locus\n";
    }
    my $key;
    my $phase;
    my $offset;
    my $top_seq;
    my $top_start;
    my $top_stop;
    my $top_loc;
    my $value;
    my $bot_seq;
    my $bot_start;
    my $bot_stop;
    my $bot_loc;
    for($phase = 0; $phase <= 5; ++$phase) {
	# top strand first
	$offset = 4 + (21 * $phase);
	# 21 mer
	$top_seq = substr($locus_for_seq,$offset,21);
	$top_start = $start + $offset;
	$top_stop = $start + $offset + 21 - 1;
	$top_loc = "$chr" . ":" . "$top_start" . "-" . "$top_stop";
	$key = "top_$phase" . "_21";
	$value = "$top_loc\t$top_seq";
	$hash{$key} = $value;
	# 20 mer
	$top_seq = substr($locus_for_seq,$offset,20);
	$top_start = $start + $offset;
	$top_stop = $start + $offset + 20 - 1;
	$top_loc = "$chr" . ":" . "$top_start" . "-" . "$top_stop";
	$key = "top_$phase" . "_20";
	$value = "$top_loc\t$top_seq";
	$hash{$key} = $value;
	# 22 mer
	$top_seq = substr($locus_for_seq,$offset,22);
	$top_start = $start + $offset;
	$top_stop = $start + $offset + 22 - 1;
	$top_loc = "$chr" . ":" . "$top_start" . "-" . "$top_stop";
	$key = "top_$phase" . "_22";
	$value = "$top_loc\t$top_seq";
	$hash{$key} = $value;
	
	# bottom strand now
	$offset = 2 + (21 * $phase);
	# 21 mer
	$bot_seq = revcomp(substr($locus_for_seq,$offset,21));
	$bot_start = $start + $offset;
	$bot_stop = $start + $offset + 21 - 1;
	$bot_loc = "$chr" . ":" . "$bot_start" . "-" . "$bot_stop";
	$key = "bot_$phase" . "_21";
	$value = "$bot_loc\t$bot_seq";
	$hash{$key} = $value;
	# 20 mer
	my $offsetx = $offset + 1;
	$bot_seq = revcomp(substr($locus_for_seq,$offsetx,20));
	$bot_start = $start + $offsetx;
	$bot_stop = $start + $offsetx + 20 - 1;
	$bot_loc = "$chr" . ":" . "$bot_start" . "-" . "$bot_stop";
	$key = "bot_$phase" . "_20";
	$value = "$bot_loc\t$bot_seq";
	$hash{$key} = $value;
	# 22 mer
	$offsetx = $offset - 1;
	$bot_seq = revcomp(substr($locus_for_seq,$offsetx,22));
	$bot_start = $start + $offsetx;
	$bot_stop = $start + $offsetx + 22 - 1;
	$bot_loc = "$chr" . ":" . "$bot_start" . "-" . "$bot_stop";
	$key = "bot_$phase" . "_22";
	$value = "$bot_loc\t$bot_seq";
	$hash{$key} = $value;
    }
    
    # test
    #while((my $x, my $y) = each %hash) {
	#print "$x\t$y\n";
    #}
    #exit;
    
    return %hash;
}

sub revcomp {
    my($sense) = @_;
    my $revcomp = reverse $sense;
    $revcomp =~ tr/ATCG/TAGC/;
    return $revcomp;
}

sub get_mir_mod {
    my($chr,$start,$stop,$strand,$type) = @_;
    my @fields = split (":",$type);
    if($strand eq "+") {
	$start = $start + $fields[0];
	$stop = $stop + $fields[1];
    } elsif ($strand eq "-") {
	$start = $start + $fields[1];
	$stop = $stop + $fields[0];
    }
    my $out = "$chr" . ":" . "$start" . "-" . "$stop";
    return $out;
}

sub get_mir_perfects {
    my($locus,$for_seq,$strand,$mir,$star) = @_;
    my %hash = ();
    my $mod;
    my $modtype;
    my $modseq;
    my $mir_chr;
    my $mir_start;
    my $mir_stop;
    my $key;
    # miRNA and variants
    if($mir =~ /^(\S+):(\d+)-(\d+)$/) {
	$mir_chr = $1;
	$mir_start = $2;
	$mir_stop = $3;
    } else {
	die "FATAL in sub-routine get_mir_perfects .. failed to parse location $mir\n";
    }
    $modtype = "0:0";
    $mod = get_mir_mod($mir_chr,$mir_start,$mir_stop,$strand,$modtype);
    $modseq = get_mir_modseq($mod,$locus,$for_seq,$strand);
    $key = "miR" . "$modtype";
    $hash{$key} = "$mod\t$modseq";

    $modtype = "0:-1";
    $mod = get_mir_mod($mir_chr,$mir_start,$mir_stop,$strand,$modtype);
    $modseq = get_mir_modseq($mod,$locus,$for_seq,$strand);
    $key = "miR" . "$modtype";
    $hash{$key} = "$mod\t$modseq";

    $modtype = "0:-2";
    $mod = get_mir_mod($mir_chr,$mir_start,$mir_stop,$strand,$modtype);
    $modseq = get_mir_modseq($mod,$locus,$for_seq,$strand);
    $key = "miR" . "$modtype";
    $hash{$key} = "$mod\t$modseq";

    $modtype = "1:1";
    $mod = get_mir_mod($mir_chr,$mir_start,$mir_stop,$strand,$modtype);
    $modseq = get_mir_modseq($mod,$locus,$for_seq,$strand);
    $key = "miR" . "$modtype";
    $hash{$key} = "$mod\t$modseq";

    $modtype = "1:0";
    $mod = get_mir_mod($mir_chr,$mir_start,$mir_stop,$strand,$modtype);
    $modseq = get_mir_modseq($mod,$locus,$for_seq,$strand);
    $key = "miR" . "$modtype";
    $hash{$key} = "$mod\t$modseq";

    $modtype = "-1:-1";
    $mod = get_mir_mod($mir_chr,$mir_start,$mir_stop,$strand,$modtype);
    $modseq = get_mir_modseq($mod,$locus,$for_seq,$strand);
    $key = "miR" . "$modtype";
    $hash{$key} = "$mod\t$modseq";

    $modtype = "-1:-2";
    $mod = get_mir_mod($mir_chr,$mir_start,$mir_stop,$strand,$modtype);
    $modseq = get_mir_modseq($mod,$locus,$for_seq,$strand);
    $key = "miR" . "$modtype";
    $hash{$key} = "$mod\t$modseq";

    $modtype = "-2:-2";
    $mod = get_mir_mod($mir_chr,$mir_start,$mir_stop,$strand,$modtype);
    $modseq = get_mir_modseq($mod,$locus,$for_seq,$strand);
    $key = "miR" . "$modtype";
    $hash{$key} = "$mod\t$modseq";

    $modtype = "-2:-3";
    $mod = get_mir_mod($mir_chr,$mir_start,$mir_stop,$strand,$modtype);
    $modseq = get_mir_modseq($mod,$locus,$for_seq,$strand);
    $key = "miR" . "$modtype";
    $hash{$key} = "$mod\t$modseq";
    
    ## now the star variants
    my $star_chr;
    my $star_start;
    my $star_stop;
    
    if($star =~ /^(\S+):(\d+)-(\d+)$/) {
	$star_chr = $1;
	$star_start = $2;
	$star_stop = $3;
    } else {
	die "FATAL in sub-routine get_mir_perfects .. failed to parse location $star\n";
    }
    
    $modtype = "0:0";
    $mod = get_mir_mod($star_chr,$star_start,$star_stop,$strand,$modtype);
    $modseq = get_mir_modseq($mod,$locus,$for_seq,$strand);
    $key = "star" . "$modtype";
    $hash{$key} = "$mod\t$modseq";

    $modtype = "0:-1";
    $mod = get_mir_mod($star_chr,$star_start,$star_stop,$strand,$modtype);
    $modseq = get_mir_modseq($mod,$locus,$for_seq,$strand);
    $key = "star" . "$modtype";
    $hash{$key} = "$mod\t$modseq";

    $modtype = "0:-2";
    $mod = get_mir_mod($star_chr,$star_start,$star_stop,$strand,$modtype);
    $modseq = get_mir_modseq($mod,$locus,$for_seq,$strand);
    $key = "star" . "$modtype";
    $hash{$key} = "$mod\t$modseq";

    $modtype = "1:1";
    $mod = get_mir_mod($star_chr,$star_start,$star_stop,$strand,$modtype);
    $modseq = get_mir_modseq($mod,$locus,$for_seq,$strand);
    $key = "star" . "$modtype";
    $hash{$key} = "$mod\t$modseq";

    $modtype = "1:0";
    $mod = get_mir_mod($star_chr,$star_start,$star_stop,$strand,$modtype);
    $modseq = get_mir_modseq($mod,$locus,$for_seq,$strand);
    $key = "star" . "$modtype";
    $hash{$key} = "$mod\t$modseq";

    $modtype = "-1:-1";
    $mod = get_mir_mod($star_chr,$star_start,$star_stop,$strand,$modtype);
    $modseq = get_mir_modseq($mod,$locus,$for_seq,$strand);
    $key = "star" . "$modtype";
    $hash{$key} = "$mod\t$modseq";
    
    # test
    #while((my $x, my $y) = each %hash) {
	#print "$x\t$y\n";
    #}
    #exit;
    
    return %hash;
    
}

sub get_mir_modseq {
    my($mod,$locus,$for_seq,$strand) = @_;
    # parse
    my $loc_start;
    my $loc_stop;
    if($locus =~ /^\S+:(\d+)-(\d+)$/) {
	$loc_start = $1;
	$loc_stop = $2;
    } else {
	die "FATAL in sub-routine get_mir_modseq .. could not parse location $locus\n";
    }
    my $mod_start;
    my $mod_stop;
    if($mod =~ /^\S+:(\d+)-(\d+)$/) {
	$mod_start = $1;
	$mod_stop = $2;
    } else {
	die "FATAL in sub-routine get_mir_modseq .. could not parse location $mod\n";
    }
    my $offset = $mod_start - $loc_start; ## zero-based
    my $length = $mod_stop - $mod_start + 1;
    my $modseq_f = substr($for_seq,$offset,$length);
    my $modseq;
    if($strand eq "-") {
	$modseq = reverse $modseq_f;
	$modseq =~ tr/ATGC/TACG/;
    } else {
	$modseq = $modseq_f;
    }
    return $modseq;
}
    
    
sub errorify_read {
    my($perfect) = @_;
    my $pick = rand();
    my $final;
    my $errors = 0;
    if($pick < $opt_e) {
	# pick a random location
	my $change_pos = int(rand(length $perfect));
	my @bases = split ('', $perfect);
	my $new_base = get_new_base($bases[$change_pos]);
	$bases[$change_pos] = $new_base;
	$final = join('', @bases);
	$errors = 1;
    } else {
	$final = $perfect;
    }
    return ($final,$errors);
}

sub get_new_base {
    my($old_base) = @_;
    my $new_base = "N";
    until (($new_base ne "N") and ($new_base ne "$old_base")) {
	my $pick = int(rand(4));
	if($pick == 0) {
	    $new_base = "A";
	} elsif ($pick == 1) {
	    $new_base = "T";
	} elsif ($pick == 2) {
	    $new_base = "G";
	} elsif ($pick == 3) {
	    $new_base = "C";
	}
    }
    return $new_base;
}
    

sub sample_mir_keys {
    my $pick = rand();
    my $key;
    if($pick < 0.6) {
	$key = "miR0:0";
    } elsif (($pick >= 0.6) and ($pick < 0.8)) {
	$key = "star0:0";
    } elsif (($pick >= 0.8) and ($pick < 0.84)) {
	$key = "miR0:-1";
    } elsif (($pick >= 0.84) and ($pick < 0.85)) {
	$key = "miR0:-2";
    } elsif (($pick >= 0.85) and ($pick < 0.89)) {
	$key = "miR1:1";
    } elsif (($pick >= 0.89) and ($pick < 0.9)) {
	$key = "miR1:0";
    } elsif (($pick >= 0.9) and ($pick < 0.92)) {
	$key = "miR-1:-1";
    } elsif (($pick >= 0.92) and ($pick < 0.925)) {
	$key = "miR-1:-2";
    } elsif (($pick >= 0.925) and ($pick < 0.945)) {
	$key = "miR-2:-2";
    } elsif (($pick >= 0.945) and ($pick < 0.95)) {
	$key = "miR-2:-3";
    } elsif (($pick >= 0.95) and ($pick < 0.96)) {
	$key = "star0:-1";
    } elsif (($pick >= 0.96) and ($pick < 0.965)) {
	$key = "star0:-2";
    } elsif (($pick >= 0.965) and ($pick < 0.975)) {
	$key = "star1:1";
    } elsif (($pick >= 0.975) and ($pick < 0.98)) {
	$key = "star1:0";
    } elsif ($pick >= 0.98) {
	$key = "star-1:-1";
    }
    return $key;
}

sub pick_an_arm {
    my($locus) = @_;
    my $chr;
    my $loc_start;
    my $loc_stop;
    if($locus =~ /^(\S+):(\d+)-(\d+)$/) {
	$chr = $1;
	$loc_start = $2;
	$loc_stop = $3;
    } else {
	die "FATAL: in sub-routine pick_an_arm : failed to parse locus name $locus\n";
    }
    my $left_start = $loc_start + 17;
    my $left_stop = $loc_start + 17 + 20;
    my $right_start = $loc_start + 17 + 20 + 48;
    my $right_stop = $loc_start + 17 + 20 + 48 + 20;
    
    my $pick = rand();
    my $mir;
    my $star;
    if($pick >= 0.5) {
	$mir = "$chr" . ":" . "$right_start" . "-" . "$right_stop";
	$star = "$chr" . ":" . "$left_start" . "-" . "$left_stop";
    } else {
	$star = "$chr" . ":" . "$right_start" . "-" . "$right_stop";
	$mir = "$chr" . ":" . "$left_start" . "-" . "$left_stop";
    }
    return ($mir,$star);
}
    
sub pick_a_strand {
    my $call = rand();
    my $strand;
    if($call >= 0.5) {
	$strand = "-";
    } else {
	$strand = "+";
    }
    return $strand;
}



sub validate_genome {
    if($opt_s) {
	# soft-masked only
	# is genome file there?
	unless(-r $opt_s) {
	    print STDERR "FATAL: Failed to find genome file $opt_s\n";
	    return 0;
	}
	my $s_fai = "$opt_s" . ".fai";
	if(-r $s_fai) {
	    return 1;
	} else {
	    print STDERR "Failed to open expected fai index file $s_fai.\nAttempting to create using samtools faidx ...";
	    system "samtools faidx $opt_s";
	    if(-r $s_fai) {
		print STDERR " Successful.\n";
		return 1;
	    } else {
		print STDERR " FAILED, FATAL.\n";
		return 0;
	    }
	}
    } elsif (($opt_d) and ($opt_n)) {
	unless(-r $opt_d) {
	    print STDERR "FATAL: Failed to find genome file $opt_d\n";
	    return 0;
	}
	unless(-r $opt_n) {
	    print STDERR "FATAL: Failed to find genome file $opt_n\n";
	    return 0;
	}
	my $un_fai = "$opt_n" . ".fai";
	my $m_fai = "$opt_d" . ".fai";
	# check -r
	unless(-r $un_fai) {
	    print STDERR "Failed to open expected fai index file $un_fai.\nAttempting to create using samtools faidx ... \n";
	    system "samtools faidx $opt_n";
	    if(-r $un_fai) {
		print STDERR "Successful.\n";
	    } else {
		print STDERR "FAILED. FATAL.\n";
		return 0;
	    }
	}
	unless(-r $m_fai) {
	    print STDERR "Failed to open expected fai index file $m_fai.\nAttempting to create using samtools faidx ... \n";
	    system "samtools faidx $opt_d";
	    if(-r $m_fai) {
		print STDERR "Successful.\n";
	    } else {
		print STDERR "FAILED. FATAL.\n";
		return 0;
	    }
	}
	
	open(UN, "$un_fai");
	open(M, "$m_fai");
	my %unm = ();
	my @fields = ();
	while (<UN>) {
	    chomp;
	    @fields = split ("\t", $_);
	    $unm{$fields[0]} = $fields[1];
	}
	close UN;
	while (<M>) {
	    chomp;
	    @fields = split ("\t", $_);
	    if(exists($unm{$fields[0]})) {
		unless($unm{$fields[0]} == $fields[1]) {
		    close M;
		    print STDERR "The genomes $opt_d and $opt_n are NOT identical according to their .fai indices. Cannot proceed.\n";
		    return 0;
		}
	    } else {
		close M;
		print STDERR "The genomes $opt_d and $opt_n are NOT identical according to their .fai indices. Cannot proceed.\n";
		return 0;
	    }
	}
	close M;
	return 1;
    }
}

sub get_genome_positions {
    my $fai;
    if($opt_s) {
	$fai = "$opt_s" . ".fai";
    } else {
	$fai = "$opt_d" . ".fai";
    }
    open(FAI, "$fai");
    my %hash = ();
    my $start;
    my $stop = 0;
    my @fields = ();
    while (<FAI>) {
	chomp;
	@fields = split ("\t", $_);
	$start = 1 + $stop;
	$stop = $start + $fields[1] - 1;
	$hash{$fields[0]} = "$start\t$stop";
	# test
	#print "$fields[0]\t$start\t$stop\n";
    }
    close FAI;
    return %hash;
}

sub get_gp_max { ## by reference, a hash
    my($hash) = @_;
    my $max = 0;
    my $key;
    my $tab;
    my @fields = ();
    while(($key,$tab) = each %$hash) {
	@fields = split ("\t", $tab);
	if($fields[1] > $max) {
	    $max = $fields[1];
	}
    }
    return $max;
}

sub get_a_locus {
    my($occupied,$genome_positions,$gp_max,$locus_size,$genome_type) = @_; ## by reference. hash, hash, scalar, scalar, scalar
    my $ok = 0;
    my $rand_loc;
    my $stop;
    my $chr;
    my $chr_start;
    my $chr_stop;
    my $tab;
    my @fields = ();
    my $true_start;
    my $true_stop;
    my $maybe;
    my $raw_seq;
    my @occ_pairs = ();
    until ($ok == 1) {
	# test
	#print "\nVariable ok is $ok at start of loop\n";
	
	# randomly select a genomic position
	$rand_loc = int(rand(($$gp_max-$$locus_size)));
	$stop = $$locus_size + $rand_loc - 1;
	
	# test 
	#print STDERR "\nrand_loc: $rand_loc stop: $stop\n";
	
	# does it fall off an end?
	while(($chr,$tab) = each %$genome_positions) {
	    @fields = split ("\t", $tab);
	    if(($rand_loc >= $fields[0]) and
	       ($rand_loc <= $fields[1])) {
		$chr_start = $chr;
	    }
	    if(($stop >= $fields[0]) and 
	       ($stop <= $fields[1])) {
		$chr_stop = $chr;
	    }
	}
	
        # test
	#print STDERR "\tchr_start: $chr_start chr_stop $chr_stop\n";
	
	unless($chr_start eq $chr_stop) {
	    next;
	}
	# convert to genome coordinates
	@fields = split("\t", $$genome_positions{$chr_start});
	$true_start = $rand_loc - $fields[0] + 1;
	$true_stop = $stop - $fields[0] + 1;
	$maybe = "$chr_start" . ":" . "$true_start" . "-" . "$true_stop";
	
	# test
	#print STDERR "\tmaybe: $maybe\n";
	
	# get the sequence, confirm it is OK
	if($opt_s) {
	    (open(FASTA, "samtools faidx $opt_s $maybe |")) || return 0;
	} elsif ($$genome_type eq "masked") {
	    (open(FASTA, "samtools faidx $opt_d $maybe |")) || return 0;
	} elsif ($$genome_type eq "notmasked") {
	    (open(FASTA, "samtools faidx $opt_n $maybe |")) || return 0;
	}
	$raw_seq = '';
	while (<FASTA>) {
	    chomp;
	    unless ($_ =~ /^>/) {
		$raw_seq .= $_;
	    }
	}
	close FASTA;
	unless(($opt_s) and ($$genome_type eq "masked")) {
	    # upper-case it, unless we have a soft-masked genome AND we want to avoid the masked parts
	    $raw_seq = uc $raw_seq;
	}
	
	# test
	#print STDERR "\traw_seq: $raw_seq\n";
	
	unless($raw_seq =~ /^[ATGC]+$/) {
	    # test
	    # print STDERR "\tFailed regex\n";
	    next;
	}
	
	# does it overlap with a previously accepted locus?
	my $fail = 0;
	if(exists($$occupied{$chr_start})) {
	    @occ_pairs = @{$$occupied{$chr_start}};
	    for(my $i = 0; $i < (scalar @occ_pairs); $i+=2) {
		if((($true_start >= $occ_pairs[$i]) and ($true_start <= $occ_pairs[($i+1)])) or
		   (($true_stop >= $occ_pairs[$i]) and ($true_stop <= $occ_pairs[($i+1)]))) {
		    $fail = 1;
		}
	    }
	}
	if($fail) {
	    # test
	    #print STDERR "\tFail after comparison to occupied .. variable fail in state $fail\n";
	    next;
	}
	
	$ok = 1;
	# test
	#print STDERR "\tvariable ok is $ok and loop should stop\n";
    }
    # add it to the occupied hash
    push(@{$occupied{$chr_start}}, $true_start);
    push(@{$occupied{$chr_start}}, $true_stop);

    # test
    #exit;
    
    # return the locus coordinates
    return ($maybe,$raw_seq);
    
}


sub get_nts {
    my($reads2get,$entries2get) = @_;
    #my $max = int(0.25 * $reads2get);
    my $logmax = log($reads2get) / log(10);
    my $interval = $logmax / $entries2get;
    my $reads_done = 0;
    my @out = ();
    my $this_log = 0;
    my $reads;
    until ($reads_done >= $reads2get) {
	$reads = 1 + int(10**$this_log);
	# exceeded?
	if(($reads_done + $reads) > $reads2get) {
	    $reads = $reads2get - $reads_done;
	}
	push(@out, $reads);
	$this_log += $interval;
	$reads_done += $reads;
	
    }
    # test
    #print "\nreads_done: $reads_done\n";

    return @out;
}
	
__END__

=head1 SYNOPSIS

sim_sRNA-seq.pl

Simulation of plant small RNA-seq data

Copyright (C) 2014 Michael J. Axtell

=head1 AUTHOR

Michael J. Axtell, Penn State University, mja18@psu.edu

=head1 LICENSE

This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 DEPENDENCIES

perl 5 <www.perl.org> at /usr/bin/perl

samtools <http://samtools.sourceforge.net/> in the PATH

=head1 INSTALL

Ensure dependicies are present (see above). Add the script sim_sRNA-seq.pl to your PATH

=head1 USAGE

    sim_sRNA-seq.pl [options] -s soft-masked-genome.fa                                                                                        

--- OR, LESS PREFERRED  ---                                                                                                                       

    sim_sRNA-seq.pl [options] -d hard-masked-genome.fa -n not-masked-genome.fa      

=head1 OPTIONS

Options:

-s : path to soft-masked genome. Lower-case letters are assumed repeat-masked. If -s is specified, neither -d nor -n can be specified.

-d : path to hard-masked genome. If -d is specified, -n must also be specified, and -s cannot be specified.

-n : path to not-masked genome. If -n is specified, -d must also be specified, and -s cannot be specified.

-v : print version number and quit.

-h : print help message and quit.

-r : desired number of reads, in millions. Default: 5.

-e : per-read probability of a single nt sequencing error (substitution). Default: 1E-4.

=head1 METHODS

=head2 Read numbers and types

30% of the simulated reads will come from roughly 100 MIRNA loci, 5% of the simulated reads will come from roughly 20 tasiRNA/phased secondary siRNA loci, and the remaining 65% from roughly 10,000 heterochromatic siRNA loci.

The abundance of reads from each locus is distributed on a log-linear scale (e.g., plotting the log10 of read number as a function of abundance rank yields a straight line).

=head2 MIRNA simulation

Valid MIRNA loci have no overlap with repeat-masked regions of the genome. The locus has a hypothetical size of 125nts, of which all must be ATGC bases (e.g. no ambiguous nts.).

The strand of the MIRNA precursor is randomly selected, as is the arm from which the mature miRNA and star come from. There is no actual hairpin sequence necessarily present at MIRNA loci .. it is only the pattern of reads that is being simulated.

The mature miRNA and mature miRNA* are defined as 'master' positions. The left-most 'master' position in a locus is a 21-mer starting at position 17 of the locus. The right-most 'master' position is a 21 mer at position 85. Assignment of the arm (i.e., whether the left-most or right-most 'master' position is the miRNA or miRNA*) is random at each locus.

Once a locus has been found and 'master' positions defined, each read is simulated according to the following probabilities. In the following list, "miR" means mature miRNA, "star" means miRNA*. The numbers after each indicate the offset at 5' and 3' ends relative to the master positions. So, "miR0:1" means the mature miRNA sequence, starting at the master 5' end, and ending 1 nt after the master 3' end.

60% miR0:0, 20% star0:0,, 4% miR0:-1, 1% miR0:-2, 4% miR1:1, 1% miR1:0, 2% miR-1:-1, 0.5% miR-1:-2, 2% miR-2:-2, 0.5% miR-2:-3, 1% star0:-1, 0.5% star0:-2, 1% star1:1, 0.5% star1:0, 2% star-1:-1.

Sampling of simulated MIRNA-derived reads continues until the required number of reads for a particular locus is recovered.

=head2 tasiRNA/phased siRNA simulation

TAS loci have no overlap with repeat-masked areas of the genome. Each locus has a nominal size of 140nts. All bases must be ATGC (non-ambiguous).

Each locus is simulated to be diced in 6 21 nt phases. At each phasing position, 21mers are the dominant size, with 20mers and 22mers being less frequent .. the 20 and 22nt variants vary in their 3' positions relative to the 'master' 21nt RNAs.

Once a locus has been identified, and all possible 20, 21, and 22mers charted, each read is simulated according to the following probabilities:

Strand of origin is 50% top, 50% bottom.

Phase position is equal chance for all (e.g. 1/6 chance for any particular phase location).

80% of the time, the 21mer is returned, 10% of the time the 20mer, and 10% the 22 mer.

=head2 Heterochromatic siRNA simulation

Heterochromatic siRNA loci can come from anywhere in the genome (repeat-masked or un-masked), and their nominal locus size is 100nts. All nts in the 100 nt locus must be ATGC (no ambiguous positions).

All possible start and stop positions for 24, 23, 22, and 21 nt RNAs are eligible from these loci. Each simulated read is identified with the following probabilities:

50% chance for top of bottom strand origin

5' position within the locus is random

90% chance of a 24 mer, 5% chance of a 23 mer, 3% chance of a 22 mer, and 2% chance of a 21 mer.

=head2 Simulation of sequencing errors

Regardless of small RNA type, once a given read has been simulated, there is a chance of introucing a single-nucleotide substitution relative to the reference genome at a randomly selected position in the read. The chance of introducing this error is given by option -e (Default is 1 in 10,000).

=head2 No overlapping loci

None of the simulated loci are allowed to have any overlap with each other.

=head1 OUTPUT

=head2 Files

Two files are created in the working directory.

simulated_sRNA-seq_reads_overview.txt : A tab-delimited text file giving the coordinates and names for each simulated locus.

simulated_sRNA-seq_reads.fasta : A FASTA file of all of the simulated reads.

=head2 Naming conventions

The FASTA header of each simulated read encode the basic information about the read. Several different fields are separated by "_" characters. For instance the following read header ">MIRNA_2_2_chr7:123063894-123063914_+_0" means ..

MIRNA: This came from a simulated MIRNA locus

2: The 2nd simulated locus of this type

2: Read number 2 from this locus

chr7:123063894-123063914: The true origin of this read.

+: The genomic strand of this read.

0: The number of sequencing errors simulated into the read.

