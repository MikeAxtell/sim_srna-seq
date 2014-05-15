#!/usr/bin/perl -w
use strict;
use Getopt::Std;

my $ver_num = "dev";

my $usage = "sim_sRNA-seq.pl version: $ver_num

Simulate plant small RNA-seq data

Usage: sim_sRNA-seq.pl [options] rep-masked-genome.fasta unmasked-genome.fasta

Options:

-v : print version number and quit
-h : print help message and quit
-r : desired total number of reads, in millions. Default: 5
-e : per-nt sequencing error rate. Default: 1E-4

";

# initialize defaults

our $opt_r = 5;
our $opt_e = 0.0001;
our $opt_h;
our $opt_v;

# get options and validate them
getopt('rehv');

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

# verify genome files
my $unmasked = pop @ARGV;
my $masked = pop @ARGV;
unless(-r $unmasked) {
    print "FATAL: Could not read unmasked genome file $unmasked\n$usage\n";
}
unless(-r $masked) {
    print "FATAL: Could not read masked genome file $masked\n$usage\n";
}

# validate genomes.
my $validated_genomes = validate_genomes($unmasked,$masked);
if($validated_genomes == -1 ) {
    print "FATAL: Genomes don't match according to their .fai files. Genomes must have the same chromosome names and lengths.\n$usage\n";
    exit;
} elsif ($validated_genomes == 0) {
    print "FATAL: .fai files for one or both of the genomes were not found. Please index the two genome versions using samtools faidx.\n$usage\n";
    exit;
} 

# build a hash of one-based cumulative positions for nts in the genome, based on order in the fai file
my %genome_positions = get_genome_positions($unmasked);
my $gp_max = get_gp_max(\%genome_positions);

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

# hash for occupied regions (those already selected)
# structure: keys: chromosome name, value: anonymous array of integers, where each pair is a coordinate pair. 
# this list is NOT sorted in any way
my %occupied = ();
## miRNA simulation
foreach my $mir_count (@mir_ns) {
    my $mir_locus = get_mir_locus(\$masked,\%occupied);  ## pass the masked genome, hash of occupied regions
}


########

sub validate_genomes {
    my($unmasked, $masked) = @_;
    my $un_fai = "$unmasked" . ".fai";
    my $m_fai = "$masked" . ".fai";
    (open(UN, "$un_fai")) || return 0;
    (open(M, "$m_fai")) || return 0;
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
		return -1;
	    }
	} else {
	    close M;
	    return -1;
	}
    }
    close M;
    return 1;
}

sub get_genome_positions {
    my($genome) = @_;
    my $fai = "$genome" . ".fai";
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

sub get_mir_locus {
    my($masked,$occupied) = @_; ## by reference ... scalar and hash
    open(MASKED, "$$masked");
    
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
	


    

