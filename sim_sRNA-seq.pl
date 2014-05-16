#!/usr/bin/perl -w
use strict;
use Getopt::Std;

my $ver_num = "dev";

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
    print "mir_locus: $mir_locus mir_count: $mir_count mir_strand: $mir_strand mature_mir: $mature_mir star_mir: $star_mir\n";
    
    ++$mir_locus_n;
    print OUTS "MIRNA_$mir_locus_n\t$mir_locus\t$mir_strand\t$mir_count\n";
	
    # simulate expression
    for(my $i = 1; $i <= $mir_count; ++$i) {
	$simulated_read_location = simulate_mirna($mature_mir,$star_mir,$mir_strand);
	
	($simulated_read,$simulated_read_errors) = errorify_read($simulated_read_location,$mir_strand);
	
	# test
	#print "\tsimulated_read_location: $simulated_read_location simuated_read: $simulated_read simulated_read_errors: $simulated_read_errors\n";
	$fasta_header = ">MIRNA_$mir_locus_n" . "__" . "$i" . "__" . "$simulated_read_location" . "__" . "$mir_strand" . "__" . "$simulated_read_errors";
	print OUTF "$fasta_header\n$simulated_read\n";
    }
}



########
sub errorify_read {
    my($location,$strand) = @_;
    # get perfect read sequence
    if($opt_s) {
	(open(FASTA, "samtools faidx $opt_s $location |")) || die "FATAL in sub-routine errorify read. Failed to open FASTA with samtools faidx\n";
    } else {
	# always can use the unmasked at this point
	(open(FASTA, "samtools faidx $opt_n $location |")) || die "FATAL in sub-routine errorify read. Failed to open FASTA with samtools faidx\n";
    }
    my $for_perfect_seq;
    while (<FASTA>) {
	chomp;
	unless($_ =~ /^>/) {
	    $for_perfect_seq .= uc $_;
	}
    }
    close FASTA;
    # revcomp if needed
    my $perfect;
    if($strand eq "-") {
	$perfect = reverse $for_perfect_seq;
	$perfect =~ tr/ATCG/TAGC/;
    } else {
	$perfect = $for_perfect_seq;
    }
    # errorify
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
    

sub simulate_mirna {
    my($mir_loc,$star_loc,$mir_strand) = @_;
    my $pick = rand();
    my $chr;
    my $start;
    my $stop;
    my $mod;
    my $modtype;
    if($pick < 0.6) {
	return $mir_loc;
    } elsif (($pick >= 0.6) and ($pick < 0.8)) {
	return $star_loc;
    } elsif (($pick >= 0.8) and ($pick < 0.95)) {
	# a variant of the mature miRNA
	# must parse the mir_loc
	if($mir_loc =~ /^(\S+):(\d+)-(\d+)$/) {
	    $chr = $1;
	    $start = $2;
	    $stop = $3;
	} else {
	    die "FATAL in sub-routine simulate_mirna .. failed to parse location $mir_loc\n";
	}
	if(($pick >= 0.8) and ($pick < 0.84)) {
	    # 0:-1
	    $modtype = "0:-1";
	} elsif (($pick >= 0.84) and ($pick < 0.85)) {
	    # 0:-2
	    $modtype = "0:-2";
	} elsif (($pick >= 0.85) and ($pick < 0.89)) {
	    # 1:1
	    $modtype = "1:1";
	} elsif (($pick >= 0.89) and ($pick < 0.9)) {
	    # 1:0
	    $modtype = "1:0";
	} elsif (($pick >= 0.9) and ($pick < 0.92)) {
	    # -1:-1
	    $modtype = "-1:-1";
	} elsif (($pick >= 0.92) and ($pick < 0.925)) {
	    # -1:-2
	    $modtype = "-1:-2";
	} elsif (($pick >= 0.925) and ($pick < 0.945)) {
	    # -2:-2
	    $modtype = "-2:-2";
	} elsif (($pick >= 0.945) and ($pick < 0.95)) {
	    # -2:-3
	    $modtype = "-2:-3";
	}
	$mod = get_mir_mod($chr,$start,$stop,$mir_strand,$modtype);
	return $mod;
    } elsif ($pick >= 0.95) {
	# a variant of the mir-star
	# must parse the star_loc
	if($star_loc =~ /^(\S+):(\d+)-(\d+)$/) {
	    $chr = $1;
	    $start = $2;
	    $stop = $3;
	} else {
	    die "FATAL in sub-routine simulate_mirna .. failed to parse location $star_loc\n";
	}
	if(($pick >= 0.95) and ($pick < 0.96)) {
	    # 0:-1
	    $modtype = "0:-1";
	} elsif (($pick >= 0.96) and ($pick < 0.965)) {
	    # 0:-2
	    $modtype = "0:-2";
	} elsif (($pick >= 0.965) and ($pick < 0.975)) {
	    # 1:1
	    $modtype = "1:1";
	} elsif (($pick >= 0.975) and ($pick < 0.98)) {
	    # 1:0
	    $modtype = "1:0";
	} elsif ($pick >= 0.98) {
	    # -1:-1
	    $modtype = "-1:-1";
	}
	$mod = get_mir_mod($chr,$start,$stop,$mir_strand,$modtype);
	return $mod;
    }
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
	


    

