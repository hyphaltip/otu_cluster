#!/usr/bin/perl -w
# $Id: needle_clusters.pl 1 2008-12-03 22:42:29Z stajich $
# Usage: needle_clusters.pl report.needle datafile.fas
# Single-Linkage clustering based on global percent identity
# Used in O'Brien et al, AEM. doi:10.1128/AEM.71.9.5544-5550.2005
# Author: Jason Stajich <jason[_AT_]open-bio.org>
# http://www.duke.edu/~jes12/

# this requires the all-vs-all needleman-wunsch report should have been already run
# see run_do_needle.sh for script to run all the pairwise alignments using
# needle (part of the EMBOSS package) http://www.emboss.org/

use strict;
use Getopt::Long;


my $Usage = <<EOF
needle_clusters.pl [options] <emboss report> <sequences>
    Options
    -h     : show brief help
    -i NUM : Minimum Percent Identity cutoff to use when 
             building single-linkage clusters (a num between 0 and 100) 
             [default 98]
    -s file: a file to provide a simple report of the sequences and
             percent identity of the alignments [default reportname.sum]
    -o file: a file to save the output of the clusters [default is STDOUT] 
            (with specifying the -o option you can run this script as
             needle_clusters.pl ... > clusters.out)
    -p file: a file to save the cluster stats wrt to the plots the seqs
             came from. 
    -w NUM : fasta sequence block width [default 60]
EOF
    ;
							

my ($minid,$Width,$output_sum, $output_clusters,$output_plots) = (98,60);
GetOptions(
	   'h|help'      => sub { print STDERR $Usage,"\n";
				  exit(0) },
	   'i|id:i'      => \$minid,
	   'o|output:s'  => \$output_clusters,
           'p|plotfile:s'=> \$output_plots,
	   's|summary:s' => \$output_sum,
	   'w|width:i'   => \$Width,
	   );

if( $minid > 100 || $minid < 0 ) { 
    die("percent id (-i option) but be between 0 and 100\n$Usage");
}

my ($reportfile,$seqfile) = @ARGV;  # get the cmdline remaining arguments
				    # from ARGV

if( ! defined $reportfile || ! defined $seqfile ) { 
    print STDERR $Usage, "\n";
    exit(0);
}

unless( defined $output_sum ) { 
    $output_sum = $reportfile.".sum";
} 

my %seqs       = &read_fasta($seqfile);
my %distmat    = &read_emboss_identity($reportfile);

my @clusters   = &single_linkage(%identities,$minid);

my $cid = 1;
open(SUM, ">$output_sum") || die("cannot open summary outputfile ($output_sum): $!");
my $clusterout;
unless( defined $output_clusters ) {
    $output_clusters = $reportfile.".clusters";
}

open($clusterout, ">$output_clusters") || 
    die("cannot open cluster outputfile ($output_clusters): $!");

my %plots;
while( my ($sequence1,$data1) = each %identities ) {
    if( defined $output_plots ) {
	if( $sequence1 =~ /^(dfmo(\d+)|(\S+))\_/ ) {
	    $plots{$2 || $3}++;
	}
    }
    while( my ($sequence2,$pid) = each %$data1 ) {
	print SUM join("\t",$sequence1, $sequence2, $pid), "\n";
    }
}

my $clusterplotfh;
my @dfmonums;
if( defined $output_plots ) {
    open($clusterplotfh, ">$output_plots") || die($!);
    # deal with either having numbers or letters, need 
    # to sort accordingly
    @dfmonums = sort { &isnum($a) && &isnum($b) ? $a <=> $b : $a cmp $b } 
    keys %plots;
    print $clusterplotfh join("\t", "CLUSTER#", @dfmonums), "\n";
}
foreach my $c ( @clusters ) {
    printf $clusterout "--- Begin cluster %d, it has %d members---\n", $cid, 
    scalar @$c; 
    my %plotdat;
    
    foreach my $gene ( @$c ) {
	if( defined $output_plots ) {
	    if( $gene =~ /^dfmo(\d+)\_/) {
		$plotdat{$1}++;
	    }
	}
	&write_fasta($clusterout,$gene,
		     $seqs{$gene} || die("cannot find $gene\n"));
    }
    print $clusterout "---end of cluster $cid---\n";
    if( %plotdat ) {
	print $clusterplotfh join("\t", "Cluster $cid",
				  map { $plotdat{$_} || 0 } @dfmonums), 
	"\n";
    }
    $cid++;
}

# --- Functions below ---

# parse Emboss needle or water output
# return a hash with key being the sequence id,
#         each value is in tern a hashref pointing to the sequences
#         which were aligned and the value for this $hash{key}->{key} will be
#         the percent identity 
# input is the name of the file
sub read_emboss_identity {
    my ($emboss_report,$no_skip_same) = @_;
    open(REPORT, "cat $emboss_report |") || die("cannot open file $emboss_report: $!");
    my (@seqs,%lengths,%gaps,%identities,%distmat,$first);
    $first = 1;
    while(<REPORT>) {
	if( /^>/ && $first ) { 
	    die("You probably switched the order of the sequence and report file, try again.\n");
	}
	if( /^\#\s+(1|2):\s+(\S+)/) {
	    $seqs[$1-1] = $2;
	} elsif( /^\#\s+Length:\s+(\d+)/ ) {
	    $lengths{$seqs[0]}->{$seqs[1]} = $1;
	} elsif( /^\#\s+Identity:\s+(\d+)\/(\d+)/ ) {
	    $identities{$seqs[0]}->{$seqs[1]} = $1;	    
	} elsif( /^\#\s+Gaps:\s+(\d+)\/(\d+)/ ) {
	    $gaps{$seqs[0]}->{$seqs[1]} = $1;
	}
	$first = 0;
    }
    close(REPORT);
    for my $seqa ( keys %lengths ) {
	for  my $seqb ( keys %{$lengths{$seqa}} ) {
	    $distmat{$seqa}->{$seqb} = 1 - ( $identities{$seqa}->{$seqb} / 
					     ($lengths{$seqa}->{$seqb} - 
					      $gaps{$seqa}->{$seqb}));
	}
    }
    return %distmat;
}

# will take as input the name of sequence file in fasta format
# will return a HASH of sequences where the key in the HASH is the
# sequence name and the data is the sequence string
sub read_fasta {
    # stolen in part from Bioperl Bio::SeqIO::fasta
    my ($seqfile) = @_;
    # store seqs by 
    my %seqs;    
    local $/ = "\n>"; # for the data input
    open(IN, $seqfile) || die("could not open file $seqfile: $!");    
    while(<IN>) {
	my ($top,$sequence) = split(/\n/,$_,2);
	my ($id,$fulldesc);
	if( $top =~ /^\s*(\S+)\s*(.*)/ ) {
	    ($id,$fulldesc) = ($1,$2);
	}	
        # FIX incase no space between > and name \AE
	if (defined $id && $id eq '') {$id=$fulldesc;} 
	$id =~ s/^>//;
	defined $sequence && $sequence =~ s/\s//g;
	$sequence =~ s/>$//;
	next unless defined $id;
	$seqs{$id} = $sequence;
    }
    close(IN);
    return %seqs;
}

# will output sequences in fasta format
# provide the filehandle, sequence id and the sequence string to 
# the write_fasta function for outputting sequence 
sub write_fasta {
    my ($fh,$id,$str) = @_;
    if(length($str) > 0) {
	$str =~ s/(.{1,$Width})/$1\n/g;
    } else {
	$str = "\n";
    }
    print $fh (">",$id,"\n",$str) or return;
}

# This function will take as input
# sequence scores - hash reference, where keys are the ids and value are the
#                   scores [percent identity for this script]
# cutoff          - the minimum cutoff value for putting sequences in the same 
#                   cluster 
sub single_linkage {
    my ($seqscores,$cutoff) = @_;
    my (@clusters,$clusterct,%seqclusterlookup);
    
    while( my ($sequence1,$data1) = each %$seqscores ) {
	# get the bin for this sequence, in the absence 
	# of seeing anything about this sequence we give it a
	# new bin number, we'll join bins if need be
	my $bin = $seqclusterlookup{$sequence1};
	unless( defined $bin ) {
	    $bin = $seqclusterlookup{$sequence1} = $clusterct++;
	    push @{$clusters[$bin]}, $sequence1;
	}
	while( my ($sequence2,$pid) = each %$data1 ) {
	    next if( $pid < $cutoff );
	    
	    my $bin2 = $seqclusterlookup{$sequence2};
	    
	    if( defined $bin2 ) {    # already seen a sequence2 before
		if( $bin != $bin2 ) { 
		    # join these two bins since they meet the threshold
		    # this means pulling everything from both bins
		    # together
		    push @{$clusters[$bin]}, @{$clusters[$bin2]};
		    $clusters[$bin2] = [];
		}
		next;
	    } else { 
		push @{$clusters[$bin]}, $sequence2;
	    }
	    $seqclusterlookup{$sequence2} = $bin2 = $bin;
	}
    }    
    # this only returns clusters which are > 0 in size (empty clusters
    # come from joining to bins together
    # it sorts them by size so largest cluster cluster comes first
    return sort { scalar @$b <=> scalar @$a } 
           grep { defined && scalar @$_ > 0} @clusters;
}

sub isnum { $_[0] =~ /^\d+$/ }
