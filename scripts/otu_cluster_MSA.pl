#!/usr/bin/perl -w
# $Id: otu_cluster_MSA.pl 1 2008-12-03 22:42:29Z stajich $
# generate cluster %id (or %dissimilarity) from a mutiple sequence alignment

use strict;

use Bio::AlignIO;

my %matrix;

my $needlefile = shift @ARGV;
my $in = Bio::AlignIO->new(-format => 'fasta',
			   -file   => $needlefile);
if( my $aln = $in->next_aln ) {
    my @seqs = $aln->each_seq;
    for my $q ( @seqs ) {
	$matrix{$q->id}->{$q->id} = 0;
	for my $h (@seqs) {
	    next if $h->id eq $q->id;
	    # cliques are distance so we do 1-%identity
	    my $pid = percent_id($h->seq, $q->seq);
	    $matrix{$q->id}->{$h->id} = sprintf("%.8f", 1-$pid);
	}
    }
} else {
    warn("no MSA (multi-fasta alignment) provided as $needlefile\n");
    exit;
}

my @names = sort keys %matrix;
print scalar @names, "\n"; # print out the 1st line of DOTUR matrix which is
                           # number of names
for my $name ( @names ) {
    my @l;
    for my $n ( @names ) {
	push @l, $matrix{$name}->{$n};
    }
    print join("\t", $name, @l),"\n";
}

sub percent_id {
    my ($qseq,$hseq) = @_;
    if( length($qseq) ne length($hseq) ) {
	warn("cannot calculate PID if they aren't aligned!\n");
	return -1;
    }
    my $len = length($qseq);
    my $gaps = 0;
    my $match = 0;
    my $mismatch = 0;
    for(my $i = 0; $i < $len; $i++)  {
	my $qchar = lc substr($qseq,$i,1);
	my $hchar = lc substr($hseq,$i,1);
	if( $qchar eq '-' || $hchar eq '-' ||
	    $qchar eq '.' || $hchar eq '.' ) {
	    $gaps++;
	} elsif( $qchar eq $hchar ) {
	    $match++;
	} else {
	    $mismatch++;
	}
    }
    return $match / ($len - $gaps);
}
