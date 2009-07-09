#!/usr/bin/perl -w
# $Id: otu_cluster.pl 1 2008-12-03 22:42:29Z stajich $
use strict;
use Bio::SearchIO;

my %matrix;
my $reportfile = shift @ARGV;

my $in = Bio::SearchIO->new(-format => 'fasta',
			    -file   => $reportfile);

while( my $r = $in->next_result ) {
    my $q = $r->query_name;
    $matrix{$q}{$q} = 0;
    while( my $hit = $r->next_hit ) {
	my $h = $hit->name;
	if( my $hsp = $hit->next_hsp ) {

	    # BEGIN copied from Bio::Search::HSP::FastaHSP::get_aln 
	    my ($hs,$qs) = ($hsp->hit_string,$hsp->query_string);

	    # we infer the end of the regional sequence where the first
	    # non space is in the homology string
	    # then we use the HSP->length to tell us how far to read
	    # to cut off the end of the sequence
	    
	    my ($start, $rest) = (0, 0);
	    if( $hsp->homology_string() =~ /^(\s+)?(.*?)\s*$/ ) {
		($start, $rest) = ($1 ? length($1) : 0, length($2));
	    }

	    $hs = substr($hs, $start,$rest);
	    $qs = substr($qs, $start,$rest);
	    # END copied from Bio::Search::HSP::FastaHSP::get_aln 

	    my $pid = &percent_id($hs,$qs);	    
	    $matrix{$q}{$h} = sprintf("%.8f",1-$pid);
	}
    }
}

my @names = sort keys %matrix;
print scalar @names, "\n"; # print out the 1st line of DOTUR matrix which is
                           # number of names
for my $name ( @names ) {
    my @l;
    for my $n ( @names ) {
	push @l, exists $matrix{$name}->{$n} ? $matrix{$name}->{$n} : 1;
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
    if( $len - $gaps ) { 
	return $match / ($len - $gaps);
    } else {
	warn("match is $match len is $len, gaps are $gaps\n");
	return 0;
    }
}
