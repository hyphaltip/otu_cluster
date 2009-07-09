#!/usr/bin/perl -w
=head1 NAME

dotur2file

=head1 USAGE

dotur2file.pl --db seqs.fa result.fn.list

=head1 DESCRIPTION

Pass a sequence database (in FastA format) and the DOTUR .list file to
build sequence files for each cluster which can be aligned.

=cut
use strict;
use Getopt::Long;
use Bio::DB::Fasta;
use Bio::SeqIO;

my $prefix = 'dotur_clust';
my $db;
GetOptions(
	   'd|db:s'     => \$db,
	   'p|prefix:s' => \$prefix,
	   );

my $in = shift @ARGV;
$db ||= shift @ARGV;

my $dbh= Bio::DB::Fasta->new($db);

open(my $fh => $in) || die "$in: $!";
while(<$fh>) {
    chomp;
    next if /^\s+$/;
    my ($cutoff,$count,@clusters) = split(/\t/,$_);
    my $n = 1;	
    my $dir = sprintf("%s_%s",$prefix,$cutoff);
    mkdir($dir);
    for my $c ( @clusters ) {
	my $out = Bio::SeqIO->new(-format => 'fasta',
				  -file   => sprintf(">%s/c%d.fa",
						     $dir,$n++));
	for my $s ( split(/,/,$c) ) {
	    my $seq = $dbh->get_Seq_by_acc($s);
	    if( ! defined $seq ) {
		warn("cannot find seq $s in db $db\n");
	    } else {
		$out->write_seq($seq);
	    }	    
	}
    }
}
