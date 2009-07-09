#!/usr/bin/perl -w
# $Id: dotur_report2table.pl 1 2008-12-03 22:42:29Z stajich $
use strict;

# USAGE: perl dotur_report2table.pl REPORT > REPORT.table

my (%groups, @group_order);
while(<>) {
    chomp;
    my ($id,$group_count,@seqsets) = split(/\t/,$_);    
    push @group_order, $id;
    my $i = 1;
    for my $set ( @seqsets ) {
	for my $s ( split(/,/,$set) ) {
	    $groups{$s}->{$id} = $i;
	}
	$i++;
    }
}

print join("\t", qw(Sequence), @group_order),"\n";

# this sequences are sorted by the 1st group ('unique')
for my $sequence ( sort { $groups{$a}->{$group_order[0]} <=> 
			      $groups{$b}->{$group_order[0]} }
			      keys %groups)  {
    # some perl magic here - sorry if you don't understand 'map'    
    print join("\t", $sequence, 
	       map { $groups{$sequence}->{$_} } @group_order), "\n";
}
