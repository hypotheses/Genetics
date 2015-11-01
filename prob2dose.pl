#!/usr/bin/perl -w
use strict;
#########################################
my ($infile, $outfile) = @ARGV;
print "Input $infile\n";
print "Output $outfile\n";
#########################################
#open(OUT,">>dose.$outfile") or die("Cannot open outfile for writing: $!\n");
open(IN,"<$infile") or die("Cannot open input file: $!\n");
while(my $line = readline(IN)) {
	next if($.==1);
	chomp($line);
	my @arr=split(/ /,$line);
	#print $arr[0].",".$arr[1].",".$arr[2].",".$arr[3].",".$arr[4];
	my $i=5;
	while(defined($arr[$i]))
		{
		my $dose = $arr[$i+1] + $arr[$i+2]*2;
		print $dose;
		$i=$i+3;
		if ( defined($arr[$i]) ) {
		    print ",";
		}
		}
		print "\n";
}
close(IN);
## close(OUT);
system("date");

