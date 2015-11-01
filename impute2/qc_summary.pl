#!/usr/bin/perl -w
use strict;
#########################################
my ($infile, $outfile) = @ARGV;
print "Input $infile\n";
print "Output $outfile\n";
#########################################

open(OUT,">dose.$outfile") or die("Cannot open outfile for writing: $!\n");
open(IN,"<$infile") or die("Cannot open input file: $!\n");
while(my $line = readline(IN)) {
	next if($.==1);
	chomp($line);
	my @arr=split(/ /,$line);
	print OUT $arr[0]."\t".$arr[1]."\t".$arr[2]."\t".$arr[3]."\t".$arr[4]."\t";
	my $i=5;
	while(defined($arr[$i]))
		{
		my $dose = $arr[$i+1] + $arr[$i+2]*2;
		print OUT $dose."\t";
		$i=$i+3
		}
		print OUT "\n";
}
close(IN);
close(OUT);
system("date");

## execute R script
open(R,">$outfile.R") or die("Cannot open R script: $!\n");
print R "tmp <- read.table('dose.$outfile',header=FALSE,sep='\t')\n";
print R "tmp1 <- tmp[,c(1,2,3,4,5)]\n";
print R "names(tmp1) <- c('CHR','SNP','POS','A1','A2')\n";
print R "tmp <- tmp[,-c(1,2,3,4,5,ncol(tmp))]\n";
print R "tmp1[,c('Theta')] <- apply(tmp,1,function(x) mean(x)/2)\n";
print R "tmp1[,c('Ratio')] <- apply(tmp,1,function(x) (sum(x**2)/ncol(tmp)-mean(x)**2)/(2*mean(x)/2)*(1-mean(x)/2))\n";
print R "write.table(tmp1,'$outfile',sep='\t',col.names=TRUE,row.names=FALSE,quote=FALSE,na='')\n";
close(R);
system("R --vanilla <$outfile.R >$outfile.Rout");
