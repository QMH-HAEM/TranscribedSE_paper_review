use strict;
use warnings;

my @SRR;
my @sample;

open IN, $ARGV[0], or die "Cannot open file. $!.";
while (<IN>){
   chomp;
   my @fields = split "\t";
   push @SRR, $fields[0];
   push @sample, $fields[1];
}
close IN;

for my $i (0..$#SRR){
   my $SRRnumber = $SRR[$i];
   my $sampleNo = $sample[$i];
   print "[", $i+1, "/", $#SRR+1, "]: Converting $SRR[$i].sra to fastq...", "\n";
   system("fasterq-dump $SRRnumber --threads 8 --mem 10000MB --split-3 --outdir ./ &");
   `mv ${SRRnumber}_1.fastq.gz ${sampleNo}_1.fastq.gz`;
   `mv ${SRRnumber}_2.fastq.gz ${sampleNo}_2.fastq.gz`;
   print "Finished converting $SRR[$i].sra to fastq.\n";
}
