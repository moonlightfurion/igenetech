open IN,"$ARGV[0]";
while(<IN>){
	chomp;
	unless($_){next;}
	my $n = (split /\s+/,$_)[0];
	system "perl microhaplotype.pl -amp amp-config.txt -bam $n/align/$n.sorted.bam -s $n";
}
close IN;
