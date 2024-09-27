open IN,"$ARGV[0]";
open OUT,">add.bed" or die $!;
while(<IN>){
	chomp;
	unless($_){next;}
	if(m/^\#/){next;}
	my @a = split(/\t/,$_);
	my @b = (split /;/,$a[2]);
	for my $pos(@b){
		my $st = $pos-1;
		print OUT "$a[1]\t$st\t$pos\n";
	}
}
close IN;
close OUT;

system "~/igenePIP_V6/software/bedtools/bin/fastaFromBed -fi ~/igenePIP_V6/database/hg19/Genome/ucsc.hg19.fasta -bed add.bed -fo add.re";

my %hash;
open IN,"add.re" or die $!;
while(<IN>){
	chomp;
	unless($_){next;}
	if(m/^\>(\S+)/){
		my ($chr,$pos) = (split /[: -]/,$1)[0,2];
		my $base = <IN>;
		chomp $base;
		$base = uc $base;
		$hash{"$chr:$pos"} = $base;
	}
}
close IN;
open IN,"$ARGV[0]";
while(<IN>){
        chomp;
        unless($_){next;}
        if(m/^\#/){print "$_\n";next;}
	my @a = split(/\t/,$_);
        my @b = (split /;/,$a[2]);
	for(my $i=0;$i<@b;$i++){
		my $info = "$a[1]:$b[$i]";
		$b[$i] = "$b[$i]-$hash{$info}";
	}
	$a[2] = join(";",@b);
	my $re = join("\t",@a);
	print "$re\n";
}
close IN;
