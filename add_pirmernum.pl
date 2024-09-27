my (%hash,%amp);
open IN,"A336V1-20211119.primer" or die $!;
while(<IN>){
	chomp;
	unless($_){next;}
	my ($chr,$st1,$ed1,$st2,$ed2,$name) = (split /\t/,$_)[0,1,2,3,4,5];
	$amp{$name} = "$chr\t$st1\t$ed1\t$st2\t$ed2";
	$st1++;
	for my $pos($st1..$ed1){
		my $info = "$chr:$pos";
		if(exists $hash{$info}{'primer'}){$hash{$info}{'primer'} = $hash{$info}{'primer'}.",$name";}
		else{$hash{$info}{'primer'} = $name;}
	}
	for my $pos($st2..$ed2){
                my $info = "$chr:$pos";
		if(exists $hash{$info}{'primer'}){$hash{$info}{'primer'} = $hash{$info}{'primer'}.",$name";}
                else{$hash{$info}{'primer'} = $name;}
        }
	my $st = $ed1+1;
	my $ed = $st2-1;
	for my $pos($st..$ed){
                my $info = "$chr:$pos";
		if(exists $hash{$info}{'insert'}){$hash{$info}{'insert'} = $hash{$info}{'insert'}.",$name";}
                else{$hash{$info}{'insert'} = $name;}
        }
}
close IN;

open IN,"$ARGV[0]";
open OUT,">pos-stat.xls" or die $!;
while(<IN>){
	chomp;
	unless($_){next;}
	if(m/^\#/){next;}
	my @a = split(/\s+/,$_);
	my @b = split(/,/,$a[2]);
	for my $pos(@b){
		my $info = "chr$a[1]:$pos";
		if(exists $hash{$info}{'insert'}){print OUT "$a[0]\t$info\tinsert\t$hash{$info}{'insert'}\n";}
		elsif(exists $hash{$info}{'primer'}){print OUT "$a[0]\t$info\tprimer\t$hash{$info}{'primer'}\n";}
		else{print OUT "$a[0]\t$info\tnone\tnone\n";}
	}
}
close IN;
close OUT;

print "Amp-id\tChr\tPos\tLink-num\tPrimer-chr\tLeft-st\tLeft-ed\tRight-st\tRight-ed\n";
my ($id,$chr,$pos,$pname);
my $lnum=0;
open IN,"pos-stat.xls" or die $!;
while(<IN>){
	chomp;
	unless($_){next;}
	my @a = split(/\t/,$_);
	unless($a[2] eq "insert"){next;}
	unless($id){
		$id = $a[0];
		($chr,$pos) = (split /\:/,$a[1])[0,1];
		$pname = $a[3];
		$lnum=1;
		next;
	}
	if($a[0] eq $id && $a[3]=~/$pname/){
		my $newpos = (split /\:/,$a[1])[1];
		$pos = $pos.";$newpos";
		$lnum++;
	}else{
		print "$id-$pname\t$chr\t$pos\t$lnum\t$amp{$pname}\n";
		$id = $a[0];
                ($chr,$pos) = (split /\:/,$a[1])[0,1];
                $pname = $a[3];
		$lnum=1;
	}
}
close IN;
print "$id-$pname\t$chr\t$pos\t$lnum\t$amp{$pname}\n";
