#!/usr/bin/perl
#date: 20230403
#author: qian
use File::Basename;
use Getopt::Long;
use Switch;
use FindBin qw($Bin);

BEGIN{
    push(@INC,$Bin);
};

my $bin = dirname $Bin;
my $usage="usage:
    $0 [Options]
Options:
    -amp  [s] amplicon file contain necessary info <requied>
    -bam  [s] bam file by bwa mem <requied>
    -s    [s] sample name <default: sample1>
    -o    [s] output dir <default: ./>
    -mq   [s] min map quality <default: 20>
    -bq   [s] min base quality <default: 20>
    -dp   [s] min depth for genotype <default: 20>
    -fre  [s] min VAF for genotype <default: 0.25>
author:
    shaoqianzhi
";

unless(@ARGV>0){
    print "$usage";
    exit;
}

GetOptions(
    "amp=s" => \$ampliconfile,
    "bam=s" => \$bamfile,
    "s=s" => \$sample,
    "o=s" => \$outdir,
    "mq=s" => \$minmapquality,
    "bq=s" => \$minbasequality,
    "dp=s" => \$mindp,
    "fre=s" => \$minfre,
);
# default
$sample||="sample1";
$outdir||="./";
unless(-e $outdir){system "mkdir -p $outdir";}
$minmapquality||=20;
$minbasequality||=20;
$mindp||=20;
$minfre||=0.25;
# start
my $time = localtime();
print "Info: $time start ana $sample\n";
# read amp
my (%hashamp,%amp);
my @allamp;
my $extendsize = 5; # 此参数为primer扩展bp数，增加兼容性
open IN,"$ampliconfile" or die $!;
# 1-T1P1  chr1    1529950-A;1529979-T;1529994-G;1529998-C 4       chr1    1529913 1529940 1530078 1530101
while(<IN>){
    chomp;
    unless($_){next;}
    if(m/^\#/){next;}
    my @a = split(/\t/,$_);
    my $chr = $a[4];
    my $st1 = $a[5]+1-$extendsize;
    my $ed1 = $a[6];
    my $st2 = $a[7];
    my $ed2 = $a[8]+$extendsize;
    for my $pos($st1..$ed1){
        my $info = "$chr:$pos";
        $hashamp{'left'}{$info}=$a[0];
    }
    for my $pos($st2..$ed2){
        my $info = "$chr:$pos";
        $hashamp{'right'}{$info}=$a[0];
    }
    $amp{$a[0]} = $a[2];
    push(@allamp,$a[0]);
    # my @allpos = split(/;/,$a[2]);
    # for my $value(@allpos){
    #    my ($pos,$ref) = (split /\-/,$value)[0,1];
    #    my $info = "$a[1]:$pos";
        #$hashamp{'pos'}{$info}=$a[0];
    #    $amp{$a[0]}{$info} = $ref;
    # }
}
close IN;
# read input bam
my ($total,$pass1,$pass2) = (0,0,0);
open IN,"samtools view -F 256 $bamfile |" or die $!;
while(<IN>){
    chomp;
    unless($_){next;}
    $total++;
    my @a = split(/\t/,$_);
    unless($a[6] eq "="){next;}
    if($a[5] eq "*" || $a[5]=~/H/){next;}
    unless($a[4]>=$minmapquality){next;}
    #check pos
    my $refplus = 0;
    while($a[5]=~/(\d+)[MD]/mg){$refplus = $refplus + $1;}
    my $refst = $a[3];
    my $refed = $a[3] + $refplus - 1;
    my $info1 = "$a[2]:$refst";
    my $info2 = "$a[2]:$refed";
    my $ampname;
    if($a[8]>0){
        unless(exists $hashamp{'left'}{$info1}){next;}
        $ampname = $hashamp{'left'}{$info1};
        #print "$hashamp{'left'}{$info}-left\t$_\n";
    }elsif($a[8]<0){
        unless(exists $hashamp{'right'}{$info2}){next;}
        $ampname = $hashamp{'right'}{$info2};
        #print "$hashamp{'right'}{$info}-right\t$_\n";
    }else{
        next;
        #die "Error: $_ insert error \n";
    }
    $pass1++;
    # check NM MD
    my ($NM,$MD);
    if($a[11]=~/NM:i:(\d+)/){$NM=$1;}
    else{die "Error: input bam format error \n";}
    if($a[12]=~/MD:Z:([^ ]+)/){$MD=$1;}
    else{die "Error: input bam format error \n";}
    # rm soft clip
    my $id = $a[0];
    my $read = $a[9];
    my $quality = $a[10];
    my $flag = $a[5];
    if($flag=~/^(\d+)S/){
        $read = substr($read,$1);
        $quality = substr($quality,$1);
    }
    if($flag=~/(\d+)S$/){
        my $p = (length $read)-$1;
        $read = substr($read,0,$p);
        $quality = substr($quality,0,$p);
    }
    # exact snp , insert and del
    my %hash;
    my ($readpos_flag,$seqpos_flag,$readpos_md,$seqpos_md,$nm) = (0,0,0,0);
    # change with flag
    while($flag=~/(\d+)([MID])/mg){
        switch($2){
            case "M"{
                $readpos_flag = $readpos_flag + $1;
                $seqpos_flag = $seqpos_flag + $1;
            }
            case "I"{
                my $newpos = $a[3]+$readpos_flag-1;
                my $newinsert = substr($read,$seqpos_flag-1,$1+1);
                my $newquality = substr($quality,$seqpos_flag-1,$1+1);
                #my $averquality = calqu($newquality);
                $hash{$newpos}{'base'} = lc($newinsert); # 插入的碱基用小写表示
                $hash{$newpos}{'quality'} = "$newquality";
                # $read = substr($read,0,$readpos_flag).substr($read,$readpos_flag+$1);
                # $quality = substr($quality,0,$readpos_flag).substr($quality,$readpos_flag+$1);
                $nm = $nm + $1;
                $seqpos_flag = $seqpos_flag + $1;
                # print "$a[2]:$newpos:$newinsert:$newquality\n";
            }
            case "D"{
                my $newpos_st = $a[3]+$readpos_flag;
                my $newpos_ed = $a[3]+$readpos_flag+$1-1;
                my $newquality = substr($quality,$seqpos_flag-2,2);
                #my $averquality = calqu($newquality);
                for my $newpos($newpos_st..$newpos_ed){
                    $hash{$newpos}{'base'} = "-"; # 待修改， 缺失的碱基用-表示
                    $hash{$newpos}{'quality'} = "$newquality";
                    # print "$a[2]:$newpos:-:$newquality\n";
                }
                #my $seq = "-" x $1;
                $nm = $nm + $1;
                $readpos_flag = $readpos_flag + $1;
            }
        }
    }
    # print "$nm\t$readpos_flag\t$read\n";
    # change with MD
    while($MD ne ""){
        if($MD=~/^(\d+)/){
            $readpos_md = $readpos_md + $1;
            $seqpos_md = $seqpos_md + $1;
            # print "S rmd is $readpos_md $seqpos_md\n";
            $MD=~s/^\d+//;
        }elsif($MD=~/^([ATCG])/){   #突变
            my $newpos = $a[3]+$readpos_md;
            my $newbase = substr($read,$seqpos_md,1);
            my $newquality = substr($quality,$seqpos_md,1);
            $hash{$newpos}{'base'} = "$newbase";
            $hash{$newpos}{'quality'} = "$newquality";
            $readpos_md = $readpos_md + 1;
            $seqpos_md = $seqpos_md + 1;
            # print "M rmd is $readpos_md $seqpos_md\n";
            $MD=~s/^([ATCG])//;
            $nm++;
            # print "$newpos:$newbase:-:$newquality\n";
        }elsif($MD=~/^\^([ATCG]+)/){   #缺失
            $MD=~s/^\^([ATCG]+)//;
            my $len = length $1;
            $readpos_md = $readpos_md + $len;
            # print "D rmd is $readpos_md $seqpos_md\n";
        }
    }
    # check nm
    unless($nm==$NM){die "Error: nm not same by $_\n";}
    # foreach my $key(keys %hash){
    #    print "$key\t$hash{$key}{'base'}\t$hash{$key}{'quality'}\n";
    # }
    # save amp info
    my $genotype;
    my @allpoint = split(/;/,$amp{$ampname});
    for my $value(@allpoint){
        #print "ana $value with $refst and $refed \n";
        my ($pos,$ref) = (split /\-/,$value)[0,1];
        my $rebase = "*";
        if($pos>=$refst && $pos<=$refed){
            if(exists $hash{$pos}){
                my $q = calqu($hash{$pos}{'quality'});
                if($q>=$minbasequality){$rebase = $hash{$pos}{'base'};}
            }else{
                $rebase = $ref;
            }
        }
        if($genotype){$genotype = $genotype.":$rebase";}
        else{$genotype = $rebase;}
        #print "out $rebase and $genotype \n";
    }
    if($genotype){
        if(exists $amp{$ampname}{$id}){$amp{$ampname}{$id}=mergeid($amp{$ampname}{$id},$genotype);}
        else{$amp{$ampname}{$id}=$genotype;}
        #print "read $ampname with $amp{$ampname}{$id}\n";
    }
}
close IN;
my $p = sprintf("%.2f",$pass1*100/$total);
print "Info: Total $total reads, $pass1 are useful, percent $p% \n";
# out detail result
open OUT,">$outdir/$sample.detail.xls" or die $!;
for my $singleamp(@allamp){
    my %hash;
    my ($totaldp,$totaluseful) = (0,0);
    foreach my $id(keys %{$amp{$singleamp}}){
        my $genotype = $amp{$singleamp}{$id};
        $totaldp++;
        if($genotype=~/\*/){next;}
        $totaluseful++;
        $hash{$genotype}++;
    }
    print OUT "$singleamp\ttotaldp=$totaldp\tusefuldp=$totaluseful\n";
    foreach my $key(sort{$hash{$b}<=>$hash{$a}} keys %hash){
        my $g = $key;
        #$g=~s/://g;
        my $p = sprintf("%.2f",$hash{$key}*100/$totaluseful);
        unless($p>=1){next;}
        print OUT "\t$g\t$hash{$key};$p%\n";
    }
}
close OUT;
# out genotype result
open IN,"$outdir/$sample.detail.xls" or die $!;
open OUT,">$outdir/$sample.genotype.xls" or die $!;
my $outname = "NA";
my $outtype = ".";
my $outtypenum = 0;
print OUT "#Amp-id\t$sample\n";
while(<IN>){
    chomp;
    unless($_){next;}
    my @a = split(/\t/,$_);
    if($a[0]){
        if($outname eq "NA"){$outname=$a[0];}
        else{
            print OUT "$outname\t$outtype\n";
            $outname=$a[0];
            $outtype = ".";
            $outtypenum = 0;
        }
    }else{
        if($outtypenum>=2){next;}
        my $type = $a[1];
        my ($dp,$fre) = (split /[; \%]/,$a[2])[0,1];
        unless($dp>=$mindp){next;}
        unless($fre>=$minfre*100){next;}
        $type=~s/://g;
        if($outtype eq "."){$outtype = $type;}
        else{$outtype = $outtype."/$type";}  # 输出的分隔符，- / ?
        $outtypenum++;
    }
}
print OUT "$outname\t$outtype\n";
close IN;
close OUT;
# end
$time = localtime();
print "Info: $time end ana $sample\n";
# zi han shu
sub calqu
{
    my $quality = shift;
    my @q = split(//,$quality);
    my $avq = 0;
    for my $sq(@q){
        $avq = $avq + ord("$sq") - 33;
    }
    my $l = @q;
    $avq = sprintf("%.2f",$avq/$l);
    return $avq;
}

sub mergeid
{
    my ($type1,$type2) = @_;
    my @a = split(/:/,$type1);
    my @b = split(/:/,$type2);
    my @re;
    for(my $i=0;$i<@a;$i++){
        if($a[$i] eq "*" && $b[$i] ne "*"){push (@re,$b[$i]);}
        elsif($a[$i] ne "*" && $b[$i] eq "*"){push (@re,$a[$i]);}
        elsif($a[$i] ne "*" && $b[$i] ne "*"){
            if($a[$i] eq $b[$i]){push (@re,$a[$i]);}
            else{push (@re,"*");}
        }else{
            push (@re,"*");
        }
    }
    my $result = join(":",@re);
    return $result;
}