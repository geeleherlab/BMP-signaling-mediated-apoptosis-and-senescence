#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $Time_Start = &sub_format_datetime(localtime($BEGIN_TIME));
my $date=`date +"%d-%m-%y"`;
print "Program Starts Time:$Time_Start\n";
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fOut,$list);
GetOptions(
				"help|?" =>\&USAGE,
				"o:s"=>\$fOut,
				"i:s"=>\$fIn,
				"l:s"=>\$list,
				) or &USAGE;
&USAGE unless ($fIn and $fOut and $list);

sub USAGE {#
	my $usage=<<"USAGE";
ProgramName:
Version:	$version
Contact:	zhangyw <yz46606\@uga.edu> 
Program Date:   $date
Description:	this program is used to convert dmr_original list to standard gabit trait form
Usage:
  Options:
  -i <dir>  input dir,directory for peak files,forced
  -l <list>  input file, sample list sample,type,replicate
                                      DMSO_1,DMSO,1
                                      DMSO_2,DMSO,2
                                      RA_1,RA,1
                                      RA_2,RA,2
  -o <dir>  output file,forced  
  -h         Help

USAGE
	print $usage;
	exit;
}
#################################read bam files
$fIn=&ABSOLUTE_DIR($fIn);
my @bam=glob"$fIn/*.bam";
my $gsize="/research/groups/geelegrp/home/yzhang24/1_RA_BMP/2_RNAChip/chip_consensus_ref/1_RA_S9_S4/list/hg38.sizes";
#################################Read sample information
my @sample;
open (IN, $list) or die $!;
while (<IN>) {
	chomp;
	next if (/^$/||$.==1);
	my @lines=split/\,/,$_;
	push @sample,$lines[0];
}
close IN;

############recall reads on reference peaks
mkdir $fOut if (! -d $fOut);
open (OUT, ">paste.sh") or die $!;
my @count;
my @type;
my @repro;
#print Dumper @sample;die;
foreach my $s ( @sample){
    foreach my $b (@bam){
        my $name=basename($b);
        if ($name =~/$s/) {
            open (OUT2, ">$fOut/$s.count.sh") or die $!;
            print OUT2 "#BSUB -P CHIP\n";
            print OUT2 "#BSUB -J $s.count\n";
            print OUT2 "#BSUB -oo sander.out -eo sander.err\n";
            print OUT2 "#BSUB -n 1\n";
            print OUT2 "#BSUB -M 80000\n";
            print OUT2 "module load  bedtools/2.17.0\n";
            print OUT2 "samtools view -F 1024 -b -q 1 $b > $s.rmdupq1.bam\n";
            print OUT2 "bedtools bamtobed -i $s.rmdupq1.bam | awk \'\{if\(\$6 ~ /-/\)\{\$2=\$3-1;\}else\{\$3=\$2+1\}; OFS=\"\\t\"; print \$0\}\'|slopBed -s -l 0 -r 200 -g $gsize -i - > $s\.bed\n";
            print OUT2 "bedtools intersect -a all_reproducible.bed -b $s\.bed -c | awk \'\{print \$NF\}\' > $s.count\n";
            close OUT2;
           `bsub < $fOut/$s.count.sh`;
		    push @type,"$s";
            push @count,"$s.count";
        }
    }
}

#####paste all counts together
print OUT "awk \'\{print \$1\"\:\"\$2\"\-\"\$3\}\' all_reproducible.bed  > All.Repro.bed.pre\n";
my $l1=join(" ",@type);
my $l2=join(" ",@count);
#print Dumper @type;
print OUT "echo \"region $l1\" > all.counts.data\n";
print OUT "paste All.Repro.bed.pre $l2 >>all.counts.data\n";
close OUT;


#######################################################################################
my $Time_End   = sub_format_datetime(localtime(time()));
print STDOUT "Program Ends Time:$Time_End\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
#######################################################################################

sub ABSOLUTE_DIR{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

#######################################################################################

sub max{#&max(lists or arry);
	#���б��е����ֵ
	my $max=shift;
	my $temp;
	while (@_) {
		$temp=shift;
		$max=$max>$temp?$max:$temp;
	}
	return $max;
}

#######################################################################################

sub min{#&min(lists or arry);
	#���б��е���Сֵ
	my $min=shift;
	my $temp;
	while (@_) {
		$temp=shift;
		$min=$min<$temp?$min:$temp;
	}
	return $min;
}

#######################################################################################

sub revcom(){#&revcom($ref_seq);
	#��ȡ�ַ������еķ��򻥲����У����ַ�����ʽ���ء�ATTCCC->GGGAAT
	my $seq=shift;
	$seq=~tr/ATCGatcg/TAGCtagc/;
	$seq=reverse $seq;
	return uc $seq;			  
}

#######################################################################################

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

#######################################################################################

sub sub_format_datetime {#Time calculation subroutine
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


