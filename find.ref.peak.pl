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
#################################read peak files
my @peaks=glob"$fIn/*.filter.narrowPeak";

#################################Read sample information
my %sample;
open (IN, $list) or die $!;
while (<IN>) {
	chomp;
	next if (/^$/||$.==1);
	my @lines=split/\,/,$_;
	$sample{$lines[1]}{$lines[2]}=$lines[0];
}
close IN;

############reproducible peaks
open (OUT, ">$fOut") or die $!;
print OUT "module load  bedtools/2.17.0\n";
my @repro;
foreach my $sam (keys %sample){
	my @keepbed;
	foreach my $a (sort keys %{$sample{$sam}}){
		my $high;
		foreach my $p (@peaks){
			if($p=~/$sample{$sam}{$a}/&& $p!~/FDR50/){
				$high=$p;
			}
		}
		#print $high;die;
		foreach my $b (sort keys %{$sample{$sam}}){
			#print $ss;die;
			if($a ne $b){
				my $low;
				my $count=0;
				foreach my $p (@peaks){
					if($p=~/$sample{$sam}{$b}/&& $p=~/FDR50/){
						$low=$p;
						$count++;
					}
				}
				if($count==1){
					if($count==(keys %{$sample{$sam}})-1){
						print OUT "bedtools intersect -a $high -b $low -u > $sample{$sam}{$a}.keep.bed\n";
						push @keepbed,"$sample{$sam}{$a}.keep.bed";
						last;
					}
					else {
						print OUT "bedtools intersect -a $high -b $low -u > $sample{$sam}{$a}.tmp.bed\n";
					}
				}
				elsif($count>=(keys %{$sample{$sam}})-1){
					print OUT "bedtools intersect -a $sample{$sam}{$a}.tmp.bed -b $low -u > $sample{$sam}{$a}.keep.bed\n";
					push @keepbed,"$sample{$sam}{$a}.keep.bed";
				}
				else {
					print OUT "bedtools intersect -a $sample{$sam}{$a}.tmp.bed -b $low -u > $sample{$sam}{$a}.tmp.bed\n";
				}
			}
		}
	}
	my $list=join(" ", @keepbed);
	print OUT "cat $list|sortBed -i | mergeBed -i stdin > $sam\_reproducible.bed\n";
	push @repro,"$sam\_reproducible.bed";
	print OUT "rm $list\n";
}

#####merge peaks
my $lst=join(" ",@repro);
print OUT "cat $lst|sortBed -i | mergeBed -i stdin > all\_reproducible.bed\n";
####remove tmp files
print OUT "rm $lst\n";



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


