use warnings;
use strict;
use Getopt::Std;
use File::Basename;

##################################################################################
# a script to clean up reads (adaptor/duplicate/contamination removal), trimming #
# external dependencies: flash, trimmomatic, bowtie2, cutadapt                   # 
# written by Sonal Singhal, sonal.singhal1 [at] gmail.com, 1 Dec 2012            #      
# updated by Ke Bi, kebi [at] berkeley.edu, 11 Oct 2013                          #	
##################################################################################




die (qq/

(version 1.11 July 31 2013)

Usage: 2-scrubReads.pl [options]

Options:
-f  CHAR   Where all raw sequence fastq files are kept
-o  CHAR   Where results will go
-a  FILE   A file containing all adapter sequences
-b  FILE   A file containing library information
-t  CHAR   Path to trimmoatic executable
-c  CHAR   Contaminant file
-e  INT    Average fragment length for your library
-m  INT    Number of bases (nt) in the begining of reads to infer duplication [20]
-w         If w is supplied, assume the dual barcode indexing. Otherwise the standard one barcode system [null]
-z         If z is supplied, use fastQC to evaluate cleaned sequence reads [null]
-i  CHAR   Instrument ID (i.e. HS, MS or HWI) [HS]
-d  CHAR   The particular library that you like to process, if "all" is specified, all libraries in the folder (-f) will be processed [all]
-g  FLOAT  Maximum allowed ratio between the number of mismatched base pairs and the overlap length [0.05]
-l  INT    Read length [100]
-n  FLOAT  Will get rid of reads for which more than n*100\% of bases are NNs [0.6]
-r  FLOAT  Will get rid of reads with any runs of bases longer than r*read_length [0.5]
-h  INT    Trimmomatic trimming cutoff [36]

Note:
Example files for -a, -b and -c are kept in https:\/\/github.com\/MVZSEQ\/XXXXXX.
It assumes a library naming convention of ONLY letters & numbers with a ABC_R[1|2].fa ending.
It assumes HiSeq reads; will work with MiSeq with a small modification.
It assumes cutadapt, cope, bowtie2, fastQC, flash in path.

\n/) if !@ARGV;

my %opts = (f=>undef, o=>undef, a=>undef, b=> undef, i=>'HS', t=>undef, c=>undef,  l=>100, n=>0.6, r=>0.5, e=>undef, g=>0.05, d=>"all", h=>36, m=>20);  
getopts('f:o:a:t:c:l:n:r:b:e:g:m:d:i:h:zw', \%opts);

#print out error messgaes
die(qq/\nHmm...I do not know where you are hiding your raw sequence files... Check "-f"? \n\n/) if (!$opts{f});
die(qq/\nHmm...I do not know where you want to save the results... Check "-o"? \n\n/) if (!$opts{o});
die(qq/\nHmm...I could not find where the adapter sequence file is.. Check "-a"? \n\n/) if (!$opts{a});
die(qq/\nHmm...I could not find where the library info file is.. Check "-b"? \n\n/) if (!$opts{b});
die(qq/\nHmm...I could not find where the contaminant file is.. Check "-c"? \n\n/) if (!$opts{c});
die(qq/\nHmm...I do not know what the instrument ID is.. Check "-i"? \n\n/) if (!$opts{i});
die(qq/\nHmm...What is the number of leading bp in a read to inder duplicates... Check "-m"? \n\n/) if (!$opts{m});
die(qq/\nHmm...I could not locate where the trimmomatic executable is.. Check "-t"? \n\n/) if (!$opts{t});
die(qq/\nMissing  version 1.10average fragment length for the library. Check "-e"? \n\n/) if (!$opts{e});
die(qq/\nYou need to specify a library name after -d, or do not use -d to enable the defualt setting [all]. \n\n/) if (!$opts{d});
die(qq/\nYou need to specify a value after -g, or do not use -g to enable the defualt setting [0.05]. \n\n/) if (!$opts{g});
die(qq/\nYou need to specify a value after -l, or do not use -l to enable the defualt setting [100]. \n\n/) if (!$opts{l});
die(qq/\nYou need to specify a value after -n, or do not use -n to enable the defualt setting [0.6]. \n\n/) if (!$opts{n});
die(qq/\nYou need to specify a value after -r, or do not use -r to enable the defualt setting [0.5]. \n\n/) if (!$opts{r});
die(qq/\nYou need to specify a value after -h, or do not use -h to enable the defualt setting [36]. \n\n/) if (!$opts{h});


my $InID = $opts{i}; 

my $dir;
if ($opts{f} =~m/\/$/) {
$dir = $opts{f};
}
else {
$dir = $opts{f} . "/";
}

my $outdir;
if ($opts{o} =~m/\/$/) {
$outdir = $opts{o};
}
else {
$outdir = $opts{o} . "/";
}


mkdir ($outdir) unless -d ($outdir);

my $tmpdir = $outdir . "tmp/";
mkdir $tmpdir unless -d $tmpdir;

my $adaptorFile = $opts{a};

die(qq/\nHmm...I do not know where the adapter sequence file is.. Check "-a"? \n\n/) unless -e $adaptorFile;

my $libInfo = $opts{b};

die(qq/\nHmm...I could not find where the library info file is.. Check "-b"? \n\n/) unless -e $libInfo;

my $contam = $opts{c};
die(qq/\nHmm...I could not find where the contaminant file is.. Check "-c"? \n\n/) unless -e $contam;

my $uniad = 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT';


my $trimmomatic = $opts{t};
die(qq/\nHmm...I could not locate where the trimmomatic executable is.. Check "-t"? \n\n/) unless -e $trimmomatic;

my $cutadapt = 'cutadapt';
my $flash = 'flash';
my $bowtie = 'bowtie2';
my $cope = 'cope';
my $readLength = $opts{l};
my $nper = $opts{n};
my $aper = $opts{r};
my $minLength = $opts{h};
my $leading = $opts{m};



my @files;
if ($opts{d} eq "all") {
@files = <$dir*_R1.fq>; 
die(qq/\nHmm...I do not know where you are hiding your raw sequence files... Check "-d" and\/or "-f"? \n\n/) if (scalar(@files) == 0);
print "\n","OK! Now processing all data files!", "\n";
}

else {
@files = <$dir*$opts{d}_R1.fq>; 
die(qq/\nHmm...I do not know where you are hiding your raw sequence files... Check "-d" and\/or "-f"? \n\n/) if (scalar(@files) == 0);
print "\n","OK! Now processing library $opts{d}!", "\n";

}


foreach my $file1 (@files) {
  open (IN,"<",$file1);
  my $firstline = <IN>;
  die(qq/\nHmm.. The instrument ID is NOT $InID... Check "-i"? \n\n/) if ($firstline !~ m/^\@$InID\S+/ );
  close IN;
  my $file2 = $file1;
  $file2 =~ s/_R1.fq/_R2.fq/;
  my $lib = $1 if basename($file1) =~ m/(\S+)_R[1|2]/i;   
  
  my $start1 = time;	
  my $dup = $outdir . $lib . '.duplicates.out';
  remove_dup($file1, $file2, $leading, $dup);
  my $time1 = int((time - $start1)/60);
  print "Found duplicates in $lib in $time1 minutes! Now this is something...\n";
  my $start2 = time;
  my $low = $outdir . $lib . '.lowComplexity.out';
  removeLowComplexity($file1,$low); removeLowComplexity($file2,$low);
  my $time2 = int((time - $start2)/60);
  print "Found low complexity reads in $lib in $time2 minutes! Whew, almost there!\n";
  
  my $start3 = time;
  my %reads = ('1' => $file1, '2' => $file2);
  my $ad = getAdaptors($lib,$adaptorFile,$libInfo);
  my @clean1;
  for my $dir (keys %reads) {
    my $out1 = trimmomatic($lib, $ad,$reads{$dir},$reads{$dir},'trim1');
    my $out2 = cutadapt($ad,$out1,$reads{$dir},'trim2');
    my $out3 = bowtie($lib,$ad,$out2,$reads{$dir},'trim3');
    my $out4 = cutadapt($ad,$out3,$reads{$dir},'trim4');
    my $final = $reads{$dir} . '_cleaned1';
    my $call = system("mv $out4 $final");
    unlink($out1,$out2,$out3);
    push(@clean1,$final);
  }	
  my $trim1 = fixMatePair($lib,\%reads,\@clean1,"trim1");	
  my $reads2 = mergeReads($lib,\%reads,$trim1,"trim2");

  my %reads2 = %{$reads2};
  my @clean2;
  for my $dir (keys %reads2) {
    my $out1 = trimmomatic($lib, $ad,$reads2{$dir},$reads2{$dir},'trim1');
    my $out2 = cutadapt($ad,$out1,$reads2{$dir},'trim2');
    my $out3 = bowtie($lib, $ad,$out2,$reads2{$dir},'trim3');
    my $out4 = cutadapt($ad,$out3,$reads2{$dir},'trim4');
    my $final = $reads2{$dir} . '_cleaned2';
    my $call = system("mv $out4 $final");
    unlink($out1,$out2,$out3);
    push(@clean2,$final);
  }
  my $trim3 = fixMatePair($lib,\%reads,\@clean2, "trim2");
  my $reads3 = mergeReads($lib,\%reads,$trim3,"trim3");
  
  my $reads4 = reallyMergeReads($lib,\%reads,$reads3,"trim4","trim5",$InID);
  my $time3 = int((time - $start3)/60);
  print "Trimmed and merged in $lib in $time3 minutes! Whew, almost there!\n";
  
  my $start4 = time;
  my $contaminants = $outdir . $lib . '.contam.out';
  removeContamination($reads4,$contam,$contaminants);
  my $time4 = int((time - $start4)/60);
  print "Removed contamination in $lib in $time4 minutes! It's going...\n";	
  
  makeFinal($outdir,$lib,$reads4,$dup,$low,$contaminants);
  
  my $call_rm1 = system("rm $dir*$lib*clean*");
  my $call_rm2 = system("rm $dir*$lib*trim*");
}

if ($opts{z}) {
  print "\n","start evaluating ... OH WOW! Your sequence data are looking much better now!" , "\n\n";
  my @clean;
  
  if ($opts{d} eq "all") {
    @clean = < $outdir*_final.txt> ;
  }
  
  else {
    @clean = <$outdir$opts{d}*_final.txt>; 
  }
  
  my $resdir =  $outdir.'evaluation/';
  mkdir $resdir unless -e $resdir;
  foreach (<@clean>) {
    my $lib = $1 if basename($_) =~ m/(\S+)_[1|2|u]_final.txt/; 
    my $call1 = system("fastqc -t 2 $_ -o $resdir");
    die(qq/\nThe program "fastQC" is not in path! \n\n/) if ($call1 == -1 );
    
    system ("rm $resdir$lib*fastqc.zip");
  }
  
}	

print "\n","Congratulations! The cleanup step is done... Now what?" , "\n\n"; 

sub makeFinal {
  my ($outdir, $lib,$reads,$dup,$low,$contam) = @_;
  
  my %junk;
  
  open(IN, "<$dup");
  while(<IN>) {
    chomp(my $line = $_);
    $junk{$line}++;
  }
  close(IN);
  open(IN, "<$low");
  while(<IN>){
    chomp(my $line = $_);
    $junk{$line}++;
  }
  close(IN);
  open(IN, "<$contam");
  while(<IN>) {
    chomp(my $line = $_);
    $junk{$line}++;
  }
  close(IN);
  
  my $new1 = $outdir . $lib . "_1_final.txt";
  my $new2 = $outdir . $lib . "_2_final.txt";
  my $newu = $outdir . $lib . "_u_final.txt";
  my %new = ('1' => $new1, '2' => $new2, 'u' => $newu);
  my %reads = %{$reads};
  foreach my $type (keys %reads) {
    open(OUT, ">$new{$type}");
    open(IN, "<$reads{$type}");
    while(<IN>) {
      chomp(my $line = $_);
      my @line =split (/\t+/, $line);
      if (scalar(@line) == 1) {
	if ($line[0] =~ m/^@($InID\S+)\/[1|2]$/) {
	  my $id = $1;
	  my $seq = <IN>; my $qualid = <IN>; my $qual = <IN>;
	  unless($junk{$id}){
	    print OUT $line, "\n", $seq,$qualid,$qual;
	  }
	}
      }	    
    }
    close(IN); close(OUT);
  }
}

sub getAdaptors {
  my ($lib,$adaptorFile,$libInfo) = @_;
  
  my %lib;
  open(IN, "<$libInfo");
    
    my $header = <IN>;

    while(<IN>) {
    chomp(my $line = $_);
    my @d = split(/\t/,$line);
    if ($d[2]) {
      $lib{$d[0]} = {'P7'=>$d[1], 'P5'=>$d[2]};
    }
    if (!$d[2]) {
      $lib{$d[0]} = {'P7'=>$d[1]};
    }
  }
  close(IN);
  
  my %P7;
  my %P5;
  
  
  open(IN, "<$adaptorFile");
  while(<IN>) {
    chomp(my $line = $_);
    if ($line =~ m/>P7_index(\d+)/) {
      my $bc = $1;
      chomp(my $seq = <IN>);
      $P7{$bc} = $seq;
    }
    if ($line =~ m/>P5_index(\d+)/) {
      my $bc2 = $1;
      chomp(my $seq2 = <IN>);
      $P5{$bc2} = $seq2;
    }
    
  } 
  close(IN);	  
  my %ad;
  if (! $opts{w}) {
    %ad = ("uni" => $uniad, "uni_rc" => rc($uniad), "index" => $P7{$lib{$lib}{'P7'}}, "index_rc" => rc($P7{$lib{$lib}{'P7'}}));
    die(qq/\nHmm...I could not find the adapter sequences for $lib ... Check the naming of the adapter sequences in "-a". \n\n/) if (!$ad{"index"});
  }
  if ($opts{w}) {
    %ad = ("uni" =>$P5{$lib{$lib}{'P5'}}, "uni_rc" => rc($P5{$lib{$lib}{'P5'}}), "index" => $P7{$lib{$lib}{'P7'}}, "index_rc" => rc($P7{$lib{$lib}{'P7'}}));
    die(qq/\nHmm...I could not find the adapter sequences for $lib ... Check the naming of the adapter sequences in "-a". \n\n/) if (!$ad{"index"} || !$ad{"uni"});
    
  } 
    
  return(\%ad);
}

sub rc {
  my ($seq) = @_;
  my $rc = $seq;
  $rc = reverse($rc);
  $rc =~ tr/ATGCatgc/TACGtacg/;
  return($rc);
}

sub removeContamination {
  my ($reads,$contam,$contaminants) = @_;
  unless (-f $contam . ".3.bt2") {
    my $bw2build = $bowtie . "-build";
    my $call1 = system("$bw2build $contam $contam");
  }                
  my %reads = %{$reads};
  my $contamout1 = $reads{'1'} . ".contam.sam1";
  my $call2 = system("$bowtie -x $contam -1 $reads{'1'} -2 $reads{'2'} --fast -S $contamout1 --sam-nohead --sam-nosq");
  die(qq/\nThe program "bowtie2" is not in path! \n\n/) if ($call2 == -1 );
  my $contamout2 = $reads{'1'} . ".contam.sam2";
  my $call3 = system("$bowtie -x $contam -U $reads{'u'} --fast -S $contamout2 --sam-nohead --sam-nosq");
  my $contamout_all = $reads{'1'} . ".contam_all.sam";
  system ("cat $contamout1 $contamout2 > $contamout_all");
  parseSAM($contaminants,$contamout_all);
  system ("rm $contamout1 $contamout2");  
}



sub parseSAM {
  my ($contaminants,$contam) = @_;
  open(OUT, ">$contaminants");
  open(IN, "<$contam");
  
  while(<IN>) {
    chomp(my $line = $_);
    my @d = split(/\t/,$line);
    if ($d[2] !~ m/\*/) {
      if ($d[5] =~ m/\d+M/) {
	my $md = $1 if $line =~ m/(MD\:Z\S+)/;
	my @a = ($md =~ m/([ATGC])/g);
	if (scalar(@a) < 2) {
	  $d[0] =~ s/\/[1|2]$//;
	  print OUT $d[0], "\n";
	}
      }
    }
  }
  close(IN); close(OUT);  
  unlink($contam);
}       

sub remove_dup {
  my ($file1, $file2, $leading, $dup) = @_;
  open (my $fh1, $file1);
  open (my $fh2, $file2);
  
  open (OUT, ">", $dup);
  my $line1=readline($fh1);
  my $line2=readline($fh2);
  
  my $count=1;
  my %hash;
  
  my $id;
  while ($line1) {
    if  ($count % 4 == 1) {
      $id = $1 if $line1 =~ m/^\@($InID\S+)\/[1|2]$/; 
    } 
    if ( $count % 4 == 2 ) {
      my $part1;
      my $part2;
      
      chomp($line1);
      chomp($line2);    
      
      $part1 = substr($line1, 0,  $leading);
      $part2 = substr($line2, 0,  $leading);
      
      my $key= $part1. "-" . $part2;
      
      if ($hash{$key}) {
	print OUT $id, "\n";
      }
      else {
	$hash{$key}++;
      }  
    }  
    $line1=readline($fh1);
    $line2=readline($fh2);
    $count++;
  }  
}



sub removeLowComplexity {
  my ($file,$low) = @_;
  open(IN, "<$file");
  open(OUT, ">>$low");
  while(<IN>) {
    chomp(my $line = $_);		
    if ($line =~ m/^@(\S+)\/[1|2]$/) {
      my $id = $1;
      chomp(my $seq = <IN>);
      my $n = int($nper*length($seq));
      my $a = int($aper*length($seq));
      my $ncounter = ($seq =~ m/N/g);
      if ($seq =~ m/[A]{$a}/i || $seq =~ m/[T]{$a}/i || $seq =~ m/[G]{$a}/i || $seq =~ m/[C]{$a}/i || $ncounter >= $n) {
	print OUT $id, "\n";
      }
    }
  }	
  close(IN); close(OUT);	
}

sub mergeReads {
  my ($lib,$orig,$reads,$base) = @_;
  my %reads = %{$reads};
  
  my $newread1 = $orig->{'1'} . '_' . $base . '_p1';
  my $newread2 = $orig->{'2'} . '_' .$base .'_p2';
  my $newreadu = $orig->{'1'} . '_' .$base .'_u';
  
  #my $call1 = system("$flash $reads{'1'} $reads{'2'} -M 100 -m 5 -x $opts{g} -f $opts{e} -o $lib");
  my $call1 = system("$flash $reads{'1'} $reads{'2'} -M 159 -m 10 -x $opts{g} -f $opts{e} -o $lib"); #use if MiSeq 150PE is run
  die(qq/\nThe program "flash" is not in path! \n\n/) if ($call1 == -1 );
  
  open (EXTEND,"<",  $lib . ".extendedFrags.fastq");
  open (NEW, ">",  $lib . ".extendedFrags.fastq1");
  while (<EXTEND>) {
    chomp(my $line = $_);	
    if ($line =~ m/^@($InID\S+)/) {
      my $new_line =  $line . "/1";
      my $seq =<EXTEND>;
      my $qualid = <EXTEND>;
      my $qual = <EXTEND>;
      print NEW $new_line , "\n" , $seq , $qualid , $qual;
     } 
  }
  close EXTEND;
  close NEW;
  my $call2 = system("cat $reads{'u'} $lib\.extendedFrags.fastq1 > $newreadu");
  my $call3 = system("mv $lib\.notCombined_1.fastq $newread1");
  my $call4 = system("mv $lib\.notCombined_2.fastq $newread2");
  my $call5 = system("rm $lib\.extendedFrags.fastq $lib\.extendedFrags.fastq1  $lib\.hist*");
  
  my %newreads = ('1' => $newread1,'2' => $newread2, 'u' => $newreadu);
  return(\%newreads);
}

sub reallyMergeReads {
  my ($lib,$orig,$reads,$base1,$base2,$InID) = @_;	
  
  my %reads = %{$reads};	
  my $newread1 = $orig->{'1'} . '_' . $base1 . '_p1';
  my $newread2 = $orig->{'2'} . '_' .$base1 .'_p2';
  my $newreadu = $orig->{'1'} . '_' .$base1 .'_u';
  my $call1 = system("cope -a $reads{'1'} -b $reads{'2'} -o $lib\.copemerged -2 $newread1 -3 $newread2 -m 0 -l 5 -c 0.9");
  die(qq/\nThe program "cope" is not in path! \n\n/) if ($call1 == -1 );
  
  my $newerread1 = $orig->{'1'} . '_' . $base2 . '_p1';
  my $newerread2 = $orig->{'2'} . '_' .$base2 .'_p2';
  my $newerreadu = $orig->{'1'} . '_' .$base2 .'_u';
  my %newerreads = ('1' => $newerread1,'2' => $newerread2, 'u' => $newerreadu);
 
  open (COPE, "<", $lib . ".copemerged");
  open (COPE_NEW, ">", $lib . ".copemerged1");  
  while (<COPE>) {
    chomp(my $line = $_);
    my @f = split /\s+/,$line;
    if (scalar(@f) == 4) {
      if ($f[1] =~ m/^($InID\S+)(\/1)(_$InID\S+)\/[1|2]$/) {
	my $full_id = $1 . $2 . $3;
	my $seq = <COPE>; my $qualid = <COPE>; my $qual = <COPE>;
	print COPE_NEW "@",$1,$2,"\n", $seq,$qualid,$qual;

      }
    }
    if (scalar(@f) == 2) {
      if ($f[0] =~ m/^@($InID\S+)(\/1)(_$InID\S+)\/[1|2]$/) {
	my $full_id = $1 . $2 . $3;
	my $seq = <COPE>; my $qualid = <COPE>; my $qual = <COPE>;
	print COPE_NEW "@",$1,$2,"\n", $seq,$qualid,$qual;	
      }
    }
  }
  close COPE;
  close COPE_NEW;
  
  my $call3 = system("cat $lib\.copemerged1 $reads{'u'} > $newerreadu");
  my $call4 = system("mv $newread1 $newerread1");
  my $call5 = system("mv $newread2 $newerread2");
  my $call6 = system("rm $lib\.copemerged $lib\.copemerged1 ");
  return(\%newerreads);
}


			
sub fixMatePair {
  my ($lib,$read,$readarray,$base) = @_;
  my @trim = @{$readarray};
  my %pair;	
  foreach my $reads (@trim) {
    open(IN, "<$reads");
    while(<IN>) {
      chomp(my $line = $_);
      if ($line =~ m/^@($InID\S+)\/[1|2]$/) {
	$pair{$1}++;
	chomp(my $seq = <IN>); chomp(my $qualid = <IN>); chomp(my $qual = <IN>);
      }
    }
    close(IN);	
  }
  my %reads = %{$read};
  my $out1 = $reads{'1'} . '_' . $base . '_p1';
  my $out2 = $reads{'2'} . '_' . $base . '_p2';
  my $outu = $reads{'1'} . '_' . $base . '_u';
  print $out1, "\n"; 
  print $out2, "\n";
  print $outu, "\n";
  open(OUT1, ">$out1"); 
  open(OUT2, ">$out2"); 
  open(OUTU, ">$outu"); 
  my %newpairs = ('1' => $out1, '2' => $out2, 'u' => $outu);
  foreach my $reads (@trim) {
    open(IN, "<$reads");
    my $file = $1 if $reads =~ m/_R(\d+)./;
    while(<IN>) {
      chomp(my $line = $_);	
      if ($line =~ m/^@($InID\S+)\/[1|2]$/) {
	my $id = $1;
	my $seq = <IN>;
	my $qualid = <IN>;
	my $qual = <IN>;
	if ($pair{$id} == 2) {
	  if ($file == 1) {
	    print OUT1 $line . "\n" . $seq . $qualid . $qual;
	  }
	  else {
	    print OUT2 $line . "\n" . $seq . $qualid . $qual;
	  }
	}
	else {
	  print OUTU $line . "\n" . $seq . $qualid . $qual;
	}
      }	
    }
    close(IN);	
  }	
  close(OUT1); close(OUT2); close(OUTU);	
  return(\%newpairs);
}

sub trimmomatic {
  my ($lib,$ad,$in,$base,$suffix) = @_;
  my $out = $base . '_' . $suffix;	
  my $adfile = $lib . "_adfile.fa";
  open(OUT, ">$adfile");
  foreach my $name (keys %{$ad}) {
    print OUT ">", $name, "\n", $ad->{$name}, "\n";
  }
  close(OUT);	
  my $call1 = system("java -classpath $trimmomatic org.usadellab.trimmomatic.TrimmomaticSE -phred33 $in $out ILLUMINACLIP:$adfile:2:40:15 SLIDINGWINDOW:4:20 MINLEN:$minLength LEADING:3 TRAILING:3");
  unlink($adfile);
  return($out);
}

sub cutadapt {
  my ($ad,$in,$base,$suffix) = @_;
  my $out  = $base . '_' . $suffix;	
  my $curRead = $in;
  my $tracker = 1;
  my %ad = %{$ad};
  foreach my $key (keys %ad) {
    my $out = $curRead . $tracker;
    my $call = system("$cutadapt -b $ad{$key} -O 4 -n 5 -e 0.15 -f fastq $curRead -o $out -m $minLength");
    die(qq/\nThe program "cutadapt" is not in path! \n\n/) if ($call == -1 );
    unlink($curRead) unless($curRead eq $in);
    $curRead = $out;
    $tracker++;
  }
  my $call2 = system("mv $curRead $out");
  return($out);
}

sub bowtie {
  my ($lib, $ad,$in,$base,$suffix) = @_;
  my $out  = $base . '_' . $suffix;	
  my $file = $lib."out.sam";
  
  my $adfile = $lib . "_adfile_norev.fa";
  open(OUT, ">$adfile");
  foreach my $name (keys %{$ad}) {
    print OUT ">", $name, "\n", $ad->{$name}, "\n" unless $name =~ m/rc/;
  }
  close(OUT);	
  
  my $bw2build = $bowtie . "-build";
  my $call1 = system("$bw2build $adfile $adfile");
  my $call2 = system("$bowtie --local -D 15 -R 2 -N 1 -L 10 -i S,1,0.75 -k 1 -x $adfile -U $in -S $file");
  die(qq/\nThe program "bowtie2" is not in path! \n\n/) if ($call2 == -1 );
  my $call3 = system("rm $adfile" . "*");
  
  open(IN, "<$file");
  open(OUT, ">$out");
  
  while(<IN>) {
    chomp(my $line1 = $_);
    my @d1 = split(/\t/,$line1);
    
    my $seq1 = $d1[9];
    my $qual1 = $d1[10];
    
    unless($line1 =~ m/^@/) {
      if ($line1 !~ m/\*/) {			
	if ($d1[5] =~ m/^(\d+)S\d+M$/) {
	  my $l = $1;
	  $seq1 = substr $seq1, 0, $l;
	  $qual1 = substr $qual1, 0, $l;
	}
	elsif ($d1[5] =~ m/^(\d+)M\d+S$/) {
	  my $start = $1;
	  $seq1 = substr $seq1, $start;
	  $qual1 = substr $qual1, $start;
	}
	else {
	  my @s;
	  while ($d1[5] =~ m/(\d+)S/g) {
	    push(@s,$1);
	  }
	  @s = sort {$a <=> $b} @s;
	  if ($s[$#s] >= $minLength) {	
	    if ($d1[5] =~ m/^(\S*)$s[$#s]/) {
	      my $match = $1;
	      my $length = $s[$#s];
	      my $start = 0;
	      while ($match =~ m/(\d+)/g) {
		$start += $1;
	      }
	      $seq1 = substr $seq1, $length;
	      $qual1 = substr $qual1, $length;	
	    }													
	  }
	  else {
	    $seq1 = 'N'; $qual1 = 'N';
	  }
	}
      }
      
      if (length($seq1) >= $minLength) {	
	print OUT "@" . $d1[0] . "\n" . $seq1 . "\n" . '+' . "\n" . $qual1 . "\n";	
      }
    }	
  }	
  unlink($file);
  close(IN); close(OUT);	
  return($out);
}
