#!/usr/bin/perl
use strict;
# This script is used to calculate the probability of a deltaSVM score higher or lower than null distribution
# This script is modified from deltasvm.pl, a simple script that calculates deltaSVM scores, which contributed by Lee et al. in their study: A method to predict the impact of regulatory variants from DNA sequence. Nat Genet. 47:955-961.

if ($#ARGV != 5) {
    print "\n";
    print "   Usage: testPosSelec.pl <ref_seq_file> <alt_seq_file> <svm_weight_file> <number_random_sample> <test> <output_file>\n\n";
    print "   ref_seq_file - reference sequences to score in FASTA format\n";
    print "   alt_seq_file - altenative sequences to score in FASTA format\n";
    print "   (NOTE: sequences should be in the same order, also have the same name)\n";
    print "   svm_weight_file - svm weights from gkm-SVM (i.e. all possible k-mer scores)\n";
    print "   number_random_sample -  how many random sequences with the same number of SNPs do you want to generate\n";
    print "   test -  highertail or lowertail test\n";
    print "   output_file - name of the output file\n\n";
    exit(0);
}

my $refseqf = $ARGV[0];
my $altseqf = $ARGV[1];
my $svmwf = $ARGV[2];
my $runNumb = $ARGV[3];
my $test = $ARGV[4];
my $outf = $ARGV[5];
my $kmerlen = 0;

my %rc;
$rc{"A"} = 'T'; $rc{"T"} = 'A'; $rc{"C"} = 'G'; $rc{"G"} = 'C';
sub revcomp {
    my $seq = shift @_;
    my $rcseq = "";
    my $l = length($seq)-1;
    while ($l >= 0) {
        $rcseq = $rcseq . $rc{substr($seq, $l, 1)};
        $l-=1;
    }

    return $rcseq;
}

my $line;

#2. read svmfile
print STDERR "reading svmwfile $svmwf..\n";
my %svmscore;
open IN, "<$svmwf";
my $bias;
chomp($bias=<IN>);
while (chomp($line=<IN>)) {
    if ($line=~/^#/) {next;}
    if ($line=~/^bias/) {next;}

    my @f = split /\s+/, $line;
    if ($#f == 1) {
        $svmscore{$f[0]} = $f[1];
        $svmscore{revcomp($f[0])} = $f[1];
    } elsif ($#f == 2) {
        $svmscore{$f[0]} = $f[2];
        $svmscore{$f[1]} = $f[2];
    } else {
        print STDERR "[ERROR] unkown format of svm weight file. quit.\n";
        exit(0);
    }
    if ($kmerlen == 0) {
        $kmerlen = length($f[0]);
    }
}
close IN;
#3. read fasta file
#Writing a DNA or protein sequence to a FASTA file in Perl is pretty easy and there are
#many ways to do so. Reading a FASTA file however, is more difficult because the sequence
#spans several lines and you don’t know you’ve reached the end of a sequence until you’ve
#reached the beginning of the next sequence
sub read_fasta_sequence {
   my ($fh, $seq_info) = @_;

   $seq_info->{seq} = undef; # clear out previous sequence

   # put the header into place
   $seq_info->{header} = $seq_info->{next_header} if $seq_info->{next_header};

   my $file_not_empty = 0;
   while (<$fh>) {
      $file_not_empty = 1;
      next if /^\s*$/;  # skip blank lines
      chomp;

      if (/^>/) { # fasta header line
         my $h = $_;
         $h =~ s/^>//;
         if ($seq_info->{header}) {
            $seq_info->{next_header} = $h;
            return $seq_info;
         }
         else { # first time through only
            $seq_info->{header} = $h;
         }
      }
      else {
         s/\s+//;  # remove any white space
         $seq_info->{seq} .= $_;
      }
   }

   if ($file_not_empty) {
      return $seq_info;
   }
   else {
      # clean everything up
      $seq_info->{header} = $seq_info->{seq} = $seq_info->{next_header} = undef;

      return;
   }
}
#4. generate random sequence based on number of snps
# make sure you call srand to seed the random number generator before you call this function.
srand(time|$$);
sub mutate {
  my($dna) = $_[0];
  my $n = $_[1]; # snp number
  # Pick random positions in the DNA
  my @position = (1 .. length($dna));
  my @randomPos;
  for (1 .. $n) {
  my $randomPos = $position[rand @position];
  push @randomPos, $randomPos;
  @position = grep { $_ != $randomPos } @position;# without replacement
  }
  # Pick a random nucleotide
  my(@nucs_A) = ('T', 'C', 'G');
  my(@nucs_T) = ('A', 'G', 'C');
  my(@nucs_C) = ('A', 'G', 'T');
  my(@nucs_G) = ('T', 'C', 'A');

  # change into new sequence
  for (my $i=0; $i < $n; $i++) {
    if ((substr($dna, $randomPos[$i], 1)=="A")) {
      substr($dna,$randomPos[$i],1,$nucs_A[rand @nucs_A]);
    }
    if ((substr($dna, $randomPos[$i], 1)=="T")) {
      substr($dna,$randomPos[$i],1,$nucs_T[rand @nucs_T]);
    }
    if ((substr($dna, $randomPos[$i], 1)=="C")) {
      substr($dna,$randomPos[$i],1,$nucs_C[rand @nucs_C]);
    }
    if ((substr($dna, $randomPos[$i], 1)=="G")) {
      substr($dna,$randomPos[$i],1,$nucs_G[rand @nucs_G]);
    }
  }
  splice( @randomPos );
  return $dna;
}

#5. scan sequence file and calculate probability
my $seqid_ref;
my $seqid_alt;
my $seq_ref;
my $seq_alt;
my $fh_ref;
my $fh_alt;
my %sequence_ref_data;
my %sequence_alt_data;
open($fh_ref, $refseqf) or die "can't open $refseqf: $!\n";
open($fh_alt, $altseqf) or die "can't open $altseqf: $!\n";
open OUT, ">$outf";

print STDERR "calculating deltaSVM using $refseqf and $altseqf..\n";
print STDERR "calculating probability of  deltaSVM higher or lower than null distribution..\n";

while (read_fasta_sequence($fh_ref, \%sequence_ref_data)) {
  $seqid_ref = $sequence_ref_data{header};
  $seq_ref = $sequence_ref_data{seq};
  read_fasta_sequence($fh_alt, \%sequence_alt_data);
  $seqid_alt = $sequence_alt_data{header};
  $seq_alt = $sequence_alt_data{seq};
  if ($seqid_ref ne $seqid_alt) {
      print STDERR "[ERROR] the sequence id in $refseqf ($seqid_ref) is different from the one in $altseqf ($seqid_alt). skip.\n";
  }
  # calculate snp number of two sequences
  my @seq_ref = split //,$seq_ref;
  my @seq_alt = split //,$seq_alt;
  my $j = 0;
  my $snpNumb=0;

  foreach  (@seq_ref) {
      if ($_ ne $seq_alt[$j]) {
          $snpNumb++;
      }
      $j++;
  }
  if ($snpNumb>1) { # at least two snps
    # calculate real deltaSVM
    my $score=0;
    my $seqlen = (length($seq_ref) > length($seq_alt)) ? length($seq_ref) : length($seq_alt);
    foreach my $i (0..($seqlen-$kmerlen)) {
      my $kmer_ref = substr($seq_ref, $i, $kmerlen);
      my $kmer_alt = substr($seq_alt, $i, $kmerlen);
      $score += ($svmscore{$kmer_alt} - $svmscore{$kmer_ref});
      #print OUT "$kmer_ref\t$kmer_alt\t$svmscore{$kmer_ref}\t$svmscore{$kmer_alt}\n";
    }
    # calculate random deltaSVM
    my $seq_random;
    my @score_random;
    my $score_random;
    for (1 .. $runNumb) {
      $score_random=0;
      $seq_random = mutate ($seq_ref,$snpNumb);
      foreach my $i (0..($seqlen-$kmerlen)) {
        my $kmer_ref = substr($seq_ref, $i, $kmerlen);
        my $kmer_alt = substr($seq_random, $i, $kmerlen);
        $score_random += ($svmscore{$kmer_alt} - $svmscore{$kmer_ref});
      }
      push @score_random, $score_random;
    }
    # calculate probabiliy of $score_abs >=/<= @score_random
    my $biggerNumb=0;
    my $lowerNumb=0;
    my $pvalue;
    if ($test eq "highertail") {
      foreach  (@score_random) {
          if ($_ >= $score) {
              $biggerNumb++;
          }
      }
      $pvalue=($biggerNumb+1)/($runNumb+1);
    }
    if ($test eq "lowertail") {
      foreach  (@score_random) {
          if ($_ <= $score) {
              $lowerNumb++;
          }
      }
      $pvalue=($lowerNumb+1)/($runNumb+1);
    }

    # output result
    #print "$seqid_ref","\n";
    #print "$score","\n";
    #print "$snpNumb","\n";
    #print "$pvalue","\n";

    print OUT "$seqid_ref","\t", "$score","\t","$snpNumb","\t","$pvalue","\n";
    splice(@score_random);
  }
}
close OUT;
close $fh_ref;
close $fh_alt;
