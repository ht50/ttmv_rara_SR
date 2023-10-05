#!/usr/bin/perl

use strict;
use warnings;
use List::Util qw(sum);
use List::Util qw(max);
use List::Util qw(min);
use Data::Dumper qw(Dumper);
use POSIX qw(ceil);
use POSIX qw(floor);
use File::Path qw(make_path);
use Getopt::Long qw(GetOptions);
use Getopt::Long qw(HelpMessage);

=pod
 
=head1 SYNOPSIS
 
  --inbam, -i	  Input bamfile (required)
  --output, -o	  Output path (required)
  --outkey, -k	  Output name stem (defaults to "", and then derived from inbam name minus bamsuffix)
  --bamx, -b	  Bamsuffix for outkey naming purposes (suffix is removed)
  --genome, -g    Genome build (defaults to "hg19"; or can be "hg38")
  --gregion, -gr  geneRegion for read extraction from bam (defaults to "", for RARA exon2+intron2+exon3 + 1kb buffers)
  --taxtype, -tt  Taxonomies for blastn search (defaults to "ttmv" but can also be "anellovirus")
  --btask, -bt    Blastn task for initial screen of reads (defaults to "megablast" but can be "blastn", "blastn-short", or "dc-megablast")
  --useunmapped, -uu	Option to include all unmapped reads and their paired reads in the blastn search (defaults to false)
  --suppvelvet, -sv	Option to supplement velvet input with unresolved soft-clipped reads to RARA geneRegion (default to true)
  --debug -d      Keep intermediate files for debugging
  --help, -h      Print this help
 
=cut
 
GetOptions(
  "inbam|i=s" => \(my $inbam = "" ),
  "output|o=s" => \(my $outpath = "" ),
  "outkey|k=s" => \(my $outkey = "" ),
  "bamx|b=s" => \(my $bamsuffix = ".bam" ),
  "genome|g=s" => \(my $genome = "hg19"),
  "gregion|gr=s" => \( my $geneRegion = "" ), 
  "taxtype|tt=s" => \( my $taxtype = "ttmv" ), # ttmv or anellovirus 
  "btask|bt=s" => \( my $btask = "megablast" ),
  "useunmapped|uu=s" => \( my $useunmapped = 0 ), 
  "suppvelvet|sv=s" => \( my $suppvelvet = 1 ),
  "debug|d" => \ (my $debug = 0 ),
  "help|h" => sub { HelpMessage(0) },
) or HelpMessage(1);
 
HelpMessage(1) unless( length($inbam)>0 );

# These commands should be in the executable path or changed to their location
my $samcmd = "samtools";	# tested with samtools version 1.17
my $blastncmd = "blastn";	# tested with BLAST+ version 2.13.0
my $vgcmd = "velvetg";		# tested with velvet version 1.2.10
my $vhcmd = "velveth";
my $bwacmd = "bwa";		# tested with bwa version 0.7.17-r1188

# These files should be present or changed to their locations
my $bwahg19 = "/reference_databases/ReferenceGenome/hg19/UCSC/hg19/Sequence/BWAIndex/genome.fa";
my $bwahg38 = "/reference_databases/ReferenceGenome/hg38/UCSC/hg38/Sequence/BWAIndex/genome.fa";
my $blastdb = "/lab-share/Path-LaMPP-e2/Public/resources-communal/blastdb/nt_viruses";
my $anelloIdFile = "anellovirus.txids";

my $taxid;
if( $taxtype eq "ttmv" ) {
  $taxid = 93678; 
} elsif( $taxtype eq "anellovirus" ) {
  $taxid = $anelloIdFile; 
} else {
  die "taxtype needs to be ttmv or anellovirus. Exiting...";
}

if( not -d $outpath ) { make_path $outpath or die "Failed to create path: $outpath"; }
if( not -d $outpath . "/results_SRnone" ) { make_path $outpath . "/results_SRnone" or die "Failed to create SRnone path"; }

my @cells;
if( $outkey eq "" ) {
  @cells = split(/\//, $inbam);
  $outkey = $cells[scalar(@cells)-1];
  $outkey =~ s/$bamsuffix//;
}
my $fbase = $outpath . "/" . $outkey;
my $fbaseSRnone = $outpath . "/results_SRnone/" . $outkey;

my $syscmd;
$syscmd = sprintf( "%s head %s | grep -w 'SN:17\\|SN:chr17' | awk '{print \$2}' | sed s/SN://g > %s_header_chr17.txt",
  $samcmd, $inbam, $fbase );
print($syscmd . "\n");
if( system($syscmd) ) { die "Failed to extract header. Exiting..."; }

my $chr17 = "";
open(FI, $fbase . "_header_chr17.txt" ) or die $!;
$chr17 = <FI>; chomp($chr17);
close FI;

if( $chr17 ne "chr17" && $chr17 ne "17" ) { die "Failed to extract chr17 header from bam. Exiting..."; }

my $bwaindex;
if( $genome eq "hg19" ) {
  if( $geneRegion eq "" ) { $geneRegion = $chr17 . ":38486109-38505716"; }
  $bwaindex = $bwahg19;
} elsif( $genome eq "hg38" ) {
  if( $geneRegion eq "" ) { $geneRegion = $chr17 . ":40329857-40349464"; }
  $bwaindex = $bwahg38; 
} else {
  die "genome needs to be hg19 or hg38. Exiting...";
}

print "************************\n" . $fbase . "\n************************\n";

$syscmd = sprintf( "time %s view %s %s | cut -f 1 | sort -u > %s_screen_reads.txt",
	$samcmd, $inbam, $geneRegion, $fbase );
print($syscmd . "\n");
if( system($syscmd) ) { die "Failed samtools extract by gene region. Exiting..."; }

if( $useunmapped ) {
  $syscmd = sprintf( "time %s view -f 4 %s | cut -f 1 | sort -u >> %s_screen_reads.txt",
	$samcmd, $inbam, $fbase );
  print($syscmd . "\n");
  if( system($syscmd) ) { die "Failed samtools extract unmapped. Exiting..."; }

  $syscmd = sprintf( "time sort -u %s_screen_reads.txt -o %s_screen_reads.txt", $fbase, $fbase );
  print($syscmd . "\n");
  if( system($syscmd) ) { die "Failed to sort into unique reads. Exiting..."; }
}

$syscmd = sprintf( "time %s view -b %s -N %s_screen_reads.txt -o %s_screen.bam",
	$samcmd, $inbam, $fbase, $fbase );
print($syscmd . "\n");
if( system($syscmd) ) { die "Failed samtools extract by read names. Exiting..."; }

$syscmd = sprintf( "time %s fasta %s_screen.bam > %s_screen.fa", $samcmd, $fbase, $fbase );
print($syscmd . "\n");
if( system($syscmd) ) { die "Failed samtools extract by read names. Exiting..."; }

my $blastout = sprintf( "%s_%s.out", $fbase, $taxtype );
if( $taxtype eq "ttmv" ) {
  $syscmd = sprintf( "time %s -query %s_screen.fa -db %s -task %s -taxids %s -outfmt 6 -out %s",
	$blastncmd, $fbase, $blastdb, $btask, $taxid, $blastout );
} else {
  $syscmd = sprintf( "time %s -query %s_screen.fa -db %s -task %s -taxidlist %s -outfmt 6 -out %s",
	$blastncmd, $fbase, $blastdb, $btask, $taxid, $blastout );
}
print($syscmd . "\n");
if( system($syscmd) ) { die "Failed blastn. Exiting..."; }

if( -z $blastout ) {
  print( "No blastn hits to " . $taxtype . ".\n" );
} else {
  $syscmd = sprintf( "time cut -f 1 %s | sed 's/\\\/1\$//g' | sed 's/\\\/2\$//g' | sort -u > %s_%s_reads.txt",
	$blastout, $fbase, $taxtype );
  print($syscmd . "\n");
  if( system($syscmd) ) { die "Failed to extract ttmv read names. Exiting..."; }

  $syscmd = sprintf( "time %s view -b %s_screen.bam -N %s_%s_reads.txt -o %s_%s.bam",
	$samcmd, $fbase, $fbase, $taxtype, $fbase, $taxtype );
  print($syscmd . "\n");
  if( system($syscmd) ) { die "Failed samtools extract by read names of ttmv blastn hits. Exiting..."; }

  $syscmd = sprintf( "time %s collate -O %s_%s.bam | %s fastq -1 %s_%s.R1.fq -2 %s_%s.R2.fq",
 	$samcmd, $fbase, $taxtype, $samcmd, $fbase, $taxtype, $fbase, $taxtype );
  print($syscmd . "\n");
  if( system($syscmd) ) { die "Failed samtools collate/fastq of ttmv blastn reads. Exiting..."; }

  $syscmd = sprintf( "%s mem -O 6 -T 20 %s %s_%s.R1.fq %s_%s.R2.fq > %s_%s_bwa.sam",
	$bwacmd, $bwaindex, $fbase, $taxtype, $fbase, $taxtype, $fbase, $taxtype );
  print($syscmd . "\n");
  if( system($syscmd) ) { die "Failed bwa mem on candidates. Exiting..."; }

  my %adets = (); my %bdets = (); my %cdets = (); my %seqs = ();
  my %seqsR1 = (); my %seqsR2 = (); my @cigars = (); my @cells2 = ();
  my $ttid; my $readname; my $cigar; my $seq; my $str; my $rnum;
  my $bkbin; my $bkbinRC; my $bk; my $bk2; my $cigstring; my $overlap;

  open( FI, $blastout ) or die $!;
  while(<FI>) {
    chomp;
    @cells = split( /\t/, $_ );
    @cells2 = split( /\//, $cells[0] );
    $bdets{ $cells[1] }{ $cells2[0] }{ "R" . $cells2[1] }{ $cells[6] . "_" . $cells[7] } =  
      $cells[1] . "|" . $cells[8] . "|" . $cells[9] . "|" . $cells[3] . "|" . $cells[4] . "|" . $cells[2] . "|" .
      $cells[10] . "|" . $cells[11];
  }
  close FI;

  open( FI, $fbase . "_" . $taxtype . "_bwa.sam" ) or die $!;
  while(<FI>) {
    if( !/^\@/ ) {
      chomp;
      my $text = $_;
      @cells = split( /\t/, $_ );
      $readname = $cells[0];
      $cigar = $cells[5];
      $seq = $cells[9];

      if( ( $cells[1] & 16 ) == 16 ) {
        $str = "-";
        $seq = reverse $seq; $seq =~ tr/ATGCatgc/TACGtacg/;
      } else {
        $str = "+";
      }

      if( ( $cells[1] & 64 ) == 64 ) {
        $rnum = "R1";
        if( ( $cells[1] & 256 ) != 256 && ( $cells[1] & 2048 ) != 2048 ) {
          $seqs{$readname."/R1"} = $seq;
        }
      } elsif( ( $cells[1] & 128 ) == 128 ) {
        $rnum = "R2";
        if( ( $cells[1] & 256 ) != 256 && ( $cells[1] & 2048 ) != 2048 ) {
          $seqs{$readname."/R2"} = $seq;
        }
      }

      my $editdist = "NA";
      my $sizeleftclip = 0; my $sizerightclip = 0;
      if( $cigar =~ /^(\d+)([SH])(.*)/ ) { $sizeleftclip = $1; }
      if( $cigar =~ /(\d+)([SH])$/ ) { $sizerightclip = $1; }
      my $loc = $cells[3];
      my $lpos = $loc; my $rpos = $loc - 1; my $sizeseq = 0;
      while( $cigar =~ /(\d+)([SHDIMN])(.*)/ ) {
        $cigar = $3;
        if( $2 ne "D" ) { $sizeseq += $1; }
        if( $2 eq "M" || $2 eq "N" || $2 eq "D" ) { $rpos += $1; }
      }
      $cigar = $cells[5];
      if( $cigar =~ /^(\d+)([SH])(.*)/ ) { $cigar = $3; }
      if( $cigar =~ /(.*)([MDI])(\d+)([SH])$/ ) { $cigar = $1 . $2; }
      if ($text =~ /NM:i:(\d+)/ ) { $editdist = $1; }
      if( $str eq "+" && $cells[2] eq "chr17" ) {
        $adets{$readname}{$rnum}{($sizeleftclip + 1) . "_" . ($sizeseq - $sizerightclip)} =
	  $cells[2] . "|" . $lpos . "|" . $rpos . "|+|" . $editdist . "|" . $cigar . "|" . $cells[5];
      } elsif( $str eq "-" && $cells[2] eq "chr17" ) { 
        $adets{$readname}{$rnum}{($sizerightclip + 1) . "_" . ($sizeseq - $sizeleftclip)} =
	  $cells[2] . "|" . $rpos . "|" . $lpos . "|-|" . $editdist . "|" . $cigar . "|" . $cells[5];
      }
    }
  }
  close FI;

  foreach $readname ( keys %adets ) {
    foreach $ttid ( keys %bdets ) {
      if( exists( $bdets{$ttid}{$readname} ) ) {
        foreach $rnum ( keys %{$adets{$readname}} ) {
          if( exists( $bdets{$ttid}{$readname}{$rnum} ) && !exists( $adets{$readname}{$rnum}{"1_0"} ) ) {
            foreach( keys %{$adets{$readname}{$rnum}} ) {
              $cdets{$ttid}{$readname}{$rnum}{$_} = $adets{$readname}{$rnum}{$_};
            }
            foreach( keys %{$bdets{$ttid}{$readname}{$rnum}} ) {
              $cdets{$ttid}{$readname}{$rnum}{$_} = $bdets{$ttid}{$readname}{$rnum}{$_};
            }
          }
        }
      }
    }
  }

  my %cigsSR = (); my %readsSR = (); my %aggSR = ();
  foreach $ttid ( keys %cdets ) {
    foreach $readname ( keys %{$cdets{$ttid}} ) {
      foreach $rnum ( ("R1", "R2") ) {
        if( exists( $cdets{$ttid}{$readname}{$rnum} ) ) {
          @cigars = keys %{$cdets{$ttid}{$readname}{$rnum}};
          if( scalar(@cigars) > 1 ) {
            my %x0; my %x1;
            foreach $cigar (@cigars) {
              @cells = split( /_/, $cigar );
              $x0{$cigar} = $cells[0]; $x1{$cigar} = $cells[1];
            }
            my %badcigs = ();
            foreach $cigar (@cigars) {
              foreach my $tmp (@cigars) {
                if( $cigar ne $tmp && $x0{$cigar} >= $x0{$tmp} && $x1{$cigar} <= $x1{$tmp} ) {
                  $badcigs{$cigar} = 1;
                }
              }
            }
            foreach $cigar ( keys %badcigs ) { delete $x0{$cigar}; delete $x1{$cigar}; }
            my @goodcigs = sort { $x0{$a} <=> $x0{$b} } keys %x0;

            @cells = split(/\|/, $cdets{$ttid}{$readname}{$rnum}{$goodcigs[0]} );
            $cigstring = $goodcigs[0] . "_" . $cells[4];
            if( scalar( @goodcigs ) > 1 ) {
              $bkbin = "";
              for( my $i=0; $i < scalar(@goodcigs)-1; $i++ ) {
                @cells = split(/\|/, $cdets{$ttid}{$readname}{$rnum}{$goodcigs[$i]} );
                @cells2 = split(/\|/, $cdets{$ttid}{$readname}{$rnum}{$goodcigs[$i+1]} );
                $bk = $cells[0] . ":" . $cells[2]; $bk2 = $cells2[0] . ":" . $cells2[1];
                $overlap = $x1{$goodcigs[$i]} - $x0{$goodcigs[$i+1]}+1;
                my $tmp; my $tmpRC;
                if( $cells[2]>$cells[1] ) {
                  $tmp = $bk . "(+)_" . $bk2 . "(+)_overlap_" . $overlap; 
                  $tmpRC = $bk2 . "(-)_" . $bk . "(-)_overlap_" . $overlap;
                } else {
                  $tmp = $bk . "(-)_" . $bk2 . "(-)_overlap_" . $overlap; 
                  $tmpRC = $bk2 . "(+)_" . $bk . "(+)_overlap_" . $overlap;
                } 
                if( $i==0 ) { $bkbin = $tmp; } else { $bkbin .= "|" . $tmp; }
                if( $i==0 ) { $bkbinRC = $tmpRC; } else { $bkbinRC = $tmpRC . "|" . $bkbinRC; }
                $cigstring .= "|" . $goodcigs[$i+1] . "_" . $cells2[4];
              }
              my $ctbk = $bkbin =~ tr/+//; my $ctbkRC = $bkbinRC =~ tr/+//;
              if( $ctbk > $ctbkRC ) {
                $cigsSR{$ttid}{$bkbin}{$readname."/".$rnum} = $cigstring . "|" . length($seqs{$readname."/".$rnum}) . "|+";
                if( !exists( $readsSR{$ttid}{$bkbin}{$readname} )) {
                  $readsSR{$ttid}{$bkbin}{$readname} = 1;
                  $aggSR{$ttid}{$bkbin} += 1;
                }  
              } else {
                $cigsSR{$ttid}{$bkbinRC}{$readname."/".$rnum} = $cigstring . "|" . length($seqs{$readname."/".$rnum}) . "|-";
                if( !exists( $readsSR{$ttid}{$bkbinRC}{$readname} )) {
                  $readsSR{$ttid}{$bkbinRC}{$readname} = 1;
                  $aggSR{$ttid}{$bkbinRC} += 1;
                }
              }
            }
          }
        }
      }
    }
  }

  my %totSR = ();
  foreach $ttid (keys %aggSR) { $totSR{$ttid} = sum(values %{$aggSR{$ttid}}); }
  open( FO, ">" . $fbase . "_results_SR_all.out" ) or die $!;
  print FO "GenBank_id\tSR\tJunction\n";
  foreach $ttid ( sort {$totSR{$b}<=>$totSR{$a}} keys %totSR ) {
    print FO $ttid . "\t" . $totSR{$ttid} . "\tall\n";
  }
  foreach $ttid ( sort {$totSR{$b}<=>$totSR{$a}} keys %totSR ) {
    foreach $bkbin ( sort {$aggSR{$ttid}{$b}<=>$aggSR{$ttid}{$a}} keys %{$aggSR{$ttid}} ) {
      print FO $ttid . "\t" . $aggSR{$ttid}{$bkbin} . "\t" . $bkbin . "\n";
    }
  }
  close FO;

  foreach $ttid ( sort {$totSR{$b}<=>$totSR{$a}} keys %totSR ) {
    open( FO, ">" . $fbase . "_results_SR_" . $ttid . ".out" ) or die $!;
    foreach $bkbin ( sort {$aggSR{$ttid}{$b}<=>$aggSR{$ttid}{$a}} keys %{$aggSR{$ttid}} ) {
      foreach $readname ( sort keys %{$cigsSR{$ttid}{$bkbin}} ) {
        print FO $readname . "\t" . $cigsSR{$ttid}{$bkbin}{$readname} . "\t" . $bkbin . "\t" . $seqs{$readname} . "\n";
      }
    }
    close FO;
  }

  my %cigsSRB = (); my %readsSRB = (); my %aggSRB = ();
  my %cigsSRN = (); my %readsSRN = (); my %aggSRN = ();
  foreach $ttid ( keys %bdets ) {
    foreach $readname ( keys %{$bdets{$ttid}} ) {
      foreach $rnum ( ("R1", "R2") ) {
        if( exists( $bdets{$ttid}{$readname}{$rnum} ) &&
            !exists( $cdets{$ttid}{$readname}{$rnum} ) ) {
          @cigars = keys %{$bdets{$ttid}{$readname}{$rnum}};
          if( scalar(@cigars) > 0 ) {
            my %x0; my %x1;
            foreach $cigar (@cigars) {
              @cells = split( /_/, $cigar );
              $x0{$cigar} = $cells[0]; $x1{$cigar} = $cells[1];
            }
            my %badcigs = ();
            foreach $cigar (@cigars) {
              foreach my $tmp (@cigars) {
                if( $cigar ne $tmp && $x0{$cigar} >= $x0{$tmp} && $x1{$cigar} <= $x1{$tmp} ) {
                  $badcigs{$cigar} = 1;
                }
              }
            }
            foreach $cigar ( keys %badcigs ) { delete $x0{$cigar}; delete $x1{$cigar}; }
            my @goodcigs = sort { $x0{$a} <=> $x0{$b} } keys %x0;

            @cells = split(/\|/, $bdets{$ttid}{$readname}{$rnum}{$goodcigs[0]} );
            $cigstring = $goodcigs[0] . "_" . $cells[4];
            if( scalar( @goodcigs ) > 1 ) {
              $bkbin = "";
              for( my $i=0; $i < scalar(@goodcigs)-1; $i++ ) {
                @cells = split(/\|/, $bdets{$ttid}{$readname}{$rnum}{$goodcigs[$i]} );
                @cells2 = split(/\|/, $bdets{$ttid}{$readname}{$rnum}{$goodcigs[$i+1]} );
                $bk = $cells[0] . ":" . $cells[2]; $bk2 = $cells2[0] . ":" . $cells2[1];
                $overlap = $x1{$goodcigs[$i]} - $x0{$goodcigs[$i+1]}+1;
                my $tmp; my $tmpRC;
                if( $cells[2]>$cells[1] ) {
                  $tmp = $bk . "(+)_" . $bk2 . "(+)_overlap_" . $overlap; 
                  $tmpRC = $bk2 . "(-)_" . $bk . "(-)_overlap_" . $overlap;
                } else {
                  $tmp = $bk . "(-)_" . $bk2 . "(-)_overlap_" . $overlap; 
                  $tmpRC = $bk2 . "(+)_" . $bk . "(+)_overlap_" . $overlap;
                } 
                if( $i==0 ) { $bkbin = $tmp; } else { $bkbin .= "|" . $tmp; }
                if( $i==0 ) { $bkbinRC = $tmpRC; } else { $bkbinRC = $tmpRC . "|" . $bkbinRC; }
                $cigstring .= "|" . $goodcigs[$i+1] . "_" . $cells2[4];
              }
              my $ctbk = $bkbin =~ tr/+//; my $ctbkRC = $bkbinRC =~ tr/+//;
              if( $ctbk > $ctbkRC ) {
                $cigsSRB{$ttid}{$bkbin}{$readname."/".$rnum} = $cigstring . "|" . length($seqs{$readname."/".$rnum}) . "|+";
                if( !exists( $readsSRB{$ttid}{$bkbin}{$readname} )) {
                  $readsSRB{$ttid}{$bkbin}{$readname} = 1;
                  $aggSRB{$ttid}{$bkbin} += 1;
                }  
              } else {
                $cigsSRB{$ttid}{$bkbinRC}{$readname."/".$rnum} = $cigstring . "|" . length($seqs{$readname."/".$rnum}) . "|-";
                if( !exists( $readsSRB{$ttid}{$bkbinRC}{$readname} )) {
                  $readsSRB{$ttid}{$bkbinRC}{$readname} = 1;
                  $aggSRB{$ttid}{$bkbinRC} += 1;
                }
              }
            } else {
              $bkbin = $ttid . ":" . $cells[1] . "-" . $cells[2];
              $cigsSRN{$ttid}{$bkbin}{$readname."/".$rnum} = $cigstring . "|" . length($seqs{$readname."/".$rnum});
              if( !exists( $readsSRN{$ttid}{$bkbin}{$readname} )) {
                $readsSRN{$ttid}{$bkbin}{$readname} = 1;
                $aggSRN{$ttid}{$bkbin} += 1;
              }
            }
          }
        }
      }
    }
  }

  my %totSRB = ();
  foreach $ttid (keys %aggSRB) { $totSRB{$ttid} = sum(values %{$aggSRB{$ttid}}); }
  open( FO, ">" . $fbase . "_results_SRnohg_all.out" ) or die $!;
  print FO "GenBank_id\tSR\tJunction\n";
  foreach $ttid ( sort {$totSRB{$b}<=>$totSRB{$a}} keys %totSRB ) {
    print FO $ttid . "\t" . $totSRB{$ttid} . "\tall\n";
  }
  foreach $ttid ( sort {$totSRB{$b}<=>$totSRB{$a}} keys %totSRB ) {
    foreach $bkbin ( sort {$aggSRB{$ttid}{$b}<=>$aggSRB{$ttid}{$a}} keys %{$aggSRB{$ttid}} ) {
      print FO $ttid . "\t" . $aggSRB{$ttid}{$bkbin} . "\t" . $bkbin . "\n";
    }
  }
  close FO;

  foreach $ttid ( sort {$totSRB{$b}<=>$totSRB{$a}} keys %totSRB ) {
    open( FO, ">" . $fbase . "_results_SRnohg_" . $ttid . ".out" ) or die $!;
    foreach $bkbin ( sort {$aggSRB{$ttid}{$b}<=>$aggSRB{$ttid}{$a}} keys %{$aggSRB{$ttid}} ) {
      foreach $readname ( sort keys %{$cigsSRB{$ttid}{$bkbin}} ) {
        print FO $readname . "\t" . $cigsSRB{$ttid}{$bkbin}{$readname} . "\t" . $bkbin . "\t" . $seqs{$readname} . "\n";
      }
    }
    close FO;
  }

  my %totSRN = ();
  foreach $ttid (keys %aggSRN) { $totSRN{$ttid} = sum(values %{$aggSRN{$ttid}}); }
  open( FO, ">" . $fbase . "_results_SRnone_all.out" ) or die $!;
  print FO "GenBank_id\tSR\tJunction\n";
  foreach $ttid ( sort {$totSRN{$b}<=>$totSRN{$a}} keys %totSRN ) {
    print FO $ttid . "\t" . $totSRN{$ttid} . "\tall\n";
  }
  foreach $ttid ( sort {$totSRN{$b}<=>$totSRN{$a}} keys %totSRN ) {
    foreach $bkbin ( sort {$aggSRN{$ttid}{$b}<=>$aggSRN{$ttid}{$a}} keys %{$aggSRN{$ttid}} ) {
      print FO $ttid . "\t" . $aggSRN{$ttid}{$bkbin} . "\t" . $bkbin . "\n";
    }
  }
  close FO;

  foreach $ttid ( sort {$totSRN{$b}<=>$totSRN{$a}} keys %totSRN ) {
    open( FO, ">" . $fbaseSRnone . "_results_SRnone_" . $ttid . ".out" ) or die $!;
    foreach $bkbin ( sort {$aggSRN{$ttid}{$b}<=>$aggSRN{$ttid}{$a}} keys %{$aggSRN{$ttid}} ) {
      foreach $readname ( sort keys %{$cigsSRN{$ttid}{$bkbin}} ) {
        print FO $readname . "\t" . $cigsSRN{$ttid}{$bkbin}{$readname} . "\t" . $bkbin . "\t" . $seqs{$readname} . "\n";
      }
    }
    close FO;
  }

  my $fbaseV = $fbase . "_vel"; 
  my %vreads = ();
  open(FI, $fbase . "_" . $taxtype . "_reads.txt") or die $!;
  while(<FI>) { chomp; $vreads{$_} = 1; }
  close FI;

  if( $suppvelvet ) {
    $syscmd = sprintf( "time %s view %s_screen.bam > %s_screen.sam", $samcmd, $fbase, $fbase );
    print($syscmd . "\n");
    if( system($syscmd) ) { die "Failed samtools bam to sam. Exiting..."; }

    my %dets = (); my %pkeys = ();
    my $pair; my $rtype; my $strand;
    open(FI, $fbase . "_screen.sam") or die $!;
    while(<FI>) {
      chomp;
      @cells = split(/\t/);
      if( ( $cells[1] & 64 ) == 64 ) {
        $pair = "R1";
      } elsif ( ( $cells[1] & 128 ) == 128 ) {
        $pair = "R2";
      } else {
      $pair = "R0";
      }
      if( ( $cells[1] & 256 ) != 256 && ( $cells[1] & 2048 ) != 2048 ) {
        $rtype = "primary";
      } elsif( ( $cells[1] & 2048 ) == 2048 ) {
        $rtype = "supplementary";
      } elsif( ( $cells[1] & 256 ) == 256 ) {
        $rtype = "secondary";
      } else {
        $rtype = "other";
      }
      if( ( $cells[1] & 16 ) == 16 ) { $strand = "-"; } else { $strand = "+"; }
      $dets{ $cells[0] }{ $pair }{ $rtype . "_" . $cells[2] . "_" . $cells[3] . "_" . $cells[5] . "_" . $strand } = $cells[9];
    }
    close FI;

    my $read; my $rkey; my $cigar; my $buffer = 20; my $mt = 30;
    foreach $read ( keys %dets ) {
      foreach $pair ( keys %{$dets{$read}} ) {
        my $alnLeft = 0; my $alnRight = 0;
        foreach $rkey ( keys %{$dets{$read}{$pair}} ) {
          @cells = split( /_/, $rkey );
          $cigar = $cells[3];
          if( ( $cigar =~ /^(\d+)M/ && $cells[4] eq "+" ) || ( $cigar =~ /M$/ && $cells[4] eq "-" ) ) { $alnLeft = 1; }
          if( ( $cigar =~ /^(\d+)M/ && $cells[4] eq "-" ) || ( $cigar =~ /M$/ && $cells[4] eq "+" ) ) { $alnRight = 1; }
          if( ( $cigar =~ /^(\d+)S(\d+)M/ && $1<=$buffer && $2 >= $mt && $cells[4] eq "+" ) ||
              ( $cigar =~ /(\d+)M(\d+)S$/ && $2<=$buffer && $1 >= $mt && $cells[4] eq "-" ) ) { $alnLeft = 1; }
          if( ( $cigar =~ /^(\d+)S(\d+)M/ && $1<=$buffer && $2 >= $mt && $cells[4] eq "-" ) ||
              ( $cigar =~ /(\d+)M(\d+)S$/ && $2<=$buffer && $1 >= $mt && $cells[4] eq "+" ) ) { $alnRight = 1; }
          if( ( $cigar =~ /^(\d+)H(\d+)M/ && $1<=$buffer && $2 >= $mt && $cells[4] eq "+" ) ||
              ( $cigar =~ /(\d+)M(\d+)H$/ && $2<=$buffer && $1 >= $mt && $cells[4] eq "-" ) ) { $alnLeft = 1; }
          if( ( $cigar =~ /^(\d+)H(\d+)M/ && $1<=$buffer && $2 >= $mt && $cells[4] eq "-" ) ||
              ( $cigar =~ /(\d+)M(\d+)H$/ && $2<=$buffer && $1 >= $mt && $cells[4] eq "+" ) ) { $alnRight = 1; }
        }
        if( $alnLeft == 0 || $alnRight == 0 ) { $vreads{$read} = 1; }
      }
    }
  }

  open(FO, ">" . $fbaseV . "_reads.txt" ) or die $!;
  foreach( keys %vreads ) { print FO $_ . "\n"; }
  close FO; 

  $syscmd = sprintf( "time %s view -b %s_screen.bam -N %s_reads.txt -o %s.bam", $samcmd, $fbase, $fbaseV, $fbaseV );
  print($syscmd . "\n");
  if( system($syscmd) ) { die "Failed samtools extract by velvet read names. Exiting..."; }

  $syscmd = sprintf( "time %s collate -O %s.bam | %s fastq -1 %s.R1.fq -2 %s.R2.fq",
 	$samcmd, $fbaseV, $samcmd, $fbaseV, $fbaseV );
  print($syscmd . "\n");
  if( system($syscmd) ) { die "Failed samtools conversion to fastq. Exiting..."; }

  $syscmd = sprintf( "time %s %s.vel 21 -fastq -shortPaired %s.R1.fq %s.R2.fq", $vhcmd, $fbaseV, $fbaseV, $fbaseV );
  print($syscmd . "\n");
  if( system($syscmd) ) { die "Failed velveth. Exiting..."; }

  my @typesV = ( "vcc100", "vcc50", "vcc20", "vcc10", "vbase", "vauto" );

  foreach my $typeV (@typesV) { 
    my $fastaV = sprintf( "%s_%s.fa", $fbaseV, $typeV );
    my $blastoutV = sprintf( "%s_%s.out", $fbaseV, $typeV );
    my $samV = sprintf( "%s_%s.sam", $fbaseV, $typeV );

    if( $typeV eq "vbase" ) {
      $syscmd = sprintf( "time %s %s.vel", $vgcmd, $fbaseV );
    } elsif( $typeV eq "vauto" ) {
      $syscmd = sprintf( "time %s %s.vel -cov_cutoff auto", $vgcmd, $fbaseV );
    } else {
      my $vcc = $typeV; $vcc =~ s/vcc//;
      $syscmd = sprintf( "time %s %s.vel -cov_cutoff %s -min_contig_lgth 200", $vgcmd, $fbaseV, $vcc );
    } 
    print($syscmd . "\n");
    if( system($syscmd) ) { die "Failed velvetg. Exiting..."; }
    $syscmd = sprintf( "mv %s.vel/contigs.fa %s", $fbaseV, $fastaV );
    print($syscmd . "\n");
    if( system($syscmd) ) { die "Failed to move velvetg output. Exiting..."; }

    if( -z $fastaV ) { next; }

    if( $taxtype eq "ttmv" ) {
      $syscmd = sprintf( "time %s -query %s_%s.fa -db %s -taxids %s -outfmt 6 -out %s",
  	$blastncmd, $fbaseV, $typeV, $blastdb, $taxid, $blastoutV );
    } else {
      $syscmd = sprintf( "time %s -query %s_%s.fa -db %s -taxidlist %s -outfmt 6 -out %s",
	$blastncmd, $fbaseV, $typeV, $blastdb, $taxid, $blastoutV );
    }
    print($syscmd . "\n");
    if( system($syscmd) ) { die "Failed blastn. Exiting..."; }

    if( -z $blastoutV ) { next; }

    $syscmd = sprintf( "%s mem -O 6 -T 20 %s %s > %s", $bwacmd, $bwaindex, $fastaV, $samV );
    print($syscmd . "\n");
    if( system($syscmd) ) { die "Failed bwa mem. Exiting..."; }

    my %adetsV = (); my %bdetsV = (); my %cdetsV = (); my %seqsV = ();
    open( FI, $blastoutV ) or die $!;
    while(<FI>) {
      chomp;
      @cells = split( /\t/, $_ );
      $bdetsV{ $cells[1] }{ $cells[0] }{ $cells[6] . "_" . $cells[7] } =  
        $cells[1] . "|" . $cells[8] . "|" . $cells[9] . "|" . $cells[3] . "|" . $cells[4] . "|" . $cells[2] . "|" .
        $cells[10] . "|" . $cells[11];
    }
    close FI;

    open( FI, $samV ) or die $!;
    while(<FI>) {
      if( !/^\@/ ) {
        chomp;
        my $text = $_;
        @cells = split( /\t/, $_ );
        $readname = $cells[0];
        $cigar = $cells[5];
        $seq = $cells[9];

        if( ( $cells[1] & 16 ) == 16 ) {
          $str = "-";
          $seq = reverse $seq; $seq =~ tr/ATGCatgc/TACGtacg/;
        } else {
          $str = "+";
        }

        if( ( $cells[1] & 256 ) != 256 && ( $cells[1] & 2048 ) != 2048 ) {
          $seqsV{$readname} = $seq;
        }

        my $editdist = "NA";
        my $sizeleftclip = 0; my $sizerightclip = 0;
        if( $cigar =~ /^(\d+)([SH])(.*)/ ) { $sizeleftclip = $1; }
        if( $cigar =~ /(\d+)([SH])$/ ) { $sizerightclip = $1; }
        my $loc = $cells[3];
        my $lpos = $loc; my $rpos = $loc - 1; my $sizeseq = 0;
        while( $cigar =~ /(\d+)([SHDIMN])(.*)/ ) {
          $cigar = $3;
          if( $2 ne "D" ) { $sizeseq += $1; }
          if( $2 eq "M" || $2 eq "N" || $2 eq "D" ) { $rpos += $1; }
        }
        $cigar = $cells[5];
        if( $cigar =~ /^(\d+)([SH])(.*)/ ) { $cigar = $3; }
        if( $cigar =~ /(.*)([MDI])(\d+)([SH])$/ ) { $cigar = $1 . $2; }
        if ($text =~ /NM:i:(\d+)/ ) { $editdist = $1; }
        if( $str eq "+" && $cells[2] eq "chr17" ) {
          $adetsV{$readname}{($sizeleftclip + 1) . "_" . ($sizeseq - $sizerightclip)} =
	    $cells[2] . "|" . $lpos . "|" . $rpos . "|+|" . $editdist . "|" . $cigar . "|" . $cells[5];
        } elsif( $str eq "-" && $cells[2] eq "chr17" ) { 
          $adetsV{$readname}{($sizerightclip + 1) . "_" . ($sizeseq - $sizeleftclip)} =
	    $cells[2] . "|" . $rpos . "|" . $lpos . "|-|" . $editdist . "|" . $cigar . "|" . $cells[5];
        }
      }
    }
    close FI;

    foreach $readname ( keys %adetsV ) {
      foreach $ttid ( keys %bdetsV ) {
        if( exists( $bdetsV{$ttid}{$readname} ) && !exists( $adetsV{$readname}{"1_0"} ) ) {
          foreach( keys %{$adetsV{$readname}} ) {
            $cdetsV{$ttid}{$readname}{$_} = $adetsV{$readname}{$_};
          }
          foreach( keys %{$bdetsV{$ttid}{$readname}} ) {
            $cdetsV{$ttid}{$readname}{$_} = $bdetsV{$ttid}{$readname}{$_};
          }
        }
      }
    }

    my %cigsSRV = (); my %readsSRV = (); my %aggSRV = ();
    my %scoreSRV = (); my %mmSRV = (); my %smSRV = ();
    foreach $ttid ( keys %cdetsV ) {
      foreach $readname ( keys %{$cdetsV{$ttid}} ) {
        if( exists( $cdetsV{$ttid}{$readname} ) ) {
          @cigars = keys %{$cdetsV{$ttid}{$readname}};
          if( scalar(@cigars) > 1 ) {
            my %x0; my %x1;
            foreach $cigar (@cigars) {
              @cells = split( /_/, $cigar );
              $x0{$cigar} = $cells[0]; $x1{$cigar} = $cells[1];
            }
            my %badcigs = ();
            foreach $cigar (@cigars) {
              foreach my $tmp (@cigars) {
                if( $cigar ne $tmp && $x0{$cigar} >= $x0{$tmp} && $x1{$cigar} <= $x1{$tmp} ) {
                  $badcigs{$cigar} = 1;
                }
              }
            }
            foreach $cigar ( keys %badcigs ) { delete $x0{$cigar}; delete $x1{$cigar}; }
            my @goodcigs = sort { $x0{$a} <=> $x0{$b} } keys %x0;

            $cigstring = ""; my $score = 0; my $mm = 0;
            if( scalar( @goodcigs ) > 1 ) {
              $bkbin = "";
              for( my $i=0; $i < scalar(@goodcigs)-1; $i++ ) {
                @cells = split(/\|/, $cdetsV{$ttid}{$readname}{$goodcigs[$i]} );
                @cells2 = split(/\|/, $cdetsV{$ttid}{$readname}{$goodcigs[$i+1]} );
                $bk = $cells[0] . ":" . $cells[2]; $bk2 = $cells2[0] . ":" . $cells2[1];
                $overlap = $x1{$goodcigs[$i]} - $x0{$goodcigs[$i+1]}+1;
                my $tmp; my $tmpRC;
                if( $cells[2]>$cells[1] ) {
                  $tmp = $bk . "(+)_" . $bk2 . "(+)_overlap_" . $overlap; 
                  $tmpRC = $bk2 . "(-)_" . $bk . "(-)_overlap_" . $overlap;
                } else {
                  $tmp = $bk . "(-)_" . $bk2 . "(-)_overlap_" . $overlap; 
                  $tmpRC = $bk2 . "(+)_" . $bk . "(+)_overlap_" . $overlap;
                } 
                if( $i==0 ) {
                  $bkbin = $tmp;
                  $bkbinRC = $tmpRC;
                  $cigstring .= $goodcigs[0] . "_" . $cells[4];
                  $mm += $cells[4];
                  @cells = split(/_/, $goodcigs[0]);
                  $score += $cells[1] - $cells[0] + 1;
                } else {
                  $bkbin .= "|" . $tmp;
                  $bkbinRC = $tmpRC . "|" . $bkbinRC;
                }
                $cigstring .= "|" . $goodcigs[$i+1] . "_" . $cells2[4];
                $mm += $cells2[4];
                @cells2 = split(/_/, $goodcigs[$i+1]);
                $score += $cells2[1] - $cells2[0] + 1;
                if( $overlap > 0 ) { $score -= $overlap; }
              }
              my $ctbk = $bkbin =~ tr/+//; my $ctbkRC = $bkbinRC =~ tr/+//;
              if( $ctbk > $ctbkRC ) {
                $cigsSRV{$ttid}{$bkbin}{$readname} = $cigstring . "|" . length($seqsV{$readname}) . "|+";
                $scoreSRV{$readname.",".$ttid.",".$bkbin} = $score;
                $mmSRV{$readname.",".$ttid.",".$bkbin} = $mm;
                $smSRV{$readname.",".$ttid.",".$bkbin} = $score - $mm; 
                if( !exists( $readsSRV{$ttid}{$bkbin}{$readname} )) {
                  $readsSRV{$ttid}{$bkbin}{$readname} = 1;
                  $aggSRV{$ttid}{$bkbin} += 1;
                }  
              } else {
                $cigsSRV{$ttid}{$bkbinRC}{$readname} = $cigstring . "|" . length($seqsV{$readname}) . "|-";
                $scoreSRV{$readname.",".$ttid.",".$bkbinRC} = $score;
                $mmSRV{$readname.",".$ttid.",".$bkbinRC} = $mm;
                $smSRV{$readname.",".$ttid.",".$bkbinRC} = $score - $mm; 
                if( !exists( $readsSRV{$ttid}{$bkbinRC}{$readname} )) {
                  $readsSRV{$ttid}{$bkbinRC}{$readname} = 1;
                  $aggSRV{$ttid}{$bkbinRC} += 1;
                }
              }
            }
          }
        }
      }
    }

    open( FO, ">" . $fbase . "_results_vel_" . $typeV . ".out" ) or die $!;
    foreach my $skey ( sort {$smSRV{$b}<=>$smSRV{$a}} keys %smSRV ) {
      @cells = split(/,/, $skey);
      $readname = $cells[0]; $ttid = $cells[1]; $bkbin = $cells[2];
      print $ttid . "\t" . $readname . "\t" . length($seqsV{$readname}) . "\t" . $scoreSRV{$skey} . "\t" . $mmSRV{$skey} .
	"\t" . $cigsSRV{$ttid}{$bkbin}{$readname} . "\t" . $bkbin . "\n";
      print FO $ttid . "\t" . $readname . "\t" . length($seqsV{$readname}) . "\t" . $scoreSRV{$skey} . "\t" . $mmSRV{$skey} .
	"\t" . $cigsSRV{$ttid}{$bkbin}{$readname} . "\t" . $bkbin . "\n";
    }
    close FO;
  }

  if( !$debug ) {
    $syscmd = sprintf( "rm -rf %s.vel", $fbaseV );
    print($syscmd . "\n");
    if( system($syscmd) ) { die "Failed to remove velvet output. Exiting..."; }

    foreach my $typeV (@typesV) { 
      if( -e $fbaseV . "_" . $typeV . ".sam" ) { 
        $syscmd = sprintf( "rm %s_%s.sam", $fbaseV, $typeV );
        print($syscmd . "\n"); system($syscmd);
      }
      if( -z $fbaseV . "_" . $typeV . ".out" ) { 
        $syscmd = sprintf( "rm %s_%s.out", $fbaseV, $typeV );
        print($syscmd . "\n"); system($syscmd);
      }
      if( -z $fbaseV . "_" . $typeV . ".fa" ) { 
        $syscmd = sprintf( "rm %s_%s.fa", $fbaseV, $typeV );
        print($syscmd . "\n"); system($syscmd);
      }
    }

    $syscmd = sprintf( "rm %s*.fq", $fbaseV );
    print($syscmd . "\n"); system($syscmd);

    $syscmd = sprintf( "rm %s_reads.txt", $fbaseV );
    print($syscmd . "\n"); system($syscmd);

    $syscmd = sprintf( "rm %s.bam", $fbaseV );
    print($syscmd . "\n"); system($syscmd);

    $syscmd = sprintf( "rm %s_%s*.fq", $fbase, $taxtype );
    print($syscmd . "\n"); system($syscmd);

    $syscmd = sprintf( "rm %s_%s.bam", $fbase, $taxtype );
    print($syscmd . "\n"); system($syscmd);

    $syscmd = sprintf( "rm %s_%s_bwa.sam", $fbase, $taxtype );
    print($syscmd . "\n"); system($syscmd);

    $syscmd = sprintf( "rm %s_%s_reads.txt", $fbase, $taxtype );
    print($syscmd . "\n"); system($syscmd);
  }

}

if( !$debug ) {
  $syscmd = sprintf( "rm %s_screen*", $fbase );
  print($syscmd . "\n"); system($syscmd); 

  $syscmd = sprintf( "rm %s_header_chr17.txt", $fbase );
  print($syscmd . "\n"); system($syscmd); 
}
