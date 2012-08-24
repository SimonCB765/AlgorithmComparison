#!/usr/bin/perl -w

# cullalignblast.perl
# creates blast database from fasta file of sequences
# runs blast of each sequence in the fasta file against the database

$filename = $ARGV[0];

$file = $filename;
$file =~ s/^.*\\//g;    # remove possible path name
$file =~ s/\.fasta//;  # remove possible extension name

# chomp($os = `uname`);

$bindir = 'C:\\Users\\Simonial\\Documents\\PhD\\Culling\\AlgorithmComparison\\PISCES\\bin';

$formatdb = "formatdb.exe";
$blastpgp = "blastpgp.exe";
$blastproc = "$bindir\\psiblastproc.pl";
$tmpdir = "C:\\Users\\Simonial\\Documents\\PhD\\Culling\\AlgorithmComparison\\PISCES\\tmp";
$z = 100000000;

# create blastdatabases and move to $tmpdir; run blast

if(! -e "$tmpdir") {
     system("mkdir $tmpdir") || die "Could not create dir $tmpdir!\n";
}

system("$formatdb -i $filename -p T -n $tmpdir\\$file");
if(! -e "$tmpdir\\$file.phr") {
     print "Failed to create database for blastpgp $file\n";
     exit;
} 

# process blast output to create list of aligned segments

open(ALIGN, "> $file.align") || die "Could not open $file.align for writing:$!";
open(LIST, "$filename") || die "Could not open $filename for reading:$!";

$index = 0;

while(<LIST>) {

   if(/^>/) {

       $head = $_;

       if($index && -e "$index.fasta") {
           close FASTA;
           &DoSearching("$index.fasta");
       }

       $index ++;

       open(FASTA, "> $index.fasta") || die "Could not open $index.fasta for writing: $!";
       print FASTA "$head";
       $flag = 1;
   } elsif($flag == 1) {
       print FASTA "$_";
   }

}

if($index && -e "$index.fasta") {
   close FASTA;
   &DoSearching("$index.fasta");
}

#system("del C:\\Users\\Simonial\\Desktop\\PISCES\\$tmpdir\\$file*.*");

close LIST;
close(ALIGN);

################################################################################
# Carry out PSI-BLAST searching against nr database and use resulted checkpoint
# files to align PDBAA to generate pdbaa.align file
################################################################################

sub DoSearching

{
   my $fastafile = shift;

   system("$blastpgp -i $fastafile -o $fastafile.tmp -e 1 -h 0.0001 -j 3 -N 18 -v 10000 -b 10000 -z $z -d $tmpdir\\$file -F F");# | $blastproc 0.0001 > $fastafile.blast");
   system($^X, "$blastproc", 0.0001, "$bindir\\$fastafile.tmp", "$bindir\\$fastafile.blast");

   unlink "$fastafile";

   if(! -e "$bindir\\$fastafile.blast" || -z "$bindir\\$fastafile.blast") {
        print("No $fastafile.blast!\n");
        return;
   }

   open(BLAST,"$bindir\\$fastafile.blast");

   while(<BLAST>) {
       chop;
       @array=split;

       if (/^Query=/) {
           $query=$array[1];
           $query =~ s/^lcl\|//;
           $qlen=$array[2];
       }

       elsif (/^Hit=/) {
           $hit=$array[1];
           $hit =~ s/^lcl\|//;
           $hlen=$array[2];
       }

       elsif (/^Stat=/) {
           $E=$array[1];
	   $E =~ s/,//g;
           if($E =~ /^e/) { $E = "1" . $E }
           $align=$array[2];
           $pc=$array[3];
           $pc =~ s/\%//g;
       }

       elsif (/^Qseq=/) {
           $qseql=$array[1];
           $qseqr=$array[2];
       }

       elsif (/^Hseq=/) {
           $hseql=$array[1];
           $hseqr=$array[2];

           if ($qlen <= $hlen) {
              printf(ALIGN "%-10s%6d%6d%6d%6d  %-10s%6d%6d%6d%6d%6d%10.1e%6d\n",
                     $query,$qlen,$qseql-1,$qseqr-$qseql+1,$qlen-$qseqr,
                     $hit,  $hlen,$hseql-1,$hseqr-$hseql+1,$hlen-$hseqr,
                     $pc,$E,$align);
           }
       }
    }

    close BLAST;

    unlink "$fastafile.blast";
}
