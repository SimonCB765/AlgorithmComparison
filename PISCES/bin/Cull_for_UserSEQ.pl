#!/usr/bin/perl -w


#######################MAKE SURE THE INPUT FILE IS DEVOID OF THE .FASTA ON THE END

my $mark = "";
$bindir = 'C:\\Users\\Simonial\\Documents\\PhD\\Culling\\AlgorithmComparison\\PISCES\\bin';

if($#ARGV < 3) {
    &PrintHelp();
    exit;
}

&SetDefault();
&ParseInput();

$file = $list;
$file =~ s/^.*\///;
$name = $list;
$name =~ s/(.*\\)*//;

system($^X, "$bindir\\PSIpdbaaalign.perl", "$file.fasta");

system("move $bindir\\$name.align $bindir\\pdbaa.align");

#system($^X, "$bindir\\Extract_Culled_SEQ.pl", $maxpc, 0, 0, "$file.fasta");

#############################################
sub PrintHelp
#############################################

{
    print "Usage: Cull_for_UserSEQ.pl
           -i input_sequence_file
              It can be a file of sequences in fasta format or output from 
              BLAST/PSI-BLAST running. If PSI-BLAST output is used, the 
              sequences will be taken from the Sbjct: lines from the hits.
           -p maxpc
              percent sequence identity threshold, the valid range is 5-100.
           -l minlen-maxlen (option)
              sequence length range, default is 20-10000
           \n";
}

#############################################
sub SetDefault
#############################################

{
        $minlen = 20;
        $maxlen = 10000;
}

#############################################
sub ParseInput
#############################################

{
        my ($i, @length);

        for $i (0 .. $#ARGV) {
                if($ARGV[$i] eq "-h") {
                        &PrintHelp();
                        exit;
                } elsif($ARGV[$i] eq "-i") {
                        $list = $ARGV[$i+1];
                } elsif($ARGV[$i] eq "-p") {
                        $maxpc = $ARGV[$i+1];
                } elsif($ARGV[$i] eq "-l") {
                        @length = split(/\-/,$ARGV[$i+1]);
                        $minlen = $length[0];
                        $maxlen = $length[1];
                }
        }

        if(!$list || !defined($maxpc)) {
                &PrintHelp();
                exit;
        }

        if($maxpc < 5 || $maxpc > 100) {
                print "Valid \$maxpc = 5 ~ 100\n";
                exit;
        }
}
