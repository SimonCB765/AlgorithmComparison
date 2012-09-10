#!/usr/bin/perl -w


# MAKE SURE THE INPUT FILE IS DEVOID OF THE .FASTA ON THE END

$bindir = 'C:\\Users\\Simonial\\Documents\\PhD\\Culling\\AlgorithmComparison\\PISCES\\bin';

&ParseInput();

$file = $list;
$file =~ s/^.*\///;
$name = $list;
$name =~ s/(.*\\)*//;

system($^X, "$bindir\\PSIpdbaaalign.perl", "$file.fasta");

system("move $bindir\\$name.align $bindir\\pdbaa.align");

#############################################
sub PrintHelp
#############################################

{
    print "Usage: Cull_for_UserSEQ.pl
           -i input_sequence_file
              A file in FASTA format is expected (with the extension removed e.g. path\to\file\Human not path\to\file\Human.fasta).
           \n";
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
                }
        }

        if(!$list{
                &PrintHelp();
                exit;
        }
}
