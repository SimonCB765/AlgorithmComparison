#!/usr/bin/perl -w

# make_filter_refined.pl   06-26-2001
#

   $E_cutoff = $ARGV[0];
   $infile = $ARGV[1];
   $output = $ARGV[2];
   %pair = ();

   # go through all created results for all rounds 
       
   $hit_flag = 0;
   $pair_flag = 0;
   $repeat = 0;
   $round = 0;
   $lastround = 0;
   %hits = ();
   %hits0 = ();

   # read information from each round's scop database searching output.

   $eof_mark = 0;
   
   open(TOREAD, "$infile");

READ:
   while(<TOREAD>) {
     chomp;
     @line = split;

     if(/^Query=/) {
	 $query = $line[1];

     } elsif(/^Results from round/) {
	 $round = $line[3];

     } elsif(/^Searching................./) {
         $lastround = $round;

     } elsif(/^  Database: /) {
	 $eof_mark = 1;
         $lastround = $round;
         last;

     } elsif(/ letters\)$/) {
         $qlength = $line[0]; $qlength =~ s/^\(//;

     } elsif(/^Sequences producing significant alignments:/) {
         $hit_flag = 1;
         $pair_flag = 0;

     } elsif(/^ \*\*\*\*\* No hits found \*\*\*\*\*\*/) {
   	 last;

     # there may be some repeats for specific hit, some of them are
     # correct, just come from different segment matchs; some of them 
     # indicate some type of error

     } elsif($hit_flag == 1 && /^>/) {
          # Check drift

          if($pair_flag == 0) {
              foreach $hit1 (keys %hits0) {
                 if($hits0{$hit1} > $E_cutoff) { next }
                 if(!$hits{$hit1}) {
                     $lastround = $round - 1;
                     # print "Drift in $round for $hit1\n";
                     last READ;
                 }
              }         

              %hits0 = %hits;
              %hits = ();
          }
           
	  $hit = $line[0]; $hit =~ s/^>//;
	  $pair_flag = 1;
	  $repeat = 0;


     } elsif($hit_flag == 1 && $pair_flag == 0) {
          if(!$line[0] || $line[0] =~ /Sequence/) { next }
          $hit1 = $line[0];
          $evalue1 = $line[-1];
          if($evalue1 =~ /^e/) { $evalue1 = "1" . $evalue1 }
          $hits{$hit1} = $evalue1;

     } elsif($pair_flag == 1 && (! /^\>/)) {
          if(/^\s+Length =/) {
      	     $hlength = $line[2];
	     $pair{$round}{$hit}{hlength} = $hlength;
       
          } elsif(/^ Score/) {
	     $repeat++;
             if($repeat > 1) { next }

	     $Score = $line[2]; 
	     $raw_score = $line[4]; 
	     $raw_score =~ s/\(//; $raw_score =~ s/\).*$//;
	     $evalue = $line[7];
	     if($evalue =~ /^e/) {
		$evalue = "1" . $evalue;
	     }

	     $pair{$round}{$hit}{Score} = $Score;
	     $pair{$round}{$hit}{raw_score} = $raw_score;
	     $pair{$round}{$hit}{evalue} = $evalue;

          } elsif(/^ Identities/) {
             if($repeat > 1) { next }

             $length = $line[2];
             $length =~ s/^.*\///;
	     $Identities = $line[3];
	     $Positives = $line[7];
	     $Identities =~ s/\(//; $Identities =~ s/\%.*$//;
	     $Positives =~ s/\(//;  $Positives =~ s/\%.*$//;

	     $pair{$round}{$hit}{Identities} = $Identities;
	     $pair{$round}{$hit}{Positives} = $Positives;
             $pair{$round}{$hit}{Len} = $length;

	     if($#line > 7) {
		$Gaps = $line[11];
		$Gaps =~ s/\(//; $Gaps =~ s/\%.*$//;
		$pair{$round}{$hit}{Gaps} = $Gaps;
	     } else {
		$pair{$round}{$hit}{Gaps} = 0;
	     }

          } elsif(/^Query:/) {
             if($repeat > 1) { next }
	     if(!defined($pair{$round}{$hit}{qstart})) {
	         $pair{$round}{$hit}{qstart} = $line[1];
	     }
             $pair{$round}{$hit}{qend} = $line[3];
	     $pair{$round}{$hit}{qbody} .= $line[2];

          } elsif(/^Sbjct:/) {
             if($repeat > 1) { next }
	     if(!defined($pair{$round}{$hit}{hstart})) {
	         $pair{$round}{$hit}{hstart} = $line[1];
	     }
             $pair{$round}{$hit}{hend} = $line[3];
	     $pair{$round}{$hit}{hbody} .= $line[2];

          } else {
	     next;

          }
      } 
	       
   }
   close(TOREAD);

   # print "lastround $lastround\n";
   open(OUTSTUFF, "> $output");

   printf(OUTSTUFF "Query= $query $qlength $query\n");
   foreach $hit (keys %{$pair{$lastround}}) {
      printf(OUTSTUFF "Hit= $hit $pair{$lastround}{$hit}{hlength} $hit\n");
      printf(OUTSTUFF "Species= NA\n");
      printf(OUTSTUFF "Stat= $pair{$lastround}{$hit}{evalue} $pair{$lastround}{$hit}{Len} $pair{$lastround}{$hit}{Identities} $pair{$lastround}{$hit}{Positives} $pair{$lastround}{$hit}{Gaps} \n");
      printf(OUTSTUFF "Qseq= $pair{$lastround}{$hit}{qstart} $pair{$lastround}{$hit}{qend} $pair{$lastround}{$hit}{qbody}\n");
      printf(OUTSTUFF "Hseq= $pair{$lastround}{$hit}{hstart} $pair{$lastround}{$hit}{hend} $pair{$lastround}{$hit}{hbody}\n");
      printf(OUTSTUFF "\n");
   }
   close(OUTSTUFF);
