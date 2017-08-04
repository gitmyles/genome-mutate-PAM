#!/usr/local/perl
# Needleman-Wunsch  Algorithm for Global Alignment of DNA Sequences
#UD CISC 636 - Dr. Li Liao - Homework #1
#Code written by Myles Johnson-Gray (mjgray@udel.edu) 
#----------------------------------------------------------------------------------------------------

use Getopt::Std; # to allow optional command line arguements

#If there exists an "-o" option on cmd line, assign the value to $opt_o.
getopt(o); 
# Usage statement for input arguments (optional command line arg  "-o 1" , input_file)
die "Usage: $0 -o 1 <input file>\n" unless @ARGV == 1;
my $input_file = @ARGV[0]; # store file path

#----------------------------------------------------------------------------------------------------

# SCORING SCHEME (determine methodology for scoring alignment)
my $MATCH    =  4; # +4 for letters that match
my $TRANSISTION = -2; # -2 for a transistion (purine->purine OR pyrimidine->pyrimidine)
my $TRANSVERSION = -3; # -3 for a transversion (purine->pyrimidine OR vice versa)
my $GAP      = -8; # -8 for any gap

#----------------------------------------------------------------------------------------------------

#PROCESS TEXT FILE (parse the input file *in FASTA format* and extract the two sequences)

sub process_file{
    #Initialize sequence strings.
    $seq1 = "";
    $seq2 = "";

    #Open input file.
    open(my $fh, '<:encoding(UTF-8)', $input_file)
      or die "Could not open file '$input_file' $!"; #exit if file can't be found
     
    $seq_count = 0; #keeps track of sequence occurances
     
     #Read each line of the file.
    while (my $row = <$fh>) {
        if (index($row, ">") == -1) { #if there is no ">" on the line, then there is a sequence to read in
            if ($seq_count < 1){ #first sequence is found
                chomp $row;
                $seq1 = $seq1 . $row; #assign sequence to seq1
            }
            else{ #second sequence
                chomp $row;
                $seq2 = $seq2 . $row; #assign sequence to seq2
            }     
            $seq_count += 1;       
        }
    }

    #Print the two sequences.
    print "Input sequences:\n";
    print "---------------------------------------------------------\n";
    print "$seq1\n";
    print "$seq2\n";
}

#----------------------------------------------------------------------------------------------------

#OUTPUT DP TABLE (called using optional cmd line arg "-o 1" ... outputs dynamic programming table)

sub output_DP{
    #Create output file that contains DP table.
    print"Printing DP table to \"align_results.txt\"\n";
    my $filename = 'align_results.txt';  #name of output file
    open(my $fh, '>', $filename); 
    print $fh "Dynamic Programming Table: <row, column> pair and resulting score. \n\n";

    #Iterate through table and display score and pointer for each row/col combination.
    for($i = 0; $i <= length($seq2); $i++) {                        
        for($j = 0; $j <= length($seq1); $j++) {                   
            $char1 =" "; #initialize char1 string
            $char2 =" "; #initialize char2 string
            
            if ($i > 0){ #first row of table isn't occupied by nuclueotide character
                $char2= substr($seq2, $i-1, 1); #assign the current character in seq2
            }
            if ($j > 0){ #first column of table isn't occupied by nucluoetide character
                $char1 = substr($seq1, $j-1, 1); #assign the current character in seq1           
            }
            
            #Output <row,column> pair with resulting score.
            #print "row($i):$char2 -> col($j):$char1     ($matrix[$i][$j]{score})\n"; #this will print to the console output (for debugging)
            print $fh "row($i):$char2 -> col($j):$char1     ($matrix[$i][$j]{score})\n";
        }
    }
    close $fh; #close output file
}

#----------------------------------------------------------------------------------------------------

#*******************MAIN  EXECUTION*******************

process_file(); #read input file into two sequences

#INITIALIZATION (initialize and fill the dynamic programming table)

@matrix; #initialize hash/table

#Prepare the first column and row of the table.
$matrix[0][0]{score}   = 0;
$matrix[0][0]{pointer} = "none";
for(my $j = 1; $j <= length($seq1); $j++) {
    $matrix[0][$j]{score}   = $GAP * $j;
    $matrix[0][$j]{pointer} = "left";
}
for (my $i = 1; $i <= length($seq2); $i++) {
    $matrix[$i][0]{score}   = $GAP * $i;
    $matrix[$i][0]{pointer} = "up";
}

# Fill the inside of the dynamic programming table.
for(my $i = 1; $i <= length($seq2); $i++) {
    for(my $j = 1; $j <= length($seq1); $j++) {
        my ($diagonal_score, $left_score, $up_score);

        #Retrive the current characters of the input sequences.
        my $letter1 = substr($seq1, $j-1, 1);
        my $letter2 = substr($seq2, $i-1, 1);                         
   
        #Calculate the scores for a match or mismatch.
        if ($letter1 eq $letter2) { #if characters are aligned...
            $diagonal_score = $matrix[$i-1][$j-1]{score} + $MATCH;
        }
        #If there is a mismatch....(transistion vs. transversion)
        elsif(($letter1 eq "A" && $letter2 eq "G") || ($letter1 eq "G" && $letter2 eq "A") || ($letter1 eq "C" && $letter2 eq "T") || ($letter1 eq "T" && $letter2 eq "C")) { #if transistion..
            $diagonal_score = $matrix[$i-1][$j-1]{score} + $TRANSISTION; 
        }
        elsif(($letter1 eq "A" && $letter2 eq "C") || ($letter1 eq "C" && $letter2 eq "A") || ($letter1 eq "A" && $letter2 eq "T") || ($letter1 eq "T" && $letter2 eq "A")
        || ($letter1 eq "G" && $letter2 eq "C") || ($letter1 eq "C" && $letter2 eq "G") || ($letter1 eq "G" && $letter2 eq "T") || ($letter1 eq "T" && $letter2 eq "G")) { #if transversion..
            $diagonal_score = $matrix[$i-1][$j-1]{score} + $TRANSVERSION; 
        }
  
        #Calculate the scores for a gap.
        $up_score   = $matrix[$i-1][$j]{score} + $GAP;
        $left_score = $matrix[$i][$j-1]{score} + $GAP;

        #Choose the best possible score.
        if ($diagonal_score >= $up_score) {
            if ($diagonal_score >= $left_score) {
                $matrix[$i][$j]{score}   = $diagonal_score; #assign diagonal score
                $matrix[$i][$j]{pointer} = "diagonal";
            }
            else {
                $matrix[$i][$j]{score}   = $left_score; #assign left score
                $matrix[$i][$j]{pointer} = "left";
            }
        } 
        else {
            if ($up_score >= $left_score) {
                $matrix[$i][$j]{score}   = $up_score; #assign up score
                $matrix[$i][$j]{pointer} = "up";
            }
            else {
                $matrix[$i][$j]{score}   = $left_score; #assign left score
                $matrix[$i][$j]{pointer} = "left";
            }
        }
    }
}

#----------------------------------------------------------------------------------------------------

# TRACEBACK (retrace backwards to find the best alignment)

#Initilaize  alignment strings.
my $align1 = "";
my $align2 = "";

#Start at the last cell of the matrix. 
my $j = length($seq1);
my $i = length($seq2);

while (1) {
    last if $matrix[$i][$j]{pointer} eq "none"; # end if you reach the first cell of the matrix
    
        if ($matrix[$i][$j]{pointer} eq "diagonal") { #if diagonal pointer...
            $align1 .= substr($seq1, $j-1, 1);
            $align2 .= substr($seq2, $i-1, 1);
            $i--;
            $j--;
        }
        elsif ($matrix[$i][$j]{pointer} eq "left") { #elif left pointer...
            $align1 .= substr($seq1, $j-1, 1);
            $align2 .= "-";
            $j--;
        }   
        elsif ($matrix[$i][$j]{pointer} eq "up") { #elif up pointer...
            $align1 .= "-";
            $align2 .= substr($seq2, $i-1, 1);
            $i--;
        }    
}

#Reverse and output the alignments.
$align1 = reverse $align1;
$align2 = reverse $align2;
print "\nReporting best alignment given scoring scheme: \n";
print "---------------------------------------------------------\n";
print "$align1\n\n";
print "$align2\n\n";

#If option "-o 1" is given, output the DP table to align_results.txt
if($opt_o == 1){
    output_DP();
}
print "---------------------------------------------------------\n";