#!/usr/local/perl
#Genome Mutation using PAM matrix
#UD CISC 636 - Dr. Li Liao - Homework #1
#Code written by Myles Johnson-Gray (mjgray@udel.edu) 
#----------------------------------------------------------------------------------------------------

use POSIX; #to use ceil(), floor()
#INITIALIZATION (open the PAM matrix and generate a random amin acid sequence)
#contributions from Dr. Li Liao (lilia0@udel.edu)

#open supplied PAM matrix
open(F, "C:/Users/Myles/Documents/School/UD/CompBio/pam1.txt");

$line = <F>; # read first line
$line =~ s/^\s+//; # remove beginning spaces.
@aa = split(/\s+/,$line); # this array has the twenty amino acids

$i=0;
while($line = <F>) {
  @col = split(/\s+/,$line);
  shift(@col);
  for($j=0; $j<=$#col; $j++) {
    $s{$aa[$i]}{$aa[$j]} = $col[$j];  # score matrix indexed by amino acid single letter code
  }
  $i++;
}
close(F); #close PAM matrix

#print "$s{A}{V}\n";  # you can access the score for a given amino acid pair
print "Initial \"ancestor\" sequence below:"; #print inital ancestor sequence
print "\n---------------------------------------------------------\n";

$i=0;
while($i<500) {
  $r = int(rand()*20); #random integer between 1 and 20
  $seq .= $aa[$r];  # concatenate aa to seq
  $i++;
}

print $seq;  # this is a random sequence of 500 aa
print "\n---------------------------------------------------------\n";

$ancestor = $seq;

$array_length = scalar @aa; #length of amino acid array (20 different amino acids)

#----------------------------------------------------------------------------------------------------

#MUTATE SEQUENCE (mutate a given amino acid sequence using a PAM matrix)

sub mutate_seq
{
  #Function inputs & variable initialization
  my $mutated_seq = $_[0];#input amino acid sequence
  my $num_loops = $_[1]; #number of times to attempt mutation
  $mut_count = 0; #count of the total mutations occured
  
  print "Preparing to mutate sequence...\n";
  
  #Attempt to preform mutation the specified number of times.
  for($h=0; $h<$num_loops; $h++){
    
    #Iterate over the input amino acid sequence.
    for($i=0; $i<length($mutated_seq); $i++){
      $mutation = int rand(10100); #generate random number to guide mutation
      $char = substr($mutated_seq, $i, 1); #holds the current amino acid character
      $score = 0; #holds current PAM matrix score
      
      #Add up scores from PAM matrix to determine the mutation.
      for($j=0; $j < $array_length; $j++){        
        $score = $s{$char}{@aa[$j]} + $score; #add amino acid pair score to current score
        #print "[$char][@aa[$j]] $score - $mutation \n"; #(For debugging. Prints aa pair, current score, and random number.) 
        
        #The outcome (mutation or not...) has been chosen.
        if($score >= $mutation){
          if($char ne @aa[$j]){ #if a mutation (amino acid change) has been chosen...
            substr($mutated_seq, $i, 1) = @aa[$j]; #mutate amino acid sequence at current spot
            $mut_count+=1; #increment mutation counter
          }
          last; #if some outcome has been chosen, break out of this loop
        }      
      }
    }
  }
  
  #Format and output information to the user.
  print "Mutation complete (mutated sequence below). \n";
  print "Total number of mutations: $mut_count";
  print "\n---------------------------------------------------------\n";  
  print $mutated_seq;
  print "\n---------------------------------------------------------\n";  
  return($mutated_seq); 
}

#----------------------------------------------------------------------------------------------------

#LOGARITHM FUNCTION (returns the log of a number with the specified base)
sub log_N
{
  my $num = shift; #number to be log'd
  my $base = shift; # base for log
  return log($num)/log($base);
}

#----------------------------------------------------------------------------------------------------
#MAIN EXECUTION: (Generate two mutants and populate the subsitution scoring matrix)

#Create two mutations from the ancestor sequence.
$mutant1 = mutate_seq($ancestor, 50);
$mutant2 = mutate_seq($ancestor, 50);
#Traverse all pairs of amino acid pairs.
for($i=0; $i<$array_length; $i++){
  for($j=0; $j<$array_length; $j++){
    #Initialize variables
    $aa_count1 = 0; #occurances of amino acid in mutant1
    $aa_count2 = 0; #occurances of amino acid in mutant2
    $align_count = 0; #number of times the amino acid pairs are aligned
    $aa_score = 0; # holds the score for the current amino acid pair

    #Iterate through both mutant sequences.
    for($k=0; $k<length($mutant1); $k++){
      
      if(substr($mutant1, $k, 1) eq @aa[$i]){ #if amino acid #1 is found in mutant1
        $aa_count1+=1; #increment counter
      }
      if(substr($mutant2, $k, 1) eq @aa[$j]){ #if amino acid #2 is found in mutant2
        $aa_count2+=1; #increment counter
      }
      
      #If the sequences are aligned...
      if((substr($mutant1, $k, 1) eq @aa[$i] && substr($mutant2, $k, 1) eq @aa[$j]) || 
      (substr($mutant1, $k, 1) eq @aa[$j] && substr($mutant2, $k, 1) eq @aa[$i])){  
        $align_count+=1; #increment the alignment calculator
      }     
    }
    
    #Calculate the PAM50 score for the amino acid pair: S(i,j) = log[f(i,j) / f(i)f(j))]
    $aa_score = ($align_count / 500) / (($aa_count1 / 500)*($aa_count2 / 500 )); 
    #$aa_score = (($align_count+.1) / (500+.1)) / (($aa_count1+.1) / (500+.1))*(($aa_count2+.1)/ (500+.1));  #PSUEDOCOUNT??
    
    if($aa_score != 0){ #check to protect against zero division
      $aa_score = log_N($aa_score, 2); #take log(base2) of the score
    }
    #Round the scores to integer values(COMMENT OUT IF YOU WANT).
    # if($number < 0.0){
    #   $aa_score = floor($aa_score); #round down
    # }
    # else{
    #   $aa_score = ceil($aa_score); #round up
    # }
        #Populate scoring matrix with amino acid pair score
    $score_matrix{@aa[$i]}{@aa[$j]} = $aa_score;
  }
}

#Display the substiution score matrix (
print"Printing scoring matrix to \"PAM50_results.txt\"\n";
my $filename = 'PAM50_results.txt';  #name of output file
open(my $fh, '>', $filename); 
print $fh "PAM 50 Substitution Scoring Matrix: <aa,aa> pair and resulting score. \n\n";

#Output <aa,aa> pair with resulting score.
for($i = 0; $i < $array_length; $i++) {                        
    for($j = 0; $j < $array_length; $j++) {                   
        #print "<@aa[$i]@aa[$j]> -> $score_matrix{@aa[$i]}{@aa[$j]}\n"; #outputs results to console...good for debugging
        print $fh "<@aa[$i]@aa[$j]> -> $score_matrix{@aa[$i]}{@aa[$j]}\n";
    }
}
close $fh; #close output file
print"\n";

# print "$score_matrix{A}{A}\n"; #You can acess the substitution score matrix in this way.