Myles Johnson-Gray
#UD CISC 636 - Dr. Li Liao - Homework #1

Both program written in PERL 5

Readme for mutateSequence.pl: (Program used to generate a PAM50 matrix from a PAM1 table.)

-Generates random amino acid sequence from PAM matrix.
-Mutate sequence based on the probabilities from the PAM matrix.
-Create two mutated sequences from the ancestor sequence and find their PAM score.
-Populate and output the PAM50 substitution scoring matrix (to PAM_50.txt).

IMPORTANT THINGS:
-The program opens the "pam1.txt" file so make sure it's in the correct directory.
-The program outputs results to PAM_50.txt

------------------------------------------------------

Readme for GrayM_align.pl:

-Reads in a txt file in FASTA format and extracts both sequences.
-Fill the dynamic programming table using the Needleman-Wunsch algorithm.
-Trace backwards to find the best alignment.
-Output best alignments given scoring scheme.
-If option "-o 1" was given, output the DP table (to align_results.txt).
-Uses Getopt::Std to allow optional command line arguements.

IMPORTANT THINGS:
-Usage statement for input arguments (optional command line arg  "-o 1" , input_file)
-The program takes an input of a text file in FASTA format (make sure this is in the correct directory).
-If the option "-o 1" is given, the DP table is printed to align_results.txt.