# Parsing-SAM-and-Extracting-Junctions
this contains the code written for the MSc Bioinformatics course - Programming and Databases for Biologists. grade: A2

this program takes two input files:
1) a SAM file containing all alignments.
2) a tab-separated file. The file has a header and contains three columns. The first is a gene ID. The second is a transcript id. The third is the location of the gene. The location is in the format TGME49_chrVIII:6,793,066..6,795,596(-). This string includes the name of the chromosome where the gene is encoded (TGME49_chrVIII), the start position (6793066), the end position (6795596), and the strand (-).

Usage: python3 myScript.py mySamFile.sam myInputTable.txt

Output:
a tab-delimited file called 2873826.txt that contains a list of all the junctions found. for each junction, the following information is reported - gene ID, start position, end position, and number of reads (alignments) found to map to that junction. the output is reported gene-wise, and there is an empty line after the set of junctions reported for each gene (for better readability).
