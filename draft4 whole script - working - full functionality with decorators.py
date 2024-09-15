#here we go. 
#we import the following modules for use through the code. they're described in a comment whenever used.
import re
import sys
import logging

#creating the command line
#this uses the inbuilt python module called sys, and a special object called argv. every argument in the command line is noted by sys and fed into the python script. this is zero indexed and the zero-th object is always the name of the python script.
sam_file = sys.argv[1]
txt_file = sys.argv[2]
#tested, working.

#creating a logger for error handling. this uses an inbuilt module called logging.
#getting the logger
logger = logging.getLogger()
#setting the level of the logger
logger.setLevel(logging.WARNING)
#creating a stream handler.
shandler = logging.StreamHandler()
#formatting the stream handler
shandler.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - {%(pathname)s:%(lineno)d} - %(message)s'))
#adding the stream handler to the logger.
logger.addHandler(shandler)

#creating a dictionary to store information about each junction. more details in comments downstream.
junction_dict = {}

#opening the filehandle.
try: #catching error in case sam file is not found.
    with open(sam_file) as in_file:

        #creating an iterator for each line.
        for line in in_file:

            #skipping the headers using an if statement. since all headers start with the symbol '@', we can use this symbol to filter them out.
            if line.startswith('@'):
                next(in_file)

            else:
                #if a line doesn't start with '@', then the header section is over and the data instances have started.

                #we split the line by the tab delimit and strip it of the 'next line' character at the end.
                line = line.rstrip().split('\t')

                #making variables for the column that we are interested in.
                #according to the assignment document, we want to consider columns - 3 (name of chromosome or RNAME), 4 (position of chromosome where the alignment starts or POS), 6 (a string describing the format - the CIGAR string), and the last column (an integer describing the number of times the given read aligned to the reference genome).
                #the given column indexes are in 'human' notation, meaning they are 1-indexed. to convert them, we subtract 1 from them, since column 2 in 'human' notation is column 1 in python notation (that is 2-1). similarly, column 1 in 'human' notation is column 0 in python notation, and so on.
                try: #catching error if the indexing is not the correct sam file indexing.
                    
                    chrom_aligned_to = line[2]
                    
                    POS = line[3]
                    try: #catching error is the POS value is wrong.
                        POS = int(POS) #we need the start position to be an integer so that mathematical calculations can be done downstream.
                    except ValueError:
                        logger.warning(f"Hey user, this is from script 2873826. The POS or position of chromosome where the alignment starts is not an integer for the query template name {line[0]}. Skipping this. Cheers!\n")
                        continue

                    CIGAR = line[5]
                    read_align_count = line[-1]

                    #print(f"{chrom_aligned_to} {POS} {CIGAR} {read_align_count}")
                    #break
                    #we have checked that this has worked as intended.

                    #since we are only interested in reads which have uniquely aligned to the reference genome, i.e. only aligned once, we use an if statement to filter these out.
                    if read_align_count.endswith('1'): #the function 'endswith' only takes strings and tuples as arguments, not integers.
                        
                        #print(read_align_count)
                        #break
                        #we have checked that this has worked as intended.

                        #we check if the read is 'split', meaning if it contains a region that has 'skips' or Ns when aligned to the reference.
                        #info - a split region in a read is called a 'junction'.
                        #to do this, the number of Ns in the CIGAR string is counted. for split reads, the count should at least be 1.
                        N_count = CIGAR.count('N')

                        #we use an if statement to filter out split reads.
                        if N_count >= 1:
                            
                            #print(line)
                            #break
                            #we have checked that this has worked as intended.

                            '''
                            how to parse the CIGAR string?
                            facts to consider - 
                            1) these strings are number:letter pairs. the letter can be M for Match, I for Insertion, D for Deletion, N for skipped region, and S for Soft Clip. the number can be any natural number (i.e. an integer, not including zero or negative numbers, and not including floats/decimal numbers). 
                            2) skipped regions cannot lie towards the edges of the read, they only lie within the read. i.e. there is always a number:letter pair to the left AND right of a number:'N' pair.
                            3) one string can contain multiple junctions.

                            the method - 
                            1) a range is created and iterated over using the N_count to extract substrings using regex. we want to extract as many substrings as there are Ns in the read.
                            2) for every number in the range, a substring is extracted from the original CIGAR string. Since we're counting the start positions and end positions of the junction, we can extract a substring that each ends at N.
                                e.g. if our CIGAR string is 5M5N6M6N7M, then our N_count will be 2, and we will extract two substrings - 5M5N and 5M5N6M6N. 
                            3) for these substrings, we extract the number:letter pairs and iterate over them to calculate the junction start and end positions.
                            4) storing these calculations is commented upon in a later comment block.
                            '''
                            
                            for number in range(1,N_count+1,): #the range start_value is 1. if only N_count was chosen as the end_value, then the last number in the range would be the just short of N_count, as python has half open notation. the default step value is 1. the end_value is modified to N_count + 1 so that with half open notation, python stops just short of (N_count + 1), which is N_count. this will return a range of numbers counting up to N_count, including N_count.

                                '''
                                to find the substring, regex is used.
                                the regex is - (([\d]+[MIDS]{1})+([\d]+N)){number}
                                let's break this apart.
                                1) ([\d]+[MIDS]{1})+
                                    any digit occuring one of more times, followed by a letter which can be M or I or D or S, repeated exactly one time.
                                    this entire pattern may repeat one of more times, until - 
                                2) ([\d]+N)
                                    any digit occuring one or more times, followed by N occuring once.
                                    this gives us back a pattern containing number(s) followed by one letter followed by number(s) followed by N.
                                3) (([\d]+[MIDS]{1})+([\d]+N)){number}
                                    putting the above pattern together, we want the entire pattern to repeat 'number' number of times. (number is the iterable variable of the range we created above.)

                                let's go back to our first example. 
                                CIGAR = '5M5N6M6N7M'
                                our range will be (1,2).
                                when the loop runs the first time, number = 1. the regex will be (([\d]+[MIDS]{1})+([\d]+N)){1}. this will return the match 5M5N.
                                when the loop runs the second time, number = 2. the regex will be (([\d]+[MIDS]{1})+([\d]+N)){2}. this will return the match 5M5N6M6N.
                                ...and so on.
                                '''

                                pattern = "(([\d]+[MIDHS]{1})+([\d]+N))" + "{" + f"{number}" + "}" #through testing, I found that I cannot directly use the 'number' variable inside a search string, so a pattern variable is built.
                                match = re.search(pattern,CIGAR)
                                #match contains the substring of the CIGAR string.
                                #so, we extract information about one junction per match.
                                if match:

                                    #number:letter pairs are extracted for each substring.
                                    pairs = re.finditer(r'([\d]+)([MNIDS]{1})',match.group())

                                    #now we calculate junction end. we just add the number for each number:letter pair to a variable.
                                    '''
                                    important note - for these calculations - I and S are ignored because they don't exist in the reference genome, and we are interested in the positions of the junctions on the reference genome.
                                    '''
                                    junction_end = 0+POS
                                    for pair in pairs: 
                                        if pair.group(2) in ['M','N','D']:
                                            junction_end=junction_end+int(pair.group(1))
                                            #this loop goes through every number:letter pair and if the letter is an M, N or D, adds the corresponding number to a variable.
                                    
                                    #we do the same for junction start. since re.finditer returns a one-shot iterator, we have to do the search again.
                                    pairs = re.finditer(r'([\d]+)([MNIDS]{1})',match.group())

                                    junction_start = 0+POS
                                    for pair in list(pairs)[:-1]: #this leaves out the last number:letter pair from the calculation; since we're calculating the start of the junction, we don't want to add in the number corresponding to the last N.
                                        if pair.group(2) in ['M','N','D']:
                                            integer = int(pair.group(1))
                                            junction_start = junction_start+integer
                                            #print(f"flag for me {junction_start}")
                                            #this is working as intended.
                                    
                                    #print(match.group())
                                    #print(POS)
                                    #print(junction_start)
                                    #print(junction_end)
                                    #we have checked that this is working as intended.
                                    '''
                                    the dictionary created at the beginning will store all the values. 
                                    we note that neither our chromosome, nor our junction_start, nor our junction_end are individually unique, it is just the combination of these values that is unique for each junction. so, for the key, we will create a string concatenating these values, as the combination is unique. this will ensure that the key remains unique.
                                    in the value for each key in the dictionary, we will store a list, like so - 
                                    ['chromosome','junction_start','junction_end','read_count']
                                    the read count is the number of reads that 'cover' the same region in the genome. so if a particular junction is covered by many reads, it has high coverage. we count this coverage as the number of reads.
                                    for each substring, we calculate the values for the junction and look them up in the dictionary. if the junction doesn't exist, we add it to the dictionary; if it does exist, then we add '1' to the read count. 
                                    '''
                                                                
                                    key = f"{chrom_aligned_to}:{junction_start}:{junction_end}"
                                    
                                    #print(key)  
                                    #this is working as intended.

                                    if key in junction_dict:
                                        list_of_values = junction_dict[key]
                                        list_of_values[3]=int(list_of_values[3])+1
                                    else:
                                        junction_dict[key]=[chrom_aligned_to,junction_start,junction_end,1] #if the junction doesn't already exist in the dictionary, then it means that our match substring is the first read to have found the junction, so the number of reads supporting the junction is 1.                           
                except IndexError:
                    logger.error(f"Hey user, this is script 2873826. The SAM file does not have the expected indexing format. Please rerun the script with the correct SAM file. Cheers!\n")
except FileNotFoundError:
    logger.error(f"Hey user, this is script 2873826. The SAM file was not found. Please rerun the script and type the file name correctly, with the address if it is not in the working directory. Cheers!\n")

#print(junction_dict)       
#print(len(junction_dict.keys())) #227
#this is working as intended.

#now that the sam file has been parsed, we move on to writing the output.
#the output file name is hardcoded. because we will refer to the gene location file and side-by-side write the output, we open two filehandles with 'with' statements as follows.
with open('2873826.txt','w') as out_file:
    try: #catching error if the txt file is not found.
        with open (txt_file) as in_file:
            next(in_file)
            for line in in_file:
                geneID, sourceID, location = line.rstrip().split('\t')
                #the location tab is in the following format - (explained using example).
                #TGME49_chrVIII:6,631,349..6,636,865(+)
                #here the first bit (until the colon) is the chromosome ID, the number until the two dots is the gene start position, the number after the two dots is the gene end position, and the bracketed sign in the end indicates the strand. here, (+) means plus strand. for the purposes of this assessment, the strand is not taken into consideration.
                #the location can be parsed using regex.
                
                #print(location)
                #this is working as intended.

                match = re.search(r'([\w]+)(:)(.+)(\.\.)(.+)(\([+-]\))',str(location))
                #the regex here is '([\w]+)(:)(.+)(\.\.)(.+)(\([+-]\))'. here ([\w]+) finds the chromosome which is composed of one or more word characters. (:) finds the semi colon. (.+) finds the first number. here, [\d] is not used because the number contains commas. (\.\.) finds the two dots. (.+) finds the second number. (\([+-]\)) finds either '(+)' or '(-)'.
                if match:
                    #print(match.groups())
                    #break
                    #the regex works as intended.
                    
                    chromosome = match.group(1)
                    #print(chromosome)

                    gene_start = match.group(3)
                    #print(gene_end)
                    #break     
                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
                    gene_end = match.group(5)
                    #print(gene_end)
                    #break

                    #for both gene_start and gene_end, we need to remove the commas and convert them into integers.
                    gene_start = int(gene_start.replace(',',''))
                    #print(gene_start)
                    #print(type(gene_start))
                    #break
                    #this works as intended.

                    gene_end = int(gene_end.replace(',',''))

                    #iterating through the dictionary
                    for value in junction_dict.values():
                        #each value is a list in the format [chromosome_aligned_to, junction_start, junction_end, read_count].
                        if chromosome == value[0] and value[1] >= gene_start and value[2] <= gene_end: #the statement checks if the value lies within the bounds of the given gene.
                            out_file.write(f"{geneID}\t{value[1]}\t{value[2]}\t{value[3]}\n")

                    #if the given gene has even a single junction, adding a new line after the junction(s) has/have been written. if gene has no junction, new line will not be added. if gene has multiple junctions, single line will be added as the loop breaks after the first execution.
                    for value in junction_dict.values():
                        if chromosome == value[0] and value[1] >= gene_start and value[2] <= gene_end:
                            out_file.write(f"\n")
                            break     
                else: #catching error if the location is written in an unexpected/inconsistent format.
                    logger.error(f'Hey user, this is script 2873826. The genomic location in the text file for gene {geneID} is in an unexpected format. Skipping this. Cheers!\n')                                                            
    except FileNotFoundError:
        logger.error(f"Hey user, this is script 2873826. The text file containing genomic locations was not found. Please rerun the script and type the file name correctly, with the address if it is not in the working directory. Cheers!\n")                