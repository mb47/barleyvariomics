import sys

#====================
# Read in variables
#====================
file=sys.argv[1]
transcriptName=sys.argv[2]
accession=sys.argv[3]
strand=sys.argv[4]

# Declare an empty dictionary
d={}
# Start a counter 
# to be able to assign number to sequences presented in the input file
# to keep track of the order of the sequences
count=0

#====================
# Read and parse the fasta file
# note: use with, so no need to close the file later
# 'r': reading mode
#====================
with open(file,'r') as f:
    for line in f:
        if line[0]=='>': #is header
            count+=1
            d[count]={} #create an empty dictionary, using the count as the key
            header=line.rstrip() #format correction: keeping the sequence header, and remove the new line symbol '\n' at the end of the header
            # if want to remove '>' as well, do: header=line[1:].rstrip()
            seq='' #create an empty string to store the sequence of the corresponding header when it spans multiple lines in fasta
        else: #is sequence
            seq+=line.rstrip() #append sequence to the string
            d[count][header]=seq #update dictionary to look like d={1:{head1:seq}}

# check data structure    
#print(d)

#====================
# writing output
#====================
# create an empty list to store all the reverse complementary sequences IN ORDER
# note: because dictionary does not have order (you can't "reverse" a dictionary)
all_rev_comp=[]
# write new fasta sequence name to stdout
print('>'+transcriptName+'_'+accession)
# iterate though keys and values in the dictionary
for k, v in d.items(): #iterate through keys (i.e., count) and values (i.e., nested dictionary {header:seq} pairs)
    for header, seq in v.items(): #iterate through the nested dictionary, for easy access to its values
        if strand=='+':
            print(seq)
        else:
            rev=seq[::-1] #print from last digit of the seq
            rev_comp=''
            for nuc in rev:
                if nuc=='A':
                    rev_comp+='T'
                elif nuc=='G':
                    rev_comp+='C'
                elif nuc=='C':
                    rev_comp+='G'
                elif nuc=='T':
                    rev_comp+='A'
            all_rev_comp.append(rev_comp)

all_rev_comp.reverse()
for i in all_rev_comp:
    print(i)
