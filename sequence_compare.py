#!/usr/bin/env python

import sys
import traceback
from Bio import SeqIO
import csv
import math
from Levenshtein import *

##Class definitions##
#This is to define the FASTAClass
class FASTAClass:
    def __init__(self, id, seq):        
        self.id = id        
        self.seq = seq

#This is to define the CompareClass
class CompareClass:
    def __init__(self, first, second, value):        
        self.first = int(first) 
        self.second = int(second)
        self.value = float(value)

##Function definitions##
#Function to split sequences.
def seqSplitter(input_file,coreNumber):
    record = {}
    #This is where the new file is read.
    sequence = open(input_file, "rU")
    record = list(SeqIO.parse(sequence, "fasta"))
    sequence.close()
    
    #List all the comparissons that need to be done based on the file index of the sequence file.
    count = 0
    compareList = [[[] for i in range(2)] for i in range(int(len(record)*(len(record)-1)/2))]
    i=0
    j=0
    for i in xrange(len(record)):
        for j in xrange(i):
            compareList[count][0] = i
            compareList[count][1] = j
            count+=1

    print "We need to do ", count, " comparissons! OR ", int(len(record)*(len(record)-1)/2)
    length = float(len(record))
    #print length*(length-1)/(2*coreNumber)
    #exit()

    #Writing all the compareList files that each core needs to handle. This is based on the number of sequences and the number of cores we are going to use.
    count=0
    writeCount=0
    for writeCount in xrange(int(coreNumber)):
        writeIndex = "Temp/"+str(writeCount+1)+".dat"
        output_handle = open(writeIndex, "wb")
        output_csv = csv.writer(output_handle, delimiter='\t')
        countI=0
        for countI in xrange(int(math.ceil(length*(length-1)/(2*coreNumber)))):
            if (count < int(length*(length-1)/2)):
                output_csv.writerow([compareList[count][0]] + [compareList[count][1]])
                #output_handle.write(compareList[count][0]+"\t"+compareList[count][1]+"\n")
                count+=1 
        output_handle.close()
        print "Writing ", writeCount+1, " file of ", int(coreNumber)

    print "Done splitting the comparisson list."
 
def seqCompare(fasta_file,line_file,coreNumber):
    #Create our hash table to add the sequences and ID table to note replicates
    sequences={}
    #Using the biopython fasta parse we can read our fasta input
    records = open(fasta_file, "rU")
    sequences = list(SeqIO.parse(records, "fasta"))
    records.close()
   
    #Create compare table of sequence indices and the Normalized Levenshtein distance
    length=float(len(sequences))
    compareList = {} #[[[float(-1)] for i in range(3)] for i in range(int(math.ceil(length*(length-1)/(2*coreNumber))))]
    #print "The size of the compare list is:\t", len(compareList), "\t", len(compareList[0])
    #Upload the comparisson indices
    line = open(line_file, "rb")
    csvread = csv.reader(line, delimiter='\t')
    #line.close()
    #print csvread.next()
    count=0 
    for row in csvread:
        #print row
        compareList[count] = CompareClass(row[0],row[1],0)
        #compareList[count][1] = row[1]
        count+=1
    line.close()
    #print "The csv file gives the first row as:\n", type(temp), "\n", compareList[0][0], "\t", compareList[0][1]
    #exit() 
    
    #Run the comparison of the sequences with the indices listed in compareList. This is going to be parallelised by splitting the list among cores on the cluster. I am going to use the Levenshtein distance to do this, and am going to store the value of distance normalized by the max length of the queries.
    count=0
    #print "Comparing:", sequences[0].id, "\t", sequences[1].id, ":\t", (1-float(distance(str(sequences[0].seq), str(sequences[1].seq)))/float(max(len(sequences[0].seq), len(sequences[1].seq))))
    #exit()
    for count in xrange(len(compareList)):
        #compareList[count][2] = float(0)
        compareList[count].value = (1-float(distance(str(sequences[compareList[count].first].seq), str(sequences[compareList[count].second].seq)))/float(max(len(sequences[compareList[count].first].seq), len(sequences[compareList[count].second].seq))))
        #print "Comparing:", sequences[int(compareList[count][0])].id, "\t", sequences[int(compareList[count][1])].id, "\t", compareList[count][2]
    #print compareList[0]
    #exit()

    #Write the list of distances to dat files with CSV writer. This can then be read by the Parse function in this script. 
    #Create a file in the Temp directory where the index list is called from.
    output_handle=open(line_file+"out","w+")
    output_csv = csv.writer(output_handle, delimiter='\t')
    countI=0
    for countI in xrange(len(compareList)):
        output_csv.writerow([compareList[countI].first] + [compareList[countI].second] + [compareList[countI].value])
            #output_handle.write(compareList[count][0]+"\t"+compareList[count][1]+"\n") 
    output_handle.close()
    print "Done. Check the output distance file in Temp folder."

#Function to read through the similarity scores and create the final file list.
def seqParse(fasta_file,line_dir,similarityCutoff,coreNumber):
    #Create hash table to read in sequences, another to write out the clean sequences and a third ID table to note replicates
    sequences={}
    sequencesClean={}
    sequencesID={}
    #Using the biopython fasta parse we can read our fasta input
    records = open(fasta_file, "rU")
    sequences = list(SeqIO.parse(records, "fasta"))
    records.close()
   
    #Create compare table of sequence indices and the Normalized Levenshtein distance
    compareList = [[[] for i in range(3)] for i in range(int(len(sequences)*(len(sequences)-1)/2))]
    #print "The size of the compare list is:\t", len(compareList), "\t", len(compareList[0])
    #Upload the comparisson indices
    #line.close()
    #print csvread.next()
    count=0
    for fileCount in xrange(coreNumber):
        readIndex = str(line_dir+str(fileCount+1)+".datout")
        line = open(readIndex, "rb")
        csvread = csv.reader(line, delimiter='\t') 
        for row in csvread:
            #print row
            compareList[count][0] = int(row[0])
            compareList[count][1] = int(row[1])
            compareList[count][2] = float(row[2])
            count+=1
        print "Reading:\t", readIndex
        line.close()
    print "Total sequences:\t", len(sequences)
    print "Total pair-wise distances calculated:\t", len(compareList)

    #Start the Clean sequence list and ID Hash table with the [0] index since it is never going to be thrown out.
    sequencesClean[0] = FASTAClass(sequences[0].id, sequences[0].seq)
    sequencesID[sequences[0].id] = sequences[0].id 
    #print "The ID file has:\t", sequencesID[sequences[0].id], "\t", sequences[0].id
    #exit()

    #Go through the comparelist and save the appropriate sequences by scanning through the file from top to bottom which is row-wise for every query.
    countI=0
    writeCount=0
    tempIndex=0
    tempCounter=1
    for countI in xrange(len(compareList)+1):
        #This is when you reach a new query_seq as seen in the first column of the compareList.
        if countI<len(compareList):
            if compareList[countI][0]!=tempIndex:
                if tempCounter==0:
                    writeCount+=1
                    sequencesClean[writeCount] = FASTAClass(sequences[tempIndex].id, sequences[tempIndex].seq)
                    sequencesID[sequences[tempIndex].id] = sequences[tempIndex].id 
                else:
                    tempCounter=0
                #This is to run the first row of a new query_seq since the count is skipped otherwise.
                tempIndex = compareList[countI][0] 
                #print "Working on:\t", sequences[tempIndex].id
                if compareList[countI][2]>=similarityCutoff:
                    tempCounter+=1
                    #If the id exists in the cleaned hash only then can I write.
                    if sequencesID[sequences[compareList[countI][1]].id]:
                        sequencesID[sequences[compareList[countI][1]].id]+="||"+sequences[compareList[countI][0]].id
            #This indicates we are on the same query_seq but with a different test_seq, which is from earlier in the list.
            elif compareList[countI][0]==tempIndex:
                if compareList[countI][2]>=similarityCutoff:
                    tempCounter+=1
                    try:
                        sequencesID[sequences[compareList[countI][1]].id]+="||"+sequences[compareList[countI][0]].id
                    except KeyError:
                        print sequences[compareList[countI][1]].id, " doesn't exist in the clean list. Skipping!"
        #This is for the last sequence since it never gets written.
        else:
            if tempCounter==0:
                writeCount+=1
                sequencesClean[writeCount] = FASTAClass(sequences[tempIndex].id, sequences[tempIndex].seq)
                sequencesID[sequences[tempIndex].id] = sequences[tempIndex].id 
        
    #Write out the clean FASTA file and the ID table in new files.
    output_handle = open("clean_"+fasta_file, "w+")
    count=0
    for count in xrange(len(sequencesClean)):
        output_handle.write(">"+sequencesClean[count].id+"\n"+str(sequencesClean[count].seq)+"\n")
    output_handle.close()
    print "Done writing the clean sequence list with ", len(sequencesClean), "sequences. Good Luck!"

    output_handle = open("id_"+fasta_file, "w+")
    for id in sequencesID:
        output_handle.write(">"+id+"\t"+sequencesID[id]+"\n")
    output_handle.close()
    print "Done writing the ID list of clean sequences. Good Luck!"

#Function to perform initial clean-up. Adapted from the Inter-webs.
def seqCleaner(fasta_file,minLength,maxLength):
    #create our hash table to add the sequences and ID table to note replicates
    sequences={}
    sequences_id={}

    #Using the biopython fasta parse we can read our fasta input
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        #Take the current sequence
        sequence=str(seq_record.seq).upper()
        #Check if the current sequence is according to the user parameters
        if (len(sequence)>=minLength and len(sequence)<=maxLength):
       # If the sequence passed in the test "is It clean?" and It isnt in the hash table , the sequence and Its id are going to be in the hash, and in the id-table
            if sequence not in sequences:
                sequences[sequence]=seq_record.id
		sequences_id[sequences[sequence]]=seq_record.id
       #If It is already in the hash table, We're just gonna concatenate the ID of the current sequence to another one that is already in the id-table
            else:
                sequences_id[sequences[sequence]]+="||"+seq_record.id
 
    #Write the clean sequences 
    #Create a file in the same directory where you ran this script
    output_file=open("clean_"+fasta_file,"w+")
    #Just Read the Hash Table and write on the file as a fasta format
    for sequence in sequences:
            output_file.write(">"+sequences[sequence]+"\n"+sequence+"\n")
    output_file.close()

    #Create a file in the same directory where you ran this script
    output_file=open("id_"+fasta_file,"w+")
    #Just Read the ID-Table and write on the file in a pseudo-fasta format
    for id in sequences_id:
            output_file.write(">"+id+"\t"+sequences_id[id]+"\n")
    output_file.close()
 
    print "Cleaning completed! Please check clean_"+fasta_file

#Funtion to query a FASTA file and eliminate all sequences present in a reference FASTA file.#
def seqQuery(queryFile,refFile):
    #create our hash table to add the sequences and ID table to note replicates
    sequencesQuery={}
    sequencesRef={}
    sequencesFinal={}

    #Using the biopython fasta parse we can read our fasta input
    records = open(queryFile, "rU")
    sequencesQuery = list(SeqIO.parse(records, "fasta"))
    records.close()

    records = open(refFile, "rU")
    sequencesRef = list(SeqIO.parse(records, "fasta"))
    records.close()
 
    #Using the biopython fasta parse we can read our fasta input
    countI=0
    cleanCount=0
    for countI in xrange(len(sequencesQuery)):
        #Compare the current sequence in the query to every sequence in the reference.
        countJ=0
        countCheck=0
        for countJ in xrange(len(sequencesRef)):
            if (str(sequencesQuery[countI].seq) == str(sequencesRef[countJ].seq)):
                countCheck+=1 
        if (countCheck==0):
            sequencesFinal[cleanCount] = FASTAClass(sequencesQuery[countI].id, sequencesQuery[countI].seq)
            cleanCount+=1
    print "A total of ", cleanCount, " sequences passed!"

    #Write the clean sequences 
    #Create a file in the same directory where you ran this script
    output_file=open("clear_"+queryFile,"w+")
    #Writing the sequences in the Final array.
    countI=0
    for countI in xrange(len(sequencesFinal)):
            output_file.write(">"+str(sequencesFinal[countI].id)+"\n"+str(sequencesFinal[countI].seq)+"\n")
    output_file.close()

    print "Query check completed! Please check clear_"+queryFile

#Function to sort through exceptions
def formatExceptionInfo(maxTBlevel=5):
    cla, exc, trbk = sys.exc_info()
    excName = cla.__name__
    try:
        excArgs = exc.__dict__["args"]
    except KeyError:
        excArgs = "<no args>"
    excTb = traceback.format_tb(trbk, maxTBlevel)
    return (excName, excArgs, excTb)
 
##Main code starts here##
#The code has 5 functions:
#    'query': Query check that eliminates sequences present in a reference FASTA file.
#    'clean': Pre-cleaner that eliminates sequences outside a window of acceptable size.
#    'split': Splitter that lists all pair-wise comparissons to be-done among the cores
#    'compare': Code that calculates Levenshtein distance of the pair lists generated.
#    'parse': Uses Levenshtein distance list to catalog all sequences below a similarity cut-off. 
userParameters=sys.argv[1:]
print userParameters
try:
    if (userParameters[0]=="help"):
        print "The code has 5 functions:"
        print "+\'query\'\tQuery check that eliminates sequences present in a reference FASTA file."
        print "+\'clean\'\tPre-cleaner that eliminates sequences outside a window of acceptable size."
        print "+\'split\'\tSplitter that lists all pair-wise comparissons to be-done among the cores."
        print "+\'compare\'\tCode that calculates Levenshtein distance of the pair lists generated."
        print "+\'parse\'\tUses Levenshtein distance list to catalog all sequences below a similarity cut-off."
    elif (userParameters[0]=="query"):
        if len(userParameters)==3:
            seqQuery(userParameters[1],userParameters[2])
        else:
            print "There is a problem with the inputs for Query-Ref check!"
    elif (userParameters[0]=="clean"):
        if len(userParameters)==4:
            seqCleaner(userParameters[1],int(userParameters[2]),int(userParameters[3]))
        else:
            print "There is a problem with the inputs for Pre-cleaning!"
    elif (userParameters[0]=="split"):
        if len(userParameters)==3:
            seqSplitter(userParameters[1],float(userParameters[2]))
        else:
            print "There is a problem with the inputs for Split!"
    elif (userParameters[0]=="compare"):
        if len(userParameters)==4:
            seqCompare(userParameters[1],userParameters[2],int(userParameters[3]))
        else:
            print "There is a problem with the inputs for Compare!"
    elif (userParameters[0]=="parse"):
        if len(userParameters)==5:
            seqParse(userParameters[1],userParameters[2],float(userParameters[3]),int(userParameters[4]))
        else:
            print "There is a problem with the inputs for Parse!" 
    else:
        print "There is a problem with inputs!"
except:
    print "Aaah! There is a problem!", sys.exc_info()[0], "/n", formatExceptionInfo() 
