"""
This program uses fastq file from Illumina sequencing as input, check average quality of all reads, trim unsatisfied
reads to desired length, map trimmed reads back to a reference genome, and call peaks based on given mapping threshold
and peak distribution width.
"""

import numpy as np
import matplotlib.pyplot as plt
from math import *

#read input files (fasta or fastaq)
def read():
    test_file = open("/Users/YanpengXi/Desktop/Python_final/Data/CA_R2.fastq","r")
    ctrl_file = open("/Users/YanpengXi/Desktop/Python_final/Data/LY_R1.fastq","r")
    chr3_file = open("/Users/YanpengXi/Desktop/Python_final/Data/chr3.fasta","r")

    test = test_file.readlines()
    ctrl = ctrl_file.readlines()
    chr3 = chr3_file.read()

    test = test[0:400] # testing first
    ctrl = ctrl[0:400]

    test_file.close()
    ctrl_file.close()

    return (test,ctrl,chr3)

# check average quality of all reads, plot quality against nucleotide position (101bp/read)
def quality(seq):
    # fastq quality score range 
    quality_str =  "#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRS"\
                   "TUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"

    quality = {'!':1,'"':2}
    for i in range (len(quality_str)):
        quality[quality_str[i]] = (i+3)

    score = []
    error = []
    for i in range(len(seq[1])):
        score.append(0)
        error.append(0)

    for i in range (0,len(seq)-3,3):
        for j in range (len(seq[i])-1):
            try:
                score[j] = score[j]+quality[seq[i][j]]
            except KeyError:
                break

    for i in range (len(score)):
        score[i] = score[i]/(len(seq)/4.0)

    for i in range (0,len(seq)-3,3):
        for j in range (len(seq[i])-1):
            try:
                error[j] = error[j]+(quality[seq[i][j]]-score[j])**2
            except KeyError:
                break
    for i in range (len(error)):
        error[i] = sqrt(error[i]/(len(seq)/4.0-1))

    # quality plot 
    x = np.arange(1,len(seq[1])+1,1)
    y = score
    plt.plot(x,y)
    plt.errorbar(x,y,error,fmt='-o')
    plt.xlabel("Nucleotide position")
    plt.ylabel("Average quality score")
    axes = plt.gca()
    axes.set_xlim([0,len(seq[1])])
    axes.set_ylim([1,94])
    plt.show()

# trim left and right borders of all reads
def trim(seq,left,right):
    for i in range (len(seq)):
        if seq[i][0]!="@" and seq[i][0]!= "+":
            seq[i] = seq[i][left:right+1]
    return(seq)


# map trimmed reads to reference genome (here only show chromosome 3 as example)
def map(test,draw):
    map_file = open("/Users/YanpengXi/Desktop/Python_final/Data/CA_R1.sam","r")
    mapped = map_file.readlines()
    map_file.close()
    Chr3 = {}
    for line in mapped:
        if line.split()[2] == "Chr3":
            chr = line.split()[2]
            location = line.split()[3]
            if Chr3.get(str(location)) == None:
                 Chr3[str(location)] = 1
            else:
                Chr3[str(location)]  = Chr3[str(location)]  +1

    location = []
    count = []
    y = []
    for element in Chr3.keys():
        if Chr3[element]<=40:
            y.append(int(element))
    y.sort()
    for element in y:
        location.append(element)
        count.append(Chr3[str(element)])

    if draw == "y":
        x = location
        y = count
        plt.plot(x,y)
        plt.xlabel("Arabidopsis Chromosome #3")
        plt.ylabel("Mapped reads")
        plt.show()

    if draw == "n":
        return (location,count)

# call peak based on given threshold and band width (no control group currently)
def call_peak(location,count,threshold,width):
    peak_location = []
    peak_count = []

    for i in range (len(location)):
        if count[i] >= threshold:
            if threshold>=count[i+width] >= threshold/2.0 and threshold>=count[i-width] >= threshold/2.0:
                peak_location.append(location[i])
                peak_count.append(count[i])

    x = peak_location
    y = peak_count
    plt.plot(x,y,'ro')
    plt.xlabel("Arabidopsis Chromosome #3")
    plt.ylabel("Called peaks")
    plt.show()

    out = open("/Users/YanpengXi/Desktop/Python_final/Data/Peaks.txt","w")
    out.write("Chromosome......Peak location (bp)\n")
    out.write("----------------------------------\n")
    for element in peak_location:
        out.write("Chromosome3......"+str(element)+"\n")
    out.close()

# introduction
def intro():
    print (""
    "||ChIP-seq Peak Caller ||\n"
    "------------------------------------\n"
    "quality()                      ####[check quality of sequence reads]\n"
    "trim(left,right)               ####[trim left and right border (bp) of all reads]\n"
    "map()                          ####[map reads to reference genome]\n"
    "call_peak(threshold, width)    ####[output write into peaks.txt]\n"
    "quit                           ####[quit program]\n"
    "-----------------------------------\n")


def main():

    test,ctrl,chr3 = read()
    user = ""
    # recurrently ask user to input a function unless "quit" is given. Typo will be noticed as "Invalid input".
    while user != "quit":
        intro()
        user = input("Please choose a function: ")
        if user == "quality()":
            quality(test)
        elif user[0:4] =="trim":
            user.replace(" ","")
            left = int(user[user.index("(")+1:user.index(",")])
            right = int(user[user.index(",")+1:user.index(")")])
            trimmed = trim(test,left,right)
            quality(trimmed)
        elif user =="map()":
            map(test,"y")
        elif user[0:4] =="call":
            user.replace(" ","")
            threshold = int(user[user.index("(")+1:user.index(",")])
            width = int(user[user.index(",")+1:user.index(")")])
            location,count = map(test,"n")
            call_peak(location,count,threshold,width)
        elif user == "quit":
            break
        else:
            print ("Invalid input!!")

main()