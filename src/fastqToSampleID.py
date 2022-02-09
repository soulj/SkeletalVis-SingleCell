#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import csv
import sys

script, sampleTable, SRR = sys.argv

def getSampleID(ID,sampleTable):
    
    ID = ID.split('.')[0]
    
    #open the sampleTable
    with open(sampleTable) as tsvfile:
            reader = csv.DictReader(tsvfile, delimiter='\t')
            
            #get all the possible file names
            for row in reader:
                fileNames=str.split(row["File"],"|")
                
                #match the fileNames and get the sampleID
                for fileName in fileNames:
                    if fileName in ID:
                        sys.stdout.write(row["Sample"])                      
            
            
if __name__ == "__main__":
    getSampleID(SRR,sampleTable)
