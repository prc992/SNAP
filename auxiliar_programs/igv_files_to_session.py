#!/usr/bin/env python

#######################################################################
#######################################################################
## Created on Jan 11th 2025 to create IGV session file from file list
#######################################################################
#######################################################################

import os
import argparse
import pandas as pd

############################################
############################################
## PARSE ARGUMENTS
############################################
############################################

Description = 'Create IGV session file from a tdf files on the current folder and a set of genes.'
Epilog = """Example usage: python igv_files_to_session.py <XML_OUT> <LIST_GENES> <GENOME>"""

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

## REQUIRED PARAMETERS
argParser.add_argument('XML_OUT', help="XML output file.")
argParser.add_argument('LIST_GENES', help="Tab-delimited file gene or locations four columns i.e. chr\tstart\tend\tname. Header isnt required.")
argParser.add_argument('GENOME', help="Full path to genome fasta file or shorthand for genome available in IGV e.g. hg19.")

args = argParser.parse_args()

############################################
############################################
## HELPER FUNCTIONS
############################################
############################################

def makedir(path):

    if not len(path) == 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

############################################
############################################
## MAIN FUNCTION
############################################
############################################

def igv_files_to_session(XMLOut,ListGenes,Genome,PathPrefix=''):

    makedir(os.path.dirname(XMLOut))

    # Define the file extension to search for
    file_extension = ".bw"
    # Get the current working directory
    current_directory = os.getcwd()
    # Find all files with the specified extension in the current directory
    fileList = [file for file in os.listdir(current_directory) if file.endswith(file_extension)]

    ## ADD RESOURCES SECTION
    XMLStr =  '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'
    XMLStr += '<Session genome="%s">\n' % (Genome)
    XMLStr += '\t<Resources>\n'
    for ifile in fileList:
        XMLStr += '\t\t<Resource path="%s" type="bw"/>\n' % (ifile)
    XMLStr += '\t</Resources>\n'

    ## ADD PANEL SECTION
    XMLStr += '\t<Panel name="DataPanel">\n'
    for ifile in fileList:
        extension = os.path.splitext(ifile)[1].lower()
        nameSample = ifile.replace(file_extension,'')
        XMLStr += '\t\t<Track altColor="0,0,178" autoScale="true" clazz="org.broad.igv.track.FeatureTrack" '
        XMLStr += 'fontSize="10" height="60" id="%s" '% (ifile)
        XMLStr += 'name="%s" renderer="BASIC_FEATURE" sortable="false" visible="true" windowFunction="count"/>\n' % (nameSample)
    XMLStr += '\t</Panel>\n'

    #Create Gene List
    GeneList = []
    gen = open(ListGenes,'r')
    while True:
        line = gen.readline()
        if line:
            strchr,start,end,name = line.strip().split('\t')
            GeneList.append((strchr,start,end,name))
        else:
            break
            fout.close()
    
    ## ADD GENE REGION SECTION
    XMLStr += '\t<PanelLayout dividerFractions="0.87"/>'
    XMLStr += '\t<Regions>\n'
    for strchr,start,end,name in GeneList:
        XMLStr += '\t\t<Region chromosome="%s" ' % (strchr)
        XMLStr += ' description="%s" ' % (name)
        XMLStr += ' end="%s" ' % (end)
        XMLStr += ' start="%s"/> \n' % (start)
    XMLStr += '\t</Regions>\n'

    ## ADD REGION OF INTEREST SECTION
    XMLStr += '\t<GeneList name="Regions of Interest">\n'
    for strchr,start,end,name in GeneList:
        XMLStr += "\t\t%s:" % (strchr)
        XMLStr += "%s-" % (start)
        XMLStr += "%s\n" % (end)

    ## ADD FRAME SECTION
  
    for strchr,start,end,name in GeneList:
        XMLStr += '\t\t<Frame chr="%s" ' % (strchr)
        XMLStr += ' end="%s" ' % (end)
        XMLStr += ' name="%s:' % (strchr)
        XMLStr += '%s-' % (start)
        XMLStr += '%s"' % (end)
        XMLStr += ' start="%s"/>\n ' % (start)
    
    XMLStr += '\t</GeneList>\n'

    XMLStr += '</Session>'
    XMLOut = open(XMLOut,'w')
    XMLOut.write(XMLStr)
    XMLOut.close()

############################################
############################################
## RUN FUNCTION
############################################
############################################

igv_files_to_session(XMLOut=args.XML_OUT,ListGenes=args.LIST_GENES,Genome=args.GENOME)

############################################
############################################
############################################
############################################
