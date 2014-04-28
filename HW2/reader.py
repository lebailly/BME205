#!/usr/bin/env python2.7

#File: lebailly/BME205/HW2/reader.py 
#Author: Chris LeBailly

"""
reader is a module which contains functions to read fasta files, read fastq
files, and to read fasta with quality files.
"""

from __future__ import division, print_function
import sys, re

#Looked at the site below when dealing with an early bug.
#http://stackoverflow.com/questions/7654971/parsing-a-fasta-file-using-a-generator-python

def read_fasta(source, alphabet):
    """
    'source' is a fasta file-like object. 'chars_not_in_seq' is a string (for 
    use in a regex) of characters outside desired sequence alphabet.

    generator function yields (name, comment, sequence), where all three are 
    strings with no line breaks.

    Warning given empty sequences and nameless sequences.
    """

    ID, seq, comment = None, [], []
    separator = re.compile('[, \t\r]')
    alpha_re = re.compile('[{}]'.format(alphabet))

    for line in source:
        line = line.rstrip()

        #Start of new sequence
        if line.startswith('>'):

            fasta_error_check(ID, seq)

            if ID: yield (ID, ' '.join(comment), ''.join(seq))

            seq, comment = [], []

            first_line = separator.split(line, maxsplit = 1)
            ID = first_line[0].lstrip('>')
            if len(first_line) == 2:
                comment.append(first_line[1])

        #Adds any semi-colon comments        
        elif line.startswith(';'):
            comment.append(line.lstrip(';'))

        #Otherwise line must be part of a sequence    
        else:
            clean_line = re.findall(alpha_re,line.upper()) #list of characters #Add .upper as option (?)
            
            seq += clean_line

    fasta_error_check(ID, seq)

    if ID: yield (ID, ' '.join(comment), ''.join(seq))

def fasta_error_check(ID, seq):
    """ Send warnings for empty sequences and nameless sequences to stderr """

    if(ID and len(seq) == 0):
        sys.stderr.write("WARNING: "+ID+" is an empty sequence.\n")
    elif(ID == '' and len(seq) != 0):
        sys.stderr.write("WARNING: Sequence with no name in input.\n")

def read_fastq(source, phred, alphabet):
    """
    'source' is a fastq file-like object. 'phred' specifies if 'source' is
    Phred+33 or Phred+64. 'chars_not_in_seq' is a string (for use in a regex) 
    of characters outside desired sequence alphabet.

    generator function yields (name, sequence, comment, quality).  Name,
    sequence and comment are strings with no line breaks.  quality is a list of
    Phred values.

    Warning given empty sequences, nameless sequences, and sequence characters
    outside of desired alphabet.  These charactersare removed from output.

    Warning given if invalid quality characters are detected.

    Error given if sequence differs from quality in length and stops reading.
    """

    #ATTEN - ADD ENUMERATE TO GIVE LINE NUMBER IN ERROR OUTPUT.

    ID, seq, qual = None, [], []
    seq_len = qual_len = 0
    done_reading_seq = phred_error = False
    comment = ''

    separator = re.compile('[, \t\r]')
    alpha_re = re.compile('[{}]'.format(alpha))

    for line in source:
        line = line.rstrip()

        #Test for ID line
        if(line.startswith('@') and seq_len == qual_len):

            fastq_error_check(ID, seq, phred_error, phred)

            #print(seq)

            if ID: yield (ID, ''.join(seq), comment, qual)

            first_line = separator.split(line, maxsplit = 1)

            ID = first_line[0].lstrip('@')
            
            if len(first_line) == 2:
                comment = first_line[1]

            seq_len, qual_len, seq, qual = 0, 0, [], []
            phred_error = False
            done_reading_seq = False
            comment = ''

        #Tests for quality not matching sequence length
        elif(qual_len > seq_len):
            raise Exception('''Sequence and quality different length in {}, 
            can not read in remaining data.\n'''.format(ID))
            

        #Test for Comment line
        elif(line.startswith('+') and not done_reading_seq):
            first_line = separator.split(line.lstrip('+'), maxsplit = 1)
            if(first_line[0] != '' and first_line[0] != ID):
                sys.stderr.write('Warning: Sequence ID @{} does not match +{}\n'
                    .format(ID,first_line[0]))
            done_reading_seq = True

        #Test for Sequence lines (no comment has yet been read)
        elif(not done_reading_seq):
            clean_line = re.findall(alpha_re, line.upper()) #Add .upper as option (?)

            seq += clean_line
            seq_len += len(clean_line)

        #Test for Quality lines
        else:
            for x in line:
                qual.append(ord(x)-phred)

                if(phred == 33 and ord(x) < 33):
                    phred_error = True
                elif(phred == 64 and ord(x) < 64):
                    phred_error = True

            qual_len += len(line)

    fastq_error_check(ID, seq, phred_error, phred)

    if ID:
        if(qual_len != seq_len):
            sys.stderr.write('''WARNING: Sequence not same length as quality in 
                {}, but last sequence in input.\n'''.format(ID))
        yield (ID, ''.join(seq), comment, qual)

def fastq_error_check(ID, seq, phred_error, phred):
    """ 
    Checks for empty or nameless sequence and gives warning if found. 
    Gives warning for bad quality character.
    """

    if(ID and len(seq) == 0):
        sys.stderr.write("WARNING: {} is an empty sequence.\n".format(ID))
    elif(ID == '' and len(seq) != 0):
        sys.stderr.write("WARNING: Sequence with no name in input.\n")
    if phred_error:
        sys.stderr.write("WARNING: Invalid Phred+{} value in {}"
            .format(phred,ID))

def read_fasta_with_quality(fasta_source, qual_source, phred, alphabet):
    """
    'fasta_source' is a fasta file-like object. 'qual_source' is a file with a 
    fasta-type ID line, followed by white-space-separated phred values. 
    'alphabet' is a string of characters 
    in the desired sequence alphabet.

    generator function yields (name, sequence, comment, quality)

    fasta sequence obtained from read_fasta function.

    Warning produced if invalid Phred value is read in.  This value is changed
    to 0 to prevent problems if these values are changed to characters.  Since
    the value is invalid it is assumed that nothing is none at that node.
    """


    ID, qual, phred_error = None, [], False
    fasta_generator = read_fasta(fasta_source, alphabet)


    for line in qual_source:
        if line.startswith('>'):
            if ID: yield (name, seq, comment, qual)
            ID, qual = line.lstrip('>'), []
            name, comment, seq = next(fasta_generator)
        else:
            phred_error = False
            for x in line.split():
                if (phred == 33 and (int(x) < 0 or int(x)>93) or 
                    phred == 64 and (int(x) < 0 or int(x)>62)):

                    qual.append(0)
                    phred_error = True
                    
                else:
                    qual.append(int(x))

        if phred_error:
            sys.stderr.write('''WARNING: Invalid Phred+{} value in {} 
                (invalid value changed to 0) \n'''.format(phred,name))

    if ID: yield (name, seq, comment, qual)
