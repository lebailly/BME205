#!/usr/bin/evn python2.7

#File: lebailly/BME205/HW9/parse.py
#Author: Chris LeBailly

"""
parse module contains several useful parsers.

fasta(source, alphabet="ACDEFGHIKLMNPQRSTVWYBZX*-", case_senitive=False) parses 
fasta from source, which is a file-like object with fasta sequences.  Yields 
(name, comment, sequence), where all three are strings with no line breaks.

fastq(source, phred, alphabet="ACDEFGHIKLMNPQRSTVWYBZX*-"): parses fastq from
source, which is a file-like object with fastq sequences. 
codon_freq_table(source) reads the frequency of codons for different species.
Source is a file-like object in a style like CodonFrequency output in GCG
Wisconsin Package.  Returns a Counter of counts for the codons of the species.

genetic_code(source, name='Standard') reads the genetic code.  Source is a file-
like object in the same format as ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt
"""

#ATTEN - Finsih module doc string.

import collections, re, sys, urllib2

def codon_freq_table(source):
    """
    Pre-condition: Source is a file-like source of a website in a style like 
    CodonFrequency output in GCG Wisconsin Package (such as http://www. kazusa.
    or.jp/codon/cgi-bin/showcodon.cgi?species=199310&aa=1&style =GCG).

    Post-condition: Returns a collections.Counter of counts for the codons.
    """

    #Counts is a Counter with codons as keys and counts as values.
    counts = collections.Counter()

    #reading_table tracks when we are reading the codon frequency table.
    reading_table = False

    for line in source:
        if(line.startswith('<PRE>')): reading_table = True
        elif(line.startswith('</PRE>')): reading_table = False
        elif(reading_table):
            split_line = line.split()

            if(len(split_line) != 0 and split_line[0] != 'AmAcid'):
                freq = float(split_line[2])
                codon = split_line[1]
                counts[codon] = freq

    return counts

def genetic_code(source, name='Standard'):
    """
    Pre-condition: source is a file-like source in the same format as 
    ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt. name is the name of the
    desired genetic code to read.  Must be the first name in the entry (not the
    second name in the entry, which some genetic codes have).

    Post-condition: Returns a dictionary with codons as keys and the amnio acids 
    they code for as values.
    """

    #Dictionary with codons as keys and the amnio acids they code for as values.
    codon_table = {}

    #reading_code tracks when we are reading the desired genetic code entry.
    reading_code = False

    for line in source:
        if(not line.startswith('--')):
            line = line.rstrip()
            if(line.startswith('  name "')):
                line = line.lstrip('  name "').rstrip('" ,')
                if(line == name): 
                    reading_code = True
                    line = source.next() #reads next line, which may be a name.
                else: reading_code = False
            if(reading_code):
                if(line.startswith('  ncbieaa')):
                    line = line.lstrip('  ncbieaa  "').rstrip('",')
                    aa_str = line
                elif(line.startswith('  -- Base1')):
                    line = line.lstrip('  -- Base1')
                    Base1 = line
                elif(line.startswith('  -- Base2')):
                    line = line.lstrip('  -- Base2')
                    Base2 = line
                elif(line.startswith('  -- Base3')):
                    line = line.lstrip('  -- Base3')
                    Base3 = line

    for index, aa in enumerate(aa_str):
        codon_table[Base1[index]+Base2[index]+Base3[index]] = aa

    return codon_table

def codon_preference_table(source): #ATTEN - merge with codon_freq_table
    """
    Pre-condition: 'source' is a file-like object with the format:
    First line is a header line.  Remaining lines are either white space
    or whitespace-separated data with the codon in the second column 
    and the counts in the third column.

    Post-condition: Returns a collections.Counter of counts for each codon.
    """

    #Counts contains codons a keys and counts as values.
    counts = collections.Counter()

    fist_line = source.readline()

    for line in source:
        split_line = line.split()

        if(len(split_line) != 0):
            freq = float(split_line[2])
            codon = split_line[1]
            counts[codon] = freq

    return counts

def codon_preference_table_with_AA(source):
    """
    Pre-condition: 'source' is a file-like object with the format:
    First line is a header line.  Remaining lines are either white space
    or whitespace-separated data with the codon in the second column 
    and the number if counts in the third column.

    Post-condition: Returns 'counts', where 'counts' has amino-acid keys and 
    Counters as values.  The inner Counter has codons as keys and counts as 
    values.
    """

    #'counts' is a dict of dict described in docstring.
    counts = {}
    codon_table = genetic_code(urllib2.urlopen(urllib2.Request(
                          'ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt')))

    for AA in 'ACEDGFIHKMLNQPSRTWVY*':
        counts[AA] = collections.Counter()

    fist_line = source.readline()

    for line in source:
        split_line = line.split()

        if(len(split_line) != 0):
            freq = float(split_line[2])
            codon = split_line[1]
            AA = codon_table[codon]
            counts[AA][codon] = freq

    return counts

def subst(source):
    """
    Pre-condition: source is a file-like object containg integer values in 
    BLOSUM format.
    
    Post-condition: Returns the substituion matrix and the alphabet as a tuple. 
    The substitution matrix is a dictionary with keys of tuples of two amnio 
    acids (such as ('M', 'A')) and a vlues of the substituion value for matching 
    those two amnio acids.  Also returns a string of the alphabet used in the 
    input substituion matrix.
    """

    #subst is the dictionary described in the docstring.
    subst = {}

    #alpha is a list of characters in the amnio acid alphabet used in source.
    alpha = None

    for line in source:
        if(not line.startswith('#')):
            if(not alpha):
                alpha = line.split()
            else:
                line_list = line.split()

                first_char = line_list.pop(0)

                for index, second_char in enumerate(alpha):
                    subst[(first_char,second_char)] = int(line_list[index])

    source.close()

    return subst, ''.join(alpha)

def fasta(source, alphabet="ACDEFGHIKLMNPQRSTVWYBZX*-", case_senitive=False):
    """
    Pre-condition: source is a file-like object containing fasta sequencs. 
    alphabet is a string containing the letters in the desired alphabet.
    If case_senitive is True, preserve case of sequences.  If it is false,
    all characters are reported as upper case characters.

    Post-condition: yields (name, comment, sequence), where all three are 
    strings with no line breaks.
    """

    #ID is the name of the sequence
    #seq is a list containing the characters of the sequence
    #comment is a list contain comments contained in the sequence

    #separator is a regular expression use separator the ID and comment.
    #alpha_re is a regular expression for the alphabet (both upper and 
    #lower case) specified.


    #ATTEN - ADD ENUMERATE TO GIVE LINE NUMBER IN ERROR OUTPUT.

    ID, seq, comment = None, [], []
    separator = re.compile('[, \t\r]')
    alpha_re = re.compile('[{}{}]'.format(alphabet, alphabet.lower()))

    for line in source:
        line = line.rstrip()

        #Start of new sequence
        if line.startswith('>'):

            _fasta_error_check(ID, seq)

            if ID: yield (ID, ' '.join(comment), ''.join(seq))

            seq, comment = [], []

            first_line = separator.split(line, maxsplit = 1)
            ID = first_line[0].lstrip('>')
            if len(first_line) == 2:
                comment.append(first_line[1])

        #Adds any semi-colon comments        
        elif line.startswith(';'): comment.append(line.lstrip(';'))

        #Otherwise line must be part of a sequence    
        else:
            if not case_senitive: clean_line = re.findall(alpha_re,line.upper())
            else: clean_line = re.findall(alpha_re, line)
            
            seq += clean_line

    _fasta_error_check(ID, seq)

    if ID: yield (ID, ' '.join(comment), ''.join(seq))

def fastq(source, phred=33, alphabet="ACDEFGHIKLMNPQRSTVWYBZX*-"):
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
    #ATTEN - ADD BETTER DOCUMENTATION!  SPELL CHECK

    ID, seq, qual = None, [], []
    seq_len = qual_len = 0
    done_reading_seq = phred_error = False
    comment = ''

    separator = re.compile('[, \t\r]')
    alpha_re = re.compile('[{}]'.format(alphabet))

    for line in source:
        line = line.rstrip()

        #Test for ID line
        if(line.startswith('@') and seq_len == qual_len):

            _fastq_error_check(ID, seq, phred_error, phred)

            if ID: yield (ID, ''.join(seq), comment, qual)

            seq_len, qual_len, seq, qual = 0, 0, [], []
            phred_error = False
            done_reading_seq = False
            comment = ''

            first_line = separator.split(line, maxsplit = 1)

            ID = first_line[0].lstrip('@')
            
            if len(first_line) == 2:
                comment = first_line[1]

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

                #ATTEN - ARE THESE BOUNDS CORRECT?
                #ATTEN - TRY TO ELINIATE PHRED BOOL
                if(phred == 33 and ord(x) < 33):
                    phred_error = True
                elif(phred == 64 and ord(x) < 64):
                    phred_error = True

            qual_len += len(line)

    _fastq_error_check(ID, seq, phred_error, phred)

    if ID:
        if(qual_len != seq_len):
            sys.stderr.write('''WARNING: Sequence not same length as quality in 
                {}, but last sequence in input.\n'''.format(ID))
        yield (ID, ''.join(seq), comment, qual)

def fasta_with_quality(fasta_source, qual_source, phred, alphabet):
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
    fasta_generator = fasta(fasta_source, alphabet)


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

def _fasta_error_check(ID, seq):
    """ Send warnings for empty sequences and nameless sequences to stderr """

    if(ID and not seq):
        sys.stderr.write("WARNING: {} is an empty sequence.\n".format(ID))
    elif(not ID and seq):
        sys.stderr.write("WARNING: Sequence with no name in input.\n")

def _fastq_error_check(ID, seq, phred_error, phred):
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