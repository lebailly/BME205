#!/usr/bin/env python2.7

#File: lebailly/BME205/HW4/Markov.py
#Author: Chris LeBailly

"""
Markov.py is a module consisting of functions related to Markov models.

functions:
get_counts(source, order, alpha, start_char, stop_char) produces a dictionary
with keys of k-mers and values of counts.

make_kmer_dict(source) take a file-like source and produces a dictionary with 
k-mer keys and count values.

compute_counts(kmer_dict, order, alpha, start_char, stop_char, pseudo_value)
produces a dict of dicts with counts from kmer_dict plus pseudo_value.

compute_log_probs(counts) converts counts to log probabilities.

coding_cost(log_probs, order, fasta, alpha, start_char, stop_char)
computes the coding cost, total characters, and total sequences in fasta.

Auxiliary functions (not needed by user):
create_pseudo_counts(order, pseudo_value, alpha, start_char, stop_char)
produces a dict of dicts with only pseudo values

print_table(table) debugging function to print dict of dicts.
"""

from __future__ import division, print_function
import argparse, collections, itertools, imp
parse = imp.load_source('parse','../parse.py')
from copy import deepcopy
from math import log

def get_counts(source, order, alpha, start_char, stop_char):
	""" 
	Counts number of k-mers in fasta source ('source'), where k = 'order'+1.
	'alpha' is a string containing the desired alphabet with start character 
    'start_char' and stop character 'stop_char'.
	
	Returns a dictionary with keys of k-mers and values of counts.
	"""

    #'counts' has k-mers for keys and counts as values
	counts = collections.Counter()

	for (fasta_id, comment, seq) in parse.fasta(source, alpha):

		#sequence is long enough start and stop don't appear in same k-mer.
		if(len(seq) >= order):

			#counts k-mers with no stop or start characters
			for start in range(len(seq)-order):
				counts[seq[start:start+order+1]] += 1

			#counts k-mers with start character and k-mers with stop characters
			for start in range(order):
				counts[start_char*(order - start)+seq[0:start+1]] += 1
				counts[seq[len(seq)-start-1:]+stop_char*(order-start)] += 1

		#sequence is short enough that start and stop appear in the same k-mer.	
		else:
			total = order+1-len(seq) #total number of start & stop chars

			#counts k-mers with start and stop characters in same k-mer
			for num_start in range(total+1):
				counts[num_start*start_char+seq+(total-num_start)*stop_char] += 1

		#counts stop character for order 0 model
		if(order == 0):
			counts['$'] += 1

	return counts

def make_kmer_dict(source):
    """
    Given a file-like object ('source') reads in whitespace-separated kmer-count 
    pairs, one per line.  Returns a dictionary with k-mer keys and count values.
    """

    #kmer_dict is a dictionary with k-mer keys and count values.
    #kmer_count is a list (created from splitting source line)
    #kmer_count[0] is the k-mer and kmer_count[1] is the count.

    kmer_dict = {}

    for line in source:
        kmer_count = line.split()

        if(len(kmer_count) == 2):
            kmer_count[1] = int(kmer_count[1])
            kmer_dict[kmer_count[0]] = kmer_count[1]
        else:
            raise Exception(
                "Kmer count input not in valid format (need 'kmer count \\n')")

    return kmer_dict

def create_pseudo_counts(order, pseudo_value, alpha, start_char, stop_char):
    """
    Expects the model order ('order'), 'pseudo_value' for each cell, 
    the alphabet which is being used ('alpha'), the start character, 
    ('start_char'), and the stop character ('stop_char').

    Creates a dictionary of dictionaries.  Generates a generic row, which is a 
    dictionary containing each letter in alpha and the stop character as keys 
    and pseudo_value as values.  The outer dictionary contains all possible 
    context (strings of length 'order' of letters in 'alpha'), including start 
    characters at the beginning of the context.  The generic row is copied and 
    stored as the value for each key (which is a context).  

    This dict of dicts is returned.
    """

    #pseduo_row is dict with alpha+stop_char as keys and pseudo_value values.
    #pseudo_counts is the dict of dict described in docstring.

    pseudo_row = {}

    for letter in alpha+stop_char:
        pseudo_row[letter] = pseudo_value

    pseudo_counts = {}

    #Generates all possible conditions
    if(order >= 0):
        for num_starts in range(order+1):
            for condition in itertools.product(alpha, repeat=order-num_starts):
                pseudo_counts[start_char*num_starts+''.join(
                    condition)] = deepcopy(pseudo_row)

    else:
        raise Exception("Error: Cannot create pseudo_counts with negative order")

    return pseudo_counts

def make_counts_table(kmer_dict, order, alpha, start_char, stop_char, 
    pseudo_value):
    """
    Expects a dictionary with kmer keys and count values (kmer_dict), the order
    of the model ('order'), the alphabet being used ('alpha'), the start 
    character ('start_char'), the stop character ('stop_char'), and the 
    pseudo_value to be filled in for each cell ('pseudo_value').

    Creates a dictionary of dictionary where the outer dictionary contains each
    possible context and the inner contains each letter in alpha and the stop 
    character.  The value is set initially at the pseudo_value.  The counts from 
    kmer_dict added to this table where the context is kmer[0:-1] and the 
    character in the row is kmer[-1].  This dict of dicts is returned.
    """

    #'counts' is a dict of dicts
    counts = create_pseudo_counts(
        order, pseudo_value, alpha, start_char, stop_char)

    for kmer, freq in kmer_dict.items():
        if((order > 0 and kmer[0:-1][-1] != stop_char and 
            kmer[-1] != start_char) or order == 0):
            if(kmer[0:-1] not in counts.keys()):
                raise Exception(
                    "K-mer context {} contains letters not in alphabet!"
                    .format(kmer[0:-1]))

            if(kmer[-1] not in counts[kmer[0:-1]].keys()):
                raise Exception("K-mer contains {} but is not in alphabet!"
                    .format(kmer[-1]))
            
            counts[kmer[0:-1]][kmer[-1]] += freq

    return counts

def compute_log_probs(counts):
    """
    Expects a dictionary of dictionary ('counts') where each inner dictionary is
    a row in the counts table.  This row is normalized (converted to 
    probabilities) such that the row sum is 1.  It is then converted to 
    -log(probability) so logs are not computed in coding_cost and to prevent 
    underflows in probability. The dictionary of dictionary with 
    log-probabilities is returned.
    """

    #log_probs is a dict of dicts
    log_probs = {}

    for context in counts.keys():
        row = {}
        total = sum(counts[context].values())
        for kmer, freq in counts[context].items():
            row[kmer] = log(total)-log(freq) 
        log_probs[context] = row

    return log_probs

def coding_cost(log_probs, order, fasta, alpha, start_char, stop_char):
    """
    Expects a dict of dicts of log-probabilities ('log_probs', outer dict is 
    context), the 'order' of the model, a 'fasta' file-like source, the alphabet
    ('alpha') used in making log_probs, the start character ('start_char'), and 
    the stop character ('stop_char').

    This function computes total coding cost (based on 'log_probs') of each 
    sequence in 'fasta'.  It uses 'order' start characters and a single stop 
    character. The function also counts to total number of characters (not 
    including any start or stop characters), and the total number of sequences.

    Returns a tuple of total cost (in log base 2), total characters, 
    and the total number of sequences.
    """
    #'seq_with_start_stop' is a string which starts with 'order' number of start
    #characters, the sequence, then one stop character.

    total_cost = total_chars = total_seqs = 0

    for name, comment, seq in parse.fasta(fasta, alpha):
        total_chars += len(seq)
        total_seqs += 1

        seq_with_start_stop = start_char*order + seq + stop_char

        for start in range(len(seq_with_start_stop)-order):
            kmer = seq_with_start_stop[start:start+order+1]
            total_cost += log_probs[kmer[0:-1]][kmer[-1]]

    return total_cost/log(2), total_chars, total_seqs

def print_table(table):
    """Debugging function to print dict of dicts of counts or log-probs"""

    #'total' is the row total

    for context, row in sorted(table.items()):
        sorted_row = sorted(row.items())
        total = sum(table[context].values())
        print('The context', context ,'has a row total of', total, 
            'and the row itself is\n', sorted_row)