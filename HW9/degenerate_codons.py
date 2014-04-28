#!/usr/bin/env python2.7

#File: lebailly/BME205/HW9/degenerate_codons.py
#Author: Chris LeBailly

"""
module degenerate_codons.py contains the degenerate_codons class and
the standard genetic code (standard_code, which is dictionary with codons as
keys and the amino acid which the codon codes for).
"""

from __future__ import division, print_function
import itertools, collections

#Ian Fiddes typed this (emailed me to save typing time)
standard_code = {'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L',
    'CTA':'L','CTG': 'L', 'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V',
    'GTC':'V','GTA':'V','GTG':'V','TCT': 'S','TCC':'S','TCA':'S','TCG':'S',
    'CCT':'P','CCC':'P','CCA':'P','CCG':'P','ACT':'T','ACC':'T','ACA':'T',
    'ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A','TAT':'Y','TAC':'Y',
    'TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','AAT':'N',
    'AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R',
    'CGG':'R','AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G',
    'GGA':'G','GGG':'G'}

class degenerate_codons(object):
	"""
	degenerate_codons class constructs dictionaries to find a mapping from 
	amino acid sets to the degenerate codons that code for the set.  The 
	compute_results method computes the imbalance and frequency for all codons.  
	output prints results to stdout (see output docstring for details).

	Methods:
	expand(dcodon): Returns a list of codons which dcodon expands to.
	imbalance(dcodon): Returns the imbalance of dcodon.
	freq(dcodon): Returns the average frequency of dcodon.
	compute_results(): Computes the imbalance and frequency for each codon in 
	each amino acid set.  Sorts results for each amino acid set by smallest
	imbalance, then by largest frequency.
	output(output_format): Prints results to stdout.  Format can be either
	'min-codons', 'full', 'all-codons', or 'minimal' (see docstring for details).
	best_codon(AAset): Given any set of amino acids and '*', returns the codon
	with minimal imbalance and highest frequency.  If none exists, adds amino 
	acids the set till one is found.
	"""

	def __init__(self, codon_counts, genetic_code=standard_code):
		"""
		Pre-condition: codon_counts is a dictionary with codons as keys and
		counts as values.  genetic_code is a dictionary with codons as keys and
		amino acids as values.  Defaults to standard_code.

		Post-condition: Creates various dictionaries needed to run 
		compute_results. Creates self.codon_expansions (which map degenerate
		codons to a list of the codons it expands to), self.amino_counts (which
		maps degenerate codons to a Counter of the amino acids is codes for), 
		and self.rev_amino_counts (which maps amino acid set to the degenerate 
		codons which code the amino acids in that set).
		"""
		
		self.codon_counts = codon_counts #See docstring for structure.
		self.total = sum(self.codon_counts.values())

		#self.genetic_code is a dictionary with codons as keys and the amino
		#acid it codes as values.
		self.genetic_code = genetic_code
		
		#self.codon_expansions has degenerate codons as keys and a list of all 
		#the codons it expands to as values.
		self.codon_expansions = {}
		
		#self.amnio_counts has degenerate codons as keys and counters as values.
		#The counters contain counts of the different amino acids the degenerate
		#codon can code for 
		self.amino_counts = {}

		#self.rev_amino_counts has an alphabetically sorted string of amino 
		#acids as keys and a list of the degenerate codons which code those 
		#amino acids as values.
		self.rev_amino_counts = {}

		self.alpha_2_base = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
				'R': 'GA', 'Y': 'TC', 'K':'GT', 'M': 'AC', 'S': 'GC', 'W': 'AT',
				'B': 'GTC', 'D': 'GAT', 'H': 'ACT', 'V': 'GCA', 'N': 'AGCT'}

		for dcodon in itertools.product(''
					.join(self.alpha_2_base.keys()), repeat=3):
			self.codon_expansions[''.join(dcodon)] = self.expand(dcodon)

		for dcodon in self.codon_expansions.keys():
			self.amino_counts[dcodon] = collections.Counter()
			for codon in self.codon_expansions[dcodon]:
				self.amino_counts[dcodon][self.genetic_code[codon]] += 1

		for dcodon in self.amino_counts.keys():
			amino_c = self.amino_counts[dcodon]
			try:
				self.rev_amino_counts[''
						.join(sorted(amino_c.keys()))].append(dcodon)
			except:
				self.rev_amino_counts[''
						.join(sorted(amino_c.keys()))] = [dcodon]

	def expand(self, dcodon):
		"""
		Pre-condition:  dcodon is a degenerate codon (a string of three letters
		using the enhanced DNA alphabet).

		Post-condition: returns a list of all codons that dcodon expands to.
		"""

		codon_list = []

		for codon in itertools.product(self.alpha_2_base[dcodon[0]],
				self.alpha_2_base[dcodon[1]], self.alpha_2_base[dcodon[2]]):
			codon_list.append(''.join(codon))

		return codon_list

	def imbalance(self, dcodon):
		"""
		Pre-condition: dcodon is a degenerate codon (a string of three letters
		using the enhanced DNA alphabet).
		
		Post-condition: Returns the imbalance of the degenerate codon.  This is
		the difference between count of the most common amino acid and least
		common amino acid (which the degenerate codon codes for) over the total
		number of amino acids which it codes for. 
		"""

		counts = self.amino_counts[dcodon].values()
		return (max(counts)-min(counts))/sum(counts)

	def freq(self, dcodon):
		"""
		Pre-condition: dcodon is a degenerate codon (a string of three letters
		using the enhanced DNA alphabet).
		
		Post-condition:  Returns the average frequency (the average of the 
		counts of the codons dcodon expands to over the total number of codons)
		"""

		ave_count = sum(self.codon_counts[codon] for codon in 
			self.codon_expansions[dcodon])/len(self.codon_expansions[dcodon])

		return ave_count/self.total

	def compute_results(self):
		"""
		Pre-condition: Class has been initialized.

		Post-condition: Computes the imbalance and frequency for each degenerate
		codon.  Sorts results for each amino acid set by smallest imbalance, 
		then by largest frequency (so the codon with the minimal imbalance and 
		highest frequency is first in this list).
		"""

		self.results = {}

		#self.results has amino acid sets as keys and lists of lists as values.
		#The inner list is of the form [degenerate codon, imbalance, frequency]
		#The outer list contains the lists for all the degenerate codon which
		#code the amino acid set (which is the key), sorted by lowest imbalance
		#then by highest frequency.

		for AAset in self.rev_amino_counts.keys():
			self.results[AAset] = []
			for dcodon in self.rev_amino_counts[AAset]:
				self.results[AAset].append([dcodon,self.imbalance(dcodon),
						self.freq(dcodon)])
			self.results[AAset].sort(key = lambda x: (x[1], -x[2]))

	def output(self, output_format):
		"""
		Pre-condition: compute_results method has been executed.  output_format
		is either 'min-codons', 'full', 'all-codons', or 'minimal'.

		Post-condition: 
		If output_format=='min-codons', 3 tab-separated fields are printed: the 
		amino acid set, the minimal imbalance, and all codons with minimal 
		imbalance that code for the amino acid set. 

		If output_format=='full', 2 space-separated fields are printed: the 
		amino acid set and a list of all each codon, its imbalance, and its
		frequency.  For example, for the amino-acid set KPQRT we get the line
		KPQRT ['MVG, 0.167, 0.0142', 'MVR, 0.167, 0.0131', 'MVA, 0.167, 0.0120']

		If output_format=='all-codons', 3 tab-separated fields are printed: the 
		amino acid set, the minimal imbalance, and all codons that code for 
		the amino acid set.

		Otherwise, 3 tab-separated fields are printed: the amino acid set, the
		minimal imbalance, and the codon with minimal imbalance and highest
		frequency. This is the 'minimal' format (default).

		For options which output multiple codons (all but minimal), codons are
		printed in the order stored in self.results (which is sorted first by 
		smallest imbalance, then by highest frequency).

		In all cases, the table is sorted alphabetically by amino-acid set.
		Imbalances are reported with 3 decimal places, frequencies with 4.
		"""

		if(output_format == 'min-codons'):
			for AA in sorted(self.results.keys()):
				min_imbalance = self.results[AA][0][1]
				min_codons = [x[0] for x in self.results[AA] 
													if x[1] == min_imbalance]
				print("{:<20}\t{:.3f}\t{}"
						.format(AA, min_imbalance, ','.join(min_codons)))
		elif(output_format == 'full'):
			for AA in sorted(self.results.keys()):
				print(AA,'['+', '.join(["'{}, {:.3f}, {:.4f}'"
					.format(x[0],x[1],x[2]) for x in self.results[AA]])+']')
		elif(output_format == 'all-codons'):
			for AA in sorted(self.results.keys()):
				min_imbalance = self.results[AA][0][1]
				all_codons = [x[0] for x in self.results[AA]]
				print("{:<20}\t{:.3f}\t{}"
						.format(AA, min_imbalance, ','.join(all_codons)))
		else:
			for AA in sorted(self.results.keys()):
				min_imbalance = self.results[AA][0][1]
				print("{:<20}\t{:.3f}\t{}"
						.format(AA, min_imbalance, self.results[AA][0][0]))

	def best_codon(self, AAset):
		"""
		Pre-condition: compute_results has been executed.  AAset is any string
		of amino acids and '*' (not necessarily a set that is coded).

		Post-condition:  Returns the codon with minimal imbalance and highest
		frequency that codes AAset.  If none exists, additional amino acids are
		added to the set until one is found.  If more than one is found when 
		adding amino acids, the one with the smallest imbalance and highest 
		frequency is picked.
		"""

		best_codon = None

		#num_extra is the number of the extra amino acids added to AAset.
		num_extra = 0

		#AAsupset is a superset of AAset.  Used to add extra amino acids if
		#AAset is not coded by any codons.

		#codon_list is a list of lists.  The inner are formated as
		#[degenerate codon, imbalance, frequency].  codon_list contains these 
		#for codons that code for AAsupset.
		codon_list = []
		
		alpha = 'ACDEFGHIKLMNPQRSTVWY*'

		while(not best_codon):
			for extra_aa in itertools.combinations(alpha, num_extra):
				AAsupset = ''.join(sorted(set(''.join([AAset,''.join(extra_aa)]))))
				try: codon_list.append(self.results[''.join(AAsupset)][0]) 
				except: pass
			if(codon_list):
				codon_list.sort(key=lambda x: (x[1], -x[2]))
				best_codon = codon_list[0][0]
			else: num_extra += 1

		return best_codon
