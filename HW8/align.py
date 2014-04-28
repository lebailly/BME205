#!/usr/bin/env python2.7

#File: lebailly/BME205/HW8/align.py
#Author: Chris LeBailly

"""
align.py contains a global_aligner class and a local_aligner class.  Both of
these use an affine gap cost algorithm.  There is also a local_linear_aligner
class which uses a linear gap cost.  There is a function, a2m_to_fasta, 
which converts a2m alignments to fasta alignments. This is used by all classes.
There is also a function, extract_local(s1_fasta_list, s2_fasta_list), which
removes leading and trailing gaps.  This is used by the local aligners.
"""

from __future__ import division, print_function
import numpy

inf = 2**31-1

class global_aligner(object):
    """ 
    global_aligner finds a global alignments of two sequences, using an affine
    gap cost.  On construction requires a substitution matrix, as well as costs 
    for opening, extending, and double gaps.

    score_a2m(s1,s2) returns the score for the global alignment of s1 to s2.
    
    align(row_seq, col_seq) computes the M, IR, and IC matrices (stored as numpy
    arrays).  Returns the score of the best global alignment.

    traceback_col_seq() returns the column sequence (aligned to the row sequence)
    in A2M format.
    """

    def __init__(self, subst, open=12, extend=1, double=3):
        """
        subst is a dictionary with keys of tuples of two amino acids (such as
        ('M','A')) and values of the substitution value for matching those two
        amino acids.

        open, extend, and double are the gap costs for opening a gap, extending
        an existing gap, and switching a gap from one sequence to another,
        respectively.  Defaults set to open=12, extend=1, double=3.
        """

        self.subst = subst
        self.open = open
        self.extend = extend
        self.double = double

    def score_a2m(self, s1, s2):
        """
        Pre-condition: s1 and s2 represent an alignment in A2M format.

        Post-condition: Returns the score of the global alignment of s1 and
        s2 with the gap costs set in construction.
        """

        score = 0

        #s1_fasta_list is a list of the s1 sequence in fasta format
        #s2_fasta_list is a list of the s2 sequence in fasta format
        s1_fasta_list, s2_fasta_list = a2m_to_fasta(s1,s2)

        if(s1_fasta_list[0] == '-' or s2_fasta_list[0] == '-'): score-=self.open
        else: score += self.subst[(s1_fasta_list[0],s2_fasta_list[0])]

        for i in xrange(1,len(s1_fasta_list)):
            if((s1_fasta_list[i] == '-' and s2_fasta_list[i-1] == '-') or 
                    (s2_fasta_list[i] == '-' and s1_fasta_list[i-1] == '-')):
                score -= self.double           
            elif((s1_fasta_list[i] == '-' and s1_fasta_list[i-1] != '-') or
                    (s2_fasta_list[i] == '-' and s2_fasta_list[i-1] != '-')): 
                score -= self.open
            elif((s1_fasta_list[i] == '-' and s1_fasta_list[i-1] == '-') or
                    (s2_fasta_list[i] == '-' and s2_fasta_list[i-1] == '-')):
                score -= self.extend
            else:
                score += self.subst[(s1_fasta_list[i],s2_fasta_list[i])]

        return score

    def align(self, row_seq, col_seq):
        """
        Pre-condition: row_seq and col_seq are sequences using the same alphabet
        as self.subst.

        Post-condition: Creates and fills in alignment matrices self.M, self.IR,
        and self.IC (which are numpy arrays with index zero corresponding to
        the start character) and stores row_seq and col_seq.  Returns the score 
        of the best global alignment.
        """

        #self.M is numpy array. self.M[r][c] is the score of the best alignment 
        #of row_seq[0]...row_seq[r-1] to col_seq[0]...col_seq[c-1] ending with
        #row_seq[r-1] aligned to col_seq[c-1] (for r,c > 0). An index of zero 
        #corresponds to the start character.

        #self.IR is numpy array. self.IR[r][c] is the score of the best alignment 
        #of row_seq[0]...row_seq[r-1] to col_seq[0]...col_seq[c-1] ending with a
        #gap aligned to col_seq[c-1] (for r,c > 0). An index of zero corresponds
        #to the start character.

        #self.IC is numpy array. self.IC[r][c] is score of the best alignment of 
        #row_seq[0]...row_seq[r-1] to col_seq[0]...col_seq[c-1] ending with
        #row_seq[r-1] aligned to a gap (for r,c > 0). An index of zero 
        #corresponds to the start character.

        self.M = numpy.zeros([len(row_seq)+1, len(col_seq)+1], dtype=int)
        self.IR = numpy.zeros([len(row_seq)+1, len(col_seq)+1], dtype=int)
        self.IC = numpy.zeros([len(row_seq)+1, len(col_seq)+1], dtype=int)

        for index, char in enumerate(row_seq):
            self.M[index+1][0] = -inf
            self.IR[index+1][0] = int(-self.open - index*self.extend)
            self.IC[index+1][0] = -inf
        for index, char in enumerate(col_seq):
            self.M[0][index+1] = -inf
            self.IR[0][index+1] = -inf
            self.IC[0][index+1] = int(-self.open - index*self.extend)

        for row_index, row_char in enumerate(row_seq):
            for col_index, col_char in enumerate(col_seq):
                self.M[row_index+1][col_index+1] = self.subst[(row_char,
                        col_char)] + max(self.M[row_index][col_index],
                        self.IR[row_index][col_index],
                        self.IC[row_index][col_index])
                self.IR[row_index+1][col_index+1] = max(
                        self.M[row_index][col_index+1] - self.open,
                        self.IR[row_index][col_index+1] - self.extend,
                        self.IC[row_index][col_index+1] - self.double)
                self.IC[row_index+1][col_index+1] = max(
                        self.M[row_index+1][col_index] - self.open,
                        self.IR[row_index+1][col_index] - self.double,
                        self.IC[row_index+1][col_index] - self.extend)

        self.row_seq = row_seq
        self.col_seq = col_seq

        return max(self.M[row_index + 1][col_index+1],
                    self.IR[row_index+1][col_index+1],
                    self.IC[row_index+1][col_index+1])


    def traceback_col_seq(self):
        """
        Pre-condition: .align method has been executed.

        Post-condition: Returns the column sequence aligned to the row sequence
        (specified when .align was executed) in A2M format.
        """

        #rev_alignment_list is a list of the alignment of self.col_seq to 
        #self.row_seq in reverse order and in A2M format.

        rev_alignment_list = []

        #row, col is the current position of the traceback in the matrices 
        #self.M, self.IR, and self.IC.

        row = len(self.row_seq)
        col = len(self.col_seq)

        #best_score is the largest of in self.M[row][col], self.IR[row][col],
        #and self.IC[row][col]. If more than one max exists, the first occurrence 
        #of the max in self.M[row][col], self.IR[row][col], & self.IC[row][col] 
        #is stored.

        #best_type indicates which matrix best_score came from.  If best_score
        #came from self.M, best_type = 0.  If best_score came from self.IR, 
        #best_type = 1.  If best_score came from self.IC, best_type = 2.

        best_type, best_score = max(enumerate([self.M[row][col],
                self.IR[row][col], self.IC[row][col]]), key=lambda x:x[1])

        while(row != 0 or col != 0):

            if best_type == 0:
                rev_alignment_list.append(self.col_seq[col-1].upper())
                subst = self.subst[(self.row_seq[row-1], self.col_seq[col-1])]

                row -= 1
                col -= 1

                if(best_score == subst + self.M[row][col]):
                    best_type = 0
                    best_score = self.M[row][col]
                elif(best_score == subst + self.IR[row][col]):
                    best_type = 1
                    best_score = self.IR[row][col]
                elif(best_score == subst + self.IC[row][col]):
                    best_type = 2
                    best_score = self.IC[row][col]

            elif best_type == 1:
                rev_alignment_list.append('-')

                row -= 1

                if(best_score + self.open == self.M[row][col]): 
                    best_type = 0
                    best_score = self.M[row][col]
                elif(best_score + self.extend == self.IR[row][col]):
                    best_type = 1
                    best_score = self.IR[row][col]
                elif(best_score + self.double == self.IC[row][col]):
                    best_type = 2
                    best_score = self.IC[row][col]

            elif best_type == 2:
                rev_alignment_list.append(self.col_seq[col-1].lower())

                col -= 1

                if(best_score + self.open == self.M[row][col]): 
                    best_type = 0
                    best_score = self.M[row][col]
                elif(best_score + self.double == self.IR[row][col]):
                    best_type = 1
                    best_score = self.IR[row][col]
                elif(best_score + self.extend == self.IC[row][col]):
                    best_type = 2
                    best_score = self.IC[row][col]

        return ''.join(rev_alignment_list[::-1])

class local_aligner(object):
    """
    local_aligner finds a local alignment of two sequences, using an affine
    gap cost.  On construction requires a substitution matrix, as well as costs 
    for opening, extending, and double gaps.

    score_a2m(s1,s2) returns the score for the local alignment of s1 to s2.
    
    align(row_seq, col_seq) computes the M, IR, and IC matrices (stored as numpy
    arrays).  Returns the score of the best local alignment.

    traceback_col_seq() returns the column sequence (aligned to the row sequence)
    in A2M format.
    """

    def __init__(self, subst, open=12, extend=1, double=3):
        """
        subst is a dictionary with keys of tuples of two amino acids (such as
        ('M','A')) and values of the substitution value for matching those two
        amino acids.
        
        open, extend, and double are the gap costs for opening a gap, extending
        an existing gap, and switching a gap from one sequence to another,
        respectively.  Defaults set to open=12, extend=1, double=3.
        """

        self.subst = subst
        self.open = open
        self.extend = extend
        self.double = double

    def score_a2m(self, s1, s2):
        """
        Pre-condition: s1 and s2 represent an alignment in A2M format.

        Post-condition: Returns the score of the local alignment of s1 and 
        s2 with the gap costs set in construction.
        """

        score = 0

        #s1_fasta_list is a list of the s1 sequence in fasta format
        #s2_fasta_list is a list of the s2 sequence in fasta format
        s1_fasta_list, s2_fasta_list = a2m_to_fasta(s1,s2)

        extract_local(s1_fasta_list, s2_fasta_list)
        
        if(len(s1_fasta_list) != 0 and len(s2_fasta_list) != 0):
            score += self.subst[(s1_fasta_list[0],s2_fasta_list[0])]

        for i in xrange(1,len(s1_fasta_list)):
            if((s1_fasta_list[i] == '-' and s2_fasta_list[i-1] == '-') or 
                    (s2_fasta_list[i] == '-' and s1_fasta_list[i-1] == '-')):
                score -= self.double           
            elif((s1_fasta_list[i] == '-' and s1_fasta_list[i-1] != '-') or
                    (s2_fasta_list[i] == '-' and s2_fasta_list[i-1] != '-')): 
                score -= self.open
            elif((s1_fasta_list[i] == '-' and s1_fasta_list[i-1] == '-') or
                    (s2_fasta_list[i] == '-' and s2_fasta_list[i-1] == '-')):
                score -= self.extend
            else:
                score += self.subst[(s1_fasta_list[i],s2_fasta_list[i])]

        return score

    def align(self, row_seq, col_seq):
        """
        Pre-condition: row_seq and col_seq are sequences using the same alphabet
        as self.subst.

        Post-condition: Creates and fills in alignment matrices self.M, self.IR,
        and self.IC (which are numpy arrays with index zero corresponding to
        the start character) and stores row_seq and col_seq.  Returns the score 
        of the best local alignment.
        """

        #self.M is numpy array. self.M[r][c] is the score of the best alignment 
        #of row_seq[0]...row_seq[r-1] to col_seq[0]...col_seq[c-1] ending with
        #row_seq[r-1] aligned to col_seq[c-1] (for r,c > 0). An index of zero 
        #corresponds to the start character.

        #self.IR is numpy array. self.IR[r][c] is the score of the best alignment 
        #of row_seq[0]...row_seq[r-1] to col_seq[0]...col_seq[c-1] ending with a
        #gap aligned to col_seq[c-1] (for r,c > 0). An index of zero corresponds
        #to the start character.

        #self.IC is numpy array. self.IC[r][c] is score of the best alignment of 
        #row_seq[0]...row_seq[r-1] to col_seq[0]...col_seq[c-1] ending with
        #row_seq[r-1] aligned to a gap (for r,c > 0). An index of zero corresponds
        #to the start character.

        self.M = numpy.zeros([len(row_seq)+1, len(col_seq)+1], dtype=int)
        self.IR = numpy.zeros([len(row_seq)+1, len(col_seq)+1], dtype=int)
        self.IC = numpy.zeros([len(row_seq)+1, len(col_seq)+1], dtype=int)

        #M[self.best_row][self.best_col] is the maximum value contained in 
        #self.M. If multiple maxes are achieved, first seen while creating 
        #self.M is used.
        
        self.best_row, self.best_col = 0,0

        for index, char in enumerate(row_seq):
            self.M[index+1][0] = -inf
            self.IR[index+1][0] = -inf
            self.IC[index+1][0] = -inf
        for index, char in enumerate(col_seq):
            self.M[0][index+1] = -inf
            self.IR[0][index+1] = -inf
            self.IC[0][index+1] = -inf

        for row_index, row_char in enumerate(row_seq):
            for col_index, col_char in enumerate(col_seq):
                self.M[row_index+1][col_index+1] = self.subst[(row_char,
                        col_char)] + max(0, self.M[row_index][col_index],
                        self.IR[row_index][col_index],
                        self.IC[row_index][col_index])
                self.IR[row_index+1][col_index+1] = max(
                        self.M[row_index][col_index+1] - self.open,
                        self.IR[row_index][col_index+1] - self.extend,
                        self.IC[row_index][col_index+1] - self.double)
                self.IC[row_index+1][col_index+1] = max(
                        self.M[row_index+1][col_index] - self.open,
                        self.IR[row_index+1][col_index] - self.double,
                        self.IC[row_index+1][col_index] - self.extend)

                if(self.M[row_index+1][col_index+1] > 
                        self.M[self.best_row][self.best_col]):
                    self.best_row = row_index + 1
                    self.best_col = col_index + 1

        self.row_seq = row_seq
        self.col_seq = col_seq

        return self.M[self.best_row][self.best_col]


    def traceback_col_seq(self):
        """
        Pre-condition: .align method has been executed.

        Post-condition: Returns the column sequence aligned to the row sequence
        (specified when .align was executed) in A2M format.
        """

        #row, col is the current position of the traceback in the matrices 
        #self.M, self.IR, and self.IC.

        row = len(self.row_seq)
        col = len(self.col_seq)

        rev_alignment_list = []

        while(row > self.best_row):
            rev_alignment_list.append('-')
            row -= 1
        while(col > self.best_col):
            rev_alignment_list.append(self.col_seq[col-1].lower())
            col -= 1

        #best_score is the largest of in self.M[row][col], self.IR[row][col],
        #and self.IC[row][col]. If more than one max exists, the first occurrence 
        #of the max in self.M[row][col], self.IR[row][col], & self.IC[row][col] 
        #is stored.

        #best_type indicates which matrix best_score came from.  If best_score
        #came from self.M, best_type = 0.  If best_score came from self.IR, 
        #best_type = 1.  If best_score came from self.IC, best_type = 2.

        #end_local_alignment is set to True when the traceback of the local
        #alignment is over (and is False while the traceback is running) 

        best_type = 0
        best_score = self.M[self.best_row][self.best_col]
        end_local_alignment = False

        while((row != 0 or col != 0) and not end_local_alignment):

            if best_type == 0:
                rev_alignment_list.append(self.col_seq[col-1].upper())
                subst = self.subst[(self.row_seq[row-1], self.col_seq[col-1])]

                row -= 1
                col -= 1

                if(best_score == subst):
                    end_local_alignment = True
                elif(best_score == subst + self.M[row][col]):
                    best_type = 0
                    best_score = self.M[row][col]
                elif(best_score == subst + self.IR[row][col]):
                    best_type = 1
                    best_score = self.IR[row][col]
                elif(best_score == subst + self.IC[row][col]):
                    best_type = 2
                    best_score = self.IC[row][col]

            elif best_type == 1:
                rev_alignment_list.append('-')

                row -= 1

                if(best_score + self.open == self.M[row][col]): 
                    best_type = 0
                    best_score = self.M[row][col]
                elif(best_score + self.extend == self.IR[row][col]):
                    best_type = 1
                    best_score = self.IR[row][col]
                elif(best_score + self.double == self.IC[row][col]):
                    best_type = 2
                    best_score = self.IC[row][col]

            elif best_type == 2:
                rev_alignment_list.append(self.col_seq[col-1].lower())

                col -= 1

                if(best_score + self.open == self.M[row][col]): 
                    best_type = 0
                    best_score = self.M[row][col]
                elif(best_score + self.double == self.IR[row][col]):
                    best_type = 1
                    best_score = self.IR[row][col]
                elif(best_score + self.extend == self.IC[row][col]):
                    best_type = 2
                    best_score = self.IC[row][col]

        while(col > 0):
            rev_alignment_list.append(self.col_seq[col-1].lower())
            col -= 1
        while(row > 0):
            rev_alignment_list.append('-')
            row -= 1

        return ''.join(rev_alignment_list[::-1])

class local_linear_aligner(object):
    """
    local_linear_aligner finds a local alignment of two sequences using
    a linear gap cost.  On construction requires a substitution matrix, as well
    as the linear gap cost.

    score_a2m(s1,s2) returns the score for the local alignment of s1 to s2.
    
    align(row_seq, col_seq) computes the B matrix.  Returns the score of the 
    best local alignment.

    traceback_col_seq() returns the column sequence (aligned to the row sequence)
    in A2M format.
    """

    def __init__(self, subst, gap=3):
        """
        subst is a dictionary with keys of tuples of two amino acids (such as
        ('M','A')) and values of the substitution value for matching those two
        amino acids.
        
        gap is the linear gap costs for any gap. Default is 3.
        """

        self.subst = subst
        self.gap = gap

    def score_a2m(self, s1, s2):
        """
        Pre-condition: s1 and s2 represent an alignment in A2M format.

        Post-condition: Returns the score of the local alignment of s1 and
        s2 with the gap costs set in construction.
        """

        score = 0

        #s1_fasta_list is a list of the s1 sequence in fasta format
        #s2_fasta_list is a list of the s2 sequence in fasta format
        s1_fasta_list, s2_fasta_list = a2m_to_fasta(s1,s2)

        extract_local(s1_fasta_list, s2_fasta_list)

        for i in xrange(len(s1_fasta_list)):
            if(s1_fasta_list[i] == '-'  or s2_fasta_list[i] == '-'):
                score -= self.gap
            else:
                score += self.subst[(s1_fasta_list[i],s2_fasta_list[i])]

        return score

    def align(self, row_seq, col_seq):
        """
        Pre-condition: row_seq and col_seq are sequences using the same alphabet
        as self.subst.

        Post-condition: Creates and fills in alignment matrix self.B (which is
        a numpy arrays with index zero corresponding to the start character) and 
        stores row_seq and col_seq. Returns the score of the best local alignment.
        """

        #self.B is numpy array. self.B[r][c] is the score of the best alignment 
        #of row_seq[0]...row_seq[r-1] to col_seq[0]...col_seq[c-1] ending with
        #row_seq[r-1] aligned to col_seq[c-1] (for r,c > 0). An index of zero 
        #corresponds to the start character.

        self.B = numpy.zeros([len(row_seq)+1, len(col_seq)+1], dtype=int)

        #M[self.best_row][self.best_col] is the maximum value contained in 
        #self.B. If multiple maxes are achieved, first seen while creating 
        #self.B is used.

        self.best_row, self.best_col = 0,0

        for index, char in enumerate(row_seq):
            self.B[index+1][0] = -inf
        for index, char in enumerate(col_seq):
            self.B[0][index+1] = -inf

        for row_index, row_char in enumerate(row_seq):
            for col_index, col_char in enumerate(col_seq):
                subst_value = self.subst[(row_char,col_char)]
                self.B[row_index+1][col_index+1] = max(
                        subst_value, 
                        subst_value + self.B[row_index][col_index],
                        self.B[row_index][col_index+1] - self.gap,
                        self.B[row_index+1][col_index] - self.gap)

                if(self.B[row_index+1][col_index+1] > 
                        self.B[self.best_row][self.best_col]):
                    self.best_row = row_index + 1
                    self.best_col = col_index + 1

        self.row_seq = row_seq
        self.col_seq = col_seq

        return self.B[self.best_row][self.best_col]


    def traceback_col_seq(self):
        """
        Pre-condition: .align method has been executed.
        Post-condition: Returns the column sequence aligned to the row sequence
        (specified when .align was executed) in A2M format.
        """

        rev_alignment_list = []

        #row, col is the current position of the traceback in the matrix self.B.

        row = len(self.row_seq)
        col = len(self.col_seq)

        while(row > self.best_row):
            rev_alignment_list.append('-')
            row -= 1
        while(col > self.best_col):
            rev_alignment_list.append(self.col_seq[col-1].lower())
            col -= 1

        #end_local_alignment is set to True when the traceback of the local
        #alignment is over (and is False while the traceback is running) 

        end_local_alignment = False

        while((row != 0 or col != 0) and not end_local_alignment):
            score = self.B[row][col]
            subst = self.subst[(self.row_seq[row-1],self.col_seq[col-1])]

            if(score == self.B[row-1][col-1] + subst):
                rev_alignment_list.append(self.col_seq[col-1].upper())
                row -= 1
                col -= 1
            elif(score == self.B[row-1][col] - self.gap):
                rev_alignment_list.append('-')
                row -= 1
            elif(score == self.B[row][col-1] - self.gap):
                rev_alignment_list.append(self.col_seq[col-1].lower())
                col -= 1
            elif(score == subst):
                rev_alignment_list.append(self.col_seq[col-1].upper())
                col -= 1
                row -= 1
                end_local_alignment = True

        while(col > 0):
            rev_alignment_list.append(self.col_seq[col-1].lower())
            col -= 1
        while(row > 0):
            rev_alignment_list.append('-')
            row -= 1

        return ''.join(rev_alignment_list[::-1])

def a2m_to_fasta(s1,s2):
    """
    Pre-condition: s1 and s2 aligned sequences in A2M format.

    Post-condition: Returns two lists of the s1 and s2 sequences in FASTA format.
    """

    #s1_fasta_list is a list of the s1 sequence in fasta format
    #s2_fasta_list is a list of the s2 sequence in fasta format 
    s1_fasta_list, s2_fasta_list = [], []
    s1_index, s2_index = 0, 0

    while(s1_index < len(s1) and s2_index < len(s2)):

        #Characters in A2M sequence align
        if((s1[s1_index] == s1[s1_index].upper() or s1[s1_index] == '-') and
           (s2[s2_index] == s2[s2_index].upper() or s2[s2_index] == '-')):
           
           s1_fasta_list.append(s1[s1_index])
           s2_fasta_list.append(s2[s2_index])

           s1_index += 1
           s2_index += 1

        #Gap in s2
        elif(s1[s1_index] == s1[s1_index].lower()):
            s1_fasta_list.append(s1[s1_index].upper())
            s2_fasta_list.append('-')
            s1_index += 1

        #Gap in s1
        elif(s2[s2_index] == s2[s2_index].lower()):
            s1_fasta_list.append('-')
            s2_fasta_list.append(s2[s2_index].upper())
            s2_index += 1

    while(s1_index < len(s1)):
        s1_fasta_list.append(s1[s1_index].upper())
        s2_fasta_list.append('-')
        s1_index += 1

    while(s2_index < len(s2)):
        s1_fasta_list.append('-')            
        s2_fasta_list.append(s2[s2_index].upper())
        s2_index += 1

    return s1_fasta_list, s2_fasta_list

def extract_local(s1_fasta_list, s2_fasta_list):
    """
    s1_fasta_list and s2_fasta_list are lists of two sequences of equal
    length in fasta format.  extract_local removes leading and trailing
    gaps, so s1_fasta_list and s2_fasta_list start and end with aligned
    characters.
    """

    while((s1_fasta_list[0] == '-' if s1_fasta_list else False)
            or (s2_fasta_list[0] == '-' if s2_fasta_list else False)):
        del s1_fasta_list[0]
        del s2_fasta_list[0]

    while((s1_fasta_list[-1] == '-' if s1_fasta_list else False) or 
            (s2_fasta_list[-1] == '-' if s2_fasta_list else False)):
        del s1_fasta_list[-1]
        del s2_fasta_list[-1]
