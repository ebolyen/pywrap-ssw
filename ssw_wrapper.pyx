#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2014--, biocore team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------


import numpy as np
cimport numpy as np

np.import_array()

from cpython cimport bool
cdef extern from "ssw.h":

    ctypedef struct cigar:
        unsigned int* seq
        int length

    ctypedef struct s_align:
        unsigned int score1
        unsigned int score2
        int ref_begin1
        int ref_end1
        int read_begin1
        int read_end1
        int ref_end2
        int* cigar
        int cigarLen

    ctypedef struct s_profile:
        pass

    cdef s_profile* ssw_init(const np.int8_t* read, 
                             const np.int32_t readLen, 
                             const np.int8_t* mat, #<THIS FUCKER
                             const np.int32_t n, 
                             const np.int8_t score_size)

    cdef void init_destroy(s_profile* p)

    cdef s_align* ssw_align(const s_profile* prof, 
                            const np.int8_t* ref, 
                            np.int32_t refLen, 
                            const np.uint8_t weight_gapO, 
                            const np.uint8_t weight_gapE, 
                            const np.uint8_t flag, 
                            const np.uint16_t filters,
                            const np.int32_t filterd,
                            const np.int32_t maskLen)

    cdef void align_destroy(s_align* a)

np_aa_table = np.array([
        23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
        23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
        23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
        23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
        23,  0, 20,  4,  3,  6, 13,  7,  8,  9, 23, 11, 10, 12,  2, 23,
        14,  5,  1, 15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23,
        23,  0, 20,  4,  3,  6, 13,  7,  8,  9, 23, 11, 10, 12,  2, 23,
        14,  5,  1, 15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23])

np_nt_table = np.array([
         4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
         4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
         4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
         4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
         4,  0,  4,  1,  4,  4,  4,  2,  4,  4,  4,  4,  4,  4,  4,  4,
         4,  4,  4,  4,  3,  0,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
         4,  0,  4,  1,  4,  4,  4,  2,  4,  4,  4,  4,  4,  4,  4,  4,
         4,  4,  4,  4,  3,  0,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4])


cdef class AlignmentStructure:
    cdef s_align *p

    cdef __constructor__(self, s_align* pointer):
        self.p = pointer
    
    def __dealloc__(self):
        if self.p is not NULL:
            print "dealloc alignment"
            align_destroy(self.p)

    def __getattr__(self, name):
        if name == 'score1':
            return self.p.score1
        if name == 'score2':
            return self.p.score2
        if name == 'ref_begin1':
            return self.p.ref_begin1
        if name == 'ref_end1':
            return self.p.ref_end1
        if name == 'read_begin1':
            return self.p.read_begin1
        if name == 'read_end1':
            return self.p.read_end1
        if name == 'ref_end2':
            return self.p.ref_end2
        if name == 'sequence':
            pass #TODO

cdef class StripedSmithWaterman:
    cdef s_profile *profile
    cdef bool is_protein
    cdef int weight_gap_open
    cdef int weight_gap_extension


    def __init__(self, read_sequence, **kwargs):
        cdef np.int8_t score_size
        score_size = 2
        cdef np.ndarray[np.int8_t, ndim=1, mode="c"] matrix
        cdef np.int32_t read_length
        cdef np.int32_t m_width
        cdef np.ndarray[np.int8_t, ndim=1, mode="c"] read_seq

        if 'score_size' in kwargs and kwargs['score_size'] is not None:
            score_size = kwargs['score_size']
        if 'protein' in kwargs and kwargs['protein'] == True:
            self.is_protein = True
            matrix = self._convert_dict2d_to_matrix(kwargs['substitution_matrix'])
        else:
            self.is_protein = False
            if 'substitution_matrix' in kwargs and kwargs['substitution_matrix'] is not None:
                matrix = self._convert_dict2d_to_matrix(kwargs['substitution_matrix'])
            else:
                #BLASTN defaults
                match = 2
                mismatch = -3
                if 'match' in kwargs and kwargs['match'] is not None:
                    match = kwargs['match']
                if 'mismatch' in kwargs and kwargs['mismatch'] is not None:
                    mismatch = kwargs['mismatch']
                matrix = self._build_match_matrix(match, mismatch)
        m_width = 576 if self.is_protein else 16
        read_seq = self._seq_converter(read_sequence, self.is_protein)
        read_length = len(read_sequence)
        self.profile = ssw_init(<np.int8_t*> read_seq.data, 
                                read_length, 
                                <np.int8_t*> matrix.data, 
                                m_width, 
                                score_size)

    def __dealloc__(self):
        if self.profile is not NULL:
            print "dealloc profile"
            init_destroy(self.profile)

    def _get_bit_flag(self):
        pass


    cdef np.ndarray[np.int8_t, ndim=1, mode="c"] _seq_converter(self, sequence, is_protein):
        cdef np.ndarray[np.int8_t, ndim=1, mode="c"] seq = np.empty(len(sequence), dtype=np.int8)
        if is_protein:
            for i, char in enumerate(sequence):
                seq[i] = np_aa_table[ord(char)]
        else:
            for i, char in enumerate(sequence):
                seq[i] = np_nt_table[ord(char)]
        print seq
        return seq

    cdef np.ndarray[np.int8_t, ndim=1, mode="c"] _build_match_matrix(self, match, mismatch):
        sequence_order = "ACGT"
        dict2d = {}
        for row in sequence_order:
            dict2d[row] = {}
            for column in sequence_order:
                dict2d[row][column] = match if row == column else mismatch
        return self._convert_dict2d_to_matrix(dict2d)

    cdef np.ndarray[np.int8_t, ndim=1, mode="c"] _convert_dict2d_to_matrix(self, dict2d):
        if self.is_protein:
            sequence_order = "ARNDCQEGHILKMFPSTWYVBZX*"                  
        else:
            sequence_order = "ACGT"
        i = 0
        length = len(sequence_order)
        cdef np.ndarray[np.int8_t, ndim=1, mode="c"] py_list_matrix = np.empty(length*length, dtype=np.int8)
        print py_list_matrix.flag['']
        for row in sequence_order:
            for column in sequence_order:
                py_list_matrix[i] = dict2d[row][column]
                i+=1
        print py_list_matrix
        return py_list_matrix

    def _handle_shared_kwargs(self, **kwargs):
        if self.weight_gap_open is None:
            #From BLASTN defaults
            self.weight_gap_open = 5
        if 'weight_gap_open' in kwargs and kwargs['weight_gap_open'] is not None:
            self.weight_gap_open = kwargs['weight_gap_open']

        if self.weight_gap_extension is None:
            #from BLASTN defaults
            self.weight_gap_extension = 2
        if 'weight_gap_extension' in kwargs and kwargs['weight_gap_extension'] is not None:
            self.weight_gap_extension = kwargs['weight_gap_extension']

    def __call__(self, reference_sequence, **kwargs):
        profile = self.profile
        cdef np.ndarray[np.int8_t, ndim=1, mode="c"] reference
        reference = self._seq_converter(reference_sequence, self.is_protein)
        cdef np.int32_t ref_length
        ref_length = len(reference_sequence)
        cdef np.uint8_t weight_gap_open
        cdef np.uint8_t weight_gap_extension
        cdef np.uint8_t bit_flag
        cdef np.uint16_t score_filter
        cdef np.int32_t distance_filter
        cdef np.int32_t mask_length

        if kwargs is not None:
            self._handle_shared_kwargs(**kwargs)

        weight_gap_open = self.weight_gap_open
        weight_gap_extension = self.weight_gap_extension
        bit_flag = 1 # self._get_bit_flag()
        score_filter = 0 # self.score_filter
        distance_filter = 0 # self.distance_filter
        mask_length = 15 # self.mask_length

        cdef s_align *align
        align = ssw_align(profile, <np.int8_t*> reference.data, ref_length, 
                          weight_gap_open, weight_gap_extension, 
                          bit_flag, score_filter, distance_filter, 
                          mask_length)
        
        alignment = AlignmentStructure()
        alignment.__constructor__(align)
        return alignment

        #return align

def striped_smith_waterman_alignment(seq1, seq2, **kwargs):
    ssw = StripedSmithWaterman(seq1, **kwargs)
    return ssw(seq2)
