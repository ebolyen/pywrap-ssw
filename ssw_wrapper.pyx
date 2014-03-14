#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2014--, biocore team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------


cimport numpy as np
import numpy as np
from cpython cimport bool
cdef extern from "ssw.h":

    ctypedef struct cigar:
        np.uint32_t *seq
        np.int32_t length

    ctypedef struct s_align:
        np.uint16_t score1
        np.uint16_t score2
        np.int32_t ref_begin1
        np.int32_t ref_end1
        np.int32_t read_begin1
        np.int32_t read_end1
        np.int32_t ref_end2
        np.uint32_t *cigar
        np.int32_t cigarLen

    ctypedef struct s_profile:
        pass

    cdef s_profile* ssw_init(const np.int8_t* read, 
                             const np.int32_t readLen, 
                             const np.int8_t* mat, 
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
        score_size = 2
        matrix = None
        m_width = None

        if 'score_size' in kwargs and kwargs['score_size'] is not None:
            score_size = kwargs['score_size']
        if 'protein' in kwargs and kwargs['protein'] == True:
            self.is_protein = True
            matrix, m_width = self._convert_dict2d_to_matrix(kwargs['substitution_matrix'])
        else:
            self.is_protein = False
            if 'substitution_matrix' in kwargs and kwargs['substitution_matrix'] is not None:
                matrix, m_width = self._convert_dict2d_to_matrix(kwargs['substitution_matrix'])
            else:
                #BLASTN defaults
                match = 2
                mismatch = -3
                if 'match' in kwargs and kwargs['match'] is not None:
                    match = kwargs['match']
                if 'mismatch' in kwargs and kwargs['mismatch'] is not None:
                    mismatch = kwargs['mismatch']
                matrix, m_width = self._build_match_matrix(match, mismatch)

        read_seq, read_length = self._seq_converter(read_sequence, self.is_protein)
        self.profile = ssw_init(read_sequence, read_length, 
                                matrix, m_width, score_size)

    def __dealloc__(self):
        if self.profile is not NULL:
            print "dealloc profile"
            init_destroy(self.profile)

    def _get_bit_flag(self):
        pass

    def _seq_converter(self, sequence, is_protein):
        seq_len = len(sequence)
        seq = np.empty(seq_len, dtype=np.int8)
        if is_protein:
            for i, char in enumerate(sequence):
                seq[i] = np_aa_table[ord(char)]
        else:
            for i, char in enumerate(sequence):
                seq[i] = np_nt_table[ord(char)]
        return (seq, seq_len)

    def _build_match_matrix(self, match, mismatch):
        sequence_order = "ACGT"
        dict2d = {}
        for row in sequence_order:
            dict2d[row] = {}
            for column in sequence_order:
                dict2d[row][column] = match if row == column else mismatch
        return self._convert_dict2d_to_matrix(dict2d)

    def _convert_dict2d_to_matrix(self, dict2d):
        py_list_matrix = ''
        if self.is_protein:
            sequence_order = "ARNDCQEGHILKMFPSTWYVBZX*"                  
        else:
            sequence_order = "ACGT"
        for row in sequence_order:
            for column in sequence_order:
                py_list_matrix += str(dict2d[row][column])
        return (py_list_matrix, len(py_list_matrix))

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
        reference = reference_sequence
        ref_length = len(reference_sequence)
        if kwargs is not None:
            self._handle_shared_kwargs(**kwargs)

        weight_gap_open = self.weight_gap_open
        weight_gap_extension = self.weight_gap_extension
        bit_flag = 1 # self._get_bit_flag()
        score_filter = 0 # self.score_filter
        distance_filter = 0 # self.distance_filter
        mask_length = 15 # self.mask_length

        cdef s_align *align
        align = ssw_align(profile, reference_sequence, ref_length, 
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
