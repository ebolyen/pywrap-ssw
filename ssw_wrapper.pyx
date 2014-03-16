#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2014--, biocore team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from collections import namedtuple

import numpy as np
cimport numpy as np

np.import_array()

from cpython cimport bool
cdef extern from "ssw.h":

    ctypedef struct cigar:
        np.uint32_t* seq
        np.int32_t length

    ctypedef struct s_align:
        np.uint16_t score1
        np.uint16_t score2
        np.int32_t ref_begin1
        np.int32_t ref_end1
        np.int32_t read_begin1
        np.int32_t read_end1
        np.int32_t ref_end2
        np.uint32_t* cigar
        np.int32_t cigarLen

    ctypedef struct s_profile:
        const np.int8_t* read
        const np.int8_t* mat
        np.int32_t readLen
        np.int32_t n
        np.uint8_t bias

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

    cdef _get_cigar(self):
        for i in range(self.p.cigarLen):
            print self.p.cigar[i]

    def __str__(self):
        return "{\n\t'%s':%d\n"\
                  "\t'%s':%d\n"\
                  "\t'%s':%d\n"\
                  "\t'%s':%d\n"\
                  "\t'%s':%d\n"\
                  "\t'%s':%d\n"\
                  "\t'%s':%d\n"\
                  "\t'%s':%d\n}" % \
                   (
                    'score1',self.p.score1,
                    'score2',self.p.score2,
                    'ref_begin1',self.p.ref_begin1,
                    'ref_end1',self.p.ref_end1,
                    'read_begin1',self.p.read_begin1,
                    'read_end1',self.p.read_end1,
                    'ref_end2',self.p.read_end1,
                    'cigarLen',self.p.cigarLen
                    )

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
        if name == 'cigar_string':
            print self._get_cigar()
        if name == 'cigarLen':
            print self.p.cigarLen
        else:
            raise AttributeError(\
                "'AlignmentStructure' object has no attribute '%s'" % name)



cdef class StripedSmithWaterman:
    cdef s_profile *profile
    cdef np.ndarray __KEEP_IT_IN_SCOPE_read
    cdef np.ndarray __KEEP_IT_IN_SCOPE_matrix
    cdef bool is_protein
    cdef np.uint8_t weight_gap_open
    cdef np.uint8_t weight_gap_extension
    cdef np.uint8_t bit_flag
    cdef np.uint16_t score_filter
    cdef np.int32_t distance_filter
    cdef np.int32_t mask_length

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
                mismatch = 3
                if 'match' in kwargs and kwargs['match'] is not None:
                    match = kwargs['match']
                if 'mismatch' in kwargs and kwargs['mismatch'] is not None:
                    mismatch = kwargs['mismatch']
                matrix = self._build_match_matrix(match, mismatch)
        m_width = 16 if self.is_protein else 4
        read_seq = self._seq_converter(read_sequence, self.is_protein)
        read_length = len(read_sequence)
        cdef s_profile* p
        self.__KEEP_IT_IN_SCOPE_read = read_seq
        self.__KEEP_IT_IN_SCOPE_matrix = matrix
        p = ssw_init(<np.int8_t*> read_seq.data, 
                                read_length, 
                                <np.int8_t*> matrix.data, 
                                m_width, 
                                score_size)


        self.profile = p
        self.test()

    def test(self):
        print <int>self.profile
        print <int> self.profile.read
        print self.profile.readLen
        print '----'
        cdef np.int8_t* ar_ay
        for i in range(self.profile.readLen):
            ar_ay = self.profile.read
            print ar_ay[i]

    def __dealloc__(self):
        if self.profile is not NULL:
            print "dealloc profile"
            init_destroy(self.profile)

    def _get_bit_flag(self):
        return self.bit_flag


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

        if self.bit_flag is None:
            self.bit_flag = 128
        if 'bit_flag' in kwargs and kwargs['bit_flag'] is not None:
            self.bit_flag = kwargs['bit_flag']

    def __call__(self, reference_sequence, **kwargs):
        self.test()

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
        bit_flag = self._get_bit_flag()
        score_filter = 0 # self.score_filter
        distance_filter = 0 # self.distance_filter
        mask_length = 15 # self.mask_length

        cdef s_align *align
        align = ssw_align(self.profile, <np.int8_t*> reference.data, ref_length, 
                          weight_gap_open, weight_gap_extension, 
                          bit_flag, score_filter, distance_filter, 
                          mask_length)
        print align.score1

        alignment = AlignmentStructure()
        alignment.__constructor__(align)
        return alignment

        #return align

def striped_smith_waterman_alignment(seq1, seq2, **kwargs):
    ssw = StripedSmithWaterman(seq1, **kwargs)
    return ssw(seq2)
