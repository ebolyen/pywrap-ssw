#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2014--, biocore team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from cpython cimport bool
import numpy as np
cimport numpy as np

cdef extern from "ssw_native.h":

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

mid_table = np.array(['M','I','D'])

cdef class AlignmentStructure:
    cdef s_align *p
    cdef str read_sequence
    cdef str reference_sequence
    cdef int index_starts_at

    def __cinit__(self, read_sequence, reference_sequence, index_starts_at):
        self.read_sequence = read_sequence
        self.reference_sequence = reference_sequence
        self.index_starts_at = index_starts_at

    cdef __constructor__(self, s_align* pointer):
        self.p = pointer
    
    def __dealloc__(self):
        if self.p is not NULL:
            align_destroy(self.p)

    cdef _get_cigar(self):
        cigar_string = ""
        for i in range(self.p.cigarLen):
            cigar_string += str(self.p.cigar[i]>>4)
            cigar_string += mid_table[self.p.cigar[i]&0xf]
        return cigar_string

    def __str__(self):
        return "{\n    '%s': %d,\n"\
                  "    '%s': %d,\n"\
                  "    '%s': %d,\n"\
                  "    '%s': %d,\n"\
                  "    '%s': %d,\n"\
                  "    '%s': %d,\n"\
                  "    '%s': %d,\n"\
                  "    '%s': '%s',\n"\
                  "    '%s': '%s',\n"\
                  "    '%s': '%s'\n}" % \
                   (
                    'optimal_alignment_score',self.optimal_alignment_score,
                    'suboptimal_alignment_score',self.suboptimal_alignment_score,
                    'target_begin',self.target_begin1,
                    'target_end_optimal',self.target_end_optimal,
                    'target_end_suboptimal',self.target_end_suboptimal,
                    'query_begin',self.query_begin1,
                    'query_end',self.query_end1,
                    'cigar',self.cigar,
                    'query_sequence', self.query_sequence,
                    'target_sequence', self.target_sequence
                    )

    def __getattr__(self, name):
        index_starts_at =  self.index_starts_at

        if name == 'optimal_alignment_score':
            return self.p.score1
        if name == 'suboptimal_alignment_score':
            return self.p.score2
        if name == 'target_begin1':
            return self.p.ref_begin1 + index_starts_at
        if name == 'target_end_optimal':
            return self.p.ref_end1 + index_starts_at
        if name == 'query_begin1':
            return self.p.read_begin1 + index_starts_at
        if name == 'query_end1':
            return self.p.read_end1 + index_starts_at
        if name == 'target_end_suboptimal':
            return self.p.ref_end2 + index_starts_at
        if name == 'cigar':
            return self._get_cigar()
        if name == 'query_sequence':
            return self.read_sequence
        if name == 'target_sequence':
            return self.reference_sequence
        else:
            raise AttributeError(\
                "'AlignmentStructure' object has no attribute '%s'" % name)

    def set_zero_based(self, is_zero_based):
        """
        """
        if is_zero_based:
            self.index_starts_at = 0
        else:
            self.index_starts_at = 1

    def is_zero_based(self):
        """
        """
        return self.index_starts_at == 0

cdef class StripedSmithWaterman:
    cdef s_profile *profile
    cdef np.uint8_t weight_gap_open
    cdef np.uint8_t weight_gap_extension
    cdef np.uint8_t bit_flag
    cdef np.uint16_t score_filter
    cdef np.int32_t distance_filter
    cdef np.int32_t mask_length
    cdef str read_sequence
    cdef int index_starts_at
    cdef bool is_protein
    cdef np.ndarray __KEEP_IT_IN_SCOPE_read
    cdef np.ndarray __KEEP_IT_IN_SCOPE_matrix

    def __init__(self, read_sequence, **kwargs):
        """
        """
        self.read_sequence = read_sequence
        self.weight_gap_open = 5 # From BLASTN defaults
        self.weight_gap_extension = 2 # From BLASTN defaults
        self.bit_flag = 1
        self.score_filter = 0 # no filter
        self.distance_filter = 0 # no filter
        temp_mask = len(self.read_sequence)/2
        self.mask_length = 15 if temp_mask < 15 else temp_mask
        # https://www.cs.utexas.edu/users/EWD/transcriptions/EWD08xx/EWD831.html
        self.index_starts_at = 0 # Dijkstra knows what's up

        self._handle_shared_kwargs(**kwargs)

        cdef np.ndarray[np.int8_t, ndim=1, mode="c"] read_seq
        read_seq = self._seq_converter(read_sequence, self.is_protein)

        cdef np.int32_t read_length
        read_length = len(read_sequence)

        cdef np.int8_t score_size
        score_size = 2
        if 'score_size' in kwargs and kwargs['score_size'] is not None:
            score_size = kwargs['score_size']

        cdef np.ndarray[np.int8_t, ndim=1, mode="c"] matrix
        if 'protein' in kwargs and kwargs['protein'] == True:
            self.is_protein = True
            matrix = self._convert_dict2d_to_matrix(kwargs['substitution_matrix'])
        else:
            self.is_protein = False
            if 'substitution_matrix' in kwargs and kwargs['substitution_matrix'] is not None:
                matrix = self._convert_dict2d_to_matrix(kwargs['substitution_matrix'])
            else:
                match = 2 # BLASTN default
                mismatch = 3 # BLASTN default
                if 'match' in kwargs and kwargs['match'] is not None:
                    match = kwargs['match']
                if 'mismatch' in kwargs and kwargs['mismatch'] is not None:
                    mismatch = kwargs['mismatch']
                matrix = self._build_match_matrix(match, mismatch)

        cdef np.int32_t m_width
        m_width = 24 if self.is_protein else 5

        cdef s_profile* p
        self.profile = ssw_init(<np.int8_t*> read_seq.data, 
                                read_length, 
                                <np.int8_t*> matrix.data, 
                                m_width, 
                                score_size)

        # A hack to keep the python GC from eating our data
        self.__KEEP_IT_IN_SCOPE_read = read_seq
        self.__KEEP_IT_IN_SCOPE_matrix = matrix

    def __call__(self, reference_sequence, **kwargs):
        """
        """
        if kwargs:
            self._handle_shared_kwargs(**kwargs)

        cdef np.ndarray[np.int8_t, ndim=1, mode="c"] reference
        reference = self._seq_converter(reference_sequence, self.is_protein)

        cdef np.int32_t ref_length
        ref_length = len(reference_sequence)

        cdef s_align *align
        align = ssw_align(self.profile, <np.int8_t*> reference.data, ref_length, 
                          self.weight_gap_open, self.weight_gap_extension, 
                          self.bit_flag, self.score_filter, self.distance_filter, 
                          self.mask_length)

        alignment = AlignmentStructure(self.read_sequence, reference_sequence, 
                                       self.index_starts_at)
        alignment.__constructor__(align) # Hack to get a pointer through
        return alignment

    def __dealloc__(self):
        if self.profile is not NULL:
            init_destroy(self.profile)

    def _handle_shared_kwargs(self, **kwargs):
        if 'weight_gap_open' in kwargs and kwargs['weight_gap_open'] is not None:
            self.weight_gap_open = kwargs['weight_gap_open']

        if 'weight_gap_extension' in kwargs and kwargs['weight_gap_extension'] is not None:
            self.weight_gap_extension = kwargs['weight_gap_extension']

        if 'bit_flag' in kwargs and kwargs['bit_flag'] is not None:
            self.bit_flag = kwargs['bit_flag']

        if 'score_filter' in kwargs and kwargs['score_filter'] is not None:
            self.score_filter = kwargs['score_filter']

        if 'distance_filter' in kwargs and kwargs['distance_filter'] is not None:
            self.distance_filter = kwargs['distance_filter']

        if 'mask_length' in kwargs and kwargs['mask_length'] is not None:
            self.mask_length = kwargs['mask_length']

        if 'zero_based' in kwargs and kwargs['zero_based'] is not None:
            if not kwargs['zero_based']:
                self.index_starts_at = 1

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
        return seq

    cdef np.ndarray[np.int8_t, ndim=1, mode="c"] _build_match_matrix(self, match, mismatch):
        sequence_order = "ACGT*"
        dict2d = {}
        for row in sequence_order:
            dict2d[row] = {}
            for column in sequence_order:
                if column == '*' or row == '*':
                    dict2d[row][column] = 0
                else:
                    dict2d[row][column] = match if row == column else -mismatch
        return self._convert_dict2d_to_matrix(dict2d)

    cdef np.ndarray[np.int8_t, ndim=1, mode="c"] _convert_dict2d_to_matrix(self, dict2d):
        if self.is_protein:
            sequence_order = "ARNDCQEGHILKMFPSTWYVBZX*"                  
        else:
            sequence_order = "ACGT*"
        i = 0
        length = len(sequence_order)
        cdef np.ndarray[np.int8_t, ndim=1, mode="c"] py_list_matrix = np.empty(length*length, dtype=np.int8)
        for row in sequence_order:
            for column in sequence_order:
                py_list_matrix[i] = dict2d[row][column]
                i+=1
        return py_list_matrix

def striped_smith_waterman_alignment(query, target, **kwargs):
    """
    """
    ssw = StripedSmithWaterman(query, **kwargs)
    return ssw(target)