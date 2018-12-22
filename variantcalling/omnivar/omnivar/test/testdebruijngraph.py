#!/usr/bin/env python

import unittest
from omnivar.debruijngraph import Node, Edge, Graph


class TestGraph(unittest.TestCase):

    def setUp(self):
        self.ref_sequence = 'agagtcgctacttatcgtgtctaactaatgggtacttcagctcaagagat'
        self.alt_sequence = 'agagtcgctacttatcgttgggtacttcagctcaagagat'
        self.graph = Graph(self.ref_sequence, kmer_size=7)
        # generate overlapping reads from sequences
        read_length = 20
        self.sample_1_reads = []
        self.sample_2_reads = []
        for i in xrange(len(self.ref_sequence) - read_length + 1):
            self.sample_1_reads.append(self.ref_sequence[i:i+read_length])
        for i in xrange(len(self.alt_sequence) - read_length + 1):
            self.sample_2_reads.append(self.alt_sequence[i:i+read_length])
        self.graph.add_sequences(self.sample_1_reads)
        self.graph.add_sequences(self.sample_2_reads)

    def testInitReferenceLessThanKmerSize(self):
        # ValueError raised if reference sequence too short
        sequence = '123456'
        with self.assertRaises(ValueError):
            graph = Graph(sequence, kmer_size=7)

    def testAddSequenceLessThanKmerSize(self):
        # Too short sequences are silently skipped
        ref_sequence = '1234567'
        sequence = '123456'
        graph = Graph(ref_sequence, kmer_size=7)
        nodes_before = len(graph.get_nodes())
        edges_before = len(graph.get_edges())

        # Confirm that sequence was never added
        graph.add_sequence(sequence)
        self.assertEqual(len(graph.get_nodes()), nodes_before)
        self.assertEqual(len(graph.get_edges()), edges_before)

    def testGetNodes(self):
        sequence = '123456'
        graph = Graph(sequence, kmer_size=5)
        self.assertEqual(len(graph.get_nodes()), 3)

    def testGetEdges(self):
        sequence = '123456'
        graph = Graph(sequence, kmer_size=5)
        self.assertEqual(len(graph.get_edges()), 2)

    def testFirstEdge(self):
        sequence = '123456'
        graph = Graph(sequence, kmer_size=5)
        self.assertEqual(graph.first_edge.get_sequence(),
                         '12345')

    def testLastEdge(self):
        sequence = '123456'
        graph = Graph(sequence, kmer_size=5)
        self.assertEqual(graph.last_edge.get_sequence(),
                         '23456')


if __name__ == '__main__':
    unittest.main()
