#!/usr/bin/env python

import unittest
from omnivar.pathfinder import Path, PathFinder
from omnivar.debruijngraph import Graph


class TestPath(unittest.TestCase):

    def setUp(self):
        self.graph = Graph('abcdefg', kmer_size=3)
        self.sample_id = 1
        first = self.graph.first_edge
        second = list(first.get_next_edges())[0]
        third = list(second.get_next_edges())[0]
        fourth = list(third.get_next_edges())[0]
        first.weight = 10
        second.weight = 10
        third.weight = 10
        fourth.weight = 10
        self.path = Path()
        self.path.add_edge(first)
        self.path.add_edge(second)
        self.path.add_edge(third)
        self.next_edge = fourth

    def testAddEdge(self):
        self.assertEqual(self.path.get_length(), 3)
        self.path.add_edge(self.next_edge)
        self.assertEqual(self.path.get_length(), 4)

    def testGetLength(self):
        self.assertEqual(self.path.get_length(), 3)

    def testGetWeights(self):
        self.assertEqual(self.path.get_weights(), [10, 10, 10])

    def testGetSequence(self):
        self.assertEqual(self.path.get_sequence(), 'abcde')

    def testClone(self):
        clone = self.path.clone()
        
        # removing edges does not affect original
        clone.edge_list.pop()
        self.assertEqual(clone.get_length(), 2)
        self.assertEqual(self.path.get_length(), 3)

    def testEndsWith(self):
        ending = self.path.get_last_edge()
        self.assertTrue(self.path.ends_with(ending))

        self.path.edge_list.pop()
        self.assertFalse(self.path.ends_with(ending))


class TestPathFinder(unittest.TestCase):

    def testFindPaths(self):
        sequence = 'abcdefgh'
        sequence_1 = sequence
        sequence_2 = 'abcefgh'
        graph = Graph(sequence, kmer_size=3)
        graph.add_sequences([sequence_1, sequence_2])

        # We should get all sequences if we have a long enough max_path_length
        paths = PathFinder(graph, minimum_weight=1).find_paths()
        self.assertEqual(len(paths), 2)
        sequences = [path.get_sequence() for path in paths]
        self.assertTrue(sequence_1 in sequences)
        self.assertTrue(sequence_2 in sequences)

    def testGetPathsNoSequences(self):
        graph = Graph('abcdefg', kmer_size=3)
        paths = PathFinder(graph, minimum_weight=1).find_paths()
        self.assertEqual(len(paths), 0)

    """
    def testGetPathFromReference(self):
        sequence = 'abcdefg'
        kmer_size = 3
        number_of_edges = len(sequence) - kmer_size + 1
        graph = Graph(sequence, kmer_size=kmer_size)
        paths = PathFinder(graph, minimum_overlap=10,
                           sample_ids=[REFERENCE_SAMPLE_ID]).find_paths()
        self.assertEqual(len(paths), 1)
        self.assertEqual(paths[0].get_sequence(), sequence)
        self.assertEqual(paths[0].get_weights(), [1]*number_of_edges)

    def testGetPathFromTwoOverlappingStrands(self):
        sequence = 'abcdefgh'
        sequence_1 = sequence[:6]
        sequence_2 = sequence[-6:]
        kmer_size = 3
        number_of_edges = len(sequence) - kmer_size + 1
        graph = Graph(sequence, kmer_size=kmer_size)
        graph.add_sequences([(1, sequence_1)], 1)
        graph.add_sequences([(2, sequence_2)], 1)
        paths = PathFinder(graph, minimum_overlap=3).find_paths()
        self.assertEqual(len(paths), 1)
        self.assertEqual(paths[0].get_sequence(), sequence)
        # six edges, with 2 overlapping
        weights = [1, 1, 2, 2, 1, 1]
        self.assertEqual(paths[0].get_weights(), weights)

    def testGetPathWithMergingStrand(self):
        sequence = 'abcdefgh'
        sequence_1 = sequence
        sequence_2 = 'xyzdefgh'
        graph = Graph(sequence, kmer_size=3)
        graph.add_sequences([(1, sequence_1)], 1)
        graph.add_sequences([(2, sequence_2)], 1)
        paths = PathFinder(graph, minimum_overlap=3).find_paths()
        self.assertEqual(len(paths), 1)
        self.assertEqual(paths[0].get_sequence(), sequence)
        # merged strand should be excluded from weights
        weights = [1, 1, 1, 1, 1, 1]
        self.assertEqual(paths[0].get_weights(), weights)

    def testGetPathWithSplittingStrand(self):
        sequence = 'abcdefgh'
        sequence_1 = sequence
        sequence_2 = 'abcdexyz'
        graph = Graph(sequence, kmer_size=3)
        graph.add_sequences([(1, sequence_1)], 1)
        graph.add_sequences([(2, sequence_2)], 1)
        paths = PathFinder(graph, minimum_overlap=3).find_paths()
        self.assertEqual(len(paths), 1)
        self.assertEqual(paths[0].get_sequence(), sequence)
        # split strand should be excluded from weights
        weights = [1, 1, 1, 1, 1, 1]
        self.assertEqual(paths[0].get_weights(), weights)

    def testGetPathThroughDoubleLoop(self):
        sequence = 'abcdefghefghijkl'
        sequences = [(1, sequence[:8]),
                     (2, sequence[2:10]),
                     (3, sequence[4:12]),
                     (4, sequence[6:14]),
                     (5, sequence[-8:])]
        graph = Graph(sequence, kmer_size=3)
        graph.add_sequences(sequences, 1)

        # weith overlap of 5, we're ok
        paths = PathFinder(graph, minimum_overlap=5).find_paths()
        self.assertEqual(len(paths), 1)
        self.assertEqual(paths[0].get_sequence(), sequence)

        # with overlap of 4, we might be fooled
        paths = PathFinder(
            graph, minimum_overlap=4, max_sequence_length=20).find_paths()
        self.assertEqual(len(paths), 3)

        # but weight might give us a clue which path to favor
        paths = sorted(paths, key=lambda p: len(p.get_sequence()))

        # * abcdefghijkl *
        #   abcdefgh
        #       efghijkl
        # abc bcd def efg fgh ghi hij jkl
        # 1   1   1   2   2   1   1   1
        self.assertEqual(
            paths[0].get_weights(), [1, 1, 1, 1, 2, 2, 1, 1, 1, 1])

        # * abcdefghefghijkl *
        #   abcdefgh
        #     cdefghef
        #       efghefgh
        #         ghefghij
        #           efghijkl
        # abc bcd cde def efg fgh ghe hef efg fgh ghi hij ijk jkl
        #  1   1   2   2   3   3   3   3   3   3   2   2   1   1
        self.assertEqual(
            paths[1].get_weights(),
            [1, 1, 2, 2, 3, 3, 3, 3, 3, 3, 2, 2, 1, 1])

        # * abcdefghefghefghijkl *
        #   abcdefgh
        #     cdefghef
        #       efghefgh
        #           efghefgh (duplicate)
        #             ghefghij
        #               efghijkl
        # abc bcd cde def efg fgh ghe hef efg fgh ghe hef efg fgh ghi hij ijk jkl
        #  1   1   2   2   3   3   2   2   2   2   2   2   3   3   2   2   1   1
        self.assertEqual(
            paths[2].get_weights(),
            [1, 1, 2, 2, 3, 3, 2, 2, 2, 2, 2, 2, 3, 3, 2, 2, 1, 1])
    """

if __name__ == '__main__':
    unittest.main()
