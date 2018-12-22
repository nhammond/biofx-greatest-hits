import argparse
import os
import unittest
import uuid
from .pathfinder import PathFinder
from .debruijngraph import Graph
from .parsers.bamparser import BamParser
from .parsers.bedparser import BedParser
from .parsers.fastaparser import FastaParser


package_root = os.path.dirname(__file__)


class OmniVar(object):

    KMER_SIZE = 41

    def __init__(self, bam_file, reference_file, bed_file, sample_id, verbosity='INFO'):
        self.reads_file = BamParser(bam_file)
        self.intervals_file = BedParser(bed_file)
        self.reference_file = FastaParser(reference_file)
        self.verbosity = verbosity

    def run(self):

        for contig, start, stop in self.intervals_file.intervals_generator():
            print "Interval %s %s:%s" % (contig, start, stop)
            reference_sequence = self.reference_file.get_sequence(contig, start, stop)
            print "...reference: %s" % str(reference_sequence)
            reads = self.reads_file.get_reads(contig, start, stop)
            graph = Graph(str(reference_sequence)[:100], kmer_size=self.KMER_SIZE)
            for read in reads:
                graph.add_sequence(str(read.query_sequence))
            pathfinder = PathFinder(graph)
            paths = pathfinder.find_paths()
            if len(paths) == 0:
                print "...no path found"
                # graph.save_image("%s_%s_%s" % (contig, start, stop))
            for path in paths:
                sequence = path.get_sequence()
                if sequence == str(reference_sequence):
                    sequence = "[reference sequence]"
                print "...length: %s sequence: %s" % (path.get_length(), sequence)
                    
            # Save image with the best path
            # graph.save_image("%s_%s_%s" % (contig, start, stop), paths=paths)
            import pdb; pdb.set_trace()


class Version(argparse.Action):
    def __call__(self, parser, namespace, values, option_string):
        with open(os.path.join(package_root,'VERSION')) as version_file:
            version = version_file.read().strip()
        print version
        exit(0)

def _parse_args():
    parser = argparse.ArgumentParser('omnivar')
    parser.add_argument('-v', '--version', nargs=0, action=Version)
    parser.add_argument('-B', '--bam', required=True, type=str)
    parser.add_argument('-r', '--reference', required=True, type=str)
    parser.add_argument('-b', '--bed', required=True, type=str)
    parser.add_argument('-s', '--sampleid', required=True, type=str)
    parser.add_argument(
        '-V', '--verbosity', help='Verbosity of terminal output',
        default='info', choices=['error', 'warning', 'info', 'debug'])
    return parser.parse_args()

def main():
    args = _parse_args()
    bam_file = args.bam
    bed_file = args.bed
    reference_file = args.reference
    sample_id = args.sampleid
    verbosity = args.verbosity.upper()
    OmniVar(bam_file, reference_file, bed_file, verbosity).run()
