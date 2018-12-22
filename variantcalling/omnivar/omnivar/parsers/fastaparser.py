from pyfaidx import Fasta


class FastaParser(object):

    def __init__(self, filename):
        self.fasta_dict = Fasta(filename)

    def get_sequence(self, contig, start, end):
        sequence = self.fasta_dict[contig]
        sequence_length = len(sequence)
        if contig not in self.fasta_dict:
            raise ValueError(
                'Contig "%s" not found in fasta.' % contig)
        if start < 0 or start >= sequence_length - 1:
            raise ValueError(
                'Start position "%s" out of bounds on contig "%s" with length "%s"'
                % (start, contig, sequence_length))
        if end < 0 or end >= sequence_length - 1:
            raise ValueError(
                'End position "%s" out of bounds on contig "%s" with length "%s"'
                % (end, contig, sequence_length))
        if end < start:
            raise ValueError(
                'End position "%s" less than start position "%s" on contig "%s"'
                % (end, start, contig))
        return sequence[start:end]
