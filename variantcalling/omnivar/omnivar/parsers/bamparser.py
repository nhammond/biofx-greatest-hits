import pysam

class BamParser(object):

    def __init__(self,filename):
        self.filename = filename
        self.samfile = pysam.AlignmentFile(filename, 'rb')

    def get_reads(self, contig, start, stop):
        return self.samfile.fetch(contig, start, stop)
