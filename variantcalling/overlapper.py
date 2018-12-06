#!/usr/bin/env python

import copy
import uuid

MIN_OVERLAP = 10
class Read:

    # Each node represents a read

    def __init__(self, readid, sequence):
        self.sequence = sequence
        self.readid = readid
        self.overlaps = set()

    def getlength(self):
        return len(self.sequence)

    def createoverlaps(self, read):
        # don't report overlaps of a read with itself
        if read is self:
            return []
        # cheks only self --> read,
        # not read --> self
        newoverlaps = set()
        for overlaplength in xrange(
                MIN_OVERLAP,
                min(read.getlength(), self.getlength())+1):
            if self.sequence[-overlaplength:] == read.sequence[:overlaplength]:
                newoverlaps.add(Overlap(self, read, overlaplength))
        self.overlaps.update(newoverlaps)
        return newoverlaps


class Overlap:
    # Overlap represents overlap from fromread to toread

    def __init__(self, fromread, toread, overlaplength):
        self.fromread = fromread
        self.toread = toread
        self.length = overlaplength

    def getsequence(self):
        return self.fromread.sequence[-self.length:]


class OverlapGraph:

    counter = 0

    def __init__(self, readidpairs, referencesequence=None):
        self.reads = set()
        self.overlaps = set()
        for readid, sequence in readidpairs:
            self.reads.add(Read(readid, sequence))
        for read1 in self.reads:
            for read2 in self.reads:
                self.overlaps.update(read1.createoverlaps(read2))
        if referencesequence:
            self.start = Read('START', referencesequence[:MIN_OVERLAP])
            self.end = Read('END', referencesequence[-MIN_OVERLAP:])
            for read in self.reads:
                self.overlaps.update(self.start.createoverlaps(read))
                self.overlaps.update(read.createoverlaps(self.end))
            self.reads.add(self.start)
            self.reads.add(self.end)
        else:
            self.start = None
            self.end = None
        print "Created %s reads with %s overlaps" % (
            len(self.reads), len(self.overlaps))

    def getcontinuoussequences(self):
        if not self.start and self.end:
            raise Exception(
                "Reference start and end not defined. Cannot get paths")
        sequence = self.start.sequence
        visitedreads = set()
        allsequences = set()
        self._crawl(self.start, sequence, visitedreads, allsequences)
        print self.counter
        return allsequences

    def _crawl(self, read, originalsequence, visitedreadsoriginalcopy, allsequences):
        if read == self.end:
            allsequences.add(originalsequence)
            self.counter += 1
            return allsequences
        for overlap in read.overlaps:
            if overlap.toread in visitedreadsoriginalcopy:
                continue
            visitedreads = copy.copy(visitedreadsoriginalcopy)
            visitedreads.add(overlap.toread)
            sequence = originalsequence + \
                       overlap.toread.sequence[overlap.length:]
            self._crawl(overlap.toread, sequence, visitedreads, allsequences)

    def saveimage(self, filenameroot):
        dotfile = filenameroot + '.dot'
        with open(dotfile, 'w') as f:
            f.write("digraph \"Graph\" {\n")
            f.write(" bgcolor=\"white\";\n")
            for read in self.reads:
                f.write("  %s [label=\"%s\"] ;\n" % (
                    read.readid, read.sequence))
            for overlap in self.overlaps:
                f.write(" %s -> %s [label=\"%d\" color=red] ;\n" %(
                    overlap.fromread.readid, overlap.toread.readid,
                    overlap.length))
            f.write("}\n")
        import pydot
        pngfile = filenameroot + '.png'
        (graph,) = pydot.graph_from_dot_file(dotfile)
        graph.write_png(pngfile)
        print "saved png file %s" % pngfile


if __name__=='__main__':
    ref = 'agagtcgctacttatcgtgtctaactaatgggtacttc' #agctcaagagat'
    # MNP complicates the graph by creating a 2nd loop to tactt
    altseq = 'agagtcgctacttatcgtgtctacttaatgggtacttcagctcaagagat'
    # splits on first step
    altseq1 = 'agagtcactacttatcgtgtctaactaatgggtacttcagctcaagagat'
    # splits on second step
    altseq2 = 'agagtcgatacttatcgtgtctaactaatgggtacttcagctcaagagat'
    # This sample is het reference. Build some reads out of both
    # sequences.

    readlength = 15
    sampleid = 1

    def getreadsandids(sequence, readlength):
        readidpairs = []
        for i in xrange(len(sequence) - readlength + 1):
            readidpairs.append((str(uuid.uuid4())[:8], sequence[i:i+readlength]))
        return readidpairs
    readidpairs = getreadsandids(ref, readlength)
    import pdb; pdb.set_trace()
#    readidpairs.extend(getreadsandids(altseq, readlength))
#    readidpairs.extend(getreadsandids(altseq1, readlength))
#    readidpairs.extend(getreadsandids(altseq2, readlength))

#    readidpairs = [(str(uuid.uuid4())[:8], 'agagtcgctacttatcgtgtct'),
#                   (str(uuid.uuid4())[:8], 'cttatcgtgtctaactaatgggtacttca'),
#                   (str(uuid.uuid4())[:8], 'aactaatgggtacttcagctcaagagat'),
#    ]
    graph = OverlapGraph(readidpairs, referencesequence=ref)
    graph.saveimage('test')
    print graph.getcontinuoussequences()
