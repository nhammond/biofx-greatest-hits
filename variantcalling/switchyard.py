#!/usr/bin/env python

"""
switchyard takes short reads from a given genomic interval
and performs local reassembly creating a read-aware debruijn
graph. Read continuity is used to select a set of reasonable
paths that traverse the interval. Paths are assumed to contain
the read entirely, although ends may be trimmed based on
base quality.

After generating a set of possible haplotypes corresponding to
these paths, switchyard re-aligns the reads using each haplotype
as reference, then scores each haplotype based on attributes of
the alignment (e.g. depth and uniformity of coverage).

The best possible genotype is then selected based on score.
"""

import copy
import uuid

REFERENCE_SAMPLE_ID = -1
REFERENCE_STRAND_ID = -1
DEFAULT_KMER_SIZE = 9


class Node:
    """Node represents kminus1mer
    """
    idcounter = 0

    def __init__(self, kminus1mer):
        self.kminus1mer = kminus1mer
        self.incomingedges = set()
        self.outgoingedges = set()
        self.nodeid = '%s_%s' % (self.idcounter, self.kminus1mer)
        Node.idcounter += 1

    def addincomingedge(self, edge):
        self.incomingedges.add(edge)

    def addoutgoingedge(self, edge):
        self.outgoingedges.add(edge)

    def getstartingstrandedges(self):
        # Strand edges that start on this node, grouped by edge
        newstarts = {}
        for edge in self.outgoingedges:
            newstarts[edge] = edge.getstartingstrandedges()
        return newstarts


class Edge:
    """Edge represents kmer
    """
    kmersize = None

    def __init__(self, fromnode, tonode):
        self.fromnode = fromnode
        self.tonode = tonode
        self.fromnode.addoutgoingedge(self)
        self.tonode.addincomingedge(self)
        self.strandedges = set()

    def getweight(self):
        return len(self.strandedges)

    def getstartingstrandedges(self):
        # Strand edges that begin on this edge
        return set(filter(
            lambda e: e.previousstrandedge is None,
            self.strandedges))

    def getsequenceletter(self):
        return self.tonode.kminus1mer[-1]


class Strand:
    """A Strand is a collection of StrandEdges that represent the path
    a single strand (or graphed sequence) takes through the graph.
    """
    idcounter = 0

    def __init__(self, strandid, sampleid):
        self.strandid = strandid
        self.sampleid = sampleid
        self.edgecount = 0
        self.firststrandedge = None
        self.laststrandedge = None
        self.strandid = self.idcounter
        Strand.idcounter += 1

    def getfirstnode(self):
        if self.firststrandedge:
            return self.firststrandedge.fromnode
        else:
            return None

    def getlastnode(self):
        if self.laststrandedge:
            return self.laststrandedge.tonode
        else:
            return None


class StrandEdge:
    """A StrandEdge strand represents the passage of 1 read along
    an edge. Collectively all the StrandEdges for 1 read
    represent the path that read takes through the graph.
    """
    def __init__(self, edge, strand, position):
        self.fromnode = edge.fromnode
        self.tonode = edge.tonode
        self.nextstrandedge = None
        self.previousstrandedge = None
        self.strand = strand
        self.edge = edge
        self.position = position

    def linktoprevious(self, previousstrandedge):
        self.previousstrandedge = previousstrandedge
        previousstrandedge.nextstrandedge = self


class Graph:
    """Graph is a DeBruijn Graph that overlays all Strands,
    including one strand representing the reference sequence
    if one is given.
    """
    def __init__(self, referencesequence, kmersize=DEFAULT_KMER_SIZE):
        self.kmersize = kmersize
        self.strands = {}
        self.nodes = {}
        self.edges = {}
        if referencesequence:
            self.referencestrand = self.createstrand(
                referencesequence, REFERENCE_STRAND_ID, REFERENCE_SAMPLE_ID)
        else:
            self.referencestrand = None
        self.linearizednodes = None
        self.linearizededges = None

    def addsequences(self, sequencelist, sampleid=REFERENCE_SAMPLE_ID):
        for (sequenceid, sequence) in sequencelist:
            strand = self.createstrand(sequence, sequenceid, sampleid)

    def getnodes(self):
        return set(self.nodes.values())

    def getedges(self):
        return set(self.edges.values())

    def getfirstnode(self):
        if self.referencestrand:
            return self.referencestrand.getfirstnode()
        else:
            raise Exception(
                "firstnode undefined because no reference strand is set.")

    def getlastnode(self):
        if self.referencestrand:
            return self.referencestrand.getlastnode()
        else:
            raise Exception(
                "lastnode undefined because no reference strand is set.")

    def createstrand(self, sequence, strandid, sampleid):
        strand = Strand(strandid, sampleid)
        self.strands[strandid] = strand
        previous = None
        position = 0
        for kmer in self.kmergenerator(sequence):
            previous = self.createstrandedgeinorder(
                kmer, strand, previous, position)
            position += 1
        strand.laststrandedge = previous
        return strand

    def createstrandedgeinorder(self, kmer, strand, previous, position):
        # Create edge between left k-minus-1mer and right k-minus-1mer nodes
        edge = self.getorcreateedge(kmer)
        strandedge = StrandEdge(edge, strand, position)
        if not strand.firststrandedge:
            strand.firststrandedge = strandedge
        edge.strandedges.add(strandedge)
        if previous:
            strandedge.linktoprevious(previous)
        strand.edgecount += 1
        return strandedge

    def getorcreateedge(self, kmer):
        try:
            edge = self.edges[kmer]
        except KeyError:
            fromnode = self.getorcreatenode(kmer[:-1])
            tonode = self.getorcreatenode(kmer[1:])
            edge = Edge(fromnode, tonode)
            self.edges[kmer] = edge
        return edge

    def getorcreatenode(self, kminus1mer):
        try:
            node = self.nodes[kminus1mer]
        except KeyError:
            node = Node(kminus1mer)
            self.nodes[kminus1mer] = node
        return node

    def _getstrandedgecounts(self):
        return [s.edgecount + 1 for s in self.strands.values()]

    def getmedianstrandedgecount(self):
        edgecounts = self._getstrandedgecounts()
        n = len(edgecounts)
        if n < 1:
            return None
        if n % 2 == 1:
            return sorted(edgecounts)[n//2]
        else:
            return sum(sorted(edgecounts)[n//2-1:n//2+1])/2.0

    def kmergenerator(self, sequence):
        # chop sequence into kmers
        for i in xrange(0, len(sequence)-(self.kmersize-1)):
            yield sequence[i:i+self.kmersize]

    def saveimage(self, filenameroot, boldnode=None):
        dotfile = filenameroot + '.dot'
        with open(dotfile, 'w') as f:
            f.write("digraph \"Graph\" {\n")
            f.write("  bgcolor=\"white\";\n")
            for node in self.getnodes():
                emphasis = ''
                kminus1mer = node.kminus1mer
                label = kminus1mer
                f.write("  %s [label=\"%s\"] ;\n" % (
                    node.nodeid, label))
            for edge in self.getedges():
                f.write(
                    '  %s -> %s [label=\"%d\" color=\"red\"] ;\n' % (
                        edge.fromnode.nodeid,
                        edge.tonode.nodeid,
                        edge.getweight(),
                    ))
            f.write("}\n")

        import pydot
        pngfile = filenameroot + '.png'
        (graph,) = pydot.graph_from_dot_file(dotfile)
        graph.write_png(pngfile)
        print "saved png file %s" % pngfile


class StrandBundle:
    def __init__(self, node, strandedges):
        # The node us usually redundant here since it should be the tonode
        # for every strandedge. But it's simpler just to require it all the
        # time.
        self.node = node
        self.strandedges = strandedges  # should be a set, may be empty

    def getdepth(self):
        return len(self.strandedges)

    def getnextedges(self, prune=True):
        # Follow strands forward one step, excluding those that end
        nextstrandedgeslist = filter(
            lambda e: e is not None,
            map(lambda e: e.nextstrandedge,
                self.strandedges))
        # Organize by edge
        nextstrandedges = {}
        for edge in self.node.outgoingedges:
            nextstrandedges[edge] = StrandBundle(edge.tonode, set())
        map(lambda e: nextstrandedges[e.edge].strandedges.add(e),
            nextstrandedgeslist)
        if prune:
            # Remove edges with no strands
            for edge, strandedges in nextstrandedges.items():
                if strandedges.getdepth() == 0:
                    del nextstrandedges[edge]
        return nextstrandedges

    def addstrandedges(self, strandedges):
        self.strandedges.update(strandedges)

    def popifoverlapexceeds(self, overlap):
        removededges = filter(
            lambda e: e.position + 1 >= overlap,
            self.strandedges)
        self.strandedges.difference_update(removededges)
        return removededges


class PathNode:
    def __init__(self, edge, depth):
        self.edge = edge
        self.depth = depth


class Path:
    class ReachedMaxPathLength(Exception):
        pass

    def __init__(self, nodelist=None, maxpathlength=1000):
        if nodelist is None:
            self.nodelist = []
        else:
            self.nodelist = copy.copy(nodelist)
        self.maxpathlength = maxpathlength

    def addedge(self, edge, depth):
        if len(self.nodelist) == self.maxpathlength:
            raise self.ReachedMaxPathLength
        self.nodelist.append(PathNode(edge, depth))

    def clone(self):
        return Path(nodelist=self.nodelist)

    def getscore(self):
        pass

    def getsequence(self):
        if len(self.nodelist) == 0:
            return ''
        sequence = self.nodelist[0].edge.fromnode.kminus1mer
        for node in self.nodelist:
            sequence += node.edge.getsequenceletter()
        return sequence


class PathFinder:

    def __init__(self, graph, minimumoverlap=None):
        self.graph = graph
        if minimumoverlap is None:
            self.minimumoverlap = self.graph.getmedianstrandedgecount() * 0.5
        else:
            self.minimumoverlap = minimumoverlap
        self.firstnode = self.graph.getfirstnode()
        self.lastnode = self.graph.getlastnode()
        self.paths = []

    def getpaths(self):
        for edge in self.firstnode.outgoingedges:
            guidestrands = StrandBundle(edge.tonode, edge.strandedges)
            overlappingstrands = StrandBundle(
                edge.tonode, set())  # Initially empty
            path = Path()
            path.addedge(edge, guidestrands.getdepth())
            self._scan(
                guidestrands, overlappingstrands, path)
        return self.paths

    def _scan(self, guidestrands, overlappingstrands, path):
        # We follow guidestrands as long as they follow the same path.
        # As we go, we collect new strands that start along the way in
        # overlappingstrands. So both guidestrands and overlappingstrands
        # lie completely in the path so far.

        # Loop until the guidestrands fork or until dead end
        while True:
            # Save the path when we hit lastnode.
            # We may still continue to find paths that
            # contain the lastnode more than once.
            if guidestrands.node is self.lastnode:
                self.paths.append(path.clone())

            nextguidestrands = guidestrands.getnextedges()
            nextedges = nextguidestrands.keys()
            nextoverlappingstrands = overlappingstrands.getnextedges(
                prune=False)
            newstarts = guidestrands.node.getstartingstrandedges()
            for edge, bundle in nextoverlappingstrands.iteritems():
                bundle.addstrandedges(newstarts[edge])

            if len(nextedges) == 0:
                # Dead end
                return

            if len(nextedges) == 1:
                # Step forward one edge
                edge = nextedges.pop()
                guidestrands = nextguidestrands[edge]
                overlappingstrands = nextoverlappingstrands[edge]
                # if any overlappingstrands have reached the
                # overlap threshold, add them to the guidestrands
                newguidestrands = overlappingstrands.popifoverlapexceeds(
                    self.minimumoverlap)
                guidestrands.addstrandedges(newguidestrands)
                try:
                    path.addedge(edge, guidestrands.getdepth())
                except ReachedMaxPathLength:
                    return
                continue

            else:
                # Guidestrands are splitting. Recurse to follow each path.
                for edge, guidestrands in nextguidestrands.iteritems():
                    overlappingstrands = nextoverlappingstrands.get(edge)
                    # if any overlappingstrands have reached the
                    # overlap threshold, add them to the guidestrands
                    newguidestrands = overlappingstrands.popifoverlapexceeds(
                        self.minimumoverlap)
                    guidestrands.addstrandedges(newguidestrands)
                    newpath = path.clone()
                    try:
                        newpath.addedge(edge, guidestrands.getdepth())
                    except newpath.ReachedMaxPathLength:
                        return
                    self._scan(guidestrands, overlappingstrands, newpath)
                return


if __name__ == '__main__':
    ref = 'agagtcgctacttatcgtgtctaactaatgggtacttcagctcaagagat'
    # MNP complicates the graph by creating a 2nd loop to tactt
    altseq = 'agagtcgctacttatcgtgtctacttaatgggtacttcagctcaagagat'
    # splits on first step
    altseq1 = 'agagtcactacttatcgtgtctaactaatgggtacttcagctcaagagat'
    # splits on second step
    altseq2 = 'agagtcgatacttatcgtgtctaactaatgggtacttcagctcaagagat'
    altseq3 = 'agcggcgcaacatgtcttgtctactataggtgttcatcaactgagaatac'
    # Merge from stray path
    altseq4 = 'cttgccgctacttatcgtgtctaactaatgggtacttcagctcaagagat'

    # This sample is het reference. Build some reads out of both
    # sequences.

    readlength = 10
    kmersize = 7
    sampleid = 1

    def getreadsandids(sequence, readlength):
        readidpairs = []
        for i in xrange(len(sequence) - readlength + 1):
            readidpairs.append((str(uuid.uuid4()), sequence[i:i+readlength]))
        return readidpairs
    readidpairs = getreadsandids(ref, readlength)
    readidpairs.extend(getreadsandids(altseq, readlength))
    readidpairs.extend(getreadsandids(altseq1, readlength))
    readidpairs.extend(getreadsandids(altseq2, readlength))
    readidpairs.extend(getreadsandids(altseq3, readlength))
    readidpairs.extend(getreadsandids(altseq4, readlength))

    graph = Graph(ref, kmersize=kmersize)
    graph.saveimage('ref')
    graph.addsequences(readidpairs, sampleid=sampleid)
    graph.saveimage('test')
    paths = PathFinder(graph).getpaths()
    for path in paths:
        print path.getsequence()
