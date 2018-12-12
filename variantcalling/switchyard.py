#!/usr/bin/env python

"""
Switchyard takes short reads from a given genomic interval
and performs local reassembly creating a read-aware debruijn
graph. Read continuity is used to select a set of reasonable
paths that traverse the interval. Paths are assumed to contain
the read entirely, although ends may be trimmed based on
base quality.

Switchyard then scores each path (haplotype) based on attributes of
the path and its read coverage (e.g. edge weight, number of duplicate
reads, path length). Genotype is then inferred using haplotype scores.
"""

import copy
import uuid

REFERENCE_SAMPLE_ID = -1
REFERENCE_STRAND_ID = -1

DEFAULT_KMER_SIZE = 7
DEFAULT_MAX_PATH_LENGTH = 1000


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

    def getstartingstrandedges(self, sampleids=None):
        # Strand edges that start on this node, grouped by edge
        newstarts = {}
        for edge in self.outgoingedges:
            newstarts[edge] = edge.getstartingstrandedges(
                sampleids=sampleids)
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

    def getweight(self, sampleids=None):
        return len(self.getstrandedges(sampleids=sampleids))

    def getstrandedges(self, sampleids=None):
        if sampleids:
            return filter(
                lambda e: e.strand.sampleid in sampleids,
                self.strandedges)
        else:
            return self.strandedges

    def getstartingstrandedges(self, sampleids=None):
        # Strand edges that begin on this edge
        startingstrandedges = set(filter(
            lambda e: e.previousstrandedge is None,
            self.strandedges))
        if sampleids:
            startingstrandedges = filter(
                lambda e: e.strand.sampleid in sampleids,
                startingstrandedges)
        return startingstrandedges

    def getsequenceletter(self):
        return self.tonode.kminus1mer[-1]


class Strand:
    """A Strand is a collection of StrandEdges that represent the path
    a single strand (or graphed sequence) takes through the graph.
    """
    def __init__(self, strandid, sampleid):
        self.strandid = strandid
        self.sampleid = sampleid
        self.edgecount = 0
        self.firststrandedge = None
        self.laststrandedge = None

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
        self.sampleids = set()
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
        self.sampleids.add(sampleid)
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

    def getweight(self):
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
                if strandedges.getweight() == 0:
                    del nextstrandedges[edge]
        return nextstrandedges

    def getweighthistory(self):
        strandedges = self.strandedges
        weighthistory = []
        while len(strandedges) > 0:
            strandedges = filter(
                lambda e: e is not None,
                map(lambda e: e.previousstrandedge,
                    strandedges))
            weighthistory.append(len(strandedges))
        return weighthistory

    def addstrandedges(self, strandedges):
        self.strandedges.update(strandedges)

    def popifoverlapexceeds(self, overlap):
        removededges = filter(
            lambda e: e.position + 1 >= overlap,
            self.strandedges)
        self.strandedges.difference_update(removededges)
        return removededges


class PathNode:
    def __init__(self, edge, weight):
        self.edge = edge
        self.weight = weight


class Path:
    class ReachedMaxPathLength(Exception):
        pass

    def __init__(self, edges=None, maxpathlength=None):
        if edges is None:
            self.nodelist = []
        else:
            self.nodelist = [PathNode(edge, 0) for edge in edges]
        self.maxpathlength = maxpathlength

    def addedge(self, edge, weight=0):
        if self.maxpathlength and len(self.nodelist) >= self.maxpathlength:
            raise self.ReachedMaxPathLength
        self.nodelist.append(PathNode(edge, weight))

    def clone(self):
        clone = Path(maxpathlength=self.maxpathlength)
        clone.nodelist = copy.copy(self.nodelist)
        return clone

    def getlength(self):
        return len(self.nodelist)

    def getscore(self):
        pass

    def getsequence(self):
        if len(self.nodelist) == 0:
            return ''
        sequence = self.nodelist[0].edge.fromnode.kminus1mer
        for node in self.nodelist:
            sequence += node.edge.getsequenceletter()
        return sequence

    def getweights(self):
        return [node.weight for node in self.nodelist]

    def adjustweightforaddedstrands(self, bundle):
        # Current step is already accounted for.
        # Add weight for history of bundle.
        counter = -1
        for weight in bundle.getweighthistory():
            self.nodelist[counter].weight += weight
            counter -= 1

    def endswith(self, path):
        if self.getlength() < path.getlength():
            return False
        for i in range(1, path.getlength()+1):
            if path.nodelist[-i].edge is not self.nodelist[-i].edge:
                return False
        return True


class PathFinder:

    def __init__(self, graph, minimumoverlap=None, anchorlength=None,
                 maxpathlength=DEFAULT_MAX_PATH_LENGTH, excludedsamples=None,
                 includedsamples=None):

        self.graph = graph
        kmersize = self.graph.kmersize
        if minimumoverlap is None:
            minimumoverlap = self.graph.getmedianstrandedgecount() * 0.5 \
                             + kmersize - 1
        if minimumoverlap < kmersize:
            raise Exception('Minimum overlap must be at least kmersize.'
                            'minimumoverlap=%s, kmersize=%s.' %
                            (minimumoverlap, kmersize))
        self.minimumoverlapedges = int(minimumoverlap - kmersize + 1)
        if anchorlength is None:
            anchorlength = self.graph.kmersize
        if anchorlength < kmersize:
            raise Exception('Anchor length must be at least kmersize.'
                            'anchorlength=%s, kmersize=%s.' %
                            (anchorlength, kmersize))
        self.anchorlength = anchorlength
        self.maxpathlength = maxpathlength
        self._initializeterminalpath()
        self.paths = []
        if excludedsamples is not None and includedsamples is not None:
            raise Exception(
                "Cannot specify both includedsamples and excludedsamples")
        if includedsamples is not None:
            self.includedsamples = includedsamples
        else:
            if excludedsamples is None:
                excludedsamples = [-1]
            self.includedsamples = self.graph.sampleids.difference(
                excludedsamples)

    def getpaths(self):
        guidestrands, overlappingstrands, path = self._getstartanchorpath()
        self._scan(guidestrands, overlappingstrands, path)
        return self.paths

    def _initializeterminalpath(self):
        # Save the end of the reference path. We will compare this
        # to the paths we traverse to know when we reach the end,
        # as determined by an overlap of anchorlength or more.
        terminaledges = []
        referencestrandedge = self.graph.referencestrand.laststrandedge
        for i in range(self.anchorlength - self.graph.kmersize + 1):
            terminaledges.insert(0, referencestrandedge.edge)
            referencestrandedge = referencestrandedge.previousstrandedge
            if referencestrandedge is None:
                raise Exception(
                    "Reference strand cannot be shorter than anchor length")
        self.terminalpath = Path(terminaledges)

    def _getstartanchorpath(self):
        anchoredges = self.anchorlength - self.graph.kmersize + 1
        referencestrandedge = self.graph.referencestrand.firststrandedge
        guidestrands, overlappingstrands, path = self._initializepath(
            referencestrandedge.edge)
        for i in range(anchoredges - 1):
            referencestrandedge = referencestrandedge.nextstrandedge
            assert referencestrandedge is not None, \
                "Reference strand cannot be shorter than anchor length"
            nextguidestrands, nextoverlappingstrands = self._getnext(
                guidestrands, overlappingstrands)
            try:
                guidestrands, overlappingstrands = self._followedge(
                    nextguidestrands, nextoverlappingstrands,
                    referencestrandedge.edge, path)
            except path.ReachedMaxPathLength:
                raise Exception(
                    "Max path length cannot be less than anchor length")
        return guidestrands, overlappingstrands, path

    def _initializepath(self, edge):
        guidestrands = StrandBundle(
            edge.tonode,
            edge.getstrandedges(sampleids=self.includedsamples))
        overlappingstrands = StrandBundle(
            edge.tonode, set())  # Initially empty
        path = Path(maxpathlength=self.maxpathlength)
        path.addedge(edge, weight=guidestrands.getweight())
        return guidestrands, overlappingstrands, path

    def _scan(self, guidestrands, overlappingstrands, path):
        # We follow guidestrands as long as they follow the same path.
        # As we go, we collect new strands that start along the way in
        # overlappingstrands. So both guidestrands and overlappingstrands
        # lie completely in the path so far.

        # Loop until the guidestrands fork or until dead end
        while True:
            # Save the path when we hit the end.
            # We may still continue to find paths that
            # contain this terminus more than once.
            if path.endswith(self.terminalpath):
                # First add overlap strands to weights, to be consistent
                # with weight at start node, which included all strands
                # without respect to minimum overlap
                completedpath = path.clone()
                completedpath.adjustweightforaddedstrands(
                    overlappingstrands)
                self.paths.append(completedpath)

            nextguidestrands, nextoverlappingstrands = self._getnext(
                guidestrands, overlappingstrands)
            nextedges = nextguidestrands.keys()

            if len(nextedges) == 0:
                # Dead end
                return

            if len(nextedges) == 1:
                # Step forward one edge
                edge = nextedges.pop()
                try:
                    guidestrands, overlappingstrands = self._followedge(
                        nextguidestrands, nextoverlappingstrands, edge, path)
                except path.ReachedMaxPathLength:
                    return
                continue

            else:
                # Guidestrands are splitting. Recurse to follow each path.
                for edge in nextguidestrands:
                    splitpath = path.clone()
                    try:
                        guidestrands, overlappingstrands = self._followedge(
                            nextguidestrands, nextoverlappingstrands,
                            edge, splitpath)
                    except path.ReachedMaxPathLength:
                        return
                    self._scan(guidestrands, overlappingstrands, splitpath)
                return

    def _followedge(self, nextguidestrands, nextoverlappingstrands,
                    edge, path):
        guidestrands = nextguidestrands[edge]
        overlappingstrands = nextoverlappingstrands[edge]
        # if any overlappingstrands have reached the
        # overlap threshold, add them to the guidestrands
        newguidestrands = overlappingstrands.popifoverlapexceeds(
            self.minimumoverlapedges)
        guidestrands.addstrandedges(newguidestrands)
        path.adjustweightforaddedstrands(
            StrandBundle(edge.tonode, newguidestrands))
        path.addedge(edge, weight=guidestrands.getweight())
        return guidestrands, overlappingstrands

    def _getnext(self, guidestrands, overlappingstrands):
        nextguidestrands = guidestrands.getnextedges()
        nextoverlappingstrands = overlappingstrands.getnextedges(
            prune=False)
        newstarts = guidestrands.node.getstartingstrandedges(
            sampleids=self.includedsamples)
        for edge, bundle in nextoverlappingstrands.iteritems():
            bundle.addstrandedges(newstarts[edge])
        return nextguidestrands, nextoverlappingstrands


if __name__ == '__main__':
    ref = 'agagtcgctacttatcgtgtctaactaatgggtacttcagctcaagagat'
    graph = Graph(ref)

    # MNP complicates the graph by creating a 2nd loop to tactt
    altseq = 'agagtcgctacttatcgtgtctacttaatgggtacttcagctcaagagat'
    # splits on first step
    altseq1 = 'agagtcactacttatcgtgtctaactaatgggtacttcagctcaagagat'
    # splits on second step
    altseq2 = 'agagtcgatacttatcgtgtctaactaatgggtacttcagctcaagagat'
    altseq3 = 'agcggcgcaacatgtcttgtctactataggtgttcatcaactgagaatac'

    def getreadsandids(sequence):
        readlength = 10
        readidpairs = []
        for i in xrange(len(sequence) - readlength + 1):
            readidpairs.append((str(uuid.uuid4()), sequence[i:i+readlength]))
        return readidpairs

    graph.addsequences(getreadsandids(ref), sampleid=1)
    graph.addsequences(getreadsandids(altseq), sampleid=1)
    graph.addsequences(getreadsandids(altseq1), sampleid=2)
    graph.addsequences(getreadsandids(altseq2), sampleid=2)
    graph.addsequences(getreadsandids(altseq3), sampleid=1)

    graph.saveimage('test')
    paths = PathFinder(graph).getpaths()
    for path in paths:
        weights = path.getweights()
        averageweight = sum(weights)/float(len(weights))
        diffsquared = [(weight-averageweight)**2 for weight in weights]
        stdev = (sum(diffsquared)/float(len(diffsquared)))**(1/2.0)
        minweight = min(weights)
        print path.getsequence()
        print "Edge weights: avg=%s, min=%s, stdev=%s" % (
            averageweight, minweight, stdev)
