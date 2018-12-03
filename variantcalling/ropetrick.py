#!/usr/bin/env python

import copy
import uuid

"""ropetrick is designed to generate a list of possible haplotypes over a
specified genome interval using short read data.
Start and end positions of the interval must be on the reference genome

1) start and end positions are on the reference genome and have good 
coverage

2) there exist one or more paths with continuous read coverage from start
to end

ropetrick works by building a DebruijnGraph from kmers of the
input sequences, which include the reference sequence on the specified interval
and a set of reads believed to fall in that interval (typically reads from a bam
file aligned with bwa or another read aligner)

Note that ropetrick does not have error tolerance, although error-tolerant
alignment is useful for identifying the reads provided as input. ropetrick 
identifies paths for which there is actual evidence in the reads provided.
"""


REFERENCE_SAMPLE_ID = -1
REFERENCE_STRAND_ID = -1
DEFAULT_KMER_SIZE = 5


class Node:
    """Node represents kminus1mer
    """
    def __init__(self, kminus1mer):
        self.kminus1mer = kminus1mer
        self.incomingedges = []
        self.outgoingedges = []

    def addincomingedge(self, edge):
        self.incomingedges.append(edge)

    def addoutgoingedge(self, edge):
        self.outgoingedges.append(edge)

class Edge:
    """Edge represents kmer
    """
    def __init__(self, fromnode, tonode):
        self.fromnode = fromnode
        self.tonode = tonode
        self.fromnode.addoutgoingedge(self)
        self.tonode.addincomingedge(self)
        self.strandedges = []

    def getweight(self):
        return len(self.strandedges)

class Strand:
    """A Strand is a collection of StrandEdges that represent the path
    a single strand (or graphed sequence) takes through the graph.
    """
    def __init__(self, strandid, sampleid):
        self.strandid = strandid
        self.sampleid = sampleid
        self.firststrandedge = None
        self.laststrandedge = None
        self.strandedges = []
        self.nodes = []

    def getfirstnode(self):
        if self.firststrandedge is None:
            return None
        return self.firststrandedge.fromnode

    def getlastnode(self):
        if self.laststrandedge is None:
            return None
        return self.laststrandedge.tonode

    def getsequence(self):
        if self.firststrandedge is None:
            return None
        else:
            seq = self.getfirstnode().kminus1mer
            seq += self.firststrandedge.getsequencetoendofstrand()
            return seq

    def addstrandedge(self, strandedge):
        # StrandEdges must be created in sequence order.
        # Set first and last edge according to order added.
        self.strandedges.append(strandedge)
        if self.firststrandedge is None:
            self.firststrandedge = strandedge
        self.laststrandedge = strandedge


class StrandEdge:
    """A StrandEdge strand represents the passage of 1 read along
    an edge. Collectively all the EdgeStrands for 1 read
    represent the path that read takes through the graph.
    """
    def __init__(self, strand, edge):
        self.edge = edge
        self.edge.strandedges.append(self)
        self.fromnode = self.edge.fromnode
        self.tonode = self.edge.tonode
        self.nextstrandedge = None
        self.previousstrandedge = None
        self.strand = strand
        self.strand.addstrandedge(self)

    def linktoprevious(self, previousstrandedge):
        self.previousstrandedge = previousstrandedge
        previousstrandedge.nextstrandedge = self

    def getsequencetoendofstrand(self):
        seq = self.getsequenceletter()
        if self.nextstrandedge:
            seq += self.nextstrandedge.getsequencetoendofstrand()
        return seq

    def getsequenceletter(self):
        return self.tonode.kminus1mer[-1]


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

    def addsequences(self, sequencelist, sampleid=REFERENCE_SAMPLE_ID):
        for (sequenceid, sequence) in sequencelist:
            strand = self.createstrand(sequence, sequenceid, sampleid)

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

    def getcontinuouspaths(self):
        ropetrick = RopeTrick(self.getfirstnode(), self.getlastnode())
        return ropetrick.getpaths()

    def createstrand(self, sequence, strandid, sampleid):
        strand = Strand(strandid, sampleid)
        self.strands[strandid] = strand
        previous = None
        for kmer in self.kmergenerator(sequence):
            previous = self.createstrandedge(kmer, strand, previous)
        return strand

    def createstrandedge(self, kmer, strand, previous):
        # Create edge between left k-minus-1mer and right k-minus-1mer nodes
        edge = self.getorcreateedge(kmer)
        strandedge = StrandEdge(strand, edge)
        if previous:
            strandedge.linktoprevious(previous)
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

    def kmergenerator(self, sequence):
        # chop sequence into kmers
        for i in xrange(0, len(sequence)-(self.kmersize-1)):
            yield sequence[i:i+self.kmersize]

    def saveimage(self, filenameroot):
        dotfile = filenameroot + '.dot'
        with open(dotfile, 'w') as f:
            f.write("digraph \"Graph\" {\n")
            f.write("  bgcolor=\"white\";\n")
            for kminus1mer, node in self.nodes.iteritems():
                if node.incomingedges:
                    label = '...'+kminus1mer[-1]
                else:
                    label = kminus1mer
		f.write("  %s [label=\"%s\"] ;\n" % (
                    kminus1mer, label))
            for kmer, node in self.nodes.iteritems():
                for edge in node.outgoingedges:
                    f.write("  %s -> %s [label=\"%d\" color=red] ;\n" % (
                        edge.fromnode.kminus1mer, edge.tonode.kminus1mer,
                        edge.getweight()))
            f.write("}\n")

        import pydot
        pngfile = filenameroot + '.png'
        (graph,) = pydot.graph_from_dot_file(dotfile)
        graph.write_png(pngfile)
        print "saved png file %s" % pngfile

class RopeTrick:

    class DeadEndPath(Exception):
        pass

    def __init__(self, initialnode, terminalnode, path=None, strandedges=None):
        # Break when we reach this node
        self.terminalnode = terminalnode

        self.initialnode = initialnode

        # RopeTrick identifies all continuous paths on the DeBruijn graph
        # from path.getcurrentnode() to terminalnode for which continuous read
        # eveidence exists, requiring that the reads used as evidence lie
        # entirely on the path.
        #
        # Rather than travel all possible paths through overlapping reads,
        # we nominate a lead strand and implicitly follow other strands that
        # align with the lead strand. When the other strands diverge from
        # the lead, we spawn a new RopeTrick object to analyze the fork.
        # When a lead strand ends, we nominate a new lead and continue.
        #
        # When evidence exists for a path that diverges from leadstrand,
        # spawn a new RopeTrick that excludes leadstrand and its followers,
        # roll back the path to where that leadstrand started, and resume
        # searching from there.
        #
        # self.strandedges corresponds to strands that are silently following
        # the lead strand. These are strands that are compatible with the
        # path thus far and have not yet provided evidence for a path other
        # than the one followed by the lead strand.
        self.strandedges = strandedges

        # self.path records all traversed nodes, plus additional info required for
        # rewinding the RopeTrick to an earlier position along the path:
        # specifically, any strandedges that were dropped because the
        # strand end was reached. These must be added back when rewinding.
        self.path = path
        if self.path:
            self.leadstrandedge = self.path.getcurrentleadstrandedge()

        # self.paths is the final result, a list of all paths with read support
        # from initial node to terminal node. That includes self.path from this
        # RopeTrick (but only if it does not dead-end) and the results from
        # recursive spawning of other RopeTricks on divergent paths.
        #
        # Don't add self.path to self.paths until we know it's
        # not a dead end
        self.paths = []

        # Between steps, self._getcurrentnode() is normally the strandedge.fromnode for
        # every strandedge in self.strandedges, including the lead strand edge

    def getpaths(self):
        if self.path is None:
            self._runfrominitialnode()
        else:
            try:
                self._findpaths()
                self.paths.append(self.path)
            except self.DeadEndPath:
                # self.path will be discarded.
                # Don't add it to self.paths, just return results that may
                # have been added from other RopeTrick objects
                pass
        return self.paths

    def _runfrominitialnode(self):
        assert len(self.initialnode.outgoingedges) > 0 and \
            len(self.initialnode.outgoingedges[0].strandedges) > 0, \
            'Start node should not be a dead end'
        for edge in  self.initialnode.outgoingedges:
            strandedges = edge.strandedges
            leadstrandedge = strandedges[0]
            path = Path([PathNode(self.initialnode, leadstrandedge)])
            # Spawn a new RopeTrick for each edge and
            # capture the results on the current RopeTrick object
            # This is a little simpler than searching one edge with the
            # current RopeTrick object since that would fork the code.
            rt = RopeTrick(self.initialnode, self.terminalnode,
                           path=path, strandedges=strandedges)
            self.paths.extend(rt.getpaths())

    def _findpaths(self):
        while self.terminalnode != self.leadstrandedge.tonode:
            self._step()
        # Add final node and end search
        self.path.append(self.terminalnode)
        return

    def _getleadstrandedge(self):
        return self.path.getcurrentleadstrandedge()

    def _step(self):
        # If lead strand has reached its end, switch strands
        if self.leadstrandedge.nextstrandedge is None:
            self._selectnewleadstrand()

        # Move everything to the next edge
        self.leadstrandedge = self.leadstrandedge.nextstrandedge
        self.strandedges = self._getnextstrandedges(self.strandedges)
        self.path.step(
            leadstrandedge=self.leadstrandedge,
            terminatedstrandedges=self._getterminatedstrandedges())

        # Also add any new strands that start along the path
        self.strandedges.extend(
            self._getstartingstrandedges(self._getcurrentnode()))

        self._processdivergentstrands()

    def _getterminatedstrandedges(self):
        return filter(
            lambda e: e.nextstrandedge==None,
            self.strandedges)

    def _getstartingstrandedges(self, node):
        # First edge of all strands beginning on node and not
        # in leadstrandhistory
        for edge in node.outgoingedges:
            return filter(
                lambda e: e.previousstrandedge is None and \
                e.strand not in self.path.leadstrandhistory,
                edge.strandedges)
                                
    def _getpreviousstrandedges(self, strandedges):
        # Get all previousstrandedges that are not None
        return filter(
            lambda e: e is not None,
            map(lambda e: e.previousstrandedge,
                strandedges))

    def _getnextstrandedges(self, strandedges):
        # Get all nextstrandedges that are not None
        return filter(
            lambda e: e is not None,
            map(lambda e: e.nextstrandedge,
                strandedges))

    def _processdivergentstrands(self):
        # Split any divergent strands into new RopeTricks,
        # one for each divergent edge
        for edge in self._getcurrentnode().outgoingedges:
            if edge.tonode == self.leadstrandedge.tonode:
                continue
            divergentnode = edge.tonode
            divergentstrandedges = filter(
                lambda e: e.tonode==divergentnode,
                self.strandedges)
            # Since the current lead strand is not on the path of these
            # divergent edges, we need to create new RopeTricks and rewind
            # to where the lead strand started, then search again from there.
            path = self.path.clone()
            rt = RopeTrick(self.initialnode, self.terminalnode,
                           path=path, strandedges=divergentstrandedges)
            self.paths.extend(rt.rewindandgetpaths())

        # Reset strandedges to contain only non-divergent strands
        self.strandedges = filter(
            lambda e: e.tonode==self.leadstrandedge.edge.tonode,
            self.strandedges)

    def _getcurrentnode(self):
        return self.path.getcurrentnode()

    def rewindandgetpaths(self):
        try:
            self._rewind()
        except self.DeadEndPath:
            # Give up, return any accumulated results
            return self.paths
        return self.getpaths()

    def _rewind(self):
        # Rewind to the start of the last leadstrand.
        # Then we can search again to find a path independent of that
        # leadstrand.
        strandtoremove = self.path.getcurrentleadstrand()
        while self.leadstrandedge is not None:
            leadstrandfirstedge = self.leadstrandedge
            self.leadstrandedge = self.leadstrandedge.previousstrandedge
            self.strandedges = self._getpreviousstrandedges(self.strandedges)
            self.strandedges.extend(self.path.getcurrentterminatedstrandedges())
            pathnode = self.path.pop()

        # While backgracking we might have picked up some strands
        # that didn't follow back to this point. Filter them out.
        self.strandedges = filter(
            lambda e: e.tonode == leadstrandfirstedge.tonode,
            self.strandedges)
        self._selectnewleadstrand()
        # replace the last pathnode to use this new lead strand
        #pathnode = self.path.pop()
        self.path.step(
            leadstrandedge=self.leadstrandedge,            
            terminatedstrandedges=self._getterminatedstrandedges())

    def _selectnewleadstrand(self):
        for strandedge in self.strandedges:
            if strandedge.strand not in self.path.leadstrandhistory and\
               strandedge.nextstrandedge is not None:
                self.leadstrandedge = strandedge
                return
        # No eligible lead strands to continue.
        raise self.DeadEndPath


class PathNode:
    def __init__(self, dbgraphnode, leadstrandedge=None,
                 terminatedstrandedges=None):
        self.node = dbgraphnode
        # leadstrandedge that leaves from the current node
        self.leadstrandedge = leadstrandedge

        # strandedges of any strands that terminated on the current node
        if terminatedstrandedges is None:
            self.terminatedstrandedges = []
        else:
            self.terminatedstrandedges = terminatedstrandedges


class Path:
    def __init__(self, pathnodeslist=None):
        if pathnodeslist:
            self.pathnodes = pathnodeslist
        else:
            self.pathnodes = []

        # initialize leadstrandhistory
        self.leadstrandhistory = set()
        for pathnode in self.pathnodes:
            if pathnode.leadstrandedge:
                self.leadstrandhistory.add(
                    pathnode.leadstrandedge.strand)

    def append(self, dbgraphnode, leadstrandedge=None,
               terminatedstrandedges=None):
        if leadstrandedge:
            self.leadstrandedges.add(leadstrandedge.strand)
        self.pathnodes.append(
            PathNode(
                dbgraphnode,
                leadstrandedge=leadstrandedge,
                terminatedstrandedges=terminatedstrandedges))
        
    def clone(self):
        path = Path(copy.copy(self.pathnodes))
        path.leadstrandhistory = self.leadstrandhistory
        return path

    def step(self, leadstrandedge=None, terminatedstrandedges=None):
        if self.pathnodes:
            assert self.pathnodes[-1].leadstrandedge.tonode == \
            leadstrandedge.fromnode, \
            'leadstrandedge does not start on current path node'
        self.leadstrandhistory.add(leadstrandedge.strand)
        self.pathnodes.append(
            PathNode(leadstrandedge.fromnode,
                 leadstrandedge=leadstrandedge,
                 terminatedstrandedges=terminatedstrandedges))

    def pop(self):
        return self.pathnodes.pop()

    def getcurrentnode(self):
        if not self.pathnodes:
            return None
        return self.pathnodes[-1].node

    def getcurrentleadstrandedge(self):
        if not self.pathnodes:
            return None
        return self.pathnodes[-1].leadstrandedge

    def getcurrentleadstrand(self):
        lastleadstrandedge = self.getcurrentleadstrandedge()
        if lastleadstrandedge is None:
            return None
        return lastleadstrandedge.strand

    def getcurrentterminatedstrandedges(self):
        if not self.pathnodes:
            return None
        return self.pathnodes[-1].terminatedstrandedges

    def getlength(self):
        return len(self.pathnodes)

    def __str__(self):
        if not self.pathnodes:
            return "empty path"
        seq = self.pathnodes[0].node.kminus1mer
        for pathnode in self.pathnodes[1:]:
            seq += pathnode.node.kminus1mer[-1]
        if self.pathnodes[-1].leadstrandedge:
            seq += self.pathnodes[-1].leadstrandedge.getsequenceletter()
        return seq

    def __repr__(self):
        return "<Path: %s>" % str(self)


if __name__=='__main__':
    ref = 'agagtcgctacttatcgtgtctaactaatgggtacttcagctcaagagat'
    # MNP complicates the graph by creating a 2nd loop to tactt
    altseq = 'agagtcgctacttatcgtgtctacttaatgggtacttcagctcaagagat'
    # splits on first step
    altseq1 = 'agagtcactacttatcgtgtctaactaatgggtacttcagctcaagagat'
    # splits on second step
    altseq2 = 'agagtcgatacttatcgtgtctaactaatgggtacttcagctcaagagat'
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

    graph = Graph(ref, kmersize=kmersize)
    graph.saveimage('ref')
    graph.addsequences(readidpairs, sampleid=sampleid)
    graph.saveimage('test')
    print graph.getcontinuouspaths()
