#!/usr/bin/env python

import copy
import uuid

"""Unfinished work. Zipper is meant to linearize complex debruijn graphs
by "unzipping" cohesive bundles of strands--splitting them off of existing
nodes and edges to create a separate path--and then "zipping" them back
together with strands that lie partially along the unzipped path (merging
synonymous nodes and edges to combine their strands).

This works for simple cases but it is not clear to me how it can be adapted
to resolve all cases.
"""

REFERENCE_SAMPLE_ID = -1
REFERENCE_STRAND_ID = -1
DEFAULT_KMER_SIZE = 9

class NoCohesiveStrands(Exception):
    pass

class Node:
    """Node represents kminus1mer
    """
    idcounter = 0

    def __init__(self, kminus1mer):
        self.kminus1mer = kminus1mer
        self.incomingedges = set()
        self.outgoingedges = set()
        # A mechanism to disable a node
        # that has been split
        self.isdead = False
        self.nodeid = '%s_%s' % (self.idcounter, self.kminus1mer)
        Node.idcounter += 1

    def addincomingedge(self, edge):
        self.incomingedges.add(edge)

    def addoutgoingedge(self, edge):
        self.outgoingedges.add(edge)

    def getincomingstrandedges(self):
        strandedges = set()
        for edge in self.incomingedges:
            strandedges.update(edge.strandedges)
        return strandedges

    def getoutgoingstrandedges(self):
        strandedges = set()
        for edge in self.outgoingedges:
            strandedges.update(edge.strandedges)
        return strandedges

    def getoriginatingstrandedges(self):
        # Strand edges that start on this node
        return list(filter(
            lambda e: len(e.previousstrandedges)==0,
            self.getoutgoingstrandedges()))

    def getterminatingstrandedges(self):
        # Strand edges that end on this node
        return filter(
            lambda e: len(e.nextstrandedges)==0,
            self.getincomingstrandedges())

    def printstrandsummary(self):
        strands = set()
        for edge in self.incomingedges:
            for strandedge in edge.strandedges:
                strands.add(strandedge.strand)
        for strand in strands:
            print strand.getsequence()
        strands = set()
        for edge in self.outgoingedges:
            for strandedge in edge.strandedges:
                strands.add(strandedge.strand)
        for strand in strands:
            print strand.getsequence()


class Edge:
    """Edge represents kmer
    """
    def __init__(self, fromnode, tonode):
        self.fromnode = fromnode
        self.tonode = tonode
        self.fromnode.addoutgoingedge(self)
        self.tonode.addincomingedge(self)
        self.strandedges = set()

    def getweight(self):
        return len(self.strandedges)

    def getnextstrandedges(self):
        return filter(
            lambda e: e is not None,
            map(lambda e: e.nextstrandedge,
                self.strandedges))

    def getpreviousstrandedges(self):
        return filter(
            lambda e: e is not None,
            map(lambda e: e.previousstrandedge,
                self.strandedges))

    def splitforward(self):
        # If strands in this edge terminate or do not cohesively
        # continue onto another edge, raise NoCohesiveStrands
        # Otherwise split away from self.tonode onto a newly
        # created node and return the next edge traversed by
        # the strands.
        print "TODO"


class Strand:
    """A Strand is a collection of StrandEdges that represent the path
    a single strand (or graphed sequence) takes through the graph.
    """
    idcounter = 0
    def __init__(self, strandid, sampleid):
        self.strandid = strandid
        self.sampleid = sampleid
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

    def getsequence(self):
        if not self.firststrandedge:
            return None
        else:
            seq = self.getfirstnode().kminus1mer
            seq += self.firststrandedge.getsequencetoendofstrand()
            return seq


class StrandEdge:
    """A StrandEdge strand represents the passage of 1 read along
    an edge. Collectively all the StrandEdges for 1 read
    represent the path that read takes through the graph.
    """
    def __init__(self, edge, strand):
        self.fromnode = edge.fromnode
        self.tonode = edge.tonode
        self.nextstrandedges = set()
        self.previousstrandedges = set()
        self.strand = strand
        self.edge = edge

    def linktoprevious(self, previousstrandedge):
        self.previousstrandedges.add(previousstrandedge)
        previousstrandedge.nextstrandedges.add(self)

    def switchtonode(self, newnode):
        self.tonode = newnode

    def switchfromnode(self, newnode):
        self.fromnode = newnode

    def getsequencetoendofstrand(self):
        seq = self.getsequenceletter()
        # If this strand is forked, all paths should give the
        # same result, so just pick one at random
        if self.nextstrandedges:
            nextstrandedge = list(self.nextstrandedges)[0]
            seq += nextstrandedge.getsequencetoendofstrand()
        return seq

    def getsequenceletter(self):
        return self.tonode.kminus1mer[-1]


class _AbstractGraph:
    """For functions shared by Graph and LinearizedGraph
    """

    def saveimage(self, filenameroot, boldnode=None):
        dotfile = filenameroot + '.dot'
        with open(dotfile, 'w') as f:
            f.write("digraph \"Graph\" {\n")
            f.write("  bgcolor=\"white\";\n")
            for node in self.getnodes():
                emphasis = ''
                if node == boldnode:
                    emphasis = '@@@'
                kminus1mer = node.kminus1mer
                label = kminus1mer
		f.write("  %s [label=\"%s%s%s\"] ;\n" % (
                    node.nodeid, emphasis, label, emphasis))
            for node in self.getnodes():
                for edge in node.outgoingedges:
                    for strandedge in edge.strandedges:
                        color = strandedge.strand.strandid % 9 + 1
                        f.write("  %s -> %s [label=\"%d\" color=\"%s\" colorscheme=\"pastel19\"] ;\n" % (
                            strandedge.fromnode.nodeid, strandedge.tonode.nodeid,
                            strandedge.strand.strandid, color))
            f.write("}\n")

        import pydot
        pngfile = filenameroot + '.png'
        (graph,) = pydot.graph_from_dot_file(dotfile)
        graph.write_png(pngfile)
        print "saved png file %s" % pngfile


class Graph(_AbstractGraph):
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
        for kmer in self.kmergenerator(sequence):
            previous = self.createstrandedgeinorder(kmer, strand, previous)
        strand.laststrandedge = previous
        return strand

    def createstrandedgeinorder(self, kmer, strand, previous):
        # Create edge between left k-minus-1mer and right k-minus-1mer nodes
        edge = self.getorcreateedge(kmer)
        strandedge = StrandEdge(edge, strand)
        if not strand.firststrandedge:
            strand.firststrandedge = strandedge
        edge.strandedges.add(strandedge)
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

    def getlinearizedgraph(self):
        graph = LinearizableGraph(self)
        graph.linearize()
        return graph


class LinearizableGraph(_AbstractGraph):
    # LinearizableGraph allows duplicate nodes and edges
    # with the same kminus1mer/kmer, so dicts that use
    # sequences as index won't work. Instead we store
    # nodes and edges in sets.

    def __init__(self, graph):
        self.strands = graph.strands
        self.nodes = graph.getnodes()
        self.edges = graph.getedges()
        # first/lastnodes are sets because the initial
        # first/lastnode may later be split into multiple
        self.firstnodes = set([])
        self.lastnodes = set([])
        firstnode = graph.getfirstnode()
        lastnode = graph.getlastnode()
        if firstnode:
            self.firstnodes.add(firstnode)
        if lastnode:
            self.firstnodes.add(lastnode)
        self.referencestrand = graph.referencestrand

    def linearize(self):
        zipper = Zipper(self)
        if not self.firstnodes and not self.lastnodes:
            raise Exception(
                'No firstnodes/lastnodes set so there is nowhere to start.')
        for node in copy.copy(self.firstnodes):
            zipper.processforward(node)
        # for node in copy.copy(self.lastnodes):
        #     zipper.processbackward(node)

    def getnodes(self):
        return self.nodes


class Zipper:

    def __init__(self, linearizablegraph):
        self.graph = linearizablegraph

    def processforward(self, node):
        assert not node.isdead, "Attempted to process a dead node"

        if not node.incomingedges:
            if not node.outgoingedges:
                # Nothing to do. We started and ended
                # on a node with no incoming or outgoing edges paths.
                return

        # Fast-forward as long as there is nothing to do
        while and len(node.outgoingedges) == 1 \
              and len(node.incomingedges) <= 1 \
              and not node.isdead:
            # Only 1 edge. Follow it to next node.
            node = list(node.outgoingedges)[0].tonode

        # When we reach a merging path, run unzip on the node, and
        # then continue scan on any split nodes produced.
        if len(node.incomingedges) > 1:
            splitnodes = self.forwardzipcycle(node)
            for splitnode in splitnodes:
                self.processforward(splitnode)
            return

        # When we reach a fork (if there is no merge on the same node),
        # processforward on each forked path.
        assert len(outgoingedges) > 1, 'Expected a fork but found none'
        for edge in copy.copy(self.outgoingedges):
            self.processforward(edge.tonode)
        return

    def forwardzipcycle(self, node):
        endpoints = []
        for edge in copy.copy(node.incomingedges):
            lastedge, stepcount = self.unzip(edge)
            endpoints.append((lastedge, stepcount))
        splitnodes = self.zip(endpoints)
        return splitnodes

    def unzip(self, edge, count=0):
        # unzip() will split "edge" off from the current "edge.tonode", and move it
        # to a newly created node that has only "edge" in its incoming edges.
        # unzip() recurses on "nextedge" (an edge that all strands in "edge" continue onto)
        # until none can befound, either because the strands end or because the strands split
        # onto different paths.

        while True:
            try:
                nextedge = edge.splitforward()
                count += 1
            except NoCohesiveStrands:
                # We've reached the end the strands in edge, or the strands have split apart.
                return edge, count

    def zip(self, endpoints):
        combinedendpoints = self.inclusivezip(endpoints)
        splitnodes = self.exclusivezip(combinedendpoints)
        return splitnodes

    def inclusivezip(self, endpoints):
        # endpoints is a list of (edge, count) tuples that represent
        # the last edge reached by each unzip operation and the number
        # of edges unzipped from the merge point to reach it.
        # Some of the endpoints may have paths to the merge node that overlap
        # partially or entirely with each other. inclusivezip() merges those
        # paths to the extent possible and returns a new set of combinedendpoints,
        # which are the most distant unzipped edges that lie on distinct paths.
        # One combinedendpoint may represent multiple endpoints if they lie on
        # a merged path.
        pass

    def exclusivezip(combinedendpoints):
        # combinedendpoints is a set of edges that represents the furthest points
        # reached 
        pass


    """
    def crawlforward(self, node):

        daughternodes = set()
        for edge in copy.copy(node.incomingedges):
            daughternode = Node(node.kminus1mer)
            edge.switchtonode(daughternode)
            self.graph.nodes.add(daughternode)
            daughternodes.add(daughternode)
        self.graph.nodes.remove(node)

        # If we split a first/last node, update the graph
        # with the new nodes
        if node in self.graph.firstnodes:
            self.graph.firstnodes.remove(node)
            self.graph.firstnodes.update(daughternodes)
        if node in self.graph.lastnodes:
            self.graph.lastnodes.remove(node)
            self.graph.lastnodes.update(daughternodes)

        # Kill this node so no other crawler will attempt to process it.
        node.isdead = True

        if not outgoingedges:
            # This path ended. Nothing else to do.
            return

        # At this point outgoingedges are contaminated with strandedges
        # from multiple daughter nodes. We need to split them.
        for edge in outgoingedges:
            daughteredges = edge.splitbyfromnode(daughternodes)
            self.graph.edges.remove(edge)
            self.graph.edges.update(daughteredges)

        #for node in daughternodes:
        #    node.printstrandsummary()
        self.graph.saveimage('after', node)
        import pdb; pdb.set_trace()

        # This crawler will run on the first nextnode
        # Start a new crawler on all the other nextnodes
        for nextnode in nextnodes:
            self.crawlforward(nextnode)
        return
    """

    def crawlbackward(self, node):
        # Recursively regress one step until done
        pass


if __name__=='__main__':
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

    #readidpairs = [
    #    (str(uuid.uuid4()), ref),
    #    (str(uuid.uuid4()), altseq3)
    #]
    graph = Graph(ref, kmersize=kmersize)
    graph.saveimage('ref')
    graph.addsequences(readidpairs, sampleid=sampleid)
    graph.saveimage('test')
    lg = LinearizableGraph(graph)
#        graph.nodes.values(), graph.edges.values(), graph.strands, referencestrand=graph.referencestrand)
    lg.saveimage('linearizable')
    lg.linearize()
    lg.saveimage('linearized')
    
