#!/usr/bin/env python

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
        # A mechanism to disable a node
        # that has been split
        self.isdead = False
        self.nodeid = self.idcounter
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

    def switchtonode(self, newnode):
        for strandedge in self.strandedges:
            strandedge.switchtonode(newnode)
        self.tonode.incomingedges.remove(self)
        newnode.incomingedges.add(self)
        self.tonode = newnode

    def switchfromnode(self, newnode):
        for strandedge in self.strandedges:
            strandedge.switchfromnode(newnode)
        self.fromnode.outgoingedges.remove(self)
        newnode.outgoingedges.add(self)
        self.fromnode = newnode

    def getnextstrandedges(self):
        nextstrandedges = set()
        map(lambda e: nextstrandedges.update(e.nextstrandedges),
            self.strandedges)
        return nextstrandedges

    def getnextstrandedgesonnextedge(self, edge):
        return set(filter(lambda e: e.edge == edge,
                          self.getnextstrandedges()))
        
    def splitbyfromnode(self, fromnodes):
        # This is called after the strandedges in this edge
        # have been split, so their fromnodes no longer match
        # the edge, and there may be different fromnodes for
        # different strandedges.
        # We need to replace this edge with daughter edges
        # that match node connections of strandedges.
        newstrandedges = self.fromnode.getoriginatingstrandedges()
        daughteredges = set()
        for fromnode in fromnodes:
            # In general nodes can have many incoming edges,
            # but these ones should only have one.
            assert len(fromnode.incomingedges) == 1, \
                'Expected just 1 incoming edge after split'
            fromedge = list(fromnode.incomingedges)[0]
            strandedges = fromedge.getnextstrandedgesonnextedge(self)
            if len(strandedges) > 0:
                # Add in strands that originate on the old fromnode.
                # Since these could correspond to any of the new nodes,
                # we add them to all.
                strandedges.update([e.strand.clone(e) for e in newstrandedges])
                for strandedge in strandedges:
                    strandedge.switchfromnode(fromnode)
                daughteredge = Edge(fromnode, self.tonode)
                daughteredge.strandedges = strandedges
                self.tonode.incomingedges.add(daughteredge)
                daughteredges.add(daughteredge)
        if len(daughteredges) == 0:
            # Special case: None of the fromnodes had an outgoing edge,
            # so all strandedges here are new. In this case they have
            # no connection to the previous sequences, and we should
            # create one daughter with just these new reads.
            severednode = Node(self.fromnode.kminus1mer)
            for strandedge in newstrandedges:
                strandedge.switchfromnode(severednode)
            daughteredge = Edge(severednode, self.tonode)
            daughteredge.strandedges = newstrandedges
            self.tonode.incomingedges.add(daughteredge)
            daughteredges.add(daughteredge)
            print "SPECIAL CASE"
            import pdb; pdb.set_trace()
        else:
            # Remove unused strandedges, since these were replaced by clones
            for strandedge in newstrandedges:
                strandedge.strand.remove(strandedge)
        self.tonode.incomingedges.remove(self)
        return daughteredges


class Strand:
    """A Strand is a collection of StrandEdges that represent the path
    a single strand (or graphed sequence) takes through the graph.
    """
    idcounter = 0
    def __init__(self, strandid, sampleid):
        self.strandid = strandid
        self.sampleid = sampleid
        self.firststrandedges = set()
        self.laststrandedges = set()
        self.strandid = self.idcounter
        Strand.idcounter += 1

    def getfirstnode(self):
        if len(self.firststrandedges) > 1:
            raise Exception('More than one first nodes exist')
        if len(self.firststrandedges) == 0:
            return None
        return list(self.firststrandedges)[0].fromnode

    def getlastnode(self):
        if len(self.laststrandedges) > 1:
            raise Exception('More than one last nodes exit')
        if len(self.laststrandedges) == 0:
            return None
        return list(self.laststrandedges)[0].tonode

    def getsequence(self):
        if not self.firststrandedges:
            return None
        else:
            firststrandedge = list(self.firststrandedges)[0]
            seq = firststrandedge.fromnode.kminus1mer
            seq += firststrandedge.getsequencetoendofstrand()
            return seq

    def clone(self, strandedge):
        clone = StrandEdge(
            strandedge.edge, strandedge.strand)
        clone.nextstrandedges = copy.copy(strandedge.nextstrandedges)
        clone.previousstrandedges = copy.copy(strandedge.previousstrandedges)
        for _next in clone.nextstrandedges:
            _next.previousstrandedges.add(clone)
        for previous in clone.previousstrandedges:
            previous.nextstrandedges.add(clone)
        if strandedge in self.firststrandedges:
            self.firststrandedges.add(clone)
        if strandedge in self.laststrandedges:
            self.laststrandedges.add(clone)
        return clone

    def remove(self, strandedge):
        for previous in copy.copy(strandedge.previousstrandedges):
            previous.nextstrandedges.remove(strandedge)
            strandedge.previousstrandedges.remove(previous)
        for _next in copy.copy(strandedge.nextstrandedges):
            _next.previousstrandedges.remove(strandedge)
            strandedge.nextstrandedges.remove(_next)
        if strandedge in self.firststrandedges:
            self.firststrandedges.remove(strandedge)
        if strandedge in self.laststrandedges:
            self.laststrandedges.remove(strandedge)


class StrandEdge:
    """A StrandEdge strand represents the passage of 1 read along
    an edge. Collectively all the EdgeStrands for 1 read
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
        self.islinearized = False  # Until we call self.linearize()
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
        return self.nodes.values()

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
        assert not self.islinearized, 'Not allowed on linearized graph'
        strand = Strand(strandid, sampleid)
        self.strands[strandid] = strand
        previous = None
        for kmer in self.kmergenerator(sequence):
            previous = self.createstrandedgeinorder(kmer, strand, previous)
        strand.laststrandedges.add(previous)
        return strand

    def createstrandedgeinorder(self, kmer, strand, previous):
        # Create edge between left k-minus-1mer and right k-minus-1mer nodes
        assert not self.islinearized, 'Not allowed on linearized graph'
        edge = self.getorcreateedge(kmer)
        strandedge = StrandEdge(edge, strand)
        if not strand.firststrandedges:
            strand.firststrandedges.add(strandedge)
        edge.strandedges.add(strandedge)
        if previous:
            strandedge.linktoprevious(previous)
        return strandedge

    def getorcreateedge(self, kmer):
        assert not self.islinearized, 'Not allowed on linearized graph'
        try:
            edge = self.edges[kmer]
        except KeyError:
            fromnode = self.getorcreatenode(kmer[:-1])
            tonode = self.getorcreatenode(kmer[1:])
            edge = Edge(fromnode, tonode)
            self.edges[kmer] = edge
        return edge

    def getorcreatenode(self, kminus1mer):
        assert not self.islinearized, 'Not allowed on linearized graph'
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
        self.nodes = set(graph.nodes.values())
        self.edges = set(graph.edges.values())
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
        crawler = Crawler(self)
        if not self.firstnodes and not self.lastnodes:
            raise Exception(
                'No firstnodes/lastnodes set so there is nowhere to start.')
        for firstnode in copy.copy(self.firstnodes):
            crawler.crawlforward(firstnode)
        for lastnode in copy.copy(self.lastnodes):
            crawler.crawlbackward(lastnode)

    def getnodes(self):
        return self.nodes


class Crawler:
    counter = 1

    def __init__(self, linearizablegraph):
        self.graph = linearizablegraph

    def crawlforward(self, node):
        while len(node.outgoingedges) == 1 \
           and len(node.incomingedges) == 1 \
           and not node.isdead:
            # Nothing to do. Fast forward
            edge = list(node.outgoingedges)[0]
            node = edge.tonode

        self.graph.saveimage('before', node)
        print "New Cycle " + node.kminus1mer
        #node.printstrandsummary()

        if node.isdead:
            raise Exception("DEADSTOP")
            # Node has been split by another crawler,
            # so stop this crawler
            return

        if not node.incomingedges:
            # Special case: Crawler was initiated on
            # a start node with no incoming edges.
            # No need to split the node. Just crawl each
            # outgoing edge.
            if not node.outgoingedges:
                # Nothing to do. We started and ended
                # on a node with no paths.
                return
            # Crawl forward on the tonode for each edge.
            outgoingedges = node.outgoingedges
            for edge in outgoingedges:
                self.crawlforward(edge.tonode)
            return

        # This is the meat of the algorithm. Split the current
        # node into one daughternode for each incoming edge,
        # and transfer the edges to their new nodes. Recurse on
        # each outgoing edge

        # Easier to save these before we split the node.
        # We'll use them later.
        outgoingedges = copy.copy(node.outgoingedges)
        nextnodes = [edge.tonode for edge in outgoingedges]

        # Break the current node into one daughternode for each
        # incoming edge. Reassign the incoming and outgoing edges
        # to the new daughter nodes
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
    #readidpairs.extend(getreadsandids(altseq, readlength))
    #readidpairs.extend(getreadsandids(altseq1, readlength))
    #readidpairs.extend(getreadsandids(altseq2, readlength))
    #readidpairs.extend(getreadsandids(altseq3, readlength))
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
    
