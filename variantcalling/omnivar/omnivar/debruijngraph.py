import pydot


DEFAULT_KMER_SIZE = 21


class Node(object):

    """A Node represents a (k-1)mer
    """

    def __init__(self, k_minus_1mer):
        self.k_minus_1mer = k_minus_1mer
        self.incoming_edges = set()
        self.outgoing_edges = set()


class Edge(object):

    """An Edge represents a kmer
    """

    def __init__(self, from_node, to_node):
        self.from_node = from_node
        self.to_node = to_node
        self.from_node.outgoing_edges.add(self)
        self.to_node.incoming_edges.add(self)
        self.weight = 0

    def get_sequence(self):
        return self.from_node.k_minus_1mer + self.get_letter()

    def get_letter(self):
        return self.to_node.k_minus_1mer[-1]

    def get_next_edges(self, minimum_weight=None):
        edges = self.to_node.outgoing_edges
        if minimum_weight is not None:
            edges = filter(lambda e: e.weight >= minimum_weight, edges)
        return edges

    def addWeight(weight=1):
        self.weight += 1


class Graph(object):

    """Graph is a DeBruijn Graph constructed from sequence data
    including short reads from one or more samples and a reference 
    sequence
    """

    def __init__(self, reference_sequence, kmer_size=DEFAULT_KMER_SIZE):
        if len(reference_sequence) < kmer_size:
            raise ValueError('Length of reference sequence is less than kmer size %s' % kmer_size)
        self.kmer_size = kmer_size
        self.nodes = {}
        self.edges = {}
        self.add_reference_sequence(reference_sequence)

    def get_nodes(self):
        return set(self.nodes.values())

    def get_edges(self):
        return set(self.edges.values())

    def add_reference_sequence(self, sequence):
        self.add_sequence(sequence, has_weight=False)
        self.first_edge = self.edges[sequence[0:self.kmer_size]]
        self.last_edge = self.edges[sequence[-self.kmer_size:]]

    def add_sequences(self, sequence_list):
        for sequence in sequence_list:
            self.add_sequence(sequence)

    def add_sequence(self, sequence, has_weight=True):
        if len(sequence) < self.kmer_size:
            return
        for kmer in self._kmer_generator(sequence):
            edge = self._get_or_create_edge(kmer)
            if has_weight:
                edge.weight += 1

    def _get_or_create_edge(self, kmer):
        try:
            edge = self.edges[kmer]
        except KeyError:
            from_node = self._get_or_create_node(kmer[:-1])
            to_node = self._get_or_create_node(kmer[1:])
            edge = Edge(from_node, to_node)
            self.edges[kmer] = edge
        return edge

    def _get_or_create_node(self, k_minus_1mer):
        try:
            node = self.nodes[k_minus_1mer]
        except KeyError:
            node = Node(k_minus_1mer)
            self.nodes[k_minus_1mer] = node
        return node

    def _kmer_generator(self, sequence):
        # chop sequence into kmers
        for i in xrange(0, len(sequence)-(self.kmer_size-1)):
            yield sequence[i:i+self.kmer_size]

    def save_image(self, filename_root):
        dotfile = filename_root + '.dot'
        with open(dotfile, 'w') as f:
            f.write("digraph \"Graph\" {\n")
            f.write("  bgcolor=\"white\";\n")            
            for node in self.get_nodes():
                f.write('  %s [label="%s"] ;\n' % (
                    node.k_minus_1mer, node.k_minus_1mer))

            for edge in self.get_edges():
                f.write(
                    '  %s -> %s [label="%d" color="red"] ;\n' % (
                        edge.from_node.k_minus_1mer,
                        edge.to_node.k_minus_1mer,
                        edge.weight,
                ))
            f.write("}\n")

        import pydot
        png_file = filename_root + '.png'
        (graph,) = pydot.graph_from_dot_file(dotfile)
        graph.write_png(png_file)
        print "saved png file %s" % png_file
