import copy
import math
from .debruijngraph import Graph


DEFAULT_MINIMUM_WEIGHT = 3


class Path(object):

    minimum_weight = DEFAULT_MINIMUM_WEIGHT

    class LowPathWeight(Exception):
        pass

    def __init__(self, edge=None):
        self.edge_list = []
        self.edge_visits = {}
        if edge:
            self.add_edge(edge)

    def get_length(self):
        return len(self.edge_list)

    def get_weights(self):
        return [edge.weight for edge in self.edge_list]

    def get_sequence(self):
        if len(self.edge_list) == 0:
            return ''
        sequence = self.edge_list[0].get_sequence()
        for edge in self.edge_list[1:]:
            sequence += edge.get_letter()
        return sequence

    def clone(self):
        clone = Path()
        clone.edge_list = copy.copy(self.edge_list)
        clone.edge_visits = copy.copy(self.edge_visits)
        return clone

    def get_last_edge(self):
        return self.edge_list[-1]

    def ends_with(self, edge):
        return edge is self.get_last_edge()

    def add_edge(self, edge):
        visits = self.edge_visits.setdefault(edge, 0) + 1
        print visits
        self.edge_visits[edge] = visits
        if float(edge.weight) / visits < self.minimum_weight:
            raise self.LowPathWeight()
        self.edge_list.append(edge)


class PathFinder(object):

    def __init__(self, graph, minimum_weight=DEFAULT_MINIMUM_WEIGHT):
        self.graph = graph
        self.paths = []
        Path.minimum_weight = minimum_weight

    def find_paths(self):
        try:
            path = self._initialize_path()
        except Path.LowPathWeight:
            return self.paths
        self._dfs(path)
        return self.paths

    def _initialize_path(self):
        path = Path()
        self._step_once(path, self.graph.first_edge)
        return path

    def _dfs(self, path):
        this_edge = path.get_last_edge()
        next_edges = this_edge.get_next_edges()
        while len(next_edges) == 1:
            next_edge = list(next_edges)[0]
            try:
                self._step_once(path, next_edge)
            except Path.LowPathWeight:
                return
            next_edges = next_edge.get_next_edges(minimum_weight=path.minimum_weight)
        if len(next_edges) == 0:
            # dead end
            return
        for next_edge in next_edges:
            new_path = path.clone()
            try:
                new_path.add_edge(next_edge)
            except Path.LowPathWeight:
                continue
            self._dfs(new_path)

    def _step_once(self, path, edge):
        path.add_edge(edge)
        if path.ends_with(self.graph.last_edge):
            self.paths.append(path)
            print path
