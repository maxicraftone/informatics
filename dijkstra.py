#!/usr/bin/python
import math
import time
import sys

class Node:
    def __init__(self, name: str) -> None:
        """
        Node class. Just holds the name of the node

        :param name: Name/ Title of the node (only for representation)
        """
        self.name = name

    def __repr__(self) -> str:
        return self.name

    def __gt__(self, node2) -> bool:
        return self.name > node2.name


class Edge:
    def __init__(self, node1: Node, node2: Node, distance: int, direction: int) -> None:
        """
        Edge class for connecting nodes through edges with distances and
        directions

        :param node1: First node
        :param node2: Second node
        :param distance: Weight of the edge
        :param direction: Direction of the connection:
                0 = both,
                1 = in direction of the first node,
                2 = in direction of the second node
        """
        self.node1 = node1
        self.node2 = node2
        self.distance = distance
        self.direction = direction

    def facing(self) -> tuple:
        """
        Returns the node(s) being faced by the edge

        :return: Nodes
        """
        # If facing both nodes
        if self.direction == 0:
            return (self.node1, self.node2)
        # If facing the first node
        elif self.direction == 1:
            return (self.node1,)
        # If facing the second node
        elif self.direction == 2:
            return (self.node2,)

    def set_node1(self, node1: Node) -> None:
        """
        Set the node1 parameter of the edge

        :param node1: Value for node1
        """
        self.node1 = node1

    def set_node2(self, node2: Node) -> None:
        """
        Set the node2 parameter of the edge

        :param node2: Value for node2
        """
        self.node2 = node2

    def set_distance(self, distance: int) -> None:
        """
        Set the distance parameter of the edge

        :param distance: Value for distance/ weight of the edge
        """
        self.distance = distance

    def set_direction(self, direction: int) -> None:
        """
        Set the direction parameter of the edge

        :param direction: Value for direction of the edge:
                0 = both directions,
                1 = in direction of the first node,
                2 = in direction of the second node
        """
        self.direction = direction

    def is_first(self, node: Node) -> bool:
        """
        Checks whether the given node is node1

        :param node: Given node
        :return: True if node == node1; False if node != node1
        """
        return node == self.node1

    def is_second(self, node: Node) -> bool:
        """
        Checks whether the given node is node2

        :param node: Given node
        :return: True if node == node2; False if node != node2
        """
        return node == self.node2

    def __repr__(self) -> str:
        # If the edge is facing both nodes
        if self.direction == 0:
            return self.node1.name + ' <-----> ' + self.node2.name
        # If the edge is facing the first node
        elif self.direction == 1:
            return self.node1.name + ' <------ ' + self.node2.name
        # If the edge is facing the second node
        elif self.direction == 2:
            return self.node1.name + ' ------> ' + self.node2.name


class Graph:
    def __init__(self) -> None:
        """
        Graph class containing nodes and edges with methods of path finding
        """
        # Define nodes and edges lists
        self.nodes = []
        self.edges = []

    def add_node(self, node: Node) -> None:
        """
        Add node to graph

        :param node: Node to be added
        """
        self.nodes.append(node)

    def add_nodes(self, nodes: list) -> None:
        """
        Add nodes to graph

        :param nodes: List of nodes to be added
        """
        self.nodes.extend(nodes)

    def set_nodes(self, nodes: list) -> None:
        """
        Set nodes of the graph

        :param nodes: List of nodes to be set
        """
        self.nodes = nodes

    def add_edge(self, edge: Edge) -> None:
        """
        Add a new edge to the graph

        :param edge: Edge to be added
        """
        self.edges.append(edge)

    def add_edges(self, edges: list) -> None:
        """
        Add edges to graph

        :param edges: List of edges to be added
        """
        self.edges.extend(edges)

    def set_edges(self, edges: list) -> None:
        """
        Set edges of the graph

        :param nodes: List of edges to be set
        """
        self.edges = edges

    def get_neighbors(self, node: Node) -> dict:
        """
        Get all neighbors of a node, with the corresponding edges

        :param node: Node to get neighbors of
        :return: Dictionary with {<node>: <edge>, ...}
        """
        neighbors = {}
        # For all nodes
        if node in self.nodes:
            # For all edges
            for e in self.edges:
                # If the node is the first in the edge declaration
                if e.is_first(node):
                    # Add the other node to the dictionary of neighbors
                    neighbors[e.node2] = e
                # If the node is the second in the edge declaration
                elif e.is_second(node):
                    # Add the other node to the dictionary of neighbors
                    neighbors[e.node1] = e
        return neighbors


    def find_path(self, start: Node, end: Node) -> list:
        """
        Find the shortest path from start node to end node using the
        Dijkstra algorithm

        :param start: Node to start from
        :param end: Node to end at
        :return: Path from start to end [[<nodes_passed>], <length>]
        """
        if start in self.nodes and end in self.nodes:
            # Init path
            path = {}
            # Set all node distances to infinity and the start distance to 0
            for node in self.nodes:
                path[node] = [math.inf, None]
            path[start] = [0, None]

            # Calc path
            ##path = self.calc_paths(start, path)
            path = self.dijkstra(start, path)

            # Create list of passed notes in the path
            path_nodes = [end]

            n = end
            while n != start:
                path_nodes.append(path[n][1])
                n = path[n][1]
                print(n)

            path_nodes.reverse()

            # Return list with passed notes and length of the path
            return [path_nodes, path[end][0]]

        elif start not in self.nodes:
            print('Start node not in graph')

        elif end not in self.nodes:
            print('End node not in graph')

    def dijkstra(self, start: Node, path: dict) -> dict:
        nodes = path
        distances = {}
        while len(nodes) > 0:
            nodes = self.sort_dict(nodes)
            u = list(nodes.keys())[0]
            nodes.pop(u)
            neighbors = self.get_neighbors(u)
            for n in neighbors:
                if n in nodes:
                    if u in distances:
                        nodes[n] = [neighbors[n].distance + distances[u][0], u]
                        distances[n] = nodes[n]
                    else:
                        nodes[n] = [neighbors[n].distance, u]
                        distances[n] = nodes[n]

        return distances

    @staticmethod
    def sort_dict(dictionary: dict) -> dict:
        """
        Function used to sort a dictionary by its values

        :param dictionary: The dictionary to be used
        :return: The sorted dictionary
        """
        return {k: v for k, v in sorted(dictionary.items(), key=lambda item: item[1])}

def benchmark(function) -> None:
    time_before = time.time_ns()
    function()
    time_delta = (time.time_ns() - time_before)
    print('Benchmark: ' + str(time_delta/1000000) + ' ms')

if __name__ == '__main__':
    graph = Graph()

    mue = Node('München')
    frb = Node('Freiburg')
    stg = Node('Stuttgart')
    fkf = Node('Frankfurt')
    kbl = Node('Koblenz')
    kln = Node('Köln')
    dsd = Node('Düsseldorf')
    brm = Node('Bremen')
    hmb = Node('Hamburg')
    kil = Node('Kiel')
    swr = Node('Schwerin')
    bln = Node('Berlin')
    drs = Node('Dresden')
    lpg = Node('Leipzig')
    eft = Node('Erfurt')
    hnv = Node('Hannover')
    mgd = Node('Magdeburg')

    #a = Node('A')
    #b = Node('B')
    #c = Node('C')
    #d = Node('D')
    #e = Node('E')
    #f = Node('F')
    #g = Node('G')
    #h = Node('H')
    #i = Node('I')
    #j = Node('J')
    #k = Node('K')
    #l = Node('L')
    #m = Node('M')
    #n = Node('N')

    graph.set_nodes([mue, frb, stg, fkf, kbl, kln, dsd, brm, hmb, kil, swr, bln, drs, lpg, eft, hnv, mgd])
    graph.set_edges([
                Edge(mue, frb, 280, 1), Edge(frb, stg, 140, 1), Edge(stg, mue, 240, 1),
                Edge(stg, fkf, 100, 1), Edge(fkf, kbl, 70, 1), Edge(kln, kbl, 70, 1),
                Edge(kln, dsd, 70, 2), Edge(kbl, dsd, 150, 1), Edge(kbl, eft, 120, 2),
                Edge(dsd, eft, 160, 1), Edge(fkf, eft, 90, 1), Edge(fkf, mue, 290, 2),
                Edge(eft, mue, 320, 1), Edge(lpg, mue, 530, 2), Edge(lpg, eft, 300, 2),
                Edge(dsd, lpg, 450, 1), Edge(dsd, brm, 230, 2), Edge(brm, eft, 330, 2),
                Edge(hmb, brm, 130, 2), Edge(kil, hmb, 130, 1), Edge(swr, kil, 320, 1),
                Edge(bln, swr, 190, 1), Edge(bln, lpg, 160, 2), Edge(drs, lpg, 70, 2),
                Edge(drs, bln, 150, 2), Edge(lpg, mgd, 150, 2), Edge(mgd, hnv, 100, 2),
                Edge(hnv, bln, 190, 0), Edge(hmb, bln, 220, 0), Edge(eft, hnv, 300, 2),
                Edge(dsd, hnv, 270, 2), Edge(hnv, brm, 130, 2), Edge(hnv, hmb, 140, 2),

                #33 Edges
    ])

    #graph.set_nodes([a, b, c, d, e, f, g, h, i, j, k, l, m, n])
    #graph.set_edges([
    #            Edge(a, d, 4, 1), Edge(a, c, 2, 0), Edge(b, d, 1, 0),
    #            Edge(b, e, 2, 2), Edge(c, d, 5, 0), Edge(c, f, 6, 0),
    #            Edge(c, g, 4, 1), Edge(d, e, 3, 0), Edge(d, k, 5, 2),
    #            Edge(e, h, 2, 2), Edge(e, i, 3, 0), Edge(f, g, 3, 2),
    #            Edge(f, l, 4, 1), Edge(g, j, 4, 2), Edge(g, k, 3, 0),
    #            Edge(h, i, 2, 0), Edge(h, n, 4, 0), Edge(j, m, 2, 0),
    #            Edge(k, m, 1, 0), Edge(k, n, 3, 0), Edge(l, m, 3, 0)])

    nodes = {}
    for node in graph.nodes:
        nodes[node.name.lower()] = node

    args = sys.argv[1:]
    if len(args) == 2:
        node1 = nodes[args[0].lower()]
        node2 = nodes[args[1].lower()]
        benchmark(lambda: print(graph.find_path(node1, node2)))

    # benchmark(lambda: print(graph.find_path(h, f)))
    #path = graph.find_path(h, f)
    #print(path)
