#!/usr/bin/python
import math

class Node:
    def __init__(self, name: str) -> None:
        """
        Initialization

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
        Initialization

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
        if self.direction == 0:
            return (self.node1, self.node2)
        elif self.direction == 1:
            return (self.node1,)
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
        if self.direction == 0:
            return self.node1.name + ' <-----> ' + self.node2.name
        elif self.direction == 1:
            return self.node1.name + ' <------ ' + self.node2.name
        elif self.direction == 2:
            return self.node1.name + ' ------> ' + self.node2.name


class Graph:
    def __init__(self) -> None:
        """
        Initialization
        """
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
        if node in self.nodes:
            for e in self.edges:
                if e.is_first(node):
                    neighbors[e.node2] = e
                elif e.is_second(node):
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
            for node in self.nodes:
                path[node] = [math.inf, None]
            path[start] = [0, None]

            # Calc path
            path = self.calc_paths(start, path)

            # Create list of passed notes in the path
            path_nodes = [end]

            n = end
            while n != start:
                path_nodes.append(path[n][1])
                n = path[n][1]

            path_nodes.reverse()

            # Return list with passed notes and length of the path
            return [path_nodes, path[end][0]]

        elif start not in self.nodes:
            print('Start node not in graph')

        elif end not in self.nodes:
            print('End node not in graph')

    def calc_paths(self, prev_node: Node, path: dict) -> dict:
        """
        Recursive function to calculate the path from the previous node to the
        nearest next one

        :param prev_node: Previous node in the path
        :param path: The path so far
        :return: Dictionary with nodes and weights {<node>: [<distance>, <prev_node>], ...}
        """
        neighbors = self.get_neighbors(prev_node)
        for n in neighbors:
            if n in neighbors[n].facing():
                if path[n][0] > neighbors[n].distance + path[prev_node][0]:
                    path[n] = [neighbors[n].distance + path[prev_node][0], prev_node]
                    path = self.sort_dict(path)

        for n in neighbors:
            if n != path[prev_node][1] and path[n][0] == neighbors[n].distance + path[prev_node][0]:
                # Calc further path with new nearest node
                path = self.calc_paths(n, path)

        return path

    @staticmethod
    def sort_dict(dictionary: dict) -> dict:
        """
        Function used to sort a dictionary by its values

        :param dictionary: The dictionary to be used
        :return: The sorted dictionary
        """
        return {k: v for k, v in sorted(dictionary.items(), key=lambda item: item[1])}


if __name__ == '__main__':
    graph = Graph()
    a = Node('A')
    b = Node('B')
    c = Node('C')
    d = Node('D')
    e = Node('E')
    f = Node('F')
    g = Node('G')
    h = Node('H')
    i = Node('I')
    j = Node('J')
    k = Node('K')
    l = Node('L')
    m = Node('M')
    n = Node('N')

    graph.set_nodes([a, b, c, d, e, f, g, h, i, j, k, l, m, n])
    graph.set_edges([
                Edge(a, d, 4, 1), Edge(a, c, 2, 0), Edge(b, d, 1, 0),
                Edge(b, e, 2, 2), Edge(c, d, 5, 0), Edge(c, f, 6, 0),
                Edge(c, g, 4, 1), Edge(d, e, 3, 0), Edge(d, k, 5, 2),
                Edge(e, h, 2, 2), Edge(e, i, 3, 0), Edge(f, g, 3, 2),
                Edge(f, l, 4, 1), Edge(g, j, 4, 2), Edge(g, k, 3, 0),
                Edge(h, i, 2, 0), Edge(h, n, 4, 0), Edge(j, m, 2, 0),
                Edge(k, m, 1, 0), Edge(k, n, 3, 0), Edge(l, m, 3, 0)])

    path = graph.find_path(l, e)
    print(path)
