#!/usr/bin/python
import shortest_path
import astar
import heapq
import math


class Node(shortest_path.Node):
    def __init__(self, name: str, x: int, y: int) -> None:
        self.name = name
        self.x = x
        self.y = y

        self.parent = None

        self.g = None
        self.h = None

    def f(self) -> int:
        self.g + self.h


class Graph(shortest_path.Graph):
    """
    def find_path(self, start: Node, end: Node) -> list:
        start.g = 0
        openlist = []
        closedlist = []
        heapq.heappush(openlist, (0, start))
        while len(openlist) != 0:
            q = heapq.heappop(openlist)[1]
            successors = self.get_neighbors(q)
            for successor, edge in successors.items():
                if successor is end:
                    return closedlist
                successor_current_cost = q.g + edge.distance
                if successor in openlist:
                    if successor.g <= successor_current_cost:
                        break
                elif successor in closedlist:
                    if successor.g <= successor_current_cost:
                        break
                    closedlist.remove(successor)
                    heapq.heappush(openlist, (successor.g + successor.h, successor))

                else:
                    heapq.heappush(openlist, (successor.g + successor.h, successor))
                    successor.h = math.sqrt((successor.x - end.x) ** 2 + (successor.y - end.y))
            print(closedlist)
            print(openlist)
            closedlist.append(q)
    """

    @staticmethod
    def heuristic_distance(node1: Node, node2: Node) -> int:
        return math.sqrt((node1.x - node2.x) ** 2 + (node1.y - node2.y) ** 2)

    @staticmethod
    def pop(l: list) -> Node:
        if len(l) == 0:
            return None
        lowest = l[0]
        for item in l:
            if item.g + item.h < lowest.g + lowest.h:
                lowest = item
        l.remove(lowest)
        return lowest

    def find_path(self, start: Node, end: Node) -> list:
        openlist = []
        closedlist = []
        start.h = Graph.heuristic_distance(start, end)
        start.g = 0
        openlist.append(start)
        while len(openlist) != 0:
            current = Graph.pop(openlist)
            if current is end:
                closedlist.append(current)
                break
            successors = self.get_neighbors(current)
            for successor, edge in successors.items():
                successor_current_cost = current.g + edge.distance
                if successor in openlist:
                    if successor.g <= successor_current_cost:
                        continue
                elif successor in closedlist:
                    if successor.g <= successor_current_cost:
                        continue
                    closedlist.remove(successor)
                    openlist.append(successor)
                else:
                    openlist.append(successor)
                    successor.h = Graph.heuristic_distance(successor, end)
                successor.g = successor_current_cost
                successor.parent = current
            closedlist.append(current)
            #print(current)
        path = []
        current = closedlist[len(closedlist)-1]
        while current is not start:
            path.append(current)
            current = current.parent
        return path


def create_nodes(x: int, y: int) -> list:
    nodes = []
    for row in range(x):
        for column in range(y):
            nodes.append(Node(str(row) + '_' + str(column), row, column))
    return nodes


def create_edges_except(nodes: list, exceptions: list) -> list:
    edges = []
    for node in nodes:
        if node.name not in exceptions:
            for element in nodes:
                if element.name not in exceptions and element is not node:
                    x_dif = node.x - element.x
                    y_dif = node.y - element.y
                    if x_dif in range(-1, 1) and y_dif in range(-1, 1):
                        edges.append(shortest_path.Edge(node, element, 1, 0))
    return edges


if __name__ == '__main__':
    graph = Graph()
    nodes = create_nodes(8, 8)
    graph.set_nodes(nodes)
    edges = create_edges_except(nodes, ['0_3', '0_4', '0_7', '1_7', '2_0', '2_3', '2_4', '2_7', '3_0', '3_3', '3_4',
                                        '3_7', '4_7', '5_2', '6_1', '6_2', '6_3', '6_4'])
    graph.set_edges(edges)
    path = graph.find_path(nodes[0], nodes[59])
    print(path)
    output = """"
  01234567
0 A..##..#
1 .......#
2 #..##..#
3 #..##..#
4 .......#
5 .#......
6 .####...
7 ...B....
        """

#      01234567
#    0 A..##..#
#    1 .......#
#    2 #..##..#
#    3 #..##..#
#    4 .......#
#    5 .#......
#    6 .####...
#    7 ...B....
