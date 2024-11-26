import numpy as np

from reference import Gene


class Network:
    """
    Represents a network of interactions, including transcription factors (TFs), miRNAs, and other genes.
    """
    def __init__(self):
        # dictionary of nodes, where each key is a node id and the value is a set of target nodes
        self.nodes = dict()                 # type: dict[str, set[str]]
        # set of transcription factor and miRNA node ids
        self.TFs = set()                    # type: set[str]
        self.miRNAs = set()                 # type: set[str]
        # dictionary of interactions, where each key is an interaction type and the value is a set of edges
        self.edge_dict = dict()             # type: dict[str, set[tuple[str, str]]]
        # set of all edges in the network, represented as tuples of nodes
        self.edges = set()                  # type: set[tuple[str, str]]
        # count of nodes regulated by either a TF or a miRNA
        self.all = 0                        # type: int

    def __eq__(self, other):
        if self.nodes.keys() != other.nodes.keys():
            return False
        for i, edges in self.edge_dict.items():
            if edges != other.edge_dict[i]:
                return False
        return True

    def __hash__(self):
        return hash(tuple(sorted(self.edges)))

    def get_interactions(self) -> set[str]:
        """ Retrieves all interaction types that have at least one edge. """
        return {interaction for interaction, edges in self.edge_dict.items() if edges}

    def count(self) -> tuple[tuple[str, int], ...]:
        """ Counts the number of edges for each interaction type. """
        return tuple([(interaction, len(edges)) for interaction, edges in sorted(self.edge_dict.items())])

    def suitable(self, interactions: list[str], n_random: int) -> bool:
        """
        Determines if the network has sufficient interaction combinations for a given set of interactions.
        :param interactions: list of interaction types to check
        :param n_random: the minimum number of random combinations required
        :return: true if there are enough combinations, false otherwise
        """
        existing_interactions = {interaction for interaction, edges in self.edge_dict.items() if len(edges) > 1}

        if not {'mirna-mirna', 'mirna-tf', 'mirna-gene'}.intersection(existing_interactions):
            print('no mirna')
            return False

        if not {'tf-mirna', 'tf-tf', 'tf-gene'}.intersection(existing_interactions):
            print('no tf')
            return False

        prod = np.prod([len(edges) for inter, edges in self.edge_dict.items() if inter in interactions and edges])

        if prod < n_random:
            print('not enough combinations')
            return False

        print('success')
        return True

    def edges_to_tsv(self, path_dict: dict[str, str]):
        """
        Writes the network's edges to TSV files for each interaction type.
        :param path_dict: a dictionary mapping interaction types to file paths
        """
        for interaction, file_path in path_dict.items():
            if interaction not in self.edge_dict.keys() or not self.edge_dict[interaction]:
                continue
            with open(file_path, 'w') as file:
                for r, t in sorted(self.edge_dict[interaction]):
                    file.write('{0}\t{1}\n'.format(r, t))

    def nodes_to_tsv(self, file_path: str):
        """
        Writes the network's nodes and their types (TF, miRNA, or other) to a TSV file.
        :param file_path: path to the output TSV file
        """
        with open(file_path, 'w') as file:
            for node in sorted(self.nodes.keys()):
                if node in self.TFs:
                    file.write('{0}\t{1}\n'.format(node, 'TF'))
                elif node in self.miRNAs:
                    file.write('{0}\t{1}\n'.format(node, 'miRNA'))
                else:
                    file.write('{0}\t{1}\n'.format(node, 'other'))

    def overview(self, file_path: str):
        """
        Writes an overview of the network to a file, including statistics about nodes, interactions, and regulators.
        :param file_path: path to the output file
        """
        actual_nodes = {node for edge in self.edges for node in edge}                       # type: set[str]

        with open(file_path, 'w') as file:
            file.write('total number of nodes by adjacency: {0:,}\n'.format(len(self.nodes)))
            file.write('number of nodes involved in edges: {0:,}\n'.format(len(actual_nodes)))
            file.write('regulators by adjacency: {0:,}\n'.format(len({n for n, targets in self.nodes.items() if targets})))
            file.write('regulators by edge list: {0:,}\n'.format(len({r for r, _ in self.edges})))
            file.write('targets by adjacency: {0:,}\n'.format(len({t for targets in self.nodes.values() for t in targets})))
            file.write('targets by edge list: {0:,}\n'.format(len({t for _, t in self.edges})))
            file.write('any regulated by TF or miRNA: {0:,}\n'.format(self.all))
            file.write('total number of TFs: {0:,}\n'.format(len(self.TFs)))
            file.write('actual TF: {0:,}\n'.format(len(actual_nodes.intersection(self.TFs))))
            file.write('total number of miRNAs: {0:,}\n'.format(len(self.miRNAs)))
            file.write('actual miRNAs: {0:,}\n'.format(len(actual_nodes.intersection(self.miRNAs))))

            for interaction, edges in sorted(self.edge_dict.items()):
                file.write('{0}: {1:,}\n'.format(interaction, len(edges)))
            file.write('total edges: {0:,}\n'.format(len(self.edges)))


class OriginalNetwork(Network):
    """ Represents the original network of interactions. """
    def __init__(self, nodes: dict[str, Gene], edge_dict: dict[str, set[tuple[str, str]]]):
        super().__init__()

        # classify nodes as TFs or miRNAs
        for gene in nodes.values():
            if gene.is_miRNA:
                self.miRNAs.add(gene.id)
            if gene.is_TF:
                self.TFs.add(gene.id)

            self.nodes[gene.id] = set()

        # filter and initialize edges
        self.edge_dict = {interaction: {(r, t) for r, t in edges if r in self.nodes and t in self.nodes}
                          for interaction, edges in edge_dict.items() if edges}
        self.edges = tuple(sorted({edge for edges in self.edge_dict.values() for edge in edges
                                   if edge[0] != edge[1]}))

        # add targets for each regulator
        for regulator, target in self.edges:
            self.nodes[regulator].add(target)

        # calculate the number of nodes regulated by TFs or miRNAs
        self.all = len({target for regulator in self.miRNAs.union(self.TFs) for target in self.nodes[regulator]})


class ReadNetwork(Network):
    """ Represents a network of interactions read from files (can be used to read in randomized networks). """
    def __init__(self, nodes: set[str], edge_paths: dict[str, str]):
        super().__init__()

        # initialize nodes and classify as TF or miRNA
        for node, kind in nodes:
            self.nodes[node] = set()

            if kind == 'TF':
                self.TFs.add(node)
            elif kind == 'miRNA':
                self.miRNAs.add(node)

        # read the edges from the edge type files
        for interaction, edge_path in edge_paths.items():
            self.edge_dict[interaction] = set()
            with open(edge_path, 'r') as file:
                for line in file:
                    if line.startswith('#'):
                        continue

                    cols = tuple(line.strip().split())     # type: tuple[str, str]

                    self.nodes[cols[0]].add(cols[1])
                    self.edge_dict[interaction].add(cols)

        # create the set of all edges (removes self-loops)
        self.edges = {(regulator, target) for regulator, targets in self.nodes.items() for target in targets}
        # calculate the number of nodes regulated by TFs or miRNAs
        self.all = len({r for n in self.miRNAs.union(self.TFs) for r in self.nodes[n]})
