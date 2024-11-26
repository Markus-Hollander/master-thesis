from collections import defaultdict
from csv import dictReader
from functools import partial
import numpy as np
import os
from shutil import copy2
import subprocess as sp
from time import time
from typing import Any

import computation
from gene_functionality import GeneExpression, ModCounter, GeneFunctionIndex
from manager import Manager
from network import OriginalNetwork, ReadNetwork
from reference import Gene, Reference
from transition import Differential
import util


def i(line: str) -> int:
    """ Converts a scientific notation string or regular string to an integer. """
    if 'e+' in line:
        a, b = line.split('e+')
        return int(float(a) * pow(10, int(b)))
    return int(line)


def reference_to_bed(gene_dict: dict[str, Gene], gene_ids: set[str], bed_path):
    """
    Writes gene and promoter intervals to a BED file.
    :param gene_dict: dictionary of genes with gene IDs as keys
    :param gene_ids: set of gene IDs to include in the BED file
    :param bed_path: bath to the result bed file
    """
    intervals = []
    with open(bed_path, 'w') as bed:
        for g in [gene for gene in gene_dict.values() if gene.id in gene_ids]:
            intervals.append((g.chromosome, g.start, g.end + 1, g.id + '/gene'))
            intervals.append((g.chromosome, g.promoter_start,  g.promoter_end + 1, g.id + '/promoter'))

        for interval in sorted(intervals, key=lambda x: (x[0], x[1], x[2])):
            bed.write('\t'.join([str(x) for x in interval]) + '\n')


def result_to_bed(p: float, results: list[tuple[str, str]], bed_path: str) -> dict[tuple[str, int, int], str]:
    """
    Processes results and writes relevant regions to a BED file.
    :param p: threshold
    :param results: a list of tuples containing cell types and file paths to result files
    :param bed_path: bath to the result bed file
    :returns: a dictionary mapping genomic regions to corresponding cell types
    """
    bin_count = dict()

    with open(bed_path, 'w') as bed_file:
        for cell, result_path in results:
            with open(result_path, 'r') as in_file:
                for line in in_file:
                    cols = line.strip().split('\t')

                    post = float(cols[8].split(';')[0].split('=')[1])

                    if post <= p:
                        continue

                    chromosome, start, end = util.process_chromosomes(cols[0]), i(cols[3]), i(cols[4]) + 1

                    bed_file.write('\t'.join([str(x) for x in [chromosome, start, end, cell]]) + '\n')
                    bin_count[(chromosome, start, end)] = cell

    return bin_count


def load_histone(file_path: str, a: str, b: str) -> tuple[dict[str, dict[str, ModCounter]], set[tuple[str, int, int]]]:
    """
    Loads histone modification data and tracks overlap regions.
    :param file_path: path to the file containing histone data
    :param a: modification type A to track
    :param b: modification type B to track
    :return: a tuple containing:
        - (dict[str, dict[str, ModCounter]]): histone modifications by region and gene.
        - (set[tuple[str, int, int]]): set of overlap regions
    """
    mods = {key: defaultdict(ModCounter)
            for key in ['combined', 'gene', 'promoter']}                        # type: dict[str, dict[str, ModCounter]]

    overlap = set()

    with open(file_path, 'r') as in_file:
        for line in in_file:
            cols = line.strip().split('\t')
            gene_id, region = cols[3].split('/')
            mod = cols[7]

            # update modification counts for region and combined categories
            for key in [region, 'combined']:
                if mod == a:
                    mods[key][gene_id].A += 1
                elif mod == b:
                    mods[key][gene_id].B += 1

                overlap.add((cols[4], int(cols[5]), int(cols[6])))

    return mods, overlap


def load_expression(gene_type_path: str, exp_threshold: float, gene_ids: set[str], file_path: str) -> dict[str, GeneExpression]:
    """
    Loads gene expression data and filters by given threshold.
    :param gene_type_path: path to the file containing gene types
    :param exp_threshold: the expression threshold for filtering results
    :param gene_ids: set of gene IDs to load
    :param file_path: path to the file containing expression data
    :return: dictionary mapping gene IDs to their corresponding gene expression data
    """
    # load gene types into dictionary
    gene_types = dict()
    with open(gene_type_path, 'r') as file:
        for line in file:
            gene_id, gene_type = line.strip().split()

            gene_types[gene_id] = gene_type

    # load expression data
    genes = dict()
    with open(file_path, 'r') as in_file:
        for line in dictReader(in_file, delimiter='\t'):
            if line['padj'] in ['NA']:
                continue

            gene_id = line['name']
            if gene_id not in gene_ids:
                continue

            # filter expression data and store gene expressions
            genes[gene_id] = GeneExpression(gene_types[gene_id], line['log2FoldChange'], line['pvalue'], line['padj'],
                                            exp_threshold)

    return genes


def load_results(cfg: dict, manager: Manager, reference: Reference, transitions: dict[tuple[str, str], Differential]):
    """
    Loads results from the analysis, including gene expression, histone data, and gene definitions.
    :param cfg: configuration dictionary containing analysis parameters
    :param manager: Manager instance to handle database operations
    :param reference: Reference instance for gene and transcript data
    :param transitions: dictionary containing transition data between cell types
    """
    t_list = list(transitions.values())

    # write gene types to file
    def write_gene_types():
        reference.write_gene_types(manager.gene_type_path)

    util.execute(write_gene_types, None, 'write gene types', 2)

    # load expression results for each transition
    func = partial(load_expression, manager.gene_type_path, float(cfg['analysis']['expression fdr']),
                   set(reference.genes.keys()))
    tasks = [t.dea_result_path for t in t_list]
    expression_results = util.multiprocess_tasks(func, tasks, 'loading {0} expression results', indent=2)

    # map gene function indices for each transition
    indices = {(t.cell_1, t.cell_2): GeneFunctionIndex(t.cell_1, t.cell_2, t.closest_prefix, t.bed_path, result)
               for t, result in zip(t_list, expression_results)}

    func = partial(reference_to_bed, reference.genes)
    tasks = [(set(indices[k].exp.keys()), t.bed_path) for k, t in transitions.items()
             if not os.path.isfile(t.bed_path)]
    util.multiprocess_tasks(func, tasks, 'generating {0} gene definition .bed files', 2)

    # process histone modification results
    histones = [histone for transition in t_list for histone in transition.histone.values()]
    histone_tasks = [(t.bed_path, h.result_bed, h.intersect_path, t.cell_1, t.cell_2, h.histone)
                     for t in t_list for h in t.histone.values()]

    func = partial(result_to_bed, float(cfg['analysis']['posterior threshold']))
    tasks = [([(h.cell_1, h.cell_1_results), (h.cell_2, h.cell_2_results)], h.result_bed) for h in histones
             if not os.path.isfile(h.result_path)]
    util.multiprocess_tasks(func, tasks, 'generating {0} result .bed files', 2)

    # intersect histone data
    util.multiprocess_tasks(util.intersect, [x[:3] for x in histone_tasks if not os.path.isfile(x[2])],
                            'generating {0} intersect .bed files', 2)

    # load histone modification data
    results = util.multiprocess_tasks(load_histone, [x[2:-1] for x in histone_tasks], 'loading {0} histone results', 2)

    def add_modifications(args):
        """ Updates gene function indices with histone modification data. """
        for task, result in args:
            _, _, _, cell_1, cell_2, histone = task
            mods, overlap = result

            indices[(cell_1, cell_2)].add_modification(histone, mods)

    # add histone modifications to gene function indices
    util.execute(add_modifications, list(zip(histone_tasks, results)), 'adding {0} histone results', 2)

    # build the databases for each gene function index
    start = time()
    print('\t\tbuilding the databases', end='\r')
    for index in indices.values():
        manager.to_db(index)
    print('\t\tbuilding the databases: done')
    util.print_time(start, 2)


def build_networks(manager: Manager, reference: Reference, specifics: set[tuple[str, tuple[str, str, bool]]]):
    """
    Builds the original and transition-specific networks based on the provided configurations.
    :param manager: manager object handling the network building process and file paths
    :param reference: reference object containing gene and edge data
    :param specifics: set of tuples defining specific transitions for network construction
    :return: tuple consisting of
        - list of transition-specific networks that were successfully built
        - dictionary of transition networks keyed by edge type and specific transition
        - set of transitions that could not be built
    """
    # 1.  BUILD ORIGINAL NETWORKS
    # 1a. complete networks
    def generate_complete_networks():
        """ Generates the complete networks for each edge type if they do not already exist. """
        for edge_type in manager.edge_types:
            node_path = manager.node_path(edge_type)
            edge_paths = {interaction: manager.edge_path(edge_type, interaction, None, -1)
                          for interaction in manager.interactions}

            # skip if all network files already exist
            if os.path.isfile(node_path) and all(os.path.isfile(edge_path) for edge_path in edge_paths.values()):
                continue

            # create and overview the original network
            network = OriginalNetwork(reference.genes, reference.edges[edge_type])
            network.overview(os.path.join(manager.network_stats_dir, 'original_edges={0}.txt'.format(edge_type)))

            # skip unsuitable networks
            if not network.suitable(manager.interactions, manager.n_random):
                continue

            # write edges and nodes to TSV files
            network.edges_to_tsv({interaction: manager.edge_path(edge_type, interaction, None, -1)
                                  for interaction in manager.interactions})
            network.nodes_to_tsv(manager.node_path(edge_type))
            manager.n_edges[(edge_type, None)] = len(network.edges)

    util.execute(generate_complete_networks, None, 'generating complete networks', 2)

    # 1b. transition specific networks
    tasks = list(specifics)

    to_remove = set()
    original_networks = dict()  # type: dict[tuple[str, tuple[str, str, bool]], OriginalNetwork]

    def generate_transition_networks():
        """ Generates transition-specific networks based on given edge types and transitions. """
        for edge_type, specific in tasks:
            gene_ids = manager.gene_ids(*specific)
            nodes = {gene_id: reference.genes[gene_id] for gene_id in gene_ids}
            network = OriginalNetwork(nodes, reference.edges[edge_type])

            network.overview(manager.network_stats_fmt.format(*specific, edge_type))

            # skip unsuitable networks
            if not network.suitable(manager.interactions, manager.n_random):
                to_remove.add((edge_type, specific))
                continue

            # write transition networks to TSV files
            network.edges_to_tsv({interaction: manager.edge_path(edge_type, interaction, specific, -1)
                                  for interaction in manager.interactions})
            network.nodes_to_tsv(manager.node_path(edge_type, specific))
            original_networks[(edge_type, specific)] = network

            manager.n_edges[(edge_type, specific)] = len(network.edges)

        return to_remove, original_networks

    to_remove, original_networks = util.execute(generate_transition_networks, None,
                                                'generating {0:,} transition networks'.format(len(tasks)), 2,
                                                removed=0)
    return list(set(specifics).difference(to_remove)), original_networks, to_remove


def degree_distributions(network: OriginalNetwork) -> list[tuple[str, int, int, str]]:
    """
    Computes the degree distribution for nodes in the network, categorized by type (TF, miRNA, genes).
    :param network: network object whose node degree distribution is to be computed
    :return: list of tuples (gene type, node degree, gene count, semicolon separated list of gene ids)
    """
    dists = {'TF': defaultdict(set),
             'miRNA': defaultdict(set),
             'genes': defaultdict(set)}

    genes = defaultdict(set)

    for node, regulated in network.nodes.items():
        if node in network.TFs:
            dists['TF'][len(regulated)].add(node)
        elif node in network.miRNAs:
            dists['miRNA'][len(regulated)].add(node)
        else:
            for r in regulated:
                genes[r].add(node)

    # aggregate the distribution for genes
    for node, regulators in genes.items():
        dists['genes'][len(regulators)].add(node)

    # format and return the degree distributions
    return [(key, b, len(gene_ids), ';'.join(sorted(gene_ids)))
            for key, dist in dists.items() for b, gene_ids in dist.items()]


def degree_correlation(manager: Manager, specific: tuple[str, str, bool], histone: str,
                       region: str,  gene_ids: str) -> tuple[int, float, str, float, float, float, str]:
    """
    Computes the degree correlation for a given set of gene IDs, histone modification, and region.
    :param manager: manager object handling the pair retrieval
    :param specific: tuple containing two cell types and a differential flag
    :param histone: histone modification name
    :param region: genomic region to analyze
    :param gene_ids: string of semicolon separated gene IDs
    :return:
    """
    pairs = manager.pairs(specific, histone, region, 'mod_diff', set(gene_ids.split(';')))

    return (len(pairs),) + computation.correlation(manager.correlation_min, pairs) + \
        computation.chi_square(manager.correlation_min, pairs)


def process_degree_distributions(manager: Manager, specifics: list[tuple[str, tuple[str, str, bool]]],
                                 networks: dict[tuple[str, tuple[str, str, bool]], OriginalNetwork],
                                 histones: dict[tuple[str, str], set[str]]):
    """
    Processes degree distributions for networks and computes degree correlations for different combinations.
    :param manager: manager object handling various network-related task
    :param specifics: list of transition specifics to process
    :param networks: networks to analyze
    :param histones: histone modification data to use in the analysis
    """
    # compute node degree distributions for each network
    tasks = [networks[specific] for specific in specifics]
    distributions = util.multiprocess_tasks(degree_distributions, tasks, 'computing node degrees for {0:,} networks', 2)
    distributions = {specific: dist for specific, dist in zip(specifics, distributions)}

    # compute degree correlations for each combination
    combinations = [specific + (histone, region, *d)
                    for specific, dist in distributions.items()
                    for histone in histones[(specific[1][0], specific[1][1])]
                    for region in manager.regions
                    for d in dist]

    f = partial(degree_correlation, manager)
    tasks = [(*comb[1:4], comb[-1]) for comb in combinations]
    correlations = util.multiprocess_tasks(f, tasks, 'computing degree correlation for {0:,} combinations', 2)

    # apply correction to p-values
    def correction():
        """ Corrects the p-values for multiple hypothesis testing using Benjamini-Hochberg procedure. """
        p_val_dict = {task: corr[3] for task, corr in zip(tasks, correlations)}

        return computation.benjamini_dict(p_val_dict, manager.correlation_q, krieger=True)

    p_values = util.execute(correction, None, 'correcting multiple hypothesis testing', 2)

    # add the degree correlation results to the manager
    def adding_degree_correlation():
        """ Adds the computed degree correlation results to the manager's database. """
        content = defaultdict(lambda: defaultdict(list))

        for comb, correlation in zip(combinations, correlations):
            edge_type, specific, histone, region, category, degree, n_genes, gene_ids = comb
            content[edge_type][(specific, histone, region)].append((category, degree, n_genes, gene_ids)
                                                                   + correlation
                                                                   + (p_values[(specific, histone, region, gene_ids)],))
        manager.add_correlation(content, degree=True)

    util.execute(adding_degree_correlation, None, 'adding degree correlation', 2)


def randomise(script_path: str, edge_path: str, rewired_path: str, n: int):
    """
    Randomizes the edges in the network by performing edge rewiring via an R script.
    :param script_path: path to the R script used for edge rewiring
    :param edge_path: path to the original edge file
    :param rewired_path: path where the rewired network will be saved
    :param n: number of rewiring iterations
    """
    try:
        cmd_args = ['Rscript', script_path, edge_path, rewired_path, str(n)]
        sp.check_call(args=' '.join(cmd_args), stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
    except sp.CalledProcessError:
        # if R script fails, copy the original edge file to the rewired path
        copy2(edge_path, rewired_path)


def compute_motifs(manager: Manager, edge_type: str,
                   specific: tuple[str, str, bool], n: int) -> list[tuple[str, str, str, str]]:
    """
    Computes motifs for TF-miRNA co-regulation based on the network data.
    :param manager: manager handling network files and interactions
    :param edge_type: type of the edge (e.g. experimental)
    :param specific: tuple containing two cell types and a differential flag
    :param n: network identifier (negative for original network, non-negative for randomized)
    :return: list of (TF name, miRNA name, motif type, semicolon-separated list of shared genes) tuples
    """
    # load nodes from the specified node file
    with open(manager.node_path(edge_type, specific), 'r') as file:
        nodes = {tuple(line.strip().split('\t')) for line in file}

    # collect edge file paths for the given interactions
    edge_paths = {interaction: manager.edge_path(edge_type, interaction, specific, n)
                  for interaction in manager.interactions}
    edge_paths = {interaction: file_path for interaction, file_path in edge_paths.items()
                  if os.path.exists(file_path)}

    # initialize the network with nodes and edges
    network = ReadNetwork(nodes, edge_paths)

    # calculate the total number of shared regulatory targets
    all_shared = len({r for n in network.miRNAs.union(network.TFs) for r in network.nodes[n]})

    # stores shared targets for tf-miRNA pairs
    shared_dict = dict()                                                    # type: dict[tuple[str, str], set[str]]

    # iterate over all pairs of TFs and miRNAs in the network
    for tf in network.TFs:
        for mirna in network.miRNAs:
            # find shared targets between the current tf and miRNA
            shared_targets = network.nodes[tf].intersection(network.nodes[mirna])

            # skip if no shared targets are found
            if not shared_targets:
                continue

            # store the shared targets in the dictionary
            shared_dict[(tf, mirna)] = shared_targets

    # compute p-values for each tf-miRNA pair with shared targets
    p_values = dict()

    for k, shared_targets in shared_dict.items():
        tf, mirna = k
        p_values[(tf, mirna)] = computation.hyper_geometric(population=all_shared,
                                                            population_successes=len(network.nodes[mirna]),
                                                            draws=len(network.nodes[tf]),
                                                            draw_successes=len(shared_targets))

    # adjust p-values for multiple testing using the benjamini-hochberg procedure
    adjusted_p_values = computation.benjamini_dict(p_values, manager.motif_q, krieger=True)

    # classify motifs based on the adjusted p-values and the regulatory relationships
    table = []
    for k, sig in adjusted_p_values.items():
        if sig:
            tf, mirna = k
            if tf not in network.nodes[mirna] and mirna not in network.nodes[tf]:
                table.append((tf, mirna, 'co-regulatory', ';'.join(shared_dict[k])))
            elif tf not in network.nodes[mirna] and mirna in network.nodes[tf]:
                table.append((tf, mirna, 'TF', ';'.join(shared_dict[k])))
            elif tf in network.nodes[mirna] and mirna not in network.nodes[tf]:
                table.append((tf, mirna, 'miRNA', ';'.join(shared_dict[k])))
            elif tf in network.nodes[mirna] and mirna in network.nodes[tf]:
                table.append((tf, mirna, 'composite', ';'.join(shared_dict[k])))

    return table


def compute_correlation(manager: Manager, edge_type: str, specific: tuple[str, str, bool], histone: str, region: str,
                        motif: str,
                        network_number: int) -> tuple[int, int, int, int, int, str, str, str, str, str, float, str]:
    """
    Computes correlations for motif interactions based on specific criteria.
    :param manager: manager object handling data retrieval and storage
    :param edge_type: type of the edge (e.g. experimental)
    :param specific: tuple containing two cell types and a differential flag
    :param histone: name of a histone modification of interest
    :param region: genomic region of interest
    :param motif: the motif type (e.g., composite, co-regulatory)
    :param network_number: network identifier (negative for original network, non-negative for randomized)
    :return: tuple consisting of counts, sorted lists of gene IDs, correlation statistics and summaries
    """
    # fetch motifs matching the criteria
    motifs = manager.motif(edge_type, specific, network_number, motif)

    # extract unique transcription factors, miRNAs, and shared targets
    tfs = {tf for tf, _, _, _ in motifs}
    miRNAs = {mirna for _, mirna, _, _ in motifs}
    targets = {target for _, _, shared_str, _ in motifs for target in shared_str.split(';')}

    # collect all unique gene identifiers, including tf and miRNA
    gene_ids = {gene_id
                for tf, mirna, shared_str, _ in motifs
                for gene_id in shared_str.split(';') + [tf, mirna]}

    # create a formatted string summarizing the motifs
    m_str = ' | '.join(['[{0}+{1}+{2}]'.format(tf, mirna, shared_str) for tf, mirna, shared_str, _ in motifs])

    # compute lengths of different groups and prepare sorted strings
    lengths = tuple(len(x) for x in [gene_ids, motifs, tfs, miRNAs, targets])   # type: tuple[int, int, int, int, int]
    strings = (';'.join(sorted(gene_ids)),
               m_str,
               ';'.join(sorted(tfs)),
               ';'.join(sorted(miRNAs)),
               ';'.join(sorted(targets)))                                       # type: tuple[str, str, str, str, str]

    # retrieve gene pairs for the correlation calculation
    pairs = manager.pairs(specific=specific, histone=histone, region=region, mod='mod_diff', gene_ids=gene_ids)

    # compute correlation for composite motifs if applicable
    if motif == 'composite':
        return lengths + strings + computation.correlation(manager.correlation_min_composite, pairs)

    # compute correlation for other motifs
    return lengths + strings + computation.correlation(manager.correlation_min, pairs)


def p_help(original: float, random: list[float], n: int) -> tuple[float, float, float, float, float]:
    """
    Calculates various p-values based on observed and random data.
    :param original: p-value of the original network
    :param random: list of p-values belonging to the randomized networks
    :param n: number of random networks generated
    :return: a tuple containing different p-value measures (smaller, greater, absolute comparisons)
    """
    def p_val(values: list[float]) -> float:
        """ Calculate p-value as the fraction of relevant random values adjusted by the number of samples. """
        return max(min((len(values) + 1) / n, 1.0), 0.0)

    # return NaN values if the original value is invalid
    if original is None or np.isnan(original):
        return np.nan, np.nan, np.nan, np.nan, np.nan

    # compute p-values for negative original values
    if original < 0:
        smaller = p_val([r for r in random if original <= r <= 0])
        greater = p_val([r for r in random if r <= original])
        other = p_val([r for r in random if r > 0])
    # compute p-values for positive original values
    else:
        smaller = p_val([r for r in random if 0.0 <= r <= original])
        greater = p_val([r for r in random if r >= original])
        other = p_val([r for r in random if r < 0])

    # compute absolute comparisons for the original value
    original = abs(original)
    absolute_smaller = p_val([r for r in random if abs(r) <= original])
    absolute_greater = p_val([r for r in random if abs(r) >= original])

    return smaller, greater, 1.0 - other, absolute_smaller, absolute_greater


def compute_p_values(manager: Manager, edge_type: str, specific: tuple[str, str, bool], histone: str,
                     region: str, motif: str) -> tuple[int, int, int, int, int, str, str, str, str, str, float, str,
                                                       str, float, float, float, float, float,
                                                       float, float, float, str]:
    """
    Computes p-values for motif interactions and associated data.
    :param manager: manager instance handling data retrieval and storage
    :param edge_type: type of the edge (e.g. experimental)
    :param specific: tuple containing two cell types and a differential flag
    :param histone: name of a histone modification of interest
    :param region: genomic region of interest
    :param motif: the motif type (e.g., composite, co-regulatory)
    :return: tuple containing gene counts, p-values, correlation metrics and chi-square statistics
    """
    # handle specific cases where motifs correspond to predefined categories
    if motif in ['all genes', 'all TF', 'all miRNA', 'other']:
        # get gene IDs based on the motif category
        if motif == 'all genes':
            gene_ids = manager.gene_ids(*specific)
        else:
            gene_ids = manager.gene_ids(*specific, gene_type=motif.replace('all ', ''))

        # compute gene pairs and chi-square statistics
        pairs = manager.pairs(specific, histone, region, 'mod_diff', gene_ids)
        chi = computation.chi_square(manager.correlation_min, pairs)

        # prepare a string representation of the pairs
        pairs_str = ';'.join([str(pair) for pair in pairs])                 # type: str

        # compute correlation statistics
        r, c = computation.correlation(manager.correlation_min, pairs)

        # return results for predefined categories
        return (len(gene_ids), 0, 0, 0, 0, ';'.join(sorted(gene_ids))) + ('n.a.',) * 4 + (r, c, pairs_str) \
               + (np.nan,) * 5 + chi

    # handle general motifs by fetching correlation data
    corr = manager.correlation(edge_type, specific, histone, region, motif)
    o_r, r_r, c = corr[:3]          #
    stats = corr[3:]                # type: tuple[int, int, int, int, int, str, str, str, str, str]

    # extract gene IDs and compute pairs
    gene_ids = set(corr[8].split(';'))
    pairs = manager.pairs(specific, histone, region, 'mod_diff', gene_ids)

    # prepare a string representation of the pairs
    pairs_str = ';'.join([str(pair) for pair in pairs])  # type: str

    # compute chi-square statistics based on motif type
    if motif == 'composite':
        chi = computation.chi_square(manager.correlation_min_composite, pairs)
    else:
        chi = computation.chi_square(manager.correlation_min, pairs)

    # return NaN p-values if random correlations are unavailable
    if not r_r:
        return stats + (o_r, c, pairs_str) + (np.nan,) * 5 + chi

    # compute p-values for the observed correlation
    return stats + (o_r, c, pairs_str) + p_help(o_r, r_r, manager.n_random) + chi


def process_randomised_networks(manager: Manager, script_path: str,
                                specifics: list[tuple[str, Any | tuple[str, str, bool]]],
                                removed: list[tuple[str, Any | tuple[str, str, bool]]],
                                histones: dict):
    """
    Processes randomised networks by generating, analyzing motifs, computing correlations, and calculating p-values.
    :param manager: manager instance handling data retrieval and storage
    :param script_path: path to the R script used for randomizing networks
    :param specifics: list of tuples containing two cell types and a differential flag
    :param removed: list of removed cell type transitions
    :param histones: mapping of cell type pairs to histone modifications
    :return:
    """
    # generate a list of all networks to process, including randomised ones
    networks = [(edge_type, interaction_type, specific, network_number)
                for interaction_type in manager.interactions
                for edge_type, specific in specifics + [(e, None) for e in manager.edge_types]
                for network_number in range(-1, manager.n_random)]

    # map each network to its corresponding file path
    edge_paths = {network: manager.edge_path(network[0], *network[1:])
                  for network in networks}  # type: dict[tuple[str, str, Any | tuple[str, str, bool], int], str]

    # prepare tasks for randomising networks that are not yet generated
    tasks = [(manager.edge_path(*network[:3], -1), file_path, manager.n_edges[(network[0], network[2])])
             for network, file_path in edge_paths.items() if network[-1] > -1]
    tasks = [(path_1, path_2, n) for path_1, path_2, n in tasks
             if os.path.exists(path_1) and not os.path.exists(path_2)]

    # parallelise the generation of randomised networks
    f = partial(randomise, script_path)
    util.multiprocess_tasks(f, tasks, 'generating {0:,} randomised networks', 2)

    # prepare tasks for computing motifs in all networks
    f = partial(compute_motifs, manager)
    tasks = list({(edge_type, specific, network_number) for edge_type, _, specific, network_number in networks
                  if manager.create_motif_table(edge_type, specific, network_number)})
    motifs = util.multiprocess_tasks(f, tasks, 'computing motifs for {0:,} networks', 2)

    # addition of motifs to the manager
    def add_motifs():
        """ Group motifs by edge type and specific cell type transition. """
        content = defaultdict(lambda: defaultdict(list))
        for key, corr in zip(tasks, motifs):
            content[key[0]][key[1]].append((key[2], corr))
        manager.add_motifs(content)
    util.execute(add_motifs, None, 'adding motifs', 2)

    # prepare tasks for computing correlations for all combinations
    f = partial(compute_correlation, manager)
    tasks = [(edge_type, specific, histone, region, motif, network_number)
             for edge_type, specific in specifics
             for histone in histones[(specific[0], specific[1])]
             for region in manager.regions
             for motif in manager.motif_names + ['all']
             for network_number in range(-1, manager.n_random)
             if not manager.correlation_table_exists(edge_type, specific, histone, region, degree=False)]
    correlation = util.multiprocess_tasks(f, tasks, 'computing correlation for {0:,} combinations', 2)

    # addition of correlations to the manager
    def add_correlation():
        """ Group correlations by edge type and specific criteria. """
        content = defaultdict(lambda: defaultdict(list))
        for key, corr in zip(tasks, correlation):
            content[key[0]][key[1:4]].append((*key[4:], *corr))
        manager.add_correlation(content)
    util.execute(add_correlation, None, 'adding correlation for {0:,} combinations'.format(len(tasks)), 2)

    # prepare tasks for computing p-values for all combinations
    f = partial(compute_p_values, manager)
    tasks = [(edge_type, specific, histone, region, motif)
             for edge_type, specific in specifics
             for histone in histones[(specific[0], specific[1])]
             for region in manager.regions
             for motif in manager.motif_names + ['all', 'all genes', 'all TF', 'all miRNA', 'other']]
    # add additional tasks for removed cell type transition
    tasks += [(edge_type, specific, histone, region, motif)
              for edge_type, specific in removed
              for histone in histones[(specific[0], specific[1])]
              for region in manager.regions
              for motif in ['all genes', 'all TF', 'all miRNA', 'other']]
    p_values = util.multiprocess_tasks(f, tasks, 'computing p-values for {0:,} combinations', 2)

    # p-value adjustment
    def adjust_p_values():
        # create a dictionary of p-values for multiple testing correction
        p_dict_list = [{task: p_row[i] for task, p_row in zip(tasks, p_values)} for i in range(13, 19)]
        adjusted = [computation.benjamini_dict(p_dict, manager.correlation_q, krieger=True) for p_dict in p_dict_list]
        # return adjusted p-values as a tuple for each task
        return [tuple(adj[task] for adj in adjusted) for task in tasks]

    adjusted = util.execute(adjust_p_values, None, 'performing multiple test correction', 2)

    # addition of p-values to the manager
    def add_p_values(args):
        """ Group p-values and adjusted values by edge type and specific criteria. """
        content = defaultdict(lambda: defaultdict(list))
        for arg in args:
            edge_type, specific, histone, region, motif = arg[0]
            p = arg[1]
            adj = arg[2]
            content[edge_type][specific[-1]].append((*specific[:-1], histone, region, motif) + p + adj)

        manager.add_p_values(content)
    util.execute(add_p_values, list(zip(tasks, p_values, adjusted)), 'adding {0:,} p-values', 2)


def process(cfg: dict, manager: Manager, reference: Reference, transitions: dict[tuple[str, str], Differential]):
    """
    Executes the full pipeline for processing networks, including differential expression, network generation, and motif
    analysis.
    :param cfg: configuration dictionary containing paths and settings
    :param manager: manager instance handling data and computations
    :param reference: reference object for annotation and data lookup
    :param transitions: transition data mapping cell pairs to differential data
    """
    # INIT
    script_path = cfg['file paths']['analysis']['rewire script']

    specifics = {(transition.cell_1, transition.cell_2, diff)
                 for transition in transitions.values()
                 for diff in manager.differential}

    histones = {(transition.cell_1, transition.cell_2): {histone for histone in transition.histone.keys()}
                for transition in transitions.values()}

    # LOAD DIFFERENTIAL EXPRESSION AND HISTONE ANALYSIS RESULTS
    util.execute(partial(load_results, cfg, manager, reference), transitions,
                 'pre-process differential expression and histone analysis results for {0} transitions', 1,
                 override=False)
    print()

    # BUILD ORIGINAL NETWORKS
    specifics, original_networks, removed = util.execute(partial(build_networks, manager, reference),
                                                         [(edge_type, specific)
                                                          for specific in specifics
                                                          for edge_type in manager.edge_types],
                                                         'building networks', 1, override=False)
    print()

    # DEGREE DISTRIBUTION
    util.execute(process_degree_distributions, (manager, specifics, original_networks, histones),
                 'process correlation depending on node degree', 1, override=False)
    print()

    # ANALYSE MOTIFS IN RANDOMISED NETWORKS
    util.execute(process_randomised_networks, (manager, script_path, specifics, removed, histones),
                 'processing motifs in randomised networks', 1, override=False)
