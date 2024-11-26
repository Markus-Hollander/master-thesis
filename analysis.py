from collections import defaultdict
from csv import dictReader
import numpy as np
import os
import pandas as pd
from statistics import stdev
from time import time
# matplotlib needs to be imported and configured first
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

import computation
from manager import Manager
from reference import Reference
from transition import Differential
import util


def nan_round(x: float, r: int):
    """
    Checks if an input is a number before rounding.
    """
    if computation.is_nan(x):
        return x
    return round(x, r)


def bubble_plot(manager: Manager):
    """
    Generates bubble plots for differential expression and histone modification data.
    :param manager: Manager object with data configuration
    """
    # initialize the overall plot style
    sns.reset_defaults()
    sns.set_style('whitegrid')
    # prepare the dictionary if it does not exist yet
    dir_path = os.path.join(manager.result_dir, 'bubble_plots')
    os.makedirs(dir_path, exist_ok=True)
    # define the plots to tbe generated
    tasks = [(('neural.progenitor.cell', 'bipolar.neuron', True), 'H3K27ac', 'promoter', 'mod_diff', 'TF'),
             (('neural.progenitor.cell', 'bipolar.neuron', False), 'H3K27ac', 'promoter', 'mod_diff', 'TF'),
             (('neural.progenitor.cell', 'bipolar.neuron', False), 'H3K4me3', 'gene', 'mod_diff', 'TF'),
             (('neural.progenitor.cell', 'bipolar.neuron', True), 'H3K4me3', 'gene', 'mod_diff', 'TF'),
             (('neural.progenitor.cell', 'bipolar.neuron', False), 'H3K4me3', 'gene', 'mod_diff', ''),
             (('neural.progenitor.cell', 'bipolar.neuron', True), 'H3K4me3', 'gene', 'mod_diff', '')]

    for task in tasks:
        gene_ids = manager.gene_ids(*task[0], task[-1])
        pairs = manager.pairs(*task[:-1], gene_ids)

        count = defaultdict(int)

        for pair in pairs:
            count[pair] += 1

        total = sum(count.values())

        scale = 10 / (total / 1000)

        exp, mod, c = zip(*[(*key, v) for key, v in count.items() if v != 0])
        c = [x * scale for x in c]
        plt.scatter(mod, exp, s=c, c='none', edgecolors='#4f81bd', linewidths=2)

        for k, v in count.items():
            y, x = k
            t = '{0:,}'.format(v)
            if (v * scale) < 500:
                if y in [-1, 0]:
                    plt.annotate(t, xy=(x, y + 0.2), color='#c0504d', weight='bold', ha='center', va='center')
                else:
                    plt.annotate(t, xy=(x, y - 0.2), color='#c0504d', weight='bold', ha='center', va='center')
            else:
                plt.annotate(t, xy=(x, y), color='#c0504d', weight='bold', ha='center', va='center')

        #
        plt.xlabel('differential modification')
        plt.ylabel('differential expression')

        plt.xticks([-1, 0, 1])
        plt.yticks([-1, 0, 1])
        plt.ylim(-1.5, 1.5)
        name = 'cell.A={0}_cell.B={1}_diff={2}_h={3}_r={4}_type={5}.png'.format(*task[0], *task[1:-2], task[-1])
        sns.despine(left=True)
        ax = plt.gca()
        ax.xaxis.grid()
        plt.savefig(os.path.join(dir_path, name), bbox_inches='tight')
        plt.clf()


def map_point(x, y, region, **kwargs):
    """
    Annotates a plot with specific points based on x, y coordinates and region.
    :param x: x-coordinates of points
    :param y: y-coordinates of points
    :param region: list of regions belonging to each point
    :param kwargs: additional arguments for customization
    """
    zip_y_annotate = list(dict.fromkeys(zip(x, y, region)))
    x_locs, x_labels = plt.xticks()
    x_labels = [t.get_text() for t in x_labels]
    x_label_loc_dict = dict(zip(x_labels, x_locs))

    x_unique = [e[0] for e in zip_y_annotate]
    y_unique = [e[1] for e in zip_y_annotate]
    r_unique = [e[2] for e in zip_y_annotate]
    x_pos_valid = [x_label_loc_dict[n] for n in x_unique]

    zip_y_annotate = list(dict.fromkeys(zip(x_pos_valid, y_unique, r_unique)))
    already_done = set()
    n = len(x_labels)

    for x, y, r in zip_y_annotate:
        if r == 'promoter':
            already_done.add(x)
            x -= 0.26 * (n / 4)
        else:
            x += 0.14 * (n / 4)

        plt.annotate('x', xy=(x, y), color='r', weight='bold')


def box_plot(df: pd.DataFrame, file_path: str, x: str, y: str, hue: str, col='', wrap=1, order=None):
    """
    Creates a boxplot from a Pandas dataframe and saves it to a file.
    :param df: pandas dataframe containing the plot data
    :param file_path: file path to save the plot
    :param x: column name for the x-axis
    :param y: column name for the y-axis
    :param hue: column name for the hue parameter
    :param col: column name to create subplot columns
    :param wrap: number of subplots per row
    :param order: order of categories on the x-axis
    """
    sns.reset_defaults()
    sns.set(style='whitegrid', font_scale=1.1)

    if df.empty:
        return
    if col:
        if order:
            g = sns.catplot(x=x, y=y, hue=hue, col=col, data=df, kind='box', col_wrap=wrap, order=order,
                            legend=False, showmeans=True)
        else:
            g = sns.catplot(x=x, y=y, hue=hue, col=col, data=df, kind='box', col_wrap=wrap, legend=False,
                            showmeans=True)
        g.map(map_point, x, 'original', 'region')
        g.set_xlabels('')
        g.set_ylabels('r')
        fig = g.fig
    else:
        box = sns.boxplot(x=x, y=y, data=df, hue=hue, linewidth=2.5)

        fig = box.figure
    fig.savefig(file_path, bbox_inches='tight', pad_inches=0)
    fig.clf()


def correlation_boxplot(manager: Manager, transitions: list[Differential]):
    """
    Creates and saves correlation boxplots for different motifs, regions, and histones.
    :param manager: Manager instance to handle data and configurations
    :param transitions: list of Differential instances representing transitions
    """
    name_dict = {'H1.hESC': 'hESC',
                 'GM23338': 'PST',
                 'common.myeloid.progenitor.CD34.positive': 'MP',
                 'CD14.positive.monocyte': 'MC',
                 'mesenchymal.stem.cell': 'MS',
                 'osteoblast': 'OB',
                 'neural.stem.progenitor.cell': 'NSP',
                 'neural.progenitor.cell': 'NP',
                 'bipolar.neuron': 'BN'}
    
    long_name_dict = {'H1.hESC': 'H1-hESC',
                      'GM23338': 'GM23338',
                      'common.myeloid.progenitor.CD34.positive': 'myeloid progenitor',
                      'mesenchymal.stem.cell': 'mysenchymal stem cell',
                      'osteoblast': 'osteoblast',
                      'neural.stem.progenitor.cell': 'neural stem progenitor',
                      'neural.progenitor.cell': 'neural progenitor',
                      'bipolar.neuron': 'bipolar neuron',
                      'CD14.positive.monocyte': 'monocyte',}

    motif_names = {'TF': 'TF-FFL',
                   'miRNA': 'miRNA-FFL',
                   'co-regulatory': 'co-regulatory',
                   'composite': 'composite-FFL',
                   'all': 'all'}

    os.makedirs(os.path.join(manager.result_dir, 'motif_by_histone_plots'), exist_ok=True)
    os.makedirs(os.path.join(manager.result_dir, 'motif_by_transition_plots'), exist_ok=True)
    os.makedirs(os.path.join(manager.result_dir, 'transition_by_motif_plots'), exist_ok=True)

    os.makedirs(os.path.join(manager.result_dir, 'transition_histone_plots'), exist_ok=True)
    os.makedirs(os.path.join(manager.result_dir, 'transition_motif_plots'), exist_ok=True)

    def transition_name(row: pd.Series) -> str:
        return '{0}\n>{1}'.format(name_dict[row['cell.A']], name_dict[row['cell.B']])
    
    def transition_long_name(row: pd.Series) -> str:
        return '{0} \n> {1}'.format(long_name_dict[row['cell.A']], long_name_dict[row['cell.B']])

    def long_motif_name(row: pd.Series) -> str:
        return motif_names[row['motif']]

    combinations = [('experimental', (t.cell_1, t.cell_2, False), histone, region, motif)
                    for t in transitions
                    for histone in t.histone.keys()
                    for region in manager.regions
                    for motif in ['co-regulatory', 'TF', 'miRNA']]

    corr_dict = {(*comb[1][:2], *comb[2:]): manager.correlation(*comb)[:2] for comb in combinations}
    corr_table = [(*key, correlations[0], corr) for key, correlations in corr_dict.items() for corr in correlations[1]]

    header = ['cell.A', 'cell.B', 'histone', 'region', 'motif', 'original', 'r']
    df = pd.DataFrame(data=corr_table, columns=header)

    grouped = df.groupby(['cell.A', 'cell.B', 'histone', 'region', 'motif'])['r'].describe().reset_index()
    grouped.to_csv(os.path.join(manager.result_dir, 'random_stats.tsv'), sep='\t')

    df['transition'] = df.apply(lambda row: transition_name(row), axis=1)
    df['alternative motif'] = df.apply(lambda row: long_motif_name(row), axis=1)
    df['transition long'] = df.apply(lambda row: transition_long_name(row), axis=1)

    for t in transitions:
        print(t.cell_1, t.cell_2)
        sub_df = df[(df['cell.A'] == t.cell_1) & (df['cell.B'] == t.cell_2)].copy()
        sub_df['motif'] = sub_df['alternative motif']

        order = ['H2AFZ', 'H3K4me2', 'H3K4me3', 'H3K27ac']
        file_path = os.path.join(manager.result_dir, 'transition_by_motif_plots',
                                 'cell.A={0}_cell.B={1}_edge=experimental_diff=False.png'.format(t.cell_1, t.cell_2))
        box_plot(sub_df, file_path, x='histone', y='r', hue='region', col='motif', wrap=1,
                 order=order)


def p_value_per_network(manager: Manager, transitions: list[Differential]):
    """
    Calculates and adjusts p-values for each transition-based network.
    :param manager: Manager instance for managing computation and storage
    :param transitions: list of Differential instances representing transitions
    """
    table_dict = dict()
    for t in transitions:
        results = manager.p_value('experimental', t.cell_1, t.cell_2, False)

        p_val_dict = {(motif, histone, region): p for histone, region, motif, p in results}

        adjusted = computation.benjamini_dict(p_val_dict, manager.correlation_q, krieger=True)

        table_dict[(t.cell_1, t.cell_2)] = {k: (p_val_dict[k], adjusted[k]) for k in p_val_dict.keys()}

    table = [['cell.A', 'cell.B', 'motif', 'histone', 'region', 'p.greater', 'significant']]

    with open(os.path.join(manager.result_dir, 'adjusted_p.tsv'), 'w') as file:
        for t, p_dict in sorted(table_dict.items()):
            for k, p_vals in sorted(p_dict.items()):
                table.append([*t, *k, *p_vals])

        file.write('\n'.join(['\t'.join([str(x) for x in line]) for line in table]))


def differential_genes(manager: Manager):
    """
    Extracts and saves lists of differentially expressed genes for all transitions.
    :param manager: Manager instance for accessing and saving data
    """
    directory = os.path.join(manager.result_dir, 'differential_genes')
    os.makedirs(directory, exist_ok=True)

    for a, b, m, genes in manager.all_genes_list():
        file_name = '_'.join([a, b, '.'.join(m.split())]) + '.txt'
        with open(os.path.join(directory, file_name), 'w') as file:
            file.write('\n'.join(genes.split(';')))


def pseudogenes(manager: Manager, reference: Reference, transitions: list[Differential]):
    """
    Analyzes and categorizes pseudogenes across transitions.
    :param manager: Manager instance for accessing configurations.
    :param reference: Reference instance containing gene information
    :param transitions: list of Differential instances for transition analysis
    """
    overview = defaultdict(lambda: defaultdict(set))

    gene_ids = reference.genes.keys()
    pseudo = {gene_id: gene.type for gene_id, gene in reference.genes.items() if 'pseudo' in gene.type}

    for t in transitions:
        k = (t.cell_1, t.cell_2)
        with open(t.dea_result_path, 'r') as file:
            for line in dictReader(file, delimiter='\t'):
                na = False
                if line['padj'] in ['NA']:
                    na = True

                gene_id = line['name']
                if gene_id not in gene_ids:
                    gene_id = None

                na_fc = False
                if line['log2FoldChange'] in ['NA']:
                    na_fc = True

                if na and na_fc and not gene_id:
                    overview[k]['NA.exp + NA.fc + NA.ID'].add(line['name'])
                elif na and na_fc and gene_id:
                    overview[k]['NA.exp + NA.fc'].add(line['name'])
                elif na and not na_fc and not gene_id:
                    overview[k]['NA.exp + NA.ID'].add(line['name'])
                elif na and not na_fc and gene_id:
                    overview[k]['NA.exp'].add(line['name'])
                elif not na and not na_fc and gene_id:
                    overview[k]['successful'].add(line['name'])

                if (na or na_fc) and gene_id:
                    if gene_id in pseudo:
                        overview[k]['NA.exp or NA.fc + pseudogene'].add(gene_id)
                    else:
                        overview[k]['NA.exp or NA.fc not pseudogene'].add(gene_id)

    with open(os.path.join(manager.result_dir, 'no_diff_exp.txt'), 'w') as file:
        file.write('pseudogenes: {0:,} ({1:.2f}%)\n\n'.format(len(pseudo), (len(pseudo) / len(gene_ids)) * 100))

        for t, d in sorted(overview.items()):
            file.write(' -- '.join(t) + '\n')

            for k, ids in sorted(d.items()):
                file.write('\t{0:<30} - {1:>6,} ({2:.2f}%)\n'.format(k, len(ids), (len(ids) / len(gene_ids)) * 100))


def motif_genes(manager: Manager):
    """
    Saves detailed gene lists for motifs including TFs, miRNAs, and target genes.
    :param manager: Manager instance for retrieving motif-related data
    """
    def fill_up(longest, l):
        l = l.split(';')
        return l + [''] * (longest - len(l))

    directory = os.path.join(manager.result_dir, 'motif_genes')
    os.makedirs(directory, exist_ok=True)

    for a, b, m, n, genes, tfs, mirnas, targets in manager.motif_genes():
        file_name = '_'.join([a, b, '.'.join(m.split())]) + '.tsv'

        genes, tfs, mirnas, targets = [fill_up(n, x) for x in [genes, tfs, mirnas, targets]]

        with open(os.path.join(directory, file_name), 'w') as file:
            file.write('\t'.join(['all', 'TFs', 'miRNAs', 'targets']) + '\n')

            for i in range(n):
                file.write('\t'.join([genes[i], tfs[i], mirnas[i], targets[i]]) + '\n')


def quant_vs_cat(manager: Manager, transition: list[Differential]):
    """
    Compares categorical and quantitative correlations across transitions.
    :param manager: Manager instance for handling data and computation
    :param transition: Differential instance representing a transition
    """
    tasks = [((t.cell_1, t.cell_2, diff), histone, region, gene_type)
             for t in transition
             for diff in manager.differential
             for histone in t.histone.keys()
             for region in manager.regions
             for gene_type in ['', 'TF', 'miRNA', 'other']]

    table = [['cell.A', 'cell.B', 'diff.', 'histone', 'region', 'type', 'r.cat', 'r.quant', 'change',
              'abs_change', 'cat_stronger']]

    for task in tasks:
        gene_ids = manager.gene_ids(*task[0], task[-1])

        cat = manager.pairs(*task[:-1], 'mod_diff', gene_ids)
        quant = manager.pairs_quant(*task[:-1], gene_ids)

        r_cat = nan_round(computation.correlation(10, cat)[0], 3)
        r_quant = nan_round(computation.correlation(10, quant)[0], 3)

        table.append([*task[0], *task[1:], r_cat, r_quant, r_cat - r_quant, abs(r_cat) - abs(r_quant),
                      abs(r_cat) >= abs(r_quant)])

    with open(os.path.join(manager.result_dir, 'cat_vs_quant.tsv'), 'w') as file:
        file.write('\n'.join(['\t'.join([str(x) for x in line]) for line in table]))

    stronger = [row[-1] for row in table[1:] if not computation.is_nan(row[6]) and not computation.is_nan(row[7])]
    print('stronger:', stronger.count(True))
    print('weaker:', stronger.count(False))
    print('total:', len(stronger))


def p_fisher_diff2(manager: Manager):
    """
    Computes Fisher transformations and correlations for motifs vs. complements.
    :param manager: Manager instance for data and configuration
    """
    corr_motifs = manager.p_comparison_all2('experimental', False)

    corr_dict = dict()
    for key, res in corr_motifs.items():
        genes, r, n = res
        cell_1, cell_2, histone, region, motif = key

        motif_gene_ids = set(genes.split(';'))
        all_gene_ids = manager.gene_ids(cell_1, cell_2, False)
        diff_gene_ids = manager.gene_ids(cell_1, cell_2, True)
        complement_all = all_gene_ids.difference(motif_gene_ids)
        complement_diff = diff_gene_ids.difference(motif_gene_ids)

        corr_all = computation.correlation(10, manager.pairs((cell_1, cell_2, False),
                                                             histone, region, 'mod_diff', complement_all))[0]
        corr_diff = computation.correlation(10, manager.pairs((cell_1, cell_2, True),
                                                              histone, region, 'mod_diff', complement_diff))[0]

        n_all = len(complement_all)
        n_diff = len(complement_diff)

        corr_dict[(*key, 'all')] = (r, n, corr_all, n_all)
        corr_dict[(*key, 'diff.')] = (r, n, corr_diff, n_diff)

    p_vals = {key: computation.fisher_correlation_comparison(*v) for key, v in corr_dict.items()}

    change = {key: v[0] - v[2] if not computation.is_nan(v[0]) and not computation.is_nan(v[2]) else np.nan
              for key, v in corr_dict.items()}

    adjusted = computation.benjamini_dict(p_vals, manager.correlation_q, krieger=True)

    with open(os.path.join(manager.result_dir, 'fisher_motif_vs_complement.tsv'), 'w') as file:
        table = [['cell.A', 'cell.B', 'histone', 'region', 'motif', 'comparison', 'r.motif', 'n.motif',
                  'r.wide', 'n.wide', 'change', 'p', 'sig.']]

        for key, corr in sorted(corr_dict.items()):
            r_1, n_1, r_2, n_2 = corr
            table.append([*key, nan_round(r_1, 3), n_1, nan_round(r_2, 3), n_2, nan_round(change[key], 3), p_vals[key],
                          adjusted[key]])

        file.write('\n'.join(['\t'.join([str(x) for x in line]) for line in table]))


def p_fisher_diff(manager: Manager):
    """
    Performs Fisher transformation-based correlation comparisons between motifs and their complements.
    :param manager: Manager instance to handle computations and data access
    """
    correlations = manager.p_comparison_diff('experimental')
    corr_motifs = manager.p_comparison_all('experimental', False)
    corr_motifs_diff = {key: (*v, *correlations[key[:-1] + ('all genes',)][2:]) for key, v in corr_motifs.items()}
    corr_motifs_all = {key: (*v, *correlations[key[:-1] + ('all genes',)][:2]) for key, v in corr_motifs.items()}

    tasks = [('fisher_diff_vs_non_diff.tsv', 'all', 'diff', correlations),
             ('fisher_motif_vs_all.tsv', 'motif', 'all', corr_motifs_all),
             ('fisher-motif_vs_diff.tsv', 'motif', 'diff', corr_motifs_diff)]

    for name, a, b, dic in tasks:
        p_vals = {key: computation.fisher_correlation_comparison(*v) for key, v in dic.items()}

        change = {key: v[2] - v[0] if not computation.is_nan(v[0]) and not computation.is_nan(v[2]) else np.nan
                  for key, v in dic.items()}

        adjusted = computation.benjamini_dict(p_vals, manager.correlation_q, krieger=True)

        with open(os.path.join(manager.result_dir, name), 'w') as file:
            table = [['cell.A', 'cell.B', 'histone', 'region', 'motif', 'r.' + a, 'n.' + a,
                      'r.' + b, 'n.' + b, 'change', 'p', 'sig.']]

            for key, corr in sorted(dic.items()):
                r_1, n_1, r_2, n_2 = corr
                table.append([*key, nan_round(r_1, 3), n_1, nan_round(r_2, 3), n_2, nan_round(change[key], 3), p_vals[key],
                              adjusted[key]])

            file.write('\n'.join(['\t'.join([str(x) for x in line]) for line in table]))


def p_steiger(manager: Manager):
    """
    Uses Steiger's statistical test for comparing two overlapping correlation coefficients from the same sample to
    compared expression-histone correlations between promoter region and gene body.
    """
    correlation_dir = dict()

    for diff in manager.differential:
        for key, mod in manager.p_comparison_region('experimental', diff).items():
            correlation_dir[(diff, *key)] = (*mod[:-1], computation.correlation(10, mod[-1])[0])

    change_dir = {key: v[2] - v[1] if not computation.is_nan(v[0]) and not computation.is_nan(v[2]) else np.nan
                  for key, v in correlation_dir.items()}

    p_val_dir = {key: computation.steiger(*v) for key, v in correlation_dir.items()}

    adjusted = computation.benjamini_dict(p_val_dir, manager.correlation_q, krieger=True)

    with open(os.path.join(manager.result_dir, 'steiger_gene_vs_promoter.tsv'), 'w') as file:
        table = [['diff', 'cell.A', 'cell.B', 'histone', 'motif', 'n', 'r.gene', 'r.prom', 'r.mod', 'change',
                  'p', 'sig.']]

        for key, corr in sorted(correlation_dir.items()):
            n, r_1, r_2, r_3 = corr
            table.append([*key, n, nan_round(r_1, 3), nan_round(r_2, 3), nan_round(r_3, 3), nan_round(change_dir[key], 3),
                          p_val_dir[key], adjusted[key]])

        file.write('\n'.join(['\t'.join([str(x) for x in line]) for line in table]))


def motif_overview(cfg: dict, manager: Manager, transitions: list[Differential]):
    """
    Generates a comprehensive overview of motifs for different transitions and edge types.
    :param cfg: run specific configuration dictionary
    :param manager: Manager instance containing motif, edge type, and differential data
    :param transitions: list of Differential objects representing transitions
    """
    start = time()
    n_random = int(cfg['analysis']['number networks'])                                      # type: int

    for edge_type in manager.edge_types:
        file_path = cfg['file paths']['analysis']['motif overview'].format(edge_type)
        table = [[''] * 3, ['cell.A', 'cell.B', 'motif']]                                   # type: list[list[str]]

        for diff in manager.differential:
            table[-2] += ['differential' if diff else 'all'] + [''] * 6
            table[-1] += ['actual.n', 'p.smaller', 'p.greater', 'random.min', 'random.median', 'random.mean', 'random.max']

        for t in transitions + [None]:
            first = True

            for motif in manager.motif_names + ['all']:
                if first:
                    cols = [t.cell_1, t.cell_2] if t else ['complete', '']
                    first = False
                else:
                    cols = ['', '']

                cols.append(motif)

                for diff in manager.differential:
                    specific = (t.cell_1, t.cell_2, diff) if t else None

                    if not manager.motif_table_exists(edge_type, specific, -1):
                        cols += [np.NaN] * 7
                        continue

                    n_actual = len(manager.motif(edge_type=edge_type, specific=specific,
                                                 network_number=-1,
                                                 motif=motif))
                    r_motifs = [len(manager.motif(edge_type=edge_type, specific=specific,
                                                  network_number=r,
                                                  motif=motif))
                                for r in range(manager.n_random)]

                    p_smaller = len([1 for r in r_motifs if r <= n_actual]) / n_random
                    p_greater = len([1 for r in r_motifs if r >= n_actual]) / n_random

                    cols += [n_actual, p_smaller, p_greater, *computation.list_stats(r_motifs)[3:]]
                table.append(cols)

        util.write_table(file_path, table, 2, range(3, 17))

    util.print_time(start, 1)


def number_of_transcripts(file_path: str):
    """
    Counts the number of unique transcripts in a FASTA file.
    :param file_path: path to a FASTA file
    """
    with open(file_path, 'r') as file:
        print('transcripts:', len({line for line in file if line.startswith('>')}))


def number_of_gene_ids_per_gene_name(file_path: str):
    """
    Analyzes the mapping of gene IDs to gene names in a tab-delimited file and computes various statistics.
    :param file_path: path to the result file
    """
    with open(file_path, 'r') as file:
        d = defaultdict(set)
        gene_lengths = []

        for line in dictReader(file, delimiter='\t'):
            d[line['gene.name']].add(line['gene.id'])
            gene_lengths.append(int(line['gene.length']))

    print('gene names:', len(d.keys()))
    print('gene IDs:', len({gene_id for ids in d.values() for gene_id in ids}))
    print('gene names with more than one ID:', len({gene_name for gene_name, ids in d.items() if len(ids) > 1}))
    max_ids = max(len(ids) for ids in d.values())
    print('most IDs per name:', max_ids)
    print('names with max IDs:', {gene_name for gene_name, ids in d.items() if len(ids) == max_ids})
    sorted_pairs = sorted([(len(ids), name) for name, ids in d.items()])
    lengths = [len(ids) for ids in d.values()]
    print('highest:', sorted_pairs[-10:])
    print('stats:', computation.list_stats(lengths, 2), stdev(lengths))
    print('gene lengths:', computation.list_stats(gene_lengths, 2), stdev(gene_lengths))
