from collections import defaultdict
import numpy as np
import os
import pandas as pd
import sqlite3
from typing import Iterable, Any

from gene_functionality import GeneFunctionIndex
import util


class Manager:
    """
    A class to manage data preprocessing, network computation, and database handling for gene expression and histone
    modification analysis.
    """
    def __init__(self, cfg: dict, reload_results: bool, recompute_networks: bool):
        """
        Initializes the Manager class and sets up directory structures based on the configuration.
        :param cfg: configuration dictionary containing paths and analysis parameters
        :param reload_results: flag to determine if existing results should be reloaded
        :param recompute_networks: flag to determine if networks should be recomputed
        """
        # base directory for processing data
        self.pre_processing_dir = cfg['file paths']['analysis']['preprocessing dir']
        # directory for storing pickle files
        self.pickle_dir = os.path.join(self.pre_processing_dir, 'pickle')
        # directories for storing network edge, node nad motif data
        self.edge_dir = os.path.join(self.pre_processing_dir, 'networks')
        self.node_dir = os.path.join(self.pre_processing_dir, 'nodes')
        self.motif_dir = os.path.join(self.pre_processing_dir, 'motifs')
        # directories for storing analysis data
        self.correlation_dir = os.path.join(self.pre_processing_dir, 'correlation')
        self.network_stats_dir = os.path.join(self.pre_processing_dir, 'network.stats')
        self.network_stats_fmt = os.path.join(self.network_stats_dir, 'cell.A={0}_cell.B={1}_diff={2}_edges={3}.txt')
        self.histogram_dir = os.path.join(self.pre_processing_dir, 'histograms')
        # directories for gene type and expression data
        self.gene_type_path = os.path.join(self.pickle_dir, 'gene_types.tsv')
        self.expression_dir = os.path.join(self.pre_processing_dir, 'results')
        # directory for results
        self.result_dir = cfg['file paths']['analysis']['result dir']

        # number of random networks to generate
        self.n_random = int(cfg['analysis']['number networks'])                             # type: int
        # FDR and p-value thresholds for correlation analysis
        self.correlation_q = float(cfg['analysis']['correlation fdr'])                      # type: float
        self.correlation_p = float(cfg['analysis']['correlation p'])                        # type: float
        # PDR and p-value thresholds for motif analysis
        self.motif_q = float(cfg['analysis']['motif fdr'])                                  # type: float
        self.motif_p = float(cfg['analysis']['motif p'])                                    # type: float

        # minimum number required for correlation analysis, general and for composite motifs
        self.correlation_min = int(cfg['analysis']['minimum number'])                       # type: int
        self.correlation_min_composite = int(cfg['analysis']['minimum number composite'])   # type: int

        # list of edge types, motif names, genomic regions, differential types, and interaction types to consider
        self.edge_types = cfg['analysis']['manager']['edge types']
        self.motif_names = cfg['analysis']['manager']['motifs']
        self.regions = cfg['analysis']['manager']['regions']
        self.differential = cfg['analysis']['manager']['differential']
        self.interactions = cfg['analysis']['manager']['interaction types']

        # dictionary to store edge counts
        self.n_edges = dict()
        # name of the expression table in the database
        self.expression_table = 'expression'

        # prepare the configured directories
        util.delete_create_dir(self.edge_dir, delete=recompute_networks)
        util.delete_create_dir(self.node_dir, delete=recompute_networks)
        util.delete_create_dir(self.motif_dir, delete=recompute_networks)
        util.delete_create_dir(self.correlation_dir, delete=recompute_networks)
        util.delete_create_dir(self.network_stats_dir, delete=recompute_networks)
        util.delete_create_dir(self.histogram_dir, delete=recompute_networks)

        util.delete_create_dir(self.expression_dir, delete=reload_results)
        util.delete_create_dir(self.pickle_dir, delete=reload_results)

        util.delete_create_dir(self.result_dir, delete=True)

    def transition_path(self, cell_1: str, cell_2: str) -> str:
        """ Constructs the file path for the transition database between two cells. """
        return os.path.join(self.expression_dir, '_'.join([cell_1, cell_2]) + '.db')

    @staticmethod
    def histone_table(histone: str, region: str) -> str:
        """ Constructs the table name for a given histone and genomic region. """
        return '_'.join([histone, region])

    @staticmethod
    def table_exists(cursor, table_name: str) -> bool:
        """ Checks if a table exists in the database. """
        try:
            cursor.execute('SELECT 1 from {0} LIMIT 1'.format(table_name))
            return True
        except sqlite3.OperationalError:
            return False

    def histone_tables_exist(self, cell_1: str, cell_2: str, histones: Iterable[str]) -> bool:
        """
        Checks if all histone modification tables exist in the database.
        :param cell_1: name of the first cell type
        :param cell_2: name of the second cell type
        :param histones: list of histones to check
        :return: True if all tables exist, False otherwise
        """
        connection = sqlite3.connect(self.transition_path(cell_1, cell_2))
        cursor = connection.cursor()

        for histone in histones:
            for region in self.regions:
                if not self.table_exists(cursor, self.histone_table(histone, region)):
                    connection.close()
                    return False

        connection.close()
        return True

    def to_db(self, index: GeneFunctionIndex):
        """
        Writes gene expression and histone modification data to the database.
        :param index: index containing expression and modification data
        """
        # establish database connection and prepare cursor
        connection = sqlite3.connect(self.transition_path(index.cell_1, index.cell_2))
        cursor = connection.cursor()
        cursor.execute('PRAGMA synchronous = OFF')
        cursor.execute('PRAGMA journal_mode = MEMORY')
        cursor.execute('BEGIN TRANSACTION')

        # SQL command to create the expression table
        create_cmd = 'CREATE TABLE {0} (id TEXT PRIMARY KEY, gene_type TEXT, is_diff INTEGER)'
        # command to insert gene expression data
        add_cmd = 'INSERT INTO {0} (id, gene_type, is_diff) VALUES (?, ?, ?)'.format(self.expression_table)

        # create and populate the expression table
        self.create_table(self.expression_table, create_cmd, cursor=cursor, drop=True)
        cursor.executemany(add_cmd.format(self.expression_table),
                           sorted([(k, e.type, e.is_diff) for k, e in index.exp.items()]))

        # commands for creating and populating histone tables
        create_cmd = 'CREATE TABLE {0} (id TEXT PRIMARY KEY, gene_type TEXT, is_diff INTEGER, fold_change REAL, ' \
                     'p_value REAL, q_value REAL, exp_cat INTEGER, mod_A INTEGER, mod_B INTEGER, unmod INTEGER, ' \
                     'mod_both INTEGER, mod_diff INTEGER, mod_diff_2 INTEGER, mod_cat INTEGER, mod_cat_2 INTEGER)'

        add_cmd = 'INSERT INTO {0} (id, gene_type, is_diff, fold_change, p_value, q_value,  exp_cat, mod_A, mod_B, ' \
                  'unmod, mod_both, mod_diff, mod_diff_2, mod_cat, mod_cat_2) VALUES ' \
                  '(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)'

        # build the tables for all genomic region and histone modification combinations
        for region, mod_index in index.indices.items():
            for histone, mod_dict in mod_index.mods.items():
                table_dict = {gene_id: (gene_id, exp.type, exp.is_diff, exp.fc, exp.p, exp.q, exp.cat)
                              for gene_id, exp in index.exp.items()}

                for gene_id, mc in mod_dict.items():
                    table_dict[gene_id] += (mc.A, mc.B, mc.unmod(), mc.both(), mc.diff(), mc.diff_2(),
                                            mc.cat(), mc.cat_2())

                table_name = self.histone_table(histone, region)
                self.create_table(table_name, create_cmd, cursor=cursor, drop=True)
                cursor.executemany(add_cmd.format(table_name), sorted(table_dict.values()))

        # commit transaction and close the connection
        connection.commit()
        connection.close()

    def gene_ids(self, cell_1: str, cell_2: str, differential: bool, gene_type='') -> set[str]:
        """
        Retrieves a set of gene IDs based on specific criteria.
        :param cell_1: name of the first cell type
        :param cell_2: name of the second cell type
        :param differential: whether to filter for differentially expressed genes
        :param gene_type: specific gene type to filter for
        :return: set of gene IDs that match the criteria
        """
        connection = sqlite3.connect(self.transition_path(cell_1, cell_2))
        cursor = connection.cursor()

        # construct the base SQL query
        cmd = 'SELECT id FROM {0}'.format(self.expression_table)
        to_add = []

        # add conditions to the query based on the parameters
        if differential:
            to_add.append('is_diff=1')
        if gene_type:
            to_add.append('gene_type="{0}"'.format(gene_type))

        if to_add:
            cmd += ' WHERE ' + ' AND '.join(to_add)

        # execute the query and fetch results
        cursor.execute(cmd)
        result = {x[0] for x in cursor.fetchall()}

        connection.close()
        return result

    def pairs(self, specific: tuple[str, str, bool], histone: str, region: str, mod: str,
              gene_ids: set[str]) -> list[tuple[int, int]]:
        """
        Retrieves pairs of expression and modification data for specified genes.
        :param specific: tuple containing two cell types and a differential flag
        :param histone: histone modification identifier
        :param region: genomic region
        :param mod: modification type to retrieve
        :param gene_ids: set of gene IDs to include
        :return: list of (expression category, modification value) pairs
        """
        cell_1, cell_2, differential = specific
        connection = sqlite3.connect(self.transition_path(cell_1, cell_2))
        cursor = connection.cursor()

        table_name = self.histone_table(histone, region)

        # construct the SQL query based on the differential flag
        if differential:
            cmd = 'SELECT id, exp_cat, {0} FROM {1} WHERE is_diff=1'.format(mod, table_name)
        else:
            cmd = 'SELECT id, exp_cat, {0} FROM {1}'.format(mod, table_name)

        # execute the query and filter results by gene IDs
        cursor.execute(cmd)
        results = cursor.fetchall()
        connection.close()

        return [(exp, mod) for gene_id, exp, mod in results if gene_id in gene_ids]

    def pairs_quant(self, specific: tuple[str, str, bool], histone: str, region: str,
                    gene_ids: set[str]) -> list[tuple[float, int]]:
        """
        Retrieves pairs of fold-change and modification differences for specified genes.
        :param specific: tuple containing two cell types and a differential flag
        :param histone: histone modification identifier
        :param region: genomic region
        :param gene_ids: set of gene IDs to include
        :return: list of (fold-change, modification difference) pairs
        """
        cell_1, cell_2, differential = specific
        connection = sqlite3.connect(self.transition_path(cell_1, cell_2))
        cursor = connection.cursor()

        table_name = self.histone_table(histone, region)

        # construct the SQL query based on the differential flag
        if differential:
            cmd = 'SELECT id, fold_change, exp_cat, mod_A, mod_B FROM {0} WHERE is_diff=1'.format(table_name)
        else:
            cmd = 'SELECT id, fold_change, exp_cat, mod_A, mod_B FROM {0} '.format(table_name)

        # execute the query, filter by gene IDs and calculate differences
        cursor.execute(cmd)
        results = cursor.fetchall()
        connection.close()

        return [(fc, B - A) if exp else (0.0, B - A) for gene_id, fc, exp, A, B in results if gene_id in gene_ids]

    def histone_pairs(self, specific: tuple[str, str, bool], histones: list[str], regions: list[str], mod: str,
                      gene_ids: set[str]) -> list[tuple[int, int]]:
        """
        Retrieves pairs of modification values for specified histones and regions.
        :param specific: tuple containing two cell types and a differential flag
        :param histones: list of histone modification identifiers to include
        :param regions: list of genomic regions to include
        :param mod: modification type to retrieve
        :param gene_ids: set of gene IDs to include
        :return: list of modification pairs
        """
        # validate the combination of histones and regions
        if len(histones) not in [1, 2] or len(regions) not in [1, 2] or (len(histones) + len(regions)) != 3:
            raise ValueError('Wrong combination of histones and regions.')

        cell_1, cell_2, differential = specific
        connection = sqlite3.connect(self.transition_path(cell_1, cell_2))
        cursor = connection.cursor()

        # retrieve modifications for each histone and region combination
        modifications = []
        for histone in histones:
            for region in regions:
                table_name = self.histone_table(histone, region)

                if differential:
                    cmd = 'SELECT id, {0} FROM {1} WHERE is_diff=1'.format(mod, table_name)
                else:
                    cmd = 'SELECT id, {0} FROM {1}'.format(mod, table_name)

                cursor.execute(cmd)
                results = cursor.fetchall()

                modifications.append({gene_id: mod for gene_id, mod in results})

        # create pairs of modifications for matching gene IDs
        connection.close()
        d_1, d_2 = modifications

        return [(d_1[gene_id], d_2[gene_id]) for gene_id in gene_ids]

    def create_table(self, table_name: str, setup_cmd: str, db_path=None, cursor=None, drop=False) -> bool:
        """
        Creates a database table if it does not already exist.
        :param table_name: name of the table to create
        :param setup_cmd: SQL command to create the table
        :param db_path: path to the database file (optional if cursor given)
        :param cursor: database cursor (optional if database path given)
        :param drop: whether to drop the table if it already exists
        :return: True if the table was created, False otherwise
        """
        if not db_path and not cursor:
            raise ValueError('Neither database path nor cursor given.')

        # open a new connection if a database path is given
        if db_path:
            connection = sqlite3.connect(db_path)
            cursor = connection.cursor()

        # drop the table if requested
        if drop:
            cursor.execute('DROP TABLE IF EXISTS {0}'.format(table_name))

        # check if the table already exists
        if self.table_exists(cursor, table_name):
            return False

        # create the table with the given command
        cursor.execute(setup_cmd.format(table_name))

        # commit and close the connection if it was opened in this function
        if db_path:
            connection.commit()
            connection.close()

        return True

    def node_path(self, edge_type: str, specific=None) -> str:
        """
        Constructs the path for saving node data.
        :param edge_type: type of the edge
        :param specific: tuple containing two cell types and a differential flag
        :return: path to the node file
        """
        if specific:
            name = 'cell.A={0}_cell.B={1}_diff={2}_edges={3}.tsv'.format(*specific, edge_type)
        else:
            name = 'complete_edges={0}.tsv'.format(edge_type)

        # ensure the node directory exists
        os.makedirs(self.node_dir, exist_ok=True)

        return os.path.join(self.node_dir, name)

    def edge_path(self, edge_type: str, interaction_type: str, specific: Any | tuple[str, str, bool],
                  network_number: int) -> str:
        """
        Constructs the path for saving edge data.
        :param edge_type: type of the edge (e.g. experimental)
        :param interaction_type: type of the interactions (e.g. tf-mirna)
        :param specific: tuple containing two cell types and a differential flag
        :param network_number: network identifier (negative for original network, non-negative for randomized)
        :return: path to the edge file
        """
        if specific:
            dir_name = 'cell.A={0}_cell.B={1}_diff={2}_edges={3}'.format(*specific, edge_type)
        else:
            dir_name = 'complete_edges={0}'.format(edge_type)

        # determine the file name based on the network number
        if network_number >= 0:
            file_name = 'random_{0}_{1}.txt'.format(network_number, interaction_type)
        else:
            file_name = 'original_{0}.tsv'.format(interaction_type)

        # ensure the edge directory exists
        dir_path = os.path.join(self.edge_dir, dir_name)
        os.makedirs(dir_path, exist_ok=True)

        return os.path.join(dir_path, file_name)

    def motif_path(self, edge_type: str, specific: Any | tuple[str, str, bool]) -> str:
        """
        Constructs the path for saving motif data.
        :param edge_type: type of the edge (e.g. experimental)
        :param specific: tuple containing two cell types and a differential flag
        :return: path to the motif file
        """
        if specific:
            name = 'cell.A={0}_cell.B={1}_diff={2}_edges={3}'.format(*specific, edge_type)
        else:
            name = 'edges={0}'.format(edge_type)

        # ensure the motif directory exists
        os.makedirs(self.motif_dir, exist_ok=True)

        return os.path.join(self.motif_dir, name + '.db')

    @staticmethod
    def motif_table(network_number: int):
        """
        Determine the motif table name based on the network number.
        :param network_number: network identifier (negative for original network, non-negative for randomized)
        :return: name of the motif table
        """
        if network_number < 0:
            return 'original'
        return 'random_{0}'.format(network_number)

    def motif_table_exists(self, edge_type: str, specific: Any | tuple[str, str, bool],
                           network_number: int) -> bool:
        """
        Checks if the motif table exists in the motif database.
        :param edge_type: type of the edge (e.g. experimental)
        :param specific: tuple containing two cell types and a differential flag
        :param network_number: network identifier (negative for original network, non-negative for randomized)
        :return: True if the table exists, False otherwise
        """
        connection = sqlite3.connect(self.motif_path(edge_type, specific))
        cursor = connection.cursor()

        exists = self.table_exists(cursor, self.motif_table(network_number))

        connection.close()
        return exists

    def create_motif_table(self, edge_type: str, specific: Any | tuple[str, str, bool], network_number: int) \
            -> bool:
        """
        Creates a motif table in the motif database.
        :param edge_type: type of the edge (e.g. experimental)
        :param specific: tuple containing two cell types and a differential flag
        :param network_number: network identifier (negative for original network, non-negative for randomized)
        :return: True if the table was successfully created, False otherwise
        """
        cmd = 'CREATE TABLE IF NOT EXISTS {0} (TF TEXT NOT NULL, miRNA TEXT NOT NULL, motif TEXT, shared TEXT, ' \
              'PRIMARY KEY (TF, miRNA))'
        table_name = self.motif_table(network_number)
        return self.create_table(table_name, cmd, db_path=self.motif_path(edge_type, specific))

    def add_motifs(self, results: dict[str, dict[Any | tuple[str, str, bool], list[tuple[int, list[tuple]]]]]):
        """
        Adds motifs to the motif database for each network.
        :param results: dictionary containing edge types, specific details, and motif data to be added
        """
        for edge_type, content in results.items():
            for specific, to_enter in content.items():
                connection = sqlite3.connect(self.motif_path(edge_type, specific))
                cursor = connection.cursor()

                # improve performance by adjusting SQLite pragmas
                cursor.execute('PRAGMA synchronous = OFF')
                cursor.execute('PRAGMA journal_mode = MEMORY')
                cursor.execute('BEGIN TRANSACTION')

                for network_number, table in to_enter:
                    cursor.executemany('INSERT INTO {0} (TF, miRNA, motif, shared) '
                                       'VALUES (?, ?, ?, ?)'.format(self.motif_table(network_number)),
                                       table)
                connection.commit()
                connection.close()

    def motif(self, edge_type: str, specific: Any | tuple[str, str, bool], network_number: int,
              motif: str) -> list[tuple[str, str, str, str]]:
        """
        Retrieves motifs from the motif database.
        :param edge_type: type of the edge (e.g. experimental)
        :param specific: tuple containing two cell types and a differential flag
        :param network_number: network identifier (negative for original network, non-negative for randomized)
        :param motif: motif name or 'all' to retrieve all motifs
        :return: database rows with the specific motifs
        """
        table_name = self.motif_table(network_number)
        connection = sqlite3.connect(self.motif_path(edge_type, specific))
        cursor = connection.cursor()

        if motif == 'all':
            cmd = 'SELECT TF, miRNA, shared, motif FROM {0}'.format(table_name)
            cursor.execute(cmd)
        else:
            cmd = 'SELECT TF, miRNA, shared, motif FROM {0} WHERE motif=?'.format(table_name)
            cursor.execute(cmd, (motif,))

        rows = cursor.fetchall()

        connection.close()
        return rows

    def correlation_path(self, edge_type: str, specific: Any | tuple[str, str, bool], degree=False) -> str:
        """
        Constructs the path for saving correlation data.
        :param edge_type: type of the edge (e.g. experimental)
        :param specific: tuple containing two cell types and a differential flag
        :param degree: True if the correlation is for node degrees, False otherwise
        :return: path to the correlation file
        """
        if specific:
            name = 'cell.A={0}_cell.B={1}_diff={2}_edges={3}.db'.format(*specific, edge_type)
        else:
            name = 'edges={0}.db'.format(edge_type)

        if degree:
            dir_path = self.histogram_dir
        else:
            dir_path = self.correlation_dir

        return os.path.join(dir_path, name)

    @staticmethod
    def correlation_table(histone: str, region: str) -> str:
        """
        Constructs the name of the correlation table.
        :param histone: histone identifier
        :param region: genomic region
        :return: table name
        """
        return '{0}_{1}'.format(histone, region)

    def correlation_table_exists(self, edge_type: str, specific: Any | tuple[str, str, bool],
                                 histone: str, region: str, degree=False) -> bool:
        """
        Checks if the correlation table exists in the database.
        :param edge_type: type of the edge (e.g. experimental)
        :param specific: tuple containing two cell types and a differential flag
        :param histone: histone identifier
        :param region: genomic region
        :param degree: True if the correlation is for node degrees, False otherwise
        :return: True if the table exists, False otherwise
        """
        table_name = self.correlation_table(histone, region)

        connection = sqlite3.connect(self.correlation_path(edge_type, specific, degree=degree))
        cursor = connection.cursor()

        exists = self.table_exists(cursor, table_name)
        connection.close()

        return exists

    def create_correlation_table(self, edge_type: str, specific: Any | tuple[str, str, bool],
                                 histone: str, region: str, degree=False) -> bool:
        """
        Creates a correlation table in the database.
        :param edge_type: type of the edge (e.g. experimental)
        :param specific: tuple containing two cell types and a differential flag
        :param histone: histone identifier
        :param region: genomic region
        :param degree: True if the correlation is for node degrees, False otherwise
        :return: True if the table was successfully created, False otherwise
        """
        if degree:
            cmd = 'CREATE TABLE IF NOT EXISTS {0} (category TEXT, degree INTEGER, n_genes INTEGER, gene_ids TEXT, ' \
                  'n_pairs INTEGER, pearson REAL, corr_comment TEXT, chi_p REAL, chi_sig INTEGER, Cramer_v REAL, ' \
                  'corrected_v REAL, chi_comment TEXT, ' \
                  'PRIMARY KEY (category, degree))'
        else:
            cmd = 'CREATE TABLE IF NOT EXISTS {0} (motif TEXT, network_number INTEGER, n_genes INTEGER, ' \
                  'n_motifs INTEGER, pearson REAL, corr_comment TEXT, n_TFs INTEGER, n_miRNAs INTEGER, ' \
                  'n_targets INTEGER, gene_list TEXT, motif_list TEXT, TF_list TEXT, miRNA_list TEXT, ' \
                  'target_list TEXT, PRIMARY KEY (motif, network_number))'

        table_name = self.correlation_table(histone, region)

        return self.create_table(table_name, cmd, db_path=self.correlation_path(edge_type, specific, degree=degree),
                                 drop=True)

    def add_correlation(self, results: dict[str, dict[tuple[tuple[str, str, bool], str, str], list[tuple]]],
                        degree=False):
        """
        Adds correlation data to the correlation database.
        :param results: dictionary containing edge types, specific details, histones, and regions, with associated data
        :param degree: True if the correlation is for node degrees, False otherwise
        """
        for edge_type, content in results.items():
            for key, table in content.items():
                specific, histone, region = key
                connection = sqlite3.connect(self.correlation_path(edge_type, specific, degree=degree))
                cursor = connection.cursor()

                # improve performance by adjusting SQLite pragmas
                cursor.execute('PRAGMA synchronous = OFF')
                cursor.execute('PRAGMA journal_mode = MEMORY')
                cursor.execute('BEGIN TRANSACTION')

                self.create_correlation_table(edge_type, specific, histone, region, degree=degree)

                if degree:
                    cmd = 'INSERT INTO {0} (category, degree, n_genes, gene_ids, n_pairs, pearson, ' \
                          'corr_comment, chi_p, Cramer_v, corrected_v, chi_comment, chi_sig)' \
                          ' VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)'.format(self.correlation_table(histone, region))
                else:
                    cmd = 'INSERT INTO {0} (motif, network_number, n_genes, n_motifs, n_TFs, n_miRNAs, n_targets, ' \
                          'gene_list, motif_list, TF_list, miRNA_list, target_list, pearson, corr_comment) ' \
                          'VALUES ' \
                          '(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)'.format(self.correlation_table(histone, region))

                cursor.executemany(cmd, table)
                connection.commit()
                connection.close()

    def correlation(self, edge_type: str, specific: tuple[str, str, bool], histone: str, region: str,
                    motif: str) -> tuple[float, list[float], str, int, int, int, int, int, str, str, str, str, str]:
        """
        Fetches the correlation data for a specified motif, histone, and region.
        :param edge_type: type of the edge (e.g. experimental)
        :param specific: tuple containing two cell types and a differential flag
        :param histone: histone identifier
        :param region: genomic region
        :param motif: motif whose correlation data is being queried
        :return: correlation of the original and randomized networks as well information about the original network
        """
        # connect to the SQLite database for the given edge type and condition
        connection = sqlite3.connect(self.correlation_path(edge_type, specific))
        cursor = connection.cursor()

        # check if the required table exists, return default if it doesn't
        if not self.table_exists(cursor, self.correlation_table(histone, region)):
            return (np.nan, [], 'table does not exist') + (0,) * 5 + ('n.a.',) * 5

        # fetch correlation data for the specified motif
        cmd = 'SELECT network_number, pearson, corr_comment, n_genes, n_motifs, n_TFs, n_miRNAs, n_targets, ' \
              'gene_list, motif_list, TF_list, miRNA_list, target_list ' \
              'FROM {0} WHERE motif="{1}"'.format(self.correlation_table(histone, region), motif)

        cursor.execute(cmd)
        results = {row[0]: row for row in cursor.fetchall()}
        connection.close()

        # return default if no data for the original network
        if -1 not in results:
            return (np.nan, [], 'no data for the original network') + (0,) * 5 + ('n.a.',) * 5

        # extract tuple data from the original network
        x = results[-1][-11:]   # type: tuple[str, int, int, int, int, int, str, str, str, str, str]

        # handle cases where the original correlation is invalid
        if results[-1][1] is None or np.isnan(results[-1][1]):
            return (np.nan, []) + x

        # handle cases where the correlation computation failed
        if results[-1][2] != 'success':
            return (np.nan, []) + x

        # store original correlation value
        original_r = results[-1][1]
        # collect valid random correlation values
        random_r = [results[i][1] for i in range(self.n_random) if results[i] and results[i][2] == 'success']
        random_r = [r for r in random_r if r is not None and not np.isnan(r)]

        # handle cases where there are no valid random correlations
        if not random_r:
            x = x[1:]           # type: tuple[int, int, int, int, int, str, str, str, str, str]
            return (original_r, [], 'no valid random correlation') + x

        # return the original correlation, random correlations, and information about the original network
        return (original_r, random_r) + x

    def p_value_path(self, edge_type: str, differential: bool) -> str:
        """
        Constructs the path for saving p-value data.
        :param edge_type: type of the edge (e.g. experimental)
        :param differential: whether to filter for differentially expressed genes
        :return: path to the p-value file
        """
        return os.path.join(self.pre_processing_dir, 'p.values_edges={0}_diff={1}.db'.format(edge_type, differential))

    @staticmethod
    def p_value_table():
        """ Returns the name of the table that stores p-values. """
        return 'p_values'

    def create_p_value_table(self, cursor) -> bool:
        """
        Creates a p-value table in the database if it doesn't already exist.
        :param cursor: the database cursor
        :return: True if the table was successfully created, False otherwise
        """
        cmd = 'CREATE TABLE IF NOT EXISTS {0} (cell_A TEXT, cell_B TEXT, histone TEXT, region TEXT, motif TEXT, ' \
              'n_genes INTEGER, n_motifs INTEGER, pearson REAL, ' \
              'p_smaller REAL, sig_smaller INTEGER, ' \
              'p_greater REAL, sig_greater INTEGER, ' \
              'p_other REAL, sig_other INTEGER, ' \
              'p_abs_smaller REAL, sig_abs_smaller INTEGER, ' \
              'p_abs_greater REAL, sig_abs_greater INTEGER, ' \
              'chi_p REAL, sig_chi INTEGER, ' \
              'corr_comment TEXT, ' \
              'Cramer_v REAL, corrected_v REAL, chi_comment TEXT, ' \
              'n_TFs INTEGER, n_miRNAs INTEGER, n_targets INTEGER, ' \
              'gene_list TEXT, motif_list TEXT, TF_list TEXT, miRNA_list TEXT, target_list TEXT, pairs TEXT, ' \
              'PRIMARY KEY (cell_A, cell_B, histone, region, motif))'

        table_name = self.p_value_table()

        return self.create_table(table_name, cmd, cursor=cursor, drop=True)

    def add_p_values(self, results: dict[str, dict[bool, list[tuple[str, str, str, str, str, int, float, float, float,
                                                                    bool, float, bool, float, bool, float, bool, float,
                                                                    bool, float, bool, float, bool, float, bool, float,
                                                                    bool, float, bool, float, bool, float, float, str,
                                                                    str]]]]):
        """
        Adds p-value results into the database for each edge type and differential condition.
        :param results: dictionary containing p-value results categorized by edge type and differential condition
        """
        for edge_type, content in results.items():
            for diff, table in content.items():
                connection = sqlite3.connect(self.p_value_path(edge_type, diff))
                cursor = connection.cursor()

                # improve performance by adjusting SQLite pragmas
                cursor.execute('PRAGMA synchronous = OFF')
                cursor.execute('PRAGMA journal_mode = MEMORY')
                cursor.execute('BEGIN TRANSACTION')

                self.create_p_value_table(cursor)
                cursor.executemany('INSERT INTO {0} (cell_A, cell_B, histone, region, motif, '
                                   'n_genes, n_motifs, n_TFs, n_miRNAs, n_targets, '
                                   'gene_list, motif_list, TF_list, miRNA_list, target_list, pearson, corr_comment, '
                                   'pairs, '
                                   'p_smaller, '
                                   'p_greater, '
                                   'p_other, '
                                   'p_abs_smaller, '
                                   'p_abs_greater, '
                                   'chi_p, Cramer_v, corrected_v, chi_comment, '
                                   'sig_smaller, '
                                   'sig_greater, '
                                   'sig_other, '
                                   'sig_abs_smaller, '
                                   'sig_abs_greater,'
                                   'sig_chi) '
                                   'VALUES ({1})'.format(self.p_value_table(), ', '.join(['?'] * 33)),
                                   table)
                connection.commit()
                connection.close()

    def p_value(self, edge_type: str, cell_1: str, cell_2: str, diff: bool) -> list[str, str, str, str]:
        """
        Fetches p-values for a specific pair of cells and their associated motif.
        :param edge_type: type of the edge (e.g. experimental)
        :param cell_1: the first cell type in the comparison
        :param cell_2: the second cell type in the comparison
        :param diff: whether to filter for differentially expressed genes
        :return: list of (histone, genomic region, motif, p-value) tuples
        """
        connection = sqlite3.connect(self.p_value_path(edge_type, diff))
        cursor = connection.cursor()

        cmd = 'SELECT histone, region, motif, p_greater ' \
              'FROM {0} WHERE cell_A="{1}" AND cell_B="{2}"'.format(self.p_value_table(), cell_1, cell_2)

        cursor.execute(cmd)
        results = cursor.fetchall()
        connection.close()

        # filter the results based on valid motifs
        return [res for res in results if res[2] in self.motif_names]

    def p_comparison_diff(self, edge_type: str) -> dict:
        """
        Compares p-values for differential conditions across multiple edge types.
        :param edge_type: type of the edge (e.g. experimental)
        """
        if len(self.differential) < 2:
            return []

        results_dicts = []
        for diff in [False, True]:
            connection = sqlite3.connect(self.p_value_path(edge_type, diff))
            cursor = connection.cursor()

            cmd = 'SELECT cell_A, cell_B, histone, region, motif, pearson, n_genes ' \
                  'FROM {0}'.format(self.p_value_table())

            cursor.execute(cmd)
            results = cursor.fetchall()

            # store results filtered by valid gene categories
            results_dicts.append({res[:5]: res[-2:] for res in results
                                  if res[-3] in ['all genes', 'other', 'all TF', 'all miRNA']})

            connection.close()

        non_diff, diff = results_dicts

        # merge the non-differential and differential results
        return {key: (*r, *diff[key]) for key, r in non_diff.items()}

    def p_comparison_all(self, edge_type: str, diff: bool) -> dict:
        """
        Fetches all p-value data for a specified edge type and differential condition.
        :param edge_type: type of the edge (e.g. experimental)
        :param diff: whether to filter for differentially expressed genes
        """
        connection = sqlite3.connect(self.p_value_path(edge_type, diff))
        cursor = connection.cursor()

        cmd = 'SELECT cell_A, cell_B, histone, region, motif, pearson, n_genes ' \
              'FROM {0}'.format(self.p_value_table())

        cursor.execute(cmd)
        results = cursor.fetchall()

        # filter results based on motif names
        results_dict = {res[:5]: res[-2:] for res in results
                        if res[-3] in self.motif_names}

        connection.close()

        return results_dict

    def p_comparison_all2(self, edge_type: str, diff: bool) -> dict:
        """
        Fetches all p-value data with additional information like gene list for a specified edge type and differential
        condition.
        :param edge_type: type of the edge (e.g. experimental)
        :param diff: whether to filter for differentially expressed genes
        """
        connection = sqlite3.connect(self.p_value_path(edge_type, diff))
        cursor = connection.cursor()

        cmd = 'SELECT cell_A, cell_B, histone, region, motif, gene_list, pearson, n_genes ' \
              'FROM {0}'.format(self.p_value_table())

        cursor.execute(cmd)
        results = cursor.fetchall()

        # filter results based on valid motifs
        results_dict = {res[:5]: res[-3:] for res in results
                        if res[-4] in self.motif_names}

        connection.close()

        return results_dict

    def p_comparison_region(self, edge_type: str, diff: bool) -> dict:
        """
        Compares p-values across regions for a specified edge type and differential condition.
        :param edge_type: type of the edge (e.g. experimental)
        :param diff: whether to filter for differentially expressed genes
        """
        connection = sqlite3.connect(self.p_value_path(edge_type, diff))
        cursor = connection.cursor()

        result_dict = defaultdict(list)

        # iterate over all regions and fetch the results for each
        for region in sorted(self.regions):

            cmd = 'SELECT cell_A, cell_B, histone, motif, n_genes, gene_list, pearson ' \
                  'FROM {0} WHERE region="{1}"'.format(self.p_value_table(), region)

            cursor.execute(cmd)
            results = cursor.fetchall()

            for res in results:
                result_dict[res[:-1]].append(res[-1])

        connection.close()

        result = dict()

        # process results for each key
        for key, v in result_dict.items():
            a, b, h, m, n, genes = key

            if not genes:
                continue
            # generate pairs for gene lists
            pairs = self.histone_pairs((a, b, diff), [h], self.regions, 'mod_diff', genes.split(';'))

            result[(a, b, h, m)] = (n, *v, pairs)

        return result

    def p_to_tsv(self):
        """ Exports the p-value data for all differential conditions to TSV files. """
        for diff in self.differential:
            connection = sqlite3.connect(self.p_value_path('experimental', diff))

            df = pd.read_sql_table(self.p_value_table(), con=connection)  # type: pd.DataFrame
            df.to_csv(os.path.join(self.result_dir, 'p_values_{0}.tsv'.format(diff)), sep='\t', index=False)

            connection.close()

    def all_genes_list(self) -> set[tuple[str, str, str, str]]:
        """
        Retrieves the list of genes marked as 'all' or 'other' in the experimental differential dataset.
        :return: set of tuples containing the gene-related information
        """
        connection = sqlite3.connect(self.p_value_path('experimental', True))
        cursor = connection.cursor()

        cmd = 'SELECT cell_A, cell_B, motif, gene_list FROM {0}'.format(self.p_value_table())
        cursor.execute(cmd)
        results = cursor.fetchall()

        return {res for res in results if 'all' in res[2] or res[2] == 'other'}

    def motif_genes(self) -> set[tuple[str, str, str, str, str, str, str, str]]:
        """
        Retrieves genes associated with specific motifs from the experimental dataset.
        :return: set of tuples containing the motif-associated genes
        """
        connection = sqlite3.connect(self.p_value_path('experimental', False))
        cursor = connection.cursor()

        cmd = 'SELECT cell_A, cell_B, motif, n_genes, gene_list, TF_list, miRNA_list, ' \
              'target_list FROM {0}'.format(self.p_value_table())
        cursor.execute(cmd)
        results = cursor.fetchall()

        return {res for res in results if res[2] in self.motif_names and res[2] != 'all'}
