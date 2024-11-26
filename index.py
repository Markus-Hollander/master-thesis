from collections import defaultdict
import os
import pandas as pd
import util


class RnaQuantRun:
    def __init__(self, cell: str, name: str, run_type: str, dirs: dict[str, str]):
        """
        :param cell: name of the cell for which the RNA-seq was conducted
        :param name: the name/identifier of the quantification run
        :param run_type: paired or single
        :param dirs: dictionary with the log and result directories
        """
        self.cell = cell                                                                    # type: str
        self.name = name                                                                    # type: str
        self.run_type = run_type                                                            # type: str
        # setup the output file and directory paths
        self.log_path = os.path.join(dirs['logs'], cell, name + '.txt')                     # type: str
        self.out_dir = os.path.join(dirs['results'], cell, name)                            # type: str
        self.quant_file = os.path.join(self.out_dir, 'quant.sf')                            # type: str
        # True if the quantification file exists, False otherwise
        self.success = os.path.isfile(self.quant_file)                                      # type: bool
        # number of files
        self.n = 0                                                                          # type: int

    def salmon_cmd_files(self) ->list[str]:
        """
        Returns the file component of the salmon quantification command for this run.
        """
        pass


class SingleRnaQuantRun(RnaQuantRun):
    def __init__(self, cell: str, name: str, run_type: str, dirs: dict[str, str], df: pd.DataFrame):
        """
        :param cell: name of the cell for which the RNA-seq was conducted
        :param name: the name/identifier of the quantification run
        :param run_type: paired or single
        :param dirs: dictionary with the log and result directories
        """
        super().__init__(cell, name, run_type, dirs)
        # list of paths to the RNA-seq read.fastq file belonging to this run
        self.files = list(sorted(df['file path'].unique()))                 # type: list[str]
        # number of RNA-seq read.fastq file belonging to this run
        self.n = len(self.files)                                            # type: int

    def salmon_cmd_files(self) -> list[str]:
        """
        Returns the file component of the salmon quantification command for this run:
        -r file_path_1 file_path_2 ...
        """
        return ['-r'] + self.files


class PairedRnaQuantRun(RnaQuantRun):
    def __init__(self, cell: str, name: str, run_type: str, dirs: dict[str, str], df: pd.DataFrame):
        """
        :param cell: name of the cell for which the RNA-seq was conducted
        :param name: the name/identifier of the quantification run
        :param run_type: paired or single
        :param dirs: dictionary with the log and result directories
        """
        super().__init__(cell, name, run_type, dirs)
        # file path index for all RNA-seq reads.fastq files belonging to this run
        path_index = {file_id: file_path
                      for file_id, file_path in {tuple(x) for x in df[['file accession', 'file path']].values}}

        # group the file paths by the sorted pair of file accession and specified paired file accession
        pair_dict = defaultdict(set)                                            # type: dict[tuple[str, str], set[str]]
        for file_id, file_path, paired_with in df[['file accession', 'file path', 'paired with']].values:
            pair_dict[tuple(sorted([file_id, paired_with]))].add(file_path)

        # set of pairs with exactly two files and not empty IDs
        complete_pairs = {pair for pair, files in pair_dict.items()
                          if len(files) == 2 and '' not in pair}                # type: set[tuple[str, str]]

        # list of paired file paths
        self.file_pairs = [(path_index[id_1], path_index[id_2])
                           for id_1, id_2 in sorted(complete_pairs)]            # type: list[tuple[str, str]]
        # number of file pairs
        self.n = len(self.file_pairs)                                           # type: int
        # number of pairs with less (or more) than two files
        self.incomplete_pairs = len(pair_dict) - len(complete_pairs)            # type: int

    def salmon_cmd_files(self) -> list[str]:
        """
        Returns the file component of the salmon quantification command for this run:
        -1 file_path_1a file_path_2a ... -2 file_path_1b file_path_2b ...
        """
        left, right = zip(*self.file_pairs)
        return ['-1'] + list(left) + ['-2'] + list(right)


class Index:
    assays = []
    output_types = []
    file_formats = []
    column = None

    def __init__(self, df: pd.DataFrame, dirs: dict[str, str]):
        """
        :param df: data frame containing the (filtered) ENCODE metadata
        :param dirs: dictionary with the relevant file directories
        """
        # extract the correctly downloaded data from the metadata frame
        df = df.loc[df['assay'].isin(self.assays)
                    & df['output type'].isin(self.output_types)
                    & df['file format'].isin(self.file_formats), :].copy()          # type: pd.DataFrame

        # add a column that contains the name of a run
        df['run name'] = df.apply(lambda row: self._get_name(row), axis=1)

        # runs grouped by the cell type
        self.runs = None                                                            # type: dict
        # set of result and log directories that need to be created
        self.dirs_to_create = set()                                                 # type: set[str]
        self.number_of_runs = 0                                                     # type: int
        self.runs_without_files = 0                                                 # type: int

        # generate a list of unique combinations of run name, cell/tissue and specified column
        runs = {tuple(x) for x in df[['run name', 'cell', self.column]].values}     # type: set[tuple[str, str, str]]

        for name, cell, x in sorted(runs):
            # extract all files belonging to this run
            run_df = df[(df['run name'] == name)                                    # type: pd.DataFrame
                        & (df[self.column] == x)
                        & (df['cell'] == cell)]

            # create the run objects
            self._add(cell, name, x, dirs, run_df)

    def _get_name(self, row: pd.Series) -> str:
        """
        Generates the run name of a file entry, which can be used to group files of an experiment run together.
        :param row: file entry of a run
        :return: run name
        """
        raise NotImplementedError

    def _add(self, cell, name, x, dirs, run_df):
        """
        Creates and adds objects representing experiment runs to the index. Should initialise the run dictionary if it
        does not exist yet.
        :param cell: cell or tissue name
        :param name: run name
        :param x: further information, e.g. run type or histone
        :param dirs: dictionary with the relevant directories
        :param run_df: data frame with file entries belonging to the run
        """
        raise NotImplementedError


class RnaIndex(Index):
    assays = ['RNA.seq']
    output_types = ['reads']
    file_formats = ['fastq']
    column = 'run type'

    def __init__(self, df: pd.DataFrame, dirs: dict[str, str], merge: bool, name_dict: dict[str, list[str]]):
        """
        :param df: data frame containing the (filtered) ENCODE metadata
        :param dirs: dictionary with the relevant file directories
        :param merge: True if technical replicates of a biological replicates should be quantified together
        :param name_dict: dictionary that contains the file name components
        """
        # number of incomplete paired-end file pairs
        self.incomplete_pairs = 0       # type: int
        self.name_dict = name_dict      # type: dict[str, list[str]]
        # True if technical replicates of a biological replicate should be merged
        self.merge = merge              # type: bool
        self.n_files = 0
        super().__init__(df, dirs)

    def _get_name(self, row: pd.Series):
        """
        Generates the run name of a file entry, which can be used to group files of an experiment run together.
        :param row: a row in the data frame
        :return: the file/directory name
        """
        if self.merge:
            name_components = [row[col] for col in self.name_dict['merge']]
        else:
            name_components = [row[col] for col in self.name_dict['no merge']]

        return '_'.join(name_components)

    def _add(self, cell, name, run_type, dirs: dict[str, str], run_df):
        """
        Creates and adds RNA quantification runs to the index. Should initialise the run dictionary if it
        does not exist yet.
        :param cell: cell or tissue name
        :param name: run name
        :param run_type: 'single' or 'paired'
        :param dirs: dictionary with the relevant directories
        :param run_df: data frame with file entries belonging to the run
        """
        # create the RnaQuantRun objects for single-end and paired-end quantification runs
        if run_type == 'single':
            run = SingleRnaQuantRun(cell, name, run_type, dirs, run_df)
        else:
            run = PairedRnaQuantRun(cell, name, run_type, dirs, run_df)
            # log the number of pairs that did not have the correct number of files
            self.incomplete_pairs += run.incomplete_pairs

        # initialise the runs dictionary if it does not exist yet
        if not self.runs:
            self.runs = defaultdict(list)

        # add the RNA quantification run to the index if it has valid files
        if run.n:
            self.runs[run.cell].append(run)
            self.number_of_runs += 1

            if isinstance(run, PairedRnaQuantRun):
                self.n_files += run.n * 2
            else:
                self.n_files += run.n
        else:
            self.runs_without_files += 1

        self.dirs_to_create.add(os.path.dirname(run.log_path))
        self.dirs_to_create.add(run.out_dir)

    def get_not_quantified_runs(self, refresh: bool) -> list[RnaQuantRun]:
        """
        :param refresh: True if all runs should be returned
        :return: list of RNA-seq quantification runs that have not been completed yet
        """
        if refresh:
            return [run for cell_list in self.runs.values() for run in cell_list]
        return [run for cell_list in self.runs.values() for run in cell_list if not run.success]

    def get_completed_quant_files(self) -> list[str]:
        """
        :return: list of RNA-seq quantification runs that have been completed
        """
        return [run.quant_file for cell_list in self.runs.values() for run in cell_list if run.success]

    def generate_design_tables(self, transitions: list):
        """
        For each differentiation transition generate the differential expression analysis design table.
        :param transitions: list of cell differentiation transitions objects that contain file paths etc.
        """
        for transition in transitions:
            with open(transition.dea_design_path, 'w') as file:
                # header
                file.write('\t'.join(['name', 'condition', 'path']) + '\n')
                # for each run belonging to the two cell types, write the relevant information
                for run in self.runs[transition.cell_1] + self.runs[transition.cell_2]:
                    file.write('\t'.join([run.name, run.cell, run.quant_file]) + '\n')


class HistoneFile:
    """
    Represents a single histone ChIP-seq .bam file.
    """
    def __init__(self, cell, name, histone, dirs: dict[str, str], bam_path: str):
        """
        :param cell: name of the cell for which the RNA-seq was conducted
        :param name: the name/identifier of the quantification run
        :param histone: histone name
        :param dirs: dictionary with output directory paths
        :param bam_path: path to a histone ChIP-seq .bam file
        """
        self.cell = cell                                                                                # type: str
        self.name = name                                                                                # type: str
        self.histone = histone                                                                          # type: str

        file_name = os.path.basename(bam_path)                                                          # type: str

        # path to the processing log of this file
        self.log_path = os.path.join(dirs['logs'], cell, histone, name + '.txt')                        # type: str
        # samtools stats for this .bam file
        self.stat_path = os.path.join(dirs['stats'], cell, histone, name + '.summary')                  # type: str
        # directory for the intermediate versions of the .bam file
        self.intermediate_dir = os.path.join(dirs['intermediate'], cell, histone)                       # type: str

        # 0. input .bam file
        self.bam_path = bam_path                                                                        # type: str
        # 1. path to the name sorted file
        self.name_sorted_path = os.path.join(self.intermediate_dir, 'name_sorted_' + file_name)         # type: str
        # 2. path to the name sorted file with marked duplicates
        self.marked_path = os.path.join(self.intermediate_dir, 'marked_' + file_name)                   # type: str
        # 3. path to the position sorted file with marked duplicates
        self.pos_sorted_path = os.path.join(self.intermediate_dir, 'pos_sorted_' + file_name)           # type: str
        # 4. final duplicate-free and position-sorted version of the input .bam file
        self.duplicate_free_path = os.path.join(dirs['results'], cell, histone, file_name)              # type: str

        # True if the quantification file exists, False otherwise
        self.duplicate_success = os.path.isfile(self.duplicate_free_path)                               # type: bool


class HistoneQuantification:
    """
    Merged file of a histone quantification run.
    """
    def __init__(self, cell, histone, runs: list[HistoneFile], dirs: dict[str, str]):
        """
        :param cell: name of the cell for which the RNA-seq was conducted
        :param histone: histone name
        :param runs: list of objects containing paths to a histone ChIP-seq .bam files
        :param dirs: dictionary with output directory paths
        """
        self.cell = cell                                                                            # type: str
        self.histone = histone                                                                      # type: str
        self.name = (cell, histone)                                                                 # tuple[str, str]
        # extract the paths to the duplicate-free, sorted versions
        self.bam_files = sorted({run.duplicate_free_path for run in runs})                          # list[str]

        # path to the processing log of this file
        self.log_path = os.path.join(dirs['logs'], cell, histone + '.txt')                          # type: str
        # samtools stats for this .bam file
        self.stat_path = os.path.join(dirs['stats'], cell, histone + '.summary')                    # type: str
        # directory for the intermediate versions of the .bam file
        self.intermediate_dir = os.path.join(dirs['intermediate'], cell)                            # type: str

        # 0. input .bam file
        self.bam_path = os.path.join(dirs['merged'], cell, histone + '_merged.bam')                 # type: str
        file_name = os.path.basename(self.bam_path)                                                 # type: str
        # 1. path to the name sorted file
        self.name_sorted_path = os.path.join(self.intermediate_dir, 'name_sorted_' + file_name)     # type: str
        # 2. path to the name sorted file with marked duplicates
        self.marked_path = os.path.join(self.intermediate_dir, 'marked_' + file_name)               # type: str
        # 3. path to the position sorted file with marked duplicates
        self.pos_sorted_path = os.path.join(self.intermediate_dir, 'pos_sorted_' + file_name)       # type: str
        # 4. final duplicate-free and position-sorted version of the input .bam file
        self.duplicate_free_path = self.bam_path                                                    # type: str

        # directory for the quantification results
        self.result_dir = os.path.join(dirs['results'], cell, histone)                              # type: str
        # out prefix needed for histoneHMM call regions
        self.out_prefix = os.path.join(self.result_dir, '_'.join([histone, cell]))                  # type: str
        # path to the histoneHMM call regions output file
        self.quant_path = os.path.join(self.out_prefix + '.txt')                                    # type: str

        # success flags for merging, duplicate removal and quantification
        self.merge_success = os.path.isfile(self.bam_path)                                          # type: bool
        self.duplicate_success = os.path.isfile(self.duplicate_free_path)                           # type: bool
        self.success = os.path.isfile(self.quant_path)                                              # type: bool


class HistoneIndex(Index):
    assays = ['ChIP.seq']
    output_types = ['alignments']
    file_formats = ['bam']
    column = 'experiment target'

    def __init__(self, df: pd.DataFrame, dirs: dict[str, str]):
        """
        Has a histone file index, as well as a histone quantification index.
        :param df: data frame containing the (filtered) ENCODE metadata
        :param dirs: dictionary with the relevant file directories
        """
        super().__init__(df, dirs)
        # histone quantification index
        self.quantification = dict()                                # type: dict[str, dict[str, HistoneQuantification]]

    def _get_name(self, row: pd.Series):
        """
        Generates the run name of a file entry, which can be used to group files of an experiment run together.
        :param row: a row in the data frame
        :return: the file/directory name
        """
        return '.'.join(os.path.basename(row['file path']).split('.')[:-1])

    def _add(self, cell, name, histone, dirs, run_df):
        """
        Creates and adds objects representing ChIP-seq histone experiment runs to the index. Should initialise the run
        dictionary if it does not exist yet.
        :param cell: cell or tissue name
        :param name: run name
        :param histone: histone target of the experiment
        :param dirs: dictionary with the relevant directories
        :param run_df: data frame with file entries belonging to the run
        """
        # initialise the run index if it does not exist yet
        if not self.runs:
            self.runs = defaultdict(util.ddl)                           # type: dict[str, dict[str, list[HistoneFile]]]

        # create an object for each file belonging to this experiment and add it to the index
        for run in [HistoneFile(cell, name, histone, dirs, file_path) for file_path in run_df['file path'].unique()]:
            self.runs[run.cell][run.histone].append(run)
            self.number_of_runs += 1

    def get_runs(self) -> list[HistoneFile]:
        """"
        :return: list of all ChIP-seq histone files in the index
        """
        return [run for run_dict in self.runs.values() for run_list in run_dict.values() for run in run_list]

    def get_quant(self, merge_successful=False, duplicate_free=False, final=False) -> list[HistoneQuantification]:
        """
        :param merge_successful: only those that could be successfully merged
        :param duplicate_free: only those for whom duplicate entries could successfully removed
        :param final: only those for whom the entire quantification pipeline was successful
        :return: list of the specified ChIP-seq histone quantification summaries in the index
        """
        # get a list of all quantification objects
        runs = [run for histone_dict in self.quantification.values() for run in histone_dict.values()]

        # return the specified subset of quantification objects
        if final:
            return [run for run in runs if run.success]
        if duplicate_free:
            return [run for run in runs if run.duplicate_success]
        if merge_successful:
            return [run for run in runs if run.merge_success]
        return runs

    def set_duplicate_removal_success(self, names: set[tuple[str, str]] | set[str], quant=False):
        """
        sets the duplicate removal success flag of the given histone quantification runs or single histone files
        to True and that of those not given to False.
        :param names: (cell, histone) of a histone quantification runs, or histone file IDs
        :param quant: True if the names belong to histone quantification runs, False for histone file IDs
        """
        runs = self.get_quant() if quant else self.get_runs()

        for run in runs:
            run.duplicate_success = run.name in names

    def add_merged_success(self, ids: set[tuple[str, str]]):
        """
        sets the merge success flag of the given histone quantification runs to True without touching the success flag
        of the others.
        :param ids: (cell, histone) of histone quantification runs
        """
        for run in self.get_quant():
            if run.name in ids:
                run.merge_success = True

    def set_merged_success(self, ids: set[tuple[str, str]]):
        """
        sets the merge success flag of the given histone quantification runs to True and that of those not given
        to False.
        :param ids: (cell, histone) of histone quantification runs
        """
        for run in self.get_quant():
            run.merge_success = run.name in ids

    def set_quant_success(self, ids: set[tuple[str, str]]):
        """
        sets the quantification success flag of the given histone quantification runs to True and that of those not
        given to False.
        :param ids: (cell, histone) of histone quantification runs
        """
        for run in self.get_quant():
            run.success = run.name in ids

    def load_quantification(self, dirs: dict[str, str]):
        """
        :param dirs: dictionary with file paths for histone quantification
        """
        # initialise the quantification dictionary
        self.quantification = defaultdict(dict)                     # type: dict[str, dict[str, HistoneQuantification]]

        # group the individual histone files by cell and histone, and then combine them
        for cell, histone_dict in self.runs.items():
            for histone, run_list in histone_dict.items():
                runs = [run for run in run_list if run.duplicate_success]                   # type: list[HistoneFile]

                self.quantification[cell][histone] = HistoneQuantification(cell, histone, runs, dirs)

        # convert to a normal dictionary
        self.quantification = dict(self.quantification)

    def get_dirs_to_create(self, merge=False) -> set[str]:
        """
        :param merge: True for getting the directories of merged files, False for single histone files
        :return: set of directory paths
        """
        runs = self.get_quant() if merge else self.get_runs()

        dirs_to_create = set()

        for run in runs:
            dirs_to_create.add(run.intermediate_dir)
            dirs_to_create.add(os.path.dirname(run.duplicate_free_path))
            dirs_to_create.add(os.path.dirname(run.log_path))
            dirs_to_create.add(os.path.dirname(run.stat_path))

        return dirs_to_create
