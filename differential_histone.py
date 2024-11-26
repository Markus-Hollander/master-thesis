from functools import partial
import multiprocessing as mp
import os
from shutil import rmtree
from typing import Callable

from index import HistoneFile, HistoneQuantification, HistoneIndex
from transition import DifferentialHistone, Differential
import util


def bam_duplicates(file_path: str) -> int:
    """
    Computes the number of duplicates in the latest version of a .bam file.
    :param file_path: path to the stat summary of a .bam file
    :return: number of duplicates in the .bam file
    """
    entries = []                                                                    # type: list[list[str]]

    # get a list of all stat entries in the file
    with open(file_path, 'r') as in_file:
        for line in in_file:
            line = line.strip()                                                     # type: str
            # a new entry starts
            if line.endswith('in total (QC-passed reads + QC-failed reads)'):
                entries.append([line])
            # a previous entry is continued
            elif entries:
                entries[-1].append(line)

    if not entries:
        raise ValueError('The .bam stat summary does not contain any entries.')

    # extract the number of duplicates from the last entry in the stat summary
    for line in entries[-1]:
        if line.endswith('duplicates'):
            return int(line.split(' + ')[0])

    raise ValueError('The last entry in the .bam stats summary does not contain duplicate information.')


def bam_sort_command(input_path: str, sorted_path: str, log_path: str, position=False) -> bool:
    """
    Sorts the given .bam file with samtools.
    :param input_path: path for the input file
    :param sorted_path: path to the output, sorted file
    :param log_path: path to the log file
    :param position: True for position sort, False for name sort
    :return: True if the command terminated without errors, False otherwise
    """
    # add -n for name sort
    cmd_args = ['samtools', 'sort'] if position else ['samtools', 'sort', '-n']                 # type: list[str]
    cmd_args += ['-o', sorted_path, input_path]                                                 # type: list[str]

    return util.base_command(log_path, cmd_args)


def bam_mark_duplicates_command(input_path: str, marked_path: str, log_path: str) -> bool:
    """
    Marks duplicates in the given .bam file.
    :param input_path: path for the input file
    :param marked_path: path to the output, marked file
    :param log_path: path to the log file
    :return: True if the command terminated without errors, False otherwise
    """
    cmd_args = ['samtools', 'fixmate', '-m', input_path, marked_path]                           # type: list[str]

    return util.base_command(log_path, cmd_args)


def bam_remove_duplicates_command(sorted_marked_path: str, final_path: str, log_path: str) -> bool:
    """
    Remove duplicates in the given .bam file.
    :param sorted_marked_path: path for the input file
    :param final_path: path to the output file without duplicates
    :param log_path: path to the log file
    :return: True if the command terminated without errors, False otherwise
    """
    cmd_args = ['samtools', 'markdup', '-r', sorted_marked_path, final_path]                    # type: list[str]

    return util.base_command(log_path, cmd_args)


def bam_stat_command(input_path: str, stat_path: str, log_path: str, overwrite=False) -> bool:
    """
    Compute stats for the given .bam file.
    :param input_path: path for the input file
    :param stat_path: path to the stat path
    :param log_path: path to the log file
    :param overwrite: True to start with an empty log, False otherwise
    :return: True if the command terminated without errors, False otherwise
    """
    cmd_args = ['samtools', 'flagstat', input_path, '>', stat_path]                             # type: list[str]

    return util.base_command(log_path, cmd_args, overwrite_log=overwrite)


def bam_index_command(input_path: str, log_path: str) -> bool:
    """
    Generate the index for the given .bam file.
    :param input_path: path for the input file
    :param log_path: path to the log file
    :return: True if the command terminated without errors, False otherwise
    """
    cmd_args = ['samtools', 'index', input_path]                                                # type: list[str]

    return util.base_command(log_path, cmd_args, overwrite_log=False)


def bam_merge_command(run: HistoneQuantification) -> tuple[str, str, bool]:
    """
    Use samtools to merge the histone ChIP-seq signal .bam files belonging to one cell-histone combination.
    :param run: object representing the histone quantification for one cell-histone combination
    :return: cell type, histone name, success flag
    """
    cmd_args = ['samtools', 'merge', run.bam_path] + list(run.bam_files)                        # type: list[str]

    return run.cell, run.histone, util.base_command(run.log_path, cmd_args, overwrite_log=True)


def remove_bam_duplicates(run: HistoneQuantification | HistoneFile) -> tuple[str, str]:
    """
    Checks if the given ChIP-seq (merged or single) signal .bam file contains duplicates. If not, just position
    sort the entries for individual files, or build the index for merged files. If yes, use samtools to mark and
    remove duplicates and position sort the entries.
    :param run: object representing a single or merged histone ChIP-seq signal .bam file
    :return: (run name, pipeline step name) if the pipeline terminated without error, ('', step name) otherwise
    """
    # generate the index for merged .bam files
    if isinstance(run, HistoneQuantification):
        if not bam_index_command(run.bam_path, run.log_path):
            return '', 'index'

    # compute the stats of the input .bam file
    if not bam_stat_command(run.bam_path, run.stat_path, run.log_path, overwrite=True):
        return '', 'stat 1'

    # check if the .bam file contains duplicates
    if bam_duplicates(run.stat_path) == 0:
        # for merged files there is nothing to do here
        if isinstance(run, HistoneQuantification):
            return run.name, 'copied'
        # if not, position sort the entries and return
        if not bam_sort_command(run.bam_path, run.duplicate_free_path, run.log_path, position=True):
            return '', 'copied'
        return run.name, 'copied'

    # to remove the duplicates, first name sort the entries
    if not bam_sort_command(run.bam_path, run.name_sorted_path, run.log_path,  position=False):
        return '', 'remove'

    # then mark the duplicates
    if not bam_mark_duplicates_command(run.name_sorted_path, run.marked_path, run.log_path):
        return '', 'remove'

    # then position sort again
    if not bam_sort_command(run.marked_path, run.pos_sorted_path, run.log_path, position=True):
        return '', 'remove'

    # then remove the duplicates
    if not bam_remove_duplicates_command(run.pos_sorted_path, run.duplicate_free_path, run.log_path):
        return '', 'remove'

    # then compute the stats of the processed .bam file
    if not bam_stat_command(run.duplicate_free_path, run.stat_path, run.log_path):
        return '', 'stat 2'

    # return an error if there are still duplicates
    if bam_duplicates(run.stat_path) != 0:
        return '', 'duplicates'

    # pipeline succeeded without errors
    return run.name, 'remove'


def histone_remove_duplicates(cfg: dict, histone_index: HistoneIndex, merge=False):
    """
    Pre-process single or merged histone ChIP-seq signal .bam files by removing duplicates. For single files,
    sort them so that they are ready for merging.
    :param cfg: configuration dictionary
    :param histone_index: contains information on all histone ChIP-seq runs
    :param merge: True to process the merged .bam files, False to process the single .bam files
    """
    # shorthand for the directory dictionaries
    sub_dir = 'quantification' if merge else 'pre-processing'                       # type: str
    path_dir = cfg['file paths']['differential histone'][sub_dir]                   # type: dict[str, str]

    # remove old files
    for key, directory in path_dir.items():
        if key in ['logs', 'results', 'stats']:
            rmtree(directory, ignore_errors=True)

    # create output directories
    for directory in histone_index.get_dirs_to_create(merge):
        os.makedirs(directory, exist_ok=True)

    # get a list of all histone ChIP-seq signal .bam files
    if merge:
        runs = histone_index.get_quant(merge_successful=True)                       # type: list[HistoneQuantification]
        print('\tprocessing {0:,} merged files'.format(len(runs)), end='\r')
    else:
        runs = histone_index.get_runs()                                             # type: list[HistoneFile]
        print('\tpre-processing {0:,} histone files'.format(len(runs)), end='\r')

    # pre-process the files in parallel
    with mp.Pool(processes=min(len(runs), 32), maxtasksperchild=1) as pool:
        results = pool.map(remove_bam_duplicates, runs)                             # type: list[tuple[str, str]]

    # extract which files could be successfully pre-processed, which were already ready and for which it failed
    copied = {name for name, flag in results if flag == 'copied' and name}          # type: set
    successful = {name for name, flag in results if name and flag == 'remove'}      # type: set
    failed = {name for name in results if not name}                                 # type: set

    if merge:
        out_str = '\tprocessed {0:,} and copied {1:,} out of {2:,} merged files, failed: {3:,}'         # type: str
    else:
        out_str = '\tpre-processed {0:,} and copied {1:,} out of {2:,} histone files, failed: {3:,}'    # type: str
    print(out_str.format(len(successful), len(copied), len(runs), len(failed)))

    # set the success flag for the pre-processed files
    if merge:
        histone_index.add_merged_success(set(list(copied) + list(successful)))
    else:
        histone_index.set_duplicate_removal_success(set(list(copied) + list(successful)))
    print('\tupdated file flags')

    # remove files of for which processing failed
    for run in runs:
        if run.name in failed:
            util.remove_file(run.duplicate_free_path)
    print('\tremoved final files of failed processing')

    # remove intermediate pre-processing files
    rmtree(path_dir['intermediate'], ignore_errors=True)

    for sub_sub_dir in ['results', 'logs', 'intermediate', 'stats']:
        util.remove_empty_sub_dirs(cfg['file paths']['differential histone'][sub_dir][sub_sub_dir])
    print('\tremoved intermediate files and empty sub-directories')


def merge_histone_alignments(cfg: dict, histone_index: HistoneIndex):
    """
    Merge single ChIP-seq histone signal .bam files of each histone-cell combination.
    :param cfg: configuration dictionary
    :param histone_index: contains information on all histone ChIP-seq runs
    """
    # shorthand for file directories
    path_dir = cfg['file paths']['differential histone']['quantification']          # type: dict[str, str]
    sub_dirs = ['stats', 'logs', 'merged']                                          # type: list[str]

    # remove previous merged files, logs and stats
    for sub_dir in sub_dirs:
        rmtree(path_dir[sub_dir], ignore_errors=True)
    print('\tremoved previous results')

    # create the necessary output directories
    for directory in histone_index.get_dirs_to_create(merge=True):
        os.makedirs(directory, exist_ok=True)
    print('\tcreated directories')

    # get all cell-histone combinations with the corresponding individual histone .bam files
    runs = histone_index.get_quant()                                                # type: list[HistoneQuantification]
    print('\tgenerating {0:,} merged files'.format(len(runs)), end='\r')

    # merge the files in parallel
    with mp.Pool(processes=min(len(runs), 32), maxtasksperchild=1) as pool:
        results = pool.map(bam_merge_command, runs)                                 # type: list[tuple[str, str, bool]]

    # get the runs that succeeded and the ones that failed
    successful = {(cell, histone) for cell, histone, flag in results if flag}       # type: set[tuple[str, str]]
    failed = {(cell, histone) for cell, histone, flag in results if not flag}       # type: set[tuple[str, str]]
    print('\tgenerated {0:,}/{1:,} merged files, failed: {2:,}'.format(len(successful), len(runs), len(failed)))

    # set the success flag for merging
    histone_index.set_merged_success(successful)
    print('\tset success flags')

    # remove files of failed merging
    for run in runs:
        if (run.cell, run.histone) in failed:
            util.remove_file(run.duplicate_free_path)

    # remove empty sub-directories
    for sub_dir in sub_dirs:
        util.remove_empty_sub_dirs(path_dir[sub_dir])
    print('\tremoved failed merge files and empty sub-directories')


def histone_call_regions_command(cfg: dict, histone_merge: bool, run: HistoneQuantification) -> tuple[str, str, bool]:
    """
    Use histoneHMM to pre-process the histone modifications for one cell-histone combination.
    :param cfg: configuration dictionary
    :param histone_merge: True if the .bam files were merged previously => do not override the log file
    :param run: merged histone ChIP-seq information for one cell-histone combination
    :return: cell, histone, and boolean representing successful execution of the command
    """
    # file paths for the histoneHMM R script, chromosome lengths and gene annotations
    script_path = cfg['file paths']['differential histone']['quantification']['histoneHMM']     # type: str
    chrom_path = cfg['file paths']['reference']['chromosome lengths']                           # type: str
    # parameters: bin size in which to divide the genome and probability cut-off
    bin_size = str(cfg['differential histone analysis']['bin size'])                            # type: str
    p = str(cfg['differential histone analysis']['probability'])                                # type: str

    cmd_args = ['Rscript', script_path, '-c', chrom_path, '-b', bin_size, '-P', p,
                '-o', run.out_prefix, '--verbose', run.duplicate_free_path]                     # type: list[str]

    return run.cell, run.histone, util.base_command(run.log_path, cmd_args, not histone_merge)


def histone_call_regions(cfg: dict, histone_index: HistoneIndex, histone_merge: bool):
    """
    Quantify the histone modifications for each cell-histone combination.
    :param cfg: configuration dictionary
    :param histone_index: contains information on all histone ChIP-seq runs
    :param histone_merge: True if the .bam files were merged previously => do not override the log file
    """
    # get all cell-histone combinations with the corresponding individual, duplicate free histone .bam files
    runs = histone_index.get_quant(duplicate_free=True)                             # type: list[HistoneQuantification]

    # remove previous results
    for run in runs:
        rmtree(run.result_dir, ignore_errors=True)
    print('\tdeleted previous results')

    # create the required result directories
    for run in runs:
        os.makedirs(run.result_dir, exist_ok=True)
        os.makedirs(os.path.dirname(run.log_path), exist_ok=True)
    print('\tcreated result directories')

    print('\tquantifying {0:,} merged histone files'.format(len(runs)), end='\r')

    command = partial(histone_call_regions_command,
                      cfg, histone_merge)             # type: Callable[[HistoneQuantification], tuple[str, str, bool]]

    # analyse the runs in parallel
    with mp.Pool(processes=min(len(runs), 32), maxtasksperchild=1) as pool:
        results = pool.map(command, runs)

    # get the runs that succeeded and the ones that failed
    successful = {(cell, histone) for cell, histone, flag in results if flag}       # type: set[tuple[str, str]]
    failed = {(cell, histone) for cell, histone, flag in results if not flag}       # type: set[tuple[str, str]]
    print('\tquantified {0:,}/{1:,} merged files, failed: {2:,}'.format(len(successful), len(runs), len(failed)))

    # set the success flag for quantification
    histone_index.set_quant_success(successful)
    print('\tset success flags')

    # remove files of failed quantification
    for run in runs:
        if (run.cell, run.histone) in failed:
            rmtree(run.result_dir, ignore_errors=True)
    print('\tremoved failed histone quantification files')


def histone_differential_command(cfg: dict, comp: DifferentialHistone) -> tuple[tuple[str, str, str], bool]:
    """
    Use histoneHMM to compute differential histone modification for a differentiation transition.
    :param cfg: dictionary with file paths and settings
    :param comp: quantified histone modifications of two cells
    :return: (cell 1, cell 2, histone), and boolean representing successful execution of the command
    """
    # file path for the histoneHMM R script
    script_path = cfg['file paths']['differential histone']['differential']['histoneHMM']           # type: str
    # parameter: probability cut-off
    p = str(cfg['differential histone analysis']['probability'])                                    # type: str

    cmd_args = ['Rscript', script_path, '-o', comp.result_dir, '-P', p, '--sample1', comp.cell_1,
                '--sample2', comp.cell_2, '--verbose', comp.quant_1_path, comp.quant_2_path]        # type: list[str]

    return (comp.cell_1, comp.cell_2, comp.histone), util.base_command(comp.log_path, cmd_args, True)


def histone_differential(cfg: dict, histone_index: HistoneIndex, transitions: dict[tuple[str, str], Differential]):
    """
    Use histoneHMM to analyse differential histone modifications of cell differentiation transitions.
    :param cfg: configuration dictionary
    :param histone_index: contains information on all histone ChIP-seq runs
    :param transitions: dict of Differential object that contains file paths etc. for this differentiation transition
    """
    # shorthand for the directory dictionary
    path_dir = cfg['file paths']['differential histone']['differential']            # type: dict[str, str]

    # remove previous results and logs
    rmtree(path_dir['results'], ignore_errors=True)
    rmtree(path_dir['logs'], ignore_errors=True)
    print('\tremoved prior results')

    # load_results the histone quantification for each cell differentiation transition
    for transition in transitions.values():
        transition.load_histone_comparisons(histone_index, cfg)
    print('\tloaded histone comparisons')

    # for each differentiation transition get the histone comparisons
    runs = [c for t in transitions.values() for c in t.histone.values()]            # type: list[DifferentialHistone]

    # create the required result directories
    for run in runs:
        os.makedirs(run.result_dir, exist_ok=True)
        os.makedirs(os.path.dirname(run.log_path), exist_ok=True)
    print('\tcreated result directories')

    # pre-load_results the differential analysis command with the configuration dictionary
    func = partial(histone_differential_command, cfg)
    print('\tanalysing {0:,} histone comparisons'.format(len(runs)), end='\r')

    # analyse the differential modification in parallel
    with mp.Pool(processes=min(len(runs), 32), maxtasksperchild=1) as pool:
        results = pool.map(func, runs)              # list[tuple[tuple[str, str, str], bool]

    # get the runs that succeeded and the ones that failed
    success = {name for name, flag in results if flag}                              # set[tuple[str, str, str]]
    failed = {name for name, flag in results if not flag}                           # set[tuple[str, str, str]]
    print('\tanalysed {0:,}/{1:,} histone comparisons, failed: {2:,}'.format(len(success), len(runs), len(failed)))

    # set the success and failure flag
    for cell_1, cell_2, histone in success:
        transitions[(cell_1, cell_2)].histone[histone].success = True

    for cell_1, cell_2, histone in failed:
        transitions[(cell_1, cell_2)].histone[histone].success = False
    print('\tset success flags')

    # remove files of failed analyses
    for run in runs:
        if (run.cell_1, run.cell_2, run.histone) in failed:
            rmtree(run.result_dir, ignore_errors=True)
    print('\tremoved failed histone comparison results')
