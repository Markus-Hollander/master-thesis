from functools import partial
import multiprocessing as mp
import os
import shlex
from shutil import rmtree
import subprocess as sp
from typing import Callable


from index import RnaQuantRun, RnaIndex
from transition import Differential
import util


def salmon_index_command(cfg: dict):
    """
    :param cfg: configuration dictionary
    Construct the transcriptome index for RNA-seq quantification with salmon.
    """
    # shortcut of the relevant configuration section
    paths = cfg['file paths']['rna seq quantification']                                             # type: dict
    # components of the command
    cmd_args = [paths['salmon'], 'index', '-t', cfg['file paths']['reference']['transcriptome fasta'],
                '-i', paths['salmon index']]
    # make sure all command components are safe for execution on the shell
    cmd_args = [shlex.quote(arg) for arg in cmd_args]                                               # type: list[str]

    # execute the command in the shell
    print('\tconstructing the transcriptome index', end='\r')
    sp.check_call(args=' '.join(cmd_args), stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
    print('\tconstructed the transcriptome index ')


def salmon_quantification_command(cfg: dict, run: RnaQuantRun) -> bool:
    """
    Constructs and executes the salmon RNA-seq quantification command for the given arguments.
    :param cfg: configuration dictionary
    :param run: RNA index row with information on result directory, log file path and run type
    :return: (True, run ID) if quantification succeeded, (False, run result directory) if it failed
    """
    # shortcuts of the relevant configuration section
    quant = cfg['rna seq quantification']
    paths = cfg['file paths']['rna seq quantification']                                         # type: dict[str, str]

    # assemble the command
    cmd_args = [paths['salmon'], 'quant', '-i', paths['salmon index'], '-l', 'A', *run.salmon_cmd_files(),
                '-o', run.out_dir, '-p', quant['number cores'], '--gcBias', '--seqBias', '--incompatPrior',
                quant['incompatible prior'], '--numBootstraps', quant['number bootstraps'],
                '--validateMapping']

    # make sure all command components are safe for execution on the shell
    cmd_args = [shlex.quote(str(arg)) for arg in cmd_args]                                      # type: list[str]

    return util.base_command(run.log_path, cmd_args, overwrite_log=True)


def salmon_rna_seq_quantification(cfg: dict, rna_index: RnaIndex, refresh: bool):
    """
    Quantify RNA-seq quantification with salmon.
    :param cfg: configuration dictionary
    :param rna_index: contains information on all RNA-seq files and quantification runs
    :param refresh: True if previous results and logs should be deleted, False otherwise
    """
    # shortcut of the relevant configuration section
    paths = cfg['file paths']['rna seq quantification']                     # type: dict[str, str]

    if refresh:
        # delete previous results and logs if specified
        rmtree(paths['results'], ignore_errors=True)
        rmtree(paths['logs'], ignore_errors=True)
        print('\tdeleted previous results and logs')

    # create the result and log (sub-) directories
    for directory in rna_index.dirs_to_create:
        os.makedirs(directory, exist_ok=True)
    print('\tcreated result and log (sub-)directories')
    print()

    # obtain a list of quantification run name, result directory, log file path and (paired or single) file list
    runs = rna_index.get_not_quantified_runs(refresh)                       # type: list[RnaQuantRun]
    # pre-load_results the quantification command with the configuration dictionary
    command = partial(salmon_quantification_command, cfg)                   # type: Callable[[RnaQuantRun], bool]

    print('\tquantifying {0} runs with {1} files'.format(rna_index.number_of_runs, rna_index.n_files), end='\r')

    # run the quantification runs with two processes in parallel
    with mp.Pool(processes=4, maxtasksperchild=1) as pool:
        results = pool.map(command, runs)                                   # type: list[bool]

    successful = [run for run, flag in zip(runs, results) if flag]          # type: list[RnaQuantRun]
    failed = [run for run, flag in zip(runs, results) if not flag]          # type: list[RnaQuantRun]

    # set the success flag
    for run in successful:
        run.success = True

    n_failed = len(failed)                                                  # type: int

    # extract the IDs of successful quantification runs, and the result directories for failed runs
    print('\tquantified: {0:5,}/{1:,}, failed files: {2}'.format(len(successful), rna_index.number_of_runs, n_failed))

    # remove results of failed quantification runs
    for run in failed:
        rmtree(run.out_dir)
    print('\tremoved results of {0} failed quantification runs'.format(n_failed))

    # remove sub directories that were created during the download steps but do not contain files
    util.remove_empty_sub_dirs(paths['results'])
    util.remove_empty_sub_dirs(paths['logs'])
    print('\tremoved empty data subdirectories')


def differential_expression_command(cfg: dict, mapping: str, transition: Differential) -> bool:
    """
    Run the differential expression analysis R-script for this transition.
    :param cfg: configuration dictionary
    :param mapping: path to the transcript -> gene mapping
    :param transition: object that contains information and file paths for a cell differentiation transition
    :return: True if the R-script terminated without error, False otherwise
    """
    cmd_args = ['Rscript', cfg['file paths']['differential expression analysis']['script'],
                shlex.quote(mapping), shlex.quote(transition.dea_design_path),
                shlex.quote(transition.dea_result_path), transition.cell_1, transition.cell_2]      # type: list[str]

    # execute the differential expression command
    if not util.base_command(transition.dea_log_path, cmd_args, overwrite_log=True):
        return False

    # because R and Excel are incredibly stupid, the header of the result table needs ot be set manually
    # and the ID columns cannot be called ID because fuck you Excel
    with open(transition.dea_result_path, 'r') as file:
        lines = [line.replace('"', '') for line in file]                                            # type: list[str]
        lines[0] = '\t'.join(['name', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj', 'weight']) + '\n'

    # write the adjusted expression file
    with open(transition.dea_result_path, 'w') as file:
        file.writelines(lines)

    return True


def transcript_mapping(fasta_path: str, mapping_path: str) -> tuple[int, int]:
    """
    Generates the transcript ID to gene ID mapping.
    :param fasta_path: path of the reference transcriptome.fasta file
    :param mapping_path: path to the mapping file
    :return: number of unique transcript IDs, number of unique gene IDs
    """
    mapping = dict()                                                                # type: dict[str, str]

    with open(fasta_path, 'r') as in_file, open(mapping_path, 'w') as out_file:
        for line in in_file:
            # skip sequence lines
            if not line.startswith('>'):
                continue

            # extract the gene ID (and remove the trailing .X number) and the transcript ID
            gene_id = line.split(' gene:')[1].split()[0].split('.')[0]              # type: str
            transcript_id = line.strip('>').split()[0]                              # type: str
            # store the mapping and write it to the output file
            mapping[transcript_id] = gene_id
            out_file.write('{0}\t{1}\n'.format(transcript_id, gene_id))

    return len(mapping.keys()), len(set(mapping.values()))


def differential_expression_analysis(cfg: dict, rna_index: RnaIndex, transitions: list[Differential], mapping: bool):
    """
    Differential gene expression analysis with DESeq2.
    :param cfg: configuration dictionary
    :param rna_index: contains information on all RNA-seq files and quantification runs
    :param transitions: list of objects that contain information such as file paths for each cell transition
    :param mapping: True if the transcript ID to gene ID mapping should be recomputed if it already exists
    """
    # shorthand
    paths = cfg['file paths']['differential expression analysis']               # type: dict[str, str]

    # remove prior results
    for directory in [paths['results'], paths['design tables'], paths['logs']]:
        rmtree(directory, ignore_errors=True)
        os.makedirs(directory, exist_ok=True)
    print('\tdeleted previous results and setup the directories')

    # generate the differential expression analysis design tables
    rna_index.generate_design_tables(transitions)
    print('\tgenerated design tables')

    # generate the transcript ID to gene ID mapping if it does not exist or should be recomputed
    if mapping or not os.path.isfile(paths['transcript mapping']):
        n, m = transcript_mapping(cfg['file paths']['reference']['transcriptome fasta'], paths['transcript mapping'])
        print('\tmapped {0:,} transcripts to {1:,} genes'.format(n, m))

    # pre-load_results the command with the config and transcript mapping path
    command = partial(differential_expression_command,
                      cfg, paths['transcript mapping'])                         # type: Callable[[Differential], bool]

    # perform the differential expression analysis in parallel
    with mp.Pool(processes=min(len(transitions), 32)) as pool:
        results = list(pool.map(command, transitions))                          # type: list[bool]
    print('\tanalysed: {0:,}/{1:,}, failed: {2:,}'.format(results.count(True), len(results), results.count(False)))

    # set the differential expression analysis success flag
    for transition, success in zip(transitions, results):
        if not success:
            util.remove_file(transition.dea_result_path)
        transition.dea_success = success
    print('\tset success flags and removed failed results')



