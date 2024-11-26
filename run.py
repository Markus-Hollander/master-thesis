import argparse as ap
import os
from time import time
import util
import warnings

import analysis
import differential_expression as dea
import differential_histone as dha
import download
from index import RnaIndex, HistoneIndex
from manager import Manager
from meta import Metadata
import processing
from reference import NewReference, Reference
from transition import Differential

# settings
warnings.simplefilter('ignore', UserWarning)
# noinspection PyBroadException


class Pipeline:
    def __init__(self, args):
        # CONFIGURATION
        # contains the path to the configuration file, as well as run flags
        self.args = args
        # load_results the configuration file
        self.cfg = util.load_yaml(args.config)              # type: dict

        # INITIALISATION
        self.meta = None                                    # type: Metadata | None
        self.transitions = dict()                           # type: dict[tuple[str, str], Differential]
        self.rna_index = None                               # type: RnaIndex | None
        self.histone_index = None                           # type: HistoneIndex | None
        self.reference = None                               # type: Reference | None

        # PIPELINES
        self._load_metadata()
        self._download_pipeline()
        self.meta.filter_download()
        self._initialise_indices()
        self._rna_seq_quantification_pipeline()
        self._differential_expression_pipeline()
        self._histone_pipeline()
        self._analysis_pipeline()

        print()

    def _load_metadata(self):
        """
        Reads in the specified ENCODE metadata.
        """
        print()
        print('METADATA')
        start = time()

        # if specified, read in the entire metadata and then filter it down and output the filtered metadata
        if self.args.meta_filter:
            self.meta = Metadata(in_path=self.cfg['file paths']['metadata']['all tsv'],
                                 download_config=self.cfg['download'],
                                 download_dir=self.cfg['file paths']['data directory'],
                                 name_components=self.cfg['file names']['download name components'],
                                 filter_dict=self.cfg['metadata filter'],
                                 adjustment_dict=self.cfg['metadata adjustment'])
            print('\tloaded and filtered entire metadata')
        # otherwise just read in the already filtered metadata
        else:
            self.meta = Metadata(in_path=self.cfg['file paths']['metadata']['filtered tsv'],
                                 download_config=self.cfg['download'])
            print('\tloaded already filtered metadata')

        # update the metadata .tsv and .db, in case the input file did not contain download information
        self.meta.db_path = self.cfg['file paths']['metadata']['filtered db']
        self.meta.tsv_path = self.cfg['file paths']['metadata']['filtered tsv']
        self.meta.update_files()
        util.print_time(start, 1)

    def _initialise_indices(self):
        """
        Initialises the RNA-seq index, the ChIP-seq histone index and the differential cell transitions if required.
        This needs to happen after the download pipeline, in order to only consider valid and successfully downloaded
        files.
        """
        print()
        print('INDEX INITIALISATION')

        # RNA INDEX
        start = time()
        # shortcut of the relevant configuration section
        paths = self.cfg['file paths']['rna seq quantification']                            # type: dict[str, str]
        # initialise the index
        self.rna_index = RnaIndex(df=self.meta.df, merge=self.args.rna_quant_merge_technical,
                                  dirs=paths, name_dict=self.cfg['file names']['quantification name components'])
        print('\tRNA-seq index initialised')
        util.print_time(start, 1)

        # HISTONE INDEX
        start = time()
        # shortcut of the relevant configuration section
        paths = self.cfg['file paths']['differential histone']['pre-processing']            # type: dict[str, str]
        quant_paths = self.cfg['file paths']['differential histone']['quantification']      # type: dict[str, str]
        # initialise the index
        self.histone_index = HistoneIndex(df=self.meta.df, dirs=paths)
        self.histone_index.load_quantification(dirs=quant_paths)
        print('\thistone index initialised')
        util.print_time(start, 1)

        # TRANSITIONS
        # for each transition, initialise the file paths
        for cell_1, cell_2 in [transition.split(' --> ') for transition in self.cfg['transitions']]:
            self.transitions[(cell_1, cell_2)] = Differential(cell_1=cell_1, cell_2=cell_2, cfg=self.cfg)

        if any([self.args.histone_differential, self.args.analysis]):
            for transition in self.transitions.values():
                transition.load_histone_comparisons(index=self.histone_index, cfg=self.cfg)
        print('\tinitialised differentiation transitions')
        util.print_time(start, 1)

        # REFERENCE
        start = time()
        if self.args.reference or not os.path.isfile(self.cfg['file paths']['reference']['genes']):
            self.reference = NewReference(self.cfg, self.args.epd, self.args.protein_coding)
            print('\tgenerated the reference genome')
        else:
            self.reference = Reference(self.cfg)
            print('\tloaded the reference genome and regulation')
        util.print_time(start, 1)

    def _download_pipeline(self):
        """
        Pipeline for downloading valid files and removing incorrectly downloaded files.
        """
        print()
        print('DOWNLOAD')
        start = time()

        if not (self.args.download or self.args.download_clean_invalid):
            self.meta.filter_download()
            print('\tfiltered for files that were successfully downloaded and valid')
            return

        # download valid files of the specified assays
        if self.args.download:
            download.download(self.meta, self.args.download_retry)
            print()

        # compute the files that are not valid downloads or were not downloaded successfully
        download.compute_invalid_downloads(self.cfg, self.meta, self.args.download_clean_invalid)

        # remove sub directories that were created during the download steps but do not contain files
        util.remove_empty_sub_dirs(self.cfg['file paths']['data directory'])
        print('\tremoved empty data subdirectories')

        util.print_time(start)

    def _rna_seq_quantification_pipeline(self):
        """
        Pipeline for quantifying raw RNA-seq data.
        """
        print()
        print('RNA-SEQ QUANTIFICATION')

        if not (self.args.rna_quant_index or self.args.rna_quant):
            print('\tnothing to do here')
            return

        # construct the salmon transcriptome index if specified
        if self.args.rna_quant_index:
            start = time()
            dea.salmon_index_command(self.cfg)
            util.print_time(start, 1)
            print()

        # quantify raw RNA-seq data with salmon if specified
        if self.args.rna_quant:
            start = time()
            dea.salmon_rna_seq_quantification(self.cfg, self.rna_index, self.args.rna_quant_fresh)
            util.print_time(start, 1)

    def _differential_expression_pipeline(self):
        """
        Pipeline for using previously computed RNA-seq quantification for differential expression analysis.
        """
        print()
        print('DIFFERENTIAL GENE EXPRESSION')
        start = time()

        if not self.args.differential_expression_analysis:
            print('\tnothing to do here')
            return

        # perform the differential gene expression analysis
        dea.differential_expression_analysis(self.cfg, self.rna_index, list(self.transitions.values()),
                                             self.args.dea_mapping)
        print()

        # output the time taken
        util.print_time(start)

    def _histone_pipeline(self):
        """
        Pipeline for differential histone modification analysis.
        """
        print()
        print('HISTONE MODIFICATIONS')

        if not any([self.args.histone_pre_process, self.args.histone_merge, self.args.histone_quant,
                    self.args.histone_merge_duplicates, self.args.histone_differential]):
            print('\tnothing to do here')
            return

        # pre-process histone .bam files, e.g. by removing duplicates and sorting
        if self.args.histone_pre_process:
            start = time()
            dha.histone_remove_duplicates(self.cfg, self.histone_index, merge=False)
            util.print_time(start, 1)
            print()

        # nothing more to do here
        if not any([self.args.histone_merge, self.args.histone_quant,
                    self.args.histone_merge_duplicates, self.args.histone_differential]):
            return

        # the histone quantification objects are needed for the next analysis steps
        self.histone_index.load_quantification(self.cfg['file paths']['differential histone']['quantification'])

        # merge the histone .bam files of each histone - cell combination
        if self.args.histone_merge:
            start = time()
            dha.merge_histone_alignments(self.cfg, self.histone_index)
            util.print_time(start, 1)
            print()

        # remove duplicates from merged .bam files and index them
        if self.args.histone_merge_duplicates:
            start = time()
            dha.histone_remove_duplicates(self.cfg, self.histone_index, merge=True)
            util.print_time(start, 1)
            print()

        # quantify the merged and indexed .bam files with histoneHMM
        if self.args.histone_quant:
            start = time()
            dha.histone_call_regions(self.cfg, self.histone_index, self.args.histone_merge)
            util.print_time(start, 1)
            print()

        # compute differential histone modifications with histoneHMM
        if self.args.histone_differential:
            start = time()
            dha.histone_differential(self.cfg, self.histone_index, self.transitions)
            util.print_time(start, 1)

    def _analysis_pipeline(self):
        print()
        print('ANALYSIS PIPELINE')

        if not any([self.args.analysis]):
            print('\tnothing to do here')
            return

        print('\tPRE-PROCESSING')
        manager = Manager(self.cfg, self.args.reload_results, self.args.recompute_networks)

        if self.args.reload_results or self.args.recompute_networks:
            processing.process(self.cfg, manager, self.reference, self.transitions)
        else:
            print('\tnothing to do here')

        print()
        print('\tOVERVIEWS')
        analysis.motif_overview(self.cfg, manager, list(self.transitions.values()))
        analysis.bubble_plot(manager)
        analysis.correlation_boxplot(manager, list(self.transitions.values()))
        analysis.p_fisher_diff(manager)
        analysis.p_steiger(manager)
        analysis.p_fisher_diff2(manager)
        analysis.quant_vs_cat(manager, list(self.transitions.values()))
        analysis.pseudogenes(manager, self.reference, list(self.transitions.values()))


def get_parameters():
    """
    Parses the command line arguments.
    :return: namespace with the command line arguments
    """
    parser = ap.ArgumentParser()
    parser.add_argument('config', type=str, help='path to the configuration file')
    # metadata run flags
    parser.add_argument('--meta-filter', action='store_true', help='filter the input metadata')
    parser.add_argument('--reference', action='store_true', help='recompute the gene definitions')
    parser.add_argument('--epd', action='store_true', help='add EPD promoter')
    parser.add_argument('--protein-coding', action='store_true', help='only consider protein coding and miRNA genes')
    # download flags
    parser.add_argument('--download', action='store_true', help='download the specified files')
    parser.add_argument('--download-retry', action='store_true', help='retry downloading files')
    parser.add_argument('--download-clean-invalid', action='store_true',
                        help='remove files that should not have been downloaded')
    # RNA-seq quantification
    parser.add_argument('--rna-quant', action='store_true', help='run RNA-seq quantification with salmon')
    parser.add_argument('--rna-quant-index', action='store_true', help='construct the quantification index')
    parser.add_argument('--rna-quant-fresh', action='store_true', help='delete previous results before quantifying')
    parser.add_argument('--rna-quant-merge-technical', action='store_true', help='merge technical replicates')
    # differential expression analysis
    parser.add_argument('--differential-expression-analysis', action='store_true',
                        help='run differential gene expression analysis')
    parser.add_argument('--dea-mapping', action='store_true', help='recompute the transcript -> gene mapping')
    # histone modification analysis
    parser.add_argument('--histone-pre-process', action='store_true', help='pre-process alignment files')
    parser.add_argument('--histone-merge', action='store_true', help='merge the pre-processed alignment .bam files')
    parser.add_argument('--histone-merge-duplicates', action='store_true',
                        help='remove duplicates in the merged .bam files')
    parser.add_argument('--histone-quant', action='store_true', help='quantify histone modifications')
    parser.add_argument('--histone-differential', action='store_true', help='differential histone quantification')

    parser.add_argument('--analysis', action='store_true', help='run the analysis pipeline')
    parser.add_argument('--recompute-networks', action='store_true', help='re-compute randomised networks')
    parser.add_argument('--reload-results', action='store_true',
                        help='reload differential expression and histone modification results')

    parser.add_argument('--beclear', action='store_true', help='run BEclear')

    return parser.parse_args()


if __name__ == '__main__':
    Pipeline(args=get_parameters())
