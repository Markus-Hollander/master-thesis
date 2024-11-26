from csv import dictReader
from statistics import stdev
from collections import namedtuple
from collections import defaultdict

from computation import list_stats
from util import process_chromosomes, process_strand, get, intersect, overlap


def p(key, val):
    """
    Utility function for printing inputs in automatically detected formats.
    """
    # single label
    if isinstance(key, str):
        # value list
        if isinstance(val, list):
            stats = list_stats(val)
            std = stdev(val) if len(val) > 1 else 0
            print('\t--> {0} n={1:,} non-zero={2:,} %={3:,} '
                  'min={4:,} median={5:,} mean={6:,.2f} max={7:,} std={8:,.2f}'.format(key, *stats, std))
        # single value
        else:
            print('\t--> {0:>7,} {1:}'.format(len(val), key))
    # multiple labels
    elif len(key) == 2:
        print('\t--> {0}) {1}: min={2:,} median={3:,} mean={4:,} max={5:,}'.format(*key, *list_stats(val)[2:]))
    elif len(key) == 3:
        print('\t--> {0}) {1:>5,} {2}: {3}'.format(key[0], len(val), *key[1:]))


def get_components(name: str) -> list[str]:
    """
    Splits a given name into components based on specific delimiters and character transitions.
    :param name: string to split into components
    :return: list of components extracted from the input name
    """
    components = ['']                                       # type: list[str]

    for c in name:
        # check if the character is a delimiter ('-' or '*') and start a new component if the current one is not empty
        if c == '-' or c == '*':
            if components[-1]:
                components.append('')
        # add to the existing but empty current component
        elif not components[-1]:
            components[-1] += c
        # append the character if both the current character and the last component are alphabetic
        elif c.isalpha() and components[-1].isalpha():
            components[-1] += c
        # start a new component if the current character is a digit and the last component is alphabetic
        elif c.isdigit() and components[-1].isalpha():
            components.append(c)
        # append the character if both the current character and the last component are numeric
        elif c.isdigit() and components[-1].isdigit():
            components[-1] += c
        # start a new component for all other cases
        else:
            components.append(c)

    # filter out any empty strings from the components list
    return [x for x in components if x]


def generate_mir_base_id(name: str) -> str:
    """
    Generates a miRBase-compliant identifier from a given name.
    :param name: input name to convert into a miRBase ID
    :return: formatted miRBase ID
    """
    # return the original name if it does not start with 'MIR' (case-insensitive)
    if not name.upper().startswith('MIR'):
        return name

    # get the components of the name
    components = get_components(name.lower())

    # determine the base prefix ('hsa-let' for 'let', otherwise 'hsa-mir')
    if 'let' in components[0]:
        components[0] = 'hsa-let'
    else:
        components[0] = 'hsa-mir'

    # identify the types of each component (True if alphabetic, False if numeric)
    component_types = [not x.isdigit() for x in components]

    # format the components based on their identified types
    if component_types == [True, False]:
        return '-'.join(components)
    elif component_types == [True, False, True]:
        return '{0}-{1}{2}'.format(*components)
    elif component_types == [True, False, False]:
        return '{0}-{1}-{2}'.format(*components)
    elif component_types == [True, False, True, False]:
        return '{0}-{1}{2}-{3}'.format(*components)
    else:
        raise ValueError('generate miRBase ID: {0}, {1}'.format(name, components))


def name_from_mir_base_id(name: str) -> str:
    """
    Converts a miRBase ID back to a simplified name format.
    :param name:  input miRBase ID to simplify
    :return: simplified name derived from the miRBase ID
    """
    # return the original name if it does not start with 'hsa-'
    if not name.startswith('hsa-'):
        return name

    # extract components by removing prefixes and case transformations
    components = get_components(name.replace('hsa-', '').replace('-3p', '').replace('-5p', '').upper())

    # change 'LET' to 'MIRLET' in the first component
    if components[0] == 'LET':
        components[0] = 'MIRLET'

    # identify the types of each component (True if alphabetic)
    component_types = [x.isalpha() for x in components]

    # reformat the name based on the identified component types
    if component_types in [[True, False], [True, False, True]]:
        return ''.join(components)
    elif component_types == [True, False, False]:
        return '{0}{1}-{2}'.format(*components)
    elif component_types == [True, False, True, False]:
        return '{0}{1}{2}-{3}'.format(*components)
    else:
        return 'ERROR'


def remove_mir_base_number(name: str) -> str:
    """
    Removes the numeric suffix from a miRBase ID to generate a base name.
    :param name: input miRBase ID to simplify
    :return: simplified name without numeric suffix
    """
    # get the components of the name
    components = get_components(name)
    # identify the types of each component (True if alphabetic)
    component_types = [x.isalpha() for x in components]

    # reformat the name based on its identified component types
    if component_types in [[True, False, True, False], [True, False, True], [True, False]]:
        return ''.join(components[:min(len(components), 3)])
    elif component_types == [True, False, False]:
        return ''.join(components[:-1])
    else:
        raise ValueError('remove number from miRBase ID: {0}, {1}'.format(name, components))


class Gene:
    """ Represents a gene with detailed attributes such as position, strand, and promoter information. """
    def __init__(self, gene_info: dict[str, str], p_start=0, p_end=0, mp_start=0, mp_end=0):
        """
        :param gene_info: dictionary with gene information
        :param p_start: default relative promoter start
        :param p_end: default relative promoter end
        :param mp_start: default relative miRNA promoter start
        :param mp_end: default relative miRNA promoter end
        """
        # basic gene information such as Ensembl ID, gene name, description and gene type
        self.id = gene_info['gene.id']                                                                      # type: str
        self.name = gene_info['gene.name']                                                                  # type: str
        self.description = gene_info['gene.description']                                                    # type: str
        self.type = gene_info['gene.type']                                                                  # type: str

        # construct the default miRBase IDs in case they are not given in the gene information
        default = '' if self.type != 'miRNA' else generate_mir_base_id(self.name)                        # type: str
        # extract the miRBase ID if given, otherwise use the default ID (for non-miRNAs)
        self.mir_base_id = get(gene_info, 'mir.base.id', default)                                           # type: str

        # gene start and end position on the respective chromosome, as well as gene length
        # coordinates are 1-based, meaning the end position is included in the length
        self.start = int(gene_info['gene.start'])                                                           # type: int
        self.end = int(gene_info['gene.end'])                                                               # type: int
        self.length = get(gene_info, 'gene.length', abs(self.start - self.end) + 1)                         # type: int

        # extract the chromosome and the strand, and make sure the strand is valid
        self.chromosome = process_chromosomes(gene_info['gene.chromosome'])                                 # type: str
        self.strand = process_strand(gene_info['gene.strand'])                                              # type: str

        if self.strand not in ['+', '-']:
            raise ValueError('Wrong strand: {0}'.format(self.strand))

        # True if the gene is on the + strand, False otherwise
        self.is_plus = self.strand == '+'                                                                   # type: bool

        # extract or set the TSS and compute the distance to the gene start (+) or end (-)
        self.tss = get(gene_info, 'gene.TSS', self.start if self.is_plus else self.end)                     # type: int
        self.tss_distance = get(gene_info, 'TSS.distance',
                                self.tss - self.start if self.is_plus else self.end - self.tss)             # type: int

        # adjust promoter start and end positions for miRNA genes
        if self.type == 'miRNA':
            p_start = mp_start
            p_end = mp_end

        # promoter start, end, length, and overlap information
        self.promoter_start = get(gene_info, 'promoter.start',
                                  self.tss + p_start if self.is_plus else self.tss - p_end)                 # type: int
        self.promoter_end = get(gene_info, 'promoter.end',
                                self.tss + p_end if self.is_plus else self.tss - p_start)                   # type: int
        self.promoter_length = get(gene_info, 'promoter.length',
                                   abs(self.promoter_start - self.promoter_end) + 1)                        # type: int
        self.promoter_overlap = get(gene_info, 'promoter.overlap',
                                    overlap(self.start, self.end, self.promoter_start, self.promoter_end))  # type: int

        # number of transcripts associated with the gene
        self.transcript_count = get(gene_info, 'transcript.count', 0)                                       # type: int
        # gene ontology (GO) identifiers
        self.GO_ids = get(gene_info, 'GO.ids', set())                                                   # type: set[str]

        # relationships with other genes
        self.contained_in = get(gene_info, 'contained.in', set())                                       # type: set[str]
        self.contains = get(gene_info, 'contains', set())                                               # type: set[str]
        self.overlaps = get(gene_info, 'overlaps', set())                                               # type: set[str]
        self.identical = get(gene_info, 'identical.range', set())                                       # type: set[str]

        # host gene information (if applicable)
        self.host_gene_id = get(gene_info, 'host.gene.id', '')                                              # type: str
        self.host_gene_name = get(gene_info, 'host.gene.name', '')                                          # type: str

        # flags for gene properties
        self.is_TF = False                                                                                  # type: bool
        self.is_miRNA = self.type == 'miRNA'                                                                # type: bool

    def change_promoter(self, tss: int, promoter_start: int, promoter_end: int):
        """
        Updates the promoter region of the gene.
        :param tss: new transcription start site
        :param promoter_start: new promoter start position
        :param promoter_end: new promoter end position
        :return:
        """
        self.tss = tss
        self.tss_distance = self.tss - self.start if self.is_plus else self.end - self.tss
        
        self.promoter_start = promoter_start
        self.promoter_end = promoter_end
        self.promoter_length = abs(self.promoter_start - self.promoter_end) + 1
        self.promoter_overlap = overlap(self.start, self.end, self.promoter_start, self.promoter_end)
    
    def change_tss(self, tss: int, p_start: int, p_end: int, mp_start: int, mp_end: int):
        """
        Updates the transcription start site and adjusts the promoter region.
        :param tss: new transcription start site
        :param p_start: relative promoter start
        :param p_end: relative promoter end
        :param mp_start: relative miRNA promoter start
        :param mp_end: relative miRNA promoter end
        """
        self.tss = tss                                                                                      # type: int
        self.tss_distance = self.tss - self.start if self.is_plus else self.end - self.tss                  # type: int

        # use miRNA-specific promoter bounds if applicable
        if self.type == 'miRNA':
            p_start = mp_start
            p_end = mp_end

        self.promoter_start = self.tss + p_start if self.is_plus else self.tss - p_end                      # type: int
        self.promoter_end = self.tss + p_end if self.is_plus else self.tss - p_start                        # type: int
        self.promoter_length = abs(self.promoter_start - self.promoter_end) + 1                             # type: int
        self.promoter_overlap = overlap(self.start, self.end, self.promoter_start, self.promoter_end)       # type: int

    def dist(self, start: int, end: int):
        """ Returns the minimum distance from the gene to a given region. """
        return min([abs(self.start - start), abs(self.end - end)])

    def tss_dist(self, tss: int):
        """ Returns the minimum distance from the gene to the given transcription start site. """
        return min([abs(x - tss) for x in [self.start, self.end]])

    def __eq__(self, other):
        return self.id == other.id

    def __hash__(self):
        return hash(self.id)


class RegulationEntry:
    """ Represents an entry for regulatory interaction between a regulator and a target. """
    def __init__(self, line: dict[str, str]):
        """
        :param line: dictionary containing regulatory entry details
        """
        # convert regulator and target names to their canonical forms
        self.regulator = name_from_mir_base_id(line['regulator'])
        self.target = name_from_mir_base_id(line['target'])

        # category of regulation (e.g., activation, repression)
        self.category = line['category']
        # flag indicating if the evidence is experimental
        self.experimental = line['evidence'] == 'Experimental'
        # source of the regulatory information
        self.source = line['source']

        # store raw regulator and target names as provided in the input
        self.raw_regulator = line['regulator']
        self.raw_target = line['target']

    def __eq__(self, other):
        if self.regulator != other.regulator:
            return False
        if self.target != other.target:
            return False
        if self.category != other.category:
            return False
        if self.experimental != other.experimental:
            return False
        if self.source != other.source:
            return False
        return True

    def __hash__(self):
        return hash((self.regulator, self.target, self.category, self.experimental, self.source))


class Reference:
    """ Represents a reference containing gene information and regulatory data. """
    def __init__(self, cfg: dict, epd=False, protein_coding=False):
        """
        :param cfg: dictionary with configuration file paths and settings
        :param epd: flag indicating whether to use EPD (default is False)
        :param protein_coding: flag indicating whether to consider only protein-coding genes (default is False)
        """
        # stores genes by their unique IDs
        self.genes = dict()                                     # type: dict[str, Gene]
        # stores sets of genes indexed by their names
        self.genes_by_name = defaultdict(set)                   # type: dict[str, set[Gene]]

        # stores regulatory edges grouped by category
        self.edges = {'all': defaultdict(set),
                      'experimental': defaultdict(set)}      # type: dict[str, dict[str, set[tuple[str, str]]]]
        # stores regulation statistics
        self.regulation_stats = dict()                          # type: dict[str, str]

        self._load(cfg, epd, protein_coding)
        self._load_regulation(cfg)

    def _load(self, cfg: dict, epd=False, protein_coding=False):
        """
        Loads gene data from the reference file specified in the configuration.

        :param cfg: dictionary with configuration file paths and settings
        :param epd: flag indicating whether to use EPD (default is False)
        :param protein_coding: flag indicating whether to consider only protein-coding genes (default is False)
        """
        with open(cfg['file paths']['reference']['genes'], 'r') as file:
            for line in dictReader(file, delimiter='\t'):
                gene = Gene(line)
                self.genes[gene.id] = gene
                self.genes_by_name[gene.name].add(gene)

    def _load_regulation(self, cfg: dict):
        """
        Loads regulation data and categorizes regulatory interactions.
        :param cfg: dictionary with configuration file paths and settings
        """
        with open(cfg['file paths']['reference']['regulation'], 'r', encoding='latin-1') as file:
            # create a set of unique RegulationEntry objects
            unique_entries = {RegulationEntry(x) for x in dictReader(file, delimiter='\t')}
            # filter out entries with errors in their regulator or target names
            unique_entries = {entry for entry in unique_entries if 'ERROR' not in [entry.regulator, entry.target]}

        # filter entries to include only those with valid regulators and targets in the reference
        entries = {e for e in unique_entries
                   if e.regulator in self.genes_by_name and e.target in self.genes_by_name}

        # initialize sets for network nodes, regulators, and targets
        nodes = set()
        all_regulators = set()
        all_targets = set()

        # mark regulators as transcription factors (TFs) based on their categories
        for entry in entries:
            for regulator in self.genes_by_name[entry.regulator]:
                if entry.category.startswith('tf'):
                    regulator.is_TF = True

        # categorize regulatory interactions and populate edge data
        for entry in entries:
            regulators = self.genes_by_name[entry.regulator]
            targets = self.genes_by_name[entry.target]

            for regulator in regulators:
                nodes.add(regulator)
                all_regulators.add(regulator)

                for target in targets:
                    nodes.add(target)
                    all_targets.add(target)

                    # determine interaction category
                    if regulator.is_TF and target.is_miRNA:
                        category = 'tf-mirna'
                    elif regulator.is_TF and target.is_TF:
                        category = 'tf-tf'
                    elif regulator.is_TF and not target.is_miRNA and not target.is_TF:
                        category = 'tf-gene'
                    elif regulator.is_miRNA and target.is_miRNA:
                        category = 'mirna-mirna'
                    elif regulator.is_miRNA and target.is_TF:
                        category = 'mirna-tf'
                    elif regulator.is_miRNA and not target.is_TF and not target.is_miRNA:
                        category = 'mirna-gene'
                    else:
                        category = 'gene-gene'

                    # add the interaction to the "all" edge set
                    self.edges['all'][category].add((regulator.id, target.id))

                    # add the interaction to the "experimental" edge set if it has experimental evidence
                    if entry.experimental:
                    if entry.experimental:
                        self.edges['experimental'][category].add((regulator.id, target.id))

        # write network statistics to a file
        with open('network_stats.txt', 'w') as file:
            file.write('number of nodes: {0:,}\n'.format(len(nodes)))
            file.write('number of regulators: {0:,}\n'.format(len(all_regulators)))
            file.write('number of targets: {0:,}\n'.format(len(all_targets)))
            file.write('number of miRNAs: {0:,}\n'.format(len({gene.id for gene in self.genes.values() if gene.is_miRNA})))
            file.write(
                'number of TFs: {0:,}\n'.format(len({gene.id for gene in self.genes.values() if gene.is_TF})))

            # write the number of edges for each category
            for category, edges in self.edges['experimental'].items():
                file.write('number of {0} edges: {1:,}\n'.format(category, len(edges)))

    def get_gene_types(self) -> list[str]:
        """ Retrieves a sorted list of all unique gene types. """
        return sorted({gene.type for gene in self.genes.values()})

    def write_gene_types(self, file_path: str):
        """
        Writes gene IDs and their types to a file.
        :param file_path: path to the file where gene types will be written
        """
        with open(file_path, 'w') as file:
            for gene in self.genes.values():
                if gene.is_TF:
                    t = 'TF'
                elif gene.is_miRNA:
                    t = 'miRNA'
                else:
                    t = 'other'

                file.write('{0}\t{1}\n'.format(gene.id, t))


MiRIAD = namedtuple('MiRIAD', ['host_gene_name', 'chromosome', 'strand', 'start', 'end'])


class FantomPromoter:
    """ Represents a Fantom promoter parsed from a formatted string. """
    def __init__(self, promoter: str):
        """
        :param promoter: string representing the promoter in the format 'name@chr:position,strand
        """
        self.input = promoter
        # process the promoter string
        promoter = promoter.split('@')[1]
        chrom, promoter = promoter.split(':')
        promoter, strand = promoter.split(',')
        start, end = promoter.split('..')

        # process chromosome name (e.g., standardize format) and strand information (e.g. +, -)
        self.chromosome = process_chromosomes(chrom)
        self.strand = process_strand(strand)
        self.start = int(start)
        self.end = int(end)


class FantomIntra:
    """ Represents an intra-chromosomal Fantom promoter group parsed from a formatted string. """
    def __init__(self, promoter: str):
        """
        :param promoter: string containing multiple promoters separated by commas
        """
        self.input = promoter
        self.names = {part.split('@')[1] for part in promoter.split(',')}


class Fantom:
    """ Represents a Fantom object that holds genomic data and promoter information. """
    def __init__(self, line: dict[str, str]):
        """
        :param line: dictionary with keys for chromosome, strand, start, end, tss, and promoter
        """
        # extract chromosome and strand information from the input
        self.input = line
        self.chromosome = process_chromosomes(line['chromosome'])
        self.strand = process_strand(line['strand'])
        self.start = int(line['start'])
        self.end = int(line['end'])

        # initialize transcription start site (TSS) and TSS availability flag
        self.tss = 0
        self.no_tss = True

        # check if TSS is valid and update attributes
        if line['tss'] not in ['Unknown', '', None]:
            self.tss = int(line['tss'])
            self.no_tss = False

        # initialize promoter and host names
        self.promoter = None                                # type: FantomPromoter | None
        self.host_names = set()
        # initialize error flag
        self.error = False

        # process the promoter field
        promoter = line['promoter']
        if promoter in ['Unknown', '', None]:
            self.error = True
            return

        # check if the promoter contains chromosome information
        if '@chr' in promoter:
            # create a FantomPromoter object for inter-chromosomal promoters
            self.promoter = FantomPromoter(promoter)
            self.intra = False
        else:
            # extract host names for intra-chromosomal promoters
            self.host_names = {part.split('@')[1] for part in promoter.split(',')}
            self.intra = True
            # set error flag if no valid host names are found
            if not self.host_names:
                self.error = True

    def __eq__(self, other):
        return self.input == other.input

    def copy(self):
        """ Creates a deep copy of the current Fantom object. """
        return Fantom(self.input)


class MiRNA:
    """ Represents a microRNA (miRNA) and its associated data, including miRBase ID and gene reference. """
    def __init__(self, mir_base_id: str, gene: Gene):
        self.mir_base_id = mir_base_id      # type: str
        self.gene = gene                    # type: Gene

        self.miriad = None                  # type: MiRIAD | None
        self.fantom = None                  # type: Fantom | None
        self.fantom_hg19 = None             # type: Fantom | None

        # initialize sets to store matching IDs for MiRIAD and Fantom
        self.miriad_match = set()           # type: set[str]
        self.fantom_match = set()           # type: set[str]

    def fantom_intra(self):
        """ Checks if the Fantom data contains intra-chromosomal promoter information. """
        return self.fantom and self.fantom.host_names

    def fantom_inter(self):
        """ Checks if the Fantom data contains inter-chromosomal promoter information. """
        return self.fantom and self.fantom.promoter

    def __eq__(self, other):
        if self.miriad != other.miriad:
            return False
        if self.fantom != other.fantom:
            return False
        return True

    def __hash__(self):
        return hash(self.mir_base_id)


class FantomProcessor:
    """ Processes Fantom data and builds mappings between miRNA IDs and Fantom objects. """
    def __init__(self, cfg):
        # initialize dictionaries for hg19 and hg38 mappings
        self.hg19 = dict()      # type: dict[str, Fantom]
        self.hg38 = dict()      # type: dict[str, Fantom]

        # load Fantom data for hg19 from the specified file
        with open(cfg['file paths']['gene definitions']['fantom mirna promoter'], 'r') as file:
            for line in dictReader(file, delimiter='\t'):
                entry = Fantom(line)

                if entry.error:
                    continue

                self.hg19[line['mir.base.id']] = entry

        # write Fantom data for hg19 to a BED file
        with open(cfg['file paths']['gene definitions']['fantom bed hg19'], 'w') as file:
            for key, entry in self.hg19.items():
                file.write('\t'.join([str(entry.chromosome), str(entry.start),
                                      str(entry.end + 1), key + '/gene']) + '\n')

                # write the TSS region if available
                if entry.tss:
                    file.write('\t'.join([str(entry.chromosome), str(entry.tss),
                                          str(entry.tss + 1), key + '/tss']) + '\n')
                # skip if the entry is intra-chromosomal
                if entry.intra:
                    continue

                file.write('\t'.join([str(entry.promoter.chromosome), str(entry.promoter.start),
                                      str(entry.promoter.end + 1), key + '/prom']) + '\n')

        # load Fantom data for hg38 from the specified BED file
        with open(cfg['file paths']['gene definitions']['fantom bed hg38'], 'r') as file:
            for line in file:
                # split the line into chromosome, positions, and tag
                chromosome, start, end, tag = line.strip().split('\t')
                start = int(start)
                end = int(end)
                mir_base_id, location = tag.split('/')

                # check if the miRBase ID already exists in the hg38 dictionary
                if mir_base_id in self.hg38:
                    entry = self.hg38[mir_base_id]
                else:
                    entry = self.hg19[mir_base_id].copy()

                # update the Fantom entry based on the location type
                if location == 'tss':
                    entry.tss = start
                elif location == 'prom':
                    entry.promoter.start = start
                    entry.promoter.end = end - 1
                else:
                    entry.start = start
                    entry.end = end - 1

                self.hg38[mir_base_id] = entry


class NewReference(Reference):
    """ Extends the Reference class to handle additional loading and processing of gene and promoter data. """
    def __init__(self, cfg: dict, epd: bool, protein_coding: bool):
        super().__init__(cfg, epd, protein_coding)

    def _load(self, cfg: dict, epd=False, protein_coding=False):
        """
        Loads all required data for processing gene and promoter information.
        :param cfg: configuration dictionary containing paths and settings for gene and promoter data
        :param epd: boolean indicating whether to load EPD promoters
        :param protein_coding: boolean indicating whether to include only protein-coding genes
        """
        self._load_gene_info(cfg, protein_coding)

        if epd:
            self._load_epd_promoters(cfg)

        self._load_overlap(cfg)
        self._load_mirna_promoters(cfg)
        self._write_definitions(cfg)

    def _load_gene_info(self, cfg: dict, protein_coding: bool):
        """
        Loads gene information, including positions, types, and annotations.
        :param cfg: configuration dictionary containing paths and settings for gene data
        :param protein_coding: boolean indicating whether to include only protein-coding genes
        """
        # retrieve promoter start and end positions from the configuration
        self.p_start = int(cfg['analysis']['promoter start'])
        self.p_end = int(cfg['analysis']['promoter end'])
        self.mp_start = int(cfg['analysis']['mirna promoter start'])
        self.mp_end = int(cfg['analysis']['mirna promoter end'])

        # load gene definitions from the specified file
        with open(cfg['file paths']['gene definitions']['ensembl definitions'], 'r') as file:
            self.genes = {line['gene.id']: Gene(line, self.p_start, self.p_end, self.mp_start, self.mp_end)
                          for line in dictReader(file, delimiter='\t')
                          if line['gene.type'] in ['protein_coding', 'miRNA']
                          or not protein_coding}                                        # type: dict[str, Gene]

        # organize genes by their names for quick lookup
        for gene in self.genes.values():
            self.genes_by_name[gene.name].add(gene)

        # load gene annotations from the specified file
        with open(cfg['file paths']['gene definitions']['ensembl annotations'], 'r') as file:
            annotations = {(line['gene.id'], line['GO.accession']) for line in dictReader(file, delimiter='\t')
                           if line['GO.accession'] and line['gene.id'] in self.genes}

            for gene_id, annotation in annotations:
                self.genes[gene_id].GO_ids.add(annotation)

    def _load_epd_promoters(self, cfg: dict):
        """
        Loads EPD promoters and filters them based on compatibility with Ensembl gene data.
        :param cfg: configuration dictionary containing paths for EPD data
        """
        # define a named tuple to represent an EPD promoter
        EpdPromoter = namedtuple('EpdPromoter', ['chromosome', 'tss'])

        # helper function to determine if TSS positions are compatible
        def f(tss_gene: int, tss_epd: int, is_plus: bool) -> bool:
            if is_plus:
                return tss_epd <= tss_gene
            return tss_epd >= tss_gene

        # load EPD promoter data from the specified file
        promoters_name = dict()
        with open(cfg['file paths']['gene definitions']['epd tss'], 'r') as file:
            for line in file:
                line = line.strip().split('\t')

                promoters_name[line[3].split('_')[0]] = EpdPromoter(chromosome=process_chromosomes(line[0]),
                                                                    tss=int(line[1]))

        # match promoters to genes by name and ensure chromosome compatibility
        promoters_id = {list(self.genes_by_name[k])[0].id: v for k, v in promoters_name.items()
                        if k in self.genes_by_name and len(self.genes_by_name[k]) == 1
                        and list(self.genes_by_name[k])[0].chromosome == v.chromosome}  # type: dict[str, EpdPromoter]

        # filter promoters based on TSS compatibility
        promoters_filtered = {k: v for k, v in promoters_id.items()
                              if f(self.genes[k].tss, v.tss, self.genes[k].is_plus)}
        promoters_rejected = {k: v for k, v in promoters_id.items()
                              if not f(self.genes[k].tss, v.tss, self.genes[k].is_plus)}

        # print diagnostic information about the loaded promoters
        print('\tEPD')
        p('promoters', promoters_name.keys())
        p('with only 1 gene on same chromosome', promoters_id.keys())
        p('distance on same chromosome', promoters_filtered.keys())
        p('miRNA promoters', {k for k in promoters_filtered.keys() if self.genes[k].type == 'miRNA'})

        print('\tdistance to Ensembl TSS')
        p('accepted +:', [v.tss - self.genes[k].tss for k, v in promoters_filtered.items()
                          if self.genes[k].is_plus])
        p('accepted -:', [self.genes[k].tss - v.tss for k, v in promoters_filtered.items()
                          if not self.genes[k].is_plus])
        p('rejected +:', [v.tss - self.genes[k].tss for k, v in promoters_rejected.items()
                          if self.genes[k].is_plus])
        p('rejected -:', [self.genes[k].tss - v.tss for k, v in promoters_rejected.items()
                          if not self.genes[k].is_plus])

        top = 30
        print('\tdistance top {0}'.format(top))
        p('accepted +:', sorted([v.tss - self.genes[k].tss for k, v in promoters_filtered.items()
                                 if self.genes[k].is_plus])[:top])
        p('accepted -:', sorted([self.genes[k].tss - v.tss for k, v in promoters_filtered.items()
                                 if not self.genes[k].is_plus])[:top])
        p('rejected +:', [v.tss - self.genes[k].tss for k, v in promoters_rejected.items()
                          if self.genes[k].is_plus and self.genes[k].tss >= v.tss])
        p('rejected -:', [self.genes[k].tss - v.tss for k, v in promoters_rejected.items()
                          if not self.genes[k].is_plus and self.genes[k].tss <= v.tss])

        # adjust TSS for genes based on filtered promoters
        for k, v in promoters_filtered.items():
            self.genes[k].change_tss(v.tss, self.p_start, self.p_end, self.mp_start, self.mp_end)

    def _load_overlap(self, cfg: dict):
        """
        Identifies overlaps, containment, and identical relationships among genes.
        :param cfg: configuration dictionary containing paths for Ensembl BED and intersect data
        """
        # define paths for BED and intersect files
        bed_path = cfg['file paths']['gene definitions']['ensembl bed']
        intersect_path = cfg['file paths']['gene definitions']['ensembl intersect']

        # write BED file for genes
        with open(bed_path, 'w') as file:
            for gene in self.genes.values():
                columns = [gene.chromosome, gene.start, gene.end, gene.id, '.', gene.strand]
                file.write('\t'.join([str(c) for c in columns]) + '\n')

        # compute intersections using an external tool
        intersect(bed_path, bed_path, intersect_path, strand=False)

        # process intersect results to classify relationships
        with open(intersect_path, 'r') as f:
            for line in f:
                _, gs_a, ge_a, id_a, _, _, _, gs_b, ge_b, id_b, _, _, _ = line.strip().split('\t')

                gs_a, ge_a, gs_b, ge_b = [int(x) for x in [gs_a, ge_a, gs_b, ge_b]]

                if id_a == id_b:
                    continue

                if [gs_a, ge_a] == [gs_b, ge_b]:
                    self.genes[id_a].identical.add(id_b)
                    self.genes[id_b].identical.add(id_a)
                elif gs_a <= gs_b <= ge_b <= ge_a:
                    self.genes[id_a].contains.add(id_b)
                    self.genes[id_b].contained_in.add(id_a)
                elif gs_b <= gs_a <= ge_a <= ge_b:
                    self.genes[id_a].contained_in.add(id_b)
                    self.genes[id_b].contains.add(id_a)
                else:
                    self.genes[id_a].overlaps.add(id_b)
                    self.genes[id_b].overlaps.add(id_a)

    def _load_mirna_promoters(self, cfg: dict):
        """
        Loads miRNA promoter information from miRIAD and FANTOM datasets.
        :param cfg: configuration dictionary containing file paths and processing parameters
        """
        # initialize miRNA gene objects from the internal gene dictionary
        mirna_g = {gene.mir_base_id: gene for gene in self.genes.values() if gene.type == 'miRNA'}
        mirna = {mir_base_id: MiRNA(mir_base_id, gene) for mir_base_id, gene in
                 mirna_g.items()}  # type: dict[str, MiRNA]

        missing_mirna_ids = set()

        # process miRIAD host gene data
        with open(cfg['file paths']['gene definitions']['miriad host genes intra'], 'r') as file:
            for line in dictReader(file, delimiter='\t'):
                # skip non-human entries or entries not present in miRNA
                if any([line['organism'] != 'human', line['mir.base.id'] not in mirna]):
                    if line['organism'] == 'human' and line['mir.base.id'] not in mirna:
                        missing_mirna_ids.add(line['mir.base.id'])
                    continue

                # identify the host gene if present in the gene list
                if line['host.gene.name'] in self.genes_by_name:
                    name = line['host.gene.name']
                else:
                    continue

                # extract the strand direction, defaulting to '+-' if not specified
                if 'direction' in line:
                    strand = process_strand(line['direction'])
                else:
                    strand = '+-'

                # create a MiRIAD entry for the miRNA
                entry = MiRIAD(host_gene_name=name, chromosome=process_chromosomes(line['chromosome']),
                               start=int(line['start']), end=int(line['end']), strand=strand)
                mirna[line['mir.base.id']].miriad = entry

        # process FANTOM data for miRNA promoters in both hg38 and hg19 assemblies
        fantom_processor = FantomProcessor(cfg)
        for mir_base_id, entry in fantom_processor.hg38.items():
            if mir_base_id in mirna:
                # filter entries that are too far apart from the miRNA gene
                if mirna[mir_base_id].gene.dist(entry.start, entry.end) > 50000:
                    continue
                mirna[mir_base_id].fantom = entry

        for mir_base_id, entry in fantom_processor.hg19.items():
            if mir_base_id in mirna:
                mirna[mir_base_id].fantom_hg19 = entry

        # validate matches between miRIAD and FANTOM data
        flags = {'fantom': defaultdict(set), 'miriad': defaultdict(set), 'both': defaultdict(set)}
        for m_id, m in mirna.items():
            f_hosts = [y for x in m.fantom.host_names for y in self.genes_by_name[x]] if m.fantom else []
            m_hosts = list(self.genes_by_name[m.miriad.host_gene_name]) if m.miriad else []

            for k, hosts in [('fantom', f_hosts), ('miriad', m_hosts), ('both', m_hosts + f_hosts)]:
                host_ids = {host.id for host in hosts if host.type != 'miRNA'}

                # categorize based on the number of host IDs
                if not host_ids:
                    flags[k][('0', 'host IDs', '< 1')].add(m.gene.id)
                    continue
                elif len(host_ids) > 1:
                    flags[k][('0', 'host IDs', '> 1')].add(m.gene.id)
                else:
                    flags[k][('0', 'host IDs', '= 1')].add(m.gene.id)

                # classify matches into tasks such as containment, overlap, etc.
                tasks = [(('1', 'contained match'), m.gene.contained_in.intersection(host_ids)),
                         (('2', 'overlap match'), m.gene.overlaps.intersection(host_ids)),
                         (('3', 'host gene on same chromosome'), {h.id for h in hosts
                                                                  if h.chromosome == m.gene.chromosome}),
                         (('4', 'no match'), set())]

                for label, matches in tasks:
                    matches = {gene_id for gene_id in matches if self.genes[gene_id].type != 'miRNA'}
                    if not matches:
                        if label[0] == '4':
                            flags[k][label + ('',)].add(m.gene.id)
                        else:
                            flags[k][label + ('< 1',)].add(m.gene.id)
                        continue
                    elif len(matches) > 1:
                        flags[k][label + ('> 1',)].add(m.gene.id)
                    else:
                        flags[k][label + ('= 1',)].add(m.gene.id)

                    if k == 'fantom':
                        m.fantom_match = matches
                    elif k == 'miriad':
                        m.miriad_match = matches
                    break

        # print statistics and debug information for miRNAs
        print('\n\tmiRNA stats')
        p('miRNAs', set(mirna_g.keys()))
        p('contained in others', {k for k, g in mirna_g.items() if g.contained_in})
        print()
        p('miRNA IDs', mirna.keys())
        p('miRIAD and FANTOM', {x.mir_base_id for x in mirna.values() if x.miriad and x.fantom})
        p('neither', {x.mir_base_id for x in mirna.values() if not x.miriad and not x.fantom})
        print('\tmiRIAD')
        p('intergenic', {key for key, val in mirna.items() if not val.miriad})
        p('intragenic', {key for key, val in mirna.items() if val.miriad})
        print('\tFANTOM')
        p('intergenic', {key for key, val in mirna.items() if val.fantom_inter()})
        p('intragenic', {key for key, val in mirna.items() if val.fantom_intra()})
        print('\texclusively in miRIAD or FANTOM')
        p('miRIAD inter', {key for key, val in mirna.items() if not val.miriad and not val.fantom})
        p('miRIAD intra', {key for key, val in mirna.items() if val.miriad and not val.fantom})
        p('FANTOM inter', {key for key, val in mirna.items() if val.fantom_inter() and not val.miriad})
        p('FANTOM intra', {key for key, val in mirna.items() if val.fantom_intra() and not val.miriad})
        print('\tmiRIAD/FANTOM')
        short = [x for x in mirna.values() if x.fantom and x.miriad]
        p('inter/inter', {x.mir_base_id for x in short if not x.miriad and x.fantom_inter()})
        p('intra/inter', {x.mir_base_id for x in short if x.miriad and x.fantom_inter()})
        p('inter/intra', {x.mir_base_id for x in short if not x.miriad and x.fantom_intra()})
        p('intra/intra', {x.mir_base_id for x in short if x.fantom_intra() and x.miriad})

        # display matching information between miRIAD and FANTOM
        for k, val in flags.items():
            print('\t' + k)
            for k, v in sorted(val.items(), key=lambda x: x[0]):
                p(k, v)

        print('\tpromoter lengths')
        p('hg38:', [abs(x.fantom.promoter.start - x.fantom.promoter.end) for x in mirna.values()
                    if x.fantom_inter()])
        p('hg19:', [abs(x.fantom_hg19.promoter.start - x.fantom_hg19.promoter.end) for x in mirna.values()
                    if x.fantom_inter()])
        print('\tdistance of closest')
        p('miRIAD:', [x.gene.dist(x.miriad.start, x.miriad.end) for x in mirna.values() if x.miriad])
        p('hg38:  ', [x.gene.dist(x.fantom.start, x.fantom.end) for x in mirna.values() if x.fantom])
        p('hg19:  ', [x.gene.dist(x.fantom_hg19.start, x.fantom_hg19.end) for x in mirna.values() if x.fantom_hg19])
        print('\tdistance of tss of closest')
        p('hg38:  ', [x.gene.tss_dist(x.fantom.tss) for x in mirna.values() if x.fantom_inter() and x.fantom.tss])
        p('hg19:  ', [x.gene.tss_dist(x.fantom_hg19.tss) for x in mirna.values() if x.fantom_inter() and x.fantom_hg19.tss])
        print('\tdistance promoter closest')
        p('hg38:  ', [x.gene.dist(x.fantom.promoter.start, x.fantom.promoter.end) for x in mirna.values()
                      if x.fantom_inter()])
        p('hg19:  ', [x.gene.dist(x.fantom_hg19.promoter.start, x.fantom_hg19.promoter.end)
                      for x in mirna.values() if x.fantom_inter()])
        print('\tdistance promoter tss')
        p('hg38 start:', [x.fantom.promoter.start - x.fantom.tss for x in mirna.values() if x.fantom_inter()])
        p('hg19 start:', [x.fantom_hg19.promoter.start - x.fantom_hg19.tss for x in mirna.values() if x.fantom_inter()])
        p('hg38 end:  ', [x.fantom.promoter.end - x.fantom.tss for x in mirna.values() if x.fantom_inter()])
        p('hg19 end:  ', [x.fantom_hg19.promoter.end - x.fantom_hg19.tss for x in mirna.values() if x.fantom_inter()])

        print('\tcongruence miRIAD and FANTOM host genes')
        p('share (all)', {x.mir_base_id for x in mirna.values() if x.miriad_match == x.fantom_match
                          and x.miriad_match and x.fantom_match})
        p('share (some)', {x.mir_base_id for x in mirna.values() if x.miriad_match.intersection(x.fantom_match)
                           and x.miriad_match and x.fantom_match})
        p('different (all)', {x.mir_base_id for x in mirna.values() if x.miriad_match and x.fantom_match
                              and not x.miriad_match.intersection(x.fantom_match)})
        p('miRIAD only', {x.mir_base_id for x in mirna.values() if x.miriad_match and not x.fantom_match})
        p('FANTOM only', {x.mir_base_id for x in mirna.values() if x.fantom_match and not x.miriad_match})

        # process FANTOM and MiRIAD matches
        flags = defaultdict(set)
        for m in mirna.values():
            host_gene = None

            if m.fantom and not m.fantom.intra:
                flags['FANTOM intergenic'].add(m.gene.id)
                m.gene.change_promoter(m.fantom.tss, m.fantom.promoter.start, m.fantom.promoter.end)
                continue
            # single FANTOM match
            if len(m.fantom_match) == 1:
                host_gene = self.genes[list(m.fantom_match)[0]]
                flags['= 1 FANTOM match'].add(m.gene.id)
            # multiple FANTOM matches
            elif len(m.fantom_match) > 1:
                in_common = m.fantom_match.intersection(m.miriad_match)
                # one match in common between FANTOM and miRIAD
                if len(in_common) == 1:
                    host_gene = self.genes[list(m.fantom_match)[0]]
                    flags['> FANTOM match, 1 in common'].add(m.gene.id)
                else:
                    if not in_common:
                        flags['> 1 FANTOM match, none in common'].add(m.gene.id)
                        in_common = m.fantom_match
                    else:
                        flags['> 1 FANTOM match, > 1 in common'].add(m.gene.id)

                    host_gene = self.genes[sorted(in_common, key=lambda g: m.gene.tss_dist(self.genes[g].tss))[0]]
            # mingle miRIAD match
            elif len(m.miriad_match) == 1:
                host_gene = self.genes[list(m.miriad_match)[0]]
                flags['= 1 miRIAD match'].add(m.gene.id)
            # multiple miRIAD matches
            elif len(m.miriad_match) > 1:
                host_gene = self.genes[sorted(m.miriad_match, key=lambda g: m.gene.tss_dist(self.genes[g].tss))[0]]
                flags['> 1 miRIAD match'].add(m.gene.id)

            # update the gene promoter based on the host gene
            if host_gene:
                m.gene.change_promoter(host_gene.tss, host_gene.promoter_start, host_gene.promoter_end)
                m.gene.host_gene_id = host_gene.id
                m.gene.host_gene_name = host_gene.name
            else:
                flags['no host gene'].add(m.gene.id)

        # output matching results
        print('\tmatches')
        for k, v in sorted(flags.items()):
            p(k, v)

        print()
        # analyze genes with and without host genes
        p('has host gene', {gene.id for gene in self.genes.values() if gene.host_gene_id})
        p('host gene contained', {gene.id for gene in self.genes.values() if gene.host_gene_id and gene.contained_in
                                  and gene.host_gene_id in gene.contained_in})
        p('contained but no host gene', {gene.id for gene in self.genes.values() if gene.host_gene_id and
                                         gene.contained_in and gene.host_gene_id not in gene.contained_in})
        p('no host gene', {gene.id for gene in self.genes.values() if gene.type == 'miRNA' and not gene.host_gene_id})

    def _write_definitions(self, cfg: dict):
        """
        Write the processed miRNA promoter definitions to the specified file path.
        :param cfg: configuration dictionary containing file paths and relevant settings
        """
        with open(cfg['file paths']['reference']['genes'], 'w') as file:
            columns = ['gene.id', 'gene.name', 'mir.base.id', 'gene.type', 'gene.description', 'gene.chromosome',
                       'gene.strand', 'gene.start', 'gene.end', 'gene.length', 'gene.TSS', 'TSS.distance',
                       'promoter.start', 'promoter.end', 'promoter.length', 'promoter.overlap', 'transcript.count',
                       'host.gene.id', 'host.gene.name', 'contained.in', 'contains', 'overlaps', 'identical.range',
                       'GO.ids']
            file.write('\t'.join(columns) + '\n')

            for gene in sorted(self.genes.values(), key=lambda g: g.id):
                columns = [gene.id, gene.name, gene.mir_base_id, gene.type, gene.description,
                           gene.chromosome, gene.strand, gene.start, gene.end, gene.length, gene.tss, gene.tss_distance,
                           gene.promoter_start, gene.promoter_end, gene.promoter_length, gene.promoter_overlap,
                           gene.transcript_count, gene.host_gene_id, gene.host_gene_name,
                           ';'.join(sorted(gene.contained_in)), ';'.join(sorted(gene.contains)),
                           ';'.join(sorted(gene.overlaps)), ';'.join(sorted(gene.identical)),
                           ';'.join(sorted(gene.GO_ids))]
                file.write('\t'.join([str(c) for c in columns]) + '\n')
