import util


class Closest:
    """
    Represents the closest feature to a genomic bin, including distance and associated pairs.
    """
    def __init__(self):
        # set of closest feature pairs
        self.pairs = set()
        # distance to the closest feature
        self.distance = 0


class BinCount:
    """
    Represents bins of genomic regions, their overlap status, and distances to the closest reference features.
    """
    def __init__(self, prefix: str, ref_path: str, bins: dict[tuple[str, int, int], str]):
        # path to save bins that do not overlap with reference features
        self.non_overlap_path = prefix + '_not_overlapping.bed'             # type: str
        # path to save closest features for non-overlapping bins
        self.closest_path = prefix + '_closest.bed'                         # type: str
        # path to the reference file
        self.ref_path = ref_path                                            # type: str

        # dictionary mapping bins to identifiers
        self.all = bins                 # type: dict[tuple[str, int, int], str]

        # set of bins that overlap with reference features
        self.overlap = set()            # type: set[tuple[str, int, int]]
        # set of bins that do not overlap with reference features
        self.not_overlapping = set()    # type: set[tuple[str, int, int]]

        # mapping of bins to Closest objects for non-overlapping bins
        self.distances = dict()         # type: dict[tuple[str, int, int], Closest]

    def add_overlap(self, overlap: set[tuple[str, int, int]], ref_path: str):
        """
        Updates overlap information and calculates closest features for non-overlapping bins.
        :param overlap: bins that overlap with reference features
        :param ref_path: path to the reference file
        """
        self.overlap = overlap
        self.not_overlapping = {key for key in self.all.keys() if key not in overlap}

        # initialize Closest objects for non-overlapping bins
        self.distances = {key: Closest() for key in self.not_overlapping}

        # write non-overlapping bins to file
        with open(self.non_overlap_path, 'w') as file:
            for bin in sorted(self.not_overlapping):
                file.write('\t'.join([bin[0], str(bin[1]), str(bin[2]), self.all[bin]]) + '\n')

        # calculate closest features using a utility function
        util.closest(ref_path, self.non_overlap_path, self.closest_path)

        # populate Closest objects with distances and feature pairs
        with open(self.closest_path, 'r') as file:
            for line in file:
                cols = line.strip().split('\t')
                bin = (cols[0], int(cols[1]), int(cols[2]))

                self.distances[bin].pairs.add(tuple(cols[7].split('/')))
                self.distances[bin].distance = int(cols[8])

    def promoter_distances(self) -> list[int]:
        """Returns distances for bins associated with promoters."""
        return [c.distance for c in self.distances.values() if 'promoter' in {x for _, x in c.pairs}]

    def gene_distances(self) -> list[int]:
        """Returns distances for bins associated with genes."""
        return [c.distance for c in self.distances.values() if 'gene' in {x for _, x in c.pairs}]

    def distances(self) -> list[int]:
        """Returns distances for all bins with associated feature pairs."""
        return [c.distance for c in self.distances.values() if c.pairs]


class GeneExpression:
    """ Represents the expression characteristics of a gene. """
    def __init__(self, gene_type: str, fc: str, p: str, q: str, t: float):
        """
        :param fc: log2 fold change
        :param q: FDR corrected p-value
        :param t: p-value threshold for differential expression
        """
        self.type = gene_type                                                                           # type: str
        self.fc = float(fc)                                                                             # type: float
        # p-value
        self.p = float(p)                                                                               # type: float
        # Benjamini-Hochberg adjusted p-value (q-value)
        self.q = float(q)                                                                               # type: float
        # True if the gene is differentially expressed, False otherwise
        self.is_diff = self.q < t                                                                       # type: bool

        # categorise the expression as either non-differentially expressed, up- or down-regulated
        if not self.is_diff:
            self.cat = 0                                                                                # type: int
        elif self.fc < 0:
            self.cat = -1
        else:
            self.cat = 1

        # categorise expression just by looking at the fold change
        if self.fc < 0:
            self.fc_cat = -1                                                                            # type: int
        elif self.fc > 0:
            self.fc_cat = 1
        else:
            self.fc_cat = 0

        # set the expression of non-differentially expressed genes to 0 and keep the fold change for the rest
        self.fc_diff = self.fc if self.is_diff else 0.0                                                 # type: float


class ModCounter:
    """ Tracks the number of modifications in a transition from cell state A to cell state B. """
    def __init__(self):
        # A = number of bins differentially modified in cell A
        self.A = 0                          # type: int
        # B = number of bins differentially modified in cell B
        self.B = 0                          # type: int

    def unmod(self) -> bool:
        """ Returns True if no modifications are present in both cell states. """
        return self.A == 0 and self.B == 0

    def both(self) -> bool:
        """ Returns True if modifications are present in both cell states. """
        return self.A > 0 and self.B > 0

    def diff(self) -> int:
        """ Returns -1 if only A is modified, 1 if only B is modified, 0 otherwise. """
        if self.unmod() or self.both():
            return 0

        return -1 if self.A else 1

    def diff_2(self) -> int:
        """ Returns -1 if only A is more modified than, 1 if B more than A, 0 otherwise. """
        if self.A > self.B:
            return -1

        if self.B > self.A:
            return 1

        return 0

    def cat(self) -> int:
        """
        Categorizes modifications (0: no modifications, 1: only A modified, 2: only B modified, 3: both modified).
        """
        if self.unmod():
            return 0
        if self.both():
            return 3
        if self.A > self.B:
            return 1
        return 2

    def cat_2(self) -> int:
        """ Categorizes modifications (0: no modifications, 1: A > B, 2: B > A, 3: both equally modified). """
        if self.unmod():
            return 0
        if self.both and self.A == self.B:
            return 3
        return 1 if self.A > self.B else 2


class ModIndex:
    """ Index of modification information for histones and genes. """
    def __init__(self):
        # {histone: {gene ID: modification information}}
        self.mods = dict()                                                      # type: dict[str, dict[str, ModCounter]]

    def add_histone(self, histone: str, counter_dict: dict[str, ModCounter]):
        """ Adds modification data for a specific histone. """
        self.mods[histone] = counter_dict

    def unmod(self, histone: str) -> set[str]:
        """ Returns genes with no modifications for the given histone. """
        return {gene_id for gene_id, c in self.mods[histone].items() if c.unmod()}

    def mod(self, histone: str) -> set[str]:
        """ Returns genes with modifications for the given histone. """
        return {gene_id for gene_id, c in self.mods[histone].items() if c.A or c.B}

    def A_only(self, histone: str) -> set[str]:
        """ Returns genes modified only in state A for the given histone. """
        return {gene_id for gene_id, c in self.mods[histone].items() if c.A and not c.B}

    def B_only(self, histone: str) -> set[str]:
        """ Returns genes modified only in state B for the given histone. """
        return {gene_id for gene_id, c in self.mods[histone].items() if c.B and not c.A}

    def A_more(self, histone: str) -> set[str]:
        """ Returns genes with more modifications in state A for the given histone. """
        return {gene_id for gene_id, c in self.mods[histone].items() if c.A > c.B}

    def B_more(self, histone: str) -> set[str]:
        """ Returns genes with more modifications in state B for the given histone. """
        return {gene_id for gene_id, c in self.mods[histone].items() if c.B > c.A}

    def both(self, histone: str) -> set[str]:
        """ Returns genes modified in both states for the given histone. """
        return {gene_id for gene_id, c in self.mods[histone].items() if c.A and c.B}

    def both_equal(self, histone: str) -> set[str]:
        """ Returns genes equally modified in both states for the given histone. """
        return {gene_id for gene_id, c in self.mods[histone].items() if c.A and c.A == c.B}


class GeneFunctionIndex:
    """
    Tracks gene functional categories, histone modifications, and expression changes.
    """
    def __init__(self, cell_1, cell_2, prefix: str, ref_path: str, exp: dict[str, GeneExpression]):
        self.cell_1 = cell_1
        self.cell_2 = cell_2
        self.bin_prefix = prefix
        self.ref_path = ref_path

        # gene expression information
        self.exp = exp                                                              # type: dict[str, GeneExpression]

        self.indices = {k: ModIndex() for k in ['promoter', 'gene', 'combined']}    # type: dict[str, ModIndex]
        self.bins = dict()                                                          # type: dict[str, BinCount]

    def add_modification(self, histone: str, indices: dict[str, dict[str, ModCounter]]):
        """
        Adds modification data for a histone.
        :param histone: histone identifier
        :param indices: modification data categorized by gene location
        """
        for location, counter_dict in indices.items():
            counter_dict = {k: counter_dict[k] if k in counter_dict.keys() else ModCounter() for k in self.exp.keys()}
            self.indices[location].add_histone(histone, counter_dict)

    def add_bin_count(self, histone: str, bin_count: dict[tuple[str, int, int], str]):
        """
        Adds bin data for a histone.
        :param histone: histone identifier
        :param bin_count: mapping of bins to identifiers
        """
        self.bins[histone] = BinCount(self.bin_prefix + '_' + histone, self.ref_path, bin_count)

    def add_closest_bins(self, histone: str, overlaps: set[tuple[str, int, int]]):
        """
        Updates overlap and closest feature data for a histone.
        :param histone: histone identifier
        :param overlaps: bins that overlap with reference features#
        """
        self.bins[histone].add_overlap(overlaps, self.ref_path)
