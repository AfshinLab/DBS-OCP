"""
Find barcodes and merge originating from the same droplet/compartment.
"""
from collections import OrderedDict, defaultdict, Counter
import dataclasses
import logging
from typing import Dict, Set, Tuple

from xopen import xopen
import numpy as np
import scipy
import scipy.stats
import scipy.sparse

from dbsocp.utils import Summary, tqdm

logger = logging.getLogger(__name__)


# TODO Parallel processing of Tabix files with pysam
# - http://databio.org/posts/tabix_files.html
# - https://github.com/databio/pararead

# Constants
MIN_OVERLAPS = 1
MAX_BARCODES_PRECENTILE = 99


def add_arguments(parser):
    parser.add_argument("input", help="Coordinate-sorted Fragment file.")
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output Fragment file with merged barcode duplicates.",
    )
    parser.add_argument(
        "-m",
        "--merges",
        help="Output TSV with barcodes that were merged in form: <old_barcode> "
        "<new_barcode>",
    )
    parser.add_argument(
        "-p",
        "--plot-similarity",
        help="Output plot of ranked jaccard similarity for overlapping barcode pairs "
        "to file",
    )
    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        default=0.5,
        help="Jaccard index threshold in range 0->1 for merging barcodes. "
        "Default: %(default)s",
    )
    parser.add_argument(
        "-s",
        "--skip-contigs",
        help="Comma separated list of contings to skip for merging",
    )
    parser.add_argument(
        "--mode",
        choices=["fragment", "cutsite"],
        default="fragment",
        help="Compare overlaps based on 'fragment' (default) or 'cutsite'.",
    )


def main(args):
    contigs = set() if args.skip_contigs is None else set(",".split(args.skip_contigs))
    run_mergedups(
        input=args.input,
        output=args.output,
        merges=args.merges,
        plot_similarity=args.plot_similarity,
        threshold=args.threshold,
        skip_contigs=contigs,
        mode=args.mode,
    )


def run_mergedups(
    input: str,
    output: str,
    merges: str,
    plot_similarity: str,
    threshold: float,
    skip_contigs: Set[str],
    mode: str,
):
    logger.info("Starting Analysis")
    summary = Summary()

    # matrix
    #                              barcodes
    #                           Bc1 Bc2 ... BcN
    #                       --------------------
    #               coord1  |   1   1   ... 0   |
    #  coordinates  coord2  |   0   1   ... 1   |
    #               coord3  |   0   0   ... 1   |
    #               coord4  |   1   1   ... 0   |
    #                       --------------------
    if mode == "fragment":
        matrix, index_barcode, barcode_counts = generate_frag_matrix(
            input, skip_contigs, summary
        )
    elif mode == "cutsite":
        matrix, index_barcode, barcode_counts = generate_cutsite_matrix(
            input, skip_contigs, summary
        )
    else:
        raise ValueError(f"Unknown mode '{mode}'.")

    # Remove coordinates that are covered by more than the MAX_BARCODES_PRECENTILE
    # percentile in number of barcodes
    bcs_per_coord = matrix.sum(axis=1)
    max_bcs = int(np.percentile(bcs_per_coord, MAX_BARCODES_PRECENTILE, axis=0))
    max_bcs = max(max_bcs, 2)
    logger.info(f"Filtering coordinate with more than {max_bcs} barcodes.")
    matrix = matrix[np.ravel(bcs_per_coord) < max_bcs, :]

    logger.info("Find overlapping barcodes")
    # Get overlapping barcodes
    # overlaps
    #                              barcodes
    #                           Bc1 Bc2 ... BcN
    #                        --------------------
    #               Bc1     |   2   2   ... 0   |
    #  barcodes     Bc2     |   2   3   ... 1   |
    #               ...     |   .   .   ... .   |
    #               BcN     |   0   1   ... 2   |
    #                        --------------------
    overlapps = matrix.transpose() * matrix

    # Only consider barcodes with more than the minimun required overlapping position.
    overlapping_bcs = overlapps > MIN_OVERLAPS

    # Get coordinates for overlaps in lower triagular along with values which gives the
    # number of overlaps
    bcs_rows, bcs_cols = overlapping_bcs.nonzero()
    del overlapping_bcs
    in_lower = bcs_rows < bcs_cols
    bcs_rows = bcs_rows[in_lower]
    bcs_cols = bcs_cols[in_lower]
    nr_overlapps = np.array(overlapps[bcs_rows, bcs_cols]).flatten()
    del overlapps

    summary["Overlapping Barcodes"] = len(nr_overlapps)

    uf, jaccard_similarity = call_duplicates(
        nr_overlapps,
        bcs_rows,
        bcs_cols,
        index_barcode,
        barcode_counts,
        threshold,
        summary,
    )

    if plot_similarity is not None and len(jaccard_similarity) > 0:
        import matplotlib.pyplot as plt

        # Remove perfect matched barcodes
        jaccard_similarity = jaccard_similarity[jaccard_similarity < 1.0]
        threshold_index = (np.abs(jaccard_similarity - threshold)).argmin()

        # Remove long tail
        jaccard_cumsum = jaccard_similarity.cumsum() / jaccard_similarity.sum()
        y = jaccard_similarity[jaccard_cumsum < 0.99]

        plt.plot(y, color="k", label="Jaccard similarity")
        plt.axvline(threshold_index, 0, 1, color="r", alpha=0.8, label="Threshold")
        plt.xlabel("Barcode pair rank")
        plt.ylabel("Jaccard similarity")
        plt.legend(loc="upper right")
        plt.savefig(plot_similarity)

    if merges is not None:
        logger.info(f"Writing merged barcodes to {merges}.")
        with open(merges, "w") as outfile:
            for component in uf.connected_components():
                for barcode in component:
                    if barcode != uf[barcode]:
                        print(barcode, uf[barcode], sep="\t", file=outfile)

    logger.info(f"Writing updated fragments to {output}.")
    write_merged_fragments(input, output, uf, summary)

    logger.info("Finished")
    summary.print_stats(name=__name__)


def write_merged_fragments(input: str, output: str, uf: "UnionFind", summary: Dict[str, int]):
    """Write merged fragments to file."""
    with open(output, "w") as outfile:
        parser = parse_fragment_file(input)
        prev_fragment = next(parser)
        for fragment in tqdm(parser, desc="Update fragments", initial=1):
            fragment.barcode = uf[fragment.barcode]

            if fragment == prev_fragment:
                prev_fragment.update(fragment)
            else:
                summary["Fragments written"] += 1
                print(prev_fragment, file=outfile)
                prev_fragment = fragment

        summary["Fragments written"] += 1
        print(prev_fragment, file=outfile)


def generate_cutsite_matrix(file: str, skip_contigs: Set[str], summary: Dict[str, int]):
    """Generate cutsite vs barcode matrix from fragment file."""
    # Contants
    BUFFER = 1000
    DOUBLE_BUFFER = 2000

    # Containers for generating matrix
    indices = []
    indptr = [0]
    barcode_index = OrderedDict()
    index_barcode = {}

    barcode_nsites = Counter()

    logger.info(f"Reading fragments from {file}")
    parser = parse_fragment_file(file)

    # Parse first frag
    first_frag = next(parser)
    summary["Fragments read"] += 1
    cutsite1, cutsite2 = first_frag.get_cutsites()

    sites_cache = defaultdict(set)
    sites_cache[cutsite1].add(first_frag.barcode)
    sites_cache[cutsite2].add(first_frag.barcode)
    next_checkpoint = cutsite1.position + DOUBLE_BUFFER
    prev_chromosome = cutsite1.chromosome

    # Parse rest
    for fragment in tqdm(parser, desc="Parsing fragments", initial=1):
        if fragment.chromosome in skip_contigs:
            continue

        cutsite1, cutsite2 = fragment.get_cutsites()
        summary["Fragments read"] += 1
        barcode = fragment.barcode
        barcode_nsites[barcode] += 2

        if cutsite1.chromosome != prev_chromosome:
            for site, barcodes in sites_cache.items():
                if len(barcodes) < 2:
                    continue

                summary["Cutsites duplicate"] += 1
                update_matrix_data(
                    barcodes, barcode_index, index_barcode, indices, indptr
                )

            summary["Cutsites"] += len(sites_cache)
            sites_cache.clear()
            next_checkpoint = cutsite1.position + DOUBLE_BUFFER
            prev_chromosome = cutsite1.chromosome

        elif cutsite1.position > next_checkpoint:
            for site in sorted(sites_cache):
                if cutsite1.position - site.position < BUFFER:
                    break

                summary["Cutsites"] += 1
                barcodes = sites_cache.pop(site)

                if len(barcodes) < 2:
                    continue

                summary["Cutsites duplicate"] += 1
                update_matrix_data(
                    barcodes, barcode_index, index_barcode, indices, indptr
                )

            next_checkpoint = cutsite1.position + DOUBLE_BUFFER

        sites_cache[cutsite1].add(barcode)
        sites_cache[cutsite2].add(barcode)

    for site, barcodes in sites_cache.items():
        if len(barcodes) < 2:
            continue

        summary["Cutsites duplicate"] += 1
        update_matrix_data(barcodes, barcode_index, index_barcode, indices, indptr)

    summary["Cutsites"] += len(sites_cache)
    sites_cache.clear()

    summary["Barcodes reads"] = len(barcode_nsites)

    logger.info("Generating Barcode vs. Fragment matrix")
    data = np.ones(len(indices))
    matrix = scipy.sparse.csr_matrix((data, indices, indptr), dtype=int)
    return matrix, index_barcode, barcode_nsites


def generate_frag_matrix(file: str, skip_contigs: Set[str], summary: Dict[str, int]):
    barcode_counts = Counter()

    # Containers for generating matrix
    indices = []
    indptr = [0]
    barcode_index = OrderedDict()
    index_barcode = {}

    prev_fragment = Fragment(None, -1, -1, None, -1)
    prev_barcodes = set()
    prev_dup = False
    logger.info(f"Reading fragments from {file}")
    for fragment in tqdm(parse_fragment_file(file), desc="Parsing fragments"):
        if fragment.chromosome in skip_contigs:
            continue

        summary["Fragments read"] += 1
        barcode = fragment.barcode
        barcode_counts[barcode] += 1

        if fragment.match_coordinates(prev_fragment):
            prev_barcodes.add(barcode)
            prev_dup = True
            continue

        if prev_dup:
            summary["Fragments duplicate"] += 1
            update_matrix_data(
                prev_barcodes, barcode_index, index_barcode, indices, indptr
            )
            prev_dup = False

        prev_fragment = fragment
        prev_barcodes = {barcode}

    if prev_dup:
        summary["Fragments duplicate"] += 1
        update_matrix_data(prev_barcodes, barcode_index, index_barcode, indices, indptr)

    summary["Barcodes reads"] = len(barcode_counts)

    logger.info("Generating Barcode vs. Fragment matrix")
    data = np.ones(len(indices))
    matrix = scipy.sparse.csr_matrix((data, indices, indptr), dtype=int)
    return matrix, index_barcode, barcode_counts


def call_duplicates(
    nr_overlapps, bcs_rows, bcs_cols, index_barcode, barcode_counts, threshold, summary
) -> "UnionFind":
    """Iterate over overlapping positions and generate duplicate calls which are used
    create a UnionFind object"""
    uf = UnionFind()
    jaccard_similarity = []
    for nr_shared, i1, i2 in tqdm(
        zip(nr_overlapps, bcs_rows, bcs_cols), desc="Find overlaps", total=len(bcs_rows)
    ):
        bc1 = index_barcode[i1]
        bc2 = index_barcode[i2]
        total = barcode_counts[bc1] + barcode_counts[bc2] - nr_shared
        jaccard_index = nr_shared / total
        jaccard_similarity.append(jaccard_index)
        logger.debug("Overlapping pair: {} {} {}".format(bc1, bc2, jaccard_index))
        if jaccard_index > threshold:
            summary["Barcodes merged"] += 1
            uf.union(bc1, bc2)
    return uf, np.array(sorted(jaccard_similarity, reverse=True))


def update_matrix_data(coord_barcodes, barcode_index, index_barcode, indices, indptr):
    for barcode in coord_barcodes:
        index = barcode_index.setdefault(barcode, len(barcode_index))
        index_barcode[index] = barcode
        indices.append(index)
    indptr.append(len(indices))


def parse_fragment_file(file: str):
    with xopen(file) as f:
        for line in f:
            chromosome, start, end, barcode, count, *_ = line.strip().split("\t")
            yield Fragment(chromosome, int(start), int(end), barcode, int(count))


@dataclasses.dataclass(eq=False)
class Fragment:
    chromosome: str
    start: int
    end: int
    barcode: str
    count: int
    __slots__ = ["chromosome", "start", "end", "barcode", "count"]

    def update(self, other: "Fragment"):
        self.count += other.count

    def match_coordinates(self, other):
        return (self.chromosome, self.start, self.end) == (
            other.chromosome,
            other.start,
            other.end,
        )

    def __eq__(self, other) -> bool:
        return (self.chromosome, self.start, self.end, self.barcode) == (
            other.chromosome,
            other.start,
            other.end,
            other.barcode,
        )

    def __str__(self):
        return (
            f"{self.chromosome}\t{self.start}\t{self.end}\t{self.barcode}\t"
            f"{self.count}"
        )

    def get_cutsites(self) -> Tuple["CutSite", "CutSite"]:
        return CutSite(self.chromosome, self.start), CutSite(self.chromosome, self.end)


@dataclasses.dataclass(frozen=True)
class CutSite:
    chromosome: str
    position: int
    __slots__ = ["chromosome", "position"]

    def __lt__(self, other) -> bool:
        return self.position < other.position


class UnionFind:
    """Union-find data structure.
    Each UnionFind instance X maintains a family of disjoint sets of
    hashable objects, supporting the following two methods:
    - X[item] returns a name for the set containing the given item.
      Each set is named by an arbitrarily-chosen one of its members; as
      long as the set remains unchanged it will keep the same name. If
      the item is not yet part of a set in X, a new singleton set is
      created for it.
    - X.union(item1, item2, ...) merges the sets containing each item
      into a single larger set.  If any item is not yet part of a set
      in X, it is added to X as one of the members of the merged set.

    Based on Josiah Carlson's code,
    http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/215912
    with significant additional changes by D. Eppstein.
    https://www.ics.uci.edu/~eppstein/PADS/UnionFind.py
    """

    def __init__(self, mapping=None):
        """Create a new  union-find structure."""
        self.parents = mapping if isinstance(mapping, dict) else {}

    def __getitem__(self, object: str) -> str:
        """Find and return the name of the set containing the object."""

        # check for previously unknown object
        if object not in self.parents:
            self.parents[object] = object
            return object

        # find path of objects leading to the root
        path = [object]
        root = self.parents[object]
        while root != path[-1]:
            path.append(root)
            root = self.parents[root]

        # compress the path and return
        for ancestor in path:
            self.parents[ancestor] = root
        return root

    def __contains__(self, item: str):
        return item in self.parents

    def __iter__(self):
        """Iterate through all items ever found or unioned by this structure."""
        return iter(self.parents)

    def items(self):
        """Iterate over tuples of items and their root"""
        for x in self:
            yield x, self[x]

    def union(self, *objects):
        """Find the sets containing the objects and merge them all."""
        roots = [self[x] for x in objects]

        # Use lexicographical ordering to set main root in set
        heaviest = sorted(roots)[0]
        for r in roots:
            if r != heaviest:
                self.parents[r] = heaviest

    def connected_components(self):
        """Iterator for sets"""
        components = defaultdict(list)
        for item, root in self.items():
            components[root].append(item)

        for component in components.values():
            yield component

    def same_component(self, *objects) -> bool:
        """Returns true if all objects are present in the same set"""
        if all(x in self for x in objects):
            return len({self[x] for x in objects}) == 1
        return False

    def update(self, other: "UnionFind"):
        """Update sets based on other UnionFind instance"""
        for x, root in other.items():
            self.union(x, root)

    @classmethod
    def from_dict(cls, mapping: Dict[str, str]):
        return cls(mapping.copy())
