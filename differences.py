from Bio import pairwise2
from typing import List, Tuple


# deletions are returned in vcf style i. e. number of residue right before deletion "-" as reference and insertion as
# alt
def diff(ref: str, alt: str, thr: int = 20) -> List[Tuple[int, int, str, str]]:
    aln = pairwise2.align.globalxx(ref, alt)[0]

    if aln.score < thr:
        return []

    aln_i: int = 0
    seq_i: int = 1  # bed-like indexing starting from 1

    out: List[Tuple[int, int, str, str]] = []

    alt_gap: bool = False

    insertion: str = ""

    alt_gap_start: int = 0

    ref_gap: bool = False
    ref_gap_start: int = 0

    seqA: str = aln.seqA
    for res1 in seqA:
        res2: str = aln.seqB[aln_i]

        if res1 != "-" and ref_gap:
            print(insertion)
            out.append((ref_gap_start, ref_gap_start, "-", insertion))
            ref_gap = False
            insertion = ""

        if res2 != "-" and alt_gap:
            out.append((alt_gap_start, seq_i - 1, insertion, "-" * len(insertion)))
            alt_gap = False
            insertion = ""

        if res1 != "-" and res2 != "-":
            if res1 != res2:
                out.append((seq_i, seq_i, res1, res2))

            seq_i += 1

        if res1 != "-" and res2 == "-":
            if not alt_gap:
                alt_gap_start = seq_i

            alt_gap = True
            insertion = insertion + res1
            seq_i += 1

        if res1 == "-" and res2 != "-":
            if not ref_gap:
                ref_gap_start = seq_i - 1

            ref_gap = True
            insertion = insertion + res2

        aln_i += 1

    return out


print(diff("ATG", "CTG"))