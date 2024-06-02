import pypangraph as pp
import pandas as pd
from Bio import SeqIO, SeqRecord, Seq
import argparse
import pathlib


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--graph",
        type=str,
        help="Path to the graph file",
        required=True,
    )
    parser.add_argument(
        "--out_fld",
        type=str,
        help="Path to the output folder",
        required=True,
    )
    args = parser.parse_args()
    return args


def create_alignment(aln):
    A, O = aln.generate_alignments()
    records = []
    for a, o in zip(A, O):
        seq = Seq.Seq(a)
        idx, num, strand = o
        record = SeqRecord.SeqRecord(seq, id=idx, description=f"{num} {strand}")
        records.append(record)
    return records


def extract_variations(aln):
    M, I, D = [], [], []

    for k, ms in aln.muts.items():
        iso, num, strand = k
        for m in ms:
            pos, alt = m
            M.append((cb, iso, num, strand, pos, alt))

    for k, ii in aln.ins.items():
        iso, num, strand = k
        for i in ii:
            P, seq = i
            pos, ins_start = P
            L = len(seq)
            I.append((cb, iso, num, strand, pos, ins_start, L, seq))

    for k, ds in aln.dels.items():
        iso, num, strand = k
        for d in ds:
            d_start, d_len = d
            D.append((cb, iso, num, strand, d_start, d_len))

    return M, I, D


if __name__ == "__main__":

    args = parse_args()

    pan = pp.Pangraph.load_json(args.graph)

    aln_fld = pathlib.Path(args.out_fld)
    aln_fld.mkdir(exist_ok=True, parents=True)

    # select core blocks
    bdf = pan.to_blockstats_df()
    core_mask = bdf["core"]
    core_blocks = bdf[core_mask].index.to_numpy()

    # for every core block create and save the alignment
    # and add mutations to dataframes
    corealn_fld = aln_fld / "core_alignments"
    corealn_fld.mkdir(exist_ok=True, parents=True)
    snps, ins, dels = [], [], []
    for cb in core_blocks:
        block = pan.blocks[cb]
        aln = block.alignment

        # create and save alignment
        aln_records = create_alignment(aln)
        aln_file = corealn_fld / f"{cb}.fa"
        SeqIO.write(aln_records, aln_file, "fasta")

        # record mutations
        M, I, D = extract_variations(aln)
        snps += M
        ins += I
        dels += D

    # finalize and save mutation dataframes
    snps = pd.DataFrame(
        snps,
        columns=[
            "block_id",
            "iso",
            "block_num",
            "strand",
            "block_aln_pos",
            "alt",
        ],
    )
    ins = pd.DataFrame(
        ins,
        columns=[
            "block_id",
            "iso",
            "block_num",
            "strand",
            "block_aln_pos",
            "ins_start",
            "ins_len",
            "seq",
        ],
    )
    dels = pd.DataFrame(
        dels,
        columns=[
            "block_id",
            "iso",
            "block_num",
            "strand",
            "block_aln_pos",
            "del_len",
        ],
    )

    snps.to_csv(aln_fld / "snps.csv")
    ins.to_csv(aln_fld / "ins.csv")
    dels.to_csv(aln_fld / "dels.csv")
