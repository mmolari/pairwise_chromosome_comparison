# %%
import pandas as pd
import numpy as np
import pypangraph as pp
import pathlib
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--graph",
        type=str,
        help="Path to the graph file",
        required=True,
    )
    parser.add_argument(
        "--core_alignments",
        type=str,
        help="Path to the core alignments folder",
        required=True,
    )
    parser.add_argument(
        "--lengths",
        type=str,
        help="Path to the lengths file",
        required=True,
    )
    parser.add_argument(
        "--block_positions",
        type=str,
        help="Path to the block positions file",
        required=True,
    )
    parser.add_argument(
        "--out_csv",
        type=str,
        help="Path to the output csv file",
        required=True,
    )
    args = parser.parse_args()
    return args


def load_dfs(args):
    pan = pp.Pangraph.load_json(args.graph)

    fld_name = pathlib.Path(args.core_alignments)
    I = pd.read_csv(fld_name / "ins.csv")
    D = pd.read_csv(fld_name / "dels.csv")
    M = pd.read_csv(fld_name / "snps.csv")

    fname = args.lengths
    Ls = pd.read_csv(fname).set_index("id")["length"].to_dict()

    fname = args.block_positions
    P = pd.read_csv(fname)
    P.set_index(["genome", "block_id", "occurrence_number"], inplace=True)

    return pan, I, D, M, Ls, P


def aln_pos_to_seq_pos(aln_pos, ins, dels):
    seq_pos = aln_pos
    for ins_loc, ins_seq in ins:
        i0, i1 = ins_loc
        if i0 < aln_pos:
            seq_pos += len(ins_seq)

    for del_loc, del_len in dels:
        if del_loc < aln_pos:
            seq_pos -= del_len

    return seq_pos


def seq_pos_to_genome_pos(block_strand, block_start, block_end, seq_pos, L):
    if block_strand == True:
        return (block_start + seq_pos) % L
    else:
        return (block_end - seq_pos) % L


if __name__ == "__main__":
    args = parse_args()
    pan, I, D, M, Ls, P = load_dfs(args)

    res = []
    for X, tp in [(M, "snp"), (I, "ins"), (D, "del")]:
        for idx, row in X.iterrows():
            bid = row["block_id"]
            main_iso = row["iso"]
            main_occ = row["block_num"]
            consensus_pos = row["block_aln_pos"]
            aln = pan.blocks[bid].alignment

            for o in aln.occs:
                iso, num, b_strand = o

                b_strand, b_start, b_end = P.loc[(iso, bid, num)]

                seq_pos = aln_pos_to_seq_pos(consensus_pos, aln.ins[o], aln.dels[o])
                genome_pos = seq_pos_to_genome_pos(
                    b_strand, b_start, b_end, seq_pos, Ls[iso]
                )

                res.append(
                    {
                        "mut_idx": idx,
                        "mut_type": tp,
                        "genome": iso,
                        "block_id": bid,
                        "occurrence_number": num,
                        "strand": b_strand,
                        "genome_pos": genome_pos,
                        "block_aln_pos": consensus_pos,
                        "is_main": main_iso == iso and main_occ == num,
                    }
                )
    res = pd.DataFrame(res)
    res.to_csv(args.out_csv, index=False)
