# %%
import numpy as np
import pandas as pd
import pypangraph as pp
from Bio import SeqIO, Seq, SeqRecord
import pathlib
import matplotlib.pyplot as plt
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--graph", type=str, help="Pangraph JSON file")
    parser.add_argument("--seq_lengths", type=str, help="Sequence lengths CSV file")
    parser.add_argument("--msu", type=str, help="Minimal synteny units CSV file")
    parser.add_argument("--out_aln_fld", type=str, help="Output alignment folder")
    parser.add_argument("--out_plot", type=str, help="Output plot file")
    parser.add_argument("--out_info", type=str, help="Output info file")
    return parser.parse_args()


def load_args():
    args = parse_args()
    pan = pp.Pangraph.load_json(args.graph)
    seq_lengths = pd.read_csv(args.seq_lengths).set_index("id")["length"].to_dict()
    msu = pd.read_csv(args.msu)
    check_signatures(msu)
    k1, k2 = pan.strains()
    mask = (msu["path"] == k2) & (msu["msu"] != 0)
    sign_dict = msu[mask].set_index("signature")
    msu.set_index(["path", "bid", "strand", "occ"], inplace=True)
    msu_dict = msu["msu"].to_dict()

    return pan, seq_lengths, msu_dict, sign_dict, msu


def check_signatures(msu):
    # check that all signatures are different
    mask = msu["msu"] != 0
    sdf = msu[mask]
    N = len(sdf)
    M = len(sdf["signature"].unique())
    assert N / 2 == M, f"Found {N/2} MSUs but only {M} unique signatures"


def extract_pathinfo(pan, path):
    blocks = pan.paths[path].block_ids
    strands = pan.paths[path].block_strands
    nums = pan.paths[path].block_nums
    return blocks, strands, nums


def msu_extremes(B, S, O, k, m, msu_dict):
    # get start and end indices of an msu in a path
    s, e = None, None
    i, N = 0, len(B)
    while (s is None) or (e is None):
        idx = (k, B[i], S[i], O[i])
        j = (i + 1) % N
        next_idx = (k, B[j], S[j], O[j])
        mi = msu_dict[idx]
        mj = msu_dict[next_idx]

        if (mi == m) and (mj != m):
            e = i
        elif (mi != m) and (mj == m):
            s = j
        i = j
    return s, e


def aln_matrix(aln1, aln2):
    align_matrix = np.array([list(aln1), list(aln2)])
    mask = np.any(align_matrix != "-", axis=0)
    align_matrix = align_matrix[:, mask]
    return align_matrix


def save_aln(A, m, k1, k2, svfld):
    svfld = pathlib.Path(svfld)
    svfld.mkdir(exist_ok=True, parents=True)
    aln1 = "".join(A[0, :])
    aln2 = "".join(A[1, :])
    rec1 = SeqRecord.SeqRecord(seq=Seq.Seq(aln1), id=f"MSU_{m}_{k1}")
    rec2 = SeqRecord.SeqRecord(seq=Seq.Seq(aln2), id=f"MSU_{m}_{k2}")
    SeqIO.write([rec1, rec2], svfld / f"MSU_{m}.fa", "fasta")


def SNPs(A):
    mask = np.all(A != "-", axis=0)
    mask &= A[0] != A[1]
    snp_idxs = np.where(mask)[0]
    return snp_idxs


def Dels(A):
    mask = A[0] == "-"
    assert np.all(A[1, mask] != "-")
    del_idxs = np.where(mask)[0]
    return del_idxs


def Ins(A):
    mask = A[1] == "-"
    assert np.all(A[0, mask] != "-")
    ins_idxs = np.where(mask)[0]
    return ins_idxs


def plot(As, start_pos, L, k1):

    kwargs = dict(
        histtype="step",
        cumulative=True,
    )

    ms = sorted(As)
    N = len(ms)
    fig, axs = plt.subplots(N, 1, figsize=(10, 2 * N))
    for y, m in enumerate(ms):
        ax = axs[y]
        A = As[m]
        bins = np.linspace(0, A.shape[1] + 2, 1000)
        snps, ins, dels = SNPs(A), Ins(A), Dels(A)
        ax.hist(snps, bins=bins, color="C0", label="SNPs", **kwargs)
        ax.hist(ins, bins=bins, color="C1", label="ins", **kwargs)
        ax.hist(dels, bins=bins, color="C2", label="dels", **kwargs)
        xticks = ax.get_xticks()
        xlabels = [f"{(int(x) + start_pos[m]) % L}" for x in xticks]
        ax.set_xticks(xticks)
        ax.set_xticklabels(xlabels)
        if start_pos[m] + A.shape[1] > L:
            ax.axvline(L - start_pos[m], color="k", lw=1, ls="--")
            ax.text(
                L - start_pos[m] + 1,
                0.0,
                f"genome start",
                ha="left",
                va="bottom",
                rotation=90,
                color="k",
            )
        ax.set_xlim(0, A.shape[1] + 2)
        ax.set_ylim(bottom=0)
        ax.set_ylabel(f"MSU {m} L={A.shape[1]/1000:.0f} kb")
        # despine
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    axs[0].legend()
    axs[-1].set_xlabel(f"Position {k1} (bp)")
    plt.tight_layout()
    return fig, axs


def MSUs_info(As):
    df = []
    for m in sorted(As):
        A = As[m]
        snps, ins, dels = SNPs(A), Ins(A), Dels(A)
        df.append(
            {
                "msu": m,
                "snps": len(snps),
                "ins": len(ins),
                "dels": len(dels),
                "len_aln": A.shape[1],
            }
        )
    return pd.DataFrame(df)


def extract_alns(pan, k1, k2, m, msu_dict, sign_dict, msu):
    B1, S1, O1 = extract_pathinfo(pan, k1)
    s, e = msu_extremes(B1, S1, O1, k1, m, msu_dict)
    i = s
    aln1, aln2 = "", ""
    while i != ((e + 1) % len(B1)):
        b1, s1, o1 = B1[i], S1[i], O1[i]
        idx1 = (k1, b1, s1, o1)
        sig = msu.loc[idx1, "signature"]
        idx2 = tuple(sign_dict.loc[sig][["path", "bid", "strand", "occ"]])
        b2, s2, o2 = idx2[1:]
        assert msu_dict[idx1] == m, f"Expected {m} but got {msu_dict[idx1]}"
        assert msu_dict[idx2] == m, f"Expected {m} but got {msu_dict[idx2]}"

        aln = pan.blocks[b1].alignment
        aln = dict(zip(*aln.generate_alignments()[::-1]))
        seq1 = Seq.Seq(aln[(k1, o1, s1)])
        seq2 = Seq.Seq(aln[(k2, o2, s2)])

        if not s1:
            seq1 = seq1.reverse_complement()
            seq2 = seq2.reverse_complement()

        aln1 += str(seq1)
        aln2 += str(seq2)

        i = (i + 1) % len(B1)
    return aln1, aln2


if __name__ == "__main__":

    args = parse_args()
    pan, Ls, msu_dict, sign_dict, msu = load_args()

    k1, k2 = pan.strains()
    B1, S1, O1 = extract_pathinfo(pan, k1)
    B2, S2, O2 = extract_pathinfo(pan, k2)

    As, start_pos = {}, {}
    msus = set(msu["msu"].unique()) - {0}
    for m in sorted(msus):
        aln1, aln2 = extract_alns(pan, k1, k2, m, msu_dict, sign_dict, msu)

        A = aln_matrix(aln1, aln2)
        save_aln(A, m, k1, k2, args.out_aln_fld)
        As[m] = A

        s, e = msu_extremes(B1, S1, O1, k1, m, msu_dict)
        start_pos[m] = pan.paths[k1].block_positions[s]

    fig, axs = plot(As, start_pos, Ls[k1], k1)
    fig.savefig(args.out_plot)
    plt.close(fig)

    df = MSUs_info(As)
    df.to_csv(args.out_info, index=False)
