# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pypangraph as pp
import segment_utils as su
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--graph", type=str, help="Pangraph JSON file")
    parser.add_argument("--seq_lengths", type=str, help="Sequence lengths CSV file")
    parser.add_argument("--msu", type=str, help="Minimal synteny units CSV file")
    parser.add_argument("--out", type=str, help="Output PDF file")
    return parser.parse_args()


def load_data(args):
    pan = pp.Pangraph.load_json(args.graph)
    seq_lengths = pd.read_csv(args.seq_lengths).set_index("id")["length"].to_dict()
    msu = pd.read_csv(args.msu)
    msu_dict = msu.set_index(["path", "bid", "strand", "occ"])["msu"].to_dict()
    sign_dict = msu.set_index(["path", "bid", "strand", "occ"])["signature"].to_dict()
    return pan, seq_lengths, msu_dict, sign_dict


def block_positions(pan):
    block_pos = {iso: {} for iso in pan.strains()}
    for iso in pan.strains():
        path = pan.paths[iso]
        B = path.block_ids
        S = path.block_strands
        O = path.block_nums
        P = path.block_positions
        N = len(B)
        for i, (b, s, o) in enumerate(zip(B, S, O)):
            block_pos[iso][(b, s, o)] = (P[i], P[(i + 1) % N])
    return block_pos


def pbc_plot(seg, c, orient, ax):
    if seg.x_runover():
        for s in seg.split_x():
            pbc_plot(s, c, orient, ax)
    elif seg.y_runover():
        for s in seg.split_y():
            pbc_plot(s, c, orient, ax)
    else:
        if seg.e1 == 0:
            seg.e1 = seg.L1
        if seg.e2 == 0:
            seg.e2 = seg.L2
        ax.plot(seg.x(), seg.y(), color=c, lw=1)


def pbc_plot_priv(seg, c, ax, kind):
    if seg.s > seg.e:
        for s in seg.split_x():
            pbc_plot_priv(s, c, ax, kind)
    else:
        x = seg.x() if kind == "x" else [0, 0]
        y = seg.x() if kind == "y" else [0, 0]
        ax.plot(x, y, color=c, lw=1)


def create_figure(pan, seq_lengths, msu_dict, sign_dict, block_pos):
    # create a figure with equal axis
    fig, axs = plt.subplots(
        2,
        2,
        figsize=(15, 15),
        sharex="col",
        sharey="row",
        gridspec_kw={
            "width_ratios": [1, 0.1],
            "height_ratios": [0.1, 1],
            "wspace": 0.05,
            "hspace": 0.05,
        },
    )

    ax = axs[1, 0]
    ax.set_aspect("equal", "box")

    x_lab, y_lab = pan.strains()

    x_path = pan.paths[x_lab]
    Bx = x_path.block_ids
    Sx = x_path.block_strands
    Ox = x_path.block_nums
    Lx = seq_lengths[x_lab]

    y_path = pan.paths[y_lab]
    By = y_path.block_ids
    Sy = y_path.block_strands
    Oy = y_path.block_nums
    Ly = seq_lengths[y_lab]

    msu_color, msu_n = {}, 0
    msu_labelled = set()
    for i, (b, s, o) in enumerate(zip(Bx, Sx, Ox)):
        x = block_pos[x_lab][(b, s, o)]
        msu_x = msu_dict[(x_lab, b, s, o)]
        sign_x = sign_dict[(x_lab, b, s, o)]
        col = "black"
        if msu_x > 0:
            if msu_x not in msu_color:
                msu_color[msu_x] = f"C{msu_n}"
                msu_n += 1
            col = msu_color[msu_x]
        for j, (c, d, p) in enumerate(zip(By, Sy, Oy)):
            if c != b:
                continue
            y = block_pos[y_lab][(c, d, p)]
            msu_y = msu_dict[(y_lab, c, d, p)]
            sign_y = sign_dict[(y_lab, c, d, p)]
            if msu_x != msu_y:
                continue
            if (msu_x != 0) and (sign_x != sign_y):
                continue
            orient = s == d
            if not msu_x in msu_labelled:
                ax.text(
                    np.mean(x),
                    np.mean(y),
                    f"{msu_x}",
                    fontsize="small",
                    color=col,
                    ha="center",
                    va="top",
                )
                msu_labelled |= {msu_x}
            seg = su.Segment(x[0], x[1], y[0], y[1], orient, Lx, Ly)
            pbc_plot(seg, col, orient, ax)

    # set minor ticks every 1e5 bp
    ax.set_xticks(np.arange(0, Lx, 1e5), minor=True)
    ax.set_yticks(np.arange(0, Ly, 1e5), minor=True)
    ax.set_xlim(0, Lx)
    ax.set_ylim(0, Ly)
    ax.grid(True, alpha=0.3)
    ax.grid(which="minor", alpha=0.1)

    private_x = set(pan.paths[x_lab].block_ids) - set(pan.paths[y_lab].block_ids)
    private_y = set(pan.paths[y_lab].block_ids) - set(pan.paths[x_lab].block_ids)
    for ax, lab, private in zip(
        [axs[0, 0], axs[1, 1]], [x_lab, y_lab], [private_x, private_y]
    ):
        path = pan.paths[lab]
        B = path.block_ids
        S = path.block_strands
        O = path.block_nums
        P = path.block_positions
        N = len(B)
        for i, (b, s, o) in enumerate(zip(B, S, O)):
            if b in private:
                pos = block_pos[lab][(b, s, o)]
                seg = su.PrivSegment(pos[0], pos[1], seq_lengths[lab])
                pbc_plot_priv(seg, "red", ax, "x" if lab == x_lab else "y")

        ax.grid(True, alpha=0.3)
        ax.grid(which="minor", alpha=0.1)

    axs[0, 0].set_yticks([])
    axs[0, 0].set_ylabel(f"{x_lab} private")
    axs[1, 1].set_xticks([])
    axs[1, 1].set_xlabel(f"{y_lab} private")
    axs[0, 1].axis("off")
    axs[1, 0].set_xlabel(f"{x_lab} (bp)")
    axs[1, 0].set_ylabel(f"{y_lab} (bp)")

    # msu_legend(fig, msu_color)

    return fig, axs


# def msu_legend(fig, msu_color):
#     ax = fig.add_subplot(111)
#     for msu, color in msu_color.items():
#         ax.plot([], [], color=color, label=f"{msu}")
#     ax.legend(loc="upper right", title="MSU", fontsize="small", title_fontsize="small")
#     ax.axis("off")


if __name__ == "__main__":
    args = parse_args()
    pan, seq_lengths, msu_dict, sign_dict = load_data(args)
    block_pos = block_positions(pan)
    fig, axs = create_figure(pan, seq_lengths, msu_dict, sign_dict, block_pos)
    fig.savefig(args.out)
    plt.close(fig)

# %%
