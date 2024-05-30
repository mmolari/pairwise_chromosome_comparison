import pypangraph as pp
import pandas as pd
from collections import defaultdict
import itertools as itt
import plotly.graph_objects as go
import plotly.subplots as sp
from dataclasses import dataclass
import copy
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--graph", type=str, help="Pangraph JSON file")
    parser.add_argument("--seq_lengths", type=str, help="Sequence lengths CSV file")
    parser.add_argument("--output", type=str, help="Output HTML file")
    return parser.parse_args()


def position_dictionary(pan):
    strains = pan.strains()
    pos = {pid: defaultdict(list) for pid in strains}
    for pid in strains:
        path = pan.paths[pid]
        B = path.block_ids
        L = len(B)
        for i, bid in enumerate(B):
            start = path.block_positions[i]
            end = path.block_positions[(i + 1) % L]
            strand = path.block_strands[i]
            occ = path.block_nums[i]
            pos[pid][bid].append((start, end, strand, occ))
    return pos


@dataclass
class Segment:
    s1: int
    e1: int
    s2: int
    e2: int
    orient: bool
    L1: int
    L2: int

    def x(self):
        return (self.s1, self.e1)

    def y(self):
        return (self.s2, self.e2) if self.orient else (self.e2, self.s2)

    def dx(self):
        return (self.e1 - self.s1) % self.L1

    def dy(self):
        return (self.e2 - self.s2) % self.L2

    def m(self):
        m = self.dy() / self.dx()
        return m if self.orient else -m

    def x_runover(self):
        return self.s1 > self.e1

    def y_runover(self):
        return self.s2 > self.e2

    def split_x(self):
        assert self.x_runover()
        Sa, Sb = copy.copy(self), copy.copy(self)
        Sa.e1 = self.L1
        Sb.s1 = 0
        if self.orient:
            Sa.e2 = Sa.s2 + (Sa.e1 - Sa.s1) * self.m()
            Sb.s2 = Sb.e2 - (Sb.e1 - Sb.s1) * self.m()
        else:
            Sa.s2 = Sa.e2 + (Sa.e1 - Sa.s1) * self.m()
            Sb.e2 = Sb.s1 - (Sb.e1 - Sb.s1) * self.m()

        return [Sa, Sb]

    def split_y(self):
        assert self.y_runover()
        Sa, Sb = copy.copy(self), copy.copy(self)
        Sa.e2 = self.L2
        Sb.s2 = 0
        if self.orient:
            Sa.e1 = Sa.s1 + (Sa.e2 - Sa.s2) / self.m()
            Sb.s1 = Sb.e1 - (Sb.e2 - Sb.s2) / self.m()
        else:
            Sa.s1 = Sa.e1 + (Sa.e2 - Sa.s2) / self.m()
            Sb.e1 = Sb.s1 - (Sb.e2 - Sb.s2) / self.m()

        return [Sa, Sb]


def display_line(bid, seg, color, lg, text, fig):
    kwargs = dict(color=color, lg=lg, text=text, fig=fig)
    if seg.s1 > seg.e1:
        for s in seg.split_x():
            display_line(bid, s, **kwargs)
    elif seg.s2 > seg.e2:
        for s in seg.split_y():
            display_line(bid, s, **kwargs)
    else:
        fig.add_trace(
            go.Scatter(
                x=seg.x(),
                y=seg.y(),
                mode="markers+lines",
                marker=dict(color=color, size=2),
                line=dict(color=color),
                name=f"{bid}",
                text=text,
                hoverinfo="text",
                legendgroup=lg,
                showlegend=False,
            ),
            row=1,
            col=2,
        )


@dataclass
class PrivSegment:
    s: int
    e: int
    L: int

    def x(self):
        return (self.s, self.e)

    def dx(self):
        return (self.e - self.s) % self.L

    def x_runover(self):
        return self.s > self.e

    def split_x(self):
        assert self.x_runover()
        Sa, Sb = copy.copy(self), copy.copy(self)
        Sa.e = self.L
        Sb.s = 0
        return [Sa, Sb]


def display_private_line(bid, seg, color, lg, text, fig, kind):
    if seg.s > seg.e:
        for s in seg.split_x():
            display_private_line(bid, s, color, lg, text, fig, kind)
    else:
        x = seg.x() if kind == "x" else [0, 0]
        y = seg.x() if kind == "y" else [0, 0]
        r, c = (2, 2) if kind == "x" else (1, 1)
        fig.add_trace(
            go.Scatter(
                x=x,
                y=y,
                mode="markers+lines",
                marker=dict(color=color, size=2),
                line=dict(color=color),
                name=f"{bid}",
                text=text,
                hoverinfo="text",
                legendgroup=lg,
                showlegend=False,
            ),
            row=r,
            col=c,
        )


def create_dotplot(pos, Ls):
    # Create the plotly figure
    fig = sp.make_subplots(
        2,
        2,
        shared_xaxes="columns",
        shared_yaxes="rows",
        column_widths=[0.1, 0.9],
        row_heights=[0.9, 0.1],
        horizontal_spacing=0.01,
        vertical_spacing=0.01,
    )

    p1, p2 = list(pos.keys())

    colors = {
        "fwd": "blue",
        "inverted": "red",
        "dupl": "gray",
        f"private {p1}": "green",
        f"private {p2}": "goldenrod",
    }

    for lg, color in colors.items():
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                marker=dict(color=color, size=5),
                name=lg,
                showlegend=True,
                legendgroup=lg,
            ),
            row=1,
            col=2,
        )

    shared = set(pos[p1].keys()) & set(pos[p2].keys())
    for bid in shared:
        l1 = pos[p1][bid]
        l2 = pos[p2][bid]

        dupl = len(l1) > 1 or len(l2) > 1
        for a1, a2 in itt.product(l1, l2):
            start1, end1, strand1, occ1 = a1
            start2, end2, strand2, occ2 = a2

            pm1 = "+" if strand1 == 1 else "-"
            pm2 = "+" if strand2 == 1 else "-"
            text = f"{bid} # {pm1}|{occ1} vs {pm2}|{occ2}"

            orient = strand1 == strand2

            if dupl:
                lg = "dupl"
            else:
                lg = "fwd" if orient else "inverted"

            seg = Segment(start1, end1, start2, end2, orient, Ls[p1], Ls[p2])
            color = colors[lg]
            display_line(bid, seg, color, lg, text, fig)

    private_1 = set(pos[p1].keys()) - set(pos[p2].keys())
    private_2 = set(pos[p2].keys()) - set(pos[p1].keys())
    for private, pid, kind in [(private_1, p1, "x"), (private_2, p2, "y")]:
        for bid in private:
            l = pos[pid][bid]
            dupl = len(l) > 1
            for a in l:
                start, end, strand, occ = a
                seg = PrivSegment(start, end, Ls[pid])
                pm = "+" if strand == 1 else "-"
                text = f"{bid} # {pm}|{occ}"

                if dupl:
                    lg = "dupl"
                else:
                    lg = f"private {pid}"

                color = colors[lg]
                display_private_line(bid, seg, color, lg, text, fig, kind)

    # Update layout for better appearance
    fig.update_layout(
        title="Genome Alignment Dotplot",
        showlegend=True,
        width=1200,
        height=1200,
    )

    fig.update_xaxes(title_text=f"{p1} genome (bp)", row=2, col=2)
    fig.update_yaxes(title_text=f"{p2} genome (bp)", row=1, col=1)

    return fig


if __name__ == "__main__":
    args = parse_args()

    pan = pp.Pangraph.load_json(args.graph)
    bdf = pan.to_blockstats_df()

    Ls = pd.read_csv(args.seq_lengths).set_index("id")["length"].to_dict()
    Ls = {k: v for k, v in Ls.items() if k in pan.strains()}
    pos = position_dictionary(pan)

    fig = create_dotplot(pos, Ls)
    fig.write_html(args.output)
