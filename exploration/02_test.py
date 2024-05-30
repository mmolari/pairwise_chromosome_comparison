# %%
import pypangraph as pp
from collections import defaultdict
import itertools as itt
import plotly.graph_objects as go
import plotly.subplots as sp
from dataclasses import dataclass

fname = "../results/CA/graph.json"
pan = pp.Pangraph.load_json(fname)
bdf = pan.to_blockstats_df()
# %%


@dataclass
class Segment:
    s1: int
    e1: int
    s2: int
    e2: int
    orient: bool
    L1: int
    L2: int


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

    def display_line(bid, s1, e1, s2, e2, strand1, strand2, color, lg, text):
        kwargs = dict(
            strand1=strand1,
            strand2=strand2,
            color=color,
            lg=lg,
            text=text,
        )
        if s1 > e1:
            return
            display_line(bid, s1, Ls[0], s2, e2, **kwargs)
            display_line(bid, 0, e1, s2, e2, **kwargs)
        elif s2 > e2:
            return
            display_line(bid, s1, e1, s2, Ls[1], **kwargs)
            display_line(bid, s1, e1, 0, e2, **kwargs)
        else:
            fig.add_trace(
                go.Scatter(
                    x=[s1, e1],
                    y=[s2, e2] if strand1 == strand2 else [e2, s2],
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

            if dupl:
                lg = "dupl"
            else:
                lg = "fwd" if strand1 == strand2 else "inverted"

            color = colors[lg]
            display_line(
                bid, start1, end1, start2, end2, strand1, strand2, color, lg, text
            )

    private_1 = set(pos[p1].keys()) - set(pos[p2].keys())
    private_2 = set(pos[p2].keys()) - set(pos[p1].keys())
    for private, pid, r, c in [(private_1, p1, 2, 2), (private_2, p2, 1, 1)]:
        for bid in private:
            l = pos[pid][bid]
            dupl = len(l) > 1
            for a in l:
                start, end, strand, occ = a
                pm = "+" if strand == 1 else "-"
                text = f"{bid} # {pm}|{occ}"

                if dupl:
                    lg = "dupl"
                else:
                    lg = f"private {pid}"

                x = [start, end] if pid == p1 else [0, 0]
                y = [start, end] if pid == p2 else [0, 0]

                color = colors[lg]
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


Ls = {"A": 4638422, "ref": 4641652}
pos = position_dictionary(pan)

fig = create_dotplot(pos, Ls)
fig.write_html("dotplot.html")

# %%
