# %%
import pypangraph as pp
from collections import defaultdict
import itertools as itt
import plotly.graph_objects as go

fname = "../results/CA/graph.json"
pan = pp.Pangraph.load_json(fname)
bdf = pan.to_blockstats_df()
# %%


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
    fig = go.Figure()

    colors = {
        "fwd": "blue",
        "inverted": "red",
        "dupl": "gray",
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
            )
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
                )
            )

    # Add scatter trace for alignments
    p1, p2 = list(pos.keys())
    for bid in pos[p1]:
        if not (bid in pos[p2]):
            continue
        l1 = pos[p1][bid]
        l2 = pos[p2][bid]

        dupl = len(l1) > 1 or len(l2) > 1
        for a1, a2 in itt.product(l1, l2):
            start1, end1, strand1, occ1 = a1
            start2, end2, strand2, occ2 = a2

            pm1 = "+" if strand1 == 1 else "-"
            pm2 = "+" if strand2 == 1 else "-"
            text = f"{bid} - {pm1}|{occ1} vs {pm2}|{occ2}"

            if dupl:
                lg = "dupl"
            else:
                lg = "fwd" if strand1 == strand2 else "inverted"

            color = colors[lg]
            display_line(
                bid, start1, end1, start2, end2, strand1, strand2, color, lg, text
            )

    # Update layout for better appearance
    fig.update_layout(
        title="Genome Alignment Dotplot",
        xaxis=dict(title=f"{p1} Coordinates"),
        yaxis=dict(title=f"{p2} Coordinates"),
        showlegend=True,
        width=1200,
        height=1200,
    )

    return fig


pos = position_dictionary(pan)

fig = create_dotplot(pos, {})
fig.write_html("dotplot.html")

# %%
