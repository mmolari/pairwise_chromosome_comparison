# %%

import pypangraph as pp
import glue_utils as gu
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--graph", type=str, required=True)
    parser.add_argument("--out", type=str, required=True)
    return parser.parse_args()


if __name__ == "__main__":

    args = parse_args()

    pan = pp.Pangraph.load_json(args.graph)
    bdf = pan.to_blockstats_df()

    paths = gu.pan_to_paths(pan)
    k1, k2 = list(paths.keys())

    glue = gu.Glue(paths)

    core_blocks = set(bdf[bdf["core"]].index)

    core_nodes = {cb: {} for cb in core_blocks}
    for k, p in paths.items():
        for n in p.nodes:
            if n.bid in core_blocks:
                core_nodes[n.bid][k] = n

    msu_id = 1
    for cb, cns in core_nodes.items():
        n1 = cns[k1]
        n2 = cns[k2]
        glue.extend(n1, n2, msu_id)
        msu_id += 1

    df = glue.to_df()
    df.to_csv(args.out, index=False)
# %%
