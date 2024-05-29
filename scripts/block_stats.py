import argparse
import pypangraph as pp


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--graph", type=str, help="input pangraph")
    parser.add_argument("--out", type=str, help="output block stats")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    pan = pp.Pangraph.load_json(args.graph)
    bdf = pan.to_blockstats_df()
    bdf = bdf.sort_values(
        ["core", "duplicated", "count", "len"], ascending=[False, True, False, False]
    )
    bdf["category"] = "core"
    mask = (~bdf["core"]) & bdf["duplicated"]
    bdf.loc[mask, "category"] = "duplicated"
    mask = (~bdf["core"]) & (~bdf["duplicated"])
    bdf.loc[mask, "category"] = "accessory"
    bdf.to_csv(args.out)
