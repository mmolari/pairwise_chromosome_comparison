import pypangraph as pp
import pandas as pd
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--graph", help="Path to the pangraph json file")
    parser.add_argument("--output", help="Path to the output csv file")
    return parser.parse_args()


def block_position_dataframe(pan):

    block_positions = []
    for path in pan.paths:
        iso = path.name
        B = path.block_ids
        S = path.block_strands
        O = path.block_nums
        P = path.block_positions

        L = len(B)
        for i in range(len(B)):
            b = B[i]
            s = S[i]
            o = O[i]
            start = P[i]
            end = P[(i + 1) % L]
            block_positions.append(
                {
                    "genome": iso,
                    "block_id": b,
                    "strand": s,
                    "occurrence_number": o,
                    "start_position": start,
                    "end_position": end,
                }
            )

    return pd.DataFrame(block_positions)


if __name__ == "__main__":
    args = parse_args()
    pan = pp.Pangraph.load_json(args.graph)
    df = block_position_dataframe(pan)
    df.to_csv(args.output, index=False)
