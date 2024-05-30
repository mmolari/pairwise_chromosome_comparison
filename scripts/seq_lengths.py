from Bio import SeqIO
import pandas as pd
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastas", type=str, nargs="+", help="Fasta files")
    parser.add_argument("--output", type=str, help="Output file")
    return parser.parse_args()


def get_seq_lengths(fasta_file):
    seq_lengths = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        L = len(record.seq)
        seq_lengths.append(
            {
                "id": record.id,
                "length": L,
                "file": fasta_file,
            }
        )
    return seq_lengths


if __name__ == "__main__":
    args = parse_args()

    seq_lengths = []
    for fasta_file in args.fastas:
        seq_lengths.extend(get_seq_lengths(fasta_file))
    df = pd.DataFrame(seq_lengths)

    df.to_csv(args.output, index=False)
