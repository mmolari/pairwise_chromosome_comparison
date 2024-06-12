import pandas as pd
from dataclasses import dataclass

verbose = False


def vprint(*args, indent=0, **kwargs):
    if verbose:
        print("--" * indent + " ", *args, **kwargs)


@dataclass(frozen=True)
class Node:
    bid: str
    strand: bool
    occ: int

    def __repr__(self) -> str:
        s = "+" if self.strand else "-"
        return f"[{self.bid}|{s}|{self.occ}]"

    def gluable(self, other: object) -> bool:
        return self.bid == other.bid and self.strand == other.strand

    def invert(self) -> "Node":
        return Node(self.bid, not self.strand, self.occ)

    def anonymous_hash(self) -> int:
        return hash((self.bid, self.strand))


@dataclass(frozen=True, unsafe_hash=False, eq=False)
class Edge:
    node1: Node
    node2: Node

    def invert(self) -> "Edge":
        return Edge(self.node2.invert(), self.node1.invert())

    def __side_eq__(self, value: object) -> bool:
        return self.node1.gluable(value.node1) and self.node2.gluable(value.node2)

    def __eq__(self, value: object) -> bool:
        return self.__side_eq__(value) or self.__side_eq__(value.invert())

    def __side_hash__(self) -> int:
        return hash((self.node1.anonymous_hash(), self.node2.anonymous_hash()))

    def __hash__(self) -> int:
        return self.__side_hash__() + self.invert().__side_hash__()


@dataclass(init=False, unsafe_hash=False, eq=False)
class Path:
    nodes: list[Node]
    node_pos: dict[Node, int]
    msu: list[int]
    signature: list[int]

    def __init__(self, nodes: list[Node]) -> None:
        self.nodes = nodes
        self.node_pos = {n: i for i, n in enumerate(nodes)}
        self.msu = [0] * len(nodes)
        self.signature = [0] * len(nodes)

    def __repr__(self) -> str:
        return "---".join([str(n) for n in self.nodes])

    @staticmethod
    def from_pangraph_path(path) -> "Path":
        B = path.block_ids
        S = path.block_strands
        O = path.block_nums
        nodes = [Node(b, s, o) for b, s, o in zip(B, S, O)]
        return Path(nodes)

    def node_idx(self, node: Node) -> int:
        return self.node_pos[node]

    def next_node_idx(self, node: Node, fwd: bool) -> int:
        pos = self.node_pos[node]
        if fwd:
            return (pos + 1) % len(self.nodes)
        else:
            return (pos - 1) % len(self.nodes)

    def next_node(self, node: Node, fwd: bool) -> Node:
        return self.nodes[self.next_node_idx(node, fwd)]

    def next_edge(self, node: Node, fwd: bool) -> Edge:
        if fwd:
            return Edge(node, self.next_node(node, fwd))
        else:
            return Edge(self.next_node(node, fwd), node)

    def get_msu(self, node: Node) -> int:
        return self.msu[self.node_pos[node]]

    def set_msu(self, node: Node, msu_id: int, signature: int) -> None:
        idx = self.node_pos[node]
        self.msu[idx] = msu_id
        self.signature[idx] = signature

    def remap_msu(self, old_id: int, new_id: int) -> None:
        for i, msu in enumerate(self.msu):
            if msu == old_id:
                self.msu[i] = new_id

    def to_df(self) -> pd.DataFrame:
        data = [n.__dict__ for n in self.nodes]
        df = pd.DataFrame(data)
        df["msu"] = self.msu
        df["signature"] = self.signature
        return df


def pan_to_paths(pan):
    res = {}
    for path in pan.paths:
        name = path.name
        res[name] = Path.from_pangraph_path(path)
    return res


class Glue:
    def __init__(self, path_dict) -> None:
        assert len(path_dict) == 2

        self.k1, self.k2 = list(path_dict.keys())
        self.path1 = path_dict[self.k1]
        self.path2 = path_dict[self.k2]

    def remap_msu(self, old_id: int, new_id: int) -> None:
        self.path1.remap_msu(old_id, new_id)
        self.path2.remap_msu(old_id, new_id)

    def extend(self, n1: Node, n2: Node, msu_id: int) -> None:
        vprint(f"Extending {n1} and {n2}", indent=1)
        assert n1.bid == n2.bid

        # if already assigned skip
        if (x := self.path1.get_msu(n1)) != 0:
            assert self.path2.get_msu(n2) == x
            vprint(f"Already assigned {n1} and {n2} to {x}", indent=1)
            return

        # if both are unassigned, extend
        vprint(f"Assigning msu of focal nodes {n1} and {n2} to {msu_id}", indent=1)
        signature = hash((n1, n2))
        self.path1.set_msu(n1, msu_id, signature)
        self.path2.set_msu(n2, msu_id, signature)
        c1 = self.glue_side(n1, n2, True, msu_id)
        c2 = self.glue_side(n1, n2, False, msu_id)

        for c in [c1, c2]:
            if c is not None:
                self.remap_msu(c[0], c[1])

    def compare(self, n1: int, n2: int, s1: bool, s2: bool) -> bool:
        e1 = self.path1.next_edge(n1, s1)
        e2 = self.path2.next_edge(n2, s2)
        vprint(f"Comparing {e1} and {e2}", indent=3)
        assert n1.bid == n2.bid, f"{n1} != {n2}"
        vprint(f"Result {e1 == e2}", indent=3)
        return e1 == e2

    def glue_side(self, n1: Node, n2: Node, fwd: bool, msu_id: int):
        vprint(f"Gluing side {fwd} of {n1} and {n2}", indent=2)

        s1 = n1.strand if fwd else not n1.strand
        s2 = n2.strand if fwd else not n2.strand

        vprint(f"Focal nodes {n1} and {n2}", indent=2)

        while self.compare(n1, n2, s1, s2):
            # next nodes
            nn1 = self.path1.next_node(n1, s1)
            nn2 = self.path2.next_node(n2, s2)
            vprint(f"Next nodes {nn1} and {nn2}", indent=3)
            # if already assigned, skip
            old_1 = self.path1.get_msu(nn1)
            old_2 = self.path2.get_msu(nn2)
            if old_1 != 0 or old_2 != 0:
                vprint(
                    f"## Already assigned {nn1} -> {old_1} or {nn2} -> {old_2}",
                    indent=3,
                )
                if old_1 == old_2:
                    vprint(f"## Reassigning {old_1} -> {msu_id}", indent=3)
                    return (old_1, msu_id)
                else:
                    return None
            vprint(f"Assigning {nn1} and {nn2} to {msu_id}", indent=3)
            signature = hash((nn1, nn2))
            self.path1.set_msu(nn1, msu_id, signature)
            self.path2.set_msu(nn2, msu_id, signature)
            n1, n2 = nn1, nn2

    def to_df(self) -> pd.DataFrame:

        df1 = self.path1.to_df()
        df1["path"] = self.k1
        df2 = self.path2.to_df()
        df2["path"] = self.k2
        df = pd.concat([df1, df2], axis=0)

        # remap msus:
        msus = sorted(set(df["msu"]) | {0})
        msu_map = {msu: i for i, msu in enumerate(msus)}
        msu_map[0] == None
        df["msu"] = df["msu"].map(msu_map)
        return df
