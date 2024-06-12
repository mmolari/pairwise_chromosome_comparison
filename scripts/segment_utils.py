import copy
from dataclasses import dataclass


@dataclass
class Segment:
    s1: int
    e1: int
    s2: int
    e2: int
    orient: bool
    L1: int
    L2: int

    def __post_init__(self):
        self.s1 %= self.L1
        self.e1 %= self.L1
        self.s2 %= self.L2
        self.e2 %= self.L2

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
        return (self.s1 > self.e1) and (self.e1 != 0)

    def y_runover(self):
        return (self.s2 > self.e2) and (self.e2 != 0)

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
