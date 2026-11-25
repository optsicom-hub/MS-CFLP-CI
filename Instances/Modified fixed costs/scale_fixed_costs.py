#!/usr/bin/env python3
"""
Generates variants of SS-CFLP-CI instances (.dzn format of wlp01...wlp19) in which the
fixed opening costs (FixedCost) are multiplied by a factor alpha. The rest of the
instance (capacities, demands, SupplyCost assignment costs and incompatibilities)
is left untouched, so the feasible region is the same as that of the original instance.

Costs are integers in the code, so each scaled FixedCost is rounded up and bounded
below by 1 (zero costs are not allowed).

For each variant, it prints the ratio between the mean fixed cost and the mean
assignment cost, to describe the cost balance in the paper.

Usage:
  python scale_fixed_costs.py wlp01.txt wlp03.txt --alpha 0.25 4 --out-dir alpha/
"""
import argparse
import re
from pathlib import Path

import numpy as np

FIXED_RE = re.compile(r"(FixedCost\s*=\s*\[)(.*?)(\]\s*;)", re.S)
SUPPLY_RE = re.compile(r"SupplyCost\s*=\s*\[\|(.*?)\|\];", re.S)


def tag(alpha):
    return f"{alpha:g}"


def read_fixed(text):
    m = FIXED_RE.search(text)
    if m is None:
        raise ValueError("FixedCost not found")
    return np.array([int(x) for x in m.group(2).replace("\n", " ").split(",") if x.strip()])


def mean_supply(text):
    rows = [r for r in SUPPLY_RE.search(text).group(1).split("|") if r.strip()]
    return float(np.mean([[float(x) for x in r.split(",")] for r in rows]))


def scale(fixed, alpha):
    scaled = np.ceil(fixed * alpha).astype(int)   
    clipped = int((scaled < 1).sum())
    return np.maximum(scaled, 1), clipped


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("instances", nargs="+")
    ap.add_argument("--alpha", nargs="+", type=float, default=[0.25, 4.0])
    ap.add_argument("--out-dir", default="alpha_instances")
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    print(f"{'instancia':<28}{'alpha':>6}{'media F':>10}{'media C':>10}{'F/C':>8}{'a 1':>6}")
    for path in args.instances:
        p = Path(path)
        text = p.read_text()
        fixed = read_fixed(text)
        mc = mean_supply(text)
        print(f"{p.stem:<28}{'1':>6}{fixed.mean():>10.1f}{mc:>10.1f}{fixed.mean() / mc:>8.2f}{'':>6}")
        for a in args.alpha:
            new, clipped = scale(fixed, a)
            block = "FixedCost = [" + ", ".join(str(int(v)) for v in new) + "];"
            new_text = FIXED_RE.sub(lambda m: block, text, count=1)
            name = f"{p.stem}_alpha{tag(a)}{p.suffix}"
            (out / name).write_text(new_text)
            chk = read_fixed((out / name).read_text())
            assert (chk == new).all() and len(chk) == len(fixed)
            assert FIXED_RE.sub("", new_text) == FIXED_RE.sub("", text)
            print(f"{name:<28}{a:>6g}{new.mean():>10.1f}{mc:>10.1f}{new.mean() / mc:>8.2f}{clipped:>6}")


if __name__ == "__main__":
    main()
