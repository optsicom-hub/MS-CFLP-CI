#!/usr/bin/env python3
"""
Generates a "remote customer" variant of an SS-CFLP-CI instance (wlpXX .dzn format).

Only the SupplyCost row of ONE customer is rewritten. The rest of the instance (the costs of the
other customers, FixedCost, Capacity, Goods and incompatibilities) remains byte-for-byte identical.

Construction:
  1. M = maximum cost over the whole SupplyCost matrix, excluding the remote customer.
  2. Each cost of the remote customer is M + min(round(|N(0, sigma)|), cap). It never falls below M
     and never exceeds M + cap (default cap = 10).
  3. One randomly chosen facility is set exactly to M, so the minimum cost of the remote customer
     is M and Dmaxmin = max_j min_i c_ij becomes exactly M. The constructive phase
     then retains all the arcs of the other customers.

For each instance, it prints M, the original and the new Dmaxmin, and checks that nothing else changes.

Usage:
  python remote_customer.py wlp07.txt wlp15.txt --out-dir remote/
  (options: --customer 0  --sigma 5  --cap 10  --seed 12345)
"""
import argparse
import re
from pathlib import Path

import numpy as np

SUPPLY_RE = re.compile(r"(SupplyCost\s*=\s*\[\|)(.*?)(\|\]\s*;)", re.S)


def parse(block):
    parts = block.split("|")              
    rows = [[int(x) for x in p.replace("\n", " ").split(",") if x.strip()] for p in parts]
    return parts, np.array(rows)


def dmaxmin(c):
    return int(c.min(axis=1).max())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("instances", nargs="+")
    ap.add_argument("--out-dir", default="remote")
    ap.add_argument("--customer", type=int, default=0, help="indice (desde 0) del cliente remoto")
    ap.add_argument("--sigma", type=float, default=5.0)
    ap.add_argument("--cap", type=int, default=10, help="maximo desplazamiento hacia arriba sobre M")
    ap.add_argument("--seed", type=int, default=12345)
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    print(f"{'instancia':<16}{'S':>6}{'W':>6}{'M':>6}{'Dmaxmin orig':>14}{'Dmaxmin nuevo':>15}{'min remoto':>12}{'max remoto':>12}")
    for path in args.instances:
        p = Path(path)
        text = p.read_text()
        m = SUPPLY_RE.search(text)
        parts, c = parse(m.group(2))
        S, W = c.shape
        j = args.customer
        if not 0 <= j < S:
            raise SystemExit(f"{p.name}: el cliente {j} no existe (S = {S})")

        others = np.delete(c, j, axis=0)
        M = int(others.max())
        d_orig = dmaxmin(c)

        off = np.minimum(np.rint(np.abs(rng.normal(0.0, args.sigma, W))).astype(int), args.cap)
        row = M + off
        row[rng.integers(W)] = M          
        new = c.copy()
        new[j] = row

        
        orig = parts[j]                    
        lead, trail = orig[:len(orig) - len(orig.lstrip())], orig[len(orig.rstrip()):]
        parts[j] = lead + ", ".join(str(int(v)) for v in row) + trail
        new_text = text[:m.start(2)] + "|".join(parts) + text[m.end(2):]
        name = f"{p.stem}_remote{p.suffix}"
        (out / name).write_text(new_text)

        m2 = SUPPLY_RE.search((out / name).read_text())
        _, c2 = parse(m2.group(2))
        assert c2.shape == c.shape and (c2 == new).all()
        assert (np.delete(c2, j, axis=0) == others).all()
        assert dmaxmin(c2) == M and c2[j].min() == M and c2[j].max() <= M + args.cap
        assert SUPPLY_RE.sub("", (out / name).read_text()) == SUPPLY_RE.sub("", text)
        old_parts = m.group(2).split("|")
        new_parts = m2.group(2).split("|")
        assert len(old_parts) == len(new_parts)
        assert all(a == b for k, (a, b) in enumerate(zip(old_parts, new_parts)) if k != j)
        print(f"{p.stem:<16}{S:>6}{W:>6}{M:>6}{d_orig:>14}{dmaxmin(c2):>15}{int(c2[j].min()):>12}{int(c2[j].max()):>12}")


if __name__ == "__main__":
    main()