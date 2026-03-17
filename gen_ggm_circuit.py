#!/usr/bin/env python3
"""
gen_ggm_circuit.py — Generate a BLIF circuit for the GGM PRF.

Usage:
    python3 gen_ggm_circuit.py --key-bits N [--depth D] [--seed S] [--out FILE]

Arguments:
    --key-bits N   Key length in bits (= node size = PRF output size). Required.
    --depth D      Input length / tree depth (default: key-bits).
    --seed S       RNG seed for the PRG graph (default: 42).
    --out FILE     Output BLIF file (default: ggm_prf_<N>bit.blif).

The circuit computes:
    GGM_PRF.eval(key, x) → value

    Inputs:  key[0..N-1]   (N-bit PRF key)
             x[0..D-1]     (D-bit input, selects tree path)
    Outputs: out[0..N-1]   (N-bit PRF output)

PRG construction (Goldreich NC0, degree-3 predicate):
    G: {0,1}^N → {0,1}^{2N}
    The bipartite graph (which input bits each output bit depends on)
    is sampled once from the seed and hardcoded into the circuit.
    Each output bit: P(x_{i0}, x_{i1}, x_{i2}) = (x_{i0} AND x_{i1}) XOR x_{i2}

GGM tree traversal (one level per x[i]):
    node_0 = key
    for i = 0..D-1:
        expanded[0..2N-1] = G(node_i)        // PRG expansion
        left[0..N-1]      = expanded[0..N-1]
        right[0..N-1]     = expanded[N..2N-1]
        node_{i+1}[j]     = MUX(x[i], left[j], right[j])  // x[i]=0→left, 1→right
    output = node_D

Circuit size:
    Per level: 2N PRG output bits (each: 1 AND + 1 XOR = 2 gates)
             + N MUX gates (each: 4 gates)
    Total: D * N * (2 + 4) = 6*D*N gates
    Plus wiring constants: O(D*N)

    Example: N=64, D=64 → ~24,576 gates
"""

import argparse
import random
import sys
import os

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description="Generate a BLIF circuit for the GGM PRF (NC0 PRG based)")
    p.add_argument("--key-bits", type=int, required=True,
                   help="Key/node/output size in bits (N)")
    p.add_argument("--depth", type=int, default=None,
                   help="Tree depth / input length (D, default = key-bits)")
    p.add_argument("--seed", type=int, default=42,
                   help="RNG seed for PRG graph sampling (default: 42)")
    p.add_argument("--out", type=str, default=None,
                   help="Output BLIF filename")
    return p.parse_args()

# ---------------------------------------------------------------------------
# PRG graph sampling
# ---------------------------------------------------------------------------
def sample_prg_graph(n_in, n_out, locality, seed):
    """
    Sample the bipartite graph for the NC0 PRG.
    Returns graph[i] = [i0, i1, i2] — the 'locality' input indices for output bit i.
    Matches the PRG_NC0 construction in prg_nc0.h.
    """
    rng = random.Random(seed)
    graph = []
    for _ in range(n_out):
        nbrs = [rng.randint(0, n_in - 1) for _ in range(locality)]
        graph.append(nbrs)
    return graph

# ---------------------------------------------------------------------------
# BLIF generation helpers
# ---------------------------------------------------------------------------
class BLIFWriter:
    def __init__(self):
        self.lines = []
        self._gate_count = 0

    def w(self, s=""):
        self.lines.append(s)

    def new_wire(self, prefix="w"):
        self._gate_count += 1
        return f"{prefix}_{self._gate_count}"

    def emit_and(self, a, b, out):
        """out = a AND b"""
        self.w(f".names {a} {b} {out}")
        self.w("11 1")
        self.w()

    def emit_xor(self, a, b, out):
        """out = a XOR b"""
        self.w(f".names {a} {b} {out}")
        self.w("10 1")
        self.w("01 1")
        self.w()

    def emit_mux(self, sel, a, b, out):
        """out = sel ? b : a   (sel=0 → a, sel=1 → b)"""
        # out = (~sel & a) | (sel & b)
        t1 = self.new_wire("mux_t1")
        t2 = self.new_wire("mux_t2")
        self.w(f".names {sel} {a} {out}")
        self.w("0 1 1")   # sel=0, a=1 → 1
        self.w()
        # Redo properly with 3-input .names
        self.lines.pop(); self.lines.pop(); self.lines.pop(); self.lines.pop()
        # Use 3-variable truth table: inputs are sel, a, b
        self.w(f".names {sel} {a} {b} {out}")
        self.w("0 1 - 1")   # sel=0, a=1 → out=1
        self.w("1 - 1 1")   # sel=1, b=1 → out=1
        self.w()

    def emit_const_0(self, out):
        self.w(f".names {out}")
        self.w()

    def emit_const_1(self, out):
        self.w(f".names {out}")
        self.w("1")
        self.w()

    def emit_buf(self, src, dst):
        """dst = src (buffer)"""
        self.w(f".names {src} {dst}")
        self.w("1 1")
        self.w()

    def total_gates(self):
        return sum(1 for l in self.lines if l.startswith(".names"))

# ---------------------------------------------------------------------------
# PRG evaluation circuit (one application of G: N bits → 2N bits)
# ---------------------------------------------------------------------------
def emit_prg(bw, graph, in_prefix, out_prefix, level):
    """
    Emit the circuit for one application of the NC0 PRG.

    graph[i] = [i0, i1, i2]: output bit i = P(in[i0], in[i1], in[i2])
    P(a,b,c) = (a AND b) XOR c

    in_prefix:  signal names in_{prefix}_b{j} for j in range(n_in)
    out_prefix: signal names out_{prefix}_b{i} for i in range(n_out)
    """
    n_out = len(graph)
    for i, (i0, i1, i2) in enumerate(graph):
        a   = f"{in_prefix}_b{i0}"
        b   = f"{in_prefix}_b{i1}"
        c   = f"{in_prefix}_b{i2}"
        ab  = f"prg_l{level}_ab_{i}"
        out = f"{out_prefix}_b{i}"
        bw.emit_and(a, b, ab)
        bw.emit_xor(ab, c, out)

# ---------------------------------------------------------------------------
# MUX layer: select left or right half of expanded based on x[level]
# ---------------------------------------------------------------------------
def emit_mux_layer(bw, n, expanded_prefix, sel_sig, out_prefix, level):
    """
    For j in range(n):
        out[j] = MUX(sel, left=expanded[j], right=expanded[n+j])
    """
    for j in range(n):
        left_sig  = f"{expanded_prefix}_b{j}"
        right_sig = f"{expanded_prefix}_b{n + j}"
        out_sig   = f"{out_prefix}_b{j}"
        bw.emit_mux(sel_sig, left_sig, right_sig, out_sig)

# ---------------------------------------------------------------------------
# Main circuit generation
# ---------------------------------------------------------------------------
def generate_ggm_circuit(n, depth, seed):
    """
    Generate a BLIF circuit for GGM_PRF.eval(key, x).

    n:     key/node/output bits
    depth: tree depth (input length)
    seed:  PRG graph seed

    Returns the BLIF string.
    """
    bw = BLIFWriter()

    # Sample the PRG graph (matches PRG_NC0 construction with locality=3)
    # G: {0,1}^n → {0,1}^{2n}
    graph = sample_prg_graph(n_in=n, n_out=2*n, locality=3, seed=seed)

    # ── Input/output declarations ─────────────────────────────────────────
    key_sigs = " ".join(f"key_b{j}" for j in range(n))
    x_sigs   = " ".join(f"x_b{i}"  for i in range(depth))
    out_sigs = " ".join(f"out_b{j}" for j in range(n))

    bw.w(f"# GGM PRF circuit — key={n} bits, depth={depth}, PRG seed={seed}")
    bw.w(f"# PRG: Goldreich NC0, locality=3, predicate P(a,b,c) = (a AND b) XOR c")
    bw.w(f"# Inputs:  key[0..{n-1}] (PRF key), x[0..{depth-1}] (PRF input)")
    bw.w(f"# Outputs: out[0..{n-1}] (PRF value = GGM leaf)")
    bw.w(f"# Gates: ~{6 * n * depth} (AND+XOR per PRG bit, 4 per MUX)")
    bw.w()
    bw.w(".model ggm_prf")
    bw.w(f".inputs {key_sigs} {x_sigs}")
    bw.w(f".outputs {out_sigs}")
    bw.w()

    # ── GGM tree traversal ────────────────────────────────────────────────
    # node_0 = key
    # For each level i = 0..depth-1:
    #   expanded = G(node_i)       [2n bits]
    #   node_{i+1}[j] = MUX(x[i], left=expanded[j], right=expanded[n+j])
    # output = node_depth

    # Wire name for current node at each level
    def node_sig(level, bit):
        if level == 0:
            return f"key_b{bit}"
        return f"node_l{level}_b{bit}"

    for level in range(depth):
        bw.w(f"# === Level {level}: expand node, select child ===")

        # Emit PRG: node_level → expanded_level (2n bits)
        # Input wire names: node_sig(level, j) for j in range(n)
        # We need a uniform prefix for the PRG input, so name it
        in_prefix  = f"node_l{level}" if level > 0 else "key"
        exp_prefix = f"exp_l{level}"

        # If level==0, input is key_b{j} directly
        if level == 0:
            # PRG reads from key_b{j} — wire names already match in_prefix=key
            pass
        # else: reads from node_l{level}_b{j} — already named correctly

        emit_prg(bw, graph, in_prefix, exp_prefix, level)
        bw.w()

        # Emit MUX layer: select left/right based on x[level]
        out_node_prefix = f"node_l{level+1}"
        if level == depth - 1:
            # Last level: output goes to final output wires
            out_node_prefix = "out"

        emit_mux_layer(bw, n, exp_prefix, f"x_b{level}", out_node_prefix, level)
        bw.w()

    bw.w(".end")

    return "\n".join(bw.lines), bw.total_gates()

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def main():
    args = parse_args()

    n     = args.key_bits
    depth = args.depth if args.depth is not None else n
    seed  = args.seed

    if n < 2:
        print("Error: --key-bits must be >= 2", file=sys.stderr)
        sys.exit(1)
    if depth < 1:
        print("Error: --depth must be >= 1", file=sys.stderr)
        sys.exit(1)

    out_file = args.out or f"ggm_prf_{n}bit_d{depth}.blif"

    print(f"Generating GGM PRF circuit:")
    print(f"  Key/output bits : {n}")
    print(f"  Depth (input)   : {depth}")
    print(f"  PRG seed        : {seed}")
    print(f"  Output file     : {out_file}")

    blif, n_gates = generate_ggm_circuit(n, depth, seed)

    with open(out_file, "w") as f:
        f.write(blif)

    n_lines = blif.count("\n")
    print(f"  Lines           : {n_lines:,}")
    print(f"  .names gates    : {n_gates:,}")
    print(f"  Estimated gates : ~{6 * n * depth:,}  (AND+XOR per PRG output + MUX)")
    print(f"Done.")

if __name__ == "__main__":
    main()
