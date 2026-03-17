# JLS22 iO — C++ Implementation

A complete implementation of the JLS22 indistinguishability obfuscation (iO)
construction (Jain, Lin, Sahai 2022), together with a GGM puncturable PRF and
a SystemVerilog circuit generator for the hardware decryption path.

Built as the cryptographic engine behind
**[Black Box Chip](https://blackboxchip.com)** — provably secure logic locking
based on *Locked Circuit Indistinguishability* (IEEE CSF 2022).

---

## Repository layout

```
include/
  field.h           Zp arithmetic and noise sampling
  lpn.h             δ-LPN encryption over Fp
  prg_nc0.h         Goldreich NC0 PRG + PRG_H
  ppe.h             Preprocessed Polynomial Encoding (JLS22 §4)
  are.h             Amortized Randomized Encoding / Yao garbling (§5)
  pre.h             PPE ∘ ARE composition (§6)
  phfe.h            DLIN-based Partially Hiding FE (§7.2)
  symbolic_garble.h Symbolic garbling — input-independent KeyGen
  fe.h              PRE ∘ PHFE sublinear FE (§7.3)
  io.h              iO via AJ15/BV15 bootstrapping (§7.1)
  ggm_prf.h         GGM puncturable PRF built from PRG_NC0

test/
  test_primitives.cpp   field, LPN, PRG          (23 tests)
  test_ppe.cpp          PPE                       (18 tests)
  test_are.cpp          ARE / Yao garbling        (19 tests)
  test_pre.cpp          PRE                       (16 tests)
  test_phfe.cpp         PHFE                      (21 tests)
  test_fe.cpp           FE                        (17 tests)
  test_io.cpp           iO                        (21 tests)
  test_ggm_prf.cpp      GGM PRF + puncturing      (22 tests)

gen_ggm_circuit.py   Generate a BLIF circuit for GGM_PRF.eval(key, x)
gen_fe_verilog.py    Generate SystemVerilog for FE.Dec(SK_C, FE.Enc(PK, x))
```

---

## Building and testing

Requires C++17 and a compiler that supports it (GCC 9+ or Clang 10+).
All headers are self-contained — no external dependencies in prototype mode.

```bash
# Compile and run a single test suite
g++ -std=c++17 -O2 -I include -o test_io test/test_io.cpp
./test_io

# Compile and run all test suites
for t in primitives ppe are pre phfe fe io; do
    g++ -std=c++17 -O2 -I include -o test_$t test/test_$t.cpp && ./test_$t
done
```

Expected output: **135/135 tests pass**.

---

## Production pairing (mcl/BN254)

`phfe.h` ships with a simulated bilinear group (field arithmetic, not
cryptographically secure) so all tests run without external dependencies.
To enable real BN254 pairings via [mcl](https://github.com/herumi/mcl):

```bash
g++ -std=c++17 -O2 -DPHFE_USE_MCL \
    -I include -I /path/to/mcl/include \
    -o test_phfe test/test_phfe.cpp \
    /path/to/mcl/lib/libmcl.a
```

With `-DPHFE_USE_MCL` the following change behaviour:
- `setup()` populates G1-encoded public key matrices
- `keygen()` populates G1-encoded secret key vectors
- `enc()` r-blinds SI: `CT_si[j] = g^{⟨r,u⟩ + SI[j]}`
- `dec()` recovers `f(PI, SI)` via pairing accumulation and discrete log

The prototype `enc_prototype()` path is always compiled in so unit tests
remain runnable without mcl.

---

## Locking a circuit with LCI

The full Indistinguishable Locking workflow from the paper is:

**Step 1 — Size your GGM PRF.**
The PRF key length `λ` is determined by the security parameter from the paper.
The depth equals the number of **input** bits of your circuit.

```bash
python3 gen_ggm_circuit.py --key-bits 64 --depth 8 --out prf.blif
```

**Step 2 — Miter the PRF circuit with your circuit.**
The miter takes your original circuit `C` and the GGM PRF circuit and wires
them together so that the locked circuit output is `C(x) XOR F_k(x)` for a
secret key `k`. Use a standard EDA miter tool (e.g. Yosys `miter` command)
or the BDD-based miter from the companion locking service.

**Step 3 — Obfuscate with iO.**
Pass the miter circuit to `iO::obfuscate`. The obfuscated circuit is
functionally equivalent to the miter but hides its structure.

```cpp
#include "io.h"

std::mt19937_64 rng(42);
iOParams params = iOParams::toy(FieldPrimes::SMALL_PRIME);
iO io(params, rng);

// Load your miter circuit
Circuit C = load_blif("miter.blif");

// Obfuscate
auto obf = io.obfuscate(C, rng);

// Evaluate on input x
std::vector<uint8_t> x = {1, 0, 1, 0, 0, 0, 0, 0};
auto y = io.eval(obf, x, rng);
```

**Step 4 — Program the hardware circuit.**
Generate the SystemVerilog decryption circuit and program the obfuscated
functional keys (`SK_C`) into it. The hardware circuit takes `x`, fresh
randomness, and `SK_C` as inputs and outputs `C(x)` without ever revealing
the key or the original circuit structure.

```bash
python3 gen_fe_verilog.py   # edit parameters at the top for your security level
# synthesise with: yosys -p "read_verilog -sv fe_dec_enc.v; synth -top fe_dec_enc"
```

The resulting synthesised netlist is your locked circuit. It is parametrised
by `SK_C` (loaded at activation time from tamper-proof memory) and computes
the correct output only when `SK_C` is present.

---

## Security notes

- **Prototype security level**: toy parameters (`K_LPN=3`, `W=64`) are for
  correctness demonstration only — not for any security claim.
- **~80-bit iO security**: set `K_LPN=160`, `N_PRIME=160`, `W=128`, `P_BITS=127` in
  `gen_fe_verilog.py`.
- **Pairing security**: the simulated pairing in prototype mode is a field
  multiplication — it provides zero hiding. Compile with `-DPHFE_USE_MCL`
  for any security-sensitive use.

---

## References

- Jain, Lin, Sahai. *Indistinguishability Obfuscation from Well-Founded
  Assumptions.* STOC 2022. (JLS22)
- El Massad et al. *Locked Circuit Indistinguishability.* IEEE CSF 2022.
- Goldreich, Goldwasser, Micali. *How to construct random functions.* JACM
  1986. (GGM PRF)

---

## License

All rights reserved. See LICENSE for details.
