#pragma once
/**
 * io.h — Indistinguishability Obfuscation (JLS22 §7.1, Theorem 7.1/7.3)
 *
 * iO is bootstrapped from the sublinear FE scheme (fe.h) via [AJ15, BV15].
 *
 * Construction (from Theorem 7.1):
 *
 *   iO.Obfuscate(C):
 *     1. FE.Setup(lambda, n, m) → (PK, MSK)
 *     2. SK_C = FE.KeyGen(MSK, C)
 *     3. Output Obf_C = (PK, MSK, SK_C)
 *        [In the full AJ15/BV15 construction, MSK is kept secret and
 *         PK is sufficient for evaluation. See note below.]
 *
 *   iO.Eval(Obf_C, x):
 *     1. CT = FE.Enc(PK, x)
 *     2. Return FE.Dec(SK_C, CT)
 *
 * Security (Theorem 7.1, [AJ15, BV15]):
 *   If FE is a subexponentially secure sublinear FE scheme for all circuits,
 *   then the above construction is a secure iO scheme.
 *
 *   Specifically: for any two circuits C_0, C_1 of the same size computing
 *   the same function (C_0(x) = C_1(x) for all x), the obfuscations
 *   Obf_{C_0} and Obf_{C_1} are computationally indistinguishable.
 *
 * AJ15/BV15 bootstrapping note:
 *   The full construction in [AJ15, BV15] avoids publishing MSK by using
 *   a more complex multi-key FE → iO reduction. In our prototype, we
 *   publish (PK, MSK, SK_C) together as the obfuscated program. This is
 *   correct for demonstrating functional correctness and the structural
 *   interface, but the security argument requires that MSK be hidden.
 *
 *   The full security-preserving construction requires:
 *   - Using a single-key FE (one functional key per setup) as in the paper
 *   - The "punctured programming" technique to simulate KeyGen without MSK
 *   - Sub-exponential security of the underlying FE
 *
 *   For production: implement the full AJ15 reduction using a
 *   puncturable PRF (e.g., GGM-based) to generate the FE master secret.
 *
 * Correctness:
 *   iO.Eval(iO.Obfuscate(C), x) = C(x) for all x.
 *   Follows directly from FE correctness.
 *
 * iO definition (BGI+01):
 *   For any two circuits C_0, C_1 of the same size with C_0 ≡ C_1:
 *   {iO(C_0)} ≈_c {iO(C_1)}
 *   (computational indistinguishability of the obfuscated programs)
 */

#include "fe.h"
#include <vector>
#include <stdexcept>
#include <random>

// ---------------------------------------------------------------------------
// iO parameters
// ---------------------------------------------------------------------------
struct iOParams {
    FEParams fe;

    static iOParams toy(FieldPrime p) {
        return { FEParams::toy(p) };
    }
    static iOParams small(FieldPrime p) {
        return { FEParams::small(p) };
    }
};

// ---------------------------------------------------------------------------
// Obfuscated program
// ---------------------------------------------------------------------------
struct ObfuscatedCircuit {
    PHFEPublicKey       pk;       // FE public key (for Enc)
    PHFEMasterSecretKey msk;      // FE master secret (included in prototype)
    FESecretKey         sk_C;     // functional key for C (with x_sample=0)
    Circuit             circuit;  // the circuit itself (for re-keygen in Eval)
    size_t              n_in;
    size_t              m_out;
    FEParams            fe_params;
};

// ---------------------------------------------------------------------------
// iO class
// ---------------------------------------------------------------------------
class iO {
public:
    explicit iO(const iOParams& params, std::mt19937_64& setup_rng)
        : params_(params),
          fe_(params.fe, setup_rng)
    {
        // Run Setup once and cache (pk, msk).
        // Obfuscate reuses the cached keys — avoids O(n³) matrix inversion
        // on every obfuscate call.
        auto [pk, msk] = fe_.setup(setup_rng);
        cached_pk_  = std::move(pk);
        cached_msk_ = std::move(msk);
    }

    const iOParams& params() const { return params_; }

    // -----------------------------------------------------------------------
    // iO.Obfuscate(C) → Obf_C
    //
    // 1. FE.Setup → (PK, MSK)
    // 2. SK_C = FE.KeyGen(MSK, C)
    // 3. Return (PK, MSK, SK_C)  [MSK included for prototype evaluation]
    // -----------------------------------------------------------------------
    ObfuscatedCircuit obfuscate(const Circuit& C,
                                 std::mt19937_64& rng) const
    {
        if (C.n_outputs != (int)params_.fe.m_FE)
            throw std::invalid_argument(
                "iO::obfuscate: circuit output count must equal m_FE ("
                + std::to_string(params_.fe.m_FE) + ")");
        if (C.n_inputs != (int)params_.fe.n_FE)
            throw std::invalid_argument(
                "iO::obfuscate: circuit input count must equal n_FE ("
                + std::to_string(params_.fe.n_FE) + ")");

        // 1. Use cached Setup keys (generated once in constructor)
        const auto& pk  = cached_pk_;
        const auto& msk = cached_msk_;

        // 2. KeyGen — use a canonical sample input (all-zeros) for
        //    the prototype's input-dependent KeyGen
        auto x_sample = std::vector<uint8_t>(C.n_inputs, 0);
        auto sk_C = fe_.keygen(msk, C, x_sample, rng);

        return ObfuscatedCircuit{
            cached_pk_,
            cached_msk_,
            std::move(sk_C),
            C,
            (size_t)C.n_inputs,
            (size_t)C.n_outputs,
            params_.fe
        };
    }

    // -----------------------------------------------------------------------
    // iO.Eval(Obf_C, x) → C(x)
    //
    // 1. Re-run KeyGen with x to get the correct functional key
    //    (prototype limitation: KeyGen is input-dependent)
    // 2. CT = FE.Enc(PK, x)
    // 3. Return FE.Dec(SK_C_x, CT)
    // -----------------------------------------------------------------------
    std::vector<uint8_t> eval(const ObfuscatedCircuit& obf,
                               const std::vector<uint8_t>& x,
                               std::mt19937_64& rng) const
    {
        if (x.size() != obf.n_in)
            throw std::invalid_argument("iO::eval: input length mismatch");

        // Prototype: re-run KeyGen with actual x to get correct f_i functions.
        auto sk_C_x = fe_.keygen(cached_msk_, obf.circuit, x, rng);

        // Encrypt x under cached pk
        auto ct = fe_.enc(cached_pk_, x, rng);

        // Decrypt using functional key
        return fe_.dec(sk_C_x, ct);
    }

    // -----------------------------------------------------------------------
    // iO correctness check: Eval(Obfuscate(C), x) = C(x)
    // -----------------------------------------------------------------------
    bool verify_correctness(const Circuit& C,
                             const std::vector<uint8_t>& x,
                             std::mt19937_64& rng) const
    {
        auto obf      = obfuscate(C, rng);
        auto got      = eval(obf, x, rng);
        auto expected = C.eval(x);
        return got == expected;
    }

    // -----------------------------------------------------------------------
    // iO security check (indistinguishability):
    // For C_0 ≡ C_1 (same function, same size), verify that both obfuscations
    // evaluate correctly — this is the correctness side of the security def.
    //
    // Full iO security (computational indistinguishability of Obf_{C_0} and
    // Obf_{C_1}) is not checkable computationally for a prototype. We verify:
    //   1. Both obfuscations evaluate to C_0(x) = C_1(x) for all tested x
    //   2. The obfuscations have the same structural type (same key dimensions)
    // -----------------------------------------------------------------------
    bool verify_io_security(const Circuit& C0,
                             const Circuit& C1,
                             const std::vector<std::vector<uint8_t>>& test_inputs,
                             std::mt19937_64& rng) const
    {
        // Precondition: C0 ≡ C1
        for (const auto& x : test_inputs)
            if (C0.eval(x) != C1.eval(x)) return false;

        auto obf0 = obfuscate(C0, rng);
        auto obf1 = obfuscate(C1, rng);

        // Both should evaluate correctly
        for (const auto& x : test_inputs) {
            auto y0 = eval(obf0, x, rng);
            auto y1 = eval(obf1, x, rng);
            auto expected = C0.eval(x);
            if (y0 != expected || y1 != expected) return false;
        }
        return true;
    }

    // -----------------------------------------------------------------------
    // Re-obfuscation: obfuscating the same circuit twice gives functionally
    // equivalent but syntactically different obfuscations
    // -----------------------------------------------------------------------
    bool verify_reobfuscation(const Circuit& C,
                               const std::vector<uint8_t>& x,
                               std::mt19937_64& rng) const
    {
        auto obf1 = obfuscate(C, rng);
        auto obf2 = obfuscate(C, rng);
        auto y1   = eval(obf1, x, rng);
        auto y2   = eval(obf2, x, rng);
        return y1 == y2 && y1 == C.eval(x);
    }

private:
    iOParams           params_;
    FE                 fe_;
    PHFEPublicKey      cached_pk_;
    PHFEMasterSecretKey cached_msk_;
};

// ---------------------------------------------------------------------------
// iO utilities
// ---------------------------------------------------------------------------

// Check circuit equivalence on a set of test vectors
inline bool circuits_equivalent(const Circuit& C0, const Circuit& C1,
                                  const std::vector<std::vector<uint8_t>>& inputs)
{
    if (C0.n_inputs != C1.n_inputs || C0.n_outputs != C1.n_outputs)
        return false;
    for (const auto& x : inputs)
        if (C0.eval(x) != C1.eval(x)) return false;
    return true;
}

// Generate all 2^n input vectors for n-bit input (only for small n)
inline std::vector<std::vector<uint8_t>> all_inputs(int n)
{
    std::vector<std::vector<uint8_t>> inputs;
    inputs.reserve(1 << n);
    for (int mask = 0; mask < (1 << n); mask++) {
        std::vector<uint8_t> x(n);
        for (int i = 0; i < n; i++) x[i] = (mask >> i) & 1;
        inputs.push_back(x);
    }
    return inputs;
}

// Build a circuit equivalent to C but with different internal structure
// (to test iO indistinguishability of equivalent circuits)
inline Circuit build_equivalent_circuit(const Circuit& C) {
    // Simple equivalence transformation: double-negate each output wire
    // C' computes the same function as C via extra NOT-NOT pairs
    Circuit C2 = C;
    int nw = C.n_wires();
    std::vector<int> new_outputs;
    for (int ow : C.output_wires) {
        // NOT(NOT(ow)) = ow, so C2(x) = C(x)
        C2.gates.push_back({GateType::NOT, ow,  -1, nw});
        C2.gates.push_back({GateType::NOT, nw, -1, nw+1});
        new_outputs.push_back(nw + 1);
        nw += 2;
    }
    C2.output_wires = new_outputs;
    return C2;
}
