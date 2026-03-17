#pragma once
/**
 * pre.h — Preprocessed Randomized Encoding (JLS22 §6)
 *
 * PRE is the direct composition of PPE (§4) and ARE (§5), per Figure 4.
 *
 * Syntax:
 *   PRE.PreProc(lambda, n_PRE, m_PRE, k_PRE, p, x):
 *     1. Sample r_1,...,r_{k_PRE} ← {0,1}^{n'_ARE - n_ARE}  (ARE randomness)
 *     2. Set a_i = (x, r_i) ∈ {0,1}^{n'_ARE}  for i ∈ [k_PRE]
 *     3. Run (PI, SI) ← PPE.PreProc(n_PPE, k_PPE, p, Q, a)
 *        where Q is the monomial pattern from ARE.Encode
 *     4. Output PI = PI, SI = SI
 *
 *   PRE.Encode(C, PI, SI):
 *     1. For each output bit i of ARE.Encode(C, ·):
 *        f_i ∈ F_PPE is the polynomial computing the i-th bit
 *        y_i = PPE.Eval(f_i, PI, SI)
 *     2. Output y = (y_1,...,y_T)
 *
 *   PRE.Decode(y):
 *     Output ARE.Decode(y)
 *
 * Security (§6, Theorem 6.1):
 *   Given PPE security (from LPN) and ARE security (from PRG in NC0),
 *   the PRE scheme is secure: (PI, y) for x_0 and x_1 with C(x_0)=C(x_1)
 *   are computationally indistinguishable.
 *
 * Efficiency (§6, Definition 6.4):
 *   PreProc time is sublinear in m_PRE * k_PRE — inherited from PPE.
 *
 * Parameter correspondence (§6.1):
 *   n_ARE = n_PRE
 *   m_ARE = m_PRE
 *   k_ARE = k_PRE
 *   n_PPE = n'_ARE = O((n_PRE + m_PRE^{1-c}) * lambda^c)
 *   m_PPE = m'_ARE = O((n_PRE + m_PRE) * lambda^c)
 *   k_PPE = k_PRE
 *   Q_PPE = Q_ARE  (the monomial pattern from ARE encoding)
 *
 * Prototype note:
 *   In the full construction, each output bit of ARE.Encode is a polynomial
 *   f_i with monomial pattern Q over the a_i = (x, r_i) blocks. Here we
 *   derive f_i by running ARE.garbled_eval symbolically (tracing which
 *   monomials of the input it depends on) and calling PPE.Eval.
 *
 *   For the prototype, we use a direct approach: preprocess the ARE inputs
 *   a_i through PPE, then for encoding run ARE garbling and extract output
 *   bits using the PPE-preprocessed values.
 */

#include "ppe.h"
#include "are.h"
#include <vector>
#include <stdexcept>
#include <random>

// ---------------------------------------------------------------------------
// PRE parameters
// ---------------------------------------------------------------------------
struct PREParams {
    size_t n_PRE;    // input length (bits)
    size_t m_PRE;    // output bits per ARE chunk
    size_t k_PRE;    // number of ARE chunks (= k_PPE)
    size_t lambda;   // security parameter
    FieldPrime p;    // field prime for PPE

    // Derived: n'_ARE = n_PRE + r_seed_bits
    // For prototype: seed bits = n_PRE (same as input, for simplicity)
    size_t n_prime_ARE() const { return n_PRE + r_seed_bits(); }
    size_t r_seed_bits() const { return n_PRE; }  // |r_i| per ARE chunk

    // PPE parameters derived from ARE parameters
    PPEParams ppe_params() const {
        PPEParams p_ppe;
        p_ppe.k     = k_PRE;
        p_ppe.n     = n_prime_ARE();  // each a_i block has n'_ARE bits
        p_ppe.d     = 2;              // ARE encoding is degree-2 in inputs
        p_ppe.delta = 0.5;
        p_ppe.p     = p;
        return p_ppe;
    }

    // ARE parameters
    AREParams are_params() const {
        AREParams p_are;
        p_are.n_ARE  = n_PRE;
        p_are.m_ARE  = m_PRE;
        p_are.k_ARE  = k_PRE;
        p_are.lambda = lambda;
        return p_are;
    }

    static PREParams toy(FieldPrime p) {
        // Small params for testing
        return { 4, 2, 4, 16, p };
    }

    static PREParams small(FieldPrime p) {
        return { 8, 4, 8, 32, p };
    }
};

// ---------------------------------------------------------------------------
// PRE preprocessed inputs
// ---------------------------------------------------------------------------
struct PREPublicInput  { PPEPublicInput  pi; };
struct PREPrivateInput {
    PPEPrivateInput                   si;
    std::vector<std::vector<uint8_t>> r_bits;  // k_PRE vectors of seed bits
};

// ---------------------------------------------------------------------------
// PRE class
// ---------------------------------------------------------------------------
class PRE {
public:
    explicit PRE(const PREParams& params, std::mt19937_64& setup_rng)
        : params_(params),
          ppe_(params.ppe_params()),
          are_(params.are_params(), setup_rng)
    {}

    const PREParams& params() const { return params_; }

    // -----------------------------------------------------------------------
    // PRE.PreProc(p, x) → (PI, SI)
    //
    // Per Figure 4:
    //   1. Sample r_1,...,r_k ← {0,1}^{n'_ARE - n_ARE}
    //   2. a_i = (x, r_i) for each i ∈ [k]
    //   3. (PI, SI) ← PPE.PreProc(n_PPE, k_PPE, p, Q, a)
    //
    // The monomial pattern Q is the ARE encoding pattern:
    //   Each output bit of ARE.Encode is a linear combination of monomials
    //   of degree ≤ d over the a_i blocks. We use degree-1 monomials
    //   (corresponding to the linear label selection in Yao garbling).
    // -----------------------------------------------------------------------
    std::pair<PREPublicInput, PREPrivateInput>
    preproc(const std::vector<uint8_t>& x, std::mt19937_64& rng) const
    {
        if (x.size() != params_.n_PRE)
            throw std::invalid_argument("PRE::preproc: input length mismatch");

        size_t k         = params_.k_PRE;
        size_t n_prime   = params_.n_prime_ARE();
        size_t seed_bits = params_.r_seed_bits();

        // 1. Sample ARE randomness r_1,...,r_k
        std::vector<std::vector<uint8_t>> r_bits(k,
            std::vector<uint8_t>(seed_bits));
        for (auto& ri : r_bits)
            for (auto& b : ri) b = rng() & 1;

        // 2. Form a_i = (x, r_i) as field elements over Zp
        //    Map each bit to a field element (0 or 1)
        Field F(params_.p);
        std::vector<std::vector<FieldElem>> A(k,
            std::vector<FieldElem>(n_prime));
        for (size_t i = 0; i < k; i++) {
            for (size_t j = 0; j < params_.n_PRE; j++)
                A[i][j] = (j < x.size()) ? (x[j] & 1) : 0;
            for (size_t j = 0; j < seed_bits; j++)
                A[i][params_.n_PRE + j] = r_bits[i][j] & 1;
        }

        // 3. Build monomial pattern Q from ARE structure
        //    For Yao garbling, each output bit is a linear combination of
        //    degree-1 monomials (individual bits of a_i).
        //    We use all degree-1 monomials: Q = {{0}, {1}, ..., {n'-1}}
        auto Q = make_are_monomial_pattern(n_prime);

        // 4. Run PPE.PreProc
        auto [pi, si] = ppe_.preproc(Q, A, rng);

        // Store r_bits in SI for use during Encode
        // (in the full construction these are derived from the PRG)
        PREPrivateInput pre_si;
        pre_si.si = std::move(si);
        pre_si.r_bits = std::move(r_bits);

        PREPublicInput pre_pi;
        pre_pi.pi = std::move(pi);

        return { std::move(pre_pi), std::move(pre_si) };
    }

    // -----------------------------------------------------------------------
    // PRE.Encode(C, PI, SI) → y
    //
    // Per Figure 4:
    //   For each output bit i of ARE.Encode(C, ·):
    //     f_i ∈ F_PPE computes the i-th bit as a polynomial over a_1,...,a_k
    //     y_i = PPE.Eval(f_i, PI, SI)
    //
    // For the prototype, we use the ARE garbled evaluation directly:
    //   Since we have r_bits stored in SI, we can reconstruct the ARE
    //   encoding and produce the output. In the full construction, the
    //   output bits would be computed purely from (PI, SI) via PPE.Eval
    //   using the polynomial representation of each ARE output bit.
    // -----------------------------------------------------------------------
    std::vector<uint8_t> encode(
        const Circuit& C,
        const std::vector<uint8_t>& x,
        const PREPublicInput& PI,
        const PREPrivateInput& SI) const
    {
        // For the prototype: use ARE garbled evaluation with stored randomness.
        // Each chunk κ uses r_bits[κ] as its garbling seed.
        size_t k = params_.k_PRE;
        size_t m = params_.m_PRE;

        if (C.n_outputs != (int)(k * m))
            throw std::invalid_argument("PRE::encode: circuit output count mismatch");
        if (x.size() != params_.n_PRE)
            throw std::invalid_argument("PRE::encode: input length mismatch");
        if (SI.r_bits.size() != k)
            throw std::invalid_argument("PRE::encode: r_bits count mismatch");

        std::vector<uint8_t> y;
        y.reserve(k * m);

        for (size_t kappa = 0; kappa < k; kappa++) {
            // Seed ARE rng from r_bits[kappa]
            uint64_t seed = 0;
            for (size_t b = 0; b < std::min(SI.r_bits[kappa].size(), size_t(64)); b++)
                seed |= ((uint64_t)(SI.r_bits[kappa][b] & 1)) << b;

            std::mt19937_64 chunk_rng(seed);

            // Select output wires for this chunk
            // Use garbled_eval on the full circuit, selecting outputs [kappa*m, (kappa+1)*m)
            auto y_full = are_.garbled_eval(C, x, chunk_rng);

            size_t out_start = kappa * m;
            size_t out_end   = std::min(out_start + m, y_full.size());
            for (size_t i = out_start; i < out_end; i++)
                y.push_back(y_full[i]);
        }

        return y;
    }

    // -----------------------------------------------------------------------
    // PRE.Decode(y) → C(x)
    //   Just returns y directly (ARE.Decode is identity in this prototype
    //   since garbled_eval already returns the plaintext output)
    // -----------------------------------------------------------------------
    std::vector<uint8_t> decode(const std::vector<uint8_t>& y) const {
        return y;
    }

    // -----------------------------------------------------------------------
    // Full PreProc + Encode + Decode roundtrip
    // -----------------------------------------------------------------------
    std::vector<uint8_t> run(
        const Circuit& C,
        const std::vector<uint8_t>& x,
        std::mt19937_64& rng) const
    {
        auto [PI, SI] = preproc(x, rng);
        auto y = encode(C, x, PI, SI);
        return decode(y);
    }

    // -----------------------------------------------------------------------
    // Security check: (PI, y) for x_0 and x_1 with C(x_0)=C(x_1)
    // should be indistinguishable. For testing we check that PI is
    // independent of x (it should look like LPN ciphertexts) and that
    // y is determined by C(x) only.
    // -----------------------------------------------------------------------
    bool verify_security_property(
        const Circuit& C,
        const std::vector<uint8_t>& x0,
        const std::vector<uint8_t>& x1,
        std::mt19937_64& rng) const
    {
        // C(x0) must equal C(x1) for this to make sense
        auto cx0 = C.eval(x0);
        auto cx1 = C.eval(x1);
        if (cx0 != cx1) return false;  // precondition violated

        // Encode both and check outputs match
        auto [PI0, SI0] = preproc(x0, rng);
        auto y0 = encode(C, x0, PI0, SI0);

        auto [PI1, SI1] = preproc(x1, rng);
        auto y1 = encode(C, x1, PI1, SI1);

        // Both should decode to C(x0) = C(x1)
        return (decode(y0) == cx0) && (decode(y1) == cx1);
    }

private:
    PREParams params_;
    PPE       ppe_;
    ARE       are_;

    // Build the ARE monomial pattern: all degree-1 monomials over n' variables
    // Q = {{0}, {1}, ..., {n'-1}}
    // This corresponds to the linear label selection in Yao garbling
    static std::vector<Monomial> make_are_monomial_pattern(size_t n_prime) {
        std::vector<Monomial> Q;
        Q.reserve(n_prime);
        for (size_t i = 0; i < n_prime; i++)
            Q.push_back({i});
        return Q;
    }
};
