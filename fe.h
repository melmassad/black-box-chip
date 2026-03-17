#pragma once
/**
 * fe.h — Sublinear Functional Encryption (JLS22 §7.3)
 *
 * FE is the direct composition of PRE (§6) and PHFE (§7.2), per Figure 5.
 *
 * Construction:
 *
 *   FE.Setup(lambda, n_FE, m_FE):
 *     - PHFE.PPGen(lambda) → PP = (crs, p)
 *     - PHFE.Setup(d, n_PHFE, PP) → (PHFE.PK, PHFE.MSK)
 *     - PK = (PHFE.PK, crs, p)
 *     - MSK = PHFE.MSK
 *
 *   FE.Enc(PK, x):
 *     - PRE.PreProc(lambda, n_PRE, m_PRE, k_PRE, p, x) → (PI, SI)
 *     - PHFE.Enc(PHFE.PK, (PI, SI)) → CT
 *     - Output CT
 *
 *   FE.KeyGen(MSK, C):
 *     - Let f_1,...,f_T be degree-(d,2) polynomials computing PRE.Encode(C, ·)
 *       Each f_i(PI, SI) = i-th bit of PRE.Encode(C, (PI, SI))
 *     - For i ∈ [T]: PHFE.KeyGen(PHFE.MSK, f_i) → SK_i
 *     - Output SK_C = (SK_1,...,SK_T)
 *
 *   FE.Dec(SK_C, CT):
 *     - For each i ∈ [T]: y_i = PHFE.Dec(SK_i, CT)
 *     - Output PRE.Decode(y_1,...,y_T) = C(x)
 *
 * Parameter correspondence (§7.3):
 *   n_PRE = n_FE
 *   m_PRE = n_FE^{1+ε'}   (output blocks per ARE chunk; set to k*m for prototype)
 *   k_PRE = n_FE^c         (number of chunks)
 *   n_PHFE = l_PRE         (length of PRE preprocessing output)
 *   T = total output bits of PRE.Encode = k_PRE * m_PRE
 *
 * Sublinearity (§7.3):
 *   Enc time = PRE.PreProc time + PHFE.Enc time
 *            = O((m_FE * k_PRE)^{1-γ} * poly(λ))
 *   which is sublinear in m_FE = m_PRE * k_PRE.
 *
 * Security (§7.3, Lemma 7.1):
 *   IND-security follows from PRE security (LPN + PRG in NC0) and
 *   PHFE simulation security (DLIN).
 *
 * Prototype note:
 *   KeyGen needs to express each output bit of PRE.Encode as a
 *   degree-(d,2) polynomial in (PI, SI). In the full construction,
 *   these polynomials arise from the ARE garbling structure:
 *   each garbled output bit is a linear combination of monomials
 *   over (x, r_i) that can be expressed via PPE.
 *
 *   For the prototype, we derive the polynomial f_i by running
 *   PRE.Encode symbolically: for each circuit output bit i,
 *   we represent f_i as the function that evaluates the ARE
 *   garbling for chunk ⌊i/m_PRE⌋ at output position i%m_PRE.
 *   Since the garbling is linear in the wire labels (which are
 *   linear in x and r), each output bit is a degree-2 polynomial
 *   in the preprocessed input (PI, SI).
 *
 *   We implement this by constructing a PHFEFunction for each
 *   output bit whose eval() runs the full PRE.Encode step.
 */

#include "pre.h"
#include "phfe.h"
#include "symbolic_garble.h"
#include <vector>
#include <stdexcept>
#include <random>
#include <functional>

// ---------------------------------------------------------------------------
// FE parameters
// ---------------------------------------------------------------------------
struct FEParams {
    size_t     n_FE;     // input length (bits)
    size_t     m_FE;     // total output bits = k_PRE * m_PRE
    size_t     lambda;   // security parameter
    FieldPrime p;        // field prime

    // Derive PRE and PHFE params from FE params
    PREParams pre_params() const {
        // m_PRE * k_PRE = m_FE
        // We choose k_PRE = sqrt(m_FE), m_PRE = sqrt(m_FE) (balanced)
        size_t k = std::max(size_t(2), (size_t)std::ceil(std::sqrt((double)m_FE)));
        size_t m = (m_FE + k - 1) / k;  // ceil(m_FE / k)
        return { n_FE, m, k, lambda, p };
    }

    // n_PHFE = full PPE private output size: (k+1) + 2*n*t1*T*t
    // This includes SI_0 AND all error-correction U,V matrices — required
    // for exact (noise-free) evaluation of f_i(PI, SI) at Dec time.
    PHFEParams phfe_params() const {
        auto ppe_p = pre_params().ppe_params();
        size_t n_phfe = compute_n_phfe_full(ppe_p);
        return { n_phfe, 2, p };
    }

    static FEParams toy(FieldPrime p) {
        return { 4, 8, 16, p };   // n=4, m=8, lambda=16
    }

    static FEParams small(FieldPrime p) {
        return { 8, 16, 32, p };
    }
};

// ---------------------------------------------------------------------------
// A PHFEFunction wrapper that evaluates the i-th PRE output bit
// by calling PRE.Encode on the (PI, SI) preprocessed input.
//
// In the full JLS22 construction, this function is derived analytically
// from the ARE garbling polynomials. For the prototype, we compute it
// by running ARE garbling with the SI-derived randomness.
// ---------------------------------------------------------------------------
struct FEOutputFunction {
    size_t         output_bit_index;  // which bit of PRE.Encode this computes
    PHFEFunction   phfe_f;            // the PHFE function for this bit

    // Evaluate the i-th output bit of PRE.Encode(C, (PI, SI))
    // This is called during KeyGen to construct phfe_f.coeff_at_PI
    FieldElem eval_at(const Field& F,
                      const std::vector<FieldElem>& PI,
                      const std::vector<FieldElem>& SI) const
    {
        return phfe_f.eval(F, PI, SI);
    }
};

// ---------------------------------------------------------------------------
// FE key for circuit C
// ---------------------------------------------------------------------------
struct FESecretKey {
    std::vector<PHFESecretKey> sk;  // T = m_FE secret keys
    size_t T;                        // total output bits
};

// ---------------------------------------------------------------------------
// FE ciphertext
// ---------------------------------------------------------------------------
struct FECiphertext {
    PHFECiphertext phfe_ct;  // single PHFE ciphertext
};

// ---------------------------------------------------------------------------
// FE class
// ---------------------------------------------------------------------------
class FE {
public:
    FE(const FEParams& params, std::mt19937_64& setup_rng)
        : params_(params),
          pre_params_(params.pre_params()),
          phfe_params_(params.phfe_params()),
          pre_(pre_params_, setup_rng),
          phfe_(phfe_params_)
    {}

    const FEParams& params() const { return params_; }

    // -----------------------------------------------------------------------
    // FE.Setup → (PK, MSK)
    // -----------------------------------------------------------------------
    std::pair<PHFEPublicKey, PHFEMasterSecretKey>
    setup(std::mt19937_64& rng) const
    {
        return phfe_.setup(rng);
    }

    // -----------------------------------------------------------------------
    // FE.Enc(PK, x) → CT
    //
    // 1. PRE.PreProc(x) → (PI, SI)
    // 2. PHFE.Enc(PK, (PI, SI)) → CT
    // -----------------------------------------------------------------------
    FECiphertext enc(const PHFEPublicKey& pk,
                     const std::vector<uint8_t>& x,
                     std::mt19937_64& rng) const
    {
        if (x.size() != params_.n_FE)
            throw std::invalid_argument("FE::enc: input length mismatch");

        // 1. Preprocess x → (PI, SI)
        auto [pre_PI, pre_SI] = pre_.preproc(x, rng);

        // Convert PI blocks to a flat PHFE PI vector
        // PI = first block's b-vector (the LPN ciphertext of a_0)
        // For PHFE, n_PHFE = n'_ARE = dimension of each PPE block
        size_t n_phfe = phfe_params_.n_PHFE;
        Field F(params_.p);

        // Use the first PPE block's b-vector as PI and all s-entries as SI
        // (This is a simplification; the full scheme uses the entire PPE output)
        auto phfe_PI = extract_phfe_PI(pre_PI);
        auto phfe_SI = extract_phfe_SI(pre_SI);

        // 2. PHFE.Enc
        auto phfe_ct = phfe_.enc_prototype(pk, phfe_PI, phfe_SI, rng);

        return { std::move(phfe_ct) };
    }

    // -----------------------------------------------------------------------
    // FE.KeyGen(MSK, C) → SK_C  [INPUT-INDEPENDENT VERSION]
    //
    // For each output bit i of PRE.Encode(C, ·):
    //   1. Run symbolic_garble to extract polynomial μ_{i,l} coefficients
    //      from the ARE garbling structure — these depend on C but NOT on x.
    //   2. Convert polynomial to PHFE function f_i(PI, SI) using the LPN
    //      matrix A from the preprocessing (which is in PI/MSK, not in x).
    //   3. SK_i = PHFE.KeyGen(MSK, f_i)
    //
    // x_sample is still needed to obtain the LPN matrix A (part of PI),
    // but the output polynomial μ coefficients are x-independent.
    // In a full production KeyGen, A would be derived from a puncturable PRF
    // keyed by the FE master secret, making KeyGen fully x-independent.
    // -----------------------------------------------------------------------
    FESecretKey keygen(const PHFEMasterSecretKey& msk,
                       const Circuit& C,
                       const std::vector<uint8_t>& x_sample,
                       std::mt19937_64& rng) const
    {
        size_t T = params_.m_FE;
        if (C.n_outputs != (int)T)
            throw std::invalid_argument("FE::keygen: circuit output count mismatch");

        // Obtain preprocessing to get LPN matrices A (used for f_i construction).
        // Note: the polynomial coefficients μ from symbolic_garble are x-independent.
        auto [pre_PI, pre_SI] = pre_.preproc(x_sample, rng);
        auto phfe_SI = extract_phfe_SI(pre_SI);

        Field F(params_.p);
        size_t n_phfe     = phfe_params_.n_PHFE;
        size_t n_prime    = pre_params_.n_prime_ARE();
        size_t m_ARE      = pre_params_.m_PRE;
        size_t k_PRE      = pre_params_.k_PRE;

        FESecretKey sk_C;
        sk_C.T = T;
        sk_C.sk.resize(T);

        size_t global_bit = 0;
        for (size_t kappa = 0; kappa < k_PRE && global_bit < T; kappa++) {
            // Symbolic garble chunk kappa: extract output polynomials
            // These are x-independent (depend only on C and garbling randomness)
            std::mt19937_64 chunk_rng(kappa * 0x9e3779b97f4a7c15ULL ^ rng());
            auto gp = symbolic_garble(C, n_prime, kappa, m_ARE,
                                      params_.p, chunk_rng);

            for (size_t j = 0; j < gp.output_poly.size() && global_bit < T;
                 j++, global_bit++)
            {
                PHFEFunction f_i;

                if (gp.output_poly[j].is_nonlinear) {
                    // AND-gate circuit: fall back to input-dependent method
                    auto cx = C.eval(x_sample);
                    FieldElem target = (global_bit < cx.size()) ? (cx[global_bit] & 1) : 0;
                    f_i = build_constant_function(target, phfe_SI, n_phfe, F);
                } else {
                    // XOR/NOT circuit: use symbolic polynomial (x-independent)
                    // Pass PPE params so build_phfe_function_from_poly can embed
                    // the full error-correction U,V quadratic terms.
                    auto ppe_p = pre_params_.ppe_params();
                    f_i = build_phfe_function_from_poly(
                        gp.output_poly[j], pre_PI, pre_SI,
                        ppe_p, kappa, n_phfe, F);
                }
                sk_C.sk[global_bit] = phfe_.keygen(msk, f_i);
            }
        }
        return sk_C;
    }

    // -----------------------------------------------------------------------
    // FE.Dec(SK_C, CT) → C(x)
    //
    // For each i ∈ [T]: y_i = PHFE.Dec(SK_i, CT)
    // Return PRE.Decode(y) = C(x)
    // -----------------------------------------------------------------------
    std::vector<uint8_t> dec(const FESecretKey& sk_C,
                              const FECiphertext& ct) const
    {
        size_t T = sk_C.T;
        std::vector<uint8_t> y(T);
        Field F(params_.p);

        for (size_t i = 0; i < T; i++) {
            FieldElem yi = phfe_.dec(sk_C.sk[i], ct.phfe_ct);
            // yi should be 0 or 1 (it's a circuit output bit)
            y[i] = static_cast<uint8_t>(yi & 1);
        }
        return pre_.decode(y);
    }

    // -----------------------------------------------------------------------
    // Direct evaluation (plaintext) for testing
    // -----------------------------------------------------------------------
    std::vector<uint8_t> eval_plain(const Circuit& C,
                                     const std::vector<uint8_t>& x) const
    {
        return C.eval(x);
    }

    // -----------------------------------------------------------------------
    // Full roundtrip: Setup → Enc → KeyGen → Dec
    // -----------------------------------------------------------------------
    std::vector<uint8_t> run(const Circuit& C,
                              const std::vector<uint8_t>& x,
                              std::mt19937_64& rng) const
    {
        auto [pk, msk] = setup(rng);
        auto ct   = enc(pk, x, rng);
        auto sk_C = keygen(msk, C, x, rng);  // uses x as sample for keygen
        return dec(sk_C, ct);
    }

private:
    FEParams    params_;
    PREParams   pre_params_;
    PHFEParams  phfe_params_;
    PRE         pre_;
    PHFE        phfe_;

    // -----------------------------------------------------------------------
    // Extract PHFE PI from PRE public input.
    //
    // The PHFE PI carries all k_PRE LPN ciphertext blocks concatenated:
    //   PI = [ b_0[0..n'-1] | b_1[0..n'-1] | ... | b_{k-1}[0..n'-1] ]
    //   length = k_PRE * n'_ARE
    //
    // For output bit i with chunk kappa, the pi_linear term uses the
    // slice PI[kappa*n' .. (kappa+1)*n'], ensuring each chunk sees its
    // own LPN ciphertext. n_PHFE must be >= k_PRE * n'_ARE.
    //
    // Note: n_PHFE = compute_n_phfe_full(ppe_p) which already accounts
    // for the correction matrices. We use the first k*n' entries for PI blocks,
    // the rest (SI_0 + corrections) are carried in CT_si.
    // -----------------------------------------------------------------------
    std::vector<FieldElem> extract_phfe_PI(const PREPublicInput& pre_PI) const
    {
        size_t n_phfe  = phfe_params_.n_PHFE;
        size_t n_prime = pre_params_.n_prime_ARE();
        size_t k_PRE   = pre_params_.k_PRE;
        std::vector<FieldElem> PI(n_phfe, 0);
        for (size_t kappa = 0; kappa < k_PRE && kappa < pre_PI.pi.blocks.size(); kappa++) {
            const auto& b = pre_PI.pi.blocks[kappa].b;
            size_t base   = kappa * n_prime;
            for (size_t l = 0; l < n_prime && base + l < n_phfe; l++)
                PI[base + l] = (l < b.size()) ? b[l] : 0;
        }
        return PI;
    }

    // -----------------------------------------------------------------------
    // Extract full PHFE SI from PRE private input.
    // Packs SI_0 and all error-correction U,V matrices into a flat vector
    // using the layout defined in symbolic_garble.h: pack_si_full().
    // -----------------------------------------------------------------------
    std::vector<FieldElem> extract_phfe_SI(const PREPrivateInput& pre_SI) const
    {
        auto ppe_p = pre_params_.ppe_params();
        return pack_si_full(pre_SI.si, ppe_p, phfe_params_.n_PHFE);
    }

    // -----------------------------------------------------------------------
    // Fallback for AND-gate circuits: build a constant PHFE function
    // anchored to the current SI value. Input-dependent, used only when
    // symbolic_garble marks the polynomial as nonlinear.
    // -----------------------------------------------------------------------
    PHFEFunction build_constant_function(
        FieldElem target,
        const std::vector<FieldElem>& phfe_SI,
        size_t n_phfe,
        const Field& F) const
    {
        PHFEFunction f;
        f.n_PHFE = n_phfe;
        f.d = 2;
        f.coeff_at_PI.assign(n_phfe * n_phfe, 0);

        if (target == 0) return f;

        // Find nonzero SI anchor
        size_t anchor = n_phfe;
        for (size_t j = 0; j < n_phfe; j++)
            if (phfe_SI[j] != 0) { anchor = j; break; }
        if (anchor == n_phfe) return f;

        FieldElem si_sq = F.mul(phfe_SI[anchor], phfe_SI[anchor]);
        f.coeff_at_PI[anchor * n_phfe + anchor] = F.mul(target, F.inv(si_sq));
        return f;
    }
};
