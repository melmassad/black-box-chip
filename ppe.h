#pragma once
/**
 * ppe.h — Preprocessed Polynomial Encoding (JLS22 §4)
 *
 * PPE allows preprocessing k input blocks x_1,...,x_k ∈ Zp^n into a public
 * part PI and private part SI such that any polynomial f with a fixed monomial
 * pattern Q can be evaluated as a degree-(O(d),2) polynomial in (PI, SI).
 *
 * Construction (§4.1):
 *
 *   PreProc(p, Q, x_1,...,x_k):
 *     1. Sample s ← Zp^k, errors e_{j,i} ~ Ber(k^{-δ}) for j∈[k], i∈[n]
 *     2. PI = (flag, {A_j, b_j = s·A_j + e_j + x_j}_{j∈[k]})
 *        where A_j ← Zp^{k×n}
 *     3. SI_0 = (1,s)^{⊗⌈d/2⌉}   [all degree-⌈d/2⌉ monomials of s]
 *     4. For each monomial Q_r ∈ Q:
 *        a. Corr_{r,j} = Mon_{Q_r}(x_j) - Mon_{Q_r}(x_j + e_j)  for j∈[k]
 *        b. Assign Corr_{r,j} into t1 matrices {M_{r,γ}} of size T×T
 *           via φ: [k] → [t1] × [T] × [T]
 *        c. Factorize each M_{r,γ} = U_{r,γ} · V_{r,γ}^T   (rank ≤ t)
 *        d. SI_r = {U_{r,γ}, V_{r,γ}}_{γ∈[t1]}
 *
 *   Eval(f, PI, SI):
 *     For each (r,j), compute w_{r,j}(PI,SI) = flag · (
 *       Mon_{Q_r}(b_j - s·A_j)  +  U_{r,φ(j)_1} · V_{r,φ(j)_1}[φ(j)_2, φ(j)_3]
 *     )
 *     Return Σ_{r,j} μ_{r,j} · w_{r,j}
 *
 * Parameters (from §4.1):
 *   t1 = ⌈k^{1-δ}⌉  (number of matrix buckets)
 *   T  = ⌈k^{δ/2}⌉  (matrix dimension)
 *   t  = k^{δ/10}    (max non-zeros per bucket, slack param t2 in paper)
 */

#include "field.h"
#include "lpn.h"
#include <vector>
#include <cmath>
#include <stdexcept>
#include <cassert>
#include <map>
#include <tuple>

// ---------------------------------------------------------------------------
// Monomial: a subset of variable indices {l_1,...,l_d} ⊆ [n]
// ---------------------------------------------------------------------------
using Monomial = std::vector<size_t>;  // sorted indices into x_j ∈ Zp^n

// Evaluate Mon_Q(v) = ∏_{l ∈ Q} v[l]  mod p
inline FieldElem eval_monomial(const Field& F,
                               const Monomial& Q,
                               const std::vector<FieldElem>& v)
{
    FieldElem result = F.one();
    for (size_t l : Q) {
        if (l >= v.size())
            throw std::out_of_range("eval_monomial: index out of range");
        result = F.mul(result, v[l]);
    }
    return result;
}

// ---------------------------------------------------------------------------
// PPE parameters
// ---------------------------------------------------------------------------
struct PPEParams {
    size_t k;       // number of input blocks
    size_t n;       // dimension of each input block x_j ∈ Zp^n
    size_t d;       // maximum monomial degree
    double delta;   // LPN noise exponent: noise rate = k^{-delta}
    FieldPrime p;   // field prime

    // Derived parameters from §4.1
    size_t t1() const { return static_cast<size_t>(std::ceil(std::pow(k, 1.0 - delta))); }
    size_t T()  const { return static_cast<size_t>(std::ceil(std::pow(k, delta / 2.0))); }
    // t is the per-bucket nonzero capacity (slack parameter t2 in paper).
    // Paper sets t = k^{delta/10}, which for small k is < 2 and causes flag=false
    // almost always. For the prototype we use max(t, T) so each T×T bucket can
    // hold up to T nonzeros — sufficient to demonstrate correctness.
    size_t t()  const {
        size_t t_theory = static_cast<size_t>(std::ceil(std::pow(k, delta / 10.0)));
        return std::max(t_theory, T());  // generous bound for small k
    }
    double noise_rate() const { return std::pow(static_cast<double>(k), -delta); }

    // Toy params for testing
    static PPEParams toy(FieldPrime p) {
        // k=64 gives t1=8, T=4, t=1, noise_rate~0.125 — tractable flag probability
        return { 64, 8, 2, 0.5, p };
    }
    static PPEParams small(FieldPrime p) {
        return { 128, 8, 2, 0.5, p };
    }
};

// ---------------------------------------------------------------------------
// The canonical injective map φ: [k] → [t1] × [T] × [T]  (§4.1)
// φ(j) = (j1, j2, j3) where j1 = j % t1, then j2,j3 from quotient
// ---------------------------------------------------------------------------
struct PhiResult { size_t j1, j2, j3; };

inline PhiResult phi_map(size_t j, size_t t1, size_t T) {
    size_t j1 = j % t1;
    size_t rem = j / t1;
    size_t j2 = rem % T;
    size_t j3 = rem / T;
    return { j1, j2, j3 };
}

// ---------------------------------------------------------------------------
// SI component for one monomial Q_r:
//   SI_r = { (U_{r,γ}, V_{r,γ}) }_{γ ∈ [t1]}
//   Each U, V is T × t (stored as vector of t column-vectors of length T)
// ---------------------------------------------------------------------------
struct SI_r {
    // U_gamma[l] = l-th column of U_{r,gamma}, length T
    // V_gamma[l] = l-th column of V_{r,gamma}, length T
    // M_{r,gamma} ≈ U_gamma * V_gamma^T
    std::vector<std::vector<FieldElem>> U;  // t1 matrices, each T×t (t columns of length T)
    std::vector<std::vector<FieldElem>> V;  // same shape
    size_t t1, T, t;

    // Access U_{gamma}[row] = column gamma*t + col?
    // We flatten: U[gamma] is a T*t vector (T rows, t cols, col-major)
    //             U[gamma][row * t + col] = U_{gamma}[row][col]
    // This matches: Corr = U_{j1} * V_{j1}^T,
    //   so Corr[j2,j3] = sum_l U_{j1}[j2,l] * V_{j1}[j3,l]
};

// ---------------------------------------------------------------------------
// Full PPE preprocessed input
// ---------------------------------------------------------------------------
struct PPEPublicInput {
    bool flag;
    // For each block j ∈ [k]: LPN encryption of x_j
    std::vector<LPNEncryption> blocks;  // size k
};

struct PPEPrivateInput {
    // SI_0 = (1,s)^{⊗⌈d/2⌉}: all degree-⌈d/2⌉ monomials in (1,s)
    // We store it as a flat vector of field elements
    std::vector<FieldElem> SI0;

    // SI_r for each monomial Q_r ∈ Q
    std::vector<SI_r> SI_monomials;  // size |Q|
};

// ---------------------------------------------------------------------------
// A polynomial f ∈ F_PPE: linear combination of monomials across blocks
//   f(x) = Σ_{r ∈ [m], j ∈ [k]} μ_{r,j} · Mon_{Q_r}(x_j)
// ---------------------------------------------------------------------------
struct PPEFunction {
    std::vector<Monomial>                       Q;      // monomial pattern, size m
    std::vector<std::vector<FieldElem>>         mu;     // mu[r][j] = coefficient, size m × k
};

// ---------------------------------------------------------------------------
// PPE class
// ---------------------------------------------------------------------------
class PPE {
public:
    explicit PPE(const PPEParams& params) : params_(params), F_(params.p) {}

    const PPEParams& params() const { return params_; }
    const Field&     field()  const { return F_; }

    // -----------------------------------------------------------------------
    // PreProc: preprocess k input blocks x_1,...,x_k
    // Returns (PI, SI) and also outputs the secret s (needed by caller for SI_0)
    // -----------------------------------------------------------------------
    std::pair<PPEPublicInput, PPEPrivateInput>
    preproc(const std::vector<Monomial>&              Q,       // monomial pattern
            const std::vector<std::vector<FieldElem>>& X,      // k blocks of n elements
            std::mt19937_64& rng) const
    {
        size_t k = params_.k;
        size_t n = params_.n;
        size_t m = Q.size();

        if (X.size() != k)
            throw std::invalid_argument("PPE::preproc: X must have k blocks");
        for (auto& xj : X)
            if (xj.size() != n)
                throw std::invalid_argument("PPE::preproc: each block must have n elements");

        // 1. Sample LPN secret s ← Zp^k
        auto s = F_.sample_vec(k, rng);

        // 2. Build PI: for each block j, LPN-encrypt x_j
        //    A_j ← Zp^{k×n}, e_j ~ Ber(k^{-δ})^n, b_j = s·A_j + e_j + x_j
        LPN lpn(LPNParams{ k, n, params_.delta, params_.p });
        PPEPublicInput PI;
        PI.flag = true;

        // Also track the noise vectors e_j for SI construction
        std::vector<std::vector<FieldElem>> E(k);  // e_j for each block

        for (size_t j = 0; j < k; j++) {
            // A_j is sampled deterministically from (params seed, j) so that
            // KeyGen and Enc always produce the same A matrix for the same block.
            // This is the prototype stand-in for PRF(MSK, j) in the full scheme.
            // The noise e_j is still fresh (from streaming rng) — it provides
            // the LPN security, while A being fixed is part of the public key.
            std::mt19937_64 det_rng(params_.p ^ (j * 0x9e3779b97f4a7c15ULL));
            auto A = F_.sample_matrix(k, n, det_rng);
            auto e = F_.sample_noise(n, params_.noise_rate(), rng);
            E[j] = e;

            std::vector<FieldElem> b(n);
            for (size_t col = 0; col < n; col++) {
                FieldElem sAcol = 0;
                for (size_t row = 0; row < k; row++)
                    sAcol = F_.add(sAcol, F_.mul(s[row], A[row][col]));
                b[col] = F_.add(F_.add(sAcol, e[col]), X[j][col]);
            }
            PI.blocks.push_back(LPNEncryption{ std::move(A), std::move(b) });
        }

        // 3. Build SI_0 = (1,s)^{⊗⌈d/2⌉}
        //    This is the tensor product of (1,s) with itself ⌈d/2⌉ times,
        //    giving all monomials of degree ≤ ⌈d/2⌉ in s.
        size_t half_d = (params_.d + 1) / 2;  // ⌈d/2⌉
        PPEPrivateInput SI;
        SI.SI0 = compute_SI0(s, half_d);

        // 4. Build SI_r for each monomial Q_r ∈ Q
        size_t t1 = params_.t1();
        size_t T  = params_.T();
        size_t t  = params_.t();

        SI.SI_monomials.resize(m);

        for (size_t r = 0; r < m; r++) {
            const Monomial& Qr = Q[r];
            SI_r& sir = SI.SI_monomials[r];
            sir.t1 = t1; sir.T = T; sir.t = t;

            // a. Compute Corr_{r,j} = Mon_{Qr}(x_j) - Mon_{Qr}(x_j + e_j)
            std::vector<FieldElem> Corr(k);
            bool flag_r = true;
            for (size_t j = 0; j < k; j++) {
                FieldElem xj_mon  = eval_monomial(F_, Qr, X[j]);
                // x_j + e_j (element-wise)
                auto xje = F_.vecadd(X[j], E[j]);
                FieldElem xje_mon = eval_monomial(F_, Qr, xje);
                Corr[j] = F_.sub(xj_mon, xje_mon);
            }

            // b. Assign Corr[j] into t1 matrices {M_{r,γ}} of size T×T
            //    M_{r,j1}[j2][j3] = Corr[j]  where (j1,j2,j3) = φ(j)
            //    We store M as vector<vector<FieldElem>> of shape t1 × (T*T)
            std::vector<std::vector<FieldElem>> M(t1,
                std::vector<FieldElem>(T * T, 0));

            for (size_t j = 0; j < k; j++) {
                auto [j1, j2, j3] = phi_map(j, t1, T);
                if (j1 >= t1 || j2 >= T || j3 >= T) {
                    PI.flag = false; flag_r = false; continue;
                }
                size_t existing = 0;
                for (size_t col = 0; col < T; col++)
                    if (M[j1][j2 * T + col] != 0 || M[j1][col * T + j3] != 0) {
                        // count nonzeros in this matrix to check overflow
                    }
                M[j1][j2 * T + j3] = F_.add(M[j1][j2 * T + j3], Corr[j]);
            }

            // Check overflow: count nonzeros per matrix ≤ t
            for (size_t gamma = 0; gamma < t1; gamma++) {
                size_t nz = 0;
                for (auto v : M[gamma]) if (v != 0) nz++;
                if (nz > t) { PI.flag = false; flag_r = false; }
            }

            // c. Factorize each M_{r,γ} = U_{r,γ} · V_{r,γ}^T (rank ≤ t)
            //    Encoding: for each nonzero M[γ][j2*T+j3], pick a column l,
            //    set U[γ][j2*t+l] = M[γ][j2*T+j3], V[γ][j3*t+l] = 1
            //    We flatten U as T*t vector (row-major: U[j2*t+l])
            sir.U.assign(t1, std::vector<FieldElem>(T * t, 0));
            sir.V.assign(t1, std::vector<FieldElem>(T * t, 0));

            for (size_t gamma = 0; gamma < t1; gamma++) {
                size_t col_idx = 0;  // next available column in U/V
                for (size_t j2 = 0; j2 < T; j2++) {
                    for (size_t j3 = 0; j3 < T; j3++) {
                        FieldElem val = M[gamma][j2 * T + j3];
                        if (val != 0 && col_idx < t) {
                            sir.U[gamma][j2 * t + col_idx] = val;
                            sir.V[gamma][j3 * t + col_idx] = 1;
                            col_idx++;
                        }
                    }
                }
            }
        }

        return { std::move(PI), std::move(SI) };
    }

    // -----------------------------------------------------------------------
    // Eval: compute f(x) from (PI, SI) using degree-(O(d),2) polynomials
    // Returns f(x_1,...,x_k) = Σ_{r,j} μ_{r,j} · Mon_{Q_r}(x_j)
    // -----------------------------------------------------------------------
    FieldElem eval(const PPEFunction& f,
                   const PPEPublicInput& PI,
                   const PPEPrivateInput& SI) const
    {
        if (!PI.flag) return F_.zero();

        size_t k  = params_.k;
        size_t t1 = params_.t1();
        size_t T  = params_.T();
        size_t t  = params_.t();
        size_t m  = f.Q.size();

        if (f.mu.size() != m)
            throw std::invalid_argument("PPE::eval: mu size mismatch");
        if (PI.blocks.size() != k)
            throw std::invalid_argument("PPE::eval: PI block count mismatch");
        if (SI.SI_monomials.size() != m)
            throw std::invalid_argument("PPE::eval: SI monomial count mismatch");

        // Recover s from SI_0 (degree-1 part is just (1, s[0], s[1], ..., s[k-1]))
        // SI_0 = tensor product (1,s)^{⊗⌈d/2⌉}, so s is at positions 1..k
        std::vector<FieldElem> s = extract_s_from_SI0(SI.SI0, k);

        FieldElem result = F_.zero();

        for (size_t r = 0; r < m; r++) {
            const Monomial& Qr  = f.Q[r];
            const SI_r& sir     = SI.SI_monomials[r];

            for (size_t j = 0; j < k; j++) {
                if (f.mu[r][j] == 0) continue;

                FieldElem w = compute_w(PI, SI, s, r, j, Qr, sir, t1, T, t);
                result = F_.add(result, F_.mul(f.mu[r][j], w));
            }
        }
        return result;
    }

    // -----------------------------------------------------------------------
    // Helper: verify that eval gives f(x) for given plaintext X
    // (used in tests to check correctness)
    // -----------------------------------------------------------------------
    FieldElem eval_plaintext(const PPEFunction& f,
                             const std::vector<std::vector<FieldElem>>& X) const
    {
        size_t k = params_.k;
        size_t m = f.Q.size();
        FieldElem result = F_.zero();
        for (size_t r = 0; r < m; r++)
            for (size_t j = 0; j < k; j++)
                result = F_.add(result,
                    F_.mul(f.mu[r][j], eval_monomial(F_, f.Q[r], X[j])));
        return result;
    }

private:
    PPEParams params_;
    Field     F_;

    // -----------------------------------------------------------------------
    // Compute SI_0 = (1,s)^{⊗h} for h = ⌈d/2⌉
    // This is the Kronecker product of (1, s[0], ..., s[k-1]) with itself h times.
    // Result is a vector of length (k+1)^h containing all degree-h monomials of s.
    // -----------------------------------------------------------------------
    std::vector<FieldElem> compute_SI0(const std::vector<FieldElem>& s,
                                       size_t h) const
    {
        size_t k = s.size();
        // base = (1, s[0], s[1], ..., s[k-1])
        std::vector<FieldElem> base(k + 1);
        base[0] = 1;
        for (size_t i = 0; i < k; i++) base[i+1] = s[i];

        // Tensor product: result = base ⊗ base ⊗ ... ⊗ base  (h times)
        std::vector<FieldElem> result = {F_.one()};
        for (size_t iter = 0; iter < h; iter++) {
            std::vector<FieldElem> next;
            next.reserve(result.size() * base.size());
            for (auto a : result)
                for (auto b : base)
                    next.push_back(F_.mul(a, b));
            result = std::move(next);
        }
        return result;
    }

    // -----------------------------------------------------------------------
    // Extract s from SI_0: SI_0[0] = 1, SI_0[1..k] = s[0..k-1]
    // (for h=1: SI_0 = (1, s[0],...,s[k-1]), length k+1)
    // (for h>1: first (k+1) entries are still (1^h, 1^{h-1}*s[0], ...))
    // We need s directly; for d ≤ 2 (h=1), it's at positions 1..k.
    // For d > 2 (h > 1), we store s separately alongside SI_0 in practice.
    // Here we encode s as positions 1..k of SI_0 when h=1, else separately.
    // -----------------------------------------------------------------------
    std::vector<FieldElem> extract_s_from_SI0(
        const std::vector<FieldElem>& SI0, size_t k) const
    {
        size_t h = (params_.d + 1) / 2;
        if (h == 1) {
            // SI0 = (1, s[0],...,s[k-1])
            if (SI0.size() < k + 1)
                throw std::runtime_error("extract_s_from_SI0: SI0 too short");
            return std::vector<FieldElem>(SI0.begin() + 1, SI0.begin() + 1 + k);
        } else {
            // For h > 1, s is still in positions 1..k of the first level
            // The tensor product layout: result[0]=1, result[1..k]=s[0..k-1]
            // (since first factor is (1,s) tensored on the right)
            // This holds because in the tensor product a ⊗ b, the first
            // (k+1) entries of a two-level tensor are b * a[0] = b * 1 = b.
            if (SI0.size() < k + 1)
                throw std::runtime_error("extract_s_from_SI0: SI0 too short for h>1");
            return std::vector<FieldElem>(SI0.begin() + 1, SI0.begin() + 1 + k);
        }
    }

    // -----------------------------------------------------------------------
    // Compute the polynomial w_{r,j}(PI, SI) from §4.1 Eval:
    //   w_{r,j} = flag · (Mon_{Qr}(b_j - s·A_j) + U_{r,j1}·V_{r,j1}[j2,j3])
    // where (j1,j2,j3) = φ(j).
    // -----------------------------------------------------------------------
    FieldElem compute_w(const PPEPublicInput& PI,
                        const PPEPrivateInput& SI,
                        const std::vector<FieldElem>& s,
                        size_t r, size_t j,
                        const Monomial& Qr,
                        const SI_r& sir,
                        size_t t1, size_t T, size_t t) const
    {
        if (!PI.flag) return F_.zero();

        const LPNEncryption& ct = PI.blocks[j];
        size_t k = params_.k;
        size_t n = params_.n;

        // Compute b_j - s·A_j  (= x_j + e_j, the noisy plaintext)
        std::vector<FieldElem> xe(n);
        for (size_t col = 0; col < n; col++) {
            FieldElem sAcol = 0;
            for (size_t row = 0; row < k; row++)
                sAcol = F_.add(sAcol, F_.mul(s[row], ct.A[row][col]));
            xe[col] = F_.sub(ct.b[col], sAcol);
        }

        // First term: Mon_{Qr}(x_j + e_j)
        FieldElem term1 = eval_monomial(F_, Qr, xe);

        // Second term: U_{r,j1}·V_{r,j1}^T[j2,j3] = Σ_l U[j1][j2,l]*V[j1][j3,l]
        auto [j1, j2, j3] = phi_map(j, t1, T);
        if (j1 >= t1 || j2 >= T || j3 >= T) return F_.zero();

        FieldElem term2 = F_.zero();
        const auto& Ug = sir.U[j1];  // T*t vector
        const auto& Vg = sir.V[j1];
        for (size_t l = 0; l < t; l++) {
            FieldElem u_val = (j2 * t + l < Ug.size()) ? Ug[j2 * t + l] : 0;
            FieldElem v_val = (j3 * t + l < Vg.size()) ? Vg[j3 * t + l] : 0;
            term2 = F_.add(term2, F_.mul(u_val, v_val));
        }

        return F_.add(term1, term2);
    }
};
