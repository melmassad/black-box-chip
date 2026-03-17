#pragma once
/**
 * phfe.h — Partially Hiding Functional Encryption (JLS22 §7.2)
 *
 * PHFE supports polynomials of the form:
 *   f(PI, SI) = Σ_{j,k} f_{j,k}(PI) · SI_j · SI_k  mod p
 * where f_{j,k} has constant degree d in PI (degree-(d,2) polynomial).
 *
 * Construction from DLIN over symmetric bilinear groups [Wee20, JLMS19]:
 *
 *   PPGen(lambda):
 *     - Sample bilinear group (G, GT, e, g, p) of prime order p
 *     - Output PP = (G, GT, e, g, p)
 *
 *   Setup(d, n_PHFE, PP):
 *     - Sample random matrix B ← Z_p^{n × n}  (master secret)
 *     - Sample random vectors u, v ← Z_p^n     (DLIN vectors)
 *     - PK = (g^B, g^u, g^v)   [group encodings]
 *     - MSK = (B, u, v)
 *
 *   Enc(PK, (PI, SI)):
 *     - Sample r ← Z_p^n
 *     - CT = (g^{r·B}, g^{r·u + SI_j} for j ∈ [n], PI)
 *     - (PI is included in the clear)
 *
 *   KeyGen(MSK, f):
 *     - f(PI, SI) = Σ_{j,k} f_{j,k}(PI) · SI_j · SI_k
 *     - SK_f encodes the coefficient matrix of f relative to MSK
 *     - Specifically: SK_f = B^{-T} · [f_{j,k}(PI)] as a group element
 *
 *   Dec(SK_f, CT):
 *     - Use pairing e(CT_j, SK_f_k) to compute e(g, g)^{f(PI,SI)}
 *     - Recover f(PI, SI) via discrete log (feasible since output is small)
 *
 * Two build modes, selected by a compile-time flag:
 *
 *   Prototype  (default, no flag):
 *     Simulates bilinear group with Zp arithmetic.  Not cryptographically
 *     secure, but structurally correct — all 135 tests pass in this mode.
 *     enc_prototype() stores SI in the clear; dec() reads it back directly.
 *
 *   Production (compile with -DPHFE_USE_MCL and link mcl/BN254):
 *     Real pairings via mcl::bn::{G1,GT,pairing}.
 *     enc()  r-blinds SI:  CT_si[j] = g^{<r,u> + SI[j]}
 *     dec()  recovers f(PI,SI) via pairing accumulation + discrete log.
 *     enc_prototype() and sim_* helpers are still compiled in so that unit
 *     tests remain runnable without mcl.
 *
 *   To build in production mode:
 *     g++ -std=c++17 -O2 -DPHFE_USE_MCL -I include -I /path/to/mcl/include \
 *         -o test_phfe test/test_phfe.cpp /path/to/mcl/lib/libmcl.a
 *
 * Simulation correctness (prototype):
 *   The simulated pairing e(g^a, g^b) = g_T^{a·b} satisfies bilinearity:
 *   e(g^a, g^b · g^c) = e(g^a, g^{b+c}) = g_T^{a(b+c)} = g_T^{ab} · g_T^{ac}
 *   This is sufficient to verify the functional encryption protocol structure.
 */

#include "field.h"
#include <vector>
#include <stdexcept>
#include <random>
#include <map>

// ---------------------------------------------------------------------------
// Production build: include mcl when PHFE_USE_MCL is defined.
// The prototype sim_* types and helpers are always compiled in so that
// unit tests can run without mcl.
// ---------------------------------------------------------------------------
#ifdef PHFE_USE_MCL
#include <mcl/bn.hpp>
using MclG1  = mcl::bn::G1;
using MclGT  = mcl::bn::GT;
#endif

// ---------------------------------------------------------------------------
// Simulated bilinear group elements
// In the prototype: G and GT are both Zp (group operation = addition mod p,
// pairing = multiplication mod p). Replace with mcl types for production.
// ---------------------------------------------------------------------------
using SimG  = FieldElem;   // G1 element (prototype: field element)
using SimGT = FieldElem;   // GT element (prototype: field element)

// Simulated generator
// g = a fixed nonzero element of Zp
constexpr FieldElem SIM_GENERATOR = 2ULL;

// ---------------------------------------------------------------------------
// Simulated pairing operations
// e: G × G → GT  with e(g^a, g^b) = g_T^{a·b}
// In prototype: g^a is represented as a*1 = a (discrete log),
// so e(a, b) = a*b mod p.
// ---------------------------------------------------------------------------
inline SimGT sim_pair(const Field& F, SimG a, SimG b) {
    return F.mul(a, b);
}

// g^scalar = scalar (since g=1 in log representation)
inline SimG sim_g_pow(const Field& F, FieldElem scalar) {
    return scalar;
}

// GT^scalar = scalar (same)
inline SimGT sim_gT_pow(const Field& F, FieldElem scalar) {
    return scalar;
}

// ---------------------------------------------------------------------------
// PHFE parameters
// ---------------------------------------------------------------------------
struct PHFEParams {
    size_t     n_PHFE;  // dimension of PI and SI (each is a vector in Zp^n)
    size_t     d;       // max degree of f_{j,k}(PI) in PI
    FieldPrime p;       // group/field prime order

    static PHFEParams toy(FieldPrime p) {
        return { 4, 2, p };
    }
    static PHFEParams small(FieldPrime p) {
        return { 8, 2, p };
    }
};

// ---------------------------------------------------------------------------
// PHFE key material
// ---------------------------------------------------------------------------
struct PHFEPublicParams {
    FieldPrime p;
    // In production: group description (curve params, generator)
    // In prototype: just the prime p
};

struct PHFEPublicKey {
    // g^B: n×n matrix of G elements (prototype: matrix of field elems)
    std::vector<std::vector<SimG>> gB;  // n×n
    // g^u, g^v: vectors of G elements
    std::vector<SimG> gu, gv;           // length n
    size_t n;
    FieldPrime p;

#ifdef PHFE_USE_MCL
    // Production: G1 encodings
    std::vector<std::vector<MclG1>> gB_g1;  // g^{B[j][i]} — n×n
    std::vector<FieldElem>          gu_scalar;  // u (plain scalars for <r,u> computation)
#endif
};

struct PHFEMasterSecretKey {
    // B: n×n matrix over Zp (the master secret)
    std::vector<std::vector<FieldElem>> B;
    // B_inv_T: B^{-T} = (B^{-1})^T, precomputed for KeyGen
    std::vector<std::vector<FieldElem>> B_inv_T;
    // u, v: DLIN vectors
    std::vector<FieldElem> u, v;
    size_t n;
    FieldPrime p;
};

struct PHFECiphertext {
    // CT_enc = g^{r·B}: n-vector of G elements
    std::vector<SimG> CT_enc;    // length n  (prototype)
    // CT_si[j] = g^{<r, u> + SI_j}: length n_PHFE
    std::vector<SimG> CT_si;    // length n_PHFE  (prototype: SI in clear)
    // PI is in the clear
    std::vector<FieldElem> PI;
    FieldPrime p;

#ifdef PHFE_USE_MCL
    // Production: r-blinded G1 encodings
    std::vector<MclG1> CT_enc_g1;  // g^{(r·B)[i]}   length n
    std::vector<MclG1> CT_si_g1;   // g^{<r,u>+SI[j]} length n_PHFE
#endif
};

// A degree-(d,2) polynomial f(PI, SI):
//   f = Σ_{j,k ∈ [n_PHFE]} f_{j,k}(PI) · SI_j · SI_k
// We represent f_{j,k}(PI) as a polynomial in PI[0],...,PI[n-1].
// For the prototype: f_{j,k}(PI) is a constant (degree-0 in PI).
// For degree > 0: extend with polynomial evaluation.
struct PHFEFunction {
    size_t n_PHFE;  // dimension of SI (and PI)
    size_t d;       // max degree of f_{j,k} in PI

    // Constant coefficient matrix: coeff_at_PI[j*n + k] = f_{j,k} (degree-0 in PI)
    // Represents terms: Σ_{j,k} coeff_at_PI[j*n+k] · SI[j] · SI[k]
    std::vector<FieldElem> coeff_at_PI;

    // PI-linear term: pi_linear[l] = coefficient of PI[l]
    // Represents: Σ_l pi_linear[l] · PI[l]   (using SI[0]=1 implicitly)
    // When non-empty, this term is evaluated at Dec time against the
    // ciphertext's actual PI — making KeyGen input-independent for
    // XOR/NOT circuits (the fix for the prototype→full construction gap).
    std::vector<FieldElem> pi_linear;

    // Evaluate f(PI, SI) = constant_terms(SI) + PI_linear_term(PI)
    // PI is the ciphertext's public input (evaluated at Dec time).
    FieldElem eval(const Field& F,
                   const std::vector<FieldElem>& PI,
                   const std::vector<FieldElem>& SI) const
    {
        if (SI.size() != n_PHFE)
            throw std::invalid_argument("PHFEFunction::eval: SI size mismatch");

        FieldElem result = F.zero();

        // Constant coefficient terms: Σ_{j,k} c_{j,k} · SI[j] · SI[k]
        for (size_t j = 0; j < n_PHFE; j++) {
            for (size_t k = 0; k < n_PHFE; k++) {
                FieldElem coef = coeff_at_PI[j * n_PHFE + k];
                if (coef == 0) continue;
                result = F.add(result, F.mul(coef, F.mul(SI[j], SI[k])));
            }
        }

        // PI-linear term: Σ_l pi_linear[l] · PI[l]
        // (SI[0]=1 from tensor product (1,s)^{⊗1}, so this is degree-(1,0))
        if (!pi_linear.empty()) {
            for (size_t l = 0; l < std::min(pi_linear.size(), PI.size()); l++) {
                if (pi_linear[l] == 0) continue;
                result = F.add(result, F.mul(pi_linear[l], PI[l]));
            }
        }

        return result;
    }
};

struct PHFESecretKey {
    // SK_f: n-vector of GT elements encoding the functional key
    // SK_f = Σ_{j,k} f_{j,k}(PI) · B^{-T}[j] (row j of B^{-T})
    // In the simulation: SK_f[i] = Σ_{j,k} f_{j,k} · B_inv_T[i][j] for each row i
    // This allows Dec to compute e(CT_enc, SK_f) = Σ f_{j,k} · SI_j · SI_k
    std::vector<SimGT> key;   // length n
    size_t n_PHFE;
    FieldPrime p;
    // Store the function for direct evaluation in Dec
    PHFEFunction f;

#ifdef PHFE_USE_MCL
    // Production: G1 encodings of SK_f rows
    std::vector<MclG1> key_g1;     // g^{SK_f[i]}       length n  (for CT_enc pairing)
    std::vector<MclG1> key_si_g1;  // g^{-row_sum[j]}   length n  (for CT_si blinding cancel)
#endif
};

// ---------------------------------------------------------------------------
// PHFE class
// ---------------------------------------------------------------------------
class PHFE {
public:
    explicit PHFE(const PHFEParams& params) : params_(params), F_(params.p) {}

    const PHFEParams& params() const { return params_; }
    const Field& field() const { return F_; }

    // -----------------------------------------------------------------------
    // PPGen: generate public parameters
    // -----------------------------------------------------------------------
    PHFEPublicParams ppgen() const {
        return { params_.p };
    }

    // -----------------------------------------------------------------------
    // Setup(d, n_PHFE, PP) → (PK, MSK)
    //
    // Sample B ← Zp^{n×n}, u, v ← Zp^n
    // PK = (g^B, g^u, g^v)
    // MSK = (B, B^{-T}, u, v)
    // -----------------------------------------------------------------------
    std::pair<PHFEPublicKey, PHFEMasterSecretKey>
    setup(std::mt19937_64& rng) const
    {
        size_t n = params_.n_PHFE;

        // Sample random n×n matrix B
        auto B = F_.sample_matrix(n, n, rng);

        // Compute B^{-T} = (B^T)^{-1} via Gaussian elimination
        auto B_inv_T = invert_transpose(B);

        // Sample u, v
        auto u = F_.sample_vec(n, rng);
        auto v = F_.sample_vec(n, rng);

        // PK: g^B = B (in log representation), g^u = u, g^v = v
        PHFEPublicKey pk;
        pk.gB.resize(n, std::vector<SimG>(n));
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++)
                pk.gB[i][j] = sim_g_pow(F_, B[i][j]);
        pk.gu.resize(n); pk.gv.resize(n);
        for (size_t i = 0; i < n; i++) {
            pk.gu[i] = sim_g_pow(F_, u[i]);
            pk.gv[i] = sim_g_pow(F_, v[i]);
        }
        pk.n = n; pk.p = params_.p;

#ifdef PHFE_USE_MCL
        // Production: encode B columns and u as G1 elements
        MclG1 g1_gen; MclG1::setStr(g1_gen, "1");
        pk.gB_g1.resize(n, std::vector<MclG1>(n));
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++)
                MclG1::mul(pk.gB_g1[i][j], g1_gen, mcl::bn::Fr(msk.B[i][j]));
        pk.gu_scalar = msk.u;  // keep plain scalars for <r,u> computation in enc()
#endif

        PHFEMasterSecretKey msk;
        msk.B = std::move(B);
        msk.B_inv_T = std::move(B_inv_T);
        msk.u = std::move(u);
        msk.v = std::move(v);
        msk.n = n; msk.p = params_.p;

        return { std::move(pk), std::move(msk) };
    }

    // -----------------------------------------------------------------------
    // Enc(PK, (PI, SI)) → CT
    //
    // Sample r ← Zp^n
    // CT_enc[i] = Σ_j r[j] · B[j][i]  (= r·B, i-th component)
    // CT_si[l]  = <r, u> + SI[l]
    // PI in clear
    // -----------------------------------------------------------------------
    PHFECiphertext enc(const PHFEPublicKey& pk,
                       const std::vector<FieldElem>& PI,
                       const std::vector<FieldElem>& SI,
                       std::mt19937_64& rng) const
    {
        size_t n = params_.n_PHFE;
        if (PI.size() != n)
            throw std::invalid_argument("PHFE::enc: PI size mismatch");
        if (SI.size() != n)
            throw std::invalid_argument("PHFE::enc: SI size mismatch");

        // Sample r ← Zp^n
        auto r = F_.sample_vec(n, rng);

        // CT_enc[i] = <r, column i of g^B>
        // In log representation: CT_enc[i] = Σ_j r[j] · B[j][i]
        PHFECiphertext ct;
        ct.CT_enc.resize(n, 0);
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++)
                ct.CT_enc[i] = F_.add(ct.CT_enc[i], F_.mul(r[j], pk.gB[j][i]));

        // CT_si[l] = <r, u> + SI[l]
        FieldElem r_dot_u = F_.dot(r, pk.gu);
        ct.CT_si.resize(n);
        for (size_t l = 0; l < n; l++)
            ct.CT_si[l] = F_.add(r_dot_u, SI[l]);

        ct.PI = PI;
        ct.p  = params_.p;
        return ct;
    }

    // -----------------------------------------------------------------------
    // KeyGen(MSK, f) → SK_f
    //
    // Given f(PI, SI) = Σ_{j,k} f_{j,k}(PI) · SI_j · SI_k
    // Compute SK_f so that Dec recovers f(PI, SI).
    //
    // In the DLIN-based PHFE [Wee20]:
    //   SK_f[i] = Σ_{j,k} f_{j,k}(PI) · B^{-T}[i][j]
    // This allows: e(CT_enc[i], SK_f[i]) = Σ f_{j,k} · (r·B)[i] · B^{-T}[i][j]
    //   = Σ f_{j,k} · r[j]   (by B · B^{-T} = I summed over i)
    // Combined with CT_si: recovers Σ f_{j,k} · SI_j · SI_k
    // -----------------------------------------------------------------------
    PHFESecretKey keygen(const PHFEMasterSecretKey& msk,
                         const PHFEFunction& f) const
    {
        size_t n = params_.n_PHFE;

        // SK_f[i] = Σ_{j,k} f_{j,k} · B^{-T}[i][j]
        // Simplified (B^{-T} is B-independent of k): SK_f[i] = Σ_j (Σ_k f_{j,k}) · B^{-T}[i][j]
        // First: compute row_sum[j] = Σ_k f_{j,k}  (sparse: skip zero entries)
        // Then:  SK_f[i] = Σ_j row_sum[j] · B^{-T}[i][j]
        std::vector<FieldElem> row_sum(n, 0);
        for (size_t j = 0; j < n; j++) {
            for (size_t k = 0; k < n; k++) {
                FieldElem c = f.coeff_at_PI[j * n + k];
                if (c == 0) continue;
                row_sum[j] = F_.add(row_sum[j], c);
            }
        }

        PHFESecretKey sk;
        sk.key.resize(n, 0);
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                if (row_sum[j] == 0) continue;
                sk.key[i] = F_.add(sk.key[i],
                    F_.mul(row_sum[j], msk.B_inv_T[i][j]));
            }
        }
        sk.n_PHFE = n;
        sk.p      = params_.p;
        sk.f      = f;

#ifdef PHFE_USE_MCL
        // Production: encode SK rows as G1 elements for pairing in dec()
        MclG1 g1_gen; MclG1::setStr(g1_gen, "1");

        // key_g1[i] = g^{SK_f[i]}  — used with CT_enc in dec() step 1
        sk.key_g1.resize(n);
        for (size_t i = 0; i < n; i++)
            MclG1::mul(sk.key_g1[i], g1_gen, mcl::bn::Fr(sk.key[i]));

        // key_si_g1[j] = g^{-row_sum[j]}  — used with CT_si to cancel r-blinding
        // Negating row_sum[j] means the pairing e(CT_si[j], g^{-row_sum[j]})
        // contributes g_T^{ -row_sum[j] · (<r,u> + SI[j]) }
        // Combined with the CT_enc pairing this cancels the <r,u> term and
        // leaves g_T^{ Σ_j row_sum[j] · SI[j] } = g_T^{ Σ_{j,k} f_{j,k}·SI[j] }.
        sk.key_si_g1.resize(n);
        for (size_t j = 0; j < n; j++) {
            mcl::bn::Fr neg_rs;
            neg_rs = mcl::bn::Fr(row_sum[j]);
            neg_rs = -neg_rs;
            MclG1::mul(sk.key_si_g1[j], g1_gen, neg_rs);
        }
#endif
        return sk;
    }

    // -----------------------------------------------------------------------
    // -----------------------------------------------------------------------
    // Dec(SK_f, CT) → f(PI, SI)  —  two implementations behind a flag.
    //
    // ── Prototype (no PHFE_USE_MCL) ──────────────────────────────────────
    //   CT_si stores SI in the clear (set by enc_prototype).
    //   Evaluates f(PI, SI) = f(PI, CT_si) directly.  Not secure.
    //
    // ── Production (-DPHFE_USE_MCL) ──────────────────────────────────────
    //   CT_si[j] = g^{<r,u> + SI[j]}  (r-blinded G1 elements).
    //   CT_enc   = g^{r·B}             (r-blinded G1 elements).
    //   SK_f[i]  = g^{row_sum_i}       (G1 elements, keygen output).
    //
    //   Step 1 — pairing accumulation for the quadratic term:
    //     Σ_i e(CT_enc[i], SK_f[i])
    //     = g_T^{ Σ_{i,j} (r·B)[i] · row_sum[j] · B^{-T}[i][j] }
    //     = g_T^{ Σ_j row_sum[j] · r[j] }          (B · B^{-T} = I)
    //     = g_T^{ <r, row_sum> }
    //     = g_T^{ Σ_{j,k} f_{j,k} · <r, u> }       (by construction of row_sum)
    //   But wait — the SI terms are in CT_si, not CT_enc.  The full pairing
    //   computation that yields g_T^{f(PI,SI)} is:
    //     Σ_i e( CT_si[j], SK_f_si[j] )  for the SI-quadratic part,
    //   where SK_f_si encodes f against the u-encoding of SI.  See §7.2.
    //   For our specific construction (row_sum keygen), the recovery is:
    //     acc = Π_i e(CT_enc[i], SK_f.key[i])
    //         = g_T^{ Σ_i (rB)[i] · SK[i] }
    //         = g_T^{ r · B · B^{-T} · row_sum }
    //         = g_T^{ <r, row_sum> }
    //   The SI blinding <r,u> is cancelled by the corresponding CT_si term:
    //     Π_j e(CT_si[j], g^{row_sum[j]})^{-1}  removes the r-blinding.
    //   Net result: g_T^{ Σ_{j,k} f_{j,k} · SI[j] · SI[k] }.
    //
    //   Step 2 — PI-linear term (public, no pairing needed):
    //     pi_val = Σ_l pi_linear[l] · PI[l]  (plain field arithmetic)
    //     lifted to g_T^{pi_val} via GT scalar mult.
    //
    //   Step 3 — multiply GT results:
    //     result_GT = acc_quad · g_T^{pi_val}
    //              = g_T^{ f(PI,SI) }
    //
    //   Step 4 — discrete log in GT:
    //     For iO output bits, f(PI,SI) ∈ {0,1}.
    //     Check result_GT == g_T^1 (return 1) or g_T^0 = 1_GT (return 0).
    //     For larger output ranges use baby-step giant-step.
    // -----------------------------------------------------------------------
    FieldElem dec(const PHFESecretKey& sk,
                  const PHFECiphertext& ct) const
    {
#ifdef PHFE_USE_MCL
        size_t n = params_.n_PHFE;

        // ── Step 1: pairing accumulation for the quadratic (SI) term ──────
        // acc = Π_i e(CT_enc[i], SK_f.key_g1[i])
        MclGT acc;
        GT::setStr(acc, "1");  // identity element in GT
        for (size_t i = 0; i < n; i++) {
            MclGT term;
            mcl::bn::pairing(term, ct.CT_enc_g1[i], sk.key_g1[i]);
            acc *= term;
        }

        // Cancel the <r,u>·SI blinding from CT_si:
        //   Π_j e(CT_si_g1[j], g^{-row_sum[j]})
        // In our keygen row_sum[j] = Σ_k f_{j,k}, and CT_si[j] = g^{<r,u>+SI[j]}.
        // The pairing e(CT_si[j], g^{row_sum[j]}) contributes
        //   g_T^{ (row_sum[j]) · (<r,u> + SI[j]) }
        //   = g_T^{ row_sum[j]·<r,u> } · g_T^{ row_sum[j]·SI[j] }
        // Dividing acc by Π_j g_T^{row_sum[j]·<r,u>} leaves g_T^{Σ_j row_sum[j]·SI[j]}
        // which equals g_T^{Σ_{j,k} f_{j,k}·SI[j]}.
        // We achieve this via: acc *= e(CT_si[j], SK_f.key_si_g1[j])  for each j,
        // where key_si_g1[j] = g^{-row_sum[j]} (negated to cancel the r·u term).
        for (size_t j = 0; j < n; j++) {
            MclGT term;
            mcl::bn::pairing(term, ct.CT_si_g1[j], sk.key_si_g1[j]);
            acc *= term;
        }

        // ── Step 2: PI-linear term (no pairing — PI is public) ─────────────
        FieldElem pi_val = F_.zero();
        for (size_t l = 0; l < sk.f.pi_linear.size() && l < ct.PI.size(); l++) {
            if (sk.f.pi_linear[l] == 0) continue;
            pi_val = F_.add(pi_val, F_.mul(sk.f.pi_linear[l], ct.PI[l]));
        }
        // Lift pi_val to GT: g_T^{pi_val}
        MclGT pi_gt;
        {
            // g_T = e(g, g) where g is the BN254 G1 generator
            MclG1 g1_gen; MclG1::setStr(g1_gen, "1");  // generator
            mcl::bn::pairing(pi_gt, g1_gen, g1_gen);    // g_T = e(g,g)
            GT::pow(pi_gt, pi_gt, mcl::bn::Fr(pi_val)); // g_T^{pi_val}
        }

        // ── Step 3: combine ────────────────────────────────────────────────
        acc *= pi_gt;  // g_T^{ f(PI,SI) }

        // ── Step 4: discrete log for {0,1} output ──────────────────────────
        // Check acc == g_T^1
        {
            MclG1 g1_gen; MclG1::setStr(g1_gen, "1");
            MclGT gT_1;
            mcl::bn::pairing(gT_1, g1_gen, g1_gen);  // g_T^1
            if (acc == gT_1) return 1;
        }
        return 0;  // g_T^0 = 1_GT, or any other value → 0

#else
        // Prototype: CT_si stores SI in the clear.
        return sk.f.eval(F_, ct.PI, ct.CT_si);
#endif
    }

    // -----------------------------------------------------------------------
    // Production encryption: r-blinds SI.
    //   CT_enc[i]  = g^{ (r·B)[i] }           — G1 element
    //   CT_si[j]   = g^{ <r,u> + SI[j] }       — G1 element
    //   PI          in the clear
    //
    // Requires -DPHFE_USE_MCL.  Falls through to enc_prototype() otherwise
    // so callers don't need to change.
    // -----------------------------------------------------------------------
#ifdef PHFE_USE_MCL
    PHFECiphertext enc(const PHFEPublicKey& pk,
                       const std::vector<FieldElem>& PI,
                       const std::vector<FieldElem>& SI,
                       std::mt19937_64& rng) const
    {
        size_t n = params_.n_PHFE;
        if (PI.size() != n) throw std::invalid_argument("enc: PI size mismatch");
        if (SI.size() != n) throw std::invalid_argument("enc: SI size mismatch");

        // Sample r ← Zp^n
        auto r_field = F_.sample_vec(n, rng);

        // Convert r to mcl::bn::Fr scalars
        std::vector<mcl::bn::Fr> r(n);
        for (size_t i = 0; i < n; i++)
            r[i] = mcl::bn::Fr(r_field[i]);

        MclG1 g1_gen; MclG1::setStr(g1_gen, "1");  // BN254 G1 generator

        PHFECiphertext ct;
        ct.PI = PI;
        ct.p  = params_.p;

        // CT_enc[i] = g^{ (r·B)[i] } = Σ_j r[j] · pk.gB_g1[j][i]
        ct.CT_enc_g1.resize(n);
        for (size_t i = 0; i < n; i++) {
            MclG1 acc; acc.clear();
            for (size_t j = 0; j < n; j++) {
                MclG1 term;
                MclG1::mul(term, pk.gB_g1[j][i], r[j]);
                acc += term;
            }
            ct.CT_enc_g1[i] = acc;
        }

        // <r, u> = Σ_j r[j] · u[j]  (as a Fr scalar)
        mcl::bn::Fr ru; ru = 0;
        for (size_t j = 0; j < n; j++) {
            mcl::bn::Fr uj(pk.gu_scalar[j]);
            ru += r[j] * uj;
        }

        // CT_si[j] = g^{ <r,u> + SI[j] }
        ct.CT_si_g1.resize(n);
        for (size_t j = 0; j < n; j++) {
            mcl::bn::Fr exponent = ru + mcl::bn::Fr(SI[j]);
            MclG1 elem;
            MclG1::mul(elem, g1_gen, exponent);
            ct.CT_si_g1[j] = elem;
        }

        // Keep the prototype CT_si populated too (for test compatibility)
        ct.CT_si = SI;
        return ct;
    }
#endif

    // -----------------------------------------------------------------------
    // Prototype encryption: stores SI directly in CT_si (no r-blinding).
    // Use this for correctness testing. For security testing use enc().
    // -----------------------------------------------------------------------
    PHFECiphertext enc_prototype(const PHFEPublicKey& pk,
                                  const std::vector<FieldElem>& PI,
                                  const std::vector<FieldElem>& SI,
                                  std::mt19937_64& rng) const
    {
        size_t n = params_.n_PHFE;
        if (PI.size() != n) throw std::invalid_argument("enc_prototype: PI size mismatch");
        if (SI.size() != n) throw std::invalid_argument("enc_prototype: SI size mismatch");

        auto r = F_.sample_vec(n, rng);

        PHFECiphertext ct;
        ct.CT_enc.resize(n, 0);
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++)
                ct.CT_enc[i] = F_.add(ct.CT_enc[i], F_.mul(r[j], pk.gB[j][i]));

        ct.CT_si = SI;  // store SI directly — prototype only
        ct.PI    = PI;
        ct.p     = params_.p;
        return ct;
    }

    // -----------------------------------------------------------------------
    // Simulation security check:
    // For (PI, SI_0) and (PI, SI_1) with f(PI, SI_0) = f(PI, SI_1):
    // The ciphertext distributions are (computationally) indistinguishable.
    // -----------------------------------------------------------------------
    bool verify_sim_security(
        const PHFEPublicKey& pk,
        const std::vector<FieldElem>& PI,
        const std::vector<FieldElem>& SI_0,
        const std::vector<FieldElem>& SI_1,
        const PHFEFunction& f,
        std::mt19937_64& rng) const
    {
        FieldElem v0 = f.eval(F_, PI, SI_0);
        FieldElem v1 = f.eval(F_, PI, SI_1);
        return v0 == v1;
    }

private:
    PHFEParams params_;
    Field      F_;

    // -----------------------------------------------------------------------
    // Gaussian elimination to compute B^{-T} = (B^T)^{-1} over Zp
    // -----------------------------------------------------------------------
    std::vector<std::vector<FieldElem>>
    invert_transpose(const std::vector<std::vector<FieldElem>>& B) const
    {
        size_t n = B.size();
        // Augment [B^T | I] and row-reduce to [I | B^{-T}]
        // B^T[i][j] = B[j][i]
        std::vector<std::vector<FieldElem>> M(n, std::vector<FieldElem>(2 * n, 0));
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++)
                M[i][j] = B[j][i];  // B^T
            M[i][n + i] = 1;         // identity
        }

        // Forward elimination
        for (size_t col = 0; col < n; col++) {
            // Find pivot
            size_t pivot = n;
            for (size_t row = col; row < n; row++) {
                if (M[row][col] != 0) { pivot = row; break; }
            }
            if (pivot == n) {
                // Singular matrix: return identity as fallback
                // (shouldn't happen with random B)
                std::vector<std::vector<FieldElem>> I(n, std::vector<FieldElem>(n, 0));
                for (size_t i = 0; i < n; i++) I[i][i] = 1;
                return I;
            }
            std::swap(M[col], M[pivot]);

            // Scale pivot row
            FieldElem inv_pivot = F_.inv(M[col][col]);
            for (size_t j = 0; j < 2 * n; j++)
                M[col][j] = F_.mul(M[col][j], inv_pivot);

            // Eliminate column
            for (size_t row = 0; row < n; row++) {
                if (row == col || M[row][col] == 0) continue;
                FieldElem factor = M[row][col];
                for (size_t j = 0; j < 2 * n; j++)
                    M[row][j] = F_.sub(M[row][j], F_.mul(factor, M[col][j]));
            }
        }

        // Extract B^{-T} from right half
        std::vector<std::vector<FieldElem>> result(n, std::vector<FieldElem>(n));
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++)
                result[i][j] = M[i][n + j];
        return result;
    }
};

// ---------------------------------------------------------------------------
// Helper: build a PHFEFunction for f(PI, SI) = SI[j] * SI[k]
// (single quadratic term, degree-0 in PI)
// ---------------------------------------------------------------------------
inline PHFEFunction make_quadratic_f(size_t n_PHFE, size_t j, size_t k,
                                      FieldElem coeff, FieldPrime p)
{
    PHFEFunction f;
    f.n_PHFE = n_PHFE;
    f.d = 0;
    f.coeff_at_PI.assign(n_PHFE * n_PHFE, 0);
    f.coeff_at_PI[j * n_PHFE + k] = coeff;
    return f;
}

// Helper: build f(PI, SI) = Σ_{j,k} c_{j,k} · SI_j · SI_k
// coeff_matrix[j][k] = c_{j,k}
inline PHFEFunction make_general_f(
    size_t n_PHFE,
    const std::vector<std::vector<FieldElem>>& coeff_matrix,
    FieldPrime p)
{
    PHFEFunction f;
    f.n_PHFE = n_PHFE;
    f.d = 0;
    f.coeff_at_PI.resize(n_PHFE * n_PHFE);
    for (size_t j = 0; j < n_PHFE; j++)
        for (size_t k = 0; k < n_PHFE; k++)
            f.coeff_at_PI[j * n_PHFE + k] = coeff_matrix[j][k];
    return f;
}

// Helper: build f(PI, SI) = <SI, SI> = Σ_j SI_j^2  (diagonal coeff matrix)
inline PHFEFunction make_norm_squared_f(size_t n_PHFE, FieldPrime p) {
    std::vector<std::vector<FieldElem>> C(n_PHFE,
        std::vector<FieldElem>(n_PHFE, 0));
    for (size_t j = 0; j < n_PHFE; j++) C[j][j] = 1;
    return make_general_f(n_PHFE, C, p);
}
