#pragma once
/**
 * symbolic_garble.h  —  Symbolic garbling + full PPE→PHFE bridge
 *
 * Closes the gap between the JLS22 prototype and the actual construction.
 *
 * The gap (what was wrong):
 *   build_phfe_function_from_poly() pre-evaluated the PI-linear term
 *   Σ_l μ_l·PI[l] at KeyGen time, using the KeyGen's own PI sample.
 *   At Dec time the ciphertext carries a *different* PI, so the evaluation
 *   was wrong for any x ≠ x_keygen.  Additionally, the PPE error-correction
 *   matrices SI_r = {U, V} were completely absent from the PHFE SI vector,
 *   so the polynomial evaluated to Σ_l μ_l·(x[l]+e[l]) instead of Σ_l μ_l·x[l].
 *
 * The fix (what this file does):
 *
 *   1. Extend SI_full to include all PPE private output:
 *        SI_full = [ SI_0  |  U-blocks  |  V-blocks ]
 *      where SI_0 = (1,s), and U/V are the error-correction rank-t factor matrices
 *      from PPE §4.1.  n_PHFE is set to hold all of this.
 *
 *   2. Store μ in pi_linear (not pre-evaluated) so eval() computes
 *      Σ_l μ_l·PI_ciphertext[l] at Dec time against the real PI.
 *
 *   3. Embed the U·V^T quadratic correction into coeff_at_PI:
 *      for each (l,γ,j2,j3,c): coeff[u_pos(l,γ,j2,c)][v_pos(l,γ,j3,c)] += μ_l
 *      so that SI[u]·SI[v] = U[γ][j2·t+c] · V[γ][j3·t+c], and summing over all
 *      (γ,j2,j3,c) gives the full PPE correction Σ_{j} Corr_j[l] = -(noise at l).
 *
 *   Result: f_i(PI_ct, SI_full_ct) = Σ_l μ_l · x[l]   exactly,
 *   with no dependence on x at KeyGen time for XOR/NOT circuits.
 */

#include "field.h"
#include "are.h"
#include "ppe.h"
#include "phfe.h"
#include "pre.h"
#include <vector>
#include <array>

// ---------------------------------------------------------------------------
// LinearPoly: degree-1 polynomial Σ_l c_l·a[l] + const  over Zp
// ---------------------------------------------------------------------------
struct LinearPoly {
    std::vector<FieldElem> coeffs;
    FieldElem              constant     = 0;
    FieldPrime             p            = 0;
    bool                   is_nonlinear = false;

    size_t n() const { return coeffs.size(); }

    FieldElem eval(const std::vector<FieldElem>& a) const {
        Field F(p);
        FieldElem r = constant;
        for (size_t l = 0; l < std::min(coeffs.size(), a.size()); l++)
            r = F.add(r, F.mul(coeffs[l], a[l]));
        return r;
    }

    LinearPoly add(const LinearPoly& o) const {
        Field F(p);
        LinearPoly r; r.p = p;
        r.is_nonlinear = is_nonlinear || o.is_nonlinear;
        r.coeffs.resize(std::max(coeffs.size(), o.coeffs.size()), 0);
        r.constant = F.add(constant, o.constant);
        for (size_t l = 0; l < coeffs.size();   l++) r.coeffs[l] = F.add(r.coeffs[l], coeffs[l]);
        for (size_t l = 0; l < o.coeffs.size(); l++) r.coeffs[l] = F.add(r.coeffs[l], o.coeffs[l]);
        return r;
    }
    LinearPoly xor_with(const LinearPoly& o) const { return add(o); }
    LinearPoly scale(FieldElem scalar) const {
        Field F(p);
        LinearPoly r; r.p = p; r.is_nonlinear = is_nonlinear;
        r.constant = F.mul(constant, scalar);
        r.coeffs.resize(coeffs.size());
        for (size_t l = 0; l < coeffs.size(); l++) r.coeffs[l] = F.mul(coeffs[l], scalar);
        return r;
    }
    static LinearPoly zero(size_t n, FieldPrime p) {
        LinearPoly z; z.p = p; z.constant = 0; z.coeffs.resize(n, 0); return z;
    }
    static LinearPoly unit(size_t l, size_t n, FieldPrime p) {
        auto z = zero(n, p); z.coeffs[l] = 1; return z;
    }
    static LinearPoly constant_poly(FieldElem c, size_t n, FieldPrime p) {
        auto z = zero(n, p); z.constant = c; return z;
    }
};

// ---------------------------------------------------------------------------
// symbolic_garble: extract μ_{i,l} polynomials from the garbled circuit.
// y_i = Σ_l μ_{i,l} · a_κ[l]  — x-independent for XOR/NOT circuits.
// ---------------------------------------------------------------------------
struct GarblePolynomials { std::vector<LinearPoly> output_poly; };

inline GarblePolynomials symbolic_garble(
    const Circuit& C, size_t n_prime, size_t kappa,
    size_t m_ARE, FieldPrime p, std::mt19937_64& rng)
{
    int nw = C.n_wires();
    Field F(p);
    // Sample labels (structural, not used symbolically)
    std::vector<std::array<Label,2>> wire_labels(nw);
    std::vector<uint8_t> perm_bits(nw);
    for (int w = 0; w < nw; w++) {
        wire_labels[w][0] = rng() & ~(1ULL<<63);
        wire_labels[w][1] = rng() & ~(1ULL<<63);
        perm_bits[w]      = rng() & 1;
    }

    std::vector<LinearPoly> wp(nw, LinearPoly::zero(n_prime, p));
    for (int j = 0; j < C.n_inputs && j < (int)n_prime; j++)
        wp[j] = LinearPoly::unit(j, n_prime, p);

    for (auto& gate : C.gates) {
        switch (gate.type) {
            case GateType::CONST0: wp[gate.out] = LinearPoly::constant_poly(0,n_prime,p); break;
            case GateType::CONST1: wp[gate.out] = LinearPoly::constant_poly(1,n_prime,p); break;
            case GateType::NOT: {
                auto neg = wp[gate.in0].scale(F.neg(1));
                neg.constant = F.add(neg.constant, 1);
                wp[gate.out] = neg;
                break;
            }
            case GateType::XOR:
                wp[gate.out] = wp[gate.in0].xor_with(wp[gate.in1]); break;
            case GateType::AND: {
                LinearPoly nl = LinearPoly::zero(n_prime, p);
                nl.is_nonlinear = true;
                wp[gate.out] = nl;
                break;
            }
        }
    }

    size_t out_start = kappa * m_ARE;
    size_t out_end   = std::min(out_start + m_ARE, (size_t)C.n_outputs);
    GarblePolynomials result;
    result.output_poly.resize(out_end - out_start);
    for (size_t i = 0; i + out_start < out_end; i++)
        result.output_poly[i] = wp[C.output_wires[out_start + i]];
    return result;
}

// ---------------------------------------------------------------------------
// SI layout for full PPE private output packed into a flat vector:
//
//   SI_full = [ SI_0 (k+1 elems) | U-blocks | V-blocks ]
//
//   U-blocks: for each monomial l ∈ [n_Q], for each γ ∈ [t1]:
//               T·t entries  →  U[γ][j2·t+c]  for j2∈[T], c∈[t]
//             total: n_Q · t1 · T · t  entries
//
//   V-blocks: same shape, same indexing, starts after all U-blocks
//
//   u_pos(l, γ, j2, c) = (k+1) + (l·t1 + γ)·(T·t) + j2·t + c
//   v_pos(l, γ, j3, c) = (k+1) + n_Q·t1·T·t + (l·t1 + γ)·(T·t) + j3·t + c
// ---------------------------------------------------------------------------

inline size_t compute_n_phfe_full(const PPEParams& pp) {
    return (pp.k + 1) + 2 * pp.n * pp.t1() * pp.T() * pp.t();
}

inline std::vector<FieldElem> pack_si_full(
    const PPEPrivateInput& pre_SI,
    const PPEParams& pp,
    size_t n_phfe_full)
{
    std::vector<FieldElem> SI(n_phfe_full, 0);
    // SI_0
    for (size_t i = 0; i < std::min(pre_SI.SI0.size(), (size_t)(pp.k+1)); i++)
        SI[i] = pre_SI.SI0[i];

    size_t n_Q   = pp.n;
    size_t t1    = pp.t1(), T = pp.T(), t = pp.t();
    size_t U_base = pp.k + 1;
    size_t V_base = U_base + n_Q * t1 * T * t;

    for (size_t l = 0; l < std::min(n_Q, pre_SI.SI_monomials.size()); l++) {
        const auto& sir = pre_SI.SI_monomials[l];
        for (size_t g = 0; g < std::min(t1, sir.U.size()); g++) {
            for (size_t j2 = 0; j2 < T; j2++) {
                for (size_t c = 0; c < t; c++) {
                    size_t flat  = j2*t + c;
                    size_t u_idx = U_base + (l*t1+g)*(T*t) + flat;
                    size_t v_idx = V_base + (l*t1+g)*(T*t) + flat;
                    if (u_idx < n_phfe_full && flat < sir.U[g].size()) SI[u_idx] = sir.U[g][flat];
                    if (v_idx < n_phfe_full && flat < sir.V[g].size()) SI[v_idx] = sir.V[g][flat];
                }
            }
        }
    }
    return SI;
}

// ---------------------------------------------------------------------------
// build_phfe_function_from_poly
//
// Translates μ polynomial into a PHFE function f_i(PI, SI_full) that
// evaluates EXACTLY to Σ_l μ_l·x[l] at Dec time — no x dependence at KeyGen.
//
// Requires:
//   n_phfe == compute_n_phfe_full(pp)
//   SI_full packed with pack_si_full()
// ---------------------------------------------------------------------------
inline PHFEFunction build_phfe_function_from_poly(
    const LinearPoly&      poly,
    const PREPublicInput&  pre_PI,
    const PREPrivateInput& pre_SI,
    const PPEParams&       pp,
    size_t                 kappa,
    size_t                 n_phfe,
    const Field&           F)
{
    PHFEFunction f;
    f.n_PHFE = n_phfe;
    f.d      = 2;
    f.coeff_at_PI.assign(n_phfe * n_phfe, 0);
    f.pi_linear.assign(n_phfe, 0);

    if (poly.is_nonlinear || pre_PI.pi.blocks.empty()) return f;

    size_t k_block = std::min(kappa, pre_PI.pi.blocks.size() - 1);
    const auto& A  = pre_PI.pi.blocks[k_block].A;
    size_t n_prime = poly.n();
    size_t k_lpn   = A.size();
    size_t n_Q     = pp.n;
    size_t t1      = pp.t1(), T = pp.T(), t = pp.t();
    size_t U_base  = pp.k + 1;
    size_t V_base  = U_base + n_Q * t1 * T * t;

    // ── Part 1a: PI-linear term ───────────────────────────────────────────
    // PI carries all k_PRE blocks concatenated: PI[kappa*n' + l] = b_kappa[l].
    // For this chunk: pi_linear[kappa*n_prime + l] = μ_l
    // so eval() computes Σ_l μ_l · PI[kappa*n_prime + l] = Σ_l μ_l · b_kappa[l].
    // This is chunk-specific — different kappas write to different PI offsets.
    size_t pi_base = kappa * n_prime;  // offset into the flat PI vector
    for (size_t l = 0; l < std::min(n_prime, poly.coeffs.size()); l++) {
        if (poly.coeffs[l] == 0) continue;
        size_t pi_idx = pi_base + l;
        if (pi_idx < n_phfe) f.pi_linear[pi_idx] = poly.coeffs[l];
    }
    // Constant term  →  c_{0,0}·SI[0]·SI[0] = const·1·1
    f.coeff_at_PI[0] = poly.constant;

    // ── Part 1b: SI-linear correction from LPN ────────────────────────────
    // -Σ_l μ_l·A[j][l] stored as coeff[0][j+1]  (SI[0]·SI[j+1] = 1·s[j])
    for (size_t j = 0; j < k_lpn && j+1 < n_phfe; j++) {
        FieldElem c = F.zero();
        for (size_t l = 0; l < std::min(n_prime, poly.coeffs.size()); l++) {
            if (poly.coeffs[l] == 0) continue;
            FieldElem Ajl = (j < A.size() && l < A[j].size()) ? A[j][l] : 0;
            c = F.sub(c, F.mul(poly.coeffs[l], Ajl));
        }
        f.coeff_at_PI[j+1] = F.add(f.coeff_at_PI[j+1], c);
        // note: coeff[0][j+1] in row-major = index 0*n_phfe + (j+1) = j+1
    }

    // ── Part 2: PPE error-correction quadratic terms ──────────────────────
    // PPE §4.1: for each block j ∈ [k], (j1,j2,j3) = φ(j, t1, T).
    // The correction term for monomial l at block j is:
    //   Σ_c U[j1][j2*t+c] · V[j1][j3*t+c]   (dot product, rank-t)
    //
    // In PHFE: each column c contributes independently, so we emit:
    //   coeff[u_pos(l,j1,j2,c)][v_pos(l,j1,j3,c)] += μ_l
    // and SI[u] · SI[v] = U[j1][j2*t+c] · V[j1][j3*t+c]  in SI_full.
    //
    // KEY: we only emit the (j2,j3) pair = φ(j)_2, φ(j)_3 for each block j.
    //      Other (j2,j3) combinations are NOT valid PPE correction terms
    //      and must not be summed over — they would add spurious garbage.
    const auto& SI_monomials = pre_SI.si.SI_monomials;
    for (size_t l = 0; l < std::min(n_prime, poly.coeffs.size()); l++) {
        if (poly.coeffs[l] == 0 || l >= SI_monomials.size()) continue;
        const auto& sir = SI_monomials[l];
        // Only block j=kappa contributes to the correction for this chunk.
        // The ARE encoding for chunk kappa uses input block a_kappa = (x, r_kappa).
        // The PPE correction Corr[j] is nonzero only for block j=kappa.
        // Looping over all j ∈ [k] would add spurious corrections from other blocks.
        for (size_t j = kappa; j == kappa; j++) {  // single iteration: j = kappa
            auto [j1, j2, j3] = phi_map(j, t1, T);
            if (j1 >= t1 || j2 >= T || j3 >= T) continue;
            if (j1 >= sir.U.size()) continue;
            for (size_t c = 0; c < t; c++) {
                size_t u_pos = U_base + (l*t1 + j1)*(T*t) + j2*t + c;
                size_t v_pos = V_base + (l*t1 + j1)*(T*t) + j3*t + c;
                if (u_pos >= n_phfe || v_pos >= n_phfe) continue;
                f.coeff_at_PI[u_pos * n_phfe + v_pos] =
                    F.add(f.coeff_at_PI[u_pos * n_phfe + v_pos],
                          poly.coeffs[l]);
            }
        }
    }

    return f;
}
