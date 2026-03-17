#pragma once
/**
 * lpn.h — Learning Parity with Noise over Fp (JLS22 §3, Definition 3.2)
 *
 * The δ-LPN assumption over Fp states that for:
 *   A  ← Z_p^{k × n}     (public matrix, k = LPN dimension, n = sample count)
 *   s  ← Z_p^{1 × k}     (secret row vector)
 *   e  ← D_{r}^{1 × n}   (noise vector, each coord nonzero w.p. r = k^{-δ})
 *
 * the distribution (A, b = s·A + e) is computationally indistinguishable
 * from (A, u) for uniform u ← Z_p^{1×n}.
 *
 * In JLS22 the LPN instance is used to build the PPE preprocessing:
 *   PI_j = (A_j, b_j = s·A_j + e_j + x_j)   for each input block x_j
 *
 * This module provides:
 *   - LPNParams:   parameter selection (k, n, delta, prime p)
 *   - LPNInstance: one LPN sample (A, b = sA + e)
 *   - LPNEncrypt:  encryption of x as (A, b = sA + e + x)  [the PPE PI step]
 *   - LPNDecrypt:  recover erroneous x+e from (A, b, s) → b - sA = x + e
 */

#include "field.h"
#include <cstdint>
#include <vector>
#include <stdexcept>
#include <cmath>

// ---------------------------------------------------------------------------
// LPN parameters (JLS22 Definition 3.2)
// ---------------------------------------------------------------------------
struct LPNParams {
    size_t      k;       // LPN secret dimension
    size_t      n;       // number of samples (columns of A)
    double      delta;   // noise rate exponent: each coord nonzero w.p. k^{-delta}
    FieldPrime  p;       // field prime

    // Noise rate r = k^{-delta}
    double noise_rate() const {
        return std::pow(static_cast<double>(k), -delta);
    }

    // Suggested parameters for security level lambda (toy values for prototype)
    static LPNParams toy(FieldPrime p) {
        // k=32, n=64, delta=0.5 -> noise rate ~1/5.6
        return { 32, 64, 0.5, p };
    }

    static LPNParams small(FieldPrime p) {
        // k=64, n=128, delta=0.5
        return { 64, 128, 0.5, p };
    }
};

// ---------------------------------------------------------------------------
// A single LPN sample: public matrix A and noisy linear combination b = sA + e
// ---------------------------------------------------------------------------
struct LPNSample {
    std::vector<std::vector<FieldElem>> A;  // k × n matrix over Zp
    std::vector<FieldElem>              b;  // length-n vector = s·A + e
};

// ---------------------------------------------------------------------------
// LPN encryption of a plaintext vector x (the PPE public-input step)
// PI_j = (A_j, b_j = s·A_j + e_j + x_j)  for x_j ∈ Z_p^n
// ---------------------------------------------------------------------------
struct LPNEncryption {
    std::vector<std::vector<FieldElem>> A;  // k × n
    std::vector<FieldElem>              b;  // length-n  = s·A + e + x
    // NOTE: s and e are secret and not stored here.
    // The encryptor must keep s separately as the decryption key.
};

// ---------------------------------------------------------------------------
// LPN class: key generation, encryption, decryption
// ---------------------------------------------------------------------------
class LPN {
public:
    explicit LPN(const LPNParams& params)
        : params_(params), F_(params.p) {}

    const LPNParams& params() const { return params_; }
    const Field&     field()  const { return F_; }

    // -----------------------------------------------------------------------
    // Key generation: sample secret s ← Z_p^k
    // -----------------------------------------------------------------------
    std::vector<FieldElem> keygen(std::mt19937_64& rng) const {
        return F_.sample_vec(params_.k, rng);
    }

    // -----------------------------------------------------------------------
    // Generate a bare LPN sample: (A, b = s·A + e)
    // A is k×n; b is length n; s is the secret.
    // -----------------------------------------------------------------------
    LPNSample sample(const std::vector<FieldElem>& s,
                     std::mt19937_64& rng) const
    {
        if (s.size() != params_.k)
            throw std::invalid_argument("LPN::sample: secret dimension mismatch");

        // Sample A ← Z_p^{k × n}  (stored as n rows of k-dim vectors,
        // i.e. A[col][row] so that s·A can be computed as dot products)
        // Convention: A is k rows × n cols, b[j] = <s, A_col_j> + e_j
        auto A = F_.sample_matrix(params_.k, params_.n, rng);
        auto e = F_.sample_noise(params_.n, params_.noise_rate(), rng);

        // b[j] = sum_i s[i] * A[i][j] + e[j]
        std::vector<FieldElem> b(params_.n, 0);
        for (size_t j = 0; j < params_.n; j++) {
            FieldElem sAj = 0;
            for (size_t i = 0; i < params_.k; i++)
                sAj = F_.add(sAj, F_.mul(s[i], A[i][j]));
            b[j] = F_.add(sAj, e[j]);
        }

        return { std::move(A), std::move(b) };
    }

    // -----------------------------------------------------------------------
    // Encrypt plaintext x ∈ Z_p^n:  (A, b = s·A + e + x)
    // This is the PPE public-input construction from JLS22 §4.1:
    //   PI_j = (A_j, b_j = <s, A_j> + e_j + x_j)
    // -----------------------------------------------------------------------
    LPNEncryption encrypt(const std::vector<FieldElem>& s,
                          const std::vector<FieldElem>& x,
                          std::mt19937_64& rng) const
    {
        if (s.size() != params_.k)
            throw std::invalid_argument("LPN::encrypt: secret dimension mismatch");
        if (x.size() != params_.n)
            throw std::invalid_argument("LPN::encrypt: plaintext length must equal n");

        auto A = F_.sample_matrix(params_.k, params_.n, rng);
        auto e = F_.sample_noise(params_.n, params_.noise_rate(), rng);

        std::vector<FieldElem> b(params_.n);
        for (size_t j = 0; j < params_.n; j++) {
            FieldElem sAj = 0;
            for (size_t i = 0; i < params_.k; i++)
                sAj = F_.add(sAj, F_.mul(s[i], A[i][j]));
            // b[j] = <s, A_j> + e_j + x_j
            b[j] = F_.add(F_.add(sAj, e[j]), x[j]);
        }

        return { std::move(A), std::move(b) };
    }

    // -----------------------------------------------------------------------
    // Partial decryption: recover x + e (erroneous plaintext)
    // Given (A, b, s): b - s·A = x + e
    // Error correction to recover exact x is handled by the PPE SI layer.
    // -----------------------------------------------------------------------
    std::vector<FieldElem> decrypt_noisy(
        const LPNEncryption& ct,
        const std::vector<FieldElem>& s) const
    {
        if (s.size() != params_.k)
            throw std::invalid_argument("LPN::decrypt_noisy: secret dimension mismatch");

        std::vector<FieldElem> xe(params_.n);
        for (size_t j = 0; j < params_.n; j++) {
            FieldElem sAj = 0;
            for (size_t i = 0; i < params_.k; i++)
                sAj = F_.add(sAj, F_.mul(s[i], ct.A[i][j]));
            // b[j] - <s, A_j> = x[j] + e[j]
            xe[j] = F_.sub(ct.b[j], sAj);
        }
        return xe;
    }

    // -----------------------------------------------------------------------
    // Homomorphic evaluation of a degree-d monomial on the ciphertext
    // For monomial Q = {l_1, ..., l_d} (subset of [n]):
    //   Mon_Q(b - sA) = Mon_Q(x + e)  (erroneous monomial evaluation)
    // This is the core homomorphic operation used in PPE.
    // For degree 1: just returns (b - sA)[l_0]
    // For degree 2: returns (b - sA)[l_0] * (b - sA)[l_1]  mod p
    // -----------------------------------------------------------------------
    FieldElem eval_monomial_noisy(
        const LPNEncryption& ct,
        const std::vector<FieldElem>& s,
        const std::vector<size_t>& monomial_indices) const
    {
        auto xe = decrypt_noisy(ct, s);
        FieldElem result = F_.one();
        for (size_t idx : monomial_indices) {
            if (idx >= params_.n)
                throw std::out_of_range("LPN::eval_monomial: index out of range");
            result = F_.mul(result, xe[idx]);
        }
        return result;
    }

private:
    LPNParams params_;
    Field     F_;
};
