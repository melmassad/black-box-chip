#pragma once
/**
 * field.h — Prime field Zp arithmetic
 *
 * Implements arithmetic over Z_p for a 64-bit prime p.
 * Uses __int128 for intermediate products to avoid overflow.
 *
 * In the JLS22 paper, all arithmetic is over Z_p where p is a
 * poly(lambda)-bit prime. For our prototype we use a 64-bit prime,
 * which is sufficient to demonstrate correctness. For production,
 * swap FieldElem/FieldPrime to NTL::ZZ_p and compile with NTL.
 *
 * NTL swap instructions:
 *   - Replace FieldPrime with NTL::ZZ (arbitrary precision)
 *   - Replace FieldElem  with NTL::ZZ_p
 *   - Replace Field class with NTL::ZZ_pContext for modulus management
 *   - All arithmetic operators already match NTL::ZZ_p's interface
 */

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>
#include <random>
#include <cassert>
#include <iostream>

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------
using FieldPrime = uint64_t;
using FieldElem  = uint64_t;   // always in [0, p)

// ---------------------------------------------------------------------------
// Field class — holds the prime and provides all arithmetic
// ---------------------------------------------------------------------------
class Field {
public:
    // Construct field with given prime p.
    // p must be < 2^63 so that intermediate __int128 products don't overflow.
    explicit Field(FieldPrime p) : p_(p) {
        if (p < 2)
            throw std::invalid_argument("Field: prime must be >= 2");
    }

    FieldPrime prime() const { return p_; }

    // Reduce x into [0, p)
    FieldElem reduce(FieldPrime x) const {
        return static_cast<FieldElem>(x % p_);
    }

    // Modular addition
    FieldElem add(FieldElem a, FieldElem b) const {
        return static_cast<FieldElem>((static_cast<__int128>(a) + b) % p_);
    }

    // Modular subtraction
    FieldElem sub(FieldElem a, FieldElem b) const {
        return static_cast<FieldElem>((static_cast<__int128>(a) + p_ - b) % p_);
    }

    // Modular multiplication
    FieldElem mul(FieldElem a, FieldElem b) const {
        return static_cast<FieldElem>(
            (static_cast<__int128>(a) * b) % p_
        );
    }

    // Modular negation
    FieldElem neg(FieldElem a) const {
        return a == 0 ? 0 : static_cast<FieldElem>(p_ - a);
    }

    // Modular exponentiation via square-and-multiply
    FieldElem pow(FieldElem base, FieldPrime exp) const {
        FieldElem result = 1;
        base = reduce(base);
        while (exp > 0) {
            if (exp & 1) result = mul(result, base);
            base = mul(base, base);
            exp >>= 1;
        }
        return result;
    }

    // Modular inverse via Fermat's little theorem (p must be prime)
    FieldElem inv(FieldElem a) const {
        if (a == 0)
            throw std::domain_error("Field::inv: zero has no inverse");
        return pow(a, p_ - 2);
    }

    // Modular division
    FieldElem div(FieldElem a, FieldElem b) const {
        return mul(a, inv(b));
    }

    // Equality
    bool eq(FieldElem a, FieldElem b) const { return a == b; }

    FieldElem zero() const { return 0; }
    FieldElem one()  const { return 1; }

    // ---------------------------------------------------------------------------
    // Vector operations (for LPN, PPE)
    // ---------------------------------------------------------------------------

    // Inner product <a, b> mod p
    FieldElem dot(const std::vector<FieldElem>& a,
                  const std::vector<FieldElem>& b) const {
        if (a.size() != b.size())
            throw std::invalid_argument("Field::dot: size mismatch");
        __int128 acc = 0;
        for (size_t i = 0; i < a.size(); i++)
            acc += static_cast<__int128>(a[i]) * b[i];
        return static_cast<FieldElem>(acc % p_);
    }

    // Matrix-vector product M * v mod p
    // M is (rows x cols) stored row-major
    std::vector<FieldElem> matvec(
        const std::vector<std::vector<FieldElem>>& M,
        const std::vector<FieldElem>& v) const
    {
        size_t rows = M.size();
        if (rows == 0) return {};
        size_t cols = M[0].size();
        if (v.size() != cols)
            throw std::invalid_argument("Field::matvec: dimension mismatch");
        std::vector<FieldElem> out(rows, 0);
        for (size_t r = 0; r < rows; r++)
            out[r] = dot(M[r], v);
        return out;
    }

    // Element-wise vector addition
    std::vector<FieldElem> vecadd(
        const std::vector<FieldElem>& a,
        const std::vector<FieldElem>& b) const
    {
        if (a.size() != b.size())
            throw std::invalid_argument("Field::vecadd: size mismatch");
        std::vector<FieldElem> out(a.size());
        for (size_t i = 0; i < a.size(); i++)
            out[i] = add(a[i], b[i]);
        return out;
    }

    // Element-wise vector subtraction
    std::vector<FieldElem> vecsub(
        const std::vector<FieldElem>& a,
        const std::vector<FieldElem>& b) const
    {
        if (a.size() != b.size())
            throw std::invalid_argument("Field::vecsub: size mismatch");
        std::vector<FieldElem> out(a.size());
        for (size_t i = 0; i < a.size(); i++)
            out[i] = sub(a[i], b[i]);
        return out;
    }

    // Scalar-vector multiplication
    std::vector<FieldElem> scalvec(
        FieldElem scalar,
        const std::vector<FieldElem>& v) const
    {
        std::vector<FieldElem> out(v.size());
        for (size_t i = 0; i < v.size(); i++)
            out[i] = mul(scalar, v[i]);
        return out;
    }

    // Zero vector of length n
    std::vector<FieldElem> zerovec(size_t n) const {
        return std::vector<FieldElem>(n, 0);
    }

    // ---------------------------------------------------------------------------
    // Sampling
    // ---------------------------------------------------------------------------

    // Sample uniformly random element in [0, p)
    // Uses rejection sampling for small primes to avoid modular bias.
    // For large primes (> 2^32) uses direct reduction (bias negligible).
    FieldElem sample_uniform(std::mt19937_64& rng) const {
        if (p_ > (1ULL << 32)) {
            // For large primes, bias from simple reduction is negligible
            return static_cast<FieldElem>(rng() % p_);
        }
        // For small primes, rejection sampling to eliminate bias
        uint64_t limit = (UINT64_MAX / p_) * p_;
        while (true) {
            uint64_t r = rng();
            if (r < limit) return static_cast<FieldElem>(r % p_);
        }
    }

    // Sample a random vector of length n
    std::vector<FieldElem> sample_vec(size_t n, std::mt19937_64& rng) const {
        std::vector<FieldElem> v(n);
        for (auto& x : v) x = sample_uniform(rng);
        return v;
    }

    // Sample a random matrix of shape (rows x cols)
    std::vector<std::vector<FieldElem>> sample_matrix(
        size_t rows, size_t cols, std::mt19937_64& rng) const
    {
        std::vector<std::vector<FieldElem>> M(rows, std::vector<FieldElem>(cols));
        for (auto& row : M)
            for (auto& x : row)
                x = sample_uniform(rng);
        return M;
    }

    // Sample Bernoulli(rate) noise: nonzero with probability rate,
    // nonzero value is uniform in [1, p)
    FieldElem sample_bernoulli(double rate, std::mt19937_64& rng) const {
        std::bernoulli_distribution bern(rate);
        if (!bern(rng)) return 0;
        // nonzero: sample uniform in [1, p)
        FieldElem v = 0;
        while (v == 0) v = sample_uniform(rng);
        return v;
    }

    // Sample noise vector: each coord nonzero with probability rate
    std::vector<FieldElem> sample_noise(
        size_t n, double rate, std::mt19937_64& rng) const
    {
        std::vector<FieldElem> e(n);
        for (auto& x : e) x = sample_bernoulli(rate, rng);
        return e;
    }

private:
    FieldPrime p_;
};

// ---------------------------------------------------------------------------
// Useful prime constants for testing
// ---------------------------------------------------------------------------
namespace FieldPrimes {
    // A well-known 64-bit prime for testing: 2^61 - 1 (Mersenne prime)
    constexpr FieldPrime MERSENNE_61 = (1ULL << 61) - 1;
    // Smaller prime for fast unit tests
    constexpr FieldPrime SMALL_PRIME = 1000000007ULL;
    // Medium prime
    constexpr FieldPrime MED_PRIME   = 4611686018427387847ULL; // ~2^62
}
