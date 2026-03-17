#pragma once
/**
 * prg_nc0.h — Pseudorandom Generator in NC⁰ (JLS22 §3, Definition 3.3)
 *
 * A PRG in NC⁰ is a Boolean function PRG: {0,1}^n → {0,1}^m where:
 *   - m = n^{1+τ} for some τ > 0  (polynomial stretch)
 *   - Each output bit depends on only d = O(1) input bits (NC⁰/locality d)
 *
 * We use Goldreich's construction [Gol00]:
 *   - Fix a random bipartite graph G: m output nodes, each connected to d=3 input nodes
 *   - Each output bit y_i = P(x_{G(i,0)}, x_{G(i,1)}, x_{G(i,2)})
 *   - Predicate P(a,b,c) = a·b ⊕ c  (degree-3, resists linear/affine attacks)
 *
 * In JLS22 this PRG is used in two places:
 *   1. ARE (§5): stretch r_i into (σ, b) for Yao garbling randomness
 *   2. ARE (§5): H stretches λ bits to 2λ+2 for garbled table computation
 *
 * Security note: Goldreich's PRG with P(a,b,c) = a·b⊕c is believed secure
 * but not proven under standard assumptions. For the prototype this suffices;
 * for production one would use a PRG based on a stream cipher (e.g., Trivium)
 * with the NC⁰ structure enforced by the circuit topology.
 */

#include <cstdint>
#include <vector>
#include <stdexcept>
#include <random>
#include <cassert>

// ---------------------------------------------------------------------------
// NC⁰ PRG parameters
// ---------------------------------------------------------------------------
struct PRGParams {
    size_t n;       // seed length (input bits)
    size_t m;       // output length (output bits), m = n^{1+tau} approximately
    size_t d;       // locality (each output bit touches d input bits), default 3

    // Compute stretch factor: m/n
    double stretch() const { return static_cast<double>(m) / n; }

    // Default small params for testing
    static PRGParams toy()   { return { 8,  32, 3 }; }
    static PRGParams small() { return { 16, 96, 3 }; }
    static PRGParams med()   { return { 32, 320, 3 }; }
};

// ---------------------------------------------------------------------------
// The predicate P: {0,1}^d → {0,1}
// Goldreich's degree-3 predicate: P(a,b,c) = a·b ⊕ c
// This has algebraic degree 2 in {a,b,c} over Z₂ (nonlinear, resists
// linear attacks). Higher-degree predicates can be swapped in here.
// ---------------------------------------------------------------------------
inline uint8_t predicate_deg3(uint8_t a, uint8_t b, uint8_t c) {
    return (a & b) ^ c;
}

// ---------------------------------------------------------------------------
// The PRG structure: sampled once at setup, then used for evaluation
// ---------------------------------------------------------------------------
class PRG_NC0 {
public:
    // Construct a PRG with given params.
    // The graph G (which input indices each output bit touches) is sampled
    // from rng at construction time and fixed for all evaluations.
    // This matches Goldreich's construction where G is a public parameter.
    PRG_NC0(const PRGParams& params, std::mt19937_64& rng)
        : params_(params)
    {
        if (params_.d < 1)
            throw std::invalid_argument("PRG_NC0: locality d must be >= 1");
        if (params_.m < params_.n)
            throw std::invalid_argument("PRG_NC0: output must be longer than input");

        // Sample the bipartite graph G:
        // For each output bit i, sample d distinct input indices from [0, n)
        std::uniform_int_distribution<size_t> idx_dist(0, params_.n - 1);
        graph_.resize(params_.m, std::vector<size_t>(params_.d));
        for (size_t i = 0; i < params_.m; i++) {
            // Sample d indices (with or without replacement — for small n
            // with replacement is fine; in practice use without replacement)
            for (size_t j = 0; j < params_.d; j++)
                graph_[i][j] = idx_dist(rng);
        }
    }

    const PRGParams& params() const { return params_; }

    // -----------------------------------------------------------------------
    // Evaluate PRG on a seed of n bits.
    // Output is m bits, each computed by the predicate on d input bits.
    // -----------------------------------------------------------------------
    std::vector<uint8_t> eval(const std::vector<uint8_t>& seed) const {
        if (seed.size() != params_.n)
            throw std::invalid_argument("PRG_NC0::eval: seed length mismatch");
        // Clamp seed bits to {0,1}
        std::vector<uint8_t> out(params_.m);
        for (size_t i = 0; i < params_.m; i++) {
            // Collect the d input bits
            const auto& nbrs = graph_[i];
            if (params_.d == 3) {
                out[i] = predicate_deg3(
                    seed[nbrs[0]] & 1,
                    seed[nbrs[1]] & 1,
                    seed[nbrs[2]] & 1
                );
            } else {
                // Generic: XOR of all inputs AND'd with first (degree-d fallback)
                // For d=1: output = x_{i0}
                // For d=2: output = x_{i0} XOR x_{i1}
                // For d>=3: first two AND'd, XOR remaining
                uint8_t acc = seed[nbrs[0]] & 1;
                for (size_t j = 1; j < params_.d - 1; j++)
                    acc &= (seed[nbrs[j]] & 1);
                acc ^= (seed[nbrs[params_.d - 1]] & 1);
                out[i] = acc;
            }
        }
        return out;
    }

    // -----------------------------------------------------------------------
    // Evaluate on a uint64 seed (convenience: unpacks bits LSB-first)
    // Only valid if n <= 64.
    // -----------------------------------------------------------------------
    std::vector<uint8_t> eval_u64(uint64_t seed_val) const {
        if (params_.n > 64)
            throw std::invalid_argument("PRG_NC0::eval_u64: n > 64");
        std::vector<uint8_t> seed(params_.n);
        for (size_t i = 0; i < params_.n; i++)
            seed[i] = static_cast<uint8_t>((seed_val >> i) & 1);
        return eval(seed);
    }

    // -----------------------------------------------------------------------
    // Pack output bits into uint64 words (for use as randomness in ARE)
    // -----------------------------------------------------------------------
    static std::vector<uint64_t> pack_bits(const std::vector<uint8_t>& bits) {
        size_t words = (bits.size() + 63) / 64;
        std::vector<uint64_t> out(words, 0);
        for (size_t i = 0; i < bits.size(); i++)
            out[i / 64] |= (static_cast<uint64_t>(bits[i] & 1) << (i % 64));
        return out;
    }

    // -----------------------------------------------------------------------
    // Accessor for the graph (needed by PPE to compute monomials symbolically)
    // graph_[i][j] = j-th input neighbour of output node i
    // -----------------------------------------------------------------------
    const std::vector<std::vector<size_t>>& graph() const { return graph_; }

    // -----------------------------------------------------------------------
    // For ARE: the PRG G that stretches (n_ARE - n_input) bits to
    // (n_ARE + m_ARE) * poly(lambda) bits of garbling randomness.
    // Returns a PRG with appropriate stretch for a circuit of size m_gates.
    // -----------------------------------------------------------------------
    static PRG_NC0 for_garbling(size_t seed_bits, size_t output_bits,
                                std::mt19937_64& rng)
    {
        PRGParams p;
        p.n = seed_bits;
        p.m = output_bits;
        p.d = 3;
        return PRG_NC0(p, rng);
    }

private:
    PRGParams params_;
    // graph_[i] = list of d input indices for output bit i
    std::vector<std::vector<size_t>> graph_;
};

// ---------------------------------------------------------------------------
// The H PRG used in ARE garbling (JLS22 §5.1):
// H: {0,1}^lambda → {0,1}^{2*lambda+2}
// This small-stretch PRG is applied per-wire in Yao garbling.
// In production this would be AES-based; here we use the same NC⁰ structure.
// ---------------------------------------------------------------------------
class PRG_H {
public:
    explicit PRG_H(size_t lambda, std::mt19937_64& rng)
        : lambda_(lambda),
          prg_(PRGParams{ lambda, 2 * lambda + 2, 3 }, rng)
    {}

    // Evaluate H on a lambda-bit label, producing (2*lambda+2) bits
    std::vector<uint8_t> eval(const std::vector<uint8_t>& label) const {
        return prg_.eval(label);
    }

    // First half: H0(label) = first (lambda+1) bits of H(label)
    std::vector<uint8_t> H0(const std::vector<uint8_t>& label) const {
        auto out = eval(label);
        out.resize(lambda_ + 1);
        return out;
    }

    // Second half: H1(label) = last (lambda+1) bits of H(label)
    std::vector<uint8_t> H1(const std::vector<uint8_t>& label) const {
        auto out = eval(label);
        std::vector<uint8_t> half(out.begin() + lambda_ + 1, out.end());
        return half;
    }

    size_t lambda() const { return lambda_; }

private:
    size_t   lambda_;
    PRG_NC0  prg_;
};
