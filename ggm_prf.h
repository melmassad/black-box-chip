#pragma once
/**
 * ggm_prf.h — GGM Pseudorandom Function Family built from PRG_NC0
 *
 * Construction [GGM84]:
 *   Given a PRG G: {0,1}^n → {0,1}^{2n} (length-doubling), define a PRF
 *   F_k: {0,1}^λ → {0,1}^n by walking a complete binary tree of depth λ:
 *
 *     F_k(x_1,...,x_λ):
 *       v_0 = k                          (root = key, length n)
 *       for i = 1 to λ:
 *         (v_L, v_R) = G(v_{i-1})        (split 2n-bit output in half)
 *         v_i = v_{x_i}                  (select left or right child)
 *       return v_λ
 *
 *   Security: if G is a secure PRG, then F is a secure PRF [GGM84].
 *
 * In JLS22 the GGM PRF is the core of the LCI locking scheme:
 *   - The PRF key k* is the Trivium keystream seed
 *   - F_{k*}(x) is XOR'd with circuit outputs to produce the locked circuit
 *   - Security reduces to PRG security of Trivium under stream-cipher assumptions
 *
 * This implementation uses PRG_NC0 as the underlying PRG.
 * The PRG is configured with m = 2n (length-doubling stretch).
 *
 * Usage:
 *   std::mt19937_64 rng(42);
 *   GGM_PRF prf(n_bits, depth, rng);       // build PRF with n-bit key, depth λ
 *   auto key   = prf.keygen(rng);           // sample k ← {0,1}^n
 *   auto value = prf.eval(key, input);      // F_k(x), input is λ bits
 *
 * Parameters:
 *   n     = node size in bits (= key size = output size)
 *   depth = input length λ (tree depth)
 */

#include "prg_nc0.h"
#include <vector>
#include <stdexcept>
#include <random>
#include <cassert>

// ---------------------------------------------------------------------------
// GGM PRF
// ---------------------------------------------------------------------------
class GGM_PRF {
public:
    // -----------------------------------------------------------------------
    // Construct a GGM PRF with:
    //   node_bits = n  (key size, output size)
    //   depth     = λ  (input length; domain is {0,1}^λ)
    //
    // The underlying PRG maps n bits → 2n bits.
    // The PRG graph is sampled once at construction and fixed.
    // -----------------------------------------------------------------------
    GGM_PRF(size_t node_bits, size_t depth, std::mt19937_64& rng)
        : node_bits_(node_bits),
          depth_(depth),
          prg_(PRGParams{ node_bits, 2 * node_bits, 3 }, rng)
    {
        if (node_bits < 2)
            throw std::invalid_argument("GGM_PRF: node_bits must be >= 2");
        if (depth < 1)
            throw std::invalid_argument("GGM_PRF: depth must be >= 1");
    }

    size_t node_bits() const { return node_bits_; }
    size_t depth()     const { return depth_; }

    // -----------------------------------------------------------------------
    // Key generation: sample k ← {0,1}^{node_bits} uniformly
    // -----------------------------------------------------------------------
    std::vector<uint8_t> keygen(std::mt19937_64& rng) const {
        std::vector<uint8_t> key(node_bits_);
        for (auto& b : key) b = rng() & 1;
        return key;
    }

    // -----------------------------------------------------------------------
    // PRF evaluation: F_k(x) where x ∈ {0,1}^depth
    //
    // Walks the GGM tree from root k, selecting left (0) or right (1)
    // child at each level according to x[i].
    // -----------------------------------------------------------------------
    std::vector<uint8_t> eval(const std::vector<uint8_t>& key,
                               const std::vector<uint8_t>& x) const
    {
        if (key.size() != node_bits_)
            throw std::invalid_argument("GGM_PRF::eval: key size mismatch");
        if (x.size() != depth_)
            throw std::invalid_argument("GGM_PRF::eval: input length must equal depth");

        std::vector<uint8_t> node = key;

        for (size_t i = 0; i < depth_; i++) {
            // Expand: G(node) → 2n bits
            auto expanded = prg_.eval(node);

            // Split into left (bits 0..n-1) and right (bits n..2n-1)
            std::vector<uint8_t> left (expanded.begin(),              expanded.begin() + node_bits_);
            std::vector<uint8_t> right(expanded.begin() + node_bits_, expanded.end());

            // Select child based on x[i]
            node = (x[i] & 1) ? right : left;
        }

        return node;
    }

    // -----------------------------------------------------------------------
    // Evaluate on an integer input (convenience for small depths ≤ 64)
    // -----------------------------------------------------------------------
    std::vector<uint8_t> eval_int(const std::vector<uint8_t>& key,
                                   uint64_t x_int) const
    {
        if (depth_ > 64)
            throw std::invalid_argument("GGM_PRF::eval_int: depth > 64");
        std::vector<uint8_t> x(depth_);
        for (size_t i = 0; i < depth_; i++)
            x[i] = (x_int >> i) & 1;
        return eval(key, x);
    }

    // -----------------------------------------------------------------------
    // Batch evaluate: compute F_k(0), F_k(1), ..., F_k(2^depth - 1)
    // by traversing the full tree (more efficient than 2^depth separate evals).
    //
    // Returns a vector of 2^depth node values, indexed by integer input.
    // Only feasible for small depth (≤ ~20).
    // -----------------------------------------------------------------------
    std::vector<std::vector<uint8_t>>
    eval_all(const std::vector<uint8_t>& key) const
    {
        if (depth_ > 20)
            throw std::invalid_argument("GGM_PRF::eval_all: depth > 20 (too large)");
        if (key.size() != node_bits_)
            throw std::invalid_argument("GGM_PRF::eval_all: key size mismatch");

        size_t n_leaves = 1ULL << depth_;

        // Build leaves indexed to match eval_int's LSB-first bit ordering:
        //   leaf index x_int: x[i] = (x_int >> i) & 1  (bit i = level i selector)
        // We compute this by iterating over all x_int values directly.
        // For tractability (depth ≤ 20), this is at most 1M evaluations.
        std::vector<std::vector<uint8_t>> leaves(n_leaves);
        for (size_t x_int = 0; x_int < n_leaves; x_int++) {
            std::vector<uint8_t> x(depth_);
            for (size_t i = 0; i < depth_; i++)
                x[i] = (x_int >> i) & 1;
            leaves[x_int] = eval(key, x);
        }
        return leaves;
    }

    // -----------------------------------------------------------------------
    // Punctured key: GGM_PRF with a hole at input x*
    //
    // A punctured key sk_{k,x*} allows evaluating F_k(x) for all x ≠ x*
    // but reveals nothing about F_k(x*). This is the key tool in the
    // AJ15/BV15 iO construction for making FE.KeyGen input-independent.
    //
    // The punctured key consists of the sibling nodes along the path to x*:
    //   For each level i, store the sibling of the node on the x* path.
    //   To evaluate F_k(x) for x ≠ x*: follow the x path until it diverges
    //   from the x* path, then use the stored sibling node.
    //
    // Returns: vector of depth sibling nodes (each node_bits bits)
    // -----------------------------------------------------------------------
    std::vector<std::vector<uint8_t>>
    puncture(const std::vector<uint8_t>& key,
             const std::vector<uint8_t>& x_star) const
    {
        if (key.size() != node_bits_)
            throw std::invalid_argument("GGM_PRF::puncture: key size mismatch");
        if (x_star.size() != depth_)
            throw std::invalid_argument("GGM_PRF::puncture: x_star length mismatch");

        // Walk the x* path, storing siblings at each level
        std::vector<std::vector<uint8_t>> siblings(depth_);
        std::vector<uint8_t> node = key;

        for (size_t i = 0; i < depth_; i++) {
            auto expanded = prg_.eval(node);
            std::vector<uint8_t> left (expanded.begin(),
                                        expanded.begin() + node_bits_);
            std::vector<uint8_t> right(expanded.begin() + node_bits_,
                                        expanded.end());
            // Sibling = the child NOT on the x* path
            siblings[i] = (x_star[i] & 1) ? left : right;
            // Continue down the x* path
            node = (x_star[i] & 1) ? right : left;
        }

        return siblings;
    }

    // -----------------------------------------------------------------------
    // Evaluate using a punctured key at input x (x must differ from x*)
    // -----------------------------------------------------------------------
    std::vector<uint8_t>
    eval_punctured(const std::vector<std::vector<uint8_t>>& siblings,
                   const std::vector<uint8_t>& x_star,
                   const std::vector<uint8_t>& x) const
    {
        if (siblings.size() != depth_)
            throw std::invalid_argument("GGM_PRF::eval_punctured: siblings size mismatch");
        if (x_star.size() != depth_ || x.size() != depth_)
            throw std::invalid_argument("GGM_PRF::eval_punctured: input length mismatch");

        // Find the first level where x diverges from x*
        size_t diverge = depth_;
        for (size_t i = 0; i < depth_; i++) {
            if ((x[i] & 1) != (x_star[i] & 1)) {
                diverge = i;
                break;
            }
        }

        if (diverge == depth_)
            throw std::invalid_argument(
                "GGM_PRF::eval_punctured: x equals x* (cannot evaluate)");

        // Start from the sibling at the divergence level
        std::vector<uint8_t> node = siblings[diverge];

        // Continue down the x path from diverge+1 to depth
        for (size_t i = diverge + 1; i < depth_; i++) {
            auto expanded = prg_.eval(node);
            std::vector<uint8_t> left (expanded.begin(),
                                        expanded.begin() + node_bits_);
            std::vector<uint8_t> right(expanded.begin() + node_bits_,
                                        expanded.end());
            node = (x[i] & 1) ? right : left;
        }

        return node;
    }

    const PRG_NC0& prg() const { return prg_; }

private:
    size_t  node_bits_;
    size_t  depth_;
    PRG_NC0 prg_;
};

// ---------------------------------------------------------------------------
// GGM PRF parameters
// ---------------------------------------------------------------------------
struct GGMParams {
    size_t node_bits;  // key/output size (= n)
    size_t depth;      // input length (= λ)

    static GGMParams toy()   { return { 8,  4 }; }   // 8-bit nodes, depth 4
    static GGMParams small() { return { 16, 8 }; }   // 16-bit nodes, depth 8
    static GGMParams med()   { return { 32, 16 }; }  // 32-bit nodes, depth 16
};
