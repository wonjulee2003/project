////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Copyright (C) 2021-2023 Crypto Lab Inc.                                    //
//                                                                            //
// - This file is part of HEaaN homomorphic encryption library.               //
// - HEaaN cannot be copied and/or distributed without the express permission //
//  of Crypto Lab Inc.                                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "PolynomialEvaluator.hpp"

using namespace HEaaN;

void levelDownIfNecessary(const HomEvaluator &eval, Ciphertext &ctxt,
                          const u64 target_level) {
    if (ctxt.getLevel() > target_level) {
        eval.levelDown(ctxt, target_level, ctxt);
    }
}

// Level down for baby-step
void levelDownForBS(const HomEvaluator &eval, std::vector<Ciphertext> &cheby_b,
                    const u64 target_idx) {
    const u64 target_level = cheby_b[target_idx].getLevel();
    for (u64 i = 1; i < target_idx; i++) {
        levelDownIfNecessary(eval, cheby_b[i], target_level);
    }
}

// Generate Chebyshev polynomials for baby-step
// cheby_b[i] : i-th degree Chebyshev polynomial
// cheby_b = [T_0, T_1, T_2, ..., T_{num_bs - 1}]
//           (T_0 is an empty ciphertext)
void genChebyshevPolyForBS(const HomEvaluator &eval, const Ciphertext &ctxt,
                           const u64 num_bs, std::vector<Ciphertext> &cheby_b) {
    cheby_b.reserve(num_bs);
    cheby_b.emplace_back(ctxt);
    cheby_b.emplace_back(ctxt);

    // Compute Chebyshev polynomials
    // For baby-step
    for (u64 i = 2; i < num_bs; i++) {
        cheby_b.emplace_back(ctxt);
        u64 half_i = i >> 1;
        if (i & 1) {
            eval.mult(cheby_b[half_i + 1], cheby_b[half_i], cheby_b[i]);
            eval.multInteger(cheby_b[i], 2, cheby_b[i]);
            eval.sub(cheby_b[i], cheby_b[1], cheby_b[i]);
        } else {
            eval.square(cheby_b[half_i], cheby_b[i]);
            eval.multInteger(cheby_b[i], 2, cheby_b[i]);
            eval.sub(cheby_b[i], 1, cheby_b[i]);
        }
    }
}

// Generate Odd degree Chebyshev polynomials for baby-step
// cheby_b[i] : i-th degree Chebyshev polynomial
// cheby_b = [T_0, T_1, T_2, ..., T_{num_bs - 1}]
//           (if i is even and i \neq 2^j then T_i is an empty ciphertext)
// For the odd degree of chebyshev, we do not need to generate even degree
// chebyshev basis So we generate odd degree and 2^j degree of chebyshev
// basises. 2^i degree to use higher degree odd basis.
void genChebyshevOddPolyForBS(const HomEvaluator &eval, const Ciphertext &ctxt,
                              const u64 num_bs,
                              std::vector<Ciphertext> &cheby_b) {
    cheby_b.reserve(num_bs);
    cheby_b.emplace_back(ctxt);
    cheby_b.emplace_back(ctxt);

    // For baby-step
    u64 half_pow = 1;
    for (u64 i = 2; i < num_bs; i++) {
        cheby_b.emplace_back(ctxt);
        if (i & 1) {
            eval.mult(cheby_b[i - half_pow], cheby_b[half_pow], cheby_b[i]);
            eval.multInteger(cheby_b[i], 2, cheby_b[i]);
            eval.sub(cheby_b[i], cheby_b[2 * half_pow - i], cheby_b[i]);
        } else if (2 * half_pow == i) {
            eval.square(cheby_b[half_pow], cheby_b[i]);
            eval.multInteger(cheby_b[i], 2, cheby_b[i]);
            eval.sub(cheby_b[i], 1, cheby_b[i]);
            half_pow *= 2;
        }
    }
}

// Generate Chebyshev polynomials for giant-step
// cheby_g[i] : (2^i * num_bs)-th degree Chebyshev Polynomial
// cheby_g = [T_{2^0 * num_bs}, T_{2^1 * num_bs}, T_{2^2 * num_bs}, ...,
//            T_{2^{log_gs} * num_bs}]
std::vector<Ciphertext> genChebyshevPolyForGS(const HomEvaluator &eval,
                      const std::vector<Ciphertext> &cheby_b,
                      const u64 log_gs) {
    std::vector<Ciphertext> cheby_g;

    const u64 num_bs = cheby_b.size();
    const u64 half_bs = num_bs >> 1;

    cheby_g.reserve(log_gs);
    for (u64 i = 0; i <= log_gs; i++)
        cheby_g.emplace_back(eval.getContext());

    if (num_bs & 1) {
        eval.mult(cheby_b[half_bs + 1], cheby_b[half_bs], cheby_g[0]);
        eval.multInteger(cheby_g[0], 2, cheby_g[0]);
        eval.sub(cheby_g[0], cheby_b[1], cheby_g[0]);
    } else {
        eval.square(cheby_b[half_bs], cheby_g[0]);
        eval.multInteger(cheby_g[0], 2, cheby_g[0]);
        eval.sub(cheby_g[0], 1, cheby_g[0]);
    }

    for (u64 i = 1; i < log_gs; i++) {
        eval.square(cheby_g[i - 1], cheby_g[i]);
        eval.multInteger(cheby_g[i], 2, cheby_g[i]);
        eval.sub(cheby_g[i], 1, cheby_g[i]);
    }

    return cheby_g;
}

// Compute baby-step without reduced coefficients.
// res = babyStep(coeffs, bs_first, bs_length)
//     = Σ coeffs[bs_first + j] * cheby_b[j]   (0 <= j < bs_length)
Ciphertext evaluateBabyStepWithoutRescale(
    const HomEvaluator &eval, const std::vector<Real> &coeffs,
    const std::vector<Ciphertext> &cheby_b, const u64 bs_first,
    const u64 bs_length, const Real multiplier) {
    Ciphertext res(eval.getContext()), tmp(eval.getContext());
    auto coeff = coeffs.begin() + static_cast<int>(bs_first + bs_length) - 1;
    auto ctxt_b = cheby_b.begin() + static_cast<int>(bs_length) - 1;

    eval.multWithoutRescale(*ctxt_b, multiplier * *coeff, res);
    coeff--;
    ctxt_b--;
    while (ctxt_b != cheby_b.begin()) {
        const Real cnst = multiplier * *coeff;
        if (std::abs(cnst) > 1.0e-8) {
            eval.multWithoutRescale(*ctxt_b, cnst, tmp);
            eval.add(res, tmp, res);
        }
        coeff--;
        ctxt_b--;
    }
    eval.add(res, multiplier * *coeff, res);

    return res;
}

// Compute baby-step without reduced coefficients.
// res = babyStep(coeffs, bs_first, bs_length)
//     = Σ coeffs[bs_first + j] * cheby_b[j]   (0 <= j < bs_length)
Ciphertext babyStepWithoutReducedCoeffs(const HomEvaluator &eval,
                                        const std::vector<Real> &coeffs,
                                        const std::vector<Ciphertext> &cheby_b,
                                        const u64 bs_first, const u64 bs_length,
                                        const Real multiplier) {
    Ciphertext res = evaluateBabyStepWithoutRescale(
        eval, coeffs, cheby_b, bs_first, bs_length, multiplier);
    eval.rescale(res);

    return res;
}

// Return reduced coefficients for level cost optimization.
// For details, see
// https://www.notion.so/Chebyshev-evaluator-Not-yet-0ffd9a2888cd4c90839cfc897ce6cf24
std::vector<Real>
makeReducedCoeffsForReducingLevelCost(const std::vector<Real> &coeffs,
                                      const u64 bs_first, const u64 bs_length,
                                      const u64 half_pow_of_two) {
    std::vector<Real> reduced_coeffs(bs_length);
    std::copy(coeffs.begin() + static_cast<long>(bs_first),
              coeffs.begin() + static_cast<long>(bs_first + bs_length),
              reduced_coeffs.begin());
    u64 pow_of_two = half_pow_of_two;
    u64 base_idx = pow_of_two;

    constexpr u64 GAP_MIN = 2;
    // base_idx = 2^k + 2^{k+1} + ... while bs_length - base_idx >= GAP_MIN
    // In this case pow_of_two >= GAP_MIN is guaranteed
    while (base_idx + GAP_MIN <= bs_length) {
        for (u64 j = 1; j < bs_length - base_idx; j++) {
            reduced_coeffs[base_idx - j] -= reduced_coeffs[base_idx + j];
            reduced_coeffs[base_idx + j] *= 2;
        }
        pow_of_two >>= 1;
        base_idx += pow_of_two;
    }

    return reduced_coeffs;
}

// Compute baby-step with reduced coefficients for reducing level cost.
// ex) bs_length = 8
// res = babyStepWith(coeffs, bs_first, bs_length)
//     = babyStepWithout(reduced_coeffs, 0, 4)
//       + cheby_b[4] * babyStepWithout(reduced_coeffs, 4, 2)
//       + cheby_b[4] * cheby_b[2] * babyStepWithout(reduced_coeffs, 6, 2)
Ciphertext babyStepWithReducedCoeffs(const HomEvaluator &eval,
                                     const std::vector<Real> &coeffs,
                                     std::vector<Ciphertext> &cheby_b,
                                     const u64 bs_first, const u64 bs_length,
                                     const Real multiplier) {
    const u64 half_pow_of_two =
        1 << (static_cast<u64>(std::floor(std::log2(cheby_b.size() - 1))));
    const std::vector<Real> reduced_coeffs =
        makeReducedCoeffsForReducingLevelCost(coeffs, bs_first, bs_length,
                                              half_pow_of_two);

    // bs_first_tmp: 2^k + 2^{k-1} + 2^{k-2} + ...
    //               while bs_length - bs_first_tmp >= GAP_MIN
    constexpr u64 GAP_MIN = 2;
    u64 pow_of_two = half_pow_of_two;
    u64 bs_first_tmp = 0;
    while (bs_first_tmp + GAP_MIN <= bs_length) {
        bs_first_tmp += pow_of_two;
        if (bs_first_tmp + GAP_MIN > bs_length) {
            bs_first_tmp -= pow_of_two;
            break;
        }
        pow_of_two >>= 1;
    }

    // Partition of [0, bs_length]
    // [0, a_0], [a_0, a_1], ...,[a_{l-1}, a_l], [a_l, bs_length]
    // where a_j = a + 2^k + 2^{k-1} + ... + 2^{k-j} and
    // a_l is the maximum value such that bs_length - a_l >= GAP_MIN
    levelDownForBS(eval, cheby_b, pow_of_two);
    Ciphertext res = babyStepWithoutReducedCoeffs(
        eval, reduced_coeffs, cheby_b, bs_first_tmp, bs_length - bs_first_tmp,
        multiplier);

    while (bs_first_tmp != 0) {
        pow_of_two <<= 1;
        bs_first_tmp -= pow_of_two;
        eval.multWithoutRescale(cheby_b[pow_of_two], res, res);
        levelDownForBS(eval, cheby_b, pow_of_two);
        Ciphertext tmp = evaluateBabyStepWithoutRescale(eval, reduced_coeffs,
                                                        cheby_b, bs_first_tmp,
                                                        pow_of_two, multiplier);
        eval.add(res, tmp, res);
        eval.rescale(res);
    }

    return res;
}

// recursive giant-step
// If gs_end - gs_first == 1
// recursiveGiantStep(gs_first, gs_end) = babyStep(gs_first * num_bs)
// If gs_end - gs_first > 1
// recursiveGiantStep(gs_first, gs_end) =
//     recursiveGiantStep(gs_first, gs_first + 2^k) +
//     cheby_g[k] * recursiveGiantStep(gs_first + 2^k, gs_end)
Ciphertext recursiveGiantStep(const HomEvaluator &eval,
                              std::vector<Real> &coeffs,
                              std::vector<Ciphertext> &cheby_b,
                              std::vector<Ciphertext> &cheby_g,
                              const u64 gs_first, const u64 gs_end,
                              const Real multiplier,
                              const bool reduce_level_cost) {
    const u64 gs_length = gs_end - gs_first;

    if (gs_length == 1) {
        const u64 bs_first = gs_first * cheby_b.size();
        const u64 bs_length =
            std::min(cheby_b.size(), coeffs.size() - bs_first);

        if (reduce_level_cost) {
            Ciphertext res = babyStepWithReducedCoeffs(
                eval, coeffs, cheby_b, bs_first, bs_length, multiplier);

            levelDownForBS(eval, cheby_b, cheby_b.size() - 1);
            return res;
        }

        return babyStepWithoutReducedCoeffs(eval, coeffs, cheby_b, bs_first,
                                            bs_length, multiplier);
    }

    // gs_idx = k such that 2^k < gs_length <= 2^{k+1}.
    const u64 gs_idx = static_cast<u64>(std::ceil(std::log2(gs_length))) - 1;
    const u64 gs_mid = gs_first + (1 << gs_idx);

    const u64 bs_first_of_mid = gs_mid * cheby_b.size();
    if (bs_first_of_mid == coeffs.size() - 1) {
        if (reduce_level_cost) {
            levelDownForBS(eval, cheby_b, cheby_b.size() - 1);
        }

        Ciphertext res =
            recursiveGiantStep(eval, coeffs, cheby_b, cheby_g, gs_first, gs_mid,
                               multiplier, false);

        const Real coeff = coeffs[bs_first_of_mid];
        if (std::abs(coeff) > 1.0e-8) {
            Ciphertext tmp(eval.getContext());
            eval.mult(cheby_g[gs_idx], multiplier * coeff, tmp);
            eval.add(res, tmp, res);
        }

        return res;
    }

    Ciphertext res = recursiveGiantStep(eval, coeffs, cheby_b, cheby_g, gs_mid,
                                        gs_end, multiplier, reduce_level_cost);

    if (gs_idx > 0) {
        const u64 target_level =
            std::min(res.getLevel(), cheby_g[gs_idx].getLevel());
        levelDownIfNecessary(eval, res, target_level);
        levelDownIfNecessary(eval, cheby_g[gs_idx], target_level);
        eval.tensor(res, cheby_g[gs_idx], res);

        const u64 gs_idx_half = gs_idx - 1;
        const u64 gs_quart = gs_first + (1 << gs_idx_half);

        Ciphertext tmp2 =
            recursiveGiantStep(eval, coeffs, cheby_b, cheby_g, gs_quart, gs_mid,
                               multiplier, false);

        levelDownIfNecessary(eval, tmp2, target_level);
        if ((gs_end & (gs_end - 1)) == 0) {
            levelDownIfNecessary(eval, cheby_g[gs_idx_half], target_level);
            eval.tensor(tmp2, cheby_g[gs_idx_half], tmp2);
        } else {
            Ciphertext tmp = cheby_g[gs_idx_half];
            eval.levelDown(tmp, target_level, tmp);
            eval.tensor(tmp2, tmp, tmp2);
        }

        eval.add(res, tmp2, res);
        eval.relinearize(res, res);
        eval.rescale(res);

        Ciphertext tmp1 =
            recursiveGiantStep(eval, coeffs, cheby_b, cheby_g, gs_first,
                               gs_quart, multiplier, false);
        eval.add(res, tmp1, res);
    } else {
        eval.mult(res, cheby_g[gs_idx], res);
        Ciphertext tmp =
            recursiveGiantStep(eval, coeffs, cheby_b, cheby_g, gs_first, gs_mid,
                               multiplier, false);
        eval.add(res, tmp, res);
    }

    return res;
}

Ciphertext transformForChebyshev(const HomEvaluator &eval,
                                 const Ciphertext &ctxt,
                                 const InputInterval &input_interval) {
    Ciphertext res{ctxt};
    const Real left_end = input_interval.left_end;
    const Real right_end = input_interval.right_end;
    if (right_end + left_end != 0)
        eval.sub(res, (right_end + left_end) / 2.0, res);
    if (right_end - left_end != 0)
        eval.mult(res, 2.0 / (right_end - left_end), res);

    return res;
}

Ciphertext evaluateChebyshevExpansion(const HomEvaluator &eval,
                                      const Ciphertext &ctxt,
                                      ChebyshevCoefficients &cheby_coeffs,
                                      const Real multiplier) {
    if (ctxt.getLevel() < cheby_coeffs.level_cost) {
        throw RuntimeException("[evaluateChebyshevExpansion] input level must "
                               "be greater than the level cost " +
                               std::to_string(cheby_coeffs.level_cost));
    }

    // Preparation
    const u64 num_bs = cheby_coeffs.num_bs;
    const u64 num_gs = cheby_coeffs.num_gs;
    const u64 log_gs = cheby_coeffs.log_gs;

    std::vector<Ciphertext> cheby_b;
    // Generate Chebyshev polynomials for baby-step
    if (cheby_coeffs.basis == PolynomialBasis::basic ||
        (num_bs & (num_bs - 1)) != 0)
        genChebyshevPolyForBS(eval, ctxt, num_bs, cheby_b);
    else if (cheby_coeffs.basis == PolynomialBasis::odd ||
             (num_bs & (num_bs - 1)) == 0)
        genChebyshevOddPolyForBS(eval, ctxt, num_bs, cheby_b);

    Ciphertext res(eval.getContext());

    if (num_gs == 1) {
        res = babyStepWithReducedCoeffs(eval, cheby_coeffs.coeffs, cheby_b, 0,
                                        cheby_coeffs.coeffs.size(), multiplier);
    } else {
        // Generate Chebyshev polynomials for giant-step
        std::vector<Ciphertext> cheby_g =
            genChebyshevPolyForGS(eval, cheby_b, log_gs);

        const bool reduce_level_cost = (num_gs & (num_gs - 1)) == 0;
        if (!reduce_level_cost) {
            levelDownForBS(eval, cheby_b, cheby_b.size() - 1);
        }

        res = recursiveGiantStep(eval, cheby_coeffs.coeffs, cheby_b, cheby_g, 0,
                                 num_gs, multiplier, reduce_level_cost);
    }

    return res;
}

void modOperate(  std::vector<Real> &coeffs, 
                  const u64 start, const u64 mid, const u64 deg){             
    u64 end = std::min(mid + deg, coeffs.size());

    for(u64 i = mid+1; i < end; i++){
        coeffs[start + deg-i+mid] -= coeffs[i]; 
        coeffs[i] *= 2;               
    }
}

Ciphertext recursiveGS( const HomEvaluator &eval,
                        std::vector<Real> &coeffs,
                        std::vector<Ciphertext> &cheby_b,
                        std::vector<Ciphertext> &cheby_g,
                        const u64 start, const int power, const u64 num_bs,
                        const Real multiplier){                       
    if(power < 0){
        Ciphertext res = babyStepWithReducedCoeffs(eval, coeffs, cheby_b, start,
                                        num_bs, multiplier);

        return res;                                        
    }
    
    u64 deg = (1 << power) * num_bs;
    u64 mid = start + deg;
    modOperate(coeffs, start, mid, deg);

    Ciphertext rest = recursiveGS(eval, coeffs, cheby_b, cheby_g,
                                        start, power-1, num_bs, multiplier);
    Ciphertext quotient = recursiveGS(eval, coeffs, cheby_b, cheby_g,
                                        mid, power-1, num_bs, multiplier);                                        

    Ciphertext res(eval.getContext());
    eval.mult(quotient, cheby_g[power], quotient);
    eval.add(rest, quotient, res);         

    return res;                                            
}

Ciphertext evaluateChebyshev( const HomEvaluator &eval,
                              const Ciphertext &ctxt,
                              ChebyshevCoefficients &cheby_coeffs,
                              const Real multiplier) {
    if (ctxt.getLevel() < cheby_coeffs.level_cost) {
        throw RuntimeException("[evaluateChebyshevExpansion] input level must "
                               "be greater than the level cost " +
                               std::to_string(cheby_coeffs.level_cost));
    }

    // Preparation
    const u64 num_bs = cheby_coeffs.num_bs;
    const u64 num_gs = cheby_coeffs.num_gs;
    const int power = cheby_coeffs.log_gs;

    std::vector<Ciphertext> cheby_b;
    // Generate Chebyshev polynomials for baby-step
    if (cheby_coeffs.basis == PolynomialBasis::basic ||
        (num_bs & (num_bs - 1)) != 0)
        genChebyshevPolyForBS(eval, ctxt, num_bs, cheby_b);
    else if (cheby_coeffs.basis == PolynomialBasis::odd ||
             (num_bs & (num_bs - 1)) == 0)
        genChebyshevOddPolyForBS(eval, ctxt, num_bs, cheby_b);

    Ciphertext res(eval.getContext());

    if (num_gs == 1) {
        res = babyStepWithReducedCoeffs(eval, cheby_coeffs.coeffs, cheby_b, 0,
                                        cheby_coeffs.coeffs.size(), multiplier);

        // res = babyStepWithoutReducedCoeffs(eval, cheby_coeffs.coeffs, cheby_b, 0,
        //                                 cheby_coeffs.coeffs.size(), multiplier);                                        
    } else {
        // Generate Chebyshev polynomials for giant-step
        std::vector<Ciphertext> cheby_g =
            genChebyshevPolyForGS(eval, cheby_b, power);

        const bool reduce_level_cost = (num_gs & (num_gs - 1)) == 0;
        if (!reduce_level_cost) {
            levelDownForBS(eval, cheby_b, cheby_b.size() - 1);
        }

        res = recursiveGS(eval, cheby_coeffs.coeffs, cheby_b, cheby_g, 0,
                                power-1, num_bs, multiplier);
    }

    return res;
}