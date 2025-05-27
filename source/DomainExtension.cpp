////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Copyright (C) 2021-2023 Crypto Lab Inc.                                    //
//                                                                            //
// - This file is part of HEaaN homomorphic encryption library.               //
// - HEaaN cannot be copied and/or distributed without the express permission //
//  of Crypto Lab Inc.                                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "DomainExtension.hpp"


Ciphertext approxDomainExtension(const HomEvaluator &eval,
                                 const Bootstrapper &btp,
                                 const Ciphertext &ctxt, const Real base_range,
                                 const Real extended_range,
                                 const Real domain_extension_rate) {
    // (yj):
    // InverseWide's target level is 4.
    // SigmoidWide's target level is 3.
    if (extended_range <= base_range)
        throw RuntimeException(
            "[approxDomainExtension] extended_range should be > base_range");

    const u64 target_level = (base_range > 1.0)
                                 ? 1 + btp.getMinLevelForBootstrap()
                                 : btp.getMinLevelForBootstrap();

    ////////////
    // SETUP  //
    ////////////

    const u64 domain_extension_order =
        static_cast<u64>(std::ceil(std::log2(extended_range / base_range) /
                                   std::log2(domain_extension_rate)));
    u64 bound = std::ceil(extended_range);
    // This yields an overflow when computing (bound * bound) with bound = 2^40.
    Real div = 1.0 / (static_cast<Real>(bound) * static_cast<Real>(bound));
    const Real div_bound = 1.0 / (1 << 20);

    Real sqrt_div = 1.0 / static_cast<Real>(bound);
    const u64 sqrt_inv_div_bound = 1 << 10;
    u64 level_cost = 3;

    /////////////////////////
    // DE - ITERATION PART //
    /////////////////////////

    Ciphertext ctxt_x{ctxt};
    Ciphertext x_stable(ctxt);
    Ciphertext tmp(eval.getContext());
    for (u64 i = 0; i < domain_extension_order; i++) {
        level_cost = (div < div_bound) ? 3 : 2;

        if (ctxt_x.getLevel() < level_cost + target_level + 1) {

            eval.multWithoutRescale(x_stable, 0.5 / static_cast<Real>(bound),
                                    ctxt_x);
            eval.rescale(ctxt_x);
            btp.bootstrap(ctxt_x, ctxt_x);
            eval.multInteger(ctxt_x, bound, ctxt_x);
            eval.conjugate(ctxt_x, tmp);
            eval.add(ctxt_x, tmp, ctxt_x);
        }

        div *= domain_extension_rate * domain_extension_rate;
        sqrt_div *= domain_extension_rate;

        Ciphertext x_copy(eval.getContext()), pow2(eval.getContext()),
            x_tmp(eval.getContext());

        x_copy = ctxt_x;
        x_tmp = x_copy;
        Real div_tmp = div;

        // This weird code is due to error expansion when you do (x+e)^3 =
        // x^3+3x^2e (x+e)^2 is better
        if (div_tmp < div_bound) {
            // (yj): pls find optimized multipliers ...
            eval.multWithoutRescale(x_copy, sqrt_div * sqrt_inv_div_bound,
                                    x_copy);
            eval.rescale(x_copy);
            div_tmp = div_bound;
        }

        Real coeff = -4.0 * div_tmp / 27.0;

        eval.square(x_copy, pow2);
        eval.multWithoutRescale(x_tmp, coeff, tmp);
        eval.rescale(tmp);
        eval.mult(pow2, tmp, tmp);
        eval.add(ctxt_x, tmp, ctxt_x);
        eval.add(x_stable, tmp, x_stable);

        bound = static_cast<u64>(
            std::ceil(static_cast<Real>(bound) / domain_extension_rate));
    }

    return x_stable;
}

void approxDomainExtensionPair(const HomEvaluator &eval,
                               const Bootstrapper &btp, const Ciphertext &ctxt1,
                               const Ciphertext &ctxt2, Ciphertext &ctxt_out1,
                               Ciphertext &ctxt_out2, const Real base_range,
                               const Real extended_range,
                               const Real domain_extension_rate) {
    // (yj):
    // InverseWide's target level is 4.
    // SigmoidWide's target level is 3.
    if (extended_range <= base_range)
        throw RuntimeException("[approxDomainExtensionPair] extended_range "
                               "should be > base_range");

    const u64 target_level = (base_range > 1.0)
                                 ? 1 + btp.getMinLevelForBootstrap()
                                 : btp.getMinLevelForBootstrap();

    ////////////
    // SETUP  //
    ////////////

    const u64 domain_extension_order =
        static_cast<u64>(std::ceil(std::log2(extended_range / base_range) /
                                   std::log2(domain_extension_rate)));
    u64 bound = std::ceil(extended_range);
    // This yields an overflow when computing (bound * bound) with bound = 2^40.
    Real div = 1.0 / (static_cast<Real>(bound) * static_cast<Real>(bound));
    const Real div_bound = 1.0 / (1 << 20);

    Real sqrt_div = 1.0 / static_cast<Real>(bound);
    const u64 sqrt_inv_div_bound = 1 << 10;
    u64 level_cost = 3;

    /////////////////////////
    // DE - ITERATION PART //
    /////////////////////////

    Ciphertext ctxt_x1{ctxt1};
    Ciphertext ctxt_x2{ctxt2};

    u64 level1 = ctxt_x1.getLevel();
    u64 level2 = ctxt_x2.getLevel();
    if (level1 < level2) {
        eval.levelDown(ctxt_x2, level1, ctxt_x2);
    }
    if (level2 < level1) {
        eval.levelDown(ctxt_x1, level2, ctxt_x1);
    }
    ctxt_out1 = ctxt_x1;
    ctxt_out2 = ctxt_x2;
    Ciphertext tmp(eval.getContext());
    Ciphertext x_copy1(eval.getContext()), x_tmp1(eval.getContext());
    Ciphertext x_copy2(eval.getContext()), x_tmp2(eval.getContext());
    Ciphertext pow2(eval.getContext());
    for (u64 i = 0; i < domain_extension_order; i++) {
        level_cost = (div < div_bound) ? 3 : 2;

        if (ctxt_x1.getLevel() < level_cost + target_level + 1) {
            eval.multImagUnit(ctxt_out2, tmp);
            eval.add(ctxt_out1, tmp, tmp);

            eval.multWithoutRescale(tmp, 0.5 / static_cast<Real>(bound),
                                    ctxt_x1);
            eval.rescale(ctxt_x1);
            btp.bootstrap(ctxt_x1, ctxt_x1, ctxt_x2);

            eval.multInteger(ctxt_x1, bound, ctxt_x1);
            eval.conjugate(ctxt_x1, tmp);
            eval.add(ctxt_x1, tmp, ctxt_x1);

            eval.multInteger(ctxt_x2, bound, ctxt_x2);
            eval.conjugate(ctxt_x2, tmp);
            eval.add(ctxt_x2, tmp, ctxt_x2);
        }

        div *= domain_extension_rate * domain_extension_rate;
        sqrt_div *= domain_extension_rate;

        x_copy1 = ctxt_x1;
        x_tmp1 = x_copy1;
        x_copy2 = ctxt_x2;
        x_tmp2 = x_copy2;

        Real div_tmp = div;

        // This weird code is due to error expansion when you do (x+e)^3 =
        // x^3+3x^2e (x+e)^2 is better
        if (div_tmp < div_bound) {
            // (yj): pls find optimized multipliers ...
            eval.multWithoutRescale(x_copy1, sqrt_div * sqrt_inv_div_bound,
                                    x_copy1);
            eval.rescale(x_copy1);

            eval.multWithoutRescale(x_copy2, sqrt_div * sqrt_inv_div_bound,
                                    x_copy2);
            eval.rescale(x_copy2);

            div_tmp = div_bound;
        }

        Real coeff = -4.0 * div_tmp / 27.0;

        eval.square(x_copy1, pow2);
        eval.multWithoutRescale(x_tmp1, coeff, tmp);
        eval.rescale(tmp);
        eval.mult(pow2, tmp, tmp);
        eval.add(ctxt_x1, tmp, ctxt_x1);
        eval.add(ctxt_out1, tmp, ctxt_out1);

        eval.square(x_copy2, pow2);
        eval.multWithoutRescale(x_tmp2, coeff, tmp);
        eval.rescale(tmp);
        eval.mult(pow2, tmp, tmp);
        eval.add(ctxt_x2, tmp, ctxt_x2);
        eval.add(ctxt_out2, tmp, ctxt_out2);

        bound = static_cast<u64>(
            std::ceil(static_cast<Real>(bound) / domain_extension_rate));
    }
}

Ciphertext approxDomainExtensionInverse(const HomEvaluator &eval,
                                        const Ciphertext &ctxt,
                                        const u64 domain_extension_order,
                                        const Real domain_extension_rate) {
    Real coeff = (4.0 / 27.0) *
                 (1.0 - 1.0 / std::pow(domain_extension_rate,
                                       2 * domain_extension_order)) /
                 (1.0 - 1.0 / std::pow(domain_extension_rate, 2));
    Real coeff2 = coeff * (4.0 / 9.0) *
                  (1.0 - 1.0 / std::pow(domain_extension_rate,
                                        2 * (domain_extension_order + 1))) /
                  (1.0 - 1.0 / std::pow(domain_extension_rate, 4));
    Real eps = 0.1;

    Ciphertext ctxt_out(eval.getContext()), ctxt_tmp(eval.getContext()),
        pow2(eval.getContext());

    eval.square(ctxt, pow2);
    eval.square(pow2, ctxt_out);
    eval.negate(ctxt_out, ctxt_out);
    eval.mult(pow2, coeff2 / (coeff + coeff2 + eps), ctxt_tmp);
    eval.add(ctxt_out, ctxt_tmp, ctxt_out);
    eval.add(ctxt_out, coeff / (coeff + coeff2 + eps), ctxt_out);

    eval.mult(ctxt, coeff + coeff2 + eps, ctxt_tmp);
    eval.mult(ctxt_tmp, pow2, ctxt_tmp);
    eval.mult(ctxt_out, ctxt_tmp, ctxt_out);
    eval.add(ctxt_out, ctxt, ctxt_out);

    return ctxt_out;
}

