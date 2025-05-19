////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Copyright (C) 2021-2023 Crypto Lab Inc.                                    //
//                                                                            //
// - This file is part of HEaaN homomorphic encryption library.               //
// - HEaaN cannot be copied and/or distributed without the express permission //
//  of Crypto Lab Inc.                                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "ApproximationMethods.hpp"

#include "SignDecompose.hpp"
#include "BasicUtils.hpp"
#include "DomainExtension.hpp"
#include "PolynomialEvaluator.hpp"
#include <cmath>


// new implimentation of SignNative using decompose functions.
// See some details on
// https://www.notion.so/Sign-function-decomposition-9b84158fe93a434db3cb27a255058b51
void approxSignDecomposition(const HomEvaluator &eval, const Bootstrapper &btp,
                             const Ciphertext &ctxt, Ciphertext &ctxt_out,
                             u64 in_prec, const int out_prec,
                             const Real scale) {
    static const u64 ONE_ITER_COST = 3; // level consumed in each iteration
    const int num_iter_shift = iterShiftWithPrec(in_prec, out_prec);

    Ciphertext sign{ctxt};
    u64 max_iter = getNumSignDecompose(in_prec);
    u64 num_iter_new = max_iter;

    if (num_iter_shift + int(max_iter) <= 0)
        num_iter_new = 1;
    else
        num_iter_new = num_iter_shift + max_iter;

    for (u64 i = 0; i < num_iter_new; i++) {
        if (sign.getLevel() < ONE_ITER_COST + btp.getMinLevelForBootstrap())
            btp.bootstrap(sign, sign);

        if (i == num_iter_new - 1)
            computeOptimalSevenDegree(eval, sign, ctxt_out, in_prec, i, scale);
        else
            computeOptimalSevenDegree(eval, sign, sign, in_prec, i, 1.0);
    }
}

void approxSignDecompositionPair(const HomEvaluator &eval,
                                 const Bootstrapper &btp,
                                 const Ciphertext &ctxt1,
                                 const Ciphertext &ctxt2, Ciphertext &ctxt_out1,
                                 Ciphertext &ctxt_out2, u64 in_prec,
                                 const int out_prec, const Real scale) {
    static const u64 ONE_ITER_COST = 3; // level consumed in each iteration
    const int num_iter_shift = iterShiftWithPrec(in_prec, out_prec);
    Ciphertext sign1{ctxt1}, sign2{ctxt2};
    makeSameLevel(eval, sign1, sign2);
    Ciphertext tmp(eval.getContext());
    u64 max_iter = getNumSignDecompose(in_prec);
    u64 num_iter_new = max_iter;

    if (num_iter_shift + int(max_iter) <= 0)
        num_iter_new = 1;
    else
        num_iter_new = num_iter_shift + max_iter;

    for (u64 i = 0; i < num_iter_new; i++) {
        if (sign1.getLevel() < ONE_ITER_COST + btp.getMinLevelForBootstrap()) {
            complexPackCtxt(eval, sign1, sign2, tmp);
            btp.bootstrap(tmp, sign1, sign2);
        }
        if (i == num_iter_new - 1) {
            computeOptimalSevenDegree(eval, sign1, ctxt_out1, in_prec, i,
                                      scale);
            computeOptimalSevenDegree(eval, sign2, ctxt_out2, in_prec, i,
                                      scale);
        } else {
            computeOptimalSevenDegree(eval, sign1, sign1, in_prec, i, 1.0);
            computeOptimalSevenDegree(eval, sign2, sign2, in_prec, i, 1.0);
        }
    }
}

void approxDiscreteEqualZeroSinc(const HomEvaluator &eval,
                                 const Bootstrapper &btp, const Ciphertext &op,
                                 Ciphertext &res, const u64 log_range,
                                 const Real prec) {
    constexpr u64 TRIGONOMETRIC_COST = 4;
    constexpr u64 DOUBLE_ANGLE_COST = 1;

    constexpr Real SINC_COEFFS[] = {1.0, 0, -1.0 / 6, 0, 1.0 / 120};
    constexpr Real COS_COEFFS[] = {1.0, 0, -1.0 / 2, 0, 1.0 / 24};

    constexpr u64 MAX_DOUBLE_ANGLE = 12;
    const u64 total_double_angle = log_range;
    const u64 current_double_angle = (total_double_angle > MAX_DOUBLE_ANGLE)
                                         ? MAX_DOUBLE_ANGLE
                                         : total_double_angle;

    const u64 log_basae_range = current_double_angle;
    const Real base_range = std::pow(2.0, log_basae_range);

    Ciphertext x{op};
    if (log_range > log_basae_range) {
        const Real range = std::pow(2.0, log_range);
        const Real range_tmp = std::pow(2.0, 4) * base_range;
        if (range > range_tmp) {
            x = approxDomainExtension(eval, btp, x, range_tmp, range, 2.25);
            x = approxDomainExtension(eval, btp, x, base_range, range_tmp, 2.0);
        } else {
            x = approxDomainExtension(eval, btp, x, base_range, range, 2.0);
        }
    }

    eval.mult(x, REAL_PI / base_range, x);

    if (x.getLevel() < TRIGONOMETRIC_COST + btp.getMinLevelForBootstrap()) {
        btp.bootstrap(x, x);
    }

    Ciphertext sinc(eval.getContext());
    Ciphertext cos(eval.getContext());
    Ciphertext tmp(eval.getContext());

    eval.mult(x, x, tmp); // tmp = x^2

    // compute sinc((pi/2^d) * x)
    eval.mult(tmp, SINC_COEFFS[4], sinc);
    eval.add(sinc, SINC_COEFFS[2], sinc);
    eval.mult(sinc, tmp, sinc);
    eval.add(sinc, SINC_COEFFS[0], sinc);

    // compute cos((pi/2^d) * x)
    eval.mult(tmp, COS_COEFFS[4], cos);
    eval.add(cos, COS_COEFFS[2], cos);
    eval.mult(cos, tmp, cos);
    eval.add(cos, COS_COEFFS[0], cos);

    // Apply double angle formula
    for (u64 i = 0; i < current_double_angle; i++) {
        if (sinc.getLevel() <
            DOUBLE_ANGLE_COST + btp.getMinLevelForBootstrap()) {
            eval.multImagUnit(cos, cos);
            eval.add(sinc, cos, sinc);
            btp.bootstrap(sinc, sinc, cos);
        }

        // compute sinc(2x)
        eval.mult(sinc, cos, sinc);

        // compute cos(2x)
        eval.square(cos, cos);
        eval.mult(cos, 2, cos);
        eval.add(cos, -1, cos);
    }

    // apply filter
    eval.sub(sinc, 0.5, res);
    int log_prec = int(std::log2(prec));
    // For interger input the value `sinc - 0.5` can be expected +-0.5.
    // So, `in_prec=1` is used in `approxSignDecomposition`.
    approxSignDecomposition(eval, btp, res, res, 2, log_prec, 0.5);
    eval.add(res, 0.5, res);
}

// See the notion page.
// https://www.notion.so/SqrtInverseForGreaterThanOne-498516cd112f4dbba7f1c796a1d8dce7?pvs=4
void approxSqrtInverseNewton(const HomEvaluator &eval, const Bootstrapper &btp,
                             const Ciphertext &ctxt, Ciphertext &ctxt_out,
                             Real initial, u64 num_iter) {
    static const u64 ONE_ITER_COST = 3;

    Ciphertext ctx{ctxt};
    Ciphertext ctxt_x(eval.getContext());
    Ciphertext ctxt_y(eval.getContext());
    Ciphertext tmp(eval.getContext());

    // y_{n+1} = 0.5 * y_n * (3 - ctxt * y_n * y_n)
    //         = y_n (1.5 + x n)
    // x_{n+1} = x_n (1.5 + x_n)^2
    // x_n = - 0.5 * ctxt * y_n * y_n = ctx * y_n * y_n
    // ctx = - 0.5 ctxt;

    // Input range are [-2^{40},2^{40}], this is square of bootstrap extend.
    // to bootstrapExtended
    Real range = 2.0 / (initial * initial) - 1;
    if (ctxt.getLevel() < 2 + btp.getMinLevelForBootstrap()) {
        if (range > std::pow(2, 20)) {
            throw RuntimeException(
                "[approxSqrtInverseNewton] The input ciphertext level >= " +
                std::to_string(2 + btp.getMinLevelForBootstrap()) +
                " is needed for such large range");
        }
    }

    // Some computations for pre-bootstrapping accroding to input level and
    // num_iters See
    // https://www.notion.so/SqrtInverseGreaterThanOne-2513684a928e4ecfb057e7faa6be2e6c?pvs=4
    // section 5.
    Real pre_boot_time_factor = 2.0;
    u64 log_range = (u64)std::floor(std::log2(range)) + 2;
    if (range > std::pow(2, 20)) {
        log_range -= 20; // bootstrapextended range issue.
        pre_boot_time_factor += 1.0;
    }

    // getLevel() can be bigger than btp.getLevelAfterFullSlotBootstrap()
    // for just FVa.
    Real iters_per_one_boot_without_pre_boot =
        std::ceil(Real(std::min(ctx.getLevel(),
                                btp.getLevelAfterFullSlotBootstrap())) /
                  2.0) -
        1.0;
    Real iters_per_one_boot_with_pre_boot =
        std::ceil(Real(btp.getLevelAfterFullSlotBootstrap()) / 2.0) - 1;
    // getLevel() can be bigger than btp.getLevelAfterFullSlotBootstrap()
    // for just FVa.
    int n_without_pre_boot =
        int(num_iter) - std::max(int(ctx.getLevel()) -
                                     int(btp.getLevelAfterFullSlotBootstrap()),
                                 0) /
                            2;
    if (std::floor(Real(n_without_pre_boot) /
                   iters_per_one_boot_without_pre_boot) >
        pre_boot_time_factor +
            std::floor(Real(num_iter) / iters_per_one_boot_with_pre_boot)) {
        eval.mult(ctx, std::pow(2, -double(log_range)), ctx);
        if (range > std::pow(2, 20)) {
            btp.bootstrapExtended(ctx, ctx);
        } else {
            btp.bootstrap(ctx, ctx);
        }

        // removing error process.
        Ciphertext ctxt_err(eval.getContext());
        eval.multInteger(ctx, 1 << (log_range - 2), ctx);
        eval.mult(ctxt, 0.25, ctxt_err);
        eval.sub(ctxt_err, ctx, ctxt_err);

        btp.bootstrap(ctxt_err, ctxt_err);
        eval.add(ctx, ctxt_err, ctx);
        eval.negate(ctx, ctx);
        eval.conjugate(ctx, tmp);
        eval.add(ctx, tmp, ctx); // ctx= -0.5 ctxt

        multConstant(eval, ctx, ctxt_x, initial * initial);
    } else {
        multConstant(eval, ctx, ctxt_x, -0.5 * initial * initial);
        eval.mult(ctx, -0.5, ctx);
    }

    for (u64 i = 0; i < num_iter; i++) {
        if (i == 0) {
            eval.add(ctxt_x, 1.5, tmp);
            eval.mult(tmp, initial, ctxt_y);
        } else {
            if (ctxt_x.getLevel() <
                ONE_ITER_COST + btp.getMinLevelForBootstrap()) {
                btp.bootstrap(ctxt_y, ctxt_y);
                eval.mult(ctx, ctxt_y, ctxt_x);
                eval.mult(ctxt_x, ctxt_y, ctxt_x);
            } else {
                eval.square(tmp, tmp);
                eval.mult(tmp, ctxt_x, ctxt_x);
            }

            eval.add(ctxt_x, 1.5, tmp);
            eval.mult(tmp, ctxt_y, ctxt_y);
        }
    }

    ctxt_out = ctxt_y;
}

void approxInverseNewton(const HomEvaluator &eval, const Bootstrapper &btp,
                         const Ciphertext &ctxt, Ciphertext &ctxt_out,
                         Real initial, u64 num_iter) {
    const u64 one_iter_cost{1};

    Ciphertext ctxt_x(eval.getContext()), ctxt_y(eval.getContext());
    Ciphertext ctxt_z(eval.getContext()), ctxt_tmp(eval.getContext());

    if (ctxt.getLevel() < btp.getLevelAfterFullSlotBootstrap())
        btp.bootstrapExtended(ctxt, ctxt_x);
    else
        ctxt_x = ctxt;

    // y_0 = initial
    // z_0 = x * y_0
    // y_1 = y_0 * (2 - z_0)
    eval.mult(ctxt_x, initial, ctxt_z);
    eval.negate(ctxt_z, ctxt_tmp);
    eval.add(ctxt_tmp, 2.0, ctxt_tmp); // tmp = 2 - z_0
    eval.mult(ctxt_tmp, initial, ctxt_y);

    // for n > 0
    // z_n = x * y_n
    //     = z_{n-1} * (2 - z_{n-1})
    // y_{n+1} = y_n * (2 - z_n)
    for (u64 iter = 1; iter < num_iter; iter++) {
        if (ctxt_y.getLevel() < one_iter_cost + btp.getMinLevelForBootstrap()) {
            btp.bootstrap(ctxt_y, ctxt_y);
            eval.mult(ctxt_x, ctxt_y, ctxt_z);
        } else {
            eval.mult(ctxt_z, ctxt_tmp, ctxt_z);
        }

        eval.negate(ctxt_z, ctxt_tmp); // tmp = 2 - z_n
        eval.add(ctxt_tmp, 2.0, ctxt_tmp);
        eval.mult(ctxt_y, ctxt_tmp, ctxt_y);
    }

    ctxt_out = ctxt_y;
}

// This method was proposed by Wilkes in 1951.
// See https://eprint.iacr.org/2019/417 for more details.
void approxSqrtWilkes(const HomEvaluator &eval, const Bootstrapper &btp,
                      const Ciphertext &ctxt, Ciphertext &ctxt_out,
                      const u64 num_iter) {
    Ciphertext tmp(eval.getContext());
    Ciphertext tmp2(eval.getContext());
    static const u64 ONE_ITER_COST = 2;

    ctxt_out = ctxt;

    if (ctxt_out.getLevel() < 1 + btp.getMinLevelForBootstrap()) {
        btp.bootstrap(ctxt_out, ctxt_out);
    }
    eval.sub(ctxt_out, 1, tmp);
    eval.mult(tmp, 0.5, tmp);

    for (u64 i = 0; i < num_iter; i++) {
        // complex bootstrap
        if (tmp.getLevel() < ONE_ITER_COST + btp.getMinLevelForBootstrap()) {
            eval.multImagUnit(tmp, tmp2);
            eval.add(ctxt_out, tmp2, tmp2);
            btp.bootstrap(tmp2, ctxt_out, tmp);
        }

        // compute a_(k+1)
        eval.negate(tmp, tmp2);
        eval.add(tmp2, 1, tmp2);
        eval.mult(ctxt_out, tmp2, ctxt_out);

        // compute h_(k+1)
        eval.sub(tmp, 1.5, tmp2);
        eval.square(tmp, tmp);
        eval.mult(tmp, tmp2, tmp);
    }
}

// bootstrapExtended is required.
void approx8thRootInverseNewton(const HomEvaluator &eval,
                                const Bootstrapper &btp, const Ciphertext &ctxt,
                                Ciphertext &ctxt_out, const Real range,
                                Real prec) {
    static const u64 ONE_ITER_COST = 4;

    Ciphertext ctxt_x(eval.getContext());
    Ciphertext ctxt_y(eval.getContext());
    Ciphertext tmp1(eval.getContext());
    Ciphertext tmp2(eval.getContext());

    // Calculating 8th root inverse via Newton method
    // y_{n+1} = 0.125 * y_n * (9 - ctxt * y_n^8 )
    //         = 1.125 * y_n + ((-0.125 * ctxt) * y_n^4) * (y_n * y_n^4)

    // Precompute x = -0.125 * ctxt
    if (ctxt.getLevel() < btp.getLevelAfterFullSlotBootstrap()) {
        Real cnst_mult = -0.125;
        if (range > MAXIMUM_FOR_BOOTSTRAP_EXTENDED) {
            DecomposedConstant dec_cnst_mult(MAXIMUM_FOR_BOOTSTRAP_EXTENDED /
                                             range);
            if (ctxt.getLevel() <
                dec_cnst_mult.exp + 2 + btp.getMinLevelForBootstrap())
                throw RuntimeException("[approx8thRootInverseNewton] Input "
                                       "Ciphertext level is not enough.");

            multDecomposedConstant(eval, ctxt, ctxt_x, dec_cnst_mult);
            btp.bootstrapExtended(ctxt_x, ctxt_x);

            cnst_mult *= range / MAXIMUM_FOR_BOOTSTRAP_EXTENDED;
        } else {
            btp.bootstrapExtended(ctxt, ctxt_x);
        }
        eval.mult(ctxt_x, cnst_mult, ctxt_x);
    } else {
        eval.mult(ctxt, -0.125, ctxt_x);
    }

    Real init = std::pow(1.0 / range, 0.125);
    u64 num_iter = 0;
    Real e_n = 1 - init;
    while (e_n >= prec) {
        e_n = e_n - (1 - e_n) * (1 - std::pow(1 - e_n, 8)) / 8;
        num_iter++;
    }
    // Excute Newton method
    for (u64 i = 0; i < num_iter; i++) {
        if (i == 0) {
            multConstant(eval, ctxt_x, tmp1, std::pow(init, 9));
            eval.add(tmp1, 1.125 * init, ctxt_y);
        } else {
            if (ctxt_y.getLevel() <
                ONE_ITER_COST + btp.getMinLevelForBootstrap())
                btp.bootstrap(ctxt_y, ctxt_y);
            eval.square(ctxt_y, tmp2);
            eval.square(tmp2, tmp2);       // tmp2 = y^4
            eval.mult(ctxt_x, tmp2, tmp1); // tmp1 = x * y^4
            eval.mult(tmp2, ctxt_y, tmp2); // tmp2 = y^5
            eval.mult(tmp1, tmp2, tmp1);

            eval.mult(ctxt_y, 1.125, tmp2);
            eval.add(tmp1, tmp2, ctxt_y);
        }
    }

    ctxt_out = ctxt_y;
}

// bootstrapExtended is required.
void approxLogEForGreaterThanOne8thRoot(const HomEvaluator &eval,
                                        const Bootstrapper &btp,
                                        const Ciphertext &op, Ciphertext &res,
                                        const Real range, const Real prec,
                                        const Real multiplier) {
    // Polynomial approximation of log(x) for x in [1/8, 1] with polynomial
    // degree 31 and baby-step 8.
    // static const ChebyshevCoefficients CHEBYSHEV_COEFFS(
    static const ChebyshevCoefficients CHEBYSHEV_COEFFS(
        {-7.8082780989e-01, 9.5518286632e-01,  -2.2808995415e-01,
         7.2613801481e-02,  -2.5990084730e-02, 9.8854881722e-03,
         -3.8322169284e-03, 1.3320451263e-03,  -6.7670146202e-04,
         5.7454951269e-04,  -2.4694765055e-04, 1.0719062667e-04,
         -4.6866237587e-05, 2.0526763511e-05,  -8.8063178062e-06,
         3.2676337823e-06,  -9.1585420985e-07, 8.2333188876e-07,
         -3.7134848934e-07, 1.6796501089e-07,  -7.6091284922e-08,
         3.4357080504e-08,  -1.5112794325e-08, 5.7063798460e-09,
         -3.3053920258e-09, 3.0309707979e-09,  -1.3918880996e-09,
         6.4013718168e-10,  -2.9480166942e-10, 1.3594371572e-10,
         -6.2760133073e-11, 3.6922291193e-11},
        8); // 8 , 4

    // Compute x = op^{-1/8}, it is insured that x is in [1/8, 1]
    Ciphertext ctxt_x(eval.getContext());
    approx8thRootInverseNewton(eval, btp, op, ctxt_x, range,
                               0.5 * prec /
                                   (multiplier > 1.0 ? multiplier : 1.0));

    // Linear transfrom from [1/8, 1] to [-1, 1]
    if (ctxt_x.getLevel() < 1 + btp.getMinLevelForBootstrap())
        btp.bootstrap(ctxt_x, ctxt_x);
    ctxt_x = transformForChebyshev(eval, ctxt_x, InputInterval(0.125, 1.0));

    // To make res to be able to do bootstrapExtended.
    if (ctxt_x.getLevel() <
        1 + CHEBYSHEV_COEFFS.level_cost + btp.getMinLevelForBootstrap())
        btp.bootstrap(ctxt_x, ctxt_x);
    // Evaluate Chebyshev expansion
    res = evaluateChebyshevExpansion(eval, ctxt_x, CHEBYSHEV_COEFFS,
                                     -8 * multiplier);
}

// bootstrap is required.
void approxLogEForLessThanOne8thRoot(const HomEvaluator &eval,
                                     const Bootstrapper &btp,
                                     const Ciphertext &op, Ciphertext &res,
                                     const Real range, const Real prec,
                                     const Real multiplier) {
    Ciphertext ctxt_x(eval.getContext());

    // from [range, 1] to [range * ceil(1/range), ceil(1/range)]
    const u64 cnst_mult_int = std::ceil(1.0 / range);
    if (op.getLevel() < 1 + btp.getMinLevelForBootstrap()) {
        btp.bootstrap(op, ctxt_x);
        eval.multInteger(ctxt_x, cnst_mult_int, ctxt_x);
    } else {
        eval.multInteger(op, cnst_mult_int, ctxt_x);
    }
    // res = multiplier * log(ceil(1/range) * x)
    //       - multiplier * log(ceil(1/range))
    approxLogEForGreaterThanOne8thRoot(eval, btp, ctxt_x, res,
                                       static_cast<Real>(cnst_mult_int), prec,
                                       multiplier);
    eval.sub(res, multiplier * std::log(cnst_mult_int), res);
}

void approxExponentialForLessThanZeroEuler(
    const HomEvaluator &eval, const Bootstrapper &btp, const Ciphertext &op,
    Ciphertext &res, const Real multiplier, Real num_iter) {

    Ciphertext x{op};
    Ciphertext tmp1(eval.getContext()), tmp2(eval.getContext());

    const Real log_scale_factor = num_iter;
    const Real native_scale_factor =
        (x.getCurrentScaleFactor() > 1)
            ? std::pow(2.0, x.getCurrentScaleFactor())
            : std::pow(2.0, 20);
    const Real scale_factor = std::round(
        native_scale_factor /
        std::round(native_scale_factor / std::pow(2.0, log_scale_factor)));
    const int iter_num = static_cast<int>(log_scale_factor);
    // (yj): iter_num may be diff. from log_scale_factor
    // log_scale_factor is chosen considering ...
    // the range of bootstrapExtended is [-2^20, 2^20].
    static const u64 ONE_ITER_COST = 3;

    eval.add(x, scale_factor, x);
    // (yj): y = D * (1 + op / 2^N)
    // ... where D is scale_factor
    // ... and N = iter_num
    for (int i = 0; i < iter_num; i++) {

        if (x.getLevel() < ONE_ITER_COST + btp.getMinLevelForBootstrap() + 2) {
            // TODO(yj): more optimiaztion is needed.
            // bootstrapExtendedAndMultiplyConstant when?
            if (i > 0) {
                eval.killImag(x, x);
                eval.killImag(tmp2, tmp2);
                eval.multImagUnit(tmp2, tmp2);
                eval.add(x, tmp2, tmp1);
                btp.bootstrapExtended(tmp1, x, tmp2);
            } else {
                btp.bootstrapExtended(x, x);
            }
        }

        eval.square(x, tmp1); // (Dy + eps)^2

        if (i > 0) {
            eval.mult(x, tmp2, tmp2);
            eval.mult(tmp2, 2.0 / scale_factor, tmp2);
            eval.sub(tmp1, tmp2, tmp1);
        }

        eval.mult(tmp1, 1.0 / scale_factor, x); // Rescale
        eval.mult(x, scale_factor, tmp2);
        eval.sub(tmp2, tmp1, tmp2); // D * eps
    }

    eval.mult(x, multiplier / scale_factor, res);
}

