////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Copyright (C) 2021-2023 Crypto Lab Inc.                                    //
//                                                                            //
// - This file is part of HEaaN homomorphic encryption library.               //
// - HEaaN cannot be copied and/or distributed without the express permission //
//  of Crypto Lab Inc.                                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "HEaaN/HEaaN.hpp"

using namespace HEaaN;

void approxSignDecomposition(const HomEvaluator &eval, const Bootstrapper &btp,
                             const Ciphertext &ctxt, Ciphertext &ctxt_out,
                             const u64 in_prec, const int out_prec,
                             const Real scale);

void approxSignDecompositionPair(const HomEvaluator &eval,
                                 const Bootstrapper &btp,
                                 const Ciphertext &ctxt1,
                                 const Ciphertext &ctxt2, Ciphertext &ctxt_out1,
                                 Ciphertext &ctxt_out2, const u64 in_prec,
                                 const int out_prec, const Real scale);

void approxDiscreteEqualZeroSinc(const HomEvaluator &eval,
                                 const Bootstrapper &btp, const Ciphertext &op,
                                 Ciphertext &res, const u64 log_range,
                                 const Real prec);

void approxSqrtInverseNewton(const HomEvaluator &eval, const Bootstrapper &btp,
                             const Ciphertext &ctxt, Ciphertext &ctxt_out,
                             Real initial, u64 num_iter);

void approxInverseNewton(const HomEvaluator &eval, const Bootstrapper &btp,
                         const Ciphertext &ctxt, Ciphertext &ctxt_out,
                         Real initial, u64 num_iter);

void approxSqrtWilkes(const HomEvaluator &eval, const Bootstrapper &btp,
                      const Ciphertext &ctxt, Ciphertext &ctxt_out,
                      const u64 num_iter);

void approx8thRootInverseNewton(const HomEvaluator &eval,
                                const Bootstrapper &btp, const Ciphertext &ctxt,
                                Ciphertext &ctxt_out, const Real range,
                                Real prec);

void approxLogEForGreaterThanOne8thRoot(const HomEvaluator &eval,
                                        const Bootstrapper &btp,
                                        const Ciphertext &op, Ciphertext &res,
                                        const Real range, const Real prec,
                                        const Real multiplier);

void approxLogEForLessThanOne8thRoot(const HomEvaluator &eval,
                                     const Bootstrapper &btp,
                                     const Ciphertext &op, Ciphertext &res,
                                     const Real range, const Real prec,
                                     const Real multiplier);

void approxExponentialForLessThanZeroEuler(
    const HomEvaluator &eval, const Bootstrapper &btp, const Ciphertext &op,
    Ciphertext &res, const Real multiplier, const Real num_iter);
