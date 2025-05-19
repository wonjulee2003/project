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

// Global constants

// The range of bootstrap extended is [-2^20, 2^20].
constexpr Real MAXIMUM_FOR_BOOTSTRAP_EXTENDED = static_cast<Real>(1 << 20);

// The magnitude of constant in `constMult` should be > 1.0e-8.
constexpr u64 LOG_MAXIMUM_FOR_BOOTSTRAP_EXTENDED = 20;
constexpr Real MINIMUM_FOR_CONST_MULT = 1.0e-8;

// constant = coeff * (MINIMUM_FOR_CONST_MULT)^(exp)
// exp is the smallest integer satisfying the magnitude of coeff is > 1.0e-8.
struct DecomposedConstant {
    explicit DecomposedConstant(const Real constant);
    explicit DecomposedConstant(const Real coefficient, const u64 exponent);

    Real cnst;
    Real coeff;
    u64 exp;
};

void multDecomposedConstant(const HomEvaluator &eval, const Ciphertext &op,
                            Ciphertext &res,
                            const DecomposedConstant &dec_cnst_mult);

void multConstant(const HomEvaluator &eval, const Ciphertext &op,
                  Ciphertext &res, const Real cnst_mult);

void complexPackCtxt(const HomEvaluator &eval, const Ciphertext &op1,
                     const Ciphertext &op2, Ciphertext &res);

void makeSameLevel(const HomEvaluator &eval, Ciphertext &op1, Ciphertext &op2);


