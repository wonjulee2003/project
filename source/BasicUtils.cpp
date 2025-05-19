////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Copyright (C) 2021-2023 Crypto Lab Inc.                                    //
//                                                                            //
// - This file is part of HEaaN homomorphic encryption library.               //
// - HEaaN cannot be copied and/or distributed without the express permission //
//  of Crypto Lab Inc.                                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "BasicUtils.hpp"

DecomposedConstant::DecomposedConstant(const Real constant) : cnst(constant) {
    coeff = std::abs(cnst);
    exp = 0;
    while (coeff < 1.0) {
        if (coeff > 1.0e-8)
            break;

        exp++;
        coeff /= MINIMUM_FOR_CONST_MULT;
    }

    coeff = (cnst > 0.0) ? coeff : -coeff;
}

DecomposedConstant::DecomposedConstant(const Real coefficient,
                                       const u64 exponent)
    : coeff(coefficient), exp(exponent) {
    cnst = coeff;
    for (u64 i = 0; i < exp; i++)
        cnst *= MINIMUM_FOR_CONST_MULT;
}

void multDecomposedConstant(const HomEvaluator &eval, const Ciphertext &op,
                            Ciphertext &res,
                            const DecomposedConstant &dec_cnst_mult) {
    eval.mult(op, dec_cnst_mult.coeff, res);
    for (u64 i = 0; i < dec_cnst_mult.exp; i++) {
        // eval.mult(res, 1.0e-8, res) is the same as multiplying zero.
        // So, multWithoutRescale is used here.
        eval.multWithoutRescale(res, MINIMUM_FOR_CONST_MULT, res);
        eval.rescale(res);
    }
}

void multConstant(const HomEvaluator &eval, const Ciphertext &op,
                  Ciphertext &res, const Real cnst_mult) {
    multDecomposedConstant(eval, op, res, DecomposedConstant(cnst_mult));
}

void complexPackCtxt(const HomEvaluator &eval, const Ciphertext &op1,
                     const Ciphertext &op2, Ciphertext &res) {
    eval.multImagUnit(op2, res);
    eval.add(res, op1, res);
}

void makeSameLevel(const HomEvaluator &eval, Ciphertext &op1, Ciphertext &op2) {
    u64 level1 = op1.getLevel();
    u64 level2 = op2.getLevel();

    if (level1 < level2) {
        eval.levelDown(op2, level1, op2);
    }
    if (level2 < level1) {
        eval.levelDown(op1, level2, op1);
    }
}

