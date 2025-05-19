////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Copyright (C) 2021-2023 CryptoLab, Inc.                                    //
//                                                                            //
// - This file is part of HEaaN homomorphic encryption library.               //
// - HEaaN cannot be copied and/or distributed without the express permission //
//  of CryptoLab, Inc.                                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "HEaaN/HEaaN.hpp"

using namespace HEaaN;

u64 getNumSignDecompose(u64 prec);

void computeOptimalSevenDegree(const HomEvaluator &eval, const Ciphertext &ctxt,
                               Ciphertext &res, const u64 in_prec,
                               const u64 index, const Real scale);

int iterShiftWithPrec(u64 in_prec, int out_prec);

