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

///
///\brief Sort the encrypted message in a ciphertext
///\param[in] eval HomEvaluator
///\param[in] btp Bootstrapper
///\param[in] op Input Ciphertext
///\param[out] res Ciphertext
///\param[in] n The number to be sorted
///\param[in] ascent Ascending or descending
///\param[in] only_last_stage Do the last stage only
///\details Sort the encrypted ciphertext `op` over the first `n` slots.
///
/// This works when `op` lies in \f$ [-0.5, 0.5] \f$.
///
/// The error for plaintext is \f$ 1.0736e-05 \f$.
///\pre Bootstrap is required.

void unitSort(const HomEvaluator &eval, const Bootstrapper &btp,
              const Ciphertext &op, Ciphertext &res, int n, bool ascent,
              bool only_last_stage = false);

