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

Ciphertext approxDomainExtension(const HomEvaluator &eval,
                                 const Bootstrapper &btp,
                                 const Ciphertext &ctxt, const Real base_range,
                                 const Real extended_range,
                                 const Real domain_extension_rate);

void approxDomainExtensionPair(const HomEvaluator &eval,
                               const Bootstrapper &btp, const Ciphertext &ctxt1,
                               const Ciphertext &ctxt2, Ciphertext &ctxt_out1,
                               Ciphertext &ctxt_out2, const Real base_range,
                               const Real extended_range,
                               const Real domain_extension_rate);

Ciphertext approxDomainExtensionInverse(const HomEvaluator &eval,
                                        const Ciphertext &ctxt,
                                        const u64 domain_extension_order,
                                        const Real domain_extension_rate);

