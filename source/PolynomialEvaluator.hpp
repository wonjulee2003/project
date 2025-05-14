#ifndef __POLYEVAL__H
#define __POLYEVAL__H

#include <iostream>
#include <cmath>
#include <string>
#include <vector>

#include "HEaaN/HEaaN.hpp"

void levelDownIfNecessary(const HEaaN::HomEvaluator &eval, HEaaN::Ciphertext &ctxt,
                          const HEaaN::u64 target_level);
                          
void levelDownForBS(const HEaaN::HomEvaluator &eval, std::vector<HEaaN::Ciphertext> &cheby_b,
                    const HEaaN::u64 target_idx);
                    
void genChebyshevPolyForBS(const HEaaN::HomEvaluator &eval, const HEaaN::Ciphertext &ctxt,
                           const HEaaN::u64 num_bs, std::vector<HEaaN::Ciphertext> &cheby_b);

void genChebyshevOddPolyForBS(const HEaaN::HomEvaluator &eval, const HEaaN::Ciphertext &ctxt,
                              const HEaaN::u64 num_bs,
                              std::vector<HEaaN::Ciphertext> &cheby_b);
                              
std::vector<HEaaN::Ciphertext> genChebyshevPolyForGS(const HEaaN::HomEvaluator &eval,
                      const std::vector<HEaaN::Ciphertext> &cheby_b,
                      const HEaaN::u64 log_gs);
                      
HEaaN::Ciphertext evaluateBabyStepWithoutRescale(
    const HEaaN::HomEvaluator &eval, const std::vector<HEaaN::Real> &coeffs,
    const std::vector<HEaaN::Ciphertext> &cheby_b, const HEaaN::u64 bs_first,
    const HEaaN::u64 bs_length, const HEaaN::Real multiplier);
    
HEaaN::Ciphertext babyStepWithoutReducedCoeffs(const HEaaN::HomEvaluator &eval,
                                        const std::vector<HEaaN::Real> &coeffs,
                                        const std::vector<HEaaN::Ciphertext> &cheby_b,
                                        const HEaaN::u64 bs_first, const HEaaN::u64 bs_length,
                                        const HEaaN::Real multiplier);

std::vector<HEaaN::Real> makeReducedCoeffsForReducingLevelCost(const std::vector<HEaaN::Real> &coeffs,
                                      const HEaaN::u64 bs_first, const HEaaN::u64 bs_length,
                                      const HEaaN::u64 half_pow_of_two);

HEaaN::Ciphertext babyStepWithReducedCoeffs(const HEaaN::HomEvaluator &eval,
                                     const std::vector<HEaaN::Real> &coeffs,
                                     std::vector<HEaaN::Ciphertext> &cheby_b,
                                     const HEaaN::u64 bs_first, 
                                     const HEaaN::u64 bs_length,
                                     const HEaaN::Real multiplier);
                                     
HEaaN::Ciphertext recursiveGiantStep(const HEaaN::HomEvaluator &eval,
                              const std::vector<HEaaN::Real> &coeffs,
                              std::vector<HEaaN::Ciphertext> &cheby_b,
                              std::vector<HEaaN::Ciphertext> &cheby_g,
                              const HEaaN::u64 gs_first, const HEaaN::u64 gs_end,
                              const HEaaN::Real multiplier,
                              const bool reduce_level_cost);
                                     
HEaaN::Ciphertext evaluateChebyshevExpansion(const HEaaN::HomEvaluator &eval,
                                      const HEaaN::Ciphertext &ctxt,
                                      const std::vector<HEaaN::Real> &cheby_coeffs,
                                      const bool basisIsOdd,
                                      const HEaaN::Real multiplier);

                                    
#endif









