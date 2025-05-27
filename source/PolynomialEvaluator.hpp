#ifndef __POLYEVAL__H
#define __POLYEVAL__H

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

#include "HEaaN/HEaaN.hpp"

#pragma once

struct Interval {
    Interval(const HEaaN::Real left, const HEaaN::Real right)
        : left_end{left}, right_end{right} {}

    HEaaN::Real left_end;
    HEaaN::Real right_end;
};

using InputInterval = Interval;
using DiffInterval = Interval;

enum PolynomialBasis { basic, odd };

struct ChebyshevCoefficients {
    ChebyshevCoefficients(std::vector<HEaaN::Real> coefficient, HEaaN::u64 num_baby_step,
                          PolynomialBasis base = PolynomialBasis::basic)
        : coeffs{std::move(coefficient)}, basis{base} {
        num_bs = (num_baby_step > 0) ? num_baby_step : coeffs.size();
        num_gs = static_cast<HEaaN::u64>(std::ceil(static_cast<HEaaN::Real>(coeffs.size()) /
                                            static_cast<HEaaN::Real>(num_bs)));
        log_bs = static_cast<HEaaN::u64>(std::ceil(std::log2(num_bs)));
        log_gs = static_cast<HEaaN::u64>(std::ceil(std::log2(num_gs)));
        level_cost = log_bs + log_gs;
    }

    std::vector<HEaaN::Real> coeffs;
    HEaaN::u64 num_bs;
    HEaaN::u64 num_gs;
    HEaaN::u64 log_bs;
    HEaaN::u64 log_gs;
    HEaaN::u64 level_cost;
    PolynomialBasis basis;
};

// method

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
                                   
HEaaN::Ciphertext transformForChebyshev(const HEaaN::HomEvaluator &eval,
                                 const HEaaN::Ciphertext &ctxt,
                                 const InputInterval &input_interval);                           

HEaaN::Ciphertext evaluateChebyshevExpansion(const HEaaN::HomEvaluator &eval,
                                      const HEaaN::Ciphertext &ctxt,
                                      const ChebyshevCoefficients &cheby_coeffs,
                                      const HEaaN::Real multiplier);
                
void modOperate(  std::vector<HEaaN::Real> &coeffs, 
                  const HEaaN::u64 start, const HEaaN::u64 mid, 
                  const HEaaN::u64 end, const HEaaN::u64 deg);                                   

HEaaN::Ciphertext recursiveGS( const HEaaN::HomEvaluator &eval,
                        std::vector<HEaaN::Real> &coeffs,
                        std::vector<HEaaN::Ciphertext> &cheby_b,
                        std::vector<HEaaN::Ciphertext> &cheby_g,
                        const HEaaN::u64 start, const int power,
                        const HEaaN::Real multiplier);

HEaaN::Ciphertext evaluateChebyshev( const HEaaN::HomEvaluator &eval,
                                     const HEaaN::Ciphertext &ctxt,
                                     ChebyshevCoefficients &cheby_coeffs,
                                     const HEaaN::Real multiplier);

#endif









