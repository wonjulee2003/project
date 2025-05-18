#ifndef __UTILS__H
#define __UTILS__H

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <random>
#include <optional>

#include "HEaaN/HEaaN.hpp"

using namespace std;
using namespace HEaaN;

int comb(int n, int r);
vector<int> PerfectMapping(int input, int bitLength, int hammingWeight);
void printMessage(const HEaaN::Message &msg, bool is_complex = false,
                         size_t start_num = 10, size_t end_num = 10,
                         const int precision = 10);
void fillRandomReal(HEaaN::Message &msg);
double randNum();
void fillReal(HEaaN::Message &msg, double val);

class Params {
public:
    int server_bin_size;
    int effective_bitLength;
    int ell;
    int hw; 
    
    Params(int server_bin_size, int effective_bitLength, int hw);
};

void approxInverseNewton(const HomEvaluator &eval,
                         const Ciphertext &ctxt, Ciphertext &ctxt_out,
                         Real initial, u64 num_iter);

void approxSqrtWilkes(const HomEvaluator &eval, 
                    //   const Bootstrapper &btp,
                      const Ciphertext &ctxt, Ciphertext &ctxt_out,
                      const u64 num_iter);

#endif
