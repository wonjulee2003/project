#ifndef __SERVER__H
#define __SERVER__H

// #include <iostream>
// #include <fstream>
// #include <sstream>

// #include <random>
// #include <vector>
// #include <iomanip>
// #include <string>

#include "HEaaN/HEaaN.hpp"
// #include "HEaaN-math/HEaaN-math.hpp"
#include "utils.hpp"
#include "PolynomialEvaluator.hpp"

class Server {
public:
    // HEaaN tools
    HEaaN::Context context;
    HEaaN::KeyPack pack;
    HEaaN::HomEvaluator evaluator;
    HEaaN::Encryptor encryptor;
    HEaaN::EnDecoder endecoder;
    // HEaaN::Bootstrapper btp;
    
    // for debugging
    HEaaN::SecretKey sk;

    // file transfer
    std::string context_string;
    std::string keypack_string;
    std::stringstream data_stream;
    std::string sk_string;
    
    // methods
    Server(const string &string1, string &string2, stringstream &&data_stream, string &string4);
    int server_hash(int mu, int bin_index, int bitLength);
    void send_result(HEaaN::Ciphertext &ctxt);
    std::stringstream server_computation(Params &params);
    void compute_sum(HEaaN::Ciphertext &ctxt);

    std::stringstream server_computation_time(Params &params);
    std::stringstream serverMultipleLabelComp(Params &params);
    std::stringstream server_computation_time0(Params &params);
};

#endif