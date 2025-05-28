#include <iostream>
#include <random>
#include <vector>
#include <tuple>
#include <sstream>

#include "HEaaN/HEaaN.hpp"
// #include "HEaaN-math/HEaaN-math.hpp"

#include "client.hpp"
#include "server.hpp"
#include "utils.hpp"
#include "PolynomialEvaluator.hpp"

inline std::string presetNamer(const HEaaN::ParameterPreset preset) {
    switch (preset) {
    case HEaaN::ParameterPreset::FVa:
        return "FVa";
    case HEaaN::ParameterPreset::FVb:
        return "FVb";
    case HEaaN::ParameterPreset::FGa:
        return "FGa";
    case HEaaN::ParameterPreset::FGb:
        return "FGb";
    case HEaaN::ParameterPreset::FTa:
        return "FTa";
    case HEaaN::ParameterPreset::FTb:
        return "FTb";
    case HEaaN::ParameterPreset::ST19:
        return "ST19";
    case HEaaN::ParameterPreset::ST14:
        return "ST14";
    case HEaaN::ParameterPreset::ST11:
        return "ST11";
    case HEaaN::ParameterPreset::ST8:
        return "ST8";
    case HEaaN::ParameterPreset::ST7:
        return "ST7";
    case HEaaN::ParameterPreset::SS7:
        return "SS7";
    case HEaaN::ParameterPreset::SD3:
        return "SD3";
    case HEaaN::ParameterPreset::CUSTOM:
        return "CUSTOM";
    case HEaaN::ParameterPreset::FVc:
        return "FVc";
    case HEaaN::ParameterPreset::FGd:
        return "FGd";
    case HEaaN::ParameterPreset::SGd0:
        return "SGd0";
    case HEaaN::ParameterPreset::FX:
        return "FX";
    default:
        throw std::invalid_argument("Not supported parameter");
    }
}
// Client Setup
//  ↓
// Encode Inputs (CW Encoding)
//  ↓
// Batch & Encrypt
//  ↓
// Send Ciphertexts →
//                      Server Setup
//                       ↓
//            Load Parameters & Keys
//                       ↓
//         Encode Server Inputs (CW)
//                       ↓
//          Homomorphic Computation
//                       ↓
//              Send Ciphertexts ←
//  ↓
// Decrypt & Decode


int main(void) {
    std::cout << "main project" << std::endl;
    Client client("SS7"); // SS7
    std::cout << "Parameter : " << presetNamer(client.preset) << std::endl;
    std::cout << getLogFullSlots(client.context) << std::endl;

    // for our deterministic hashing(assignment) scheme, effective_bitLength > 15
    // 2^15 slots in FGb parameters
    int effective_bitLength = 19; 
    // large hamming weight parameters causes depletion due to extensive multiplication.

    int hw = 7;

    // ======

    int client_set = 4096;
    int server_set = 1<<20;
    int N = 1 << getLogFullSlots(client.context);

    int num_bins = 1, num_balls = 1;
    // int gamma = ceil(client_set *1.27 / N);
    int gamma = 1;

    num_bins = N*gamma;
    num_balls = 3*server_set;

    // mu = ceil((float)num_balls/num_bins + 2 * sqrt((float)num_balls*log2(num_bins)/num_bins));
    // server_bin_size used for PEPSI

    // int server_bin_size = ceil((float)num_balls/num_bins + 2 * sqrt((float)num_balls*log2(num_bins)/num_bins));
    int server_bin_size = 526;

    Params params(server_bin_size, effective_bitLength, hw);
    
    for(int i = 0; i < 2; i++){
        std::cout << std::endl << "=========== " << i+1 << " times ===========" << std::endl;
 
        auto [context_string, keypack_string, data_stream, sk_string] = client.client_preprocessing(params);
        
        std::cout << context_string << ", " << keypack_string << std::endl;

        Server server(context_string, keypack_string, std::move(data_stream), sk_string);
        // std::stringstream final_stream = std::move(server.server_computation(params));
        // std::stringstream final_stream = std::move(server.serverMultipleLabelComp(params));
        std::stringstream final_stream = std::move(server.server_computation_time(params));

        client.load_decrypt(final_stream);
    }

    return 0;
}
