#ifndef __CLIENT__H
#define __CLIENT__H

#include <iostream>
#include <random>
#include <vector>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <stdexcept>
#include <tuple>
#include <fstream>
#include <string>

#include "HEaaN/HEaaN.hpp"
#include "utils.hpp"

class Client {
public:
    // HEaaN tools
    HEaaN::ParameterPreset preset;
    HEaaN::Context context;
    HEaaN::SecretKey sk;
    HEaaN::KeyPack pack;
    HEaaN::KeyGenerator keygen;
    HEaaN::Encryptor encryptor;
    // HEaaN::EnDecoder endecoder;

    // file transfer
    std::string context_string;
    std::string keypack_string;
    std::stringstream data_stream;
    
    // methods
    Client(const string &name);
    HEaaN::ParameterPreset get_preset_from_string(const string& name);
    void save_parameters();
    int client_hash(int bin_index, int bitLength);
    void save_ciphertext(std::vector <HEaaN::Ciphertext> input);
    std::tuple<std::string, std::string, std::stringstream, std::string> client_preprocessing(Params &params);
    bool load_decrypt(stringstream &final_stream);
};


#endif