#include "client.hpp"

using namespace std;
using namespace HEaaN;

ParameterPreset Client::get_preset_from_string(const string& name) {
    static const std::unordered_map<string, ParameterPreset> preset_map = {
        {"FVa", ParameterPreset::FVa},
        {"FVb", ParameterPreset::FVb},
        {"FGa", ParameterPreset::FGa},
        {"FGb", ParameterPreset::FGb},
        {"FTa", ParameterPreset::FTa},
        {"FTb", ParameterPreset::FTb},
        {"ST19", ParameterPreset::ST19},
        {"ST14", ParameterPreset::ST14},
        {"ST11", ParameterPreset::ST11},
        {"ST8", ParameterPreset::ST8},
        {"ST7", ParameterPreset::ST7},
        {"SS7", ParameterPreset::SS7},
        {"SD3", ParameterPreset::SD3},
        {"CUSTOM", ParameterPreset::CUSTOM},
        {"FX", ParameterPreset::FX},
        {"FVc", ParameterPreset::FVc},
        {"FGd", ParameterPreset::FGd},
        {"SGd0", ParameterPreset::SGd0},
    };

    auto it = preset_map.find(name);
    if (it != preset_map.end()) {
        return it->second;
    }
    throw std::invalid_argument("Unknown preset name: " + name);
}

//Constructor
Client::Client(const string &name)
    : preset(get_preset_from_string(name)), context(makeContext(preset)), sk(context), pack(context),
    keygen(context,sk,pack), encryptor(context)
    // endecoder(context)
    {}

// save context and keypack and send over the server
void Client::save_parameters() {
    context_string = "context.bin";
    std::cout << "Trying to open file: [" << context_string << "]" << std::endl;
    saveContextToFile(context, context_string);
    std::ifstream f(context_string);
    if (!f.is_open()) {
            std::cerr << "File not found or not readable: " << context_string << std::endl;
    }
    
    keypack_string = "keypack";
    std::cout << "Trying to open file: [" << keypack_string << "]" << std::endl;
    pack.save(keypack_string);   
}

int Client::client_hash(int bin_index, int bitLength){
    return (bin_index+1) % (1 << bitLength);
}

void Client::save_ciphertext(vector <Ciphertext> input){
    data_stream.str("");
    for (int i = 0; i < input.size(); i++){
        input[i].save(data_stream);
    }
}

tuple<string, string, stringstream, string> Client::client_preprocessing(Params &params){
    // set up timer for experimentation
    cout << "Client preprocessing begins" << endl;
    int server_bin_size = params.server_bin_size, effective_bitLength = params.effective_bitLength, 
    ell = params.ell, hw = params.hw;

    // generate procedure keys
    cout << "Generating encryption, multiplication, and rotation keys ..." << endl;
    keygen.genEncKey();
    keygen.genMultKey();
    keygen.genRotKeyBundle();
    cout << "Done." << endl;

    const auto log_slots = getLogFullSlots(context);

    vector<Ciphertext> input;
    vector<vector<int>> encoding_vector(ell,vector<int>(1<<log_slots, 0));

    for (int bin_index = 0; bin_index < 1<<log_slots; bin_index++) {
        int hashed_val = client_hash(bin_index, effective_bitLength);
        vector<int> bit_val = PerfectMapping(hashed_val, ell, hw);
        for (int bit_index = 0; bit_index < ell; bit_index++) {
            encoding_vector[bit_index][bin_index] = bit_val[bit_index];
        }
    }

    // copy over to message object and encrypt
    for (size_t l = 0; l < ell; l++) {
        HEaaN::Complex complex{0};
        Message msg(log_slots, complex);
        for (size_t i = 0; i < msg.getSize(); ++i) {
            msg[i].real(encoding_vector[l][i]);
            if (i == msg.getSize()-1){
                std::cout << std::endl << "Message : " << msg.getSize() << std::endl;
                printMessage(msg);
                std::cout << std::endl;
            }
        }

        Ciphertext ctxt(context);
        std::cout << "Encrypt ... ";
        encryptor.encrypt(msg, pack, ctxt); 
        std::cout << "done" << std::endl;

        input.push_back(ctxt);
    } 
    
    cout << "Saving client ciphertexts ...";
    save_ciphertext(input);
    cout << "done" << endl;

    // for dubugging purposes ONLY.
    string sk_string = "secretkey.bin";
    sk.save(sk_string);

    save_parameters();
    return {context_string, keypack_string, move(data_stream), sk_string};
}

// Load and decrypt.
bool Client::load_decrypt(stringstream &final_stream) {
    Ciphertext ctxt(context);
    ctxt.load(final_stream);

    Decryptor decryptor(context);
    Message dmsg;
    decryptor.decrypt(ctxt, sk, dmsg);

    std::cout << std::endl << "(decrypted) Average : " << dmsg[0].real() << std::endl;
    std::cout << "Result Ciphertext - level " << ctxt.getLevel() << std::endl;

    return true;
}
