#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>

#include "HEaaN/HEaaN.hpp"

#include "utils.hpp"
#include "PolynomialEvaluator.hpp"

using namespace HEaaN;

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

int main(void){
    ParameterPreset preset = ParameterPreset::SS7;
    Context context = makeContext(preset);
    std::cout << "Parameter : " << presetNamer(preset) << std::endl
              << std::endl;

    const auto log_slots = getLogFullSlots(context);
    const auto num_slots = UINT64_C(1) << log_slots;

    SecretKey sk(context);
    KeyPack pack(context);
    KeyGenerator keygen(context, sk, pack);

    std::cout << "Generate encryption key ... ";
    keygen.genEncKey();
    std::cout << "done" << std::endl;

    std::cout << "Generate multiplication key ... ";
    keygen.genMultKey();

    // std::cout << "Generating BTS key ..." << std::endl;

    // keygen.genConjKey();
    // keygen.genRotKeysForBootstrap(log_slots);

    std::cout << "done" << std::endl;

    Timer_micro key_timer;
    long key_time;

    Encryptor enc(context);
    Decryptor dec(context);
    EnDecoder endecoder(context);

    HomEvaluator eval(context, pack);
    // Bootstrapper btp(eval);

    // std::cout << std::endl << "Number of slots = " << num_slots << std::endl;

    Message msg(log_slots), msg_zero(log_slots);
    for (size_t i = 0; i < num_slots; ++i) {
        msg_zero[i].real(0.0);
        msg_zero[i].imag(0.0);
    }

    std::cout << "Construct message" << std::endl;
    key_timer.start();

    for (int i = 0; i < num_slots; ++i) {
        msg[i].real(i);
        msg[i].imag(0.0);
    }

    key_time = key_timer.end_and_get();
    std::cout << "done ..." << std::endl;
    std::cout << key_time << " usec" << std::endl;
    // printMessage(msg);
    std::cout << std::endl;


    std::cout << "Encrypt msg (make ctxt)" << std::endl;
    Ciphertext ctxt(context);
    key_timer.start();
    
    enc.encrypt(msg, pack, ctxt);
    
    key_time = key_timer.end_and_get();
    std::cout << "done ..." << std::endl;
    std::cout << key_time << " usec" << std::endl << std::endl;


    std::cout << "Encode msg (make ptxt)" << std::endl;
    key_timer.start();
    Plaintext ptxt(context);

    ptxt = endecoder.encode(msg, 
                            ctxt.getLevel(), ctxt.getRescaleCounter());

    // ptxt = endecoder.encode(msg, 7, 0);

    key_time = key_timer.end_and_get();
    std::cout << "done ..." << std::endl;
    std::cout << key_time << " usec" << std::endl << std::endl;                


    // std::cout << "Check NTT for ptxt" << std::endl;
    // std::cout << ptxt.getMx().isNTT() << std::endl << std::endl;
        

    std::cout << "Encrypt ptxt (make ctxt)" << std::endl;
    key_timer.start();
    
    enc.encrypt(ptxt, pack, ctxt);
    
    key_time = key_timer.end_and_get();
    std::cout << "done ..." << std::endl;
    std::cout << key_time << " usec" << std::endl << std::endl;

    
    // std::cout << "Check NTT for ctxt" << std::endl;
    // std::cout << ctxt.getPoly(0).getMx().isNTT() << std::endl;




    Ciphertext test1(context), test2(context);
    enc.encrypt(msg_zero, pack, test1);
    enc.encrypt(msg_zero, pack, test2);


    std::cout << "Add ctxt and ctxt" << std::endl;
    key_timer.start();

    eval.add(ctxt, test1, test1);

    key_time = key_timer.end_and_get();
    std::cout << "done ..." << std::endl;
    std::cout << key_time << " usec" << std::endl << std::endl;      
     
    
    std::cout << "mult ptxt and ctxt" << std::endl;
    key_timer.start();

    eval.mult(ctxt, ptxt, test1);

    key_time = key_timer.end_and_get();
    std::cout << "done ..." << std::endl;
    std::cout << key_time << " usec" << std::endl << std::endl; 
    
    
    std::cout << "mult ctxt and ctxt" << std::endl;
    key_timer.start();

    eval.mult(ctxt, test2, test1);

    key_time = key_timer.end_and_get();
    std::cout << "done ..." << std::endl;
    std::cout << key_time << " usec" << std::endl << std::endl;     

    // ============================================

    Message msg1(log_slots), zero_msg(log_slots);
    fillReal(msg1, 0.0);

    std::cout << "Outside" << std::endl;
    key_timer.start();

    for(int i = 0; i < 100; i++){
        eval.add(msg1, 0, zero_msg);
    }

    key_time = key_timer.end_and_get();
    std::cout << "done ..." << std::endl;
    std::cout << key_time << " usec" << std::endl << std::endl;    

    std::cout << "Inside" << std::endl;
    key_timer.start();
    
    for (int i = 0; i < 100; i++){
        HEaaN::Complex complex{0};
        Message msg(log_slots, complex);
    }

    key_time = key_timer.end_and_get();
    std::cout << "done ..." << std::endl;
    std::cout << key_time << " usec" << std::endl << std::endl;    
    
    
    
    // std::cout << "Sum" << std::endl;
    // key_timer.start();

    // double n = 0;
    // for(int i = 0; i < msg.getSize() ; i++){
    //     n += msg[i].real();
    // }

    // key_time = key_timer.end_and_get();
    // std::cout << "done ..." << std::endl;
    // std::cout << key_time << " usec" << std::endl << std::endl;    


    std::cout << sizeof(Message) << std::endl;
    std::cout << sizeof(Plaintext) << std::endl;
    std::cout << sizeof(Ciphertext) << std::endl;

    return 0;
}