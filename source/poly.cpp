#include <iostream>
#include <vector>

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
    ParameterPreset preset = ParameterPreset::FGb;
    Context context = makeContext(preset);
    std::cout << "Parameter : " << presetNamer(preset) << std::endl
              << std::endl;

    const auto log_slots = getLogFullSlots(context);
    const auto num_slots = UINT64_C(1) << log_slots;

    // This is the approximation polynomial of degree 7
    // for sigmoid function 1 / (1 + exp(-x)).
    std::vector<Real> chebyshev_coef_deg_3 = {-1.5, 2.75, -1.5, 0.25};
    static const std::vector<Real> chebyshev_coef_deg_5{
	-28.75, 50.875, -30.0, 9.0625, 
	-1.25, 0.0625};


    SecretKey sk(context);
    KeyPack pack(context);
    KeyGenerator keygen(context, sk, pack);

    std::cout << "Generate encryption key ... ";
    keygen.genEncKey();
    std::cout << "done" << std::endl;

    std::cout << "Generate multiplication key ... ";
    keygen.genMultKey();
    std::cout << "done" << std::endl;

    Encryptor enc(context);
    Decryptor dec(context);
    HomEvaluator eval(context, pack);

    Message msg(log_slots);

    std::cout << std::endl << "Number of slots = " << num_slots << std::endl;

    for (size_t i = 0; i < num_slots; ++i) {
        msg[i].real(static_cast<HEaaN::Real>(i) /
                    static_cast<HEaaN::Real>(num_slots));
        msg[i].imag(0.0);
    }
    msg[0].real(2.0);
    printMessage(msg);

    Ciphertext ctxt(context);
    enc.encrypt(msg, pack, ctxt);

    std::cout << "Level before computing : " << ctxt.getLevel() << std::endl;

    eval.add(evaluateChebyshevExpansion(eval, ctxt, 
    chebyshev_coef_deg_5, false, 1), 0, ctxt);


    std::cout << "done ..." << std::endl;
    std::cout << "Level after computing : " << ctxt.getLevel() << std::endl;

    HEaaN::Message dmsg;
    std::cout << "Decrypt ... ";
    dec.decrypt(ctxt, sk, dmsg);
    std::cout << "done" << std::endl;

    std::cout << std::endl << "Result vector : " << std::endl;
    printMessage(dmsg);

    return 0;
}