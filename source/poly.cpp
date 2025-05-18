#include <iostream>
#include <vector>
#include <cmath>

#include "HEaaN/HEaaN.hpp"
// #include "HEaaN-math/HEaaN-math.hpp"

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

    std::vector<Real> chebyCoefDeg3 = {
	5.875, -10.5, 6.0, -1.5, 
	0.125};
    u64 num_baby_step;

    if(chebyCoefDeg3.size() <= 3) num_baby_step = chebyCoefDeg3.size();
    else num_baby_step = pow(2,floor(ceil(log2(chebyCoefDeg3.size()))/2));

    ChebyshevCoefficients chebyshev_coef_deg_3(chebyCoefDeg3, num_baby_step);

    std::cout << "num_baby_step : " << num_baby_step << std::endl;
    std::cout << "num_bs : " << chebyshev_coef_deg_3.num_bs << std::endl;
    std::cout << "num_gs : " << chebyshev_coef_deg_3.num_gs << std::endl;
    std::cout << "log_gs : " << chebyshev_coef_deg_3.log_gs << std::endl;

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
    msg[0].real(0.0); msg[1].real(1.0);
    msg[2].real(2.0); msg[3].real(3.0);
    msg[4].real(4.0); msg[5].real(5.0);
    printMessage(msg);

    Ciphertext ctxt(context), ctxt_out(context);
    enc.encrypt(msg, pack, ctxt);

    std::cout << "Level before computing : " << ctxt.getLevel() << std::endl;

    eval.add(evaluateChebyshev(eval, ctxt, 
    chebyshev_coef_deg_3, 1), 0, ctxt);

    std::cout << "done ..." << std::endl;
    std::cout << "Level after computing : " << ctxt.getLevel() << std::endl;

    HEaaN::Message dmsg;
    std::cout << "Decrypt ... ";
    dec.decrypt(ctxt, sk, dmsg);
    std::cout << "done" << std::endl;

    std::cout << std::endl << "Result vector : " << std::endl;
    printMessage(dmsg);

    // evalCheby 할 때 재귀 부분에서 필요 이상으로 많은 연산 수행 중 
    // 이거 줄여야 한다. (완료) 

    // sqrt test

    std::cout << std::endl << "Sqrt test" << std::endl;

    for (size_t i = 0; i < num_slots; ++i) {
        msg[i].real(0.0);
        msg[i].imag(0.0);
    }
    msg[0].real(1.0); 
    msg[1].real(4.0);
    msg[2].real(6.25); 
    printMessage(msg);

    enc.encrypt(msg, pack, ctxt);
    std::cout << "Level before mult : " << ctxt.getLevel() << std::endl;

    eval.mult(ctxt, 1.0/9, ctxt);

    std::cout << "Level before computing : " << ctxt.getLevel() << std::endl;

    approxSqrtWilkes(eval, ctxt, ctxt_out, 3);

    std::cout << "Level before mult : " << ctxt_out.getLevel() << std::endl;
    eval.mult(ctxt_out, 3, ctxt_out);

    std::cout << "Level after computing : " << ctxt_out.getLevel() << std::endl;

    std::cout << "Decrypt ... ";
    dec.decrypt(ctxt_out, sk, dmsg);
    std::cout << "done" << std::endl;

    std::cout << std::endl << "Result vector : " << std::endl;
    printMessage(dmsg);

    return 0;
}