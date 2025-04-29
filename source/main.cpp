#include <iostream>
#include <random>

#include "HEaaN/HEaaN.hpp"

inline double randNum() {
    static std::default_random_engine gen{std::random_device()()};
    std::uniform_real_distribution<double> dist(-1.0L, 1.0L);
    return dist(gen);
}

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
    default:
        throw std::invalid_argument("Not supported parameter");
    }
}

int main() {
    std::cout << randNum() << std::endl;
    HEaaN::ParameterPreset preset = HEaaN::ParameterPreset::SS7;
    HEaaN::Context context = makeContext(preset);
    std::cout << "Parameter : " << presetNamer(preset) << std::endl
              << std::endl;
}