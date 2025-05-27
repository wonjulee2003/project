////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Copyright (C) 2021-2023 Crypto Lab Inc.                                    //
//                                                                            //
// - This file is part of HEaaN homomorphic encryption library.               //
// - HEaaN cannot be copied and/or distributed without the express permission //
//  of Crypto Lab Inc.                                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <iostream>

class ProgressBar {
public:
    static ProgressBar *getInstance() {
        static ProgressBar instance;
        return &instance;
    }

    inline static bool getSwitch() { return getInstance()->onoff_; }
    inline static void setSwitch(const bool onoff) {
        getInstance()->onoff_ = onoff;
    }

    inline static void printBar(double percent, const char *funct_name) {
        constexpr char BAR = '=';
        constexpr char BAR_END = '>';
        constexpr char BLANK = ' ';
        constexpr int LEN = 40;
        constexpr double TICK = 100.0 / LEN;

        if (getSwitch()) {
            printf("\r%s [", funct_name);
            const int bar_count = static_cast<int>(percent / TICK);
            for (int i = 0; i < LEN; i++) {
                if (bar_count > i) {
                    if (bar_count <= i + 1) {
                        printf("%c", BAR_END);
                    } else {
                        printf("%c", BAR);
                    }
                } else {
                    printf("%c", BLANK);
                }
            }
            printf("] %0.2lf%%", percent);
            fflush(stdout);
        }
    }

    inline static void printBarEnd() {
        if (getSwitch()) {
            printf("\n\n");
        }
    }

private:
    ProgressBar() = default;
    ~ProgressBar() = default;

    bool onoff_ = false;
};
