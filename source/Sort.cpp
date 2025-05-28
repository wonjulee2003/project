////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Copyright (C) 2021-2023 Crypto Lab Inc.                                    //
//                                                                            //
// - This file is part of HEaaN homomorphic encryption library.               //
// - HEaaN cannot be copied and/or distributed without the express permission //
//  of Crypto Lab Inc.                                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "Sort.hpp"

#include "ApproximationMethods.hpp"
#include "ProgressBar.hpp"

#include <algorithm>
#include <vector>

Message genMasks(u64 log_slots, u64 dist) {
    Message mask(log_slots);
    std::generate_n(mask.begin(), mask.getSize(), [dist, i = 0]() mutable {
        u64 j = i++ / dist;
        return (j & U64ONE) ? REAL_ZERO : REAL_ONE;
    });
    return mask;
}

void compForSort(const HomEvaluator &eval, const Bootstrapper &btp,
                 const Ciphertext &op, const Message &mask, int dist,
                 Ciphertext &res) {
    Ciphertext rot(eval.getContext());
    eval.leftRotate(op, dist, rot);
    if (rot.getLevel() < 1 + btp.getMinLevelForBootstrap()){
        btp.bootstrap(rot, rot);}
    eval.mult(rot, mask, rot);
    eval.sub(op, rot, res);
    approxSignDecomposition(eval, btp, res, res, 13, -18, 0.5);
    eval.add(res, 0.5, res);
}

void assembleSlot(const HomEvaluator &eval, std::vector<Ciphertext> &op_arr,
                  u64 dist, Ciphertext &res) {
    res = op_arr[0];
    for (u64 i = 1; i < op_arr.size(); i++) {
        eval.rightRotate(op_arr[i], i * dist, op_arr[i]);
        eval.add(res, op_arr[i], res);
    }
}

void flipCtxt(const HomEvaluator &eval, Ciphertext &res, int dist) {
    u64 log_slots = res.getLogSlots();
    Message mask(log_slots);
    std::generate_n(mask.begin(), mask.getSize(), [dist, i = 0]() mutable {
        u64 j = i++ / (2 * dist);
        return (j & U64ONE) ? REAL_ONE : -REAL_ONE;
    });

    mask.to(res.getDevice());
    eval.mult(res, mask, res);
}

///////////////////////////////////////////////////////////////////////////////
/////////////////////////////// two Sorter ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void twoSubSort(const HomEvaluator &eval, const Ciphertext &op1,
                const Ciphertext &op2, const Ciphertext &op_comp, bool ascent,
                std::vector<Ciphertext> &maxmin_arr) {
    if (ascent) {
        eval.sub(op1, op2, maxmin_arr[1]);
        eval.mult(op_comp, maxmin_arr[1], maxmin_arr[1]);
        eval.add(maxmin_arr[1], op2, maxmin_arr[1]);
        eval.add(op1, op2, maxmin_arr[0]);
        eval.sub(maxmin_arr[0], maxmin_arr[1], maxmin_arr[0]);
    } else {
        eval.sub(op1, op2, maxmin_arr[0]);
        eval.mult(op_comp, maxmin_arr[0], maxmin_arr[0]);
        eval.add(maxmin_arr[0], op2, maxmin_arr[0]);
        eval.add(op1, op2, maxmin_arr[1]);
        eval.sub(maxmin_arr[1], maxmin_arr[0], maxmin_arr[1]);
    }
}

void twoSort(const HomEvaluator &eval, const Bootstrapper &btp, Ciphertext &res,
             int dist1, bool ascent) {
    u64 log_slots = res.getLogSlots();
    Message mask = genMasks(log_slots, dist1);
    mask.to(res.getDevice());

    if (res.getLevel() < 3 + btp.getMinLevelForBootstrap()){
        btp.bootstrap(res, res);}

    Ciphertext ctxt_comp(eval.getContext());
    compForSort(eval, btp, res, mask, dist1, ctxt_comp);

    if (ctxt_comp.getLevel() < 1 + btp.getMinLevelForBootstrap()){
        btp.bootstrap(ctxt_comp, ctxt_comp);}
    Ciphertext ctxt1(eval.getContext()), ctxt2(eval.getContext());
    eval.mult(res, mask, ctxt1);
    eval.leftRotate(res, dist1, res);
    eval.mult(res, mask, ctxt2);

    std::vector<Ciphertext> ctxt_sort_arr;
    ctxt_sort_arr.reserve(2);
    ctxt_sort_arr.emplace_back(eval.getContext());
    ctxt_sort_arr.emplace_back(eval.getContext());
    twoSubSort(eval, ctxt1, ctxt2, ctxt_comp, ascent, ctxt_sort_arr);
    assembleSlot(eval, ctxt_sort_arr, dist1, res);
}

///////////////////////////////////////////////////////////////////////////////
/////////////////////////////// main Sorter ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void unitSort(const HomEvaluator &eval, const Bootstrapper &btp,
              const Ciphertext &op, Ciphertext &res, int n, bool ascent,
              bool only_last_stage) {
    res = op;

    int start_idx = 0;
    int log_n = std::ceil(std::log2(n));
    if (only_last_stage)
        start_idx = log_n - 1;

    double count = 0.0;
    const double max_count_inv = 100.0 * 2 / ((log_n + 1) * log_n);
    for (int i = start_idx; i < log_n; i++) {
        int j = i;
        bool flip = !(i == log_n - 1);

        if (flip) {
            if (res.getLevel() < 1 + btp.getMinLevelForBootstrap()){
                btp.bootstrap(res, res);}
            flipCtxt(eval, res, 1 << i);
        }

        while (j >= 0) {
            twoSort(eval, btp, res, 1 << j, ascent);
            count += 1;
            ProgressBar::printBar(count * max_count_inv, "unitSort");
            j--;
        }

        if (flip) {
            if (res.getLevel() < 1 + btp.getMinLevelForBootstrap()){
                btp.bootstrap(res, res);}
            flipCtxt(eval, res, 1 << i);
        }
    }
    ProgressBar::printBarEnd();
}

