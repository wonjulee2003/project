#include "server.hpp"

using namespace std;
using namespace HEaaN;
// using namespace HEaaN::Math;

int Server::server_hash(int mu, int bin_index, int bitLength){
    return ((bin_index+1)+mu) % (1 << bitLength);
}

// constructor
Server::Server(const string &string1, string &string2, stringstream &&source_stream, string &string4)
    : // context(makeContextFromFile(string1)), pack(context, string2), evaluator(context, pack), encryptor(context), btp(evaluator),
       context(makeContextFromFile(string1)), pack(context, string2), evaluator(context, pack), encryptor(context), endecoder(context),  
    //   context(makeContextFromFile(string1)), pack(context, string2), evaluator(context, pack), encryptor(context),
    // sk passed over for debugging purposes
    sk(context, string4),
    context_string{string1}, keypack_string{string2}, data_stream{source_stream.str()}, sk_string{string4}
{
    std::cout << getDegree(context) << std::endl;
}

void Server::compute_sum(Ciphertext &ctxt) {
    const auto log_slots = getLogFullSlots(context);
    size_t rotation_required = std::ceil(log_slots);
        HEaaN::Ciphertext ctxt_rotated(context);
        for (size_t i = 0; i < rotation_required; ++i) {
            evaluator.leftRotate(ctxt, UINT64_C(1) << i, ctxt_rotated);
            evaluator.add(ctxt, ctxt_rotated, ctxt);
        }
}

void Server::send_result(Ciphertext &ctxt) {
    data_stream.str("");
    ctxt.save(data_stream);
}


stringstream Server::server_computation(Params &params){
    cout << "Server computation begins" << endl;

    int server_bin_size = params.server_bin_size, effective_bitLength = params.effective_bitLength, 
    ell = params.ell, hw = params.hw;
    
    const auto log_slots = getLogFullSlots(context);

    vector<Ciphertext> inputs;

    //loading ciphertexts
    Ciphertext ct_temp(context);
    for (int i=0;i<params.ell;i++){
        ct_temp.load(data_stream);
        inputs.push_back(ct_temp);
    }
    data_stream.str("");
    std::cout << "Saved number of ciphertexts : " << inputs.size() << std::endl;


    // create labels for Circuit + Labelled PSI
    // using single ciphertext for simulation
    Message label_msg(log_slots);
    double val = 5;
    fillReal(label_msg, val);

    // subtraction vector of [hw]
    vector<Message> subtraction(hw-1);
    for (double i=1;i<hw;i++){
        Message msg(log_slots);
        fillReal(msg, i);
        subtraction[i-1] = msg;
    }
    vector<Message> factorial(hw-1);
    for (double i=2;i<hw+1;i++) {
        Message msg(log_slots);
        fillReal(msg, 1.0/i);
        factorial[i-2] = msg;
    }

    // store final results here
    vector<Ciphertext>mu_result;
    vector<Ciphertext>intersection;

    Message zero_msg(log_slots); 
    fillReal(zero_msg, 0.0);

    vector<Ciphertext> intersection_cheby;

    const std::vector<Real> chebyCoefDeg3 = ChebyCoefTable::ChebyCoefDeg[3];
    u64 num_baby_step;

    if(chebyCoefDeg3.size() <= 3) num_baby_step = chebyCoefDeg3.size();
    else num_baby_step = pow(2,floor(ceil(log2(chebyCoefDeg3.size()))/2));

    // beginning of main loop
    for (int i=0; i<server_bin_size; i++) {
        vector<vector<int>> operand_vec(ell,vector<int>(1<<log_slots, 0));

        for (int bin_index = 0; bin_index < 1<<log_slots; bin_index++) {
            int hashed_val = server_hash(i, bin_index, effective_bitLength);
            vector<int> bit_val = PerfectMapping(hashed_val, ell, hw);
            for (int bit_index = 0; bit_index < ell; bit_index++) {
                operand_vec[bit_index][bin_index] = bit_val[bit_index];
            }
        } 
        
        Ciphertext ctxt(context);
        // Message zero_msg(log_slots); 
        fillReal(zero_msg, 0.0);
        encryptor.encrypt(zero_msg, pack, ctxt);
        

        for (int l=0; l<ell; l++){
            HEaaN::Complex complex{0};
            Message msg(log_slots, complex);
            Ciphertext temp(context);
            for (int k=0; k<msg.getSize(); k++) {
                msg[k].real(operand_vec[l][k]);
                // if (k == msg.getSize()-1 && i == 4){
                //     std::cout << std::endl << "Server Message : " << std::endl;
                //     printMessage(msg);
                //     std::cout << std::endl;
                // }
            }

            evaluator.mult(inputs[l], msg, temp);
            // std::cout << "Result Ciphertext - level " << temp.getLevel() << std::endl;

            // sum all individual ciphertext vectors in temp_ct to obtain h' (ctxt)
            evaluator.add(ctxt, temp, ctxt);
        }

        // add Chebyshev Basis Evaluation Function HERE
        // input : ciphertext ctxt, poly coefficients (potentially using a lookup table for small domains)
        // output : ciphertext ctxt (or other ciphertext storage)

        Ciphertext ctxt_test(context);
        ChebyshevCoefficients chebyshev_coef_deg_3(chebyCoefDeg3, num_baby_step);

        if (i==0){
            std::cout << "Ciphertext - initial level " << ctxt.getLevel() << std::endl;
        }
        
        evaluator.add(evaluateChebyshev(evaluator, ctxt, 
            chebyshev_coef_deg_3, 1.0/6), 0, ctxt_test);

        if (i==0){
            std::cout << "Chebyshev Ciphertext - post-comp level " << ctxt_test.getLevel() << std::endl;
        }

        intersection_cheby.push_back(ctxt_test);

        // Alternative: Direct Product Evaluation
        Ciphertext temp(context);
        Ciphertext result(context);

        evaluator.mult(ctxt, 1.0, result);

        if (i==0){
            std::cout << "Result Ciphertext - initial level " << result.getLevel() << std::endl;
        }

        for (int k=0;k<hw-1;k++){
            evaluator.sub(ctxt, subtraction[k], temp);
            evaluator.mult(result, temp, result);
            evaluator.mult(result, factorial[k], result);
        }
        // store final value in mu_result[i] and multiply by the label
        intersection.push_back(result);
        evaluator.mult(result, label_msg, result);

        if (i==0){
            std::cout << "Result Ciphertext - post-comp level " << result.getLevel() << std::endl;
        }

        mu_result.push_back(result);


        // bootstrapping? 
        // Using FGb parameters, we reach level 6 for result ciphertexts using hamming weight of 3. 
        // level 2 for hamming weight of 5
    }

    // ================================

    // sum each index of mu_result[i]
    Ciphertext size(context);
    // Message zero_msg(log_slots); 
    fillReal(zero_msg, 0.0);
    encryptor.encrypt(zero_msg, pack, size);

    for (int i = 0; i < mu_result.size(); i++){
        evaluator.add(size, intersection[i], size);
    }

    Ciphertext final(context);
    encryptor.encrypt(zero_msg, pack, final);

    for (int i = 0; i < mu_result.size(); i++){
        evaluator.add(final, mu_result[i], final);
    }
    std::cout << "Result Ciphertext - level " << final.getLevel() << std::endl;

    // ================================

    // compute the size of intersection
    std::cout << "Compute intersection size ... ";
    compute_sum(size);
    std::cout << "done" << std::endl;

    // DEBUGGING
    Decryptor decryptor(context);
    Message dmsg(log_slots);
    decryptor.decrypt(size, sk, dmsg);
    printMessage(dmsg);

    // compute average
    std::cout << "Compute average ... ";
    compute_sum(final);
    evaluator.mult(final, 1.0 / (1 << log_slots), final);
    std::cout << "done" << std::endl;

    
    // compute standard deviation on vector (create a function)


    // // encrypt and send over (create a function)
    send_result(final);
    cout << "Server side completed." << endl;

    return move(data_stream);
}

stringstream Server::server_computation_time(Params &params){
    cout << "Server computation <time> begins" << endl;

    int server_bin_size = params.server_bin_size, effective_bitLength = params.effective_bitLength, 
    ell = params.ell, hw = params.hw;

    const auto log_slots = getLogFullSlots(context);

    vector<Ciphertext> inputs;

    std::cout << "Loading Client ctxt(communication cost)" << std::endl;

    std::ofstream file;
    file.open("/home/jimineum/ckks_psi/project/load_ctxt.txt", ios::app);

    if(!file.is_open()){
        std::cout << "Cannot open the file" << std::endl;
        exit(1);
    }

    Timer key_timer;
    key_timer.start();
    
    //loading ciphertexts
    Ciphertext ct_temp(context);
    for (int i=0;i<params.ell;i++){
        ct_temp.load(data_stream);
        inputs.push_back(ct_temp);
    }
    long key_time = key_timer.end_and_get();
    std::cout << key_time << " ms" << std::endl;

    file << key_time << std::endl;
    file.close();



    data_stream.str("");
    std::cout << "Saved number of ciphertexts : " << inputs.size() << std::endl;

    // create labels for Circuit + Labelled PSI
    // using single ciphertext for simulation
    // Message label_msg(log_slots);
    // double val = 5;
    // fillReal(label_msg, val);

    // find fac for Chebyshev
    Real fac = 1;
    for(Real i = 2; i <= hw; i++) fac *= i;

    Ciphertext init(context);
    // Message zero_msg(log_slots); 
    // fillReal(zero_msg, 0.0);

    vector<Ciphertext> intersection_cheby(server_bin_size, init);

    // const std::vector<Real> chebyCoefDeg = ChebyCoefTable::ChebyCoefDeg[hw];
    // u64 num_baby_step;

    // if(hw+1 <= 3) num_baby_step = hw+1;
    // else num_baby_step = pow(2,floor(ceil(log2(hw+1))/2));


    // std::cout << "Server data preprocessing begins" << std::endl;

    // file.open("/home/jimineum/ckks_psi/project/server_dataprocessing.txt", ios::app);

    // if(!file.is_open()){
    //     std::cout << "Cannot open the file" << std::endl;
    //     exit(1);
    // }

    // key_timer.start();

    // vector<vector<vector<int>>> operand_vec(server_bin_size, 
    //     vector<vector<int>>(ell,vector<int>(1<<log_slots, 0)));

    // for (int i=0; i<server_bin_size; i++) {
    //     for (int bin_index = 0; bin_index < 1<<log_slots; bin_index++) {
    //         int hashed_val = server_hash(i, bin_index, effective_bitLength);
    //         vector<int> bit_val = PerfectMapping(hashed_val, ell, hw);
    //         for (int bit_index = 0; bit_index < ell; bit_index++) {
    //             operand_vec[i][bit_index][bin_index] = bit_val[bit_index];
    //         }
    //     } 
    // }

    // key_time = key_timer.end_and_get();
    // std::cout << key_time << " ms" << std::endl;

    // file << key_time << std::endl;
    // file.close();

    int gamma = 1;

    std::cout << "Server Computation begins" << std::endl;

    file.open("/home/jimineum/ckks_psi/project/server_comp.txt", ios::app);

    if(!file.is_open()){
        std::cout << "Cannot open the file" << std::endl;
        exit(1);
    }

    key_timer.start();

        // beginning of main loop

        #pragma omp parallel for
        for (int i=0; i<server_bin_size; i++) {
            vector<vector<int>> operand_vec(ell,vector<int>(1<<log_slots, 0));

            // std::cout << "Server preprocessing" << std::endl;
            Timer key_timer1;
            // key_timer1.start();

            for (int bin_index = 0; bin_index < 1<<log_slots; bin_index++) {
                int hashed_val = server_hash(i, bin_index, effective_bitLength);
                vector<int> bit_val = PerfectMapping(hashed_val, ell, hw);
                for (int bit_index = 0; bit_index < ell; bit_index++) {
                    operand_vec[bit_index][bin_index] = bit_val[bit_index];
                }
            } 
            
            // key_time = key_timer1.end_and_get();
            // std::cout << key_time << " ms" << std::endl << std::endl;

            // std::cout << "Multiply and Add" << std::endl;
            // key_timer1.start();

            Ciphertext ctxt(context);
            Message zero_msg(log_slots); 
            fillReal(zero_msg, 0.0);
            encryptor.encrypt(zero_msg, pack, ctxt);
            

            for (int l=0; l<ell; l++){
                // Timer tempTime;
                // long temp_time;

                // if(l == 0){
                //     std::cout << "running time for constructing message" << std::endl;
                //     tempTime.start();
                // }

                HEaaN::Complex complex{0};
                Message msg(log_slots, complex);
                Ciphertext temp(context);
                for (int k=0; k<msg.getSize(); k++) {
                    // msg[k].real(operand_vec[i][l][k]); // separate

                    msg[k].real(operand_vec[l][k]); // combine
                    // if (k == msg.getSize()-1 && i == 4){
                    //     std::cout << std::endl << "Server Message : " << std::endl;
                    //     printMessage(msg);
                    //     std::cout << std::endl;
                    // }
                }

                // if(l == 0){
                //     temp_time = tempTime.end_and_get();
                //     std::cout << temp_time << " ms" << std::endl;
                // }

                // tempTime.start();

                evaluator.mult(inputs[l], msg, temp);

                // temp_time = tempTime.end_and_get();
                // std::cout << temp_time << " ms" << std::endl;

                // sum all individual ciphertext vectors in temp_ct to obtain h' (ctxt)
                
                // tempTime.start();

                evaluator.add(ctxt, temp, ctxt);

                // temp_time = tempTime.end_and_get();
                // std::cout << temp_time << " ms" << std::endl;
            }

            // key_time = key_timer1.end_and_get();
            // std::cout << key_time << " ms" << std::endl << std::endl;

            // add Chebyshev Basis Evaluation Function HERE
            // input : ciphertext ctxt, poly coefficients (potentially using a lookup table for small domains)
            // output : ciphertext ctxt (or other ciphertext storage)


            std::cout << "Poly Eval" << std::endl;
            key_timer1.start();

            Ciphertext ctxt_test(context);
            const std::vector<Real> chebyCoefDeg = ChebyCoefTable::ChebyCoefDeg[hw];
            u64 num_baby_step;

            if(hw+1 <= 3) num_baby_step = hw+1;
            else num_baby_step = pow(2,floor(ceil(log2(hw+1))/2));

            ChebyshevCoefficients chebyshev_coef_deg(chebyCoefDeg, num_baby_step);

            // if (i==0){
            //     std::cout << "Ciphertext - initial level " << ctxt.getLevel() << std::endl;
            // }
            
            evaluator.add(evaluateChebyshev(evaluator, ctxt, 
                chebyshev_coef_deg, 1.0/fac), 0, ctxt_test);

            // if (i==0){
            //     std::cout << "Chebyshev Ciphertext - post-comp level " << ctxt_test.getLevel() << std::endl;
            // }

            evaluator.add(ctxt_test,0,intersection_cheby[i]);

            key_time = key_timer1.end_and_get();
            std::cout << key_time << " ms" << std::endl << std::endl;


            // Alternative: Multiply like Tree

            // std::cout << "Poly Eval" << std::endl;
            // key_timer1.start();

            // vector<Ciphertext> subtraction(hw, init);
            // for(int i = 0; i < hw; i++){
            //     Ciphertext temp(context);
            //     evaluator.sub(ctxt, i, temp);
            //     evaluator.mult(temp, 1.0/(hw-i), temp);
            //     subtraction[i] = temp;
            // }

            // Ciphertext result(context);
            // result = multiplyLikeTree(evaluator, subtraction, 0, hw-1);
            
            // // divide by 1/hw!
            // evaluator.add(result, 0, intersection_cheby[i]);

            // key_time = key_timer1.end_and_get();
            // std::cout << key_time << " ms" << std::endl << std::endl;

            // bootstrapping? 
            // Using FGb parameters, we reach level 6 for result ciphertexts using hamming weight of 3. 
            // level 2 for hamming weight of 5
        }
        // ================================

        // // sum each index of mu_result[i]
        // Ciphertext size(context);
        // // Message zero_msg(log_slots); 
        // fillReal(zero_msg, 0.0);
        // encryptor.encrypt(zero_msg, pack, size);

        // for (int i = 0; i < mu_result.size(); i++){
        //     evaluator.add(size, intersection[i], size);
        // }

        // Ciphertext final(context);
        // encryptor.encrypt(zero_msg, pack, final);

        // for (int i = 0; i < mu_result.size(); i++){
        //     evaluator.add(final, mu_result[i], final);
        // }
        // std::cout << "Result Ciphertext - level " << final.getLevel() << std::endl;

        // ================================

        // Check the case that only find intersection

        // std::cout << "Combine" << std::endl;
        // Timer key_timer1;
        // key_timer1.start();

        Ciphertext size(context);
        Message zero_msg(log_slots); 
        fillReal(zero_msg, 0.0);
        encryptor.encrypt(zero_msg, pack, size);

        // #pragma omp parallel for
        for (int i = 0; i < intersection_cheby.size(); i++){
            // #pragma omp critical

            evaluator.add(size, intersection_cheby[i], size);
        }

        // key_time = key_timer1.end_and_get();
        // std::cout << key_time << " ms" << std::endl << std::endl;

    key_time = key_timer.end_and_get();
    std::cout << key_time << " ms" << std::endl;
    file << key_time << std::endl;
    file.close();


    // compute the size of intersection
    // std::cout << "Compute intersection size ... ";
    // compute_sum(size);
    // std::cout << "done" << std::endl;

    // DEBUGGING
    Decryptor decryptor(context);
    Message dmsg(log_slots);
    decryptor.decrypt(size, sk, dmsg);
    printMessage(dmsg);

    // compute average
    // std::cout << "Compute average ... ";
    // compute_sum(final);
    // evaluator.mult(final, 1.0 / (1 << log_slots), final);
    // std::cout << "done" << std::endl;

    
    // compute standard deviation on vector (create a function)


    // // encrypt and send over (create a function)
    // send_result(final);
    send_result(size);
    cout << "Server side completed." << endl;

    return move(data_stream);
}


stringstream Server::server_computation_time0(Params &params){
    cout << "Server computation <time> begins" << endl;

    int server_bin_size = params.server_bin_size, effective_bitLength = params.effective_bitLength, 
    ell = params.ell, hw = params.hw;

    const auto log_slots = getLogFullSlots(context);

    Timer_micro key_timer;
    long key_time;

    std::cout << "Loading Client ctxt(communication cost)" << std::endl;
    key_timer.start();

    // std::ofstream file;
    // file.open("/home/jimineum/ckks_psi/project/load_ctxt.txt", ios::app);

    // if(!file.is_open()){
    //     std::cout << "Cannot open the file" << std::endl;
    //     exit(1);
    // }

    
    //loading ciphertexts

    vector<Ciphertext> inputs;
    Ciphertext ct_temp(context);

    for (int i=0;i<params.ell;i++){
        ct_temp.load(data_stream);
        
        // down level for efficient computing
        evaluator.levelDown(ct_temp, ct_temp.getLevel()-2, ct_temp);

        inputs.push_back(ct_temp);
    }

    // ct_temp level = 5
    // std::cout << "ct_temp level : " << ct_temp.getLevel() << std::endl;

    key_time = key_timer.end_and_get();
    std::cout <<  "done ..." << std::endl;
    std::cout << key_time << " usec" << std::endl;

    // file << key_time << std::endl;
    // file.close();


    data_stream.str("");
    std::cout << "Saved number of ciphertexts : " << inputs.size() << std::endl;


    std::cout << std::endl << "Server preprocessing" << std::endl;
    key_timer.start();

    // init_ct level = 4
    Ciphertext init_ct(context);
    evaluator.levelDown(init_ct, init_ct.getLevel()-3, init_ct);

    Message zero_msg(log_slots); fillReal(zero_msg, 0.0);

    // init_pt level = 5
    Plaintext init_pt = endecoder.encode(zero_msg, ct_temp.getLevel(), 0);
    encryptor.encrypt(init_pt, pack, init_ct);    

    // Ciphertext ctxt(context), temp(context);
    // evaluator.levelDown(ctxt, init_ct.getLevel(), ctxt);
    // evaluator.levelDown(temp, init_ct.getLevel()-1, temp);

    // Message zero_msg(log_slots); 
    // fillReal(zero_msg, 0.0);

    int gamma = 1;
    vector<Ciphertext> intersection(server_bin_size, init_ct);
    
    std::cout << sizeof(intersection) << std::endl;

    // Server preprocessing

    std::cout << "size : " << server_bin_size << std::endl;

    vector<vector<Plaintext>> encode_message(server_bin_size, 
                                    vector<Plaintext>(ell, init_pt));

    std::cout << "begin" << std::endl;

    #pragma omp parallel for
    for (int i=0; i<server_bin_size; i++) {
        vector<vector<int>> operand_vec(ell,vector<int>(1<<log_slots, 0));

        // #pragma omp parallel for
        for (int bin_index = 0; bin_index < 1<<log_slots; bin_index++) {
            int hashed_val = server_hash(i, bin_index, effective_bitLength);
            vector<int> bit_val = PerfectMapping(hashed_val, ell, hw);

            for (int bit_index = 0; bit_index < ell; bit_index++) {
                operand_vec[bit_index][bin_index] = bit_val[bit_index];
            }
        } 

        HEaaN::Complex complex{0};
        Message msg(log_slots, complex);

        for (int l=0; l<ell; l++){

            // #pragma omp parallel for
            int n = 0;
            for (int k=0; k<msg.getSize(); k++) {
                int num = operand_vec[l][k];

                msg[k].real(num);
                n += num; 
            }

            if (n != 0){
                evaluator.add(
                    endecoder.encode(msg, init_ct.getLevel(), init_ct.getRescaleCounter()),
                    0, encode_message[i][l]
                );
            }
        }
    
    }

    key_time = key_timer.end_and_get();
    std::cout <<  "done ..." << std::endl;
    std::cout << key_time << " usec" << std::endl;


    std::cout << std::endl << "Server Computation begins" << std::endl;

    // file.open("/home/jimineum/ckks_psi/project/server_comp.txt", ios::app);

    // if(!file.is_open()){
    //     std::cout << "Cannot open the file" << std::endl;
    //     exit(1);
    // }

    key_timer.start();

        // beginning of main loop

        #pragma omp parallel for
        for (int i=0; i<server_bin_size; i++) {
            // vector<vector<int>> operand_vec(ell,vector<int>(1<<log_slots, 0));

            // std::cout << "Server preprocessing" << std::endl;
            // Timer_micro key_timer1;
            // key_timer1.start();

            // for (int bin_index = 0; bin_index < 1<<log_slots; bin_index++) {
            //     int hashed_val = server_hash(i, bin_index, effective_bitLength);
            //     vector<int> bit_val = PerfectMapping(hashed_val, ell, hw);
            //     for (int bit_index = 0; bit_index < ell; bit_index++) {
            //         operand_vec[bit_index][bin_index] = bit_val[bit_index];
            //     }
            // } 

            // std::vector<Plaintext> encode_message;

            // HEaaN::Complex complex{0};
            // Message msg(log_slots, complex);
            
            // for (int l=0; l<ell; l++){
            //     for (int k=0; k<msg.getSize(); k++) {
            //         msg[k].real(operand_vec[l][k]); 
            //     }

            //     encode_message.push_back(endecoder.encode( msg, 
            //                             init.getLevel(), init.getRescaleCounter()));
            // }

            // key_time = key_timer1.end_and_get();
            // std::cout << key_time << " usec" << std::endl << std::endl;


            // std::cout << "Multiply and Add" << std::endl;
            // key_timer1.start();

            Ciphertext ctxt(context), ct_temp(context), temp(context);
            evaluator.add(init_ct, init_ct, ctxt);

            // Timer_micro tt;
            // long tt_time;

            for (int l=0; l<ell; l++){
                // Timer tempTime;
                // long temp_time;

                // tempTime.start();

                // Ciphertext temp(context);

                // std::cout << "1 mult" << std::endl;
                // tt.start();

                evaluator.mult(inputs[l], encode_message[i][l], ct_temp);
                // std::cout << "ct_temp level : " << ct_temp.getLevel() << std::endl;
                // std::cout << "ct_temp level : " << ctxt.getLevel() << std::endl;
                
                // tt_time = tt.end_and_get();
                // std::cout << tt_time << " usec" << std::endl;

                // sum all individual ciphertext vectors in temp_ct to obtain h' (ctxt)
                
                // tempTime.start();

                // std::cout << "1 add" << std::endl;
                // tt.start();

                evaluator.add(ctxt, ct_temp, ctxt);

                // tt_time = tt.end_and_get();
                // std::cout << tt_time << " usec" << std::endl;

                // temp_time = tempTime.end_and_get();
                // std::cout << temp_time << " ms" << std::endl;
            }

            // key_time = key_timer1.end_and_get();
            // std::cout << key_time << " usec" << std::endl << std::endl;

            // add Chebyshev Basis Evaluation Function HERE
            // input : ciphertext ctxt, poly coefficients (potentially using a lookup table for small domains)
            // output : ciphertext ctxt (or other ciphertext storage)

            // Alternative: Multiply like Tree

            // std::cout << "Poly Eval" << std::endl;
            // key_timer1.start();

            vector<Ciphertext> subtraction(hw, init_ct);
            for(int i = 0; i < hw; i++){
                // Ciphertext temp(context);
                evaluator.sub(ctxt, i, ct_temp);
                evaluator.mult(ct_temp, 1.0/(hw-i), temp);
                evaluator.add(temp, 0, subtraction[i]);
            }

            evaluator.add(multiplyLikeTree(evaluator, subtraction, 0, hw-1), 0, intersection[i]);
            
            // divide by 1/hw!
            // evaluator.add(result, 0, intersection[i]);

            // key_time = key_timer1.end_and_get();
            // std::cout << key_time << " usec" << std::endl << std::endl;

            

            // bootstrapping? 
            // Using FGb parameters, we reach level 6 for result ciphertexts using hamming weight of 3. 
            // level 2 for hamming weight of 5
        }
        // ================================

        // Check the case that only find intersection

        // Ciphertext size(context);
        // Message zero_msg(log_slots); 
        // fillReal(zero_msg, 0.0);
        // encryptor.encrypt(zero_msg, pack, size);
        // evaluator.levelDown(size, intersection[0].getLevel(), size);

        evaluator.levelDown(init_ct, intersection[0].getLevel(), init_ct);

        // std::cout << "result level : "<< intersection[0].getLevel() << std::endl;


        // Timer_micro key_timer1;

        // std::cout << "Combine" << std::endl;
        // key_timer1.start();

        // #pragma omp parallel for
        for (int i = 0; i < server_bin_size; i++){
            // #pragma omp critical

            evaluator.add(intersection[i], init_ct, init_ct);
        }

        // key_time = key_timer1.end_and_get();
        // std::cout << key_time << " usec" << std::endl << std::endl;

    key_time = key_timer.end_and_get();
    std::cout << key_time << " usec" << std::endl;


    // file << key_time << std::endl;
    // file.close();


    // compute the size of intersection
    // std::cout << "Compute intersection size ... ";
    // compute_sum(size);
    // std::cout << "done" << std::endl;

    // DEBUGGING
    Decryptor decryptor(context);
    Message dmsg(log_slots);
    decryptor.decrypt(init_ct, sk, dmsg);
    printMessage(dmsg);

    // compute average
    // std::cout << "Compute average ... ";
    // compute_sum(final);
    // evaluator.mult(final, 1.0 / (1 << log_slots), final);
    // std::cout << "done" << std::endl;

    
    // compute standard deviation on vector (create a function)


    // // encrypt and send over (create a function)
    // send_result(final);
    send_result(init_ct);
    cout << "Server side completed." << endl;

    return move(data_stream);
}



stringstream Server::serverMultipleLabelComp(Params &params){
    cout << "Server Multiple Label computation begins" << endl;

    int server_bin_size = params.server_bin_size, effective_bitLength = params.effective_bitLength, 
    ell = params.ell, hw = params.hw;
    
    const auto log_slots = getLogFullSlots(context);

    vector<Ciphertext> inputs;

    //loading ciphertexts
    Ciphertext ct_temp(context);
    for (int i=0;i<params.ell;i++){
        ct_temp.load(data_stream);
        inputs.push_back(ct_temp);
    }
    data_stream.str("");
    std::cout << "Saved number of ciphertexts : " << inputs.size() << std::endl;


    // create labels for Circuit + "Mutiple" Labelled PSI
    // using single ciphertext for simulation

    // Assume that the second label consists of 3 types. 
    // Test the case that we choose 2 types among the 3 labels.
    
    // like choosing Korean and Japanese among Korean, Japanese, and Chinese.

    vector<string> label1((1<<log_slots));
    Message mask(log_slots);
    Message label_val(log_slots);
    // int elt_idx[SIZE];

    for(int elt = 0; elt < (1<<log_slots); elt++){
        if(elt < 10) label1[elt] = (elt % 2 == 0) ? "Korean" : "Chinese";
        else if(elt < 20) label1[elt] = (elt % 2) ? "Japanese" : "Chinese";
        else label1[elt] = "Chinese";
    }

    for(int bin_index = 0; bin_index < (1<<log_slots); bin_index++){
        if(label1[bin_index].compare("Korean") == 0) {
            mask[bin_index].real(1.0);
            mask[bin_index].imag(0.0);

            label_val[bin_index].real(10.0);
            label_val[bin_index].imag(0.0);
        }
        else if(label1[bin_index].compare("Japanese") == 0) {
            mask[bin_index].real(1.0);
            mask[bin_index].imag(0.0);

            label_val[bin_index].real(5.0);
            label_val[bin_index].imag(0.0);
        }
        else {
            mask[bin_index].real(0.0);
            mask[bin_index].imag(0.0);

            label_val[bin_index].real(0.0);
            label_val[bin_index].imag(0.0);
        }
    }

    
    // subtraction vector of [hw]
    vector<Message> subtraction(hw-1);
    for (double i=1;i<hw;i++){
        Message msg(log_slots);
        fillReal(msg, i);
        subtraction[i-1] = msg;
    }
    vector<Message> factorial(hw-1);
    for (double i=2;i<hw+1;i++) {
        Message msg(log_slots);
        fillReal(msg, 1.0/i);
        factorial[i-2] = msg;
    }

    // store final results here
    vector<Ciphertext> intersection;
    vector<Ciphertext> intersection_cheby;

    Message zero_msg(log_slots); 
    fillReal(zero_msg, 0.0);

    // prepare for cheby setting
    const std::vector<Real> chebyCoefDeg3 = ChebyCoefTable::ChebyCoefDeg[3];
    u64 num_baby_step;

    if(chebyCoefDeg3.size() <= 3) num_baby_step = chebyCoefDeg3.size();
    else num_baby_step = pow(2,floor(ceil(log2(chebyCoefDeg3.size()))/2));

    // beginning of main loop
    for (int i=0; i<server_bin_size; i++) {
        vector<vector<int>> operand_vec(ell,vector<int>((1<<log_slots), 0));

        for (int bin_index = 0; bin_index < (1<<log_slots); bin_index++) {
            // In this case, bin_index is an element, or id 
            // that is compared to the element of client.

            // Also, server_hash is permutatoin based hashing that just 
            // compresses a size of the input. 
            // (Do not change the location!)

            int hashed_val = server_hash(i, bin_index, effective_bitLength);

            vector<int> bit_val = PerfectMapping(hashed_val, ell, hw);
            for (int bit_index = 0; bit_index < ell; bit_index++) {
                operand_vec[bit_index][bin_index] = bit_val[bit_index];
            }
        } 
        
        Ciphertext ctxt(context);
        encryptor.encrypt(zero_msg, pack, ctxt);
        

        for (int l=0; l<ell; l++){
            HEaaN::Complex complex{0};
            Message msg(log_slots, complex);
            Ciphertext temp(context);

            for (int k=0; k<msg.getSize(); k++) {
                msg[k].real(operand_vec[l][k]);
                // if (k == msg.getSize()-1 && i == 4){
                //     std::cout << std::endl << "Server Message : " << std::endl;
                //     printMessage(msg);
                //     std::cout << std::endl;
                // }
            }

            evaluator.mult(inputs[l], msg, temp);
            // std::cout << "Result Ciphertext - level " << temp.getLevel() << std::endl;

            // sum all individual ciphertext vectors in temp_ct to obtain h' (ctxt)
            evaluator.add(ctxt, temp, ctxt);
        }

        // add Chebyshev Basis Evaluation Function HERE
        // input : ciphertext ctxt, poly coefficients (potentially using a lookup table for small domains)
        // output : ciphertext ctxt (or other ciphertext storage)

        Ciphertext ctxt_test(context);
        ChebyshevCoefficients chebyshev_coef_deg_3(chebyCoefDeg3, num_baby_step);

        if (i==0){
            std::cout << "Ciphertext - initial level " << ctxt.getLevel() << std::endl;
        }
        
        evaluator.add(evaluateChebyshev(evaluator, ctxt, 
            chebyshev_coef_deg_3, 1.0/6), 0, ctxt_test);

        if (i==0){
            std::cout << "Chebyshev Ciphertext - post-comp level " << ctxt_test.getLevel() << std::endl;
        }

        intersection_cheby.push_back(ctxt_test);


        // Alternative: Direct Product Evaluation

        // Ciphertext temp(context);
        // Ciphertext result(context);

        // evaluator.add(ctxt, 0, result);

        // // if (i==0){
        // //     std::cout << "Result Ciphertext - initial level " << result.getLevel() << std::endl;
        // // }

        // for (int k=0;k<hw-1;k++){
        //     evaluator.sub(ctxt, subtraction[k], temp);
        //     evaluator.mult(result, temp, result);
        //     evaluator.mult(result, factorial[k], result);
        // }
        // // store final value in mu_result[i] and multiply by the label
        // intersection.push_back(result);

        // if (i==0){
        //     std::cout << "Result Ciphertext - post-comp level " << result.getLevel() << std::endl;
        // }
    }

    // sum each elt of intersection to get the total intersection.
    Ciphertext size(context);
    encryptor.encrypt(zero_msg, pack, size);

    for (int i = 0; i < server_bin_size; i++){
        evaluator.add(size, intersection_cheby[i], size);
    }

    // multiply mask(ptxt) and compute the size of intersection and average

    Ciphertext temp(context);
    evaluator.mult(size, mask, temp);

    std::cout << "Compute intersection size ... ";
    compute_sum(temp);
    std::cout << "done" << std::endl;

    // DEBUGGING
    Decryptor decryptor(context);
    Message dmsg(log_slots);
    decryptor.decrypt(temp, sk, dmsg);
    printMessage(dmsg);

    // compute average
    std::cout << "Compute average ... ";

    Ciphertext final(context), squareEx(context);
    evaluator.mult(size, label_val, final);
    evaluator.mult(final, final, squareEx);
    
    Ciphertext tempInv(context);
    approxInverseNewtonX(evaluator, temp, tempInv, 0.05, 4);

    compute_sum(final);

    std::cout << std::endl;   
    std::cout << "Use approxInverse" << std::endl; 

    evaluator.mult(final, tempInv, temp);
    std::cout << "Result level : " << temp.getLevel() << std::endl;

    decryptor.decrypt(temp, sk, dmsg);
    printMessage(dmsg);

    std::cout << std::endl;

    // evaluator.mult(final, 1.0 / 10, final);
    // std::cout << "done" << std::endl;

    // Compute E(X^2) to find the standard deviation

    Ciphertext temp2(context);
    compute_sum(squareEx);

    std::cout << "Compute E(X^2) using approxInverse" << std::endl;

    evaluator.mult(squareEx, tempInv, temp2);
    std::cout << "level after computing E(X^2) : " << temp2.getLevel() << std::endl;

    evaluator.mult(temp, temp, temp);
    evaluator.sub(temp2, temp, temp); // E(X^2) - E(X)^2

    // Boostrapping
    // btp.bootstrap(temp, temp, false);    
    std::cout << "level after bootstrapping : " << temp.getLevel() << std::endl;

    std::cout << "Message before computing sqrt" << std::endl;

    decryptor.decrypt(temp, sk, dmsg);
    printMessage(dmsg);

    std::cout << std::endl;

    evaluator.mult(temp, 1.0/9, temp);
    approxSqrtWilkesX(evaluator, temp, temp, 3);
    evaluator.mult(temp, 3, temp);

    std::cout << "level after computing the standard deviation: " << temp2.getLevel() << std::endl;

    decryptor.decrypt(temp, sk, dmsg);
    printMessage(dmsg);

    std::cout << std::endl;
    

    // // encrypt and send over (create a function)
    send_result(final);
    cout << "Server side completed." << endl;

    return move(data_stream);
}
