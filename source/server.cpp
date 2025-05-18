#include "server.hpp"

using namespace std;
using namespace HEaaN;
// using namespace HEaaN::Math;

int Server::server_hash(int mu, int bin_index, int bitLength){
    return ((bin_index+1)+mu) % (1 << bitLength);
}

// constructor
Server::Server(const string &string1, string &string2, stringstream &&source_stream, string &string4)
    : context(makeContextFromFile(string1)), pack(context, string2), evaluator(context, pack), encryptor(context),
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
        Message zero_msg(log_slots); 
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


    // sum each index of mu_result[i]
    Ciphertext size(context);
    Message zero_msg(log_slots); 
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
    const std::vector<Real> chebyCoefDeg3 = {
	-1.5, 2.75, -1.5, 0.25};
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

        Ciphertext temp(context);
        Ciphertext result(context);

        evaluator.add(ctxt, 0, result);

        // if (i==0){
        //     std::cout << "Result Ciphertext - initial level " << result.getLevel() << std::endl;
        // }

        for (int k=0;k<hw-1;k++){
            evaluator.sub(ctxt, subtraction[k], temp);
            evaluator.mult(result, temp, result);
            evaluator.mult(result, factorial[k], result);
        }
        // store final value in mu_result[i] and multiply by the label
        intersection.push_back(result);

        if (i==0){
            std::cout << "Result Ciphertext - post-comp level " << result.getLevel() << std::endl;
        }


        // bootstrapping? 
        // Using FGb parameters, we reach level 6 for result ciphertexts using hamming weight of 3. 
        // level 2 for hamming weight of 5
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
    approxInverseNewton(evaluator, temp, tempInv, 0.05, 4);

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
    
    std::cout << "Message before computing sqrt" << std::endl;

    decryptor.decrypt(temp, sk, dmsg);
    printMessage(dmsg);

    std::cout << std::endl;

    evaluator.mult(temp, 1.0/9, temp);
    approxSqrtWilkes(evaluator, temp, temp2, 3);
    evaluator.mult(temp2, 3, temp2);

    std::cout << "level after computing the standard deviation: " << temp2.getLevel() << std::endl;

    decryptor.decrypt(temp2, sk, dmsg);
    printMessage(dmsg);

    std::cout << std::endl;
    
    // compute standard deviation on vector (create a function)


    // // encrypt and send over (create a function)
    send_result(final);
    cout << "Server side completed." << endl;

    return move(data_stream);
}
