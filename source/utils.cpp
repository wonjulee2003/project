#include "utils.hpp"


using namespace std;
using namespace HEaaN;

int comb(int n, int k) {
    if (k > n) {
        return 0;
    }
    int res = 1;
    for (int i {1}; i <= k; i++) {
        res *= n--;
        res /= i;
    }
    return res;
}

float log2_comb(int n, int k) {
    if (k > n) {
        return 0;
    }
    float r = 0.0;
    for (uint64_t d = 1; d <= k; ++d) {
        r += log2(n--)-log2(d);
    }
    return r;
}

vector<int> PerfectMapping(int input, int bitLength, int hammingWeight){
    int r = input, h = hammingWeight;
    vector<int> res(bitLength, 0);
    int mod_size = comb(bitLength, h);
    if (r >= mod_size){
            cout << "Overflow" << r << endl;
        r %= mod_size;
    }
    for (int n = bitLength-1; n>=0; n--){
        if (r >= comb(n, h)){
            res[n] = 1;
            r -= comb(n, h);
            h --;
        }
        if(h==0) break;
    }
    return res;
}

void printMessage(const HEaaN::Message &msg, bool is_complex,
                         size_t start_num, size_t end_num,
                         const int precision) {
    const size_t msg_size = msg.getSize();
    std::cout.precision(precision);
    std::cout << "[ ";
    for (size_t i = 0; i < start_num; ++i) {
        if (is_complex)
            std::cout << msg[i] << ", ";
        else
            std::cout << msg[i].real() << ", ";
    }
    std::cout << "..., ";
    for (size_t i = end_num; i > 1; --i) {
        if (is_complex)
            std::cout << msg[msg_size - i] << ", ";
        else
            std::cout << msg[msg_size - i].real() << ", ";
    }
    if (is_complex)
        std::cout << msg[msg_size - 1] << " ]" << std::endl;
    else
        std::cout << msg[msg_size - 1].real() << " ]" << std::endl;
}

double randNum() {
    static std::default_random_engine gen{std::random_device()()};
    std::uniform_real_distribution<double> dist(-1.0L, 1.0L);
    return dist(gen);
}


void fillRandomReal(HEaaN::Message &msg) {
    for (size_t idx = 0; idx < msg.getSize(); idx ++) {
        msg[idx].real(randNum());
        msg[idx].imag(0.0);
    }
}

void fillReal(HEaaN::Message &msg, double val) {
    for (size_t idx = 0; idx < msg.getSize(); idx ++) {
        msg[idx].real(val);
        msg[idx].imag(0.0);
    }
}

Params::Params(int server_bin_size, int effective_bitLength, int hw)
    : server_bin_size{server_bin_size}, effective_bitLength{effective_bitLength}, ell{hw}, hw{hw}
     {
        while (log2_comb(ell, hw) < effective_bitLength){
            ell++;
        }
        // bitLength = effective_bitLength + log_modulus_degree
     }

void approxInverseNewton(const HomEvaluator &eval,
                         const Ciphertext &ctxt, Ciphertext &ctxt_out,
                         Real initial, u64 num_iter) {
    const u64 one_iter_cost{1};

    Ciphertext ctxt_x(eval.getContext()), ctxt_y(eval.getContext());
    Ciphertext ctxt_z(eval.getContext()), ctxt_tmp(eval.getContext());

    ctxt_x = ctxt;

    // y_0 = initial
    // z_0 = x * y_0
    // y_1 = y_0 * (2 - z_0)
    eval.mult(ctxt_x, initial, ctxt_z);
    eval.negate(ctxt_z, ctxt_tmp);
    eval.add(ctxt_tmp, 2.0, ctxt_tmp); // tmp = 2 - z_0
    eval.mult(ctxt_tmp, initial, ctxt_y);

    // for n > 0
    // z_n = x * y_n
    //     = z_{n-1} * (2 - z_{n-1})
    // y_{n+1} = y_n * (2 - z_n)
    for (u64 iter = 1; iter < num_iter; iter++) {
        eval.mult(ctxt_z, ctxt_tmp, ctxt_z);

        eval.negate(ctxt_z, ctxt_tmp); // tmp = 2 - z_n
        eval.add(ctxt_tmp, 2.0, ctxt_tmp);
        eval.mult(ctxt_y, ctxt_tmp, ctxt_y);
    }

    ctxt_out = ctxt_y;
}

void approxSqrtWilkes(const HomEvaluator &eval, 
                    //   const Bootstrapper &btp,
                      const Ciphertext &ctxt, Ciphertext &ctxt_out,
                      const u64 num_iter) {
    Ciphertext tmp(eval.getContext());
    Ciphertext tmp2(eval.getContext());
    static const u64 ONE_ITER_COST = 2;

    ctxt_out = ctxt;

    // if (ctxt_out.getLevel() < 1 + btp.getMinLevelForBootstrap()) {
    //     btp.bootstrap(ctxt_out, ctxt_out);
    // }
    eval.sub(ctxt_out, 1, tmp);
    eval.mult(tmp, 0.5, tmp);

    for (u64 i = 0; i < num_iter; i++) {
        // complex bootstrap
        // if (tmp.getLevel() < ONE_ITER_COST + btp.getMinLevelForBootstrap()) {
        //     eval.multImagUnit(tmp, tmp2);
        //     eval.add(ctxt_out, tmp2, tmp2);
        //     btp.bootstrap(tmp2, ctxt_out, tmp);
        // }

        // compute a_(k+1)
        eval.negate(tmp, tmp2);
        eval.add(tmp2, 1, tmp2);
        eval.mult(ctxt_out, tmp2, ctxt_out);

        // compute h_(k+1)
        eval.sub(tmp, 1.5, tmp2);
        eval.square(tmp, tmp);
        eval.mult(tmp, tmp2, tmp);
    }
}