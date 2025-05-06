#include "utils.hpp"


using namespace std;

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