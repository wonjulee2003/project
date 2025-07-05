#include <iostream>
#include <math.h>
#include <random>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>

#include <stdlib.h>
#include <time.h>

typedef long long int lld;
int C[1001][1001];

/*
b : bin size
h : hamming weight
lambda : bit length after Perm-based hashing 
l : bit length after Constant Weight Encoding

rehash 몇 번을 하면 좋을지 parameter 설정해야함. 
*/

lld MakeRandomVal(int lower, int upper);
void assignVal(int *arr_ptr, int size);
void printArr(int *arr_ptr, int size);
lld HashFunc(lld input, lld coef, lld mod, lld div);
int PermBasedHashing(int input, int leftBitSize, int rightBitSize, lld coef);
void CuckooHashing(int *arr_ptr, std::unordered_map<int,int> &idx_ptr, int *bin_ptr, int elementBitSize, int binBitSize, const int MAXLOOP, int size, int binSize);
void CuckooHashingWithFixedCoef(std::vector<lld> coefs, int *arr_ptr, std::unordered_map<int,int> &idx_ptr, int *bin_ptr, 
                    int elementBitSize, int binBitSize, const int MAXLOOP, int size, 
                    int binSize);

int PerfectMapping(int input, int bitLength, int hammingWeight);
int comb(int n, int r);

bool test = false;
int rehashCnt = 0;

int main(void){
    // 1.

    // int arr[12] = {
    //     804165511, 1039738916, 498641203, 8771098, 
    //     850867609, 737815609, 242774061, 947127015,
    //     453724639, 756848357, 135097, 2348
    // };

    // int bin[16] = {};

    // int leftBitSize = 4, rightBitSize = 28;
    // lld temp = 2 * MakeRandomVal(0, 1 << (rightBitSize - 1) - 1) + 1;

    // for(int i = 0; i < 12; i++){


    //     PermBasedHashing(arr[i], 4, 28);
    // }

    // 2.

    // int arr[12] = {
    //     1298,
    //     11320,
    //     22179,
    //     26560,
    //     29566,
    //     629,
    //     12460,
    //     32671,
    //     1332,
    //     24473,
    //     1960,
    //     16310
    // };

    // 3.

    // int arr[12] = {
    //     1462,
    //     24452,
    //     30651,
    //     17310,
    //     21006,
    //     358,
    //     4242,
    //     17379,
    //     7712,
    //     11982,
    //     27332,
    //     18513
    // };

    // int bin[16] = {};
    // std::unordered_map<int,int> elt_idx;

    // int elementBitSize = 32, binBitSize = 4;
    // const int MAXLOOP = 11; // ln(arrSize) / ln(1.27)
    // int size = sizeof(arr)/sizeof(int);

    // std::cout << "Start" << std::endl << std::endl;

    // CuckooHashing(arr, elt_idx, bin, elementBitSize, binBitSize, MAXLOOP, size);

    // for(int i = 0; i < 16; i++){
    //     std::cout << "index : " << i << ", value : " << bin[i] << std::endl;
    // }


    // =========================================
    
    // 2. To check cuckoo hahsing works well. (random arr and random coefs)

    // 12 / 16
    // int inputSize = 4096;
    // int binSize = 8192;
    // int iterNum = 1<<30;

    // int arr[inputSize];
    // int bin[binSize];

    // for(int i = 0; i < binSize; i++){
    //     bin[i] = 0;
    // }

    // std::unordered_map<int,int> elt_idx;

    // int elementBitSize = 32, binBitSize = (int)log2(binSize);
    // const int MAXLOOP = 35; // ln(arrSize) / ln(1.27) = 11
    // int size = sizeof(arr)/sizeof(int);

    // std::cout << "Start" << std::endl << std::endl;

    // int j;
    // for(j = 1; j <= iterNum; j++){ 
    //     // j = 26 일 때 rehash 가 필요한 경우 발생
    //     // lld temp = 2 * MakeRandomVal(0, (1 << (elementBitSize - 1)) - 1) + 1 로 할 경우

    //     std::cout << j << "th execution" << std::endl << std::endl;
    //     for(int i = 0; i < inputSize; i++){
    //         int temp = MakeRandomVal(1,1<<31-1);
    //         arr[i] = temp;
    //         // std::cout << temp << std::endl;
    //     }

    //     CuckooHashing(arr, elt_idx, bin, elementBitSize, binBitSize, MAXLOOP, size, binSize);
        
    //     if(test){
    //         break;
    //     }

    //     for(int i = 0; i < inputSize; i++){
    //         arr[i] = 0;
    //     }

    //     for(int i = 0; i < binSize; i++){
    //         bin[i] = 0;
    //     }

    //     elt_idx.clear();

    //     // std::cout << "===================" << std::endl << std::endl;
    // }

    // std::cout << j-1 << " times succeed !" << std::endl;
    // std::cout << "Rehash count : " << rehashCnt << std::endl;

    // ============================================

    // 3. To check rehash occurs for fixed coefs. (random arr and random coefs)

    int inputSize = 1024; // 4096, 8192
    int binSize = 2048;
    int iterNum = 1<<20;

    int arr[inputSize];
    int bin[binSize];

    for(int i = 0; i < binSize; i++){
        bin[i] = 0;
    }

    std::unordered_map<int,int> elt_idx;

    int elementBitSize = 32, binBitSize = (int)log2(binSize);
    const int MAXLOOP = 35; // ln(arrSize) / ln(1.27) = 11

    std::vector<lld> coefs;
    while(coefs.size() != 3){
        lld temp = 2 * MakeRandomVal(0, (1 << (elementBitSize - 1)) - 1) + 1;
        // lld temp = 2 * MakeRandomVal(0, 1 << (elementBitSize - 1) - 1) + 1;

        if(find(coefs.begin(), coefs.end(), temp) == coefs.end()){
            coefs.push_back(temp);
        }
    }

    std::cout << "Start" << std::endl;
    std::cout << "Test " << iterNum << " times!" << std::endl << std::endl;

    int j;
    for(j = 1; j <= iterNum; j++){ 
        // j = 26 일 때 rehash 가 필요한 경우 발생
        // lld temp = 2 * MakeRandomVal(0, (1 << (elementBitSize - 1)) - 1) + 1 로 할 경우

        std::cout << j << "th execution" << std::endl << std::endl;
        assignVal(arr, inputSize);
        // printArr(arr, inputSize);
        CuckooHashingWithFixedCoef(coefs, arr, elt_idx, bin, elementBitSize, binBitSize, MAXLOOP, inputSize, binSize);
        
        if(test){
            break;
        }

        for(int i = 0; i < inputSize; i++){
            arr[i] = 0;
        }

        for(int i = 0; i < binSize; i++){
            bin[i] = 0;
        }

        elt_idx.clear();

        // std::cout << "===================" << std::endl << std::endl;
    }

    std::cout << j-1 << " times succeed !" << std::endl;
    printArr(arr, inputSize);

    // srand(time(NULL));

    // for(int i = 0; i < 12; i++){
    //     std::cout << rand() << std::endl;
    // }
    
    // return 0;
}

lld MakeRandomVal(int lower, int upper){
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_int_distribution<lld> dis(lower,upper);

    return dis(gen);
}

void assignVal(int *arr_ptr, int size){
    std::unordered_set<int> cnt;

    for(int i = 0; i < size; i++){
        while(true){
            int temp = MakeRandomVal(1,1<<31-1);

            if(cnt.find(temp) == cnt.end()){
                arr_ptr[i] = temp;
                cnt.insert(temp);
                break;
            }
        }        
    }

}

void printArr(int *arr_ptr, int size){
    for(int i = 0; i  < size; i++){
        std::cout << arr_ptr[i] << std::endl;
    }
    std::cout << std::endl;
}

lld HashFunc(lld input, lld coef, lld mod, lld div){
    // mod = 1 << bitlength of input
    // div = 1 << (bitlength of input - bitlength of ouput)

    return ((coef*input) % mod) / div;
}

int PermBasedHashing(int input, int leftBitSize, int rightBitSize, lld coef){ 
    // permutation based hashing for reducing the input size
    // input : x1 + x2 bit -> output : x2 bit
    
    // Assume that the right bit size is bigger than the left bit size 

    lld mod = (lld)1 << rightBitSize;
    lld div = (lld)1 << (rightBitSize - leftBitSize);

    int leftval = input >> rightBitSize;
    lld rightval = (input << leftBitSize) >> leftBitSize;

    return leftval ^ HashFunc(rightval, coef, mod, div);
}

void CuckooHashing( int *arr_ptr, std::unordered_map<int,int> &idx_ptr, int *bin_ptr, 
                    int elementBitSize, int binBitSize, const int MAXLOOP, int size, 
                    int binSize){
    // Using 3 hash functions

    // Assume that the binSize is a power of 2.
    // and also assume that elements of input are not 0.
    
    lld mod = (lld)1 << elementBitSize;
    lld div = (lld)1 << (elementBitSize - binBitSize);

    int cnt = 1;

    while(true){
        std::vector<lld> coefs;

        while(coefs.size() != 3){
            lld temp = 2 * MakeRandomVal(0, (1 << (elementBitSize - 1)) - 1) + 1;
            // lld temp = 2 * MakeRandomVal(0, 1 << (elementBitSize - 1) - 1) + 1;

            if(find(coefs.begin(), coefs.end(), temp) == coefs.end()){
                coefs.push_back(temp);
            }
        }
        
        int i = 0;

        while(i < size){
            int element = arr_ptr[i];
            int idx = idx_ptr[element];
            bool flag = true;

            // insert 

            for(int j = 0; j < MAXLOOP; j++){
                lld coef = coefs[idx];
                lld hash_val = HashFunc(element, coef, mod, div);

                // std::cout << element << " " << hash_val << std::endl;

                if(bin_ptr[hash_val] == 0){
                    bin_ptr[hash_val] = element;
                    idx_ptr[element] = idx;
                    i++;
                    flag = false;

                    break;
                }
                else{
                    int temp_val = bin_ptr[hash_val];
                    int temp_idx = idx_ptr[temp_val];
                    
                    bin_ptr[hash_val] = element;
                    idx_ptr[element] = idx;

                    element = temp_val;
                    idx = (temp_idx + 1) % 3;
                }
            }

            // if flag = true, then we require rehash
            if(flag) break; 

            // std::cout << "Insert one element" << std::endl;

        }

        if(i == size) break;
        else{
            std::cout << "We need rehash!" << std::endl;
            
            for(int i = 0; i < binSize; i++){
                bin_ptr[i] = 0;
            }

            idx_ptr.clear();

            if(cnt == 10){
                test = true;
                break;
            }
            cnt++;
            rehashCnt++;
        }
    }
    
    return ;
}

void CuckooHashingWithFixedCoef(std::vector<lld> coefs, int *arr_ptr, std::unordered_map<int,int> &idx_ptr, int *bin_ptr, 
                    int elementBitSize, int binBitSize, const int MAXLOOP, int size, 
                    int binSize){
    // Using 3 hash functions

    // Assume that the binSize is a power of 2.
    // and also assume that elements of input are not 0.
    
    lld mod = (lld)1 << elementBitSize;
    lld div = (lld)1 << (elementBitSize - binBitSize);

    int i = 0;

    while(i < size){
        int element = arr_ptr[i];
        int idx = idx_ptr[element];
        bool flag = true;

        // insert 

        for(int j = 0; j < MAXLOOP; j++){
            lld coef = coefs[idx];
            lld hash_val = HashFunc(element, coef, mod, div);

            // std::cout << element << " " << hash_val << std::endl;

            if(bin_ptr[hash_val] == 0){
                bin_ptr[hash_val] = element;
                idx_ptr[element] = idx;
                i++;
                flag = false;

                break;
            }
            else{
                int temp_val = bin_ptr[hash_val];
                int temp_idx = idx_ptr[temp_val];
                
                bin_ptr[hash_val] = element;
                idx_ptr[element] = idx;

                element = temp_val;
                idx = (temp_idx + 1) % 3;
            }
        }

        // if flag = true, then we require rehash
        if(flag) break; 

        // std::cout << "Insert one element" << std::endl;

    }

    if(i < size) test = true;
            
    return ;
}

int PerfectMapping(int input, int bitLength, int hammingWeight){ // Constant-Weight Encoding
    int r = input, h = hammingWeight, ret = 0;
    int idx = 1 << (bitLength -1);

    for(int n = bitLength - 1; n >= 0; n--){      
        if(n < h){
            ret |= idx;
            h--;
        }else if(r >= comb(n,h)){
            ret |= idx;
            r -= comb(n,h);
            h--;
        }

        if(h == 0) break;
        idx >>= 1;
    }

    return ret;
}

int comb(int n, int r){
    // We can assume that n >= r.
     
    if(n == r || r == 0){
        return 1;
    }
    else if(C[n][r] != 0){
        return C[n][r];
    }

    int ret = comb(n-1,r) + comb(n-1,r-1);
    C[n][r] = ret;
    return ret;
}


// bool Lookup(std::vector<lld> coefs, int *bin_ptr, lld mod, lld div, int element){
//     To remove the repeated values

//     for(lld coef : coefs){
//         if(bin_ptr[HashFunc(coef,element,mod,div)] == 0){
//             return true;
//         }
//     }
//     return false;
// }

// bool Insert(std::vector<lld> coefs, int *bin_ptr, lld mod, lld div, int element, int MAXLOOP){
//     if (Lookup(coefs, bin_ptr, mod, div, element)){
        

//     }  
// }

