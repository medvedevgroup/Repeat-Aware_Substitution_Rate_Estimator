#include <iostream>
#include <random>
#include <string>
#include <set>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <map>
#include <unordered_map>
using namespace std;


#include "mpreal.h"
using mpfr::mpreal;
// Required precision of computations in decimal digits
const int digits = 50;  // Play with it to check different precisions



set<string> kspectrum(string s, int k){
    set<string> kmers;
    for(int i=0; i<s.length()-k+1; i++){
        string kmer = s.substr(i, k);
        kmers.insert(kmer);
    }
    return kmers;
}


mpreal Polynomial (const unordered_map<int,int>& coeffs_Map, mpreal q){
    mpreal::set_default_prec(mpfr::digits2bits(digits));
    mpreal poly = 0;
    for (const auto& [exp, mul] : coeffs_Map) {
        poly += mul * pow(q , exp);
    }
    return poly;
}

mpreal Polynomial_prime (const unordered_map<int,int>& coeffs_Map, mpreal q){
    mpreal::set_default_prec(mpfr::digits2bits(digits));
    mpreal poly_prime = 0;
    for (const auto& [exp, mul] : coeffs_Map) {
        poly_prime += exp * mul * pow(q , exp -1);
    }
    return poly_prime;
}

mpreal Newton_Method (const unordered_map<int, int>& coeffs_Map, const set<string>& ksp_set, int k, double epsilon, double I_obs){
    mpreal::set_default_prec(mpfr::digits2bits(digits));
    //set<string> ksp_set = kspectrum(seq, k);
    
    // for (const auto& pair : coeffs_Map) {
    //     std::cout << "Key: " << pair.first << ", Value: " << pair.second << std::endl;
    // }
    //cout << "a: " << (ksp_set.size() - I_obs) << endl;
    mpreal q_hat = 1 - pow(0.09, k); 
    mpreal f = Polynomial(coeffs_Map, q_hat) - (ksp_set.size() - I_obs);      // P(q) - a
    while (abs(f) > epsilon) {
        f = Polynomial(coeffs_Map, q_hat) - (ksp_set.size() - I_obs);      // P(q) - a
        mpreal fp = Polynomial_prime (coeffs_Map, q_hat);   // P'(q)
        q_hat -= f / fp;
        //cout << "f: " << f << endl;
        //cout << "fp: " << fp << endl;
        //cout << "q_hat: " << q_hat << endl;
    }


    return q_hat;
}

