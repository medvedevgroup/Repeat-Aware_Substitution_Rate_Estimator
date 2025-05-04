#include <iostream>
#include <random>
#include <string>
#include <set>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <map>
#include <unordered_map>
#include "Newton_method.h"
#include "sketch.h"
#include <cmath>
using namespace std;


#include "mpreal.h"
using mpfr::mpreal;
// Required precision of computations in decimal digits
const int digits = 50;  // Play with it to check different precisions
unordered_map<string, set<int>> occurrence; // kmer, occ
set<string> ksp_set1; // kspec of origin string s
set<string> ksp_set2; // kspec of mutated string t
unordered_map<int, int> coeffs_Map; // a * x^b, b->a



set<string> kspectrum_update(string s, int k){
    set<string> kmers;
    for(int i=0; i<s.length()-k+1; i++){
        string kmer = s.substr(i, k);
        kmers.insert(kmer);
        occurrence[kmer].insert(i);
    }
    return kmers;
}


int intersect_size(const set<string>& set1, const set<string>& set2){
    std::vector<string> intersect;

    std::set_intersection(set1.begin(), set1.end(),
                          set2.begin(), set2.end(),
                          std::back_inserter(intersect));
    return intersect.size();
}

// string readSubstring(string filename, int start ,int n){
//     ifstream file(filename); 
//     string line;
//     string sequence;

//     while (getline(file, line)) {
//         if (line.empty()) continue;
//         if (line[0] == '>') { 
//             continue;
//         } else {
//             sequence += line; 
//         }
//     }

//     file.close();
//     string seq = sequence.substr(start,n);
//     transform(seq.begin(), seq.end(), seq.begin(),::toupper);
//     return seq;
// }

std::optional<unordered_map<int, int>> parse_csv(const std::string& filename) {
    std::ifstream file(filename);
    std::unordered_map<int, int> occ_map;
    
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file '" << filename << "'.\n";
        return std::nullopt;
    }
    
    std::string line;
    int line_number = 0;
    
    try {
        while (std::getline(file, line)) {
            line_number++;
            
            // Skip empty lines
            if (line.empty()) {
                continue;
            }
            
            std::istringstream iss(line);
            int occ, count;
            char comma;
            
            if (iss >> occ >> comma >> count) {
                if (comma == ',') {
                    occ_map[occ] = count;
                } else {
                    std::cerr << "Error: Invalid format at line " << line_number 
                              << ". Expected comma between values.\n";
                    return std::nullopt;
                }
            } else {
                std::cerr << "Error: Could not parse values at line " << line_number << ".\n";
                return std::nullopt;
            }
        }
        
        file.close();
        
        if (occ_map.empty()) {
            std::cerr << "Error: No valid occurrence data found in file '" << filename << "'.\n";
            return std::nullopt;
        }
        
        return occ_map;
        
    } catch (const std::exception& e) {
        std::cerr << "Error processing file '" << filename << "': " << e.what() << "\n";
        file.close();
        return std::nullopt;
    }
}

std::optional<set<string>> parse_sequence2kspec(string filename, bool mode, int k){ 
    // if mode is true, then we call `kspectrum_update` to parse distribution information
    // otherwise, call `kspectrum`, just get kmers
    ifstream file(filename);
    set<string> kspec;
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file '" << filename << "'.\n";
        return std::nullopt;
    }
    
    string line;
    string sequence;
    bool found_header = false;
    
    try {
        while (std::getline(file, line)) {
            // Skip empty lines
            if (line.empty()) {
                continue;
            }
            
            // Check if this is a header line
            if (line[0] == '>') {
                // skip, we only extract one
                if (found_header && !sequence.empty()) {
                    break;
                }
                
                found_header = true;
                continue;
            }
            
            // collect sequence data
            if (found_header) {
                sequence += line;
            }
        }
        
        file.close();
        
        // Check l(s) < k?
        if (sequence.empty()) {
            std::cerr << "Error: No sequence found in file '" << filename << "'.\n";
            return std::nullopt;
        }
        
        if (sequence.length() < k) {
            std::cerr << "Error: Sequence length (" << sequence.length() 
                      << ") is less than k (" << k << ") in file '" << filename << "'.\n";
            return std::nullopt;
        }
        
        if (mode){
            kspec = kspectrum_update(sequence, k);
        }
        else{
            kspec = kspectrum(sequence, k);
        }
        return kspec;
        
    } catch (const std::exception& e) {
        std::cerr << "Error processing file '" << filename << "': " << e.what() << "\n";
        file.close();
        return std::nullopt;
    }
}

std::optional<set<string>> parse_kmerfile2kspect(const std::string& filename, int k) {
    std::ifstream file(filename);
    std::set<std::string> kmers;
    
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file '" << filename << "'.\n";
        return std::nullopt;
    }
    
    std::string line;
    std::string current_kmer;
    bool reading_kmer = false;
    int line_number = 0;
    
    try {
        while (std::getline(file, line)) {
            line_number++;
            
            if (line.empty()) {
                continue;
            }
            
            if (line[0] == '>') {
                if (reading_kmer && current_kmer.length() != k) {
                    std::cerr << "Error: K-mer at line " << line_number 
                              << " has length " << current_kmer.length() 
                              << ", expected " << k << ".\n";
                    return std::nullopt;
                }
                kmers.insert(current_kmer);
                current_kmer.clear();
                reading_kmer = true;
                continue;
            }
            
            // read kmer
            if (reading_kmer) {
                current_kmer += line;     
            }
        }
        if (!current_kmer.empty()) {
            kmers.insert(current_kmer);
        }
        file.close();
        
        if (kmers.empty()) {
            std::cerr << "Error: No valid k-mers found in file '" << filename << "'.\n";
            return std::nullopt;
        }
        
        return kmers;
        
    } catch (const std::exception& e) {
        std::cerr << "Error processing file '" << filename << "': " << e.what() << "\n";
        file.close();
        return std::nullopt;
    }
}



void print_usage() {
    std::cout << "Usage: ./Mutation_rate_estimator [OPTIONS]\n"
              << "Options:\n"
              << "  --mode [sequence|mixture|kmer]    Operation mode\n"
              << "  --input1 Fasta FILE                     First input fasta\n"
              << "  --input2 Fasta FILE                     Second input fasta\n"
              << "  --dist FILE                       Distribution file (required only for kmer mode) in csv form\n"
              << "  --k INT                           k-mer length\n"
              << "  --theta FLOAT                     Sketching parameter (default: 1.0, i.e., without sketching)\n"
              << "  --e FLOAT                         Error_bound for Newton's Method (default: 1e-5)\n"
              << "  -h, --help                        Display this help message\n";
}

struct Parameters {
    std::string mode;
    std::string input1;
    std::string input2;
    std::string dist;
    int k = 0;
    double theta = 1.0;
    double error_bound = 1e-5;
    bool valid = true;
    bool help_printed = false;
};

Parameters parse_arguments(int argc, char* argv[]) {
    Parameters params;
    
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        
        if (arg == "--help" || arg == "-h") {
            print_usage();
            params.valid = false;
            params.help_printed = true;
            return params;
        } else if (arg == "--mode") {
            if (i + 1 < argc) {
                params.mode = argv[++i];
                if (params.mode != "sequence" && params.mode != "mixture" && params.mode != "kmer") {
                    std::cerr << "Error: Invalid mode. Use 'sequence', 'mixture', or 'kmer'.\n";
                    params.valid = false;
                }
            } else {
                std::cerr << "Error: --mode requires an argument.\n";
                params.valid = false;
            }
        } else if (arg == "--input1") {
            if (i + 1 < argc) {
                params.input1 = argv[++i];
            } else {
                std::cerr << "Error: --input1 requires a filename.\n";
                params.valid = false;
            }
        } else if (arg == "--input2") {
            if (i + 1 < argc) {
                params.input2 = argv[++i];
            } else {
                std::cerr << "Error: --input2 requires a filename.\n";
                params.valid = false;
            }
        } else if (arg == "--dist") {
            if (i + 1 < argc) {
                params.dist = argv[++i];
            } else {
                std::cerr << "Error: --dist requires a filename.\n";
                params.valid = false;
            }
        } else if (arg == "--k") {
            if (i + 1 < argc) {
                try {
                    params.k = std::stoi(argv[++i]);
                    if (params.k <= 0) {
                        std::cerr << "Error: k must be a positive integer.\n";
                        params.valid = false;
                    }
                } catch (const std::exception&) {
                    std::cerr << "Error: --k requires a valid integer.\n";
                    params.valid = false;
                }
            } else {
                std::cerr << "Error: --k requires an integer argument.\n";
                params.valid = false;
            }
        } else if (arg == "--theta") {
            if (i + 1 < argc) {
                try {
                    params.theta = std::stod(argv[++i]);
                    if (params.theta <= 0 || params.theta > 1) {
                        std::cerr << "Error: theta must be between 0 and 1.\n";
                        params.valid = false;
                    }
                } catch (const std::exception&) {
                    std::cerr << "Error: --theta requires a valid float value.\n";
                    params.valid = false;
                }
            } else {
                std::cerr << "Error: --theta requires a float argument.\n";
                params.valid = false;
            }
        } else if (arg == "--e") {
            if (i + 1 < argc) {
                try {
                    params.error_bound = std::stod(argv[++i]);
                    if (params.error_bound <= 0) {
                        std::cerr << "Error: error bound must be positive.\n";
                        params.valid = false;
                    }
                } catch (const std::exception&) {
                    std::cerr << "Error: --e requires a valid float value.\n";
                    params.valid = false;
                }
            } else {
                std::cerr << "Error: --e requires a float argument.\n";
                params.valid = false;
            }
        } else {
            std::cerr << "Warning: Unknown option: " << arg << "\n";
        }
    }
    
    // Validate required parameters

    if (params.k == 0) {
        std::cerr << "Error: user must provide k.\n";
        params.valid = false;
    }
    
    if (params.mode.empty()) {
        std::cerr << "Error: --mode is required.\n";
        params.valid = false;
    }
    
    if (params.input1.empty()) {
        std::cerr << "Error: --input1 is required.\n";
        params.valid = false;
    }
    else{
        if (params.mode != "kmer"){
            auto kspec1 = parse_sequence2kspec(params.input1, true, params.k);
            if (!kspec1) {
                std::cerr << "Failed to process the first input file.\n";
                params.valid = false;
            }
            else{
                ksp_set1 = *kspec1;
            }
        }
        else{
            auto kspec1 = parse_kmerfile2kspect(params.input1, params.k);
            if (!kspec1) {
                std::cerr << "Failed to process the first input file.\n";
                params.valid = false;
            }
            else{
                ksp_set1 = *kspec1;
            } 
        }
    }
    
    if (params.input2.empty()) {
        std::cerr << "Error: --input2 is required.\n";
        params.valid = false;
    }
    else{
        if (params.mode == "sequence"){
            auto kspec2 = parse_sequence2kspec(params.input2, false, params.k);
            if (!kspec2) {
                std::cerr << "Failed to process the second input file.\n";
                params.valid = false;
            }
            else{
                ksp_set2 = *kspec2;
            }
        }
        else{
            auto kspec2 = parse_kmerfile2kspect(params.input2, params.k);
            if (!kspec2) {
                std::cerr << "Failed to process the second input file.\n";
                params.valid = false;
            }
            else{
                ksp_set2 = *kspec2;
            } 
        }
    }
    
    // require dist file for kmer mode
    if (params.mode == "kmer") {
        if (params.dist.empty()){
            std::cerr << "Error: --dist is required for kmer mode.\n";
            params.valid = false;
        }
        else{
            auto occ_map = parse_csv(params.dist);
            if (!occ_map) {
                std::cerr << "Failed to process the occurrence distribution file.\n";
                params.valid = false;
            }
            else{
                coeffs_Map = *occ_map;
            }
        }
    }
    
    return params;
}


int main (int argc, char* argv[]){
    mpreal::set_default_prec(mpfr::digits2bits(digits));
    


    Parameters params = parse_arguments(argc, argv);
    
    if (!params.valid) {
        if (!params.help_printed){
            print_usage();
        }
        return 1;
    }

    int k = params.k;
    double epsilon = params.error_bound;
   
    

    std::cout << "Running Mutation_rate_estimator with parameters:\n";
    std::cout << "  Mode: " << params.mode << "\n";
    std::cout << "  Input1: " << params.input1 << "\n";
    std::cout << "  Input2: " << params.input2 << "\n";
    
    // dist file need, if kmer mode
    if (!params.dist.empty()) {
        std::cout << "  Dist: " << params.dist << "\n";
    }
    
    std::cout << "  k: " << params.k << "\n";
    std::cout << "  Theta: " << params.theta << "\n";
    std::cout << "  Error_bound for Newton's Method: " << epsilon << "\n";
    
    if (params.mode != "kmer"){
        for (const auto& pair : occurrence){
            int copy = pair.second.size();
            coeffs_Map[copy]++;
        }
    }
    
    int I_obs;
    double theta = params.theta;
    if (theta < 1.0){
        std::random_device rd;
        uint32_t seed = rd();
        set<string> kmer_set_sketch1 = kspectrum_sketch(ksp_set1, theta, seed);
        set<string> kmer_set_sketch2 = kspectrum_sketch(ksp_set2, theta, seed);
        int I_obs_sketch = round(static_cast<double>(intersect_size(kmer_set_sketch1, kmer_set_sketch2)) / (theta));
        mpreal q_hat_sketch = 0;
        mpreal r_hat_sketch = 0;
        if (I_obs_sketch < ksp_set1.size() && I_obs_sketch > 0 ){
            q_hat_sketch = Newton_Method (coeffs_Map, ksp_set1, k, epsilon, I_obs_sketch);
            r_hat_sketch = 1.0 - pow(1.0 - q_hat_sketch, 1.0/k);
        }
        else if (I_obs_sketch == 0){
            q_hat_sketch = 1.0;
            r_hat_sketch = 1.0;
        }

        cout << "r_sketch: " << r_hat_sketch << endl;
    }
    else{
        I_obs = intersect_size(ksp_set1, ksp_set2);
        mpreal q_hat = Newton_Method (coeffs_Map, ksp_set1, k, epsilon, I_obs);
        mpreal r_hat = 1.0 - pow(1.0 - q_hat, 1.0/k);
        cout << "r_strong: " << r_hat << endl;
    }



}