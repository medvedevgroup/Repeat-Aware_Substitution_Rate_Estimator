#ifndef SKETCH_H
#define SKETCH_H
#include <iostream>
#include <random>
#include <string>
#include <set>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <map>
#include <unordered_map>
#include "MurmurHash3.h"
using namespace std;

set<string> kspectrum_sketch(const set<string>& kmer_set, double theta, uint32_t seed);
// int intersect_size_sketch(const set<string>& set1, string t, int k, double theta, uint32_t seed);
#endif
