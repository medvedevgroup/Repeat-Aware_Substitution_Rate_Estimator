#ifndef NEWTON_METHOD_H
#define NEWTON_METHOD_H
#include "mpreal.h"
#include <string>
#include <set>
#include <unordered_map>
using mpfr::mpreal;
using namespace std;
mpreal Newton_Method (const unordered_map<int, int>& coeffs_Map, const set<string>& ksp_set, int k, double epsilon, double I_obs);
std::set<std::string> kspectrum(std::string s, int k);
#endif