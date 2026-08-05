#ifndef GLOBALS_H
#define GLOBALS_H

#include <random>
using namespace std;

extern string name_inst;
extern int seed;
extern std::mt19937 rng;
extern double pct_open;
extern double pct_close;
extern int subproblem_size;
extern int sc_norm;
extern int fixed_norm;
extern int capacity_norm;
extern int strategy;
extern double fo_literature;
extern double t_ours_to_match;
//extern string results_path;
extern int divisor_time;
#endif // GLOBALS_H