#include "globals.h"
#include <chrono>
using namespace std;

//int seed = std::chrono::system_clock::now().time_since_epoch().count();
string name_inst = " ";
int seed = 12345;
std::mt19937 rng(seed);
double pct_open = 0.0;
double pct_close = 0.0;
int subproblem_size = 0;
int sc_norm = 0;
int fixed_norm = 0;
int capacity_norm = 0;
int strategy = 2;
double t_ours_to_match = -1.0;
//string results_path = "convergence_results.txt";
int divisor_time = 1;
