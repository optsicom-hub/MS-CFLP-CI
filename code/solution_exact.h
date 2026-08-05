#ifndef SOLUTION_EXACT_H
#define SOLUTION_EXACT_H


#include <climits>
#include <set>
#include "problem.h"

using namespace std;

#define LITTLE_P

typedef struct {
    int cost;
    int Dmax_sol;
    long long run_time;
    long long target_time;
    //int sol_count;
    int **X;
    double gap;
    int solcount;

    int *fopen;
} solutionex_t;



solutionex_t *create_solutionex(int, int);
void copysolution_ex(const problem_t *, solutionex_t *, solutionex_t *);
void delete_solutionex(solutionex_t *solutionex, int num_S);
void free_solution(solutionex_t *);

#endif //SOLUTION_EXACT_H


