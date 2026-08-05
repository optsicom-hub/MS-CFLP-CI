#include "solution_exact.h"
#include "problem.h"

solutionex_t *create_solutionex(int num_S, int num_W)
{
    auto *solutionex = new solutionex_t;
    solutionex->cost = 0;
    solutionex->Dmax_sol = 0;
    solutionex->run_time = 0;
    solutionex->target_time = 0;
    solutionex->gap = -1.0;
    solutionex->X = new int*[num_S];
    solutionex->fopen = new int[num_W];
    solutionex->solcount = 0;

    for (int i = 0; i < num_S; i++)
    {
        solutionex->X[i] = new int[num_W];

        for (int j = 0; j < num_W; j++) {
            solutionex->X[i][j] = 0;
        }
    }

    for (int j = 0; j < num_W; j++) {
        solutionex->fopen[j] = 0;
    }
    return solutionex;
}

void copysolution_ex(const problem_t *problem, solutionex_t *solex_1, solutionex_t * solex_2) {
    solex_1->cost = solex_2->cost;
    solex_1->Dmax_sol = solex_2->Dmax_sol;
    solex_1->run_time = solex_2->run_time;
    solex_1->target_time = solex_2->target_time;
    solex_1->gap = solex_2->gap;
    int S = problem->num_S;
    int W = problem->num_W;
    solex_1->solcount = solex_2->solcount;

    for (int i = 0; i < S; i++) {
        for (int j = 0; j < W; j++) {
            solex_1->X[i][j] = solex_2->X[i][j];
        }
    }

    for (int j = 0; j < W; j++) {
        solex_1->fopen[j] = solex_2->fopen[j];
    }
}

void delete_solutionex(solutionex_t *solutionex, int num_S)
{
    if (solutionex == nullptr) {
        return;
    }

    if (solutionex->X != nullptr) {
        for (int i = 0; i < num_S; i++)
        {
            if (solutionex->X[i] != nullptr) {
                delete[] solutionex->X[i];
                solutionex->X[i] = nullptr;
            }
        }
    }

    if (solutionex->X != nullptr) {
        delete[] solutionex->X;
        solutionex->X = nullptr;
    }

    if (solutionex->fopen != nullptr) {
        delete[] solutionex->fopen;
        solutionex->fopen = nullptr;
    }

    delete solutionex;
}

void free_solution(solutionex_t *solutionex)
{
    if (solutionex != NULL) {
        delete []solutionex;
    }
}