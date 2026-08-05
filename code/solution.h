#ifndef SOLUTION_H
#define SOLUTION_H
#include <sstream>
#include "problem.h"
#include <vector>

#include "solution_exact.h"

using namespace std;

#define LITTLE_P

typedef struct {
    int num_S;
    int num_W;
    double cost;
    int **assignment_f;
    int **assignment_c;
    int *fopen;
    int *f_assigned;
    int *g_assigned;
    int *c_done;
    int n_open;
    int c_incompleted;
    int **facility_client;
    int **discards;
} solution_t;



solution_t *create_solution(int, int);
void printsol (solution_t *);
void free_solution(solution_t *);
void init_solution(const problem_t *problem, solution_t *solution);
int GetsolutionCost(problem_t *, solution_t *, vector<int>, vector<int>, solutionex_t *);
void delete_solution(solution_t *solution);
void saveSolution(solution_t *solucion, const string &file_name, int S, int W);
void copy_solution(solution_t *sol1, const solution_t *sol2);
void DropCustomer(const problem_t *problem, solution_t *solution, int cliente, int facility_old);
void InsertCustomer_2q(const problem_t *problem, solution_t *sol, int cliente, int facility, int suma_valor ,int done);
vector<int> ClosestFac_norm(const problem_t*problem, solution_t *solution, vector<int> v_cli);
void GenerateSubproblem(const problem_t *problem, problem_t *subprob, solution_t *sol, int ncli, int nsupp, vector<int> v_cli, vector<int> v_supp);
#endif
