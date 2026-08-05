#ifndef PROBLEM_H
#define PROBLEM_H

using namespace std;

#define LITTLE_P

typedef struct {
    int num_S;
    int num_W;
    int *capacity;
    double *fixedcost;
    int *goods;
    double **supplycost;
    int **inc;
    int num_Inc;
    int **f_best_c;
    double **c_best_f;
    int **c_best_f_param;
    int **c_best;
    int **c_best_sc;
    int **c_worst;

} problem_t;

typedef struct {
    int num_S;
    int num_W;
    int *capacity;
    double *fixedcost;
    int *goods;
    double **supplycost;
    int **inc;
    int num_Inc;
    int **f_best_c;
    double **c_best_f;
    int **c_best_f_param;
    int **c_best;
    int **c_best_sc;
    int **c_worst;

} subproblem_t;

problem_t *create_problem(int, int);
void free_problem(problem_t *);
int CalcCmaxMin(const problem_t *);
void printProblem(const problem_t *problem);

#endif /* !PROBLEM_H */