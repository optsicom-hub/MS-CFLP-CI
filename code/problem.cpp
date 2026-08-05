#include <stdio.h>
#include "problem.h"
#include <vector>
#include <iostream>

#include "solution.h"

using namespace std;

problem_t *create_problem(int num_S, int num_W)
{

    auto *problem = new problem_t;

    problem->num_S = num_S;
    problem->num_W = num_W;
    problem->capacity = new int[num_W];
    problem->fixedcost = new double[num_W];
    problem->goods = new int[num_S];
    problem->supplycost = new double*[num_S];
    problem->inc = new int*[num_S];
    problem->f_best_c = new int*[num_W];
    problem->c_best_f = new double*[num_S];
    problem->c_best_f_param = new int*[num_S];
    problem->c_best = new int*[num_S];
    problem->c_best_sc = new int*[num_S];
    problem->c_worst = new int*[num_S];

    for (int i = 0; i < num_S; i++)
    {
        problem->inc[i] = new int[num_S];
        problem->supplycost[i] = new double[num_W];
        problem->c_best_f[i] = new double[3];
        problem->c_best_f_param[i] = new int[2];
        problem->c_best[i] = new int[num_W];
        problem->c_best_sc[i] = new int[num_W];
        problem->c_worst[i] = new int[num_W];
    }
    for (int i = 0; i < num_W; i++)
    {
        problem->f_best_c[i] = new int[num_S];
    }


    return problem;
}

void free_problem(problem_t *problem)
{
    if (problem != NULL) {
        delete[] problem;
    }
}

int CalcCmaxMin(const problem_t *problem) {

        if (problem->num_S == 0 || problem->num_W == 0) {
            return 0;
        }

        int max_of_mins = std::numeric_limits<int>::min();

        for (int i = 0; i < problem->num_S; i++) {

            int min_cost = *std::min_element(
                problem->supplycost[i],
                problem->supplycost[i] + problem->num_W
            );

            max_of_mins = std::max(max_of_mins, min_cost);
        }

        return max_of_mins;
}


void printProblem(const problem_t *problem) {
    cout << "Instance" <<endl;

    cout << "S " << problem->num_S << endl;
    cout << "W " << problem->num_W << endl;

    cout << "Capacity: " << endl;
    for (int j = 0 ; j < problem->num_W; j++) {
        cout << problem->capacity[j] <<  " ";
    }
    cout <<endl;

    cout << "Fixed Cost: " << endl;
    for (int j = 0 ; j < problem->num_W; j++) {
        cout << problem->fixedcost[j] << " ";
    }
    cout <<endl;

    cout << "Goods: " << endl;
    for (int i = 0 ; i < problem->num_S; i++) {
        cout << problem->goods[i] << " ";
    }
    cout <<endl;

    cout << "Supply Cost: " << endl;
    for (int i = 0 ; i < problem->num_S; i++) {
        for (int j = 0; j < problem->num_W; j++) {
            cout << problem->supplycost[i][j] << " ";
        }
        cout << endl;
    }

    cout << "Incompatibilities: " << endl;
    for (int i = 0 ; i < problem->num_S; i++) {
        for (int j = 0; j < problem->num_S; j++) {
            cout << problem->inc[i][j] << " ";
        }
        cout << endl;
    }

    cout << "Num inc: " << problem->num_Inc <<endl;


}


