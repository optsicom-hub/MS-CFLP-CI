#include "solution.h"
#include "problem.h"
#include "globals.h"
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <vector>
#include <fstream>

#include "exact.h"

using namespace std;

solution_t *create_solution(int num_S, int num_W)
{

    auto *solution = new solution_t;

    solution->num_S = num_S;
    solution->num_W = num_W;
    solution->assignment_f = new int*[num_S];
    solution->assignment_c = new int*[num_S];
    solution->facility_client = new int*[num_W];
    solution->fopen = new int[num_W]();
    solution->f_assigned = new int[num_W]();
    solution->g_assigned = new int[num_S]();
    solution->c_done = new int[num_S]();
    solution->cost = 0.0;
    solution->n_open = 0;
    solution->c_incompleted = 0;

    for (int i = 0; i < num_S; i++)
    {
        solution->assignment_f[i] = new int[1]();
        solution->assignment_c[i] = new int[1]();
    }

   for (int i = 0; i < num_W; i++) {
        solution->facility_client[i] = new int[1]();
   }

    return solution;
}

void init_solution(const problem_t *problem, solution_t *ini_sol) {
    int S = problem->num_S;
    int W = problem->num_W;
    ini_sol->num_S = S;
    ini_sol->num_W = W;

    for (int j = 0; j < W; j++)
    {
        ini_sol->fopen[j] = 0;
        ini_sol->f_assigned[j] = 0;
        ini_sol->facility_client[j][0] = 0;
    }

    for (int j = 0; j < S; j++)
    {
        ini_sol->c_done[j] = 0;
        ini_sol->g_assigned[j] = 0;

    }

    for (int i = 0; i < S; i++) {
        ini_sol->assignment_f[i][0] = 0;
        ini_sol->assignment_c[i][0] = 0;
    }

    ini_sol->cost =0;
    ini_sol->n_open = 0;
    ini_sol->c_incompleted = 0;

}

int GetsolutionCost(problem_t *problem, solution_t *solution, vector<int> v_cli, vector<int> v_supp, solutionex_t *exact_sol) {
    int fixed_cost = 0;
    int sc_cost = 0;

    for (int j = 0; j < problem->num_W; j++) {
        int f = v_supp[j];

        for (int i = 1; i < solution->facility_client[f][0]+1; i++) {
            int c = solution->facility_client[f][i];

            auto pos_cli = find(v_cli.begin(), v_cli.end(), c);
            int pos = -1;

            if (pos_cli != v_cli.end()) {
                pos = distance(v_cli.begin(), pos_cli);
            }

            for (int k = 1; k < solution->assignment_f[c][0]+1; k++) {
                if (solution->assignment_f[c][k] == f) {
                    exact_sol->X[i][j] = solution->assignment_c[c][k];
                    sc_cost += solution->assignment_c[c][k]*problem->supplycost[pos][j];
                    break;
                }
            }

        }

    }

    for (int j = 0; j < v_supp.size(); j++) {
        int f = v_supp[j];
        if (solution->fopen[f] == 1) {
            fixed_cost += problem->fixedcost[j];
        }
    }

    return fixed_cost + sc_cost;

}

void printsol (solution_t * sol) {
	cout << "Assignment: " << endl;
    for (int i = 0; i < sol->num_S; i++) {
        cout << "i = " << i << ": ";
        for (int j = 0; j < sol->assignment_f[i][0]+1;j++) {
            cout << sol->assignment_f[i][j] << " ";
        }
        cout << " | ";
        for (int j = 0; j < sol->assignment_f[i][0]+1;j++) {
            cout << sol->assignment_c[i][j] << " ";
        }
        cout << endl;
    }

    cout << "F open: " << endl;
    for (int i = 0; i < sol->num_W; i++) {
        cout << "i = " << i << ": ";
        cout << sol->fopen[i] << endl;
    }

    cout << "F assigned: " << endl;
    for (int i = 0; i < sol->num_W; i++) {
        cout << "i = " << i << ": ";
        cout << sol->f_assigned[i] << endl;
    }

    cout << "G assigned: " << endl;
    for (int i = 0; i < sol->num_S; i++) {
        cout << "i = " << i << ": ";
        cout << sol->g_assigned[i] << endl;
    }

    cout << "C done: " << endl;
    for (int i = 0; i < sol->num_S; i++) {
        cout << "i = " << i << ": ";
        cout << sol->c_done[i] << endl;
    }

    cout << "N open: " << sol->n_open << endl;
    cout << "C incompleted: " << sol->c_incompleted << endl;

    cout << "Facility client: " << endl;
    for (int i = 0; i < sol->num_W; i++) {
        cout << "i = " << i << ": ";
        for (int j = 0; j < sol->facility_client[i][0]+1; j++) {
            cout << sol->facility_client[i][j] << " ";
        }
        cout << endl;
    }
    cout << "/////////////////////////////"<<endl;
}

void free_solution(solution_t *solution)
{
    if (solution != NULL) {
        delete []solution;
    }
}

void delete_solution(solution_t *solution)
{
    if (solution == nullptr) return;

    for (int i = 0; i < solution->num_S; i++)
    {
        delete[] solution->assignment_f[i];
        delete[] solution->assignment_c[i];
    }

    for (int i = 0; i < solution->num_W; i++) {
        delete[] solution->facility_client[i];
    }

    delete[] solution->assignment_f;
    delete[] solution->assignment_c;

    delete[] solution->facility_client;


    delete[] solution->fopen;
    delete[] solution->f_assigned;
    delete[] solution->g_assigned;
    delete[] solution->c_done;

    delete solution;
}

void saveSolution(solution_t *solucion, const string &file_name, int S, int W) {

    ofstream ofs(file_name);
    if (!ofs) {
        cerr << "Error opening " << file_name << "." << endl;
        return;
    }

    ofs << "[" << endl;
    int c;
    int copy;
    int copy2;
    for (unsigned i = 0; i < S; ++i) {
        ofs << "  (";
        for (unsigned j = 0; j < W; ++j) {
            c = 0;
            int l = solucion->assignment_f[i][0];
            for (unsigned k = 1; k < l+1; k++) {
                if (solucion->assignment_f[i][k] == j) {
                    copy = solucion->assignment_c[i][k];
                    c++;
                }else {
                    copy2 = 0;
                }
            }
            if (c > 0) {
                ofs << copy;
            }else {
                ofs << copy2;
            }
            if (j < W - 1) {
                ofs << ", ";
            }
        }
        ofs << ")";
        ofs << endl;
    }
    ofs << "]" << endl;

    ofs.close();
}

void copy_solution(solution_t *sol1, const solution_t *sol2) {
    sol1->cost = sol2->cost;
    sol1->cost = sol2->cost;
    sol1->num_S = sol2->num_S;
    sol1->num_W = sol2->num_W;
    sol1->n_open = sol2->n_open;
    sol1->c_incompleted = sol2->c_incompleted;
    int S = sol2->num_S;
    int W = sol2->num_W;

    for (int i = 0; i < W; i++) {
        int longitud = sol2->facility_client[i][0];

        delete[] sol1->facility_client[i];
        sol1->facility_client[i] = new int[longitud + 1];

        std::memcpy(sol1->facility_client[i], sol2->facility_client[i], (longitud + 1) * sizeof(int));

        sol1->fopen[i] = sol2->fopen[i];
        sol1->f_assigned[i] = sol2->f_assigned[i];
    }

    for (int i = 0; i < S; i++) {
        int l = sol2->assignment_f[i][0];
        delete[] sol1->assignment_f[i];
        sol1->assignment_f[i] = new int[l + 1];
        std::memcpy(sol1->assignment_f[i], sol2->assignment_f[i], (l + 1) * sizeof(int));

        delete[] sol1->assignment_c[i];
        sol1->assignment_c[i] = new int[l + 1];
        std::memcpy(sol1->assignment_c[i], sol2->assignment_c[i], (l + 1) * sizeof(int));

        sol1->c_done[i] = sol2->c_done[i];
        sol1->g_assigned[i] = sol2->g_assigned[i];
    }

}

void DropCustomer(const problem_t *problem, solution_t *solution, int cliente, int facility_old) {
    int pos_facility;
    for (int i = 1; i < solution->assignment_f[cliente][0]+1; i++) {
        if (solution->assignment_f[cliente][i] == facility_old) {
            pos_facility = i;
            break;
        }
    }

    int q = solution->assignment_c[cliente][pos_facility];

    solution->g_assigned[cliente] = solution->g_assigned[cliente] - q;
    solution->f_assigned[facility_old] = solution->f_assigned[facility_old] - q;

    if (solution->c_done[cliente] == 2 && q > 0) {
        solution->c_incompleted--;
    }

    if (solution->g_assigned[cliente] == 0) {
        solution->c_done[cliente] = 0;
    }else {
        solution->c_done[cliente] = 1;
    }


    solution->cost = solution->cost - q*problem->supplycost[cliente][facility_old];
    int open = solution->fopen[facility_old];
    if (open == 1 && solution->f_assigned[facility_old] == 0) {
        solution->fopen[facility_old] = 0;
        solution->n_open--;
        solution->cost = solution->cost - problem->fixedcost[facility_old];
    }

    int l = solution->assignment_f[cliente][0];
    vector<int> assignment_c;
    vector<int> assignment_f;

    assignment_c.push_back(l-1);
    assignment_f.push_back(l-1);

    for (int i = 1; i < l+1; i++) {
        if (solution->assignment_f[cliente][i] != facility_old) {
            assignment_c.push_back(solution->assignment_c[cliente][i]);
            assignment_f.push_back(solution->assignment_f[cliente][i]);
        }
    }
    delete [] solution->assignment_c[cliente];
    solution->assignment_c[cliente] = new int [l];
    delete [] solution->assignment_f[cliente];
    solution->assignment_f[cliente] = new int [l];

    for (int i = 0; i < l; i++) {
        solution->assignment_c[cliente][i] = assignment_c[i];
        solution->assignment_f[cliente][i] = assignment_f[i];
    }

    int l2 = solution->facility_client[facility_old][0];
    vector <int> v_facilities;

    v_facilities.push_back(l2-1);
    for (int i = 1; i < l2+1; i++) {
        if (solution->facility_client[facility_old][i] != cliente) {
            int r = solution->facility_client[facility_old][i];
            v_facilities.push_back(r);
        }
    }
    delete [] solution->facility_client[facility_old];
    solution->facility_client[facility_old] = new int [l2];

    for (int j = 0; j < l2; j++) {
        solution->facility_client[facility_old][j] = v_facilities[j];
    }
}

void InsertCustomer_2q(const problem_t *problem, solution_t *sol, int cliente, int facility, int suma_valor ,int done) {
    int l, pos2;
    int f_open = sol->fopen[facility];

if ((sol->f_assigned[facility] + suma_valor <= problem->capacity[facility]) && (done < 2) && (suma_valor != 0)) {

    sol->cost = sol->cost + suma_valor*problem->supplycost[cliente][facility];
    sol->fopen[facility] = 1;
    sol->f_assigned[facility] = sol->f_assigned[facility] + suma_valor;
    sol->g_assigned[cliente] = sol->g_assigned[cliente] + suma_valor;

    if (sol->g_assigned[cliente] == problem->goods[cliente]) {
        sol->c_done[cliente] = 2;
        if (done < 2 && sol->c_done[cliente] == 2) {
            sol->c_incompleted++;
        }
    }else if (sol->g_assigned[cliente] > 0) {
        sol->c_done[cliente] = 1;
    }

    if (f_open == 0 && sol->fopen[facility] == 1) {
        sol->cost = sol->cost + problem->fixedcost[facility];
        sol->n_open++;
    }

    l = sol->assignment_f[cliente][0];
    std::vector<int> matrix_f;
    std::vector<int> matrix_c;

    for (int i = 0; i < l+1; i++) {
        matrix_f.push_back(sol->assignment_f[cliente][i]);
        matrix_c.push_back(sol->assignment_c[cliente][i]);
    }
    auto start = matrix_f.begin()+1;
    auto end = start + l;

    auto pos = std::find(start, end, facility);
    if (pos != end) {
        pos2 = distance(start, pos);
        matrix_c[pos2+1] = matrix_c[pos2+1] + suma_valor;

        delete [] sol->assignment_c[cliente];
        sol->assignment_c[cliente] = new int [l+1];

        for (int i = 0; i < l+1; i++) {
            sol->assignment_c[cliente][i] = matrix_c[i];
        }
    }else {
        matrix_f.push_back(facility);
        matrix_c.push_back(suma_valor);

        delete [] sol->assignment_f[cliente];
        delete [] sol->assignment_c[cliente];
        sol->assignment_f[cliente] = new int [l+2];
        sol->assignment_c[cliente] = new int [l+2];

        sol->assignment_f[cliente][0] = matrix_f[0]+1;
        sol->assignment_c[cliente][0] = matrix_c[0]+1;
        for (int i = 1; i < l+2; i++) {
            sol->assignment_f[cliente][i] = matrix_f[i];
            sol->assignment_c[cliente][i] = matrix_c[i];
        }
    }

    int longitud = sol->facility_client[facility][0];
    std::vector<int> v_temporal(longitud+2);

    for (int i = 0; i < longitud+1; i++) {
            v_temporal[i] = sol->facility_client[facility][i];
    }
    v_temporal[longitud+1] = -1;


    if (longitud != 0) {
        auto start = v_temporal.begin()+1;
        auto end = v_temporal.end();
        pos = find(start, end, cliente);
        if (pos == end) {
            v_temporal[0] = v_temporal[0]+1;

            v_temporal[longitud+1] = cliente;

            delete [] sol->facility_client[facility];

            sol->facility_client[facility] = new int [longitud+2];

            for (int i = 0; i < longitud+2; i++) {
                sol->facility_client[facility][i] = v_temporal[i];
            }


        }
    }else {

        v_temporal[0] = v_temporal[0]+1;

        v_temporal[longitud+1] = cliente;

        delete [] sol->facility_client[facility];

        sol->facility_client[facility] = new int [longitud+2];

        for (int i = 0; i < longitud+2; i++) {
            sol->facility_client[facility][i] = v_temporal[i];
        }
    }

    }

}

vector<int> ClosestFac_norm(const problem_t*problem, solution_t *solution, vector<int> v_cli) {

    int n = v_cli.size();
    int cli, q;

    vector<vector<int>> SumCostFac;
    vector<int> vector_q;

    for (int j = 0; j < problem->num_W; j++) {
        int suma = 0;

        if (solution->fopen[j]==0) {
            for (int i = 0; i < n; i++) {
                int cliente = v_cli[i];
                suma += problem->supplycost[cliente][j];
            }

            suma = ceil(1.0*suma/sc_norm + (1.0*problem->fixedcost[j]/fixed_norm)/(1.0*problem->capacity[j]/capacity_norm));
            SumCostFac.push_back({j,suma});

        }

    }

    std::sort(SumCostFac.begin(), SumCostFac.end(), [](const std::vector<int>& a, const std::vector<int>& b) {
        return a[1] < b[1];
    });

    vector<int> ReturnFac;

    for (int i = 0; i < SumCostFac.size(); i++) {
        ReturnFac.push_back(SumCostFac[i][0]);
    }

    return ReturnFac;
}

void GenerateSubproblem(const problem_t *problem, problem_t *subprob, solution_t *sol, int ncli, int nsupp, vector<int> v_cli, vector<int> v_supp) {
    subprob->num_S = ncli;
    subprob->num_W = nsupp;

    for (int j = 0; j < nsupp; j++) {
        int f = v_supp[j];
        subprob->capacity[j] = problem->capacity[f];
        subprob->fixedcost[j] = problem->fixedcost[f];
    }

    for (int i = 0; i < ncli; i++) {
        int c = v_cli[i];

        int suma = 0;

        for (int j = 1; j < sol->assignment_f[c][0]+1; j++) {
            int f = sol->assignment_f[c][j];

            auto it = find(v_supp.begin(), v_supp.end(), f);

            if (it != v_supp.end()) {
                suma += sol->assignment_c[c][j];
            }


        }

        subprob->goods[i] = suma;

    }

    for (int i = 0; i < ncli; i++) {
        int c = v_cli[i];
        for (int j = 0; j < nsupp; j++) {
            int f = v_supp[j];
            subprob->supplycost[i][j] = problem->supplycost[c][f];
        }
    }

    int cont = 0;
    for (int i1 = 0; i1 < ncli; i1++) {
        int c1 = v_cli[i1];
        for (int i2 = 0; i2 < ncli; i2++) {
            int c2 = v_cli[i2];
            subprob->inc[i1][i2] = problem->inc[c1][c2];

            if (i1 > i2 && subprob->inc[i1][i2]>0) {
                cont++;
            }
        }
    }

    subprob->num_Inc = cont;

    for (int i = 0; i < subprob->num_S; i++){

        for (int j = 0; j < subprob->num_W; j++){
            subprob->c_best[i][j] = 0;
            subprob->c_best_sc[i][j] = 0;
            subprob->c_worst[i][j] = 0;
        }

        subprob->c_best_f[i][0] = 0;
        subprob->c_best_f[i][1] = 0;
        subprob->c_best_f[i][2] = 0;

        subprob->c_best_f_param[i][0] = 0;
        subprob->c_best_f_param[i][1] = 0;

    }
    for (int j = 0; j < subprob->num_W; j++)
    {
        for (int i = 0; i < subprob->num_S; i++) {
            subprob->f_best_c[j][i] = 0;
        }

    }

}

vector<pair<int,double>> ClosestFac_norm2(const problem_t* problem, solution_t* solution, vector<int> v_cli) {

    int n = v_cli.size();

    vector<pair<int,double>> ReturnFac;

    for (int j = 0; j < problem->num_W; j++) {
        if (solution->fopen[j] == 0) {
            double suma = 0.0;
            for (int i = 0; i < n; i++) {
                int cliente = v_cli[i];
                suma += problem->supplycost[cliente][j];
            }

            suma = suma/sc_norm + (1.0*problem->fixedcost[j]/fixed_norm)/(1.0*problem->capacity[j]/capacity_norm);
            ReturnFac.push_back({j, suma});
        }
    }

    return ReturnFac;
}
