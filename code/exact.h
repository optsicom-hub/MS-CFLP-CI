#ifndef EXACT_H_
#define EXACT_H_

#include "problem.h"
#include "solution.h"
#include "solution_exact.h"
#include <gurobi_c++.h>
#include <map>
#include <unordered_map>
#include <chrono>
#include <fstream>
#include <string>
#include <iomanip> // SEARCH_PROFILE

using namespace std;

struct PairHash {
    size_t operator()(const std::pair<int,int>& p) const {
        return std::hash<long long>()(((long long)p.first << 32) | (unsigned int)p.second);
    }
};

// Fires on every new MIP incumbent; records the first moment the objective
// meets the fo_literature target into the global t_ours_to_match.
class ConvergenceCallback : public GRBCallback {
public:
    std::chrono::steady_clock::time_point start_time;
    explicit ConvergenceCallback(const std::chrono::steady_clock::time_point& t0)
        : start_time(t0) {}
protected:
    void callback() override;
};

class IntervalConvergenceCallback : public GRBCallback {
public:
    std::ofstream& file;
    const std::chrono::steady_clock::time_point& global_start;
    double last_obj;
    std::string phase;
    double fo_literature;
    double& t_ours_to_match;

    IntervalConvergenceCallback(std::ofstream& f, const std::chrono::steady_clock::time_point& start, const std::string& p, double fo_lit, double& t_match)
        : file(f), global_start(start), last_obj(1e30), phase(p), fo_literature(fo_lit), t_ours_to_match(t_match) {}

protected:
    void callback() override {
        if (where == GRB_CB_MIPSOL) {
            double obj = getDoubleInfo(GRB_CB_MIPSOL_OBJ);
            if (obj < last_obj) {
                last_obj = obj;
                auto now = std::chrono::steady_clock::now();
                double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - global_start).count() / 1000.0;
                file << elapsed << " " << phase << " " << obj << "\n";
                file.flush();
            }
            if (obj <= fo_literature && t_ours_to_match < 0.0) {
                auto now = std::chrono::steady_clock::now();
                t_ours_to_match = std::chrono::duration_cast<std::chrono::milliseconds>(now - global_start).count() / 1000.0;
            }
        }
    }
};

// SEARCH_PROFILE: logs every improving incumbent of the constructive phase
// as "<time> CONST <obj>" (time in seconds since t0, 3 decimals).
class SearchProfileCallback : public GRBCallback { // SEARCH_PROFILE
public:
    std::ofstream& file;
    const std::chrono::steady_clock::time_point& t0;
    double best_obj;

    SearchProfileCallback(std::ofstream& f, const std::chrono::steady_clock::time_point& start)
        : file(f), t0(start), best_obj(1e30) {}

protected:
    void callback() override {
        if (where == GRB_CB_MIPSOL) {
            double obj = getDoubleInfo(GRB_CB_MIPSOL_OBJ);
            if (obj < best_obj) {
                best_obj = obj;
                double t = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
                file << std::fixed << std::setprecision(3) << t << " CONST " << obj << "\n";
                file.flush();
            }
        }
    }
};

// LOG_CONSTRUCTIVE: one line per run in resumen_fo.csv (defined in exact.cpp)
void WriteFoCsv(const std::string& name_inst, double z_const, double z_final,
    double t_best_const, double t_best_fo, double total_time,
    long long n_iter_fo, long long n_mejoras_fo);

class Exact {
public:
    Exact();
    ~Exact();

    void addConstraintsToGurobiModel2Opt(GRBModel& model,
                                     std::map<std::pair<int,int>, GRBVar>& X,
                                     GRBVar* Y,
                                     std::map<std::pair<int,int>, GRBVar>& Z,
                                     const problem_t *problem, int Dmax);

    void addConstraintsToGurobiModel3Opt(GRBModel& model,
                                 GRBVar** X,
                                 GRBVar** Z,
                                 GRBVar* Y,
                                 const problem_t *problem);

    void buildGurobiModel_ex(GRBModel& model,
                             std::unordered_map<std::pair<int,int>, GRBVar, PairHash>& X2,
                             std::unordered_map<std::pair<int,int>, GRBVar, PairHash>& Z2,
                             GRBVar* Y,
                             const problem_t* problem,
                             int DmaxParam,
                             std::map<int, GRBConstr>& cap_constrs,
                             std::map<int, GRBConstr>& goods_constrs);

    void extendGurobiModel_ex(GRBModel& model,
                              std::unordered_map<std::pair<int,int>, GRBVar, PairHash>& X2,
                              std::unordered_map<std::pair<int,int>, GRBVar, PairHash>& Z2,
                              GRBVar* Y,
                              const problem_t* problem,
                              int d,
                              std::map<int, GRBConstr>& cap_constrs,
                              std::map<int, GRBConstr>& goods_constrs);

    void runGurobi_ex(const problem_t*problem, solutionex_t *sol, string name_inst, int exact_mode, string ini_sol, double maxnode, double gap, const long long tlimit, int upperbound, int DmaxParam, const std::chrono::steady_clock::time_point& global_start, GRBEnv* external_env = nullptr);

    void RunGurobiEx_Interval(const problem_t *problem, solution_t *sol, solutionex_t *sol_ex, solutionex_t *best_ex, long long time_limit, const std::chrono::steady_clock::time_point& start_time,
        int DmaxIni, string name_inst, double upp_bound);

    void ConversionSolexToSol(const problem_t *problem, solution_t *sol, solutionex_t *sol_ex);

    vector<int> KnapsackClosedFac(const problem_t* problem, solution_t* solution,
                               const vector<int>& v_cust_subproblem,
                               int n_supp_close,
                               const vector<int>& v_supp_open);

};
#endif //EXACT_H


