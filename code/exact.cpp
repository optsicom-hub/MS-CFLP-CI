#include <algorithm>
#include <iostream>
#include <random>
#include <vector>
#include <string>
#include <gurobi_c++.h>
#include <fstream>
#include <cmath>
#include <map>
#include <regex>
#include <chrono>
#include <numeric>
#include <filesystem>
#include <iomanip>
#include <sstream> // SEARCH_PROFILE
#include <functional> // LOG_CONSTRUCTIVE
#include "exact.h"
#include "solution_exact.h"
#include "problem.h"
#include "solution.h"
#include "globals.h"
using namespace std;
Exact::Exact()
{
}

Exact::~Exact(void)
{
}


void Exact::addConstraintsToGurobiModel3Opt(GRBModel& model,
                                 GRBVar** X,
                                 GRBVar** Z,
                                 GRBVar* Y,
                                 const problem_t *problem)
{
    try {

        // Capacity
        for (int j = 0; j < problem->num_W; ++j) {
            GRBLinExpr lhs = 0;

            std::vector<GRBVar> vars;
            std::vector<double> coeffs;

            for (int i = 0; i < problem->num_S; ++i) {
                vars.push_back(X[i][j]);
                coeffs.push_back(1.0);
            }

            if (!vars.empty()) {
                lhs.addTerms(&coeffs[0], &vars[0], vars.size());
            }

            lhs -= problem->capacity[j] * Y[j];

            model.addConstr(lhs <= 0);
        }

        // Goods/Demand
        for (int i = 0; i < problem->num_S; ++i) {
            GRBLinExpr lhs = 0;

            std::vector<GRBVar> vars;
            std::vector<double> coeffs;

            for (int j = 0; j < problem->num_W; ++j) {
                vars.push_back(X[i][j]);
                coeffs.push_back(1.0);
            }

            if (!vars.empty()) {
                lhs.addTerms(&coeffs[0], &vars[0], vars.size());
            }

            model.addConstr(lhs == problem->goods[i]);
        }

        // Upper bound
        for (int i = 0; i < problem->num_S; ++i) {
            for (int j = 0; j < problem->num_W; ++j) {
                GRBLinExpr lhs = X[i][j] - problem->goods[i] * Z[i][j];
                model.addConstr(lhs <= 0);
            }
        }

    	// Incompatibilities
    	std::vector<std::pair<int, int>> incompatible_pairs;
    	for (int i1 = 0; i1 < problem->num_S; ++i1) {
    		for (int i2 = i1 + 1; i2 < problem->num_S; ++i2) {
    			if (problem->inc[i1][i2] == 1) {
    				incompatible_pairs.push_back({i1, i2});
    			}
    		}
    	}

    	for (int j = 0; j < problem->num_W; ++j) {
    		std::vector<GRBLinExpr> lhs_exprs;
    		lhs_exprs.reserve(incompatible_pairs.size());

    		for (const auto& pair : incompatible_pairs) {
    			GRBLinExpr lhs = Z[pair.first][j] + Z[pair.second][j];
    			lhs_exprs.push_back(lhs);
    		}

    		for (const auto& lhs : lhs_exprs) {
    			model.addConstr(lhs <= Y[j]);
    		}
    	}


    } catch (GRBException& e) {
        std::cerr << "Error Gurobi: " << e.getMessage() << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}

void Exact::runGurobi_ex(const problem_t*problem, solutionex_t *sol, string name_inst, int exact_mode, string ini_sol, double maxnode, double gap, const long long tlimit, int upperbound, int DmaxParam, const std::chrono::steady_clock::time_point& global_start, GRBEnv* external_env) {
    try {
    	auto start_time = std::chrono::steady_clock::now();
    	auto current_time = std::chrono::steady_clock::now();
    	auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count();

		// Set exact_mode = 1

    	vector<vector<int>> supplycost;
        vector<int> fixedcost;
        vector<int> capacity;
        vector<int> goods;
        vector<vector<int>> incomp;

		int DmaxReturn = 0;

    	vector<vector<int>> ini_x = vector(problem->num_S, std::vector<int>(problem->num_W, 0));
    	vector<int> ini_y = vector(problem->num_W, 0);

    	GRBEnv local_env(true);
    	GRBEnv* env_ptr;
    	if (external_env != nullptr) {
    		env_ptr = external_env;
    	} else {
    		if (exact_mode == 2) {
    			local_env.set(GRB_IntParam_OutputFlag, 1);
    		} else {
    			local_env.set(GRB_IntParam_OutputFlag, 0);
    		}
    		//local_env.set(GRB_IntParam_OutputFlag, 1);
    		local_env.set(GRB_IntParam_Threads, 1);
    		local_env.set(GRB_IntParam_Seed, 12345);
    		local_env.start();
    		env_ptr = &local_env;
    	}

        GRBModel model = GRBModel(*env_ptr);
        model.set(GRB_StringAttr_ModelName, "MyModel");



    	GRBVar** X1;
    	GRBVar** Z1;
    	std::map<std::pair<int,int>, GRBVar> X2;
    	std::map<std::pair<int,int>, GRBVar> Z2;

    	// Create arrays for variable attributes
    	int num_Y = problem->num_W;
    	double* lb_Y = new double[num_Y];
    	double* ub_Y = new double[num_Y];
    	double* obj_Y = new double[num_Y];
    	char* vtype_Y = new char[num_Y];

    	for (int j = 0; j < num_Y; ++j) {
    		lb_Y[j] = 0.0;
    		ub_Y[j] = 1.0;
    		obj_Y[j] = 0.0;
    		vtype_Y[j] = GRB_BINARY;
    	}

    	// Single batch call to create ALL variables at once
    	GRBVar* Y = model.addVars(lb_Y, ub_Y, obj_Y, vtype_Y, nullptr, num_Y);

    	for (int j = 0; j < num_Y; ++j) {
    		Y[j].set(GRB_DoubleAttr_Start, ini_y[j]);
    	}

    	delete[] lb_Y;
    	delete[] ub_Y;
    	delete[] obj_Y;
    	delete[] vtype_Y;


    	if (exact_mode == 1) {

    		int total_vars = problem->num_S * problem->num_W;
    		double* lb_X1 = new double[total_vars];
    		double* ub_X1 = new double[total_vars];
    		char* vtype_X1 = new char[total_vars];


    		int idx = 0;
    		for (int i = 0; i < problem->num_S; ++i) {
    			for (int j = 0; j < problem->num_W; ++j) {
    				lb_X1[idx] = 0.0;
    				ub_X1[idx] = problem->goods[i];
    				vtype_X1[idx] = GRB_INTEGER;
    				idx++;
    			}
    		}

    		GRBVar* X1_flat = model.addVars(lb_X1, ub_X1, nullptr, vtype_X1, nullptr, total_vars);

    		double* lb_Z1 = new double[total_vars];
    		double* ub_Z1 = new double[total_vars];
    		char* vtype_Z1 = new char[total_vars];

    		for (int k = 0; k < total_vars; ++k) {
    			lb_Z1[k] = 0.0;
    			ub_Z1[k] = 1.0;
    			vtype_Z1[k] = GRB_BINARY;
    		}

    		GRBVar* Z1_flat = model.addVars(lb_Z1, ub_Z1, nullptr, vtype_Z1, nullptr, total_vars);

    		X1 = new GRBVar*[problem->num_S];
    		Z1 = new GRBVar*[problem->num_S];

    		idx = 0;
    		for (int i = 0; i < problem->num_S; ++i) {
    			X1[i] = new GRBVar[problem->num_W];
    			Z1[i] = new GRBVar[problem->num_W];

    			for (int j = 0; j < problem->num_W; ++j) {
    				X1[i][j] = X1_flat[idx];
    				Z1[i][j] = Z1_flat[idx];
    				idx++;
    			}
    		}

    		delete[] lb_X1;
    		delete[] ub_X1;
    		delete[] vtype_X1;
    		delete[] lb_Z1;
    		delete[] ub_Z1;
    		delete[] vtype_Z1;
    		delete[] X1_flat;
    		delete[] Z1_flat;

    		GRBLinExpr obj = 0;
    		std::vector<GRBVar> all_vars;
    		std::vector<double> all_coeffs;

    		for (int i = 0; i < problem->num_S; ++i) {
    			for (int j = 0; j < problem->num_W; ++j) {
    				all_vars.push_back(X1[i][j]);
    				all_coeffs.push_back(static_cast<double>(problem->supplycost[i][j]));
    			}
    		}

    		for (int j = 0; j < problem->num_W; ++j) {
    			all_vars.push_back(Y[j]);
    			all_coeffs.push_back(static_cast<double>(problem->fixedcost[j]));
    		}

    		if (!all_vars.empty()) {
    			obj.addTerms(&all_coeffs[0], &all_vars[0], all_vars.size());
    		}

    		model.setObjective(obj, GRB_MINIMIZE);

    		addConstraintsToGurobiModel3Opt(model, X1, Z1, Y, problem);

    	}

    	current_time = std::chrono::steady_clock::now();
    	elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count();

    	model.update();

    	model.set(GRB_DoubleParam_MIPGap, gap);

    	if (maxnode != -1.0) {
    		model.set(GRB_DoubleParam_NodeLimit, maxnode); // we don't use it right now
    	}

    	std::string baseFileName;

    	if (upperbound != -1.0) {
    		model.set(GRB_DoubleParam_Cutoff, upperbound);
    	}
    	size_t lastDotPos = name_inst.rfind('.');
    	if (lastDotPos != std::string::npos && lastDotPos < name_inst.length() - 1) {
    		baseFileName = name_inst.substr(0, lastDotPos);
    	} else {
    		baseFileName = name_inst;
    	}

    	model.set(GRB_StringParam_LogFile, "gurobi_output_" + baseFileName + ".log");
    	if (exact_mode == 2) {
    		model.set(GRB_StringParam_LogFile, "gurobi_output_" + baseFileName + ".log"); // Name of the log file
    	}


    	model.set(GRB_IntParam_LogToConsole, 1);
    	current_time = std::chrono::steady_clock::now();
    	elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count();

    	sol->lower_bound = -1.0; // EXACT_FULL
    	sol->mip_gap = -1.0;     // EXACT_FULL

    	if (tlimit != -1 && tlimit-elapsed_time/1000.0 > 0) {
    		model.set(GRB_DoubleParam_TimeLimit, tlimit-elapsed_time/1000.0);

    		/*if (exact_mode == 2) {
    			model.set(GRB_IntParam_MIPFocus, 1);
    			model.set(GRB_IntParam_Presolve, 2);
    		}*/

    		current_time = std::chrono::steady_clock::now();
    		elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count();
    		model.optimize();

    		// EXACT_FULL: lower bound always read; MIP gap only if there is a solution (-1 otherwise)
    		{
    			try { sol->lower_bound = model.get(GRB_DoubleAttr_ObjBound); } catch (...) { sol->lower_bound = -1.0; } // EXACT_FULL
    			int solcount_ef = 0;
    			try { solcount_ef = model.get(GRB_IntAttr_SolCount); } catch (...) { solcount_ef = 0; } // EXACT_FULL
    			if (solcount_ef > 0) {
    				try { sol->mip_gap = model.get(GRB_DoubleAttr_MIPGap); } catch (...) { sol->mip_gap = -1.0; }         // EXACT_FULL
    			}
    		}

    		vector<double> ret_vec;
    	}

    	if (model.get(GRB_IntAttr_SolCount) == 0) {
    		sol->cost = -1;

    		sol->Dmax_sol = -1;
    		sol->target_time = 0;
			sol->gap = -1.0;

    		for (int i = 0; i < problem->num_S; i++) {
    			for (int j = 0; j < problem->num_W; j++ ) {
    				sol->X[i][j] = 0;
    			}
    		}

    		for (int j = 0; j < problem->num_W; j++ ) {
    			sol->fopen[j] = 0;
    		}
    	}else {
    		sol->cost = model.get(GRB_DoubleAttr_ObjVal);
    		if (exact_mode == 1) {
    			for (int i = 0; i < problem->num_S; i++) {
    				for (int j = 0; j < problem->num_W; j++ ) {
						sol->X[i][j] = abs(X1[i][j].get(GRB_DoubleAttr_X));
    				}
    			}
    		}else {
    			for (int i = 0; i < problem->num_S; i++) {
    				for (int j = 0; j < problem->num_W; j++ ) {
    					if (X2.count({i, j})) {
    						sol->X[i][j] = abs(X2[{i,j}].get(GRB_DoubleAttr_X));
    					}else {
    						sol->X[i][j] = 0;
    					}

    				}
    			}
    		}
    		sol->Dmax_sol = DmaxReturn;

    		for (int j = 0; j < problem->num_W; j++) {
    			for (int i = 0; i < problem->num_S; i++) {
    				if (sol->X[i][j]>0) {
    					sol->fopen[j] = 1;
    					break;
    				}
    			}
    		}


    		int status = model.get(GRB_IntAttr_Status);

    		if (status == GRB_OPTIMAL || status == GRB_SUBOPTIMAL) {
    			sol->gap = model.get(GRB_DoubleAttr_MIPGap);
    		} else {
    			sol->gap = -1.0;
    		}


    	}

    	sol->solcount = model.get(GRB_IntAttr_SolCount);

    } catch (GRBException e) {
        std::cerr << "Error Gurobi: " << e.getMessage() << std::endl;
        std::cerr << "Error Code: " << e.getErrorCode() << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

}

void Exact::ConversionSolexToSol(const problem_t *problem, solution_t *sol, solutionex_t *sol_ex) {
	sol->cost = sol_ex->cost;

	sol->num_W = problem->num_W;
	sol->num_S = problem->num_S;

	for (int j = 0; j < problem->num_W; j++) {
		sol->fopen[j] = sol_ex->fopen[j];
		sol->f_assigned[j] = 0;
	}

	sol->c_incompleted = 0;
	for (int i = 0; i < problem->num_S; i++) {

		int contador = 0;

		for (int k = 0; k < problem->num_W; k++) {
			if(sol_ex->X[i][k]>0) {
				contador++;
			}
		}

		if (sol->assignment_f[i] != nullptr) {
			delete[] sol->assignment_f[i];
			sol->assignment_f[i] = new int[contador + 1];
			sol->assignment_f[i][0] = 0;
		}

		if (sol->assignment_c[i] != nullptr) {
			delete[] sol->assignment_c[i];
			sol->assignment_c[i] = new int[contador + 1];
			sol->assignment_c[i][0] = 0;
		}


		sol->g_assigned[i] = 0;
		for (int j = 0; j < problem->num_W; j++) {
			if(sol_ex->X[i][j]>0) {
				int q = sol_ex->X[i][j];
				sol->assignment_f[i][0]++;
				sol->assignment_c[i][0]++;

				int pos = sol->assignment_f[i][0];
				sol->assignment_f[i][pos] = j;
				sol->assignment_c[i][pos] = q;

				sol->g_assigned[i] += q;
				sol->f_assigned[j] += q;
			}
		}

		if (sol->g_assigned[i] == 0) {
			sol->c_done[i] = 0;
		}else if(sol->g_assigned[i] > 0 && sol->g_assigned[i] < problem->goods[i]) {
			sol->c_done[i] = 1;
		}else {
			sol->c_done[i] = 2;
			sol->c_incompleted++;
		}


	}
	sol->n_open = 0;
	for (int j = 0; j < problem->num_W; j++) {
		if (sol->f_assigned[j] > 0) {
			sol->n_open++;
		}

		int longitud = 0;
		vector<int> v_cli;
		for (int i = 0; i < problem->num_S; i++) {
			if (sol_ex->X[i][j]>0) {
				longitud++;
				v_cli.push_back(i);
			}
		}

		if (sol->facility_client[j] != nullptr) {
			delete[] sol->facility_client[j];
			sol->facility_client[j] = new int[longitud + 1];

			sol->facility_client[j][0] = longitud;
			for (int k = 1; k < longitud + 1; k++) {
				sol->facility_client[j][k] = v_cli[k-1];
			}
		}


	}


}

void Exact::buildGurobiModel_ex(GRBModel& model,
                                std::unordered_map<std::pair<int,int>, GRBVar, PairHash>& X2,
                                std::unordered_map<std::pair<int,int>, GRBVar, PairHash>& Z2,
                                GRBVar* Y,
                                const problem_t* problem,
                                int DmaxParam,
                                std::map<int, GRBConstr>& cap_constrs,
                                std::map<int, GRBConstr>& goods_constrs) {
	std::vector<std::pair<int, int>> valid_indices;
	std::vector<double> x_upper_bounds;

	for (int i = 0; i < problem->num_S; ++i) {
		for (int j = 0; j < problem->num_W; ++j) {
			if (problem->supplycost[i][j] <= DmaxParam) {
				valid_indices.push_back({i, j});
				x_upper_bounds.push_back(problem->goods[i]);
			}
		}
	}

	int num_valid = valid_indices.size();

	double* lb_X2 = new double[num_valid];
	double* ub_X2 = new double[num_valid];
	char* vtype_X2 = new char[num_valid];

	for (int k = 0; k < num_valid; ++k) {
		lb_X2[k] = 0.0;
		ub_X2[k] = x_upper_bounds[k];
		vtype_X2[k] = GRB_INTEGER;
	}

	GRBVar* X2_array = model.addVars(lb_X2, ub_X2, nullptr, vtype_X2, nullptr, num_valid);

	double* lb_Z2 = new double[num_valid];
	double* ub_Z2 = new double[num_valid];
	char* vtype_Z2 = new char[num_valid];

	for (int k = 0; k < num_valid; ++k) {
		lb_Z2[k] = 0.0;
		ub_Z2[k] = 1.0;
		vtype_Z2[k] = GRB_BINARY;
	}

	GRBVar* Z2_array = model.addVars(lb_Z2, ub_Z2, nullptr, vtype_Z2, nullptr, num_valid);

	for (int k = 0; k < num_valid; ++k) {
		X2[valid_indices[k]] = X2_array[k];
		Z2[valid_indices[k]] = Z2_array[k];
	}

	delete[] lb_X2;
	delete[] ub_X2;
	delete[] vtype_X2;
	delete[] lb_Z2;
	delete[] ub_Z2;
	delete[] vtype_Z2;
	delete[] X2_array;
	delete[] Z2_array;

	GRBLinExpr obj = 0;
	std::vector<GRBVar> all_vars;
	std::vector<double> all_coeffs;

	for (const auto& entry : X2) {
		all_vars.push_back(entry.second);
		all_coeffs.push_back(static_cast<double>(problem->supplycost[entry.first.first][entry.first.second]));
	}

	for (int j = 0; j < problem->num_W; ++j) {
		all_vars.push_back(Y[j]);
		all_coeffs.push_back(static_cast<double>(problem->fixedcost[j]));
	}

	if (!all_vars.empty()) {
		obj.addTerms(&all_coeffs[0], &all_vars[0], all_vars.size());
	}

	model.setObjective(obj, GRB_MINIMIZE);

	// Capacity constraints (stored for incremental extension)
	std::vector<std::vector<GRBVar>> vars_by_warehouse(problem->num_W);
	for (const auto& entry : X2)
		vars_by_warehouse[entry.first.second].push_back(entry.second);
	for (int j = 0; j < problem->num_W; ++j) {
		if (!vars_by_warehouse[j].empty()) {
			GRBLinExpr lhs = 0;
			std::vector<double> coeffs(vars_by_warehouse[j].size(), 1.0);
			lhs.addTerms(&coeffs[0], &vars_by_warehouse[j][0], vars_by_warehouse[j].size());
			lhs -= problem->capacity[j] * Y[j];
			cap_constrs[j] = model.addConstr(lhs <= 0);
		}
	}

	// Goods/demand constraints (stored for incremental extension)
	std::vector<std::vector<GRBVar>> vars_by_supplier(problem->num_S);
	for (const auto& entry : X2)
		vars_by_supplier[entry.first.first].push_back(entry.second);
	for (int i = 0; i < problem->num_S; ++i) {
		if (!vars_by_supplier[i].empty()) {
			GRBLinExpr lhs = 0;
			std::vector<double> coeffs(vars_by_supplier[i].size(), 1.0);
			lhs.addTerms(&coeffs[0], &vars_by_supplier[i][0], vars_by_supplier[i].size());
			goods_constrs[i] = model.addConstr(lhs == problem->goods[i]);
		}
	}

	// Upper bound constraints
	for (const auto& entry : X2) {
		int i = entry.first.first;
		int j = entry.first.second;
		auto it_Z = Z2.find({i, j});
		if (it_Z != Z2.end())
			model.addConstr(entry.second - problem->goods[i] * it_Z->second <= 0);
	}

	// Incompatibility constraints
	std::vector<std::pair<int,int>> incompatible_pairs;
	for (int i1 = 0; i1 < problem->num_S; ++i1)
		for (int i2 = i1+1; i2 < problem->num_S; ++i2)
			if (problem->inc[i1][i2] == 1)
				incompatible_pairs.push_back({i1, i2});

	std::vector<std::vector<int>> incomp_neighbors(problem->num_S);
	for (const auto& p : incompatible_pairs) {
		incomp_neighbors[p.first].push_back(p.second);
		incomp_neighbors[p.second].push_back(p.first);
	}

	std::vector<std::vector<std::pair<int,GRBVar>>> Z_by_warehouse(problem->num_W);
	for (const auto& entry : Z2)
		Z_by_warehouse[entry.first.second].push_back({entry.first.first, entry.second});

	for (int j = 0; j < problem->num_W; ++j) {
		const auto& suppliers_at_j = Z_by_warehouse[j];
		if (suppliers_at_j.empty()) continue;
		std::map<int,GRBVar> supplier_to_var;
		for (const auto& sv : suppliers_at_j) supplier_to_var[sv.first] = sv.second;
		for (size_t idx1 = 0; idx1 < suppliers_at_j.size(); ++idx1) {
			int i1 = suppliers_at_j[idx1].first;
			GRBVar z1 = suppliers_at_j[idx1].second;
			for (int i2 : incomp_neighbors[i1]) {
				if (i2 > i1) {
					auto it = supplier_to_var.find(i2);
					if (it != supplier_to_var.end())
						model.addConstr(z1 + it->second <= Y[j]);
				}
			}
		}
	}
}

void Exact::extendGurobiModel_ex(GRBModel& model,
                                  std::unordered_map<std::pair<int,int>, GRBVar, PairHash>& X2,
                                  std::unordered_map<std::pair<int,int>, GRBVar, PairHash>& Z2,
                                  GRBVar* Y,
                                  const problem_t* problem,
                                  int d,
                                  std::map<int, GRBConstr>& cap_constrs,
                                  std::map<int, GRBConstr>& goods_constrs) {
	// Precompute incompatibility neighbors
	std::vector<std::vector<int>> incomp_neighbors(problem->num_S);
	for (int i1 = 0; i1 < problem->num_S; ++i1)
		for (int i2 = i1+1; i2 < problem->num_S; ++i2)
			if (problem->inc[i1][i2] == 1) {
				incomp_neighbors[i1].push_back(i2);
				incomp_neighbors[i2].push_back(i1);
			}

	for (int i = 0; i < problem->num_S; ++i) {
		for (int j = 0; j < problem->num_W; ++j) {
			if (problem->supplycost[i][j] != d) continue;

			// New variables with objective coefficient set directly
			GRBVar x_new = model.addVar(0.0, problem->goods[i],
			                            static_cast<double>(problem->supplycost[i][j]),
			                            GRB_INTEGER);
			GRBVar z_new = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);

			X2[{i, j}] = x_new;
			Z2[{i, j}] = z_new;

			// Capacity constraint: add x_new or create if warehouse j had no vars before
			if (cap_constrs.count(j)) {
				model.chgCoeff(cap_constrs[j], x_new, 1.0);
			} else {
				GRBLinExpr lhs = x_new - problem->capacity[j] * Y[j];
				cap_constrs[j] = model.addConstr(lhs <= 0);
			}

			// Goods constraint: add x_new or create if supplier i had no vars before
			if (goods_constrs.count(i)) {
				model.chgCoeff(goods_constrs[i], x_new, 1.0);
			} else {
				goods_constrs[i] = model.addConstr(x_new == problem->goods[i]);
			}

			// Upper bound
			model.addConstr(x_new - problem->goods[i] * z_new <= 0);

			// Incompatibility constraints with already-existing Z variables at warehouse j
			for (int i2 : incomp_neighbors[i]) {
				if (Z2.count({i2, j}))
					model.addConstr(z_new + Z2[{i2, j}] <= Y[j]);
			}
		}
	}
}

// Instance name: file name only, without path and without ".txt"
static std::string InstShortName(const std::string& name_inst) {
	std::string inst_short = std::filesystem::path(name_inst).filename().string();
	const std::string ext = ".txt";
	if (inst_short.size() >= ext.size() &&
		inst_short.compare(inst_short.size() - ext.size(), ext.size(), ext) == 0) {
		inst_short.erase(inst_short.size() - ext.size());
	}
	return inst_short;
}

// LOG_CONSTRUCTIVE: Gurobi results of the last optimize() of the constructive phase
// (the model actually solved, with d_final). Unknown values stay "NA".
struct ConstructiveLog {
	std::string status = "NA", solcount = "NA", objval = "NA", objbound = "NA",
				mipgap = "NA", runtime = "NA", grb_timelimit = "NA";
	bool has_obj = false, has_bound = false;
	double z = 0.0, lb = 0.0;
	void reset() { *this = ConstructiveLog(); }
};

static std::string FmtNum(double v) {
	std::ostringstream o;
	o.imbue(std::locale::classic());
	o << std::setprecision(10) << v;
	return o.str();
}

// LOG_CONSTRUCTIVE: append one line to resumen_constructivo.csv (header if the file does not exist)
static void WriteConstructiveCsv(const problem_t* problem, const std::string& name_inst, int dmaxmin,
	int d_final, size_t num_arcs, int num_vars, int num_binvars, int num_constrs,
	const ConstructiveLog& cl, long long time_limit_const, double construction_secs) {

	const std::string csv_file = "resumen_constructivo.csv";
	const long long total_arcs = static_cast<long long>(problem->num_S) * problem->num_W;
	auto ival = [](long long v) { return v < 0 ? std::string("NA") : std::to_string(v); };

	// lb_valid = false when ObjBound is unknown, 0 or a huge negative number (no useful bound)
	const bool lb_valid = cl.has_bound && cl.lb != 0.0 && cl.lb > -1e20;
	std::string gap_pct = "NA";
	if (cl.has_obj && cl.has_bound && cl.z != 0.0) {
		gap_pct = FmtNum(100.0 * (cl.z - cl.lb) / cl.z); // same definition as Gurobi's MIPGap
	}

	std::error_code ec;
	bool write_header = !std::filesystem::exists(csv_file, ec) ||
						std::filesystem::file_size(csv_file, ec) == 0;
	if (ec) write_header = true;

	std::ofstream f(csv_file, std::ios::app);
	f.imbue(std::locale::classic());
	if (write_header) {
		f << "instancia,W,S,media_supplycost,dmaxmin,d_final,total_arcos,num_arcos,pct_arcos,"
			 "num_vars,num_bin_vars,num_constrs,status,solcount,z,lb,lb_valid,grb_mipgap,gap_pct,"
			 "runtime_s,tiempo_construccion_s,time_limit_const,grb_timelimit\n";
	}
	f << InstShortName(name_inst) << ","
	  << problem->num_W << ","
	  << problem->num_S << ","
	  << FmtNum(CalcMeanCost(problem)) << ","
	  << dmaxmin << ","
	  << ival(d_final) << ","
	  << total_arcs << ","
	  << (d_final < 0 ? std::string("NA") : std::to_string(num_arcs)) << ","
	  << (d_final < 0 || total_arcs == 0 ? std::string("NA") : FmtNum(100.0 * num_arcs / total_arcs)) << ","
	  << ival(num_vars) << ","
	  << ival(num_binvars) << ","
	  << ival(num_constrs) << ","
	  << cl.status << ","
	  << cl.solcount << ","
	  << cl.objval << ","
	  << cl.objbound << ","
	  << (lb_valid ? "true" : "false") << ","
	  << cl.mipgap << ","
	  << gap_pct << ","
	  << cl.runtime << ","
	  << FmtNum(construction_secs) << ","
	  << time_limit_const << ","
	  << cl.grb_timelimit << "\n";
}

// LOG_CONSTRUCTIVE: append one line to resumen_fo.csv (header if the file does not exist).
// z_const < 0 means the constructive found no solution -> z/improvement/t_best_const are NA.
void WriteFoCsv(const std::string& name_inst, double z_const, double z_final,
	double t_best_const, double t_best_fo, double total_time,
	long long n_iter_fo, long long n_mejoras_fo) {

	const std::string csv_file = "resumen_fo.csv";
	const bool has_const = z_const >= 0.0;
	auto na_or = [&](double v) { return has_const ? FmtNum(v) : std::string("NA"); };

	std::error_code ec;
	bool write_header = !std::filesystem::exists(csv_file, ec) ||
						std::filesystem::file_size(csv_file, ec) == 0;
	if (ec) write_header = true;

	std::ofstream f(csv_file, std::ios::app);
	f.imbue(std::locale::classic());
	if (write_header) {
		f << "instancia,z_const,z_final,mejora_pct,t_best_const,t_best_fo,tiempo_total_s,n_iter_fo,n_mejoras_fo\n";
	}
	f << InstShortName(name_inst) << ","
	  << na_or(z_const) << ","
	  << na_or(z_final) << ","
	  << (has_const && z_const != 0.0 ? FmtNum(100.0 * (z_const - z_final) / z_const) : std::string("NA")) << ","
	  << na_or(t_best_const) << ","
	  << na_or(t_best_fo) << ","
	  << FmtNum(total_time) << ","
	  << n_iter_fo << ","
	  << n_mejoras_fo << "\n";
}

void Exact::RunGurobiEx_Interval(const problem_t *problem, solution_t *sol, solutionex_t *sol_ex, solutionex_t *best_ex, long long time_limit, const std::chrono::steady_clock::time_point& start_time,
	int DmaxIni, string name_inst, double upp_bound) {

	// Reduced-model size log (resumen_modelo.txt): values of the final (feasible) iteration.
	// They are (re)set at each d right after model.update(), so on exit they describe the
	// last model built/optimized; resumen_modelo.txt is still only written when resumen_found.
	auto construction_start = std::chrono::steady_clock::now();
	bool resumen_found = false;
	int resumen_d = -1;
	size_t resumen_arcs = 0;
	int resumen_numvars = -1;
	int resumen_numbinvars = -1; // LOG_CONSTRUCTIVE
	int resumen_numconstrs = -1;

	// LOG_CONSTRUCTIVE: the CSV line is written when this guard is destroyed, i.e. on every exit
	// path of this function (normal end, early return on time limit, caught GRBException/bad_alloc).
	ConstructiveLog clog;
	struct OnExit {
		std::function<void()> f;
		~OnExit() { try { f(); } catch (...) {} }
	} clog_guard{[&]() {
		if (!log_constructive) return;
		double construction_secs = std::chrono::duration_cast<std::chrono::milliseconds>(
			std::chrono::steady_clock::now() - construction_start).count() / 1000.0;
		WriteConstructiveCsv(problem, name_inst, DmaxIni, resumen_d, resumen_arcs, resumen_numvars,
			resumen_numbinvars, resumen_numconstrs, clog, time_limit, construction_secs);
	}};

	try {

		auto current_time = std::chrono::steady_clock::now();
		auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count();
		long long elapsed_gurobi;

		std::string baseFileName;
		size_t lastDotPos = name_inst.rfind('.');
		if (lastDotPos != std::string::npos && lastDotPos < name_inst.length() - 1) {
			baseFileName = name_inst.substr(0, lastDotPos);
		} else {
			baseFileName = name_inst;
		}

		// Time constructive - Search Profile Experiment
		//std::ofstream conv_file("convergence_" + baseFileName + ".txt");

		GRBEnv env(true);
		env.set(GRB_IntParam_OutputFlag, 0);
		env.set(GRB_IntParam_LogToConsole, 0);
		env.set(GRB_IntParam_Threads, 1);
		env.set(GRB_IntParam_Seed, 12345);
		env.start();

		// SEARCH_PROFILE: open (truncate) once at the start of the constructive phase
		std::ofstream sp_file; // SEARCH_PROFILE
		if (search_profile) { // SEARCH_PROFILE
			sp_file.open("convergence_" + std::filesystem::path(name_inst).stem().string() + ".txt", std::ios::trunc);
		}
		SearchProfileCallback sp_cb(sp_file, start_time); // SEARCH_PROFILE: same start_time as the algorithm

		GRBModel model(env);
		if (search_profile) model.setCallback(&sp_cb); // SEARCH_PROFILE

		int num_Y = problem->num_W;
		double* lb_Y = new double[num_Y];
		double* ub_Y = new double[num_Y];
		double* obj_Y = new double[num_Y];
		char* vtype_Y = new char[num_Y];
		for (int j = 0; j < num_Y; ++j) {
			lb_Y[j] = 0.0;
			ub_Y[j] = 1.0;
			obj_Y[j] = 0.0;
			vtype_Y[j] = GRB_BINARY;
		}
		GRBVar* Y = model.addVars(lb_Y, ub_Y, obj_Y, vtype_Y, nullptr, num_Y);
		delete[] lb_Y;
		delete[] ub_Y;
		delete[] obj_Y;
		delete[] vtype_Y;

		for (int j = 0; j < num_Y; ++j) {
			Y[j].set(GRB_DoubleAttr_Start, 0.0);
		}

		std::unordered_map<std::pair<int,int>, GRBVar, PairHash> X2;
		std::unordered_map<std::pair<int,int>, GRBVar, PairHash> Z2;
		std::map<int, GRBConstr> cap_constrs;
		std::map<int, GRBConstr> goods_constrs;
		buildGurobiModel_ex(model, X2, Z2, Y, problem, DmaxIni, cap_constrs, goods_constrs);

		for (int d = DmaxIni; d <= 200; d++) { // 200 is large enough to ensure that a feasible problem will be found in our instances
			auto start_gurobi = std::chrono::steady_clock::now();

			current_time = std::chrono::steady_clock::now();
			elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count();

			if (d > DmaxIni) {
				extendGurobiModel_ex(model, X2, Z2, Y, problem, d, cap_constrs, goods_constrs);
			}

			model.update();

			// Size of the model about to be optimized (only the last d survives on exit)
			resumen_d = d;
			resumen_arcs = X2.size();
			clog.reset(); // LOG_CONSTRUCTIVE: discard results of previous (infeasible) d
			resumen_numvars = model.get(GRB_IntAttr_NumVars);
			resumen_numbinvars = model.get(GRB_IntAttr_NumBinVars); // LOG_CONSTRUCTIVE
			resumen_numconstrs = model.get(GRB_IntAttr_NumConstrs);

			clog.grb_timelimit = FmtNum(time_limit - elapsed_time / 1000.0); // LOG_CONSTRUCTIVE
			model.set(GRB_DoubleParam_TimeLimit, time_limit - elapsed_time / 1000.0);
			//cout << "T restant " << time_limit - elapsed_time / 1000.0<<endl;
			model.set(GRB_IntParam_MIPFocus, 1);
			model.set(GRB_IntParam_Presolve, 2);
			model.set(GRB_DoubleParam_MIPGap, 0.0);
			if (upp_bound != -1.0) {
				model.set(GRB_DoubleParam_Cutoff, upp_bound);
			}

			if (search_profile) sp_cb.best_obj = 1e30; // SEARCH_PROFILE: reset best value for each resolve (each d)

			model.optimize();

			// LOG_CONSTRUCTIVE: each attribute in its own try/catch, a missing value stays "NA"
			if (log_constructive) {
				try { clog.status = std::to_string(model.get(GRB_IntAttr_Status)); } catch (GRBException&) {}
				int cl_solcount = -1;
				try { cl_solcount = model.get(GRB_IntAttr_SolCount); clog.solcount = std::to_string(cl_solcount); } catch (GRBException&) {}
				try { clog.lb = model.get(GRB_DoubleAttr_ObjBound); clog.objbound = FmtNum(clog.lb); clog.has_bound = true; } catch (GRBException&) {}
				try { clog.runtime = FmtNum(model.get(GRB_DoubleAttr_Runtime)); } catch (GRBException&) {}
				if (cl_solcount > 0) {
					try { clog.z = model.get(GRB_DoubleAttr_ObjVal); clog.objval = FmtNum(clog.z); clog.has_obj = true; } catch (GRBException&) {}
					try { clog.mipgap = FmtNum(model.get(GRB_DoubleAttr_MIPGap)); } catch (GRBException&) {}
				}
			}

			// SEARCH_PROFILE: bound/gap/status per d -> convergence_bound_<inst>.txt and console.
			// Each attribute is read in its own try/catch so a missing value never reaches the outer catch.
			if (search_profile) { // SEARCH_PROFILE
				try {
					int sp_status = model.get(GRB_IntAttr_Status);
					int sp_solcount = model.get(GRB_IntAttr_SolCount);
					auto sp_fmt = [](double v) { std::ostringstream o; o << std::setprecision(10) << v; return o.str(); };
					std::string sp_obj = "NA", sp_bound = "NA", sp_gap = "NA";
					try { sp_bound = sp_fmt(model.get(GRB_DoubleAttr_ObjBound)); } catch (GRBException&) {}
					if (sp_solcount > 0) {
						try { sp_obj = sp_fmt(model.get(GRB_DoubleAttr_ObjVal)); } catch (GRBException&) {}
						try { sp_gap = sp_fmt(model.get(GRB_DoubleAttr_MIPGap)); } catch (GRBException&) {}
					}
					std::string sp_line = std::to_string(d) + " " + std::to_string(sp_status) + " " + sp_obj + " " + sp_bound + " " + sp_gap;
					std::ofstream sp_bfile("convergence_bound_" + std::filesystem::path(name_inst).stem().string() + ".txt", std::ios::app);
					sp_bfile << sp_line << "\n";
					sp_bfile.flush();
					std::cout << "[SP-BOUND] " << sp_line << std::endl;
				} catch (GRBException&) {}
			}

			current_time = std::chrono::steady_clock::now();
			elapsed_gurobi = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_gurobi).count();

			int solcount = model.get(GRB_IntAttr_SolCount);
			sol_ex->solcount = solcount;

			if (solcount > 0) {
				sol_ex->cost = model.get(GRB_DoubleAttr_ObjVal);
				sol_ex->Dmax_sol = 0;
				for (int i = 0; i < problem->num_S; i++) {
					for (int j = 0; j < problem->num_W; j++) {
						if (X2.count({i, j})) {
							sol_ex->X[i][j] = abs(X2[{i,j}].get(GRB_DoubleAttr_X));
						} else {
							sol_ex->X[i][j] = 0;
						}
					}
				}
				for (int j = 0; j < problem->num_W; j++) {
					sol_ex->fopen[j] = 0;
					for (int i = 0; i < problem->num_S; i++) {
						if (sol_ex->X[i][j] > 0) {
							sol_ex->fopen[j] = 1;
							break;
						}
					}
				}
				int status = model.get(GRB_IntAttr_Status);
				if (status == GRB_OPTIMAL || status == GRB_SUBOPTIMAL) {
					sol_ex->gap = model.get(GRB_DoubleAttr_MIPGap);
				} else {
					sol_ex->gap = -1.0;
				}
				sol_ex->run_time = elapsed_gurobi / 1000.0;
				copysolution_ex(problem, best_ex, sol_ex);
				if (search_profile) { // SEARCH_PROFILE
					std::cout << "[SP] CONST obj=" << best_ex->cost << std::endl;
				}
				best_ex->target_time = elapsed_time / 1000.0;

				resumen_found = true; // resumen_d/arcs/numvars/numconstrs already set after model.update()
				break;
			}

			if (solcount == 0) {
				sol_ex->cost = -1;
				sol_ex->Dmax_sol = -1;
				sol_ex->target_time = 0;
				sol_ex->gap = -1.0;
				elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count();
				if (elapsed_time / 1000.0 >= time_limit) {
					return;
				}
			}

			current_time = std::chrono::steady_clock::now();
			elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count();
		}

		ConversionSolexToSol(problem, sol, best_ex);

		if (resumen_found) {
			double construction_secs = std::chrono::duration_cast<std::chrono::milliseconds>(
				std::chrono::steady_clock::now() - construction_start).count() / 1000.0;

			const std::string resumen_file = "resumen_modelo.txt";

			const std::string inst_short = InstShortName(name_inst);

			double total_arcs = static_cast<double>(problem->num_S) * problem->num_W;
			double arcs_frac = total_arcs > 0 ? resumen_arcs / total_arcs : 0.0;

			std::error_code ec;
			bool write_header = !std::filesystem::exists(resumen_file, ec) ||
								std::filesystem::file_size(resumen_file, ec) == 0;
			if (ec) write_header = true;

			std::ofstream resumen(resumen_file, std::ios::app);
			resumen.imbue(std::locale::classic());
			if (write_header) {
				resumen << "instancia,d_final,num_arcos,frac_arcos,num_vars,num_constrs,tiempo_construccion_s,fo_inicial\n";
			}
			resumen << inst_short << ","
					<< resumen_d << ","
					<< resumen_arcs << ","
					<< std::fixed << std::setprecision(4) << arcs_frac << ","
					<< resumen_numvars << ","
					<< resumen_numconstrs << ","
					<< std::setprecision(3) << construction_secs << ","
					<< std::defaultfloat << std::setprecision(10) << best_ex->cost << "\n";
		} // ofstream closes here

		} catch (GRBException e) {
			std::cerr << "Error Gurobi: " << e.getMessage() << std::endl;
			std::cerr << "Error Code: " << e.getErrorCode() << std::endl;
		} catch (const std::exception& e) {
			std::cerr << "Error: " << e.what() << std::endl;
		}
}

vector<int> Exact::KnapsackClosedFac(const problem_t* problem, solution_t* solution,
                               const vector<int>& v_cust_subproblem,
                               int n_supp_close,
                               const vector<int>& v_supp_open) {
    try {
        vector<pair<int,double>> closest = ClosestFac_norm2(problem, solution, v_cust_subproblem);

        int total_demand = 0;
        for (int c : v_cust_subproblem) {
            for (int j = 1; j < solution->assignment_f[c][0] + 1; j++) {
                int f = solution->assignment_f[c][j];
                auto it = find(v_supp_open.begin(), v_supp_open.end(), f);
                if (it != v_supp_open.end()) {
                    total_demand += solution->assignment_c[c][j];
                }
            }
        }

        int n_closed = closest.size();

        GRBEnv env_knap = GRBEnv(true);
        env_knap.set(GRB_IntParam_OutputFlag, 0);
        env_knap.start();
        GRBModel model = GRBModel(env_knap);

        vector<double> lb_y(n_closed, 0.0);
        vector<double> ub_y(n_closed, 1.0);
        vector<double> obj_y(n_closed);
        vector<char> vtype_y(n_closed, GRB_BINARY);

        for (int idx = 0; idx < n_closed; idx++) {
            obj_y[idx] = closest[idx].second;
        }

        GRBVar* y = model.addVars(lb_y.data(), ub_y.data(), obj_y.data(),
                                   vtype_y.data(), nullptr, n_closed);

        model.update();

        GRBLinExpr cap_expr = 0;
        for (int idx = 0; idx < n_closed; idx++) {
            cap_expr += problem->capacity[closest[idx].first] * y[idx];
        }
        model.addConstr(cap_expr >= total_demand, "capacity");

        model.set(GRB_DoubleParam_TimeLimit, 5.0);
        model.optimize();

        vector<int> selected;
        if (model.get(GRB_IntAttr_SolCount) > 0) {
            for (int idx = 0; idx < n_closed; idx++) {
                if (y[idx].get(GRB_DoubleAttr_X) > 0.5) {
                    selected.push_back(closest[idx].first);
                }
            }
        } else {
            vector<int> fallback = ClosestFac_norm(problem, solution, v_cust_subproblem);
            for (int i = 0; i < min(n_supp_close, (int)fallback.size()); i++) {
                selected.push_back(fallback[i]);
            }
        }

        delete[] y;
        return selected;

    } catch (GRBException& e) {
        vector<int> fallback = ClosestFac_norm(problem, solution, v_cust_subproblem);
        vector<int> selected;
        for (int i = 0; i < min(n_supp_close, (int)fallback.size()); i++) {
            selected.push_back(fallback[i]);
        }
        return selected;
    }
}
