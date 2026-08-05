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

    	if (exact_mode == 2) {
    		model.set(GRB_StringParam_LogFile, "gurobi_output_" + baseFileName + ".log"); // Name of the log file
    	}


    	model.set(GRB_IntParam_LogToConsole, 1);
    	current_time = std::chrono::steady_clock::now();
    	elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count();

    	if (tlimit != -1 && tlimit-elapsed_time/1000.0 > 0) {
    		model.set(GRB_DoubleParam_TimeLimit, tlimit-elapsed_time/1000.0);

    		/*if (exact_mode == 2) {
    			model.set(GRB_IntParam_MIPFocus, 1);
    			model.set(GRB_IntParam_Presolve, 2);
    		}*/

    		current_time = std::chrono::steady_clock::now();
    		elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count();
    		model.optimize();

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

void Exact::RunGurobiEx_Interval(const problem_t *problem, solution_t *sol, solutionex_t *sol_ex, solutionex_t *best_ex, long long time_limit, const std::chrono::steady_clock::time_point& start_time,
	int DmaxIni, string name_inst, double upp_bound) {

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

		GRBModel model(env);

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

		for (int d = DmaxIni; d <= 60; d++) { // 60 is large enough to ensure that a feasible problem will be found
			auto start_gurobi = std::chrono::steady_clock::now();

			current_time = std::chrono::steady_clock::now();
			elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count();

			if (d > DmaxIni) {
				extendGurobiModel_ex(model, X2, Z2, Y, problem, d, cap_constrs, goods_constrs);
			}

			model.update();
			model.set(GRB_DoubleParam_TimeLimit, time_limit - elapsed_time / 1000.0);
			//cout << "T restant " << time_limit - elapsed_time / 1000.0<<endl;
			model.set(GRB_IntParam_MIPFocus, 1);
			model.set(GRB_IntParam_Presolve, 2);
			model.set(GRB_DoubleParam_MIPGap, 0.0);
			if (upp_bound != -1.0) {
				model.set(GRB_DoubleParam_Cutoff, upp_bound);
			}

			model.optimize();

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
				best_ex->target_time = elapsed_time / 1000.0;
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
		} catch (GRBException e) {
			std::cerr << "Error Gurobi: " << e.getMessage() << std::endl;
			std::cerr << "Error Code: " << e.getErrorCode() << std::endl;
		} catch (const std::exception& e) {
			std::cerr << "Error: " << e.what() << std::endl;
		}
}