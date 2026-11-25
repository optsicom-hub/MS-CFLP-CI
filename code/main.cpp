#include "problem.h"
#include "solution.h"
#include "solution_exact.h"
#include "globals.h"
#include "exact.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <vector>
#include <fstream>
#include <random>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <thread>

using std::filesystem::path;
using namespace std;
namespace fs = std::filesystem;

problem_t *readinstanceOpt(problem_t *problem, string fichero);
int minimum(int a, int b);
int maximum(int a, int b);
void saveToFileString(const std::string &filename, std::vector<string> data);

int main(int argc, char **argv) {

	//search_profile = false; // SEARCH_PROFILE: set to false to disable

	auto *problem = new problem_t;
	string name_inst = argv[1];

	for (int i = 1; i < argc; ++i) {
		std::string arg_str = argv[i];

		if (arg_str == "--pct_open") {
			if (i + 1 < argc) {
				try {
					pct_open = std::stod(argv[++i]);
				} catch (const std::exception& e) {
					std::cerr << "ERROR: invalid value for --pct_open: " << argv[i] << " (" << e.what() << ")" << std::endl;
					return 1;
				}
			} else {
				std::cerr << "ERROR: missing value for --pct_open." << std::endl;
				return 1;
			}
		}else if (arg_str == "--pct_close") {
			if (i + 1 < argc) {
				try {
					pct_close = std::stod(argv[++i]);
				} catch (const std::exception& e) {
					std::cerr << "ERROR: invalid value for --pct_close: " << argv[i] << " (" << e.what() << ")" << std::endl;
					return 1;
				}
			} else {
				std::cerr << "ERROR: missing value for --pct_close." << std::endl;
				return 1;
			}
		}else if (arg_str == "--subproblem_size") {
			if (i + 1 < argc) {
				try {
					subproblem_size = std::atoi(argv[++i]);
				} catch (const std::exception& e) {
					std::cerr << "ERROR: invalid value for --subproblem_size: " << argv[i] << " (" << e.what() << ")" << std::endl;
					return 1;
				}
			} else {
				std::cerr << "ERROR: missing value for --subproblem_size." << std::endl;
				return 1;
			}
		}else if (arg_str == "--strategy") {
			if (i + 1 < argc) {
				try {
					strategy = std::atoi(argv[++i]);
				} catch (const std::exception& e) {
					std::cerr << "ERROR: invalid value for --strategy: " << argv[i] << " (" << e.what() << ")" << std::endl;
					return 1;
				}
			} else {
				std::cerr << "ERROR: missing value for --strategy." << std::endl;
				return 1;
			}
		}else if (arg_str == "--seed") {
			if (i + 1 < argc) {
				try {
					seed = std::atoi(argv[++i]);
				} catch (const std::exception& e) {
					std::cerr << "ERROR: invalid value for --seed: " << argv[i] << " (" << e.what() << ")" << std::endl;
					return 1;
				}
			} else {
				std::cerr << "ERROR: missing value for --seed." << std::endl;
				return 1;
			}
		}else if (arg_str == "--divisor_time") {
			if (i + 1 < argc) {
				try {
					divisor_time = std::atoi(argv[++i]);
				} catch (const std::exception& e) {
					std::cerr << "ERROR: invalid value for --divisor_time: " << argv[i] << " (" << e.what() << ")" << std::endl;
					return 1;
				}
			} else {
				std::cerr << "ERROR: missing value for --divisor_time." << std::endl;
				return 1;
			}
		}else if (arg_str == "--log_constructive") { // LOG_CONSTRUCTIVE
			if (i + 1 < argc) {
				try {
					log_constructive = std::stoi(argv[++i]) != 0;
				} catch (const std::exception& e) {
					std::cerr << "ERROR: invalid value for --log_constructive: " << argv[i] << " (" << e.what() << ")" << std::endl;
					return 1;
				}
			} else {
				std::cerr << "ERROR: missing value for --log_constructive." << std::endl;
				return 1;
			}
		}else if (arg_str == "--search_profile") { // SEARCH_PROFILE
			if (i + 1 < argc) {
				try {
					search_profile = std::stoi(argv[++i]) != 0;
				} catch (const std::exception& e) {
					std::cerr << "ERROR: invalid value for --search_profile: " << argv[i] << " (" << e.what() << ")" << std::endl;
					return 1;
				}
			} else {
				std::cerr << "ERROR: missing value for --search_profile." << std::endl;
				return 1;
			}
		}
	}


	problem = readinstanceOpt(problem, name_inst);

	// Needed to obtain Pass Mark results
	const double ratio_PASS = 2650.0/2998.0;

	// Select depending on the execution
	//const double tmax = double(problem->num_W)*ratio_PASS;
	//const double tmax = 10.0 * sqrt(double(problem->num_W)) * ratio_PASS;
	//const double tmax = 10*sqrt(problem->num_W);
	// double tmax = 7200.0;
	const double tmax = problem->num_W;
	//const double tmax = 10*sqrt(problem->num_W);

	const long long time_limit1 = tmax/divisor_time;
	Exact ex;

	int time_limit_ls;
	int time_limit_ls_ex0;
	long long elapsed_time_ex0;

	double total_time_secs = 0.0;

	for (int execution = 0; execution < 1; execution++) {
		auto *solution = create_solution(problem->num_S, problem->num_W);

		auto *best_solution = create_solution(problem->num_S, problem->num_W);
		init_solution(problem,best_solution);

		auto *solution_exact = create_solutionex(problem->num_S, problem->num_W);

		auto *best_exact = create_solutionex(problem->num_S, problem->num_W);

		auto start_time = std::chrono::steady_clock::now();
		long long elapsed_time = 0;
		auto current_time = std::chrono::steady_clock::now();

		// Instant (in seconds since start_time) at which the best objective value found so far was obtained
		double best_found_time = 0.0;
		double t_best_const = 0.0; // LOG_CONSTRUCTIVE: best_found_time right after the constructive phase

			int Dmaxheur = 0;

			int Dmaxmin;

			if (execution == 0) {

				Dmaxmin = CalcCmaxMin(problem);

				/*{
					double meanCost = CalcMeanCost(problem);
					std::ios_base::fmtflags oldFlags = cout.flags();
					std::streamsize oldPrecision = cout.precision();
					cout << "Mean supply cost: " << std::fixed << std::setprecision(4) << meanCost << endl;
					cout << "Dmaxmin: " << Dmaxmin << endl;
					cout.flags(oldFlags);
					cout.precision(oldPrecision);

				}*/

				current_time = std::chrono::steady_clock::now();
				elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count();

				best_exact->cost = -1;
				best_exact->run_time = 0;

				// target_time is expressed in seconds elsewhere (see RunGurobiEx_Interval), keep this placeholder consistent
				best_exact->target_time = elapsed_time / 1000.0;
				best_exact->Dmax_sol = Dmaxheur;

				//comment to execute MILP in 7200 seconds
				ex.RunGurobiEx_Interval(problem, solution, solution_exact, best_exact, time_limit1-elapsed_time/1000.0, start_time, Dmaxmin, name_inst, -1.0);

				// RunGurobiEx_Interval sets best_exact->target_time to the instant (in seconds since
				// start_time) at which the MILP solution that becomes the initial best_solution was found.
				if (best_exact->cost != -1) {
					best_found_time = static_cast<double>(best_exact->target_time);
				}
				t_best_const = best_found_time; // LOG_CONSTRUCTIVE

				current_time = std::chrono::steady_clock::now();
				double elapsed_interval = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count() / 1000.0;

				// uncomment to execute MILP in 7200 seconds
				/*ex.runGurobi_ex(problem, solution_exact, name_inst,1," ", -1.0, 0.0, 7200, -1.0, 0, start_time);
				{
					double elapsed_milp = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time).count() / 1000.0;
					std::ofstream f_milp("results_milp_W.txt", std::ios::app);
					f_milp << std::fixed << std::setprecision(4); // EXACT_FULL
					f_milp << name_inst << " " << solution_exact->cost << " " << elapsed_milp
					       << " " << solution_exact->lower_bound << " " << solution_exact->mip_gap << "\n"; // EXACT_FULL
				}
				 exit(0);*/
			}

		if (execution == 0) {
			current_time = std::chrono::steady_clock::now();
			elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count();
			elapsed_time_ex0 = elapsed_time;
		}else {
			elapsed_time = elapsed_time_ex0;
		}

		copy_solution(best_solution, solution);

		delete_solutionex(solution_exact,problem->num_S);

		vector<vector<int> > count_closestfac;

		int fac;

		for (int i = 0; i < problem->num_W; i++) {
			count_closestfac.push_back({i, 0});
		}

		for (int i = 0; i < problem->num_S; i++) {
			fac = problem->c_best[i][0];
			count_closestfac[fac][1]++;
		}
		for (int i = 0; i < problem->num_W; i++) {
			count_closestfac[i][1]++;
		}

		std::vector<int> v_facilities;
		for (const auto &row: count_closestfac) {
			v_facilities.push_back(row[1]);
		}

		std::discrete_distribution<int> dist_best_option(v_facilities.begin(), v_facilities.end());

		if (execution == 0) {
			current_time = std::chrono::steady_clock::now();
			elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count();
		}

		long long elapsed_gurobi;

		if (execution == 0) {
			time_limit_ls = maximum(tmax-time_limit1, tmax-elapsed_time/1000.0);
			time_limit_ls_ex0 = time_limit_ls;

		}else {
			time_limit_ls = time_limit_ls_ex0;
		}

			pct_close = 1.0-pct_open;

			copy_solution(best_solution,solution);

			int n_supp_open = 0;
			int n_supp_close = 0;
			auto start_ls = std::chrono::steady_clock::now();

			long long elapsed_ls = 0;

			std::ofstream sp_file; // SEARCH_PROFILE
			if (search_profile) { // SEARCH_PROFILE
				sp_file.open("convergence_" + fs::path(name_inst).stem().string() + ".txt", std::ios::app);
			}

			long long n_iter_fo = 0, n_mejoras_fo = 0; // LOG_CONSTRUCTIVE: F&O iterations / accepted improvements

			while(elapsed_ls/1000.0 < time_limit_ls) {
				n_iter_fo++; // LOG_CONSTRUCTIVE

				n_supp_open = minimum(ceil(pct_open*subproblem_size),best_solution->n_open);
				n_supp_close = minimum(ceil(pct_close*subproblem_size),problem->num_W-best_solution->n_open);

				if (n_supp_open + n_supp_close > subproblem_size) {
					n_supp_close--;
				}
				vector<int> v_supp;
				vector<int> v_supp_old;

				std::vector<int> v_num_customers;

				for (int j = 0; j < problem->num_W; j++) {
					v_num_customers.push_back(best_solution->facility_client[j][0]);

				}

				std::discrete_distribution<int> dist_num_cust(v_num_customers.begin(), v_num_customers.end());
				while(v_supp.size()<n_supp_open) {
					int f = dist_num_cust(rng);

					if (best_solution->fopen[f] == 1) {
						// check that it is not in the vector
						auto it_found = std::find(v_supp.begin(), v_supp.end(), f);
						if (it_found == v_supp.end()) {
							v_supp.push_back(f);
						}
					}
				}

				vector<int> v_cust_subproblem;

				for (int f : v_supp) {
					for (int z = 1; z < best_solution->facility_client[f][0]+1; z++) {
						int c = best_solution->facility_client[f][z];

						auto it_found = std::find(v_cust_subproblem.begin(), v_cust_subproblem.end(), c);

						if (it_found == v_cust_subproblem.end()) {
							v_cust_subproblem.push_back(c);
						}
					}
				}


				sort(v_cust_subproblem.begin(), v_cust_subproblem.end());
				int ncli = static_cast<int>(v_cust_subproblem.size());


				vector<int> closest_suppliers;
				if (strategy == 1) {
					closest_suppliers = ClosestFac_norm(problem, best_solution, v_cust_subproblem);
				}

				if (strategy < 3) {
					int iter = 0;

					while(v_supp.size()<n_supp_close+n_supp_open) {
						int f;

						if (strategy == 1) {
							f = closest_suppliers[iter];
						}else if(strategy == 2) {
							f = dist_best_option(rng);
						}

						if (best_solution->fopen[f] == 0) {
							auto it_found = std::find(v_supp.begin(), v_supp.end(), f);

							if (it_found == v_supp.end()) {
								v_supp.push_back(f);
							}

						}
						if (strategy == 1) {
							iter++;
						}

					}
				}else {
					vector<int> knap_selected = ex.KnapsackClosedFac(
					problem, best_solution,
					v_cust_subproblem,
					n_supp_close,
					v_supp
					);

					for (int f : knap_selected) {
						v_supp.push_back(f);
					}
				}


				sort(v_supp.begin(), v_supp.end());

				//int nsupp = n_supp_open + n_supp_close;
				int nsupp = (strategy == 3) ? static_cast<int>(v_supp.size()) : n_supp_open + n_supp_close;

				auto *subproblem = new problem_t;
				subproblem = create_problem(ncli, nsupp);

				auto *solution_exact = create_solutionex(ncli, nsupp);

				GenerateSubproblem(problem, subproblem, best_solution, ncli, nsupp, v_cust_subproblem, v_supp);

				int coste_pre = GetsolutionCost(subproblem, best_solution, v_cust_subproblem, v_supp, solution_exact);

				current_time = std::chrono::steady_clock::now();
				elapsed_ls = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_ls).count();

				auto start_gurobi = std::chrono::steady_clock::now();

				if (time_limit_ls-elapsed_ls/1000.0>0) {
					ex.runGurobi_ex(subproblem, solution_exact, name_inst,1, " ", -1.0, 0.0, time_limit_ls-elapsed_ls/1000.0, coste_pre-1, 0, start_time);
				}
				current_time = std::chrono::steady_clock::now();
				elapsed_gurobi = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_gurobi).count();

				if (solution_exact->cost < coste_pre && solution_exact->solcount > 0) {
					n_mejoras_fo++; // LOG_CONSTRUCTIVE
					for (int c : v_cust_subproblem) {

						for (int k = best_solution->assignment_f[c][0]; k > 0; k--) {

							auto it = find(v_supp.begin(), v_supp.end(), best_solution->assignment_f[c][k]);

							if (it != v_supp.end()){
								DropCustomer(problem, best_solution, c, best_solution->assignment_f[c][k]);
							}
						}
					}

					for (int i = 0; i < subproblem->num_S; i++) {
						int c = v_cust_subproblem[i];
						for (int j = 0; j < subproblem->num_W; j++) {
							int f = v_supp[j];
							if (solution_exact->X[i][j] > 0) {
								InsertCustomer_2q(problem, best_solution, c, f, solution_exact->X[i][j], best_solution->c_done[c]);
							}
						}
					}

					if (search_profile) { // SEARCH_PROFILE
						std::cout << "[SP] FO coste=" << best_solution->cost << std::endl;
					}

					// Record the instant (seconds since start_time) at which the solution improved
					current_time = std::chrono::steady_clock::now();
					best_found_time = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count() / 1000.0;

					if (search_profile) { // SEARCH_PROFILE
						double sp_t = std::chrono::duration<double>(std::chrono::steady_clock::now() - start_time).count();
						sp_file << std::fixed << std::setprecision(3) << sp_t << " FO " << best_solution->cost << "\n";
						sp_file.flush();
					}

				}

				current_time = std::chrono::steady_clock::now();
				elapsed_ls = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_ls).count();

				delete_solutionex(solution_exact, ncli);

			}

		// Uncomment to save final solution
		//string sol_filename = "solution" + to_string(execution) + name_inst;
		//saveSolution(best_solution, sol_filename, problem->num_S, problem->num_W);

		int fo_final = best_solution->cost;
		total_time_secs = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time).count() / 1000.0;

		// LOG_CONSTRUCTIVE: z_const = fo_inicial of the constructive (-1 if none), t_best_fo = best_found_time
		// (same value written to summary.txt), independent of search_profile
		if (log_constructive) {
			WriteFoCsv(name_inst, best_exact->cost, fo_final, t_best_const, best_found_time, total_time_secs,
				n_iter_fo, n_mejoras_fo);
		}

		delete_solution(solution);
		delete_solution(best_solution);
		delete_solutionex(best_exact, problem->num_S);

		// Add execution to summary file: instance, total runtime, final objective value and the
		// instant (seconds since start_time) at which that final objective value was obtained
		vector<string> data = {name_inst, to_string(total_time_secs), to_string(fo_final), to_string(best_found_time)};//,
		saveToFileString("summary.txt", data);

		cout << fo_final <<endl;
	}



	return 0;

}


problem_t *readinstanceOpt(problem_t *problem_, const string fichero)
{
    int W, S;

    ifstream in(fichero);
    if (!in) {
        fprintf(stdout, "Error Instance Reading");
        return nullptr;
    }


    in.sync_with_stdio(false);

    char dummy[200];
    char d;

    in.getline(dummy, 200, '=');
    in >> W;
    in.getline(dummy, 200, '=');
    in >> S;

    problem_ = create_problem(S, W);
    // Capacity
    in.getline(dummy, 200, '[');
    capacity_norm = 0;
    for (int j = 0; j < W; j++) {
        in >> problem_->capacity[j] >> d;
        if (problem_->capacity[j] > capacity_norm) {
            capacity_norm = problem_->capacity[j];
        }
    }
    // Fixed Costs
    in.getline(dummy, 200, '[');
    fixed_norm = 0;
    for (int j = 0; j < W; j++) {
        in >> problem_->fixedcost[j] >> d;
        if (problem_->fixedcost[j] > fixed_norm) {
            fixed_norm = problem_->fixedcost[j];
        }
    }
    // Goods
    in.getline(dummy, 200, '[');
    for (int j = 0; j < S; j++) {
        in >> problem_->goods[j] >> d;
    }
    // Supply Costs
    in.getline(dummy, 200, '[');
    in.getline(dummy, 200, '|');

    sc_norm = 0;
    for (int i = 0; i < S; i++) {
        for (int j = 0; j < W; j++) {
            in >> problem_->supplycost[i][j] >> d;
            if (problem_->supplycost[i][j] > sc_norm) {
                sc_norm = problem_->supplycost[i][j];
            }
        }
    }
    // Incompatibilities
    in.getline(dummy, 200, '=');
    in >> problem_->num_Inc >> d;
    in.getline(dummy, 200, '|');


    // Initialize incompatibility matrix
    for (int i = 0; i < S; i++) {
        memset(problem_->inc[i], 0, S * sizeof(int));
    }

    int fila, col;
    for (int i = 0; i < problem_->num_Inc; i++) {
        in >> fila >> d >> col >> d;
        problem_->inc[fila-1][col-1] = 1;
        problem_->inc[col-1][fila-1] = 1;
    }
    in.close();

    // Compute f_best_c (best suppliers for each warehouse)
    std::vector<std::pair<int, int>> temp_pairs;
    temp_pairs.reserve(S);

    for (int i = 0; i < W; i++) {
        temp_pairs.clear();

        for (int j = 0; j < S; j++) {
            temp_pairs.push_back({problem_->supplycost[j][i], j});
        }

        std::sort(temp_pairs.begin(), temp_pairs.end());

        for (int f = 0; f < S; f++) {
            problem_->f_best_c[i][f] = temp_pairs[f].second;
        }
    }


    for (int i = 0; i < S; i++) {
        double min_scost;
        int facility = 0;

        // Compute cost for first warehouse
    	min_scost = problem_->supplycost[i][0] + problem_->fixedcost[0] / problem_->capacity[0];


        // Find minimum across all warehouses
        for (int j = 1; j < W; j++) {
            double cost;
        	cost = problem_->supplycost[i][j] + problem_->fixedcost[j] / problem_->capacity[j];

            if (cost < min_scost) {
                min_scost = cost;
                facility = j;
            }
        }

        problem_->c_best_f[i][0] = i;
        problem_->c_best_f[i][1] = facility;
        problem_->c_best_f[i][2] = min_scost;
    }

    std::sort(problem_->c_best_f, problem_->c_best_f + S, [](const double* a, const double* b) {
        return a[2] > b[2];
    });


    std::vector<std::pair<double, int>> v_facilities;
    std::vector<std::pair<double, int>> v_facilities_sc;
    v_facilities.reserve(W);
    v_facilities_sc.reserve(W);

    for (int i = 0; i < S; i++) {
        v_facilities.clear();
        v_facilities_sc.clear();

        for (int j = 0; j < W; j++) {
            double cost;
        	cost = problem_->supplycost[i][j] + problem_->fixedcost[j] / problem_->capacity[j];


            v_facilities.push_back({cost, j});
            v_facilities_sc.push_back({problem_->supplycost[i][j], j});
        }

        std::sort(v_facilities.begin(), v_facilities.end());
        std::sort(v_facilities_sc.begin(), v_facilities_sc.end());

        for (int k = 0; k < W; k++) {
            problem_->c_best[i][k] = v_facilities[k].second;
            problem_->c_best_sc[i][k] = v_facilities_sc[k].second;
        }
    }

    std::vector<std::pair<double, int>> v_worst;
    v_worst.reserve(W);

    for (int i = 0; i < S; i++) {
        v_worst.clear();

        for (int j = 0; j < W; j++) {
            v_worst.push_back({problem_->supplycost[i][j] + problem_->fixedcost[j], j});
        }

        std::sort(v_worst.begin(), v_worst.end(), std::greater<std::pair<double, int>>());

        for (int k = 0; k < W; k++) {
            problem_->c_worst[i][k] = v_worst[k].second;
        }
    }

    for (int i = 0; i < S; i++) {
        problem_->c_best_f_param[i][0] = i;

        // Find position of supplier i in c_best_f
        for (int j = 0; j < S; j++) {
            if (problem_->c_best_f[j][0] == i) {
                problem_->c_best_f_param[i][1] = j;
                break;
            }
        }
    }


    return problem_;
}

void saveToFileString(const std::string &filename, std::vector<string> data) {

	std::ofstream outFile(filename, std::ios::app);
	if (outFile.is_open()) {
		outFile << std::fixed;

		for (const auto &value : data) {

			std::ostringstream oss;
			oss << value;

			std::string strValue = oss.str();
			std::replace(strValue.begin(), strValue.end(), '.', ',');

			outFile << strValue << " ";
		}

		outFile << "\n"; 
		outFile.close();
	} else {
		std::cerr << "Error." << std::endl;
	}
}

int minimum(int a, int b) {
	if (a < b) {
		return a;
	}else {
		return b;
	}
}

int maximum(int a, int b) {
	if (a < b) {
		return b;
	}else {
		return a;
	}
}
