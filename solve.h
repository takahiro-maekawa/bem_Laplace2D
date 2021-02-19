#ifndef INCLUDED_SOLVE
#define INCLUDED_SOLVE

#include <vector>
#include "geometry.h"

using namespace std;

void solve_equations_and_get_result(MatrixXd A, VectorXd b, VectorXd& phai);

VectorXd solve_problem_and_get_beta_vector(vector<Element> Elements, vector<Source> Sources);

vector<double> insert_result(vector<vector<double>> positions, vector<Source> Sources, VectorXd beta_vector);

#endif
