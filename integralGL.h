#ifndef INCLUDED_INTEGRALGL
#define INCLUDED_INTEGRALGL
#include <vector>

using namespace std;

void return_eta_weights(vector<double>& etas, vector<double>& weights, int n_pt);

void insert_1d_position(vector<double> &rs_elem_node, double r1, double r2, vector<double> etas);
#endif
