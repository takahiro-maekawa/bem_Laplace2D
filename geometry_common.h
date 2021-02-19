#ifndef INCLUDED_GEOMETRY_COMMON
#define INCLUDED_GEOMETRY_COMMON

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Core>
#include "geometry.h"

using namespace std;

void make_Elements(vector<Element>& Elements, vector<Node> Nodes, vector<int> boundary_types, vector<double> boundary_values);

#endif
