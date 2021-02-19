#ifndef INCLUDED_CORE_MATRIX
#define INCLUDED_CORE_MATRIX

#include <vector>
#include "geometry.h"

void make_coeff_matrix_proto(MatrixXd& K_mat, MatrixXd& G_mat, MatrixXd& R_mat, vector<Element> Elements, vector<Source> Sources);
void transform_matrix(MatrixXd& A, VectorXd& b, VectorXd& b_value_vec, vector<int>& Dirichlet_indexes, vector<int>& Neumann_indexes,
  MatrixXd K_mat, MatrixXd G_mat, vector<Element> Elements);

#endif
