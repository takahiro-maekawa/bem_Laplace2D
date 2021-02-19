#include <vector>
#include "geometry_common.h"
#include "core_matrix.h"
#include "general_func.h"

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>

using namespace std;
using namespace Eigen;

//called in solve_equations_and_get_result
//得られた方程式の解やディレクレ条件から境界要素における電位ポテンシャルを取得
VectorXd get_phai_vector_from_solution(VectorXd solution_vec, VectorXd b_value_vec, vector<int> Dirichlet_indexes, vector<int> Neumann_indexes){
  int n_elem = solution_vec.size();
  int n_elem_Dir = Dirichlet_indexes.size();
  VectorXd phai; phai.resize(n_elem);

  for (int i=0; i<n_elem; ++i){
    if (i < n_elem_Dir){
      phai(Dirichlet_indexes[i]) = b_value_vec(i);
    }else{
      phai(Neumann_indexes[i-n_elem_Dir]) = solution_vec(i);
    }
  }

  return phai;
}

//called by solve_problem_and_get_beta_vector
//calls     get_phai_vector_from_solution
//連立方程式を解いた後、境界要素における値を計算、作るまでには[core_matrix.cpp]を利用
void solve_equations_and_get_result(MatrixXd A, VectorXd b, VectorXd& phai,
   VectorXd b_value_vec, vector<int> Dirichlet_indexes, vector<int> Neumann_indexes){

  VectorXd solution_vec_ = A.fullPivLu().solve(b);

  VectorXd solution_vec;
  phai = get_phai_vector_from_solution(solution_vec_, b_value_vec, Dirichlet_indexes, Neumann_indexes);

  //phai = {2,4};
  cout << "phai on boundary\n";
  for (int i = 0; i<24; ++i){
    cout<<i<<"\t"<<phai(i) <<"\n";
  }
  cout<<"\n";
}

//called by main()
//calls     make_coeff_matrix_proto[core_matrix], transform_matrix[core_matrix],
//          solve_equations_and_get_result
//連立方程式の構築~境界要素におけるポテンシャルを取得
VectorXd solve_problem_and_get_beta_vector(vector<Element>Elements, vector<Source>Sources){
  VectorXd beta_vector;

  vector<int> Dirichlet_indexes, Neumann_indexes;
  MatrixXd K_mat, G_mat, R_mat, A;
  VectorXd b, b_value_vec, phai;

  // in core_matrix.cpp
  make_coeff_matrix_proto(K_mat, G_mat, R_mat, Elements, Sources);
  transform_matrix(A, b, b_value_vec, Dirichlet_indexes, Neumann_indexes, K_mat, G_mat, Elements);

  solve_equations_and_get_result(A, b, phai, b_value_vec, Dirichlet_indexes, Neumann_indexes);

  beta_vector = R_mat * phai;

  return beta_vector;
}

//called in insert_result
//calls     green_func[general_func]
//グリーン関数が格納された行列を別途計算
MatrixXd make_green_func_mat_for_insertion(vector<vector<double>> positions, vector<Source> Sources){
  int n_pos = positions.size();
  int n_src = Sources.size();
  MatrixXd green_func_mat;
  green_func_mat.resize(n_pos, n_src);
  for(int i=0; i<n_pos; ++i){
    for(int j=0; j<n_src; ++j){
      green_func_mat(i,j) = green_func(Sources[j].position[0], Sources[j].position[1], positions[i][0], positions[i][1]);
    }
  }
  return green_func_mat;
}

//called by *
//calls     make_green_func_mat_for_insertion
// 受信点に対応する結果を出力（実質ドット積を計算するだけ）
vector<double> insert_result(vector<vector<double>> positions, vector<Source> Sources, VectorXd beta_vector){
  int n_pos = positions.size();
  MatrixXd green_func_mat;
  VectorXd v_;

  green_func_mat = make_green_func_mat_for_insertion(positions, Sources);
  v_ = green_func_mat * beta_vector;

  vector<double> phai_for_pts(n_pos);
  for(int i = 0; i < n_pos; ++i){
    phai_for_pts[i] = v_(i);
  }

  return phai_for_pts;
}
