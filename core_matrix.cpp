#include <vector>
#include <random>
#include <iostream>
#include "geometry.h"
#include "core_matrix.h"
#include "general_func.h"
#include "integralGL.h"

using namespace std;

// new added
// called in make_Gmat
//対角成分に微小値を足す
void avoid_diverge_hack(MatrixXd& M){
  int n = M.innerSize();
  double eps =1e-10;

  for (int i=0; i<n; ++i){
    M(i,i) = M(i,i) + eps;
  }
}

//called in make_coeff_matrix_proto
//[L]の計算 (文献参照)
MatrixXd make_Lmat(vector<Element> Elements, int n_elem){

  MatrixXd L_mat(n_elem,n_elem);

  for (int i =0; i < n_elem; ++i){
    for (int j =0; j < n_elem; ++j){
      if (i == j){
        L_mat(i,j) = Elements[i].S;
      }else{
        L_mat(i,j) = 0.0;
      }
    }
  }
  return L_mat;
}

//called in make_coeff_matrix_proto
//calls return_eta_weights[integralGL], insert_1d_position[integralGL], green_func[general_func], avoid_diverge_hack
//[G]の計算　(文献参照)
MatrixXd make_Gmat(vector<Element> Elements, vector<Source> Sources, int n_elem, int n_source){
  MatrixXd G_mat(n_source, n_elem);
  vector<double> etas, weights;
  int npts_for_integral = 3;

  return_eta_weights(etas, weights, npts_for_integral);

  for (int j =0; j < n_source; ++j){
    double xtx = Sources[j].position[0]; double ytx = Sources[j].position[1];

    for (int m =0; m < n_elem; ++m){
      vector<double> xs_elem_node(npts_for_integral), ys_elem_node(npts_for_integral);
      insert_1d_position(xs_elem_node, Elements[m].Nodes[0].position[0], Elements[m].Nodes[1].position[0], etas);
      insert_1d_position(ys_elem_node, Elements[m].Nodes[0].position[1], Elements[m].Nodes[1].position[1], etas);

      G_mat(j,m) = 0;

      for (int ipt=0; ipt < npts_for_integral; ++ipt){
        double xm_ipt = xs_elem_node[ipt]; double ym_ipt = ys_elem_node[ipt];
        double v = green_func(xtx, ytx, xm_ipt, ym_ipt);

        G_mat(j,m) += v * weights[ipt];
      }
      G_mat(j,m) = G_mat(j,m) * Elements[m].S / 2.0;
    }
  }

  avoid_diverge_hack(G_mat);
  return G_mat;
}

//called in make_coeff_matrix_proto
//calls return_eta_weights[integralGL], insert_1d_position[integralGL], green_func[general_func], green_func_grad[general_func]
//[F]の計算　(文献参照)
MatrixXd make_Fmat(vector<Element> Elements, vector<Source> Sources, int n_elem, int n_source){
  MatrixXd F_mat(n_source, n_source);
  vector<double> etas, weights;
  int npts_for_integral = 3;

  return_eta_weights(etas, weights, npts_for_integral);

  for (int j =0; j < n_source; ++j){
    double xtx1 = Sources[j].position[0]; double ytx1 = Sources[j].position[1];

    for (int k =0; k < n_source; ++k){
      double xtx2 = Sources[k].position[0]; double ytx2 = Sources[k].position[1];

      F_mat(j,k) = 0;

      for (int m =0; m < n_elem; ++m){
        double nxm =  Elements[m].normal_vector[0]; double nym =  Elements[m].normal_vector[1];

        vector<double> xs_elem_node(npts_for_integral), ys_elem_node(npts_for_integral);
        insert_1d_position(xs_elem_node, Elements[m].Nodes[0].position[0], Elements[m].Nodes[1].position[0], etas);
        insert_1d_position(ys_elem_node, Elements[m].Nodes[0].position[1], Elements[m].Nodes[1].position[1], etas);

        double integral_v = 0;

        for (int ipt=0; ipt < npts_for_integral; ++ipt){
          double xm_ipt = xs_elem_node[ipt]; double ym_ipt = ys_elem_node[ipt];

          double v1 = green_func(xtx1, ytx1, xm_ipt, ym_ipt) * green_func_grad(xtx2, ytx2, xm_ipt, ym_ipt, nxm, nym);
          double v2 = green_func(xtx2, ytx2, xm_ipt, ym_ipt) * green_func_grad(xtx1, ytx1, xm_ipt, ym_ipt, nxm, nym);

          integral_v += (v1 + v2) * weights[ipt];

        }
        integral_v = integral_v * Elements[m].S / 2.0;

        F_mat(j,k) += integral_v;
      }

      F_mat(j,k) = F_mat(j,k) / 2.0;
    }
  }
  return F_mat;
}

//called in  Culc_K_mat_from_G_L_F
//[R]の計算　(文献参照)
MatrixXd make_R_mat_from_G_L(MatrixXd G_mat, MatrixXd L_mat){
  MatrixXd R_mat;
  MatrixXd G_mat_T = G_mat.transpose();
  MatrixXd G_mat_T_inv = G_mat_T.inverse();

  R_mat = G_mat_T_inv * L_mat;
  return R_mat;
}

//called in   Culc_K_mat_from_G_L_F
//calls       make_R_mat_from_G_L
//[K]の計算　(文献参照)
void Culc_K_mat_from_G_L_F(MatrixXd& K_mat, MatrixXd G_mat, MatrixXd L_mat,MatrixXd F_mat){

  MatrixXd R_mat = make_R_mat_from_G_L(G_mat, L_mat);

  //K_mat = R_mat_T * F_mat; written in the document but maybe wrong
  K_mat = F_mat* R_mat;
}

//called by main()
//calls make_Lmat, make_Gmat, make_Fmat, Culc_K_mat_from_G_L_F
//計算に必要な行列を用意
void make_coeff_matrix_proto(MatrixXd& K_mat, MatrixXd& G_mat, MatrixXd& R_mat, vector<Element> Elements, vector<Source> Sources){
  int n_elem = Elements.size();
  int n_source = Sources.size();

  MatrixXd L_mat = make_Lmat(Elements, n_elem);
  G_mat = make_Gmat(Elements, Sources, n_elem, n_source);
  MatrixXd F_mat = make_Fmat(Elements, Sources, n_elem,n_source);
  R_mat = make_R_mat_from_G_L(G_mat, L_mat);

  Culc_K_mat_from_G_L_F(K_mat,G_mat,L_mat,F_mat);

}

//called in   transform_matrix
//ディレクレ条件及びノイマン条件に対応する要素のインデックスを取得
void make_boundary_indexes(vector<int> &Dirichlet_indexes, vector<int> &Neumann_indexes,
                vector<Element> Elements){
  int n_elem = Elements.size();

  for (int i_elem=0; i_elem < n_elem; ++i_elem){
    if(Elements[i_elem].boundary_type == 0){
      Dirichlet_indexes.push_back(i_elem);
    }else{
      Neumann_indexes.push_back(i_elem);
    }
  }

}

//called in   transform_K_G_to_A_B
//行列[K], [G]を入れ替えたりして、別の行列を作成する関数
void transform_matrix_swap_boundary(MatrixXd& M, MatrixXd mat1, MatrixXd mat2,
   vector<int> Dirichlet_indexes, vector<int> Neumann_indexes){
  int n_elem =  Dirichlet_indexes.size() + Neumann_indexes.size();
  int n_Dir_elem =  Dirichlet_indexes.size();

  for (int i = 0; i<n_elem; ++i){
    for (int j = 0; j<n_elem; ++j){
      if(j<n_Dir_elem){
        int index = Dirichlet_indexes[j];

        M(i,j) = mat1(i, index);

      }else{
        int index = Neumann_indexes[j - n_Dir_elem];
        M(i,j) =-1 * mat2(i, index);
      }
    }
  }

}

//called in   transform_matrix
//calls       transform_matrix_swap_boundary
//[K], [G]から連立方程式に利用できる行列を作る
void transform_K_G_to_A_B(MatrixXd& A, MatrixXd& B, MatrixXd K_mat,
   MatrixXd G_mat, vector<int> Dirichlet_indexes, vector<int> Neumann_indexes){

  transform_matrix_swap_boundary(A, G_mat, K_mat, Dirichlet_indexes, Neumann_indexes);
  transform_matrix_swap_boundary(B, K_mat, G_mat, Dirichlet_indexes, Neumann_indexes);

}

//called in  transform_matrix
//ディレクレ、ノイマン条件に対応する値を格納
VectorXd make_boundary_value_v(MatrixXd B, vector<Element> Elements, vector<int> Dirichlet_indexes, vector<int> Neumann_indexes){
  VectorXd b_value_vec(Elements.size());
  int n_elem =  Elements.size();
  int index = 0;
  for(int i=0; i<Dirichlet_indexes.size(); ++i){
    b_value_vec(index) = Elements[Dirichlet_indexes[i]].boundary_value;
    index += 1;
  }

  for(int i=0; i<Neumann_indexes.size(); ++i){
    b_value_vec(index) = Elements[Neumann_indexes[i]].boundary_value;
    index += 1;
  }

  return b_value_vec;

}

//called by main()
//calls       make_boundary_indexes, transform_K_G_to_A_B, make_b_value_v
//用意した行列から連立方程式を作成する関数
void transform_matrix(MatrixXd& A, VectorXd& b, VectorXd& b_value_vec, vector<int>& Dirichlet_indexes, vector<int>& Neumann_indexes,
  MatrixXd K_mat, MatrixXd G_mat, vector<Element> Elements){
  //vector<int> Dirichlet_indexes, Neumann_indexes;

  int n_elem = Elements.size();

  A.resize(n_elem, n_elem);
  MatrixXd B = A;

  make_boundary_indexes(Dirichlet_indexes, Neumann_indexes, Elements);
  transform_K_G_to_A_B(A, B, K_mat, G_mat, Dirichlet_indexes, Neumann_indexes);

  b_value_vec = make_boundary_value_v(B, Elements, Dirichlet_indexes, Neumann_indexes);
  b = B*b_value_vec;
}
