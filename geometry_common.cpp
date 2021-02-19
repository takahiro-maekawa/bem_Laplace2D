#include <vector>
#include <iostream>
#include <cmath>

#include "geometry_common.h"
#include "geometry.h"

//called in set_Nodes_to_Elements
//要素端の添字を取得([[0,1], [1,2]...[n-1,0]])
vector<vector<int>> prepare_indexes(int n_elem){
  vector<vector<int>> indexes(n_elem, vector<int>(2, 0));
  for (int i=0; i<n_elem; i++){
    for (int j=0; j<2; j++){
      indexes[i][j] = i+j;
    }
  }
  indexes[n_elem-1][1] = 0;

  return indexes;
}

//called in make_Elements
//calls     prepare_indexes
//要素クラスに要素端情報を挿入
void set_Nodes_to_Elements(vector<Element>& Elements, vector<Node> Nodes, int n_elem){
  
  vector<vector<int>> indexes = prepare_indexes(n_elem);

  for(int i=0; i<n_elem; ++i){
    for (int j=0; j<2; ++j){
      Elements[i].Nodes.resize(2);
      Elements[i].Nodes[j] = Nodes[indexes[i][j]];
    }
  }
}

//called in make_Elements
//境界条件を要素クラスに挿入
void set_boundaries_to_Elements(vector<Element>& Elements, vector<int> boundary_types,
  vector<double> boundary_values, int n_elem){
  for(int i=0; i<n_elem; ++i){
    Elements[i].boundary_type=boundary_types[i];
    Elements[i].boundary_value=boundary_values[i];
  }
}

//called in make_Elements
//idを要素クラスに挿入
void set_id_to_Elements(vector<Element>& Elements, int n_elem){
  for (int i = 0; i<n_elem; ++i){
    Elements[i].id = i;
  }
}

//called in make_Elements
//要素端の重心を要素の位置として定義
void set_position_of_Elements(vector<Element>& Elements, int n_elem){
  //n_node = Elements
  //cout << "n_elem\t" << n_elem << "\n";

  for(int i = 0; i<n_elem; ++i){
    int n_node = Elements[i].Nodes.size();
    int n_dim = Elements[i].Nodes[0].position.size();

    for(int i_dim=0; i_dim < n_dim; ++i_dim){
      double pos = 0;

      for (int i_node=0; i_node<n_node; ++i_node){
        pos = pos + Elements[i].Nodes[i_node].position[i_dim];
      }

      pos = pos / n_node;

      Elements[i].position[i_dim] = pos;

    }
  }
}

//called in make_Elements
//線要素に垂直なベクトルを取得
void set_l_and_normal_vector_of_Elements(vector<Element>& Elements, int n_elem){
  for(int i = 0; i<n_elem; ++i){
    double dx = Elements[i].Nodes[1].position[0] - Elements[i].Nodes[0].position[0];
    double dy = Elements[i].Nodes[1].position[1] - Elements[i].Nodes[0].position[1];

    double l = sqrt(dx*dx + dy*dy);
    Elements[i].normal_vector[0] = dy / l;
    Elements[i].normal_vector[1] = -dx / l;
    Elements[i].S = l;
  }
}

//called iby[geometry_sample]
//calls     set_Nodes_to_Elements, set_boundaries_to_Elements, set_id_to_Elements,
//          set_position_of_Elements, set_l_and_normal_vector_of_Elements
//初期化
void make_Elements(vector<Element>& Elements, vector<Node> Nodes, vector<int> boundary_types, vector<double> boundary_values){

  int n_elem =Elements.size();

  set_Nodes_to_Elements(Elements, Nodes, n_elem);

  set_boundaries_to_Elements(Elements, boundary_types, boundary_values, n_elem);

  set_id_to_Elements(Elements, n_elem);

  set_position_of_Elements(Elements, n_elem);

  set_l_and_normal_vector_of_Elements(Elements, n_elem);

}
