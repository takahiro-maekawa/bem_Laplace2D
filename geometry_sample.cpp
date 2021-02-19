#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include "geometry_sample.h"
#include "geometry_common.h"

#define PI 3.14159265359

using namespace std;

//called by prepare_Nodes, Source_geometry_include_sample
//四分円境界に対応する位置を取得
vector<vector<double>> quarter_positions(double a, double b){
  vector<vector<double>> positions(24, vector<double>(2));

  int index = 0;
  double theta;
  //内側の弧
  for(int i=0; i<5; ++i){
    theta = 90. * (1. - (double)i / 4.) / 180. * PI;
    positions[index]= {a * cos(theta), a * sin(theta)};
    index += 1;
  }

  //四分円の下側境界
  for(int i=0; i<5; ++i){
    positions[index] = {a + (b - a)/ 6. * ((double)i + 1) , 0.0};
    index += 1;
  }

  //外側の弧
  for(int i=0; i<9; ++i){
    theta = 90. * (double)i / 8. / 180. * PI;
    positions[index]= {b * cos(theta), b * sin(theta)};
    index += 1;
  }

  //四分円の左側境界
  for(int i=0; i<5; ++i){
    positions[index]= {0, b + (a - b)/ 6. * ((double)i + 1)};
    index += 1;
  }

  return positions;
}


//called in  Element_geometry_include_sample
//calls      quarter_positions
//要素端を用意
void prepare_Nodes(vector<Node> & Nodes, double a, double b){
  Nodes.resize(24);

  vector<vector<double>> positions = quarter_positions(a, b);

  for (int i = 0; i<Nodes.size(); ++i){
    Nodes[i].position = positions[i];
    Nodes[i].id = i;
  }
}

//called in  Element_geometry_include_sample
//境界条件の設定
void prepare_boundary_info(vector<int> &boundary_types, vector<double> &boundary_values, double phai_a, double phai_b){
  boundary_types.resize(24);
  boundary_values.resize(24);

  int index = 0;

  for(int i=0; i<4; ++i){
    boundary_types[index] = 0;
    boundary_values[index] = phai_a;
    index += 1;
  }

  for(int i=0; i<6; ++i){
    boundary_types[index] = 1;
    boundary_values[index] = 0.;
    index += 1;
  }
  for(int i=0; i<8; ++i){
    boundary_types[index] = 0;
    boundary_values[index] = phai_b;
    index += 1;
  }
  for(int i=0; i<6; ++i){
    boundary_types[index] = 1;
    boundary_values[index] = 0.;
    index += 1;
  }
}

//called by   main()
//calls       prepare_Nodes, prepare_boundary_info, make_Elements
//            all above in [geometry_common.cpp]
//サンプル問題の要素を用意
vector<Element> Element_geometry_include_sample(double a, double b, double phai_a, double phai_b){
  vector<Element> Elements;

  vector<Node> Nodes;
  vector<int> boundary_types;
  vector<double> boundary_values;

  prepare_Nodes(Nodes, a, b);
  prepare_boundary_info(boundary_types, boundary_values, phai_a, phai_b);

  Elements.resize(24);

  make_Elements(Elements, Nodes, boundary_types, boundary_values);

  return Elements;
}

//called by   main()
//calls   quarter_positions
//サンプル問題の外部電荷を用意
vector<Source> Source_geometry_include_sample(double a, double b, double offset_c, double offset_r){
  vector<Source> Sources;
  Sources.resize(24);

  vector<int> v(24);// = {0, 1, 2, 3, 4, 5};
  for (int i=0; i<24; ++i){
    v[i] = i;
  }

  vector<vector<double>> positions = quarter_positions(a + sqrt(offset_c) * offset_c - offset_r, b + offset_c +  offset_r);
  for (int i = 0; i<Sources.size(); ++i){
    Sources[i].position ={positions[i][0] - offset_c, positions[i][1] - offset_c};
  }
  Sources[5].position = {Sources[5].position[0] - 0.2, Sources[5].position[1] - 0.15 };

  return Sources;
}
