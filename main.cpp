#include <vector>
#include <fstream>
#include "geometry_sample.h"
#include "solve.h"

//　受信点(内部)を用意
vector<vector<double>> make_positions_sample(double r_min, double r_max, int n_r, int n_theta){
  vector<vector<double>> positions(n_r * n_theta, vector<double>(2));
  int index = 0;
  double r, theta;

  for (int i = 0; i< n_r; ++i){
    r = r_min + (r_max - r_min)/(n_r-1) * i;

    for (int j = 0; j< n_theta; ++j){
      theta = 0 + (PI / 2. - 0)/(n_theta-1) * j;
      positions[index] = {r * cos(theta), r * sin(theta)};
      index = index + 1;
    }
  }

  return positions;
}

//　結果をcsvに出力
void write_result_csv_sample(vector<vector<double>> positions,
   vector<double> results, int n_r, int n_theta, string filepath){
    int n_result = results.size();
    vector<vector<double>> XX(n_r, vector<double>(n_theta)), YY(n_r, vector<double>(n_theta)),  ZZ(n_r, vector<double>(n_theta));

    ofstream writing_file;
    writing_file.open(filepath, std::ios::out);

    writing_file<< "i,j,X,Y,Z\n";
    for (int index = 0; index < n_result; ++index){
      int i = index / n_theta;
      int j = index % n_theta;

      XX[i][j] = positions[index][0];
      YY[i][j] = positions[index][1];
      ZZ[i][j] = results[index];

      writing_file<< i <<","<<j<<","<<XX[i][j] <<","<<YY[i][j] <<","<<ZZ[i][j]<<"\n";
    }

}

//calls Element_geometry_include_sample[geometry_sample], Source_geometry_include_sample[geometry_sample],
//      solve_problem_and_get_beta_vector[solve]
//      make_positions_sample, insert_result[solve], write_result_csv_sample
int main(){
  // in geometry_sample.cpp
  vector<Element> Elements = Element_geometry_include_sample(2, 5, 4, 15);

  // in geometry_sample.cpp
  vector<Source> Sources = Source_geometry_include_sample(2, 5, 2, 2);

  //in solve.cpp
  VectorXd beta_vector = solve_problem_and_get_beta_vector(Elements, Sources);

  //sample receiver points
  double r_min = 2; double r_max = 5; int n_r = 100; int n_theta = 40;
  vector<vector<double>> positions = make_positions_sample(r_min, r_max, n_r, n_theta);
  vector<double> results = insert_result(positions, Sources, beta_vector);
  write_result_csv_sample(positions, results, n_r, n_theta, "test.csv");

  return 1;
}
