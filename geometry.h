#ifndef INCLUDED_GEOMETRY
#define INCLUDED_GEOMETRY

#define PI 3.14159265359

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Core>

using namespace std;
using namespace Eigen;

class Node{
public:
  int id;
  vector<double> position = {0.0, 0.0};

};

class Element{
public:
  int id;
  vector<Node> Nodes;
  vector<double> position = {0.0, 0.0};
  vector<double> normal_vector = {0.0, 1.0};
  int boundary_type = 0; //ディリクレなら0, ノイマンなら1
  double boundary_value = 0.0;
  double S = 1.;

};

class Source{
public:
  int id;
  vector<double> position = {0.0, 0.0};

};

#endif
