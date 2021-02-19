#include <cmath>
#include <iostream>
#include "general_func.h"
#define PI 3.14159265359

double get_x(double xtx, double xrx){
  return xtx - xrx;
}

double get_y(double ytx, double yrx){
  return ytx - yrx;
}

//2点間の距離
double get_r(double xtx, double ytx, double xrx, double yrx){
  double x = get_x(xtx, xrx);
  double y = get_y(ytx, yrx);

  return sqrt(x*x + y*y);

}

//２次元ラプラス問題の解
double green_func(double xtx, double ytx, double xrx, double yrx){
  double r = get_r(xtx, ytx, xrx, yrx);

  //std::cout << "get_r(xtx, ytx, xrx, yrx)\t"<<get_r(0, 0, 3, 4)<<"\n";
  return -1. / 2./ PI * log(r);
}

//２次元ラプラス問題の勾配
double green_func_grad(double xtx, double ytx, double xrx, double yrx, double nx, double ny){
  double v,r;
  v = -nx * get_x(xtx, xrx) - ny * get_y(ytx, yrx);
  r = get_r(xtx, ytx, xrx, yrx);

  v = v/2./PI/r/r;

  return v;
}
