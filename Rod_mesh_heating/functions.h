#ifndef FUNCTION_H_INCLUDED_
#define FUNCTION_H_INCLUDED_
#include <vector>
using std::vector;

vector<double> gauss(vector< vector<double> > matrix, vector<double> col, int n);
void init_vertical(vector< vector<vector<double> > > &h_new_temps, vector< vector<double> > &v_known, double r, double tcu, double tcd);
void init_horizontal(vector< vector<vector<double> > > &h_new_temps, vector< vector<double> > &v_known, double r);
void print_picture(vector< vector<double> > v_temps, vector< vector<double> > h_temps, vector<int> v_intersect, vector<int> h_intersect);
void update_known_temps(vector<double> &known, vector<double> &last_temps, vector<double> new_temps, int length);
void disconnect_points(vector<int> v_intersect, vector<int> h_intersect, vector< vector<int> > &intersect_array, vector<int> non_intersecting_points);

#endif