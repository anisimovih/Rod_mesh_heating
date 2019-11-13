#pragma once
#include <vector>
#include "../Task.h"
#include "../Point.h"
using std::vector;
using std::string;
using std::map;

class RodMeshHeating : public Task{
public:
	RodMeshHeating();
    ~RodMeshHeating();
    void parseParam(map<string,string> mapParam);
    double getTime();
    map<string, Point> getPoints();
    void nextStep();

	void print_picture(); // ��������� ����� ��� ������, ����� ����� ������.


private:
	double t;  // ��������� �����
	double dX; // ��� �� �������
	int v_len, h_len;  // ����� ����� � ������������ � ������������� ��������
	int v_num, h_num;  // ����� ������������ � ������������� ��������
	double at;  // ����������� ���������������������� (����)
	double dt;  // ��� �� �������
	double r;  // �������������� ����������� (at*dt/dX^2)
	double Tcu; // ����������� ������.
	double Tcd; // ����������� �����.
	double Ts; // ����������� ����������� �������.
	vector<int> non_intersecting_points; // � ������ ������ ����������� ��������� (������ ���� �� 1, ����� �������, ������ ����, �� ���� ��� ����� � �����).
	map<string, Point> result;  // ��������� � �������, ��� �������� ����������� ����������.
	vector<int> v_intersect{}; // ���������� ������������ �������
	vector<int> h_intersect{}; // ���������� �������������� �������
	vector< vector<int> > intersect_2; // ���������� ����������� ��������
	vector< vector<vector<double> > > h_unknown;  // ������� �������������� ��������
	vector< vector<double> > h_last_temps;  // ���������� ����������� � ������ �������������� ��������
	vector< vector<double> > h_known;  // ��������� ������ �������������� ��������
	vector< vector<double> > h_new_temps;  // ����� �������� ���������� � ������ �������������� ��������
	vector< vector<vector<double> > > v_unknown;  // ������� �������������� ��������
	vector< vector<double> > v_last_temps;  // ���������� ����������� � ������ �������������� ��������
	vector< vector<double> > v_known;  // ��������� ������ �������������� ��������
	vector< vector<double> > v_new_temps;  // ����� �������� ���������� � ������ �������������� ��������


	vector<int> strToVec(string str);
	vector<double> gauss(vector< vector<double> > matrix, vector<double> col, int n);
	void init_vertical();
	void init_horizontal();
	void update_known_temps(vector<double> &known, vector<double> &last_temps, vector<double> new_temps, int length);
	void disconnect_points();
};

