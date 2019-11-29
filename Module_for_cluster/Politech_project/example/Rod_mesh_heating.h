#pragma once
#include <vector>
#include <deque>
#include "../Task.h"
#include "../Point.h"
using std::vector;
using std::deque;
using std::string;
using std::map;

class RodMeshHeating : public Task {
public:
	RodMeshHeating();
	~RodMeshHeating();
	void parseParam(map<string, string> mapParam);
	double getTime();
	map<string, Point> getPoints();
	void nextStep();

	void print_picture(); // Временный метод для дебага, потом нужно убрать.


private:
	double t;  // Суммарное время
	double dX; // Шаг по стержнюa
	int v_len, h_len;  // Число шагов у вертикальных и горионтальных стержней
	int v_num, h_num;  // Число вертикальных и горионтальных стержней
	double at_v, at_h;  // Коэффициент температуропроводности
	double dt;  // Шаг по времени
	double r_v, r_h;  // Результирующий коэффициент (at*dt/dX^2)
	double Tcu;  // Температура сверху.
	double Tcd;  // Температура снизу.
	double Ts;  // Изначальная температура системы.
	double T_error;  // Погрешность на концах горионтальных стержней.
	double Tm;  // Температура плавления.
	double Lp;  // Мощность лампочки.
	double radius_v, radius_h;  // Радиус стержня.
	double density_v, density_h;  // Плотность стержня.
	double specific_heat_v, specific_heat_h;  // Удельная теплоемкость стержня.
	double absorption_coefficient_v, absorption_coefficient_h; // Коэффициент поглащения стержня.
	vector<int> non_intersecting_points; // В данных точках пересечение выключено (Отсчет идет от 1, слева направо, сверху вниз, то есть как буквы в книге).
	map<string, Point> result;  // Контейнер с данными, для отправки клиентскому приложению.
	vector<int> v_intersect{}; // Координаты вертикальных сержней
	vector<int> h_intersect{}; // Координаты горизонтальных сержней
	vector< vector<int> > intersect_2; // Координаты пересеченйи стержней
	vector< deque<deque<double> > > h_unknown;  // Матрицы горизонтальных стержней
	vector< deque<double> > h_last_temps;  // Предыдущие температуры в точках горизонтальных стержней
	vector< deque<double> > h_known;  // Известные данные горизонтальных стержней
	vector< vector<double> > h_new_temps;  // Новые значения температур в точках горизонтальных стержней

	vector<int> extend_len; // Новые длины стержней.
	vector<int> extend;  // Удлинение стержней

	
	vector< vector<vector<double> > > v_unknown;  // Матрицы горизонтальных стержней
	vector< vector<double> > v_last_temps;  // Предыдущие температуры в точках горизонтальных стержней
	vector< vector<double> > v_known;  // Известные данные горизонтальных стержней
	vector< vector<double> > v_new_temps;  // Новые значения температур в точках горизонтальных стержней


	vector<int> strToVec(string str);
	void initial_values_calc();
	template < typename a_type, typename y_type>
	vector<double> gauss(a_type matrix, y_type col, int n);
	void infinite_slau(vector<double> &h_new_temps, deque<deque<double>> &unknown, deque<double> &known, deque<double> &last, double r, double error, int &extended, int &extended_len);
	template < typename T_1, typename T_2>
	void init(vector< T_2 > &unknown, vector< T_1 > &last_temps, vector< T_1 > &known, vector< vector<double> > &new_temps, int num, int len, double r, bool vertical = false);
	template < typename known_type, typename known_type_2>
	void update_known_temps(known_type &known, known_type_2 &last_temps, vector<double> new_temps, int length);
	void disconnect_points();
	void lamp();
};

