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

	void print_picture(); // Временный метод для дебага, потом нужно убрать.


private:
	double t;  // Суммарное время
	double dX; // Шаг по стержню
	int v_len, h_len;  // Число шагов у вертикальных и горионтальных стержней
	int v_num, h_num;  // Число вертикальных и горионтальных стержней
	double at;  // Коэффициент температуропроводности (Медь)
	double dt;  // Шаг по времени
	double r;  // Результирующий коэффициент (at*dt/dX^2)
	double Tcu; // Температура сверху.
	double Tcd; // Температура снизу.
	double Ts; // Изначальная температура системы.
	vector<int> non_intersecting_points; // В данных точках пересечение выключено (Отсчет идет от 1, слева направо, сверху вниз, то есть как буквы в книге).
	map<string, Point> result;  // Контейнер с данными, для отправки клиентскому приложению.
	vector<int> v_intersect{}; // Координаты вертикальных сержней
	vector<int> h_intersect{}; // Координаты горизонтальных сержней
	vector< vector<int> > intersect_2; // Координаты пересеченйи стержней
	vector< vector<vector<double> > > h_unknown;  // Матрицы горизонтальных стержней
	vector< vector<double> > h_last_temps;  // Предыдущие температуры в точках горизонтальных стержней
	vector< vector<double> > h_known;  // Известные данные горизонтальных стержней
	vector< vector<double> > h_new_temps;  // Новые значения температур в точках горизонтальных стержней
	vector< vector<vector<double> > > v_unknown;  // Матрицы горизонтальных стержней
	vector< vector<double> > v_last_temps;  // Предыдущие температуры в точках горизонтальных стержней
	vector< vector<double> > v_known;  // Известные данные горизонтальных стержней
	vector< vector<double> > v_new_temps;  // Новые значения температур в точках горизонтальных стержней


	vector<int> strToVec(string str);
	vector<double> gauss(vector< vector<double> > matrix, vector<double> col, int n);
	void init_vertical();
	void init_horizontal();
	void update_known_temps(vector<double> &known, vector<double> &last_temps, vector<double> new_temps, int length);
	void disconnect_points();
};

