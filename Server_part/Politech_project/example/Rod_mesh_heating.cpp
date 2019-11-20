#include "Rod_mesh_heating.h"
#include <iostream>
#include <iomanip>      // std::setw
#include <cmath>		// sqrt
using std::pair;
using std::to_string;
using std::cout;
using std::endl;
using std::stod;
using std::setw;
using std::abs;


RodMeshHeating::RodMeshHeating() :Task() {
}


RodMeshHeating::~RodMeshHeating() {
	result.clear();
}


map<string, Point> RodMeshHeating::getPoints() {
	result.clear();

	for (int j = 0; j < v_num; j++)
	{
		for (int i = 0; i < v_len; i++)
		{
			Point point(v_intersect[j] * dX, i * dX, v_new_temps[j][i]);
			result.insert(pair<string, Point>(to_string(v_intersect[j] * dX) + to_string(i * dX), point));
		}
	}

	for (int j = 0; j < h_num; j++)
	{
		for (int i = 0; i < h_len; i++)
		{
			Point point(i * dX, h_intersect[j] * dX, h_new_temps[j][i]);
			result.insert(pair<string, Point>(to_string(i * dX) + to_string(h_intersect[j] * dX), point));
		}
	}

	return result;
}


double RodMeshHeating::getTime() {
	return t;
}


void RodMeshHeating::nextStep() {
	// Рассчет температур вертикальных стержней.
	for (int j = 0; j < v_num; j++)
	{
		v_new_temps[j] = gauss<vector<vector<double>>, vector<double>>(v_unknown[j], v_known[j], v_len);
	}
	// Перисвоение новых значений горизонтальным стержням(в точках пересечения).
	for (int j = 0; j < h_num; j++)
	{
		for (int k = 0; k < v_num; k++)
		{
			if (intersect_2[k + j * v_num][0] != -1)
			{
				h_known[j][intersect_2[k + j * v_num][1]] = v_new_temps[k][intersect_2[k + j * v_num][0]];
			}
		}
	}

	// Рассчет температур горизонтальных стержней.
	for (int j = 0; j < h_num; j++)
	{
		//h_new_temps[j] = gauss<deque<deque<double>>, deque<double>>(h_unknown[j], h_known[j], h_len);
		h_new_temps[j] = infinite_slau(h_unknown[j], h_known[j], h_len, r_h, Ts, T_error);
	}
	// Перисвоение новых значений вертикальным стержням(в точках пересечения).
	for (int j = 0; j < h_num; j++)
	{
		for (int k = 0; k < v_num; k++)
		{
			if (intersect_2[k + j * v_num][0] != -1)
			{
				v_new_temps[k][intersect_2[k + j * v_num][0]] = h_new_temps[j][intersect_2[k + j * v_num][1]];
			}
		}
	}

	lamp();  // Добавление воздействия лампочки.

	// Запись найденных значений вертикальных стержней в матрицу.
	for (int j = 0; j < v_num; j++)
	{
		update_known_temps<vector<double>>(v_known[j], v_last_temps[j], v_new_temps[j], v_len);
	}

	// Запись найденных значений горизонтальных стержней в матрицу.
	for (int j = 0; j < h_num; j++)
	{
		update_known_temps<deque<double>>(h_known[j], h_last_temps[j], h_new_temps[j], h_len);
	}

	t += dt;
}


void RodMeshHeating::parseParam(map<string, string> mapParam) {
	dX = stod(mapParam["dx"]);
	dt = stod(mapParam["dt"]);
	Lp = stod(mapParam["Мощность лампочки"]);
	Ts = stod(mapParam["Изначальная температура стержней"]);
	Tcu = stod(mapParam["Температура тела сверху"]);
	Tcd = stod(mapParam["Температура тела снизу"]);
	T_error = stod(mapParam["Погрешность на концах горионтальных стержней"]);
	v_len = stod(mapParam["Длина вертикальных стержней"]) / dX;
	h_len = stod(mapParam["Длина горизонтальных стержней"]) / dX;
	v_num = stod(mapParam["Число вертикальных стержней"]);
	h_num = stod(mapParam["Число горизонтальных стержней"]);
	at_v = stod(mapParam["Коэффициент температуропроводности вертикальных стержней"]);
	at_h = stod(mapParam["Коэффициент температуропроводности горизонтальных стержней"]);
	radius_v = stod(mapParam["Радиус вертикальных стержней"]);
	radius_h = stod(mapParam["Радиус горизонтальных стержней"]);
	density_v = stod(mapParam["Плотность вертикальных стержней"]);
	density_h = stod(mapParam["Плотность горизонтальных стержней"]);
	specific_heat_v = stod(mapParam["Теплоемкость вертикальных стержней"]);
	specific_heat_h = stod(mapParam["Теплоемкость горизонтальных стержней"]);
	absorption_coefficient_v = stod(mapParam["Коэффициент поглощения вертикальных стержней"]);
	absorption_coefficient_h = stod(mapParam["Коэффициент поглощения горизонтальных стержней"]);
	non_intersecting_points = strToVec(mapParam["Точки отсутствия пересечений"]);


	initial_values_calc();
}


// Парсинг строки в вектор.
vector<int> RodMeshHeating::strToVec(string str)
{
	vector<int> arr = {};

	string num = "";

	for (int i = 0; str[i] != '\0'; i++)
		if (str[i] != ',')
			num += str[i];
		else
		{
			arr.push_back(stoi(num));
			i++;
			num = "";
		}
	arr.push_back(stoi(num));

	return arr;
}


//  Рассет начальных значений, в зависимоти от переданных параметров.
void RodMeshHeating::initial_values_calc()
{
	r_v = at_v * dt / dX / dX;
	r_h = at_h * dt / dX / dX;

	// Рассчет оптимальных координат стержней, в зависимости от их длин.
	for (int i = 1; i < v_num + 1; i++)
	{
		v_intersect.push_back(int(h_len / (v_num + 1) * i) - 1);
	}

	for (int i = 1; i < h_num + 1; i++)
	{
		h_intersect.push_back(int(v_len / (h_num + 1)  * i) - 1);
	}


	init_vertical();
	init_horizontal();
	vector< vector<int> > intersect_2_new(v_num * h_num, vector<int>(2)); // Пересечения стержней
	intersect_2 = intersect_2_new;
	disconnect_points();
}


// Решение уравнений.
template < typename a_type, typename y_type>
vector<double> RodMeshHeating::gauss(a_type coefficient_matrix, y_type equations_values, int dimension)
{
	double max;
	int k, index;
	const double eps = 0.00001;  // точность
	vector<double> x(dimension);
	k = 0;

	while (k < dimension)
	{
		// Нормализация уравнений
		for (int i = k; i < dimension; i++)
		{
			double temp = coefficient_matrix[i][k];
			if (abs(temp) < eps) continue; // для нулевого коэффициента пропустить
			for (int j = 0; j < dimension; j++)
				coefficient_matrix[i][j] = coefficient_matrix[i][j] / temp;
			equations_values[i] = equations_values[i] / temp;
			if (i == k)  continue; // уравнение не вычитать само из себя
			for (int j = 0; j < dimension; j++)
				coefficient_matrix[i][j] = coefficient_matrix[i][j] - coefficient_matrix[k][j];
			equations_values[i] = equations_values[i] - equations_values[k];
		}
		k++;
	}
	// обратная подстановка
	for (k = dimension - 1; k >= 0; k--)
	{
		x[k] = equations_values[k];
		for (int i = 0; i < k; i++)
			equations_values[i] = equations_values[i] - coefficient_matrix[i][k] * x[k];
	}

	return x;
}

// Решение бесконечной системы уравнений(горизонтальные стержни) уравнений.
vector<double> RodMeshHeating::infinite_slau(deque<deque<double>> unknown, deque<double> known, int num, double r, double borders, double error)
{
	vector<double> result;
	vector<double> final_result(num);
	int count = 0, size = num;

	result = gauss<deque<deque<double>>, deque<double>>(unknown, known, num);
	while (result[num - 1] > borders + error || result[0] > borders + error)
	{
		count++;
		num += 2;

		deque<double> new_first_row(num, 0);
		new_first_row[0] = r + 1;
		new_first_row[1] = -r;
		unknown.push_front(new_first_row);
		known.push_front(10);

		unknown[1][0] = 2 * r + 1;
		unknown[1].push_front(-r);
		unknown[1].push_back(0);
		for (size_t i = 2; i < num - 2; i++)
		{
			unknown[i].push_front(0);
			unknown[i].push_back(0);
		}
		unknown[num - 2].push_front(0);
		unknown[num - 2][num - 2] = 2 * r + 1;
		unknown[num - 2].push_back(-r);

		deque<double> new_last_row(unknown.size() + 1, 0);
		new_last_row[num - 2] = -r;
		new_last_row[num - 1] = r + 1;
		unknown.push_back(new_last_row);
		known.push_back(10);

		result = gauss<deque<deque<double>>, deque<double>>(unknown, known, num);
	}

	int k = 0;
	for (int i = count; i < count + size; i++)
	{
		final_result[k] = result[i];
		k++;
	}

	return final_result;
}


// Инициализация данных у вертикальных стержней.
void RodMeshHeating::init_vertical()
{
	vector<deque<double>> abs;
	deque<double> abc;
	vector< vector<vector<double> > > v_unknown_new(v_num, vector< vector<double> >(v_len, vector<double>(v_len, 0)));  // Матрицы горизонтальных стержней
	vector< vector<double> > v_last_temps_new(v_num, vector<double>(v_len, Ts));  // Предыдущие температуры в точках горизонтальных стержней
	vector< vector<double> > v_known_new(v_num, vector<double>(v_len, Ts));  // Известные данные горизонтальных стержней
	vector< vector<double> > v_new_temps_new(v_num, vector<double>(v_len));  // Температуры в точках горизонтальных стержней

	v_unknown = v_unknown_new;
	v_last_temps = v_last_temps_new;
	v_known = v_known_new;
	v_new_temps = v_new_temps_new;

	for (size_t i = 0; i < v_new_temps.size(); i++)
	{
		for (size_t j = 0; j < v_new_temps[0].size(); j++)
		{
			if (j != 0)
				v_unknown[i][j][j - 1] = -r_v;
			v_unknown[i][j][j] = r_v * 2 + 1;
			if (j != v_new_temps[0].size() - 1)
				v_unknown[i][j][j + 1] = -r_v;
		}

		v_known[i][0] += Tcu * r_v; // Добавляем граничные условия
		v_known[i][v_known[0].size() - 1] += Tcd * r_v; // Добавляем граничные условия
	}
}


// Инициализация данных у горизонтальных стержней.
void RodMeshHeating::init_horizontal()
{
	vector< deque<deque<double> > > h_unknown_new(h_num, deque< deque<double> >(h_len, deque<double>(h_len, 0)));  // Матрицы горизонтальных стержней
	vector< vector<double> > h_last_temps_new(h_num, vector<double>(h_len, Ts));  // Предыдущие температуры в точках горизонтальных стержней
	vector< deque<double> > h_known_new(h_num, deque<double>(h_len, Ts));  // Известные данные горизонтальных стержней
	vector< vector<double> > h_new_temps_new(h_num, vector<double>(h_len));  // Температуры в точках горизонтальных стержней

	h_unknown = h_unknown_new;
	h_last_temps = h_last_temps_new;
	h_known = h_known_new;
	h_new_temps = h_new_temps_new;

	for (size_t i = 0; i < h_new_temps.size(); i++)
	{
		for (size_t j = 0; j < h_new_temps[0].size(); j++)
		{
			if (j != 0)
				h_unknown[i][j][j - 1] = -r_h;
			h_unknown[i][j][j] = r_h * 2 + 1;
			if (j != h_new_temps[0].size() - 1)
				h_unknown[i][j][j + 1] = -r_h;
		}

		// Т.к. у боковых точек только 1 соседняя
		h_unknown[i][0][0] = r_h + 1;
		h_unknown[i][h_new_temps[0].size() - 1][h_new_temps[0].size() - 1] = r_h + 1;
	}
}


// Вывод картинки в консоль.
void RodMeshHeating::print_picture()
{
	for (size_t i = 0; i < v_new_temps[0].size(); i++) // столбцы
	{
		for (size_t j = 0; j < h_new_temps[0].size(); j++)  // строки
		{
			bool printed = false;
			for (size_t k = 0; k < v_intersect.size(); k++)
			{
				if (j == v_intersect[k])
				{
					cout << setw(10) << v_new_temps[k][i];
					printed = true;
					break;
				}
			}
			if (printed == false)
			{
				for (size_t k = 0; k < h_intersect.size(); k++)
				{
					if (i == h_intersect[k])
					{
						cout << setw(10) << h_new_temps[k][j];
						printed = true;
						break;
					}
				}
			}
			if (printed == false)
				cout << setw(10) << " ";
		}
		cout << endl;
	}
	cout << endl << endl;
}


// Обновление полученных температур.
template < typename known_type>
void RodMeshHeating::update_known_temps(known_type &known, vector<double> &last_temps, vector<double> new_temps, int length)
{
	for (int i = 0; i < length; i++)
	{
		known[i] -= last_temps[i];
		known[i] += new_temps[i];
		last_temps[i] = new_temps[i];
	}
}


//Формирование массива с координатами пересеченйи стержней.
void RodMeshHeating::disconnect_points()
{
	int non_intersecting_points_num = non_intersecting_points.size();

	for (int j = 0; j < h_num; j++)
	{
		for (int i = 0; i < v_num; i++)
		{
			bool non_intersect = false;
			for (int k = 0; k < non_intersecting_points_num; k++)
			{
				if (i + j * v_num == non_intersecting_points[k] - 1)
					non_intersect = true;
			}
			if (non_intersect == false)
			{
				intersect_2[i + j * v_num][0] = h_intersect[j];
				intersect_2[i + j * v_num][1] = v_intersect[i];
			}
			else
			{
				intersect_2[i + j * v_num][0] = -1;
				intersect_2[i + j * v_num][1] = -1;
			}
		}
	}
}


void RodMeshHeating::lamp()
{
	double dT = Lp / (2 * 3.14 * 3.14 * radius_v * radius_v * radius_v * density_v * dX * dt * specific_heat_v) * absorption_coefficient_v;

	for (int i = 0; i < h_num; i++)
	{
		for (int j = 0; j < h_len; j++)
		{
			double x = abs(j - (h_len - 1) / 2);
			double y = abs(h_intersect[i] - (h_len - 1) / 2);
			double length = sqrt(x * x + y * y);
			if (length != 0)
				h_new_temps[i][j] += dT / length;
		}
	}

	dT = Lp / (2 * 3.14 * 3.14 * radius_h * radius_h * radius_h * density_h * dX * dt * specific_heat_h) * absorption_coefficient_h;

	for (int i = 0; i < v_num; i++)
	{
		for (int j = 0; j < v_len; j++)
		{
			double x = abs(v_intersect[i] - (v_len - 1) / 2);
			double y = abs(j - (v_len - 1) / 2);
			double length = sqrt(x * x + y * y);
			if (length != 0)
				v_new_temps[i][j] += dT / length;
		}
	}
}
