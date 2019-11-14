#include "Rod_mesh_heating.h"
#include <iostream>
#include <vector>
#include <iomanip>      // std::setw
using std::map;
using std::pair;
using std::to_string;
using std::string;
using std::cout;
using std::endl;
using std::stod;
using std::setw;
using std::vector;


RodMeshHeating::RodMeshHeating():Task(){
}


RodMeshHeating::~RodMeshHeating(){
    result.clear();
}


map<string, Point> RodMeshHeating::getPoints(){
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


double RodMeshHeating::getTime(){
    return t;
}


void RodMeshHeating::nextStep(){
	// Рассчет температур вертикальных стержней.
	for (int j = 0; j < v_num; j++)
	{
		v_new_temps[j] = gauss(v_unknown[j], v_known[j], v_len);
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
		h_new_temps[j] = gauss(h_unknown[j], h_known[j], h_len);
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

	// Обновление известных значений горизонтальных стержней.
	for (int j = 0; j < v_num; j++)
	{
		update_known_temps(v_known[j], v_last_temps[j], v_new_temps[j], v_len);
	}
	// Обновление известных значений горизонтальных стержней.
	for (int j = 0; j < h_num; j++)
	{
		update_known_temps(h_known[j], h_last_temps[j], h_new_temps[j], h_len);
	}
	t += dt;
}


void RodMeshHeating::parseParam(map<string, string> mapParam){
	dX = stod(mapParam["d_x"]);
	v_len = stod(mapParam["v_len"]) / dX;
	h_len = stod(mapParam["h_len"]) / dX;
	v_num = stod(mapParam["v_num"]);
	h_num = stod(mapParam["h_num"]);
	at = stod(mapParam["at"]);
	dt = stod(mapParam["dt"]);
	Tcu = stod(mapParam["Tcu"]);
	Tcd = stod(mapParam["Tcd"]);
	Ts = stod(mapParam["Ts"]);
	non_intersecting_points = strToVec(mapParam["non_intersect"]);

	r = at * dt / dX / dX;

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


// Решение уравнений.
vector<double> RodMeshHeating::gauss(vector< vector<double> > matrix, vector<double> col, int n)
{
	double max;
	int k, index;
	const double eps = 0.0000000000000001;  // точность
	vector<double> x(n);
	vector<double> y(n);
	vector< vector<double> > a(n, vector<double>(n));
	k = 0;

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			a[i][j] = matrix[i][j];
		}
	}

	for (int i = 0; i < n; i++)
	{
		y[i] = col[i];
	}

	while (k < n)
	{
		// Поиск строки с максимальным a[i][k]
		max = abs(a[k][k]);
		index = k;
		for (int i = k + 1; i < n; i++)
		{
			if (abs(a[i][k]) > max)
			{
				max = abs(a[i][k]);
				index = i;
			}
		}
		// Перестановка строк
		if (max < eps)
		{
			// нет ненулевых диагональных элементов, по идее не должны сюда попадать.
			cout << "Решение получить невозможно из-за нулевого столбца ";
			cout << index << " матрицы A" << endl;
			return x;
		}
		for (int j = 0; j < n; j++)
		{
			double temp = a[k][j];
			a[k][j] = a[index][j];
			a[index][j] = temp;
		}
		double temp = y[k];
		y[k] = y[index];
		y[index] = temp;
		// Нормализация уравнений
		for (int i = k; i < n; i++)
		{
			double temp = a[i][k];
			if (abs(temp) < eps) continue; // для нулевого коэффициента пропустить
			for (int j = 0; j < n; j++)
				a[i][j] = a[i][j] / temp;
			y[i] = y[i] / temp;
			if (i == k)  continue; // уравнение не вычитать само из себя
			for (int j = 0; j < n; j++)
				a[i][j] = a[i][j] - a[k][j];
			y[i] = y[i] - y[k];
		}
		k++;
	}
	// обратная подстановка
	for (k = n - 1; k >= 0; k--)
	{
		x[k] = y[k];
		for (int i = 0; i < k; i++)
			y[i] = y[i] - a[i][k] * x[k];
	}

	return x;
}


// Инициализация данных у вертикальных стержней.
void RodMeshHeating::init_vertical()
{
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
				v_unknown[i][j][j - 1] = -r;
			v_unknown[i][j][j] = r * 2 + 1;
			if (j != v_new_temps[0].size() - 1)
				v_unknown[i][j][j + 1] = -r;
		}

		v_known[i][0] += Tcu * r; // Добавляем граничные условия
		v_known[i][v_known[0].size() - 1] += Tcd * r; // Добавляем граничные условия
	}
}


// Инициализация данных у горизонтальных стержней.
void RodMeshHeating::init_horizontal()
{
	vector< vector<vector<double> > > h_unknown_new(h_num, vector< vector<double> >(h_len, vector<double>(h_len, 0)));  // Матрицы горизонтальных стержней
	vector< vector<double> > h_last_temps_new(h_num, vector<double>(h_len, Ts));  // Предыдущие температуры в точках горизонтальных стержней
	vector< vector<double> > h_known_new(h_num, vector<double>(h_len, Ts));  // Известные данные горизонтальных стержней
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
				h_unknown[i][j][j - 1] = -r;
			h_unknown[i][j][j] = r * 2 + 1;
			if (j != h_new_temps[0].size() - 1)
				h_unknown[i][j][j + 1] = -r;
		}

		// Т.к. у боковых точек только 1 соседняя
		h_unknown[i][0][0] = r + 1;
		h_unknown[i][h_new_temps[0].size() - 1][h_new_temps[0].size() - 1] = r + 1;
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
void RodMeshHeating::update_known_temps(vector<double> &known, vector<double> &last_temps, vector<double> new_temps, int length)
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
