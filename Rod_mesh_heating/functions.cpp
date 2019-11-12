﻿#include "functions.h"
#include <iostream>
#include <iomanip>      // std::setw
#include <vector>
#include "functions.h" 
using std::cout;
using std::endl;
using std::vector;
using std::setw;

vector<double> gauss(vector< vector<double> > matrix, vector<double> col, int n)
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
			// нет ненулевых диагональных элементов
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

void init_vertical(vector< vector<vector<double> > > &h_new_temps, vector< vector<double> > &v_known, double r, double tcu, double tcd)
{
	for (size_t i = 0; i < h_new_temps.size(); i++)
	{
		for (size_t j = 0; j < h_new_temps[0].size(); j++)
		{
			if (j != 0)
				h_new_temps[i][j][j - 1] = -r;
			h_new_temps[i][j][j] = r * 2 + 1;
			if (j != h_new_temps[0].size() - 1)
				h_new_temps[i][j][j + 1] = -r;
		}

		v_known[i][0] += tcu * r; // Добавляем граничные условия
		v_known[i][v_known[0].size() - 1] += tcd * r; // Добавляем граничные условия
	}
}

void init_horizontal(vector< vector<vector<double> > > &h_new_temps, vector< vector<double> > &v_known, double r)
{
	for (size_t i = 0; i < h_new_temps.size(); i++)
	{
		for (size_t j = 0; j < h_new_temps[0].size(); j++)
		{
			if (j != 0)
				h_new_temps[i][j][j - 1] = -r;
			h_new_temps[i][j][j] = r * 2 + 1;
			if (j != h_new_temps[0].size() - 1)
				h_new_temps[i][j][j + 1] = -r;
		}

		// Т.к. у боковых точек только 1 соседняя
		h_new_temps[i][0][0] = r + 1;
		h_new_temps[i][h_new_temps[0].size() - 1][h_new_temps[0].size() - 1] = r + 1;
	}
}

void print_picture(vector< vector<double> > v_temps, vector< vector<double> > h_temps, vector<int> v_intersect, vector<int> h_intersect)
{
	for (size_t i = 0; i < v_temps[0].size(); i++) // столбцы
	{
		for (size_t j = 0; j < h_temps[0].size(); j++)  // строки
		{
			bool printed = false;
			for (size_t k = 0; k < v_intersect.size(); k++)
			{
				if (j == v_intersect[k])
				{
					cout << setw(10) << v_temps[k][i];
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
						cout << setw(10) << h_temps[k][j];
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
}

void update_known_temps(vector<double> &known, vector<double> &last_temps, vector<double> new_temps, int length)
{
	for (int i = 0; i < length; i++)
	{
		known[i] -= last_temps[i];
		known[i] += new_temps[i];
		last_temps[i] = new_temps[i];
	}
}

void disconnect_points(vector<int> v_intersect, vector<int> h_intersect, vector< vector<int> > &intersect_array, vector<int> non_intersecting_points)
{
	int h_num = h_intersect.size();
	int v_num = v_intersect.size();
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
				intersect_array[i + j * v_num][0] = h_intersect[j];
				intersect_array[i + j * v_num][1] = v_intersect[i];
			}
			else
			{
				intersect_array[i + j * v_num][0] = -1;
				intersect_array[i + j * v_num][1] = -1;
			}
		}
	}
}