#include "pch.h" // Предкомпилированные заголовки Visual Studio, по факту не нужны.
#include <iostream>
#include <iomanip>      // std::setw
#include <vector>


using namespace std;


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


int main()
{
	setlocale(LC_ALL, "Russian");
	int v_len = 17, h_len = 13; // Число шагов у вертикальных и горионтальных стержней
	int v_num = 4, h_num = 5; // Число вертикальных и горионтальных стержней
	//double at = 0.00008418; // Коэффициент температуропроводности (Алюминий)
	double at = 0.000111; // Коэффициент температуропроводности (Медь)
	double dt = 0.001; // Шаг по времени
	double lambda = 0.001; // Шаг по стержню
	//double r = at * dt / lambda / lambda; // Результирующий коэффициент
	//double r = 0.000001; // Минимальный коэффициент, влияющий на результат
	double r = 0.5;
	int steps_num = 50; // Число шагов вычислений.
	double tcu = 100; // Температура сверху.
	double tcd = 100; // Температура снизу.
	double ts = 10; // Изначальная температура системы.

	vector< vector<vector<double> > > h_unknown(h_num, vector< vector<double> >(h_len, vector<double>(h_len, 0)));  // Матрицы горизонтальных стержней
	vector< vector<double> > h_last_temps(h_num, vector<double>(h_len, ts));  // Предыдущие температуры в точках горизонтальных стержней
	vector< vector<double> > h_known(h_num, vector<double>(h_len, ts));  // Известные данные горизонтальных стержней
	vector< vector<double> > h_new_temps(h_num, vector<double>(h_len));  // Температуры в точках горизонтальных стержней

	vector< vector<vector<double> > > v_unknown(v_num, vector< vector<double> >(v_len, vector<double>(v_len, 0)));  // Матрицы горизонтальных стержней
	vector< vector<double> > v_last_temps(v_num, vector<double>(v_len, ts));  // Предыдущие температуры в точках горизонтальных стержней
	vector< vector<double> > v_known(v_num, vector<double>(v_len, ts));  // Известные данные горизонтальных стержней
	vector< vector<double> > v_new_temps(v_num, vector<double>(v_len));  // Температуры в точках горизонтальных стержней

	vector<int> v_intersect{ 2, 5, 7, 10 }; // Координаты вертикальных сержней
	vector<int> h_intersect{ 2, 5, 8, 11, 14 }; // Координаты горизонтальных сержней

	vector< vector<int> > intersect_2(v_num * h_num, vector<int>(2)); // Пересечения стержней
	vector<int> non_intersecting_points_v{ 1, 5, 9, 13 }; // В данных точках пересечение выключено (Отсчет идет от 1, слева направо, сверху вниз, то есть как буквы в книге).

	init_vertical(v_unknown, v_known, r, tcu, tcd); // Инициализация матриц вертикальных стержней
	init_horizontal(h_unknown, h_known, r); // Инициализация матриц горизонтальных стержней

	disconnect_points(v_intersect, h_intersect, intersect_2, non_intersecting_points_v);

	// Рассчет значений
	for (int i = 0; i < steps_num; i++)
	{
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

		// Отрисвока картинки на каждом шаге.
		//print_picture(v_new_temps, h_new_temps, v_len, h_len, v_num, h_num, v_intersect, h_intersect);
	}

	print_picture(v_new_temps, h_new_temps, v_intersect, h_intersect);
	
	return 0;
}