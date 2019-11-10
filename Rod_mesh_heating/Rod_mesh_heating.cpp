#include "pch.h"
#include <iostream>
#include <iomanip>      // std::setw


using namespace std;

double *gauss(double **matrix, double *col, int n)
{
	double *x, **a, *y, max;
	int k, index;
	const double eps = 0.0000000000000001;  // точность
	x = new double[n];
	y = new double[n];
	a = new double*[n];
	k = 0;


	for (int i = 0; i < n; i++)
	{
		a[i] = new double[n];
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
			return 0;
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


void init_horizontal(double ***z, double **y, double **new_temps, double **last_temsp, int number, int size, double ts, double r)
{
	for (int i = 0; i < number; i++)
	{
		y[i] = new double[size];
		z[i] = new double*[size];
		last_temsp[i] = new double[size];
		for (int j = 0; j < size; j++)
		{
			y[i][j] = ts;

			z[i][j] = new double[size];
			for (int k = 0; k < size; k++)
			{
				z[i][j][k] = 0;
			}

			z[i][j][j - 1] = -r;
			z[i][j][j] = r * 2 + 1;
			z[i][j][j + 1] = -r;

			last_temsp[i][j] = ts;
		}

		// Т.к. у боковых точек только 1 соседняя
		z[i][0][0] = r + 1;
		z[i][size - 1][size - 1] = r + 1;

		new_temps[i] = new double[size];
	}
}


void init_vertical(double ***z, double **y, double **new_temps, double **last_temsp, int number, int size, int tcu, int tcd, double ts, double r)
{
	for (int i = 0; i < number; i++)
	{
		y[i] = new double[size];
		z[i] = new double*[size];
		last_temsp[i] = new double[size];
		for (int j = 0; j < size; j++)
		{
			y[i][j] = ts;

			z[i][j] = new double[size];
			for (int k = 0; k < size; k++)
			{
				z[i][j][k] = 0;
			}

			z[i][j][j - 1] = -r;
			z[i][j][j] = r * 2 + 1;
			z[i][j][j + 1] = -r;

			last_temsp[i][j] = ts;
		}

		y[i][0] += tcu * r; // Добавляем граничные условия
		y[i][size - 1] += tcd * r; // Добавляем граничные условия

		new_temps[i] = new double[size];
	}
}


void print_picture(double **v_temps, double **h_temps, int v_len, int h_len, int v_num, int h_num, int intersect[2][5])
{
	for (int i = 0; i < v_len; i++) // столбцы
	{
		for (int j = 0; j < h_len; j++)  // строки
		{
			// TODO: Нужно заменить нормальным циклом.
			if (j == intersect[0][0])
				cout << setw(10) << v_temps[0][i];
			else if (j == intersect[0][1])
				cout << setw(10) << v_temps[1][i];
			else if (j == intersect[0][2])
				cout << setw(10) << v_temps[2][i];
			else if (j == intersect[0][3])
				cout << setw(10) << v_temps[3][i];

			else if (i == intersect[1][0])
				cout << setw(10) << h_temps[0][j];
			else if (i == intersect[1][1])
				cout << setw(10) << h_temps[1][j];
			else if (i == intersect[1][2])
				cout << setw(10) << h_temps[2][j];
			else if (i == intersect[1][3])
				cout << setw(10) << h_temps[3][j];
			else if (i == intersect[1][4])
				cout << setw(10) << h_temps[4][j];
			else
				cout << setw(10) << " ";
		}
		cout << endl;
	}
	cout << endl << endl;
}


void update_known_temps(double *known, double *last_temps, double *new_temps, int length)
{
	for (int i = 0; i < length; i++)
	{
		known[i] -= last_temps[i];
		known[i] += new_temps[i];
		last_temps[i] = new_temps[i];
	}
}


int main()
{
	setlocale(LC_ALL, "Russian");
	int v_len = 17, h_len = 13; // Число шагов у вертикальных и горионтальных стержней
	int v_num = 4, h_num = 5; // Число вертикальных и горионтальных стержней
	double **h_new_temps = new double*[h_num];  // Температуры в точках горизонтальных стержней
	double **h_last_temps = new double*[h_num];  // Предыдущие температуры в точках горизонтальных стержней
	double **h_known = new double*[h_num];  // Известные данные горизонтальных стержней
	double ***h_unknown = new double**[h_num];  // Матрицы горизонтальных стержней

	double **v_new_temps = new double*[v_num];  // Температуры в точках горизонтальных стержней
	double **v_last_temps = new double*[v_num];  // Предыдущие температуры в точках горизонтальных стержней
	double **v_known = new double*[v_num];  // Известные данные горизонтальных стержней
	double ***v_unknown = new double**[v_num];  // Матрицы горизонтальных стержней

	//Точки пересечения стержней
	int intersect[2][5] = {
		{2, 5, 7, 10, 0},// горизонталь
		{2, 5, 8, 11, 14} // вертикаль
	};

	//double at = 0.00008418; // Коэффициент температуропроводности (Алюминий)
	double at = 0.000111; // Коэффициент температуропроводности (Медь)
	double dt = 0.001; // Шаг по времени
	double lambda = 0.001; // Шаг по стержню
	double r = at * dt / lambda / lambda; // Результирующий коэффициент
	//double r = 0.000001; // Минимальный коэффициент, влияющий на результат
	//double r = 1;
	int steps_num = 50; // Число шагов вычислений.
	double tcu = 100; // Температура сверху.
	double tcd = 100; // Температура снизу.
	double ts = 10; // Изначальная температура системы.

	init_vertical(v_unknown, v_known, v_new_temps, v_last_temps, v_num, v_len, tcu, tcd, ts, r); // Инициализация матриц вертикальных стержней
	init_horizontal(h_unknown, h_known, h_new_temps, h_last_temps, h_num, h_len, ts, r); // Инициализация матриц горизонтальных стержней

	// Рассчет значений
	for (int i = 0; i < steps_num; i++)
	{
		// Рассчет температур вертикальных стержней.
		for (int j = 0; j < v_num; j++)
		{
			v_new_temps[j] = gauss(v_unknown[j], v_known[j], v_len);
		}

		// Перисвоение новых значений горизонтальным стержням(в точках пересечения).
		for (int k = 0; k < h_num; k++)
		{
			for (int j = 0; j < v_num; j++)
			{
				h_known[k][intersect[0][j]] = v_new_temps[j][intersect[1][k]];
			}
		}

		// Рассчет температур горизонтальных стержней.
		for (int j = 0; j < h_num; j++)
		{
			h_new_temps[j] = gauss(h_unknown[j], h_known[j], h_len);
		}

		// Перисвоение новых значений вертикальным стержням(в точках пересечения).
		for (int k = 0; k < v_num; k++)
		{
			for (int j = 0; j < h_num; j++)
			{
				v_new_temps[k][intersect[1][j]] = h_new_temps[j][intersect[0][k]];
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
		//print_picture(v_new_temps, h_new_temps, v_len, h_len, v_num, h_num, intersect);
	}

	print_picture(v_new_temps, h_new_temps, v_len, h_len, v_num, h_num, intersect);

	return 0;
}