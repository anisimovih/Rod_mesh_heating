#include "pch.h" // Предкомпилированные заголовки Visual Studio, по факту не нужны.
#include <iostream>
#include <vector>
#include "functions.h" 
using std::vector;


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

	disconnect_points(v_intersect, h_intersect, intersect_2, non_intersecting_points_v); // Создание массива пересейчений с учетом указанных непересекающихс точек.

	// Рассчет значений
	for (int i = 0; i < steps_num; i++)
	{
		// Рассчет температур вертикальных стержней.
		for (int j = 0; j < v_num; j++)
		{
			v_new_temps[j] = gauss(v_unknown[j], v_known[j], v_len);
		}

		for (int j = 0; j < h_num; j++)
		{
			for (int k = 0; k < v_num; k++)
			{
				if (intersect_2[k + j * v_num][0] != -1)
				{
					// Перисвоение новых значений горизонтальным стержням(в точках пересечения).
					h_known[j][intersect_2[k + j * v_num][1]] = v_new_temps[k][intersect_2[k + j * v_num][0]];
				}
			}

			// Рассчет температур горизонтальных стержней.
			h_new_temps[j] = gauss(h_unknown[j], h_known[j], h_len);
			for (int k = 0; k < v_num; k++)
			{
				if (intersect_2[k + j * v_num][0] != -1)
				{
					// Перисвоение новых значений вертикальным стержням(в точках пересечения).
					v_new_temps[k][intersect_2[k + j * v_num][0]] = h_new_temps[j][intersect_2[k + j * v_num][1]];
				}
			}
			// Обновление известных значений горизонтальных стержней.
			update_known_temps(h_known[j], h_last_temps[j], h_new_temps[j], h_len);
		}

		// Обновление известных значений горизонтальных стержней.
		for (int j = 0; j < v_num; j++)
		{
			update_known_temps(v_known[j], v_last_temps[j], v_new_temps[j], v_len);
		}

		// Отрисвока картинки на каждом шаге.
		//print_picture(v_new_temps, h_new_temps, v_len, h_len, v_num, h_num, v_intersect, h_intersect);
	}

	print_picture(v_new_temps, h_new_temps, v_intersect, h_intersect);
	
	return 0;
}