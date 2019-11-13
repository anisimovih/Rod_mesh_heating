#include "./example/Rod_mesh_heating.h"
#include <iostream>
using std::map;
using std::string;

int main()
{
	RodMeshHeating *test = new RodMeshHeating;
	map <string, string> params_map = { { "d_x", "0.05" },
										{ "v_len", "1" },
										{ "h_len", "1" },
										{ "v_num", "4" },
										{ "h_num", "5" },
										{ "at", "0.000111" },
										{ "dt", "1" },
										{ "Tcu", "100" },
										{ "Tcd", "100" },
										{ "Ts", "10" },
										{ "non_intersect", "1, 5, 9, 13"}
									  };
	test->parseParam(params_map);
	for (int i = 0; i < 50; i++)
	{
		test->nextStep();
	}
	// test->print_picture();  // Отрисовка картинки в консоли.
	test->getPoints();
	test->getTime();
}
