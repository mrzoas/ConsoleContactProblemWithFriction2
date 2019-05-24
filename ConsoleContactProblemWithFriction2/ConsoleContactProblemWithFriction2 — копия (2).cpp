// ConsoleContactProblemWithFriction2.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
// Эта программа начата 2019.02.25 на основе программы, использующей CUDA. Хотя к моменту создания проекта я уже полностью переписал код,
// решающий задачу методом Ньютона. Теперь же удалил весь код для CUDA
// Необходимо добавить в пограмму проксимальную регуляризацию. Попробовать разные способы для исключения узлов Дирихле.
// Сравнить с результатами другой программы.

// Версия без ругуляризации скопирована в каталоге
// Это попытка сделать версию с регуляризацией 2019.02.26

#include "pch.h"
#include <iostream>


//#include <stdio.h>
#include <math.h>
#include <fstream>
#include <chrono>

#define IDX2C(i,j,ld) (((j)*(ld))+(i)) // Используется для привычной индексации элементов матрицы


static void makeATable(double lambda, double mu, double h, double tA[12][7][4]);

int main()
{
	static double E;		// Модуль Юнга
	static double nu;		// коэффициент Пуассона
	static double f;		// коэффициент трения
	static double h;		// Шаг разбиения
	static double epsilon;	// Точность расчетов
	static long long r;		// Параметр метода МФЛ
	static int n;			// Число узлов вдоль одной из сторон единичного квадрата

	// Ввод параметров
	std::cout << "Input E:    ";	  std::cin >> E;
	std::cout << "Input n:    ";	  std::cin >> n;
	nu = 0.34;
	f = 0.5;

	r = 1e12;	//1e+12
	epsilon = 1e-5;	// 1e-5

	h = 1.0 / (n - 1);
	long long nFE = n * n + n;	// Количество конечных элементов (КЭ)
	long long nFE2 = 2 * nFE;	// Удвоенное количество КЭ, каждому соответствует две компоненты вектора перемещения
	int nw = n;					// Количество узлов в одном теле по горизонтали
	int nh = (n + 1) / 2;		// Количество узлов в одном теле по вертикали

	int iA, jA;		// Индексы для прохода по матрице A
	int ix, iy;		// Индексы для прохода по узлам
	int nb;			// Индекс тела (0 - верхнее и 1 - нижнее)
	int k;			// Индекс направления (0 - горизонтальное и 1 - вертикальное)

	double lambda, mu;			// Параметры Ламе
	double koeffA[12][7][4];	// Массив с коэффициентами
	mu = E / (2 + 2 * nu);
	lambda = nu * E / ((1 + nu) * (1 - 2 * nu));
	makeATable(lambda, mu, h, koeffA);

	// Заполнение матрицы жесткости A
	double* A = new double[nFE2 * nFE2];
	for (iA = 0; iA < nFE2; iA++)
	{
		for (jA = 0; jA < nFE2; jA++)
		{
			A[IDX2C(iA, jA, nFE2)] = 0;
		}
	}
	int nt = -1;	// Тип узла
	int iA0, iA1;	// Номер компоненты вектора перемещений. Номера 1 и 2 компонент v на конкретном конечном элементе
	for (k = 0; k <= 1; k++)
	{// Перебираем компоненты вектора перемещения
		for (nb = 0; nb <= 1; nb++)
		{// Перебираем элементы верхнего тела, затем нижнего. По большей части их матрицы жесткости похожи
			for (iy = 0; iy < nh; iy++)
			{// Проходим по строчкам тела
				for (ix = 0; ix < nw; ix++)
				{// Проходим по точкам тела в строке iy

					// Определяем номер строки в матрице A для конкретной компоненты вектора перемещений
					iA = k * nFE + nb * nw * nh + iy * nw + ix;
					iA0 = 0 * nFE + nb * nw * nh + iy * nw + ix;
					iA1 = 1 * nFE + nb * nw * nh + iy * nw + ix;

					if ((ix > 0) && (ix < nw - 1) && (iy > 0) && (iy < nh - 1)) nt = 0;	// Внутренние узлы тела
					else if ((iy == 0) && (ix == 0)) nt = 1;							// Верхний левый угол
					else if ((iy == 0) && (ix > 0) && (ix < nw - 1)) nt = 2;			// Внутренние узлы верхней границы
					else if ((iy == 0) && (ix == nw - 1)) nt = 3;						// Верхний правый угол
					else if ((iy > 0) && (iy < nh - 1) && (ix == 0)) nt = 4;			// Внутренние узлы левой границы
					else if ((iy > 0) && (iy < nh - 1) && (ix == nw - 1)) nt = 5;		// Внутренние узлы правой границы
					else if ((iy == nh - 1) && (ix == 0)) nt = 6;						// Нижний левый угол
					else if ((iy == nh - 1) && (ix > 0) && (ix < nw - 1)) nt = 7;		// Внутренние узлы нижней границы
					else if ((iy == nh - 1) && (ix == nw - 1)) nt = 8;					// Нижний правый угол
					else nt = -1;

					if (koeffA[nt][0][0 + k * 2] != 0) A[IDX2C(iA, iA0,				nFE2)] = koeffA[nt][0][0 + k * 2];
					if (koeffA[nt][1][0 + k * 2] != 0) A[IDX2C(iA, iA0 + 1,			nFE2)] = koeffA[nt][1][0 + k * 2];
					if (koeffA[nt][2][0 + k * 2] != 0) A[IDX2C(iA, iA0 + 1 - nw,	nFE2)] = koeffA[nt][2][0 + k * 2];
					if (koeffA[nt][3][0 + k * 2] != 0) A[IDX2C(iA, iA0 - nw,		nFE2)] = koeffA[nt][3][0 + k * 2];
					if (koeffA[nt][4][0 + k * 2] != 0) A[IDX2C(iA, iA0 - 1,			nFE2)] = koeffA[nt][4][0 + k * 2];
					if (koeffA[nt][5][0 + k * 2] != 0) A[IDX2C(iA, iA0 - 1 + nw,	nFE2)] = koeffA[nt][5][0 + k * 2];
					if (koeffA[nt][6][0 + k * 2] != 0) A[IDX2C(iA, iA0 + nw,		nFE2)] = koeffA[nt][6][0 + k * 2];

					if (koeffA[nt][0][1 + k * 2] != 0) A[IDX2C(iA, iA1,				nFE2)] = koeffA[nt][0][1 + k * 2];
					if (koeffA[nt][1][1 + k * 2] != 0) A[IDX2C(iA, iA1 + 1,			nFE2)] = koeffA[nt][1][1 + k * 2];
					if (koeffA[nt][2][1 + k * 2] != 0) A[IDX2C(iA, iA1 + 1 - nw,	nFE2)] = koeffA[nt][2][1 + k * 2];
					if (koeffA[nt][3][1 + k * 2] != 0) A[IDX2C(iA, iA1 - nw,		nFE2)] = koeffA[nt][3][1 + k * 2];
					if (koeffA[nt][4][1 + k * 2] != 0) A[IDX2C(iA, iA1 - 1,			nFE2)] = koeffA[nt][4][1 + k * 2];
					if (koeffA[nt][5][1 + k * 2] != 0) A[IDX2C(iA, iA1 - 1 + nw,	nFE2)] = koeffA[nt][5][1 + k * 2];
					if (koeffA[nt][6][1 + k * 2] != 0) A[IDX2C(iA, iA1 + nw,		nFE2)] = koeffA[nt][6][1 + k * 2];

				}
			}
		}
	}

	// Третий способ
	// Отметим в специальном массиве переменные, на которые наложено условие Дирихле
	//bool * indexDirihle = new bool[nFE2];
	//for (int i = 0; i < nFE2; i++)
	//{
	//	indexDirihle[i] = false;
	//}
	//int countDirihle = 0; // Количество переменных, на которые наложено условие Дирихле
	// Конец кусочка третьего способа

	for (iy = 0; iy < nh - 1; iy++)
	{// Правая стенка верхнего тела не должна перемещаться вдоль горизонтальной оси

		// Первый спобоб - Большое число на диагонали
		//A[IDX2C(iy * nw + nw - 1, iy * nw + nw - 1, nFE2)] = 1.0e+100; // Вообще и 1.0e+10 вполне подходило

		// Второй способ - зануление строки и столбца и единица на диагонали
		for (int i = 0; i < nFE2; i++)
		{
			A[IDX2C(0 * nFE + 0 * nw * nh + nw * iy + nw - 1, i, nFE2)] = 0;
			A[IDX2C(i, 0 * nFE + 0 * nw * nh + nw * iy + nw - 1, nFE2)] = 0;

			//A[IDX2C(1 * nFE + 0 * nw * nh + nw * iy + nw - 1, i, nFE2)] = 0; // полностью
			//A[IDX2C(i, 1 * nFE + 0 * nw * nh + nw * iy + nw - 1, nFE2)] = 0; // полностью
		}
		A[IDX2C(0 * nFE + iy * nw + nw - 1, 0 * nFE + iy * nw + nw - 1, nFE2)] = 1;

		//A[IDX2C(1 * nFE + iy * nw + nw - 1, 1 * nFE + iy * nw + nw - 1, nFE2)] = 1; // полностью

		// Третий способ - вычеркивание строк и столбцов из матрицы (для этого сперва запоминаем их индексы)
		//indexDirihle[iy * nw + nw - 1] = true;
		//countDirihle++;
	}
	for (ix = 0; ix < nw; ix++)
	{// Нижняя сторона нижнего тела полностью закреплена

		// Первый спобоб
		//A[IDX2C(0 * nFE + 1 * nw * nh + nw * (nh - 1) + ix, 0 * nFE + 1 * nw * nh + nw * (nh - 1) + ix, nFE2)] = 1.0e+100; // Большое число
		//A[IDX2C(1 * nFE + 1 * nw * nh + nw * (nh - 1) + ix, 1 * nFE + 1 * nw * nh + nw * (nh - 1) + ix, nFE2)] = 1.0e+100; // Большое число

		// Второй способ
		for (int i = 0; i < nFE2; i++)
		{
			A[IDX2C(0 * nFE + 1 * nw * nh + nw * (nh - 1) + ix, i, nFE2)] = 0;
			A[IDX2C(i, 0 * nFE + 1 * nw * nh + nw * (nh - 1) + ix, nFE2)] = 0;
			A[IDX2C(1 * nFE + 1 * nw * nh + nw * (nh - 1) + ix, i, nFE2)] = 0;
			A[IDX2C(i, 1 * nFE + 1 * nw * nh + nw * (nh - 1) + ix, nFE2)] = 0;
		}
		A[IDX2C(0 * nFE + 1 * nw * nh + nw * (nh - 1) + ix, 0 * nFE + 1 * nw * nh + nw * (nh - 1) + ix, nFE2)] = 1;
		A[IDX2C(1 * nFE + 1 * nw * nh + nw * (nh - 1) + ix, 1 * nFE + 1 * nw * nh + nw * (nh - 1) + ix, nFE2)] = 1;

		// Третий способ
		//indexDirihle[0 * nFE + 1 * nw * nh + nw * (nh - 1) + ix] = true;
		//countDirihle++;
		//indexDirihle[1 * nFE + 1 * nw * nh + nw * (nh - 1) + ix] = true;
		//countDirihle++;
	}

	// Модифицируем матрицу A в соответствии с заменой переменных
	for (ix = 0; ix < nw; ix++)
	{// Проходим по конечным элементам нижнего тела, лежащим в зоне контакта. К их строкам в матрице A нужно прибавить соответствующие строки
		iA = 1 * nFE + 1 * nw * nh + ix;
		for (jA = 0; jA < nFE2; jA++)
		{
			A[IDX2C(iA, jA, nFE2)] += A[IDX2C(iA - nw, jA, nFE2)]; // Добавлем к коэффициенту нижнего тела коэффициент верхнего
		}
		for (jA = 0; jA < nFE2; jA++)
		{
			A[IDX2C(jA, iA, nFE2)] += A[IDX2C(jA, iA - nw, nFE2)]; // Сразу же (после прохода по всей строке) проделаем то же самое для столбца
		}
	}
	for (ix = 0; ix < nw; ix++)
	{
		iA = 1 * nFE + 1 * nw * (nh - 1) + ix;
		for (jA = 0; jA < nFE2; jA++)
		{
			A[IDX2C(iA, jA, nFE2)] *= -1;
			A[IDX2C(jA, iA, nFE2)] *= -1; // Меняем знак у коэффициента верхнего тела
		}
	}// Это только замена вертикальной составляющей. Для задачи с трением нужно это же проделать и для горизонтальной

	// Замена для горизонтальной составляющей
	for (ix = 0; ix < nw; ix++)
	{// Проходим по конечным элементам нижнего тела, лежащим в зоне контакта. К их строкам в матрице A нужно прибавить соответствующие строки
		iA = 0 * nFE + 1 * nw * nh + ix;
		for (jA = 0; jA < nFE2; jA++)
		{
			A[IDX2C(iA, jA, nFE2)] += A[IDX2C(iA - nw, jA, nFE2)]; // Добавлем к коэффициенту нижнего тела коэффициент верхнего
		}
		for (jA = 0; jA < nFE2; jA++)
		{
			A[IDX2C(jA, iA, nFE2)] += A[IDX2C(jA, iA - nw, nFE2)]; // Сразу же (после прохода по всей строке) проделаем то же самое для столбца
		}
	}
	for (ix = 0; ix < nw; ix++)
	{
		iA = 0 * nFE + 1 * nw * (nh - 1) + ix;
		for (jA = 0; jA < nFE2; jA++)
		{
			A[IDX2C(iA, jA, nFE2)] *= -1;
			A[IDX2C(jA, iA, nFE2)] *= -1; // Меняем знак у коэффициента верхнего тела
		}
	}

	// Формируем вектор нагрузки FP
	double* FP = new double[nFE2];
	for (iA = 0; iA < nFE2; iA++)
	{
		FP[iA] = 0;
	}
	for (int k = 0; k <= 1; k++)
	{// Направление нагрузки по горизонтали/вертикали
		for (int nb = 0; nb <= 1; nb++)
		{// Перебираем элементы верхнего тела, затем нижнего
			for (int iy = 0; iy < nh; iy++)
			{// Проходим по строчкам тела
				for (int ix = 0; ix < nw; ix++)
				{// Проходим по точкам тела в строке iy
				 // Определяем тип узла для правильного вычисления интеграла
					if ((ix > 0) && (ix < nw - 1) && (iy > 0) && (iy < nh - 1)) nt = 0; // Внутренние узлы тела
					else if ((iy == 0) && (ix == 0)) nt = 1; // Верхний левый угол
					else if ((iy == 0) && (ix > 0) && (ix < nw - 1)) nt = 2; // Внутренние узлы верхней границы
					else if ((iy == 0) && (ix == nw - 1)) nt = 3; // Верхний правый угол
					else if ((iy > 0) && (iy < nh - 1) && (ix == 0)) nt = 4; // Внутренние узлы левой границы
					else if ((iy > 0) && (iy < nh - 1) && (ix == nw - 1)) nt = 5; // Внутренние узлы правой границы
					else if ((iy == nh - 1) && (ix == 0)) nt = 6; // Нижний левый угол
					else if ((iy == nh - 1) && (ix > 0) && (ix < nw - 1)) nt = 7; // Внутренние узлы нижней границы
					else if ((iy == nh - 1) && (ix == nw - 1)) nt = 8; // Нижний правый угол
					else nt = -1;

					double ku, pointFP = 0;

					// Для объемных сил
					if (nt == 0) ku = h * h;
					else if ((nt == 1) || nt == 8) ku = h * h / 6;
					else if ((nt == 3) || nt == 6) ku = h * h / 3;
					else ku = h * h / 2;

					pointFP = 0 * ku; // Здесь определяем с помощью условий, для каких точек задать какую силу P

									  // Для поверхностных сил
					if (nt == 1 || nt == 3 || nt == 6 || nt == 8) ku = h / 2;
					else ku = h;

					// Здесь определяем с помощью условий, для каких точек задать какую силу F
					if ((nb == 0) && (k == 1) && (iy == 0) && (ix < nw / 3)) pointFP += -60 * ku; //-60


					// Записываем полученное значение в вектор FP
					FP[k * nFE + nb * nw * nh + iy * nw + ix] = pointFP;
				}
			}
		}
	}


	// Третий способ
	//int * newIndex = new int[nFE2]; // Массив показывает, какие индексы стали у узлов после удаления узлов Дирихле
	//int count = 0; // Количество узлов, если выкинуть узлы Дирихле
	//double* newA = new double[(nFE2 - countDirihle) * (nFE2 - countDirihle)];
	//for (int i = 0; i < nFE2; i++)
	//{
	//	if (indexDirihle[i] == false)
	//	{
	//		newIndex[i] = count;
	//		FP[count] = FP[i];
	//		for (int j = i, jj = count; j < nFE2; j++)
	//		{
	//			if (indexDirihle[j] == false)
	//			{
	//				newA[IDX2C(count, jj, nFE2 - countDirihle)] = A[IDX2C(i, j, nFE2)];
	//				newA[IDX2C(jj, count, nFE2 - countDirihle)] = A[IDX2C(j, i, nFE2)];
	//				jj++;
	//			}
	//		}
	//		count++;
	//	}
	//	else
	//	{
	//		newIndex[i] = NULL;
	//	}
	//}
	//for (int i = 0; i < (nFE2 - countDirihle) * (nFE2 - countDirihle); i++)
	//{
	//	A[i] = newA[i];
	//}
	//nFE2 = count;



	//std::fstream ofs("output.txt", std::ios_base::out);

	//for (iA = 0; iA < nFE2; iA++)
	//{
	//	for (jA = 0; jA < nFE2; jA++)
	//	{
	//		ofs << A[IDX2C(iA, jA, nFE2)] << "\t";
	//	}
	//	ofs << "\n";
	//}
	//ofs << "\n\n\n\n\n";
	//ofs.close();
	// Конец кусочка с третьим способом



	double delta = 0;
	int iter1 = 0, iter2 = 0, iter3 = 0;
	double* v = new double[nFE2];
	double* v2 = new double[nFE2];
	double* vOld = new double[nFE2];
	double* hg = new double[nFE2]; // шаг умножить на трение в этой точке
	double* hlrt = new double[nFE2]; // шаг умножить на положительную срезку
	for (iA = 0; iA < nFE2; iA++)
	{
		v[iA] = 0;
		v2[iA] = 0;
		hg[iA] = 0;
		hlrt[iA] = 0;
	}
	double* Hessian = new double[nFE2 * nFE2];
	double* l = new double[nw]; // двойственная переменная
	double* t = new double[nw]; // скачок вертикальной составляющей (замена переменных)
	double* s = new double[nw]; // скачок горизонтальной составляющей (замена переменных)
	double* g = new double[nw]; // трение
	for (ix = 0; ix < nw; ix++)
	{
		l[ix] = 0;
		t[ix] = 0;
		s[ix] = 0;
		g[ix] = 0;
	}

	// для регуляризации
	double* vk = new double[nFE2];
	double* tk = new double[nw]; // скачок вертикальной составляющей (замена переменных)
	double* sk = new double[nw]; // скачок горизонтальной составляющей (замена переменных)
	for (iA = 0; iA < nFE2; iA++)
	{
		vk[iA] = 0;
	}
	for (ix = 0; ix < nw; ix++)
	{
		tk[ix] = 0;
		sk[ix] = 0;
	}
	double ku; // Какой-то коэффициент для регуляризации



	// Засекаем время начала работы алгоритма
	auto start_time = std::chrono::steady_clock::now();
	auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time);
	int min, sec, msec;


	do
	{// Метод последовательных приближений
		do
		{// Метод Удзавы
			do
			{// Метод поточечной релаксации
				for (iA = 0; iA < nFE2; iA++)
				{
					vOld[iA] = v[iA];
				}

				
				for (int iy = 0; iy < nh; iy++)
				{// Проходим по строчкам тела
					for (int ix = 0; ix < nw; ix++)
					{

						if ((iy > 0) && (iy < nh - 1) && (ix > 0) & (ix < nw - 1)) nt = 0;
						else if ((iy == 0) && (ix == 0)) nt = 1;
						else if ((iy == 0) && (ix > 0) && (ix < nw - 1)) nt = 2;
						else if ((iy == 0) && (ix == nw-1)) nt = 3;
						else if ((iy > 0) && (iy < nh-1) && (ix == 0)) nt = 4;
						else if ((iy > 0) && (iy < nh-1) && (ix == nw - 1)) nt = 5;
						else if ((iy == nh - 1) && (ix == 0)) nt = 6;
						else if ((iy == nh - 1) && (ix > 0) && (ix < nw - 1)) nt = 7;
						else if ((iy == nh - 1) && (ix == nw - 1)) nt = 8;
						else nt = -1;

						// Для объемных сил
						if (nt == 0) ku = h * h;
						else if ((nt == 1) || nt == 8) ku = h * h / 6;
						else if ((nt == 3) || nt == 6) ku = h * h / 3;
						else ku = h * h / 2;
						
						vk[0 * nFE + 0 * nw * nh + iy * nw + ix] = ku * v[0 * nFE + 0 * nw * nh + iy * nw + ix];
						vk[0 * nFE + 1 * nw * nh + iy * nw + ix] = ku * v[0 * nFE + 1 * nw * nh + iy * nw + ix];
						vk[1 * nFE + 0 * nw * nh + iy * nw + ix] = ku * v[1 * nFE + 0 * nw * nh + iy * nw + ix];
						vk[1 * nFE + 1 * nw * nh + iy * nw + ix] = ku * v[1 * nFE + 1 * nw * nh + iy * nw + ix];
					}
					if ((ix > 0) && (ix < nw - 1)) ku = h;
					else ku = h / 2;
					tk[ix] = ku * t[ix];
					sk[ix] = ku * s[ix];
				}

				for (int k = 0; k <= 1; k++)
				{// По горизонтали/вертикали
					for (int nb = 0; nb <= 1; nb++)
					{// Перебираем элементы верхнего тела, затем нижнего
						for (int iy = 0; iy < nh; iy++)
						{// Проходим по строчкам тела
							for (int ix = 0; ix < nw; ix++)
							{
								iA = k * nFE + nb * nw * nh + iy * nw + ix;
								iA0 = 0 * nFE + nb * nw * nh + iy * nw + ix;
								iA1 = 1 * nFE + nb * nw * nh + iy * nw + ix;

								if (nb == 1 && iy == nh - 1)
								{
									continue; // нижняя часть нижнего тела закреплена, так что пропускаем
								}
								else if ((nb == 0) && (iy == nh - 1))
								{
									if (k == 0)
									{ // s - горизональный скачок
										double psi = - (FP[iA] + vk[iA]);
										for (jA = 0; jA < nFE2; jA++)
										{
											psi += A[IDX2C(iA, jA, nFE2)] * v[jA];
										}
										psi -= A[IDX2C(iA, iA, nFE2)] * v[iA]; // Вычитаем лишнее слагаемое (думаю, так быстрей, чем выполняться проверку в цикле)

										if (psi + h * g[ix] < 0)
											v[iA] = (-psi - h * g[ix]) / A[IDX2C(iA, iA, nFE2)];
										else if (psi - h * g[ix] > 0)
											v[iA] = (-psi + h * g[ix]) / A[IDX2C(iA, iA, nFE2)];
										else
											v[iA] = 0;
									}
									else if (k == 1)
									{ // t - вертикальный скачок
										double psi = FP[iA] + vk[iA];
										for (jA = 0; jA < nFE2; jA++)
										{
											psi -= A[IDX2C(iA, jA, nFE2)] * v[jA];
										}
										psi += A[IDX2C(iA, iA, nFE2)] * v[iA]; // Прибавляем лишнее вычитаемое (думаю, так быстрей, чем выполняться проверку в цикле)

										psi /= A[IDX2C(iA, iA, nFE2)];

										if (l[ix] + r * psi <= 0)
											v[iA] = psi;
										else
										{
											psi = FP[iA] + vk[iA];
											for (jA = 0; jA < nFE2; jA++)
											{
												if (jA == iA) continue;

												psi -= A[IDX2C(iA, jA, nFE2)] * v[jA];
											}
											psi -= h * l[ix];
											psi /= (A[IDX2C(iA, iA, nFE2)] + h * r);
											v[iA] = psi;
										}
									}
								}
								else
								{
									v[iA] = FP[iA] + vk[iA];

									if ((nb == 0 && iy == nh - 1) || (nb == 0 && iy == nh - 2) || (nb == 1 && iy == 0))
									{
										for (jA = 0; jA < nFE2; jA++)
										{
											if (jA == iA) continue;

											v[iA] -= A[IDX2C(iA, jA, nFE2)] * v[jA];
										}
									}
									else
									{

										if (k == 0) v[iA] -= A[IDX2C(iA, iA1, nFE2)] * v[iA1];
										else v[iA] -= A[IDX2C(iA, iA0, nFE2)] * v[iA0];

										if (ix != (nw - 1))
										{
											v[iA] -= A[IDX2C(iA, iA0 + 1, nFE2)] * v[iA0 + 1];
											v[iA] -= A[IDX2C(iA, iA1 + 1, nFE2)] * v[iA1 + 1];
										}

										if ((ix != (nw - 1)) && (iy != 0))
										{
											v[iA] -= A[IDX2C(iA, iA0 + 1 - nw, nFE2)] * v[iA0 + 1 - nw];
											v[iA] -= A[IDX2C(iA, iA1 + 1 - nw, nFE2)] * v[iA1 + 1 - nw];
										}

										if (iy != 0)
										{
											v[iA] -= A[IDX2C(iA, iA0 - nw, nFE2)] * v[iA0 - nw];
											v[iA] -= A[IDX2C(iA, iA1 - nw, nFE2)] * v[iA1 - nw];
										}

										if (ix != 0)
										{
											v[iA] -= A[IDX2C(iA, iA0 - 1, nFE2)] * v[iA0 - 1];
											v[iA] -= A[IDX2C(iA, iA1 - 1, nFE2)] * v[iA1 - 1];
										}

										if ((ix != 0) && (iy != (nh - 1)))
										{
											v[iA] -= A[IDX2C(iA, iA0 - 1 + nw, nFE2)] * v[iA0 - 1 + nw];
											v[iA] -= A[IDX2C(iA, iA1 - 1 + nw, nFE2)] * v[iA1 - 1 + nw];
										}

										if ((iy != (nh - 1)))
										{
											v[iA] -= A[IDX2C(iA, iA0 + nw, nFE2)] * v[iA0 + nw];
											v[iA] -= A[IDX2C(iA, iA1 + nw, nFE2)] * v[iA1 + nw];
										}
									}

									v[iA] /= A[IDX2C(iA, iA, nFE2)];
								}


								// если третий способ
								// v2[indexDirihle[index]] = hlrt[indexDirihle[index]] - FP[indexDirihle[index]];
							}
						}
					}
				}


				for (int ix = 0; ix < nw; ix++)
				{// Дублируем элементы в переменную t
					t[ix] = v[1 * nFE + (nh - 1) * nw + ix];

					// если третий способ
					// t[ix] = v[indexDirihle[1 * nFE + (nh - 1) * nw + ix]];
				}

				delta = 0;
				for (int k = 0; k <= 1; k++)
				{// Проход по компонентам
					for (int nb = 0; nb <= 1; nb++)
					{// Перебираем элементы верхнего тела, затем нижнего
						for (int iy = 0; iy < nh; iy++)
						{// Проходим по строчкам тела
							for (int ix = 0; ix < nw; ix++)
							{// Проходим по точкам тела в строке iy
								int index = k * nFE + nb * nw * nh + iy * nw + ix;
								delta = fmax(delta, abs(v[index] - vOld[index]));

								// если третий способ
								// delta = fmax(delta, abs(v[indexDirihle[index]] - vOld[indexDirihle[index]]));
							}
						}
					}
				}
				iter1++;
				elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time);
				msec = elapsed_ms.count();
				min = msec / 60000;
				sec = (msec - min * 60000) / 1000;
				msec = msec % 1000;

				if (iter1 % 1000 == 0) printf("      %i:%.2i:%.3i      %.2i    %.10f    \n", min, sec, msec, iter1, delta);



				/*if (iter3 > 0)
				{
					std::cout << std::endl;
					for (int ix = 0; ix < nw; ix++)
					{
						printf("%.3i    %.10f     %.10f\n", ix, l[ix], t[ix]);
					}
					std::cout << std::endl;
					system("pause");
				}*/



			} while (delta > epsilon * h * 1e-4);
			printf("      %i:%.2i:%.3i      %.2i    %.10f    \n", min, sec, msec, iter1, delta);


			delta = 0;
			double lPointOld;
			for (int ix = 0; ix < nw; ix++)
			{
				lPointOld = l[ix];

				l[ix] = fmax(0, l[ix] + r * t[ix]);
				delta = fmax(delta, abs(l[ix] - lPointOld));
			}
			iter2++;
			printf("      -------------------------------------------------------\n");
			printf("   %i                        %.10f\n", iter2, delta);
			printf("   ==========================================================\n");
		} while (delta >= epsilon);

		delta = 0;
		double gPointOld;
		for (int ix = 0; ix < nw; ix++)
		{
			gPointOld = g[ix];

			g[ix] = f * abs(l[ix]);

			delta = fmax(delta, abs(g[ix] - gPointOld));
		}
		iter3++;
		printf("%i                             %.10f\n", iter3, delta);
		printf("#############################################################\n");

	} while (delta >= epsilon);

	std::cout << std::endl;
	for (int ix = 0; ix < nw; ix++)
	{
		printf("%.3i    %.10f     %.10f\n", ix, l[ix], t[ix]);
	}
	std::cout << std::endl;






	// Сохранение полученного решения в файл
	std::ofstream out;          // поток для записи
	out.open("output_ut.txt"); // окрываем файл для записи
	if (out.is_open())
	{
		for (int ix = 0; ix < nw; ix++)
		{
			out << "(" << ix * h + v[0 * nFE + 1 * nw * nh + ix] - v[0 * nFE + (nh - 1) * nw + ix] << "," << v[1 * nFE + 1 * nw * nh + ix] - v[1 * nFE + (nh - 1) * nw + ix] << ")" << std::endl;
		}

		for (int ix = 0; ix < nw; ix++)
		{
			out << "(" << ix * h + v[0 * nFE + 1 * nw * nh + ix] << "," << v[1 * nFE + 1 * nw * nh + ix] << ")" << std::endl;
		}
	}


	return 0;
}


static void makeATable(double lambda, double mu, double h, double tA[12][7][4])
{
	double tI[4][6][7] = {
	{
		{ +0.5, -0.5, +0.0, +0.0, +0.0, +0.0, +0.0 },
		{ +0.0, +0.0, +0.0, +0.0, +0.0, +0.0, +0.0 },
		{ +0.5, +0.0, +0.0, +0.0, -0.5, +0.0, +0.0 },
		{ +0.5, +0.0, +0.0, +0.0, -0.5, +0.0, +0.0 },
		{ +0.0, +0.0, +0.0, +0.0, +0.0, +0.0, +0.0 },
		{ +0.5, -0.5, +0.0, +0.0, +0.0, +0.0, +0.0 }
	},
	{
		{ +0.0, +0.0, +0.0, +0.0, +0.0, +0.0, +0.0 },
		{ +0.5, +0.0, +0.0, -0.5, +0.0, +0.0, +0.0 },
		{ +0.5, +0.0, +0.0, -0.5, +0.0, +0.0, +0.0 },
		{ +0.0, +0.0, +0.0, +0.0, +0.0, +0.0, +0.0 },
		{ +0.5, +0.0, +0.0, +0.0, +0.0, +0.0, -0.5 },
		{ +0.5, +0.0, +0.0, +0.0, +0.0, +0.0, -0.5 }
	},
	{
		{ +0.0, +0.5, -0.5, +0.0, +0.0, +0.0, +0.0 },
		{ +0.0, +0.0, +0.0, +0.0, +0.0, +0.0, +0.0 },
		{ -0.5, +0.0, +0.0, +0.5, +0.0, +0.0, +0.0 },
		{ +0.0, +0.0, +0.0, +0.0, +0.5, -0.5, +0.0 },
		{ +0.0, +0.0, +0.0, +0.0, +0.0, +0.0, +0.0 },
		{ -0.5, +0.0, +0.0, +0.0, +0.0, +0.0, +0.5 }
	},
	{
		{ +0.0, +0.0, +0.0, +0.0, +0.0, +0.0, +0.0 },
		{ +0.0, +0.0, -0.5, +0.5, +0.0, +0.0, +0.0 },
		{ -0.5, +0.0, +0.0, +0.0, +0.5, +0.0, +0.0 },
		{ +0.0, +0.0, +0.0, +0.0, +0.0, +0.0, +0.0 },
		{ +0.0, +0.0, +0.0, +0.0, +0.0, -0.5, +0.5 },
		{ -0.5, +0.5, +0.0, +0.0, +0.0, +0.0, +0.0 }
	} };

	double tC[9][7][4];
	// Первый индекс - тип узла. Зависит от его расположения относительно окружающих КЭ.
	// 0 - большинство узлов, со всех сторон окружены другими КЭ
	// 1 - верхний левый узел, рядом с ним только три КЭ, причем общий носитель - только один треугольный элемент
	// 2 - верхние не угловые узлы, у них есть пять соседей, пересечение с ними может быть в одном из трех треугольников
	// ...

	for (int i = 0; i <= 6; i++)
	{
		for (int j = 0; j <= 3; j++)
		{
			tC[0][i][j] = tI[j][0][i] + tI[j][1][i] + tI[j][2][i] + tI[j][3][i] + tI[j][4][i] + tI[j][5][i];
			tC[1][i][j] = tI[j][5][i];
			tC[2][i][j] = tI[j][3][i] + tI[j][4][i] + tI[j][5][i];
			tC[3][i][j] = tI[j][3][i] + tI[j][4][i];
			tC[4][i][j] = tI[j][0][i] + tI[j][1][i] + tI[j][5][i];
			tC[5][i][j] = tI[j][2][i] + tI[j][3][i] + tI[j][4][i];
			tC[6][i][j] = tI[j][0][i] + tI[j][1][i];
			tC[7][i][j] = tI[j][0][i] + tI[j][1][i] + tI[j][2][i];
			tC[8][i][j] = tI[j][2][i];
		}
	}

	double unciahh = h * h / 12;
	double tU[6][7]= {
		{unciahh, unciahh / 2, unciahh / 2,           0,           0,           0,           0},
		{unciahh,           0, unciahh / 2, unciahh / 2,           0,           0,           0},
		{unciahh,           0,           0, unciahh / 2, unciahh / 2,           0,           0},
		{unciahh,           0,           0,           0, unciahh / 2, unciahh / 2,           0},
		{unciahh,           0,           0,           0,           0, unciahh / 2, unciahh / 2},
		{unciahh, unciahh / 2,           0,           0,           0,           0, unciahh / 2},
	};
	double tAu[9][7];
	for (int i = 0; i <= 6; i++)
	{
		tAu[0][i] = tU[0][i] + tU[1][i] + tU[2][i] + tU[3][i] + tU[4][i] + tU[5][i];
		tAu[1][i] = tU[5][i];
		tAu[2][i] = tU[3][i] + tU[4][i] + tU[5][i];
		tAu[3][i] = tU[3][i] + tU[4][i];
		tAu[4][i] = tU[0][i] + tU[1][i] + tU[5][i];
		tAu[5][i] = tU[2][i] + tU[3][i] + tU[4][i];
		tAu[6][i] = tU[0][i] + tU[1][i];
		tAu[7][i] = tU[0][i] + tU[1][i] + tU[2][i];
		tAu[8][i] = tU[2][i];
	}

	for (int i = 0; i <= 8; i++)
	{
		for (int j = 0; j <= 6; j++)
		{
			// Коэффициенты при x1 и x2 при рассчете компоненты x1
			tA[i][j][0] = (2 * mu + lambda) * tC[i][j][0] + mu * tC[i][j][1] / 2 +tAu[i][j];
			tA[i][j][1] = lambda * tC[i][j][2] + mu * tC[i][j][3] / 2;
			// Коэффициенты при x1 и x2 при рассчете компоненты x2
			tA[i][j][2] = mu * tC[i][j][2] / 2 + lambda * tC[i][j][3];
			tA[i][j][3] = mu * tC[i][j][0] / 2 + (2 * mu + lambda) * tC[i][j][1] +tAu[i][j];
		}
	}
}