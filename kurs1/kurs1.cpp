// kurs.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <clocale>
#include <cmath>
#include <cstdlib>
#include <stdlib.h>
#include <string>
#include <vector>

// Function to derive/integrate
double function(double x)
{
    return std::sin(x) * x * x;
}

// Approximating derivative using Gaussian method
double gaussian_derive(double (*func)(double), double x0, double step)
{
    return (func(x0 + step) - func(x0 - step)) / (2 * step);
}

// Calculate Legendre polynomial
double legendre(int n, double x)
{
    if (n == 0)
    {
        return 1.0;
    }
    else if (n == 1)
    {
        return x;
    }
    else
    {
        return ((2.0 * n - 1.0) * x * legendre(n - 1, x) - (n - 1) * legendre(n - 2, x)) / n;
    }
}

// Calculate derivative of Legendre polynomial
double dlegendre(int n, double x)
{
    if (n == 0)
    {
        return 0.0;
    }
    else if (n == 1)
    {
        return 1.0;
    }
    else
    {
        return (n / (x * x - 1.0)) * (x * legendre(n, x) - legendre(n - 1, x));
    }
}

// Calculate integral using Gauss-Legendre method
double gauss_legendre_integrate(double (*f)(double), double a, double b, int n)
{
    double* x = new double[n];
    double* w = new double[n];

    int    m = (n + 1) / 2;
    double xm = 0.5 * (b + a);
    double xl = 0.5 * (b - a);

    for (int i = 1; i <= m; i++)
    {
        double z = cos(3.141592 * (i - 0.25) / (n + 0.5));
        double p0 = 1.0;
        double p1 = z;

        for (int j = 2; j <= n; j++)
        {
            double pj = ((2.0 * j - 1.0) * z * p1 - (j - 1.0) * p0) / j;
            p0 = p1;
            p1 = pj;
        }

        x[i - 1] = xm - xl * z;
        x[n - i] = xm + xl * z;
        w[i - 1] = 2.0 * xl / ((1.0 - z * z) * pow(dlegendre(n, z), 2));
        w[n - i] = w[i - 1];
    }

    double integral = 0.0;
    for (int i = 0; i < n; i++)
    {
        integral += w[i] * (*f)(x[i]);
    }

    delete[] x;
    delete[] w;

    return integral;
}

int main(void)
{
    std::setlocale(LC_ALL, "");

    std::string function_str = "x^2 * sin(x)";

    double step, x;

    std::cout << "Введите x для вычисление производной от x: ";
    std::cin >> x;

    std::cout << "Введите размер шага для вычисления производной (меньше - точнее): ";
    std::cin >> step;

    std::cout << std::endl
        << "Функция - " << function_str << "; Производная x: " << gaussian_derive(function, x, step) << std::endl
        << std::endl;

    double a, b;
    int    n;

    std::cout << "Введите нижнее значение производной (a): ";
    std::cin >> a;

    std::cout << "Введите верхнее значение производной (b): ";
    std::cin >> b;

    std::cout << "Введите количество итераций для вычисления интеграла (больше - точнее): ";
    std::cin >> n;

    std::cout << std::endl
        << "Функция - " << function_str << "; Интеграл от a до b: " << gauss_legendre_integrate(function, a, b, n)
        << std::endl;

    return EXIT_SUCCESS;
}
// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
