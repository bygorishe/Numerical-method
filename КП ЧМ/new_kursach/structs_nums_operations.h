#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <windows.h>
#include <iostream>
#include <set>
#include <vector>
#include <fstream> 
#include <algorithm>

//typedef double real;

using namespace std;

struct knot {
    double x = 0, y = 0;
};

struct triangle {
    int knot_nums[4]{};  //4 - номер подобласти
};

struct local {
    // сортируем и потом юзаем в построение глобальной
    vector<int> knot_nums;
    vector<vector<double>> A;
    vector<double> b;
};

struct SLAE {
    vector<int> jg, ig;
    vector<double> ggl, ggu, b, di;
};

int num_knots, num_lambda, triangle_num, num_bounds_1, num_bounds_2, num_bounds_3, n;

vector<double> operator + (vector<double> vector_1, const vector<double>& vector_2) {
    size_t size = vector_1.size();
    for (size_t i = 0; i < size; ++i)
        vector_1[i] -= vector_2[i];
    return vector_1;
}

vector<double>& operator += (vector<double>& vector_1, const vector<double>& vector_2) {
    size_t size = vector_1.size();
    for (size_t i = 0; i < size; ++i)
        vector_1[i] += vector_2[i];
    return vector_1;
}

vector<double> operator * (const double& w, vector<double> vector) {
    size_t size = vector.size();
    for (size_t i = 0; i < size; ++i)
        vector[i] *= w;
    return vector;
}

double scalar(vector<double> x1, vector<double> x2) {
    double sum = 0.;
    for (int i = 0; i < num_knots; i++)
        sum += x1[i] * x2[i];
    return sum;
}