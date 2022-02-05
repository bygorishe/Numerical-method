#pragma once
#include "structs_nums_operations.h"

void read_knots(vector<knot>& knots) {
    ifstream input_knots("knots.txt");
    input_knots >> num_knots;
    knots.resize(num_knots);
    for (int i = 0; i < num_knots; i++)
        input_knots >> knots[i].x >> knots[i].y;
}

void create_lambda_f_gamma(vector<int>& lambda, vector<int>& f, vector<int>& gamma) {
    int temp;
    ifstream lambda_file("coef.txt");
    lambda_file >> num_lambda;
    lambda.resize(num_lambda);
    for (int i = 0; i < num_lambda; i++) {
        lambda_file >> temp;
        lambda[i] = temp;
        // cout << lambda[i];
    }
    f.resize(num_lambda);
    for (int i = 0; i < num_lambda; i++) {
        lambda_file >> temp;
        f[i] = temp;
        // cout << f[i];
    }
    gamma.resize(num_lambda);
    for (int i = 0; i < num_lambda; i++) {
        lambda_file >> temp;
        gamma[i] = temp;
        // cout << gamma[i];
    }
}

void create_time_vector(vector<double>& time) {
    double temp;
    ifstream time_file("time.txt");
    for (int i = 0; i < 4; i++) {
        time_file >> temp;
        time[i] = temp;
        // cout << time[i];
    }    
}

void create_triangles(vector<triangle>& triangle_list) {
    ifstream f("triangles.txt");
    f >> triangle_num;
    triangle_list.resize(triangle_num);
    int number;
    for (int i = 0; i < triangle_num; i++) {
        triangle triad;
        for (int j = 0; j < 3; j++) {
            f >> number;
            triad.knot_nums[j] = number - 1;
        }
        f >> number;
        triad.knot_nums[3] = number - 1; //подобласть
        triangle_list[i] = triad;
    }
}

double func_lambda(int number_f) {
    switch (number_f) {
    case 0: {
        return 0;
        break;
    }
    case 1: {
        return 1;
        break;
    }
    }
}

double func_f(int number_f, knot knots, double time) {
    switch (number_f) {
    case 0: {
        return 0;
        break;
    }
    case 1: {
        return 4 * time * time * time * knots.x;
        break;
    }
    case 2: {
        return (knots.x + knots.y);
        break;
    }
    }
}

double func_gamma(int number_f) {
    switch (number_f) {
    case 0: {
        return 0;
        break;
    }
    case 1: {
        return 1;
        break;
    }
    case 2: {
        return 2;
        break;
    }
    }
}