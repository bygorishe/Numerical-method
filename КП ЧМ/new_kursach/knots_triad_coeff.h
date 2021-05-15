#pragma once
#include "structs_nums_operations.h"

knot* read_knots() {
    ifstream input_knots("knots.txt");
    if (!input_knots.is_open()) {
        input_knots.close();
    }
    input_knots >> num_knots;
    knot* knots = new knot[num_knots];
    for (int i = 0; i < num_knots; i++)
        input_knots >> knots[i].x >> knots[i].y;
    return knots;
}

void create_lambda_f_gamma(vector<int>& lambda, vector<int>& f, vector<int>& gamma) {
    double temp;
    ifstream lambda_file("lambda.txt");
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

void create_triangles(vector<triangle*>& triangle_list) {
    ifstream f("triangles.txt");
    f >> triangle_num;
    triangle_list.resize(triangle_num);
    int number;
    for (int i = 0; i < triangle_num; i++) {
        triangle* triad = new triangle;
        for (int j = 0; j < 3; j++) {
            f >> number;
            triad->knot_nums[j] = number - 1;
        }
        f >> number;
        triad->knot_nums[3] = number - 1; //подобласть
        triangle_list[i] = triad;
    }
}

double func_f(int number_f, int i, knot*& knots) {
    switch (number_f) {
    case 1: {
        return knots[i].x * knots[i].x;
        //return knots[i].y * knots[i].y;
        break;
    }
    case 2: {
        return (knots[i].x * knots[i].y);
        //return (knots[i].x + knots[i].y);
        break;
    }
    case 3: {
        return (knots[i].x * knots[i].y) * (knots[i].x * knots[i].y) + 10.;
        break;
    }
    case 4: {
        return pow(knots[i].x, 3) - 6 * knots[i].x;
        break;
    }
    case 5: {
        return 3 * cos(knots[i].x + knots[i].y);
        break;
    }
    case 6: {
        return 0.;
        break;
    }
    case 7: {
        return 1.;
        break;
    }
    case 8: {
        return knots[i].x * knots[i].x - 2.;
        break;
    }
    case 9: {
        return -2.;
        break;
    }
    case 10: {
        return knots[i].x;
        break;
    }
    case 11: {
        return pow(knots[i].x, 4) - 12 * knots[i].x * knots[i].x;
        break;
    }
    }
}

double func_lambda(int number_f, knot*& knots) {
    switch (number_f) {
    case 1: {
        return 0.;
        break;
    }
    case 2: {
        return 1.;
        break;
    }
    case 3: {
        return 5.;
        break;
    }
    }
}

double func_gamma(int number_f, knot*& knots) {
    switch (number_f) {
    case 1: {
        return 0.;
        break;
    }
    case 2: {
        return 5.;
        break;
    }
    case 3: {
        return 1.;
        break;
    }
    }
}