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
    double temp;
    ifstream lambda_file("lambda.txt");
    lambda_file >> num_lambda;
    lambda.resize(num_lambda);
    for (int i = 0; i < num_lambda; i++) {
        lambda_file >> temp;
        lambda[i] = temp;
    }
    f.resize(num_lambda);
    for (int i = 0; i < num_lambda; i++) {
        lambda_file >> temp;
        f[i] = temp;
    }
    gamma.resize(num_lambda);
    for (int i = 0; i < num_lambda; i++) {
        lambda_file >> temp;
        gamma[i] = temp;
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

double func_f(int number_f, knot& knots) {
    switch (number_f) {
    case 1: {
        return knots.x * knots.x;
        //return knots[i].y * knots[i].y;
        break;
    }
    case 2: {
        return (knots.x * knots.y);
        //return (knots[i].x + knots[i].y);
        break;
    }
    case 3: {
        return (knots.x * knots.y) * (knots.x * knots.y) + 10.;
        break;
    }
    case 4: {
        return pow(knots.x, 3) - 6 * knots.x;
        break;
    }
    case 5: {
        return 3 * cos(knots.x + knots.y);
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
        return knots.x * knots.x - 2.;
        break;
    }
    case 9: {
        return -2.;
        break;
    }
    case 10: {
        return knots.x;
        break;
    }
    case 11: {
        return pow(knots.x, 4) - 12 * knots.x * knots.x;
        break;
    }
    }
}

double func_lambda(int number_f) {
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

double func_gamma(int number_f) {
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