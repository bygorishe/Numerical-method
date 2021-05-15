#pragma once
#include "knots_triad_coeff.h"

void read_bounds_2_3(vector<local*>& vector_bounds, ifstream& input, int& num_bounds) {
    input >> num_bounds;
    vector_bounds.resize(num_bounds);
    int num;
    for (int i = 0; i < num_bounds; i++) {
        local* local_bound = new local;
        local_bound->knot_nums.resize(4);
        for (int j = 0; j < 2; j++) {//первые два - номера узлов ребра, 3,4 - значение
            input >> num;                     //дл€ 3 кр 3 и 4 - номер функции ub, дл€ 2 кр - тетта
            local_bound->knot_nums[j] = num - 1;
        }
        input >> local_bound->knot_nums[2];
        input >> local_bound->knot_nums[3];
        vector_bounds[i] = local_bound;
    }
}

void read_bounds_1(vector<local*>& vector_bounds, ifstream& input, int& num_bounds) {
    input >> num_bounds;
    vector_bounds.resize(num_bounds);
    int num;
    for (int i = 0; i < num_bounds; i++) {
        local* local_bound = new local;
        local_bound->knot_nums.resize(3);
        for (int j = 0; j < 2; j++) {//первые два - номера узлов ребра
            input >> num;
            local_bound->knot_nums[j] = num - 1;
        }
        input >> num;
        local_bound->knot_nums[2] = num; //номер уравнени€
        local_bound->b.resize(2);
        vector_bounds[i] = local_bound;
    }
}

void build_bound(vector<local*>& vector_bound, int flag, knot*& knots, double betta) {
    int iA = 0;
    for (vector <local*>::iterator iter = vector_bound.begin(); iter != vector_bound.end(); iter++, iA++) {
        local* bounds = *iter;

        bounds->b.resize(2);

        double h = sqrt((knots[bounds->knot_nums[1]].x - knots[bounds->knot_nums[0]].x) * (knots[bounds->knot_nums[1]].x - knots[bounds->knot_nums[0]].x)
            + (knots[bounds->knot_nums[1]].y - knots[bounds->knot_nums[0]].y) * (knots[bounds->knot_nums[1]].y - knots[bounds->knot_nums[0]].y));

        if (flag == 3) { //только дл€ 3 краевых
            bounds->A.resize(2);
            for (int i = 0; i < 2; i++)
                bounds->A[i].resize(2);

            //bounds->knot_nums[2] - B и тд
            bounds->A[0][0] = bounds->A[1][1] = betta * h / 3;
            bounds->A[1][0] = bounds->A[0][1] = betta * h / 6;


            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++)
                    cout << bounds->A[i][j] << " ";
                cout << endl;
            }

            bounds->b[0] = betta * h * (2 * bounds->knot_nums[2] + bounds->knot_nums[3]) / 6;
            bounds->b[1] = betta * h * (bounds->knot_nums[2] + 2 * bounds->knot_nums[3]) / 6;

        }
        else {//вторые
            bounds->b[0] = h * (2 * bounds->knot_nums[2] + bounds->knot_nums[3]) / 6;
            bounds->b[1] = h * (bounds->knot_nums[2] + 2 * bounds->knot_nums[3]) / 6;
        }
        vector_bound[iA] = bounds;
    }
}

double ub(knot*& knots, int i, int num) {
    switch (num) {
    case 12: {
        return knots[i].x;
        break;
    }
    case 13: {
        //return (knots[i].x * knots[i].y);
        return (knots[i].x * knots[i].x);
        break;
    }
    case 3: {
        return cos(knots[i].x + knots[i].y);
        break;
    }
    case 14: {
        return pow(knots[i].x, 3);
        break;
    }
    case 1: {
        return pow(knots[i].x, 4);
        break;
    }
    case 5: {
        return 1;
        break;
    }
    default: {
        return 1;
        break;
    }
    }
}

void use_bounds(vector<local*>& vector_bound, vector<set<int>>& L, SLAE& slae, int bound_num, knot*& knots) {
    if (bound_num == 3) { // 3 краевые услови€
        for (vector <local*>::iterator iter = vector_bound.begin(); iter != vector_bound.end(); iter++) {
            local* bound_iter = *iter;
            //local* A_iter = Local_A[0];

            //заносим высе диагональные элементы
            for (int k = 0; k < 2; k++) {
                slae.di[bound_iter->knot_nums[k]] += bound_iter->A[k][k];
                //cout << "di" << bound_iter->knot_nums[k] << "=" << slae.di[bound_iter->knot_nums[k]] << endl;
                slae.b[bound_iter->knot_nums[k]] += ub(knots, bound_iter->knot_nums[k], 1); // последнее - номер фйункции
                //cout << "b" << bound_iter->knot_nums[k] << "=" << slae.b[bound_iter->knot_nums[k]] << endl;
            }
            //начинаем цикл по строкам нижнего
            for (int i = 0; i < 2; i++) {
                //устанавливаем начальное значение нижней границы поиска
                int ibeg = slae.ig[bound_iter->knot_nums[i]];
                //cout << "knot_nums: " << bound_iter->knot_nums[i] << endl;
                //cout << "ibeg=" << ibeg << endl;

                for (int j = 0; j < i; j++) { // do j=1,i-1
                    int iend = slae.ig[bound_iter->knot_nums[i] + 1] - 1;
                    //cout << "iend=" << iend << endl;

                    //ploblems
                    //while на  for
                    while (slae.jg[ibeg] != bound_iter->knot_nums[j]) {
                        int ind = (ibeg + iend) / 2;
                        if (slae.jg[ind] < bound_iter->knot_nums[j])
                            ibeg = ind + 1;
                        else
                            iend = ind;
                    }
                    slae.ggu[ibeg] += bound_iter->A[j][i];
                    slae.ggl[ibeg] += bound_iter->A[i][j];
                    ibeg++;
                }
            }
        }
    }
    if (bound_num == 2) { // 2 краевые услови€
        for (vector <local*>::iterator iter = vector_bound.begin(); iter != vector_bound.end(); iter++) {
            local* bound_iter = *iter;
            //local* A_iter = Local_A[0];

            //заносим высе диагональные элементы
            for (int k = 0; k < 2; k++) {
                slae.b[bound_iter->knot_nums[k]] += bound_iter->b[k];
                //cout << "b" << bound_iter->knot_nums[k] << "=" << slae.b[bound_iter->knot_nums[k]] << endl;
            }
        }
    }
    /*if (bound_num == 1) { // 1 краевые услови€
        for (vector <local*>::iterator iter = vector_bound.begin(); iter != vector_bound.end(); iter++) {
            local* bound_iter = *iter;
            //local* A_iter = Local_A[0];

            //заносим высе диагональные элементы
            for (int k = 0; k < 2; k++) {
                slae.b[bound_iter->knot_nums[k]] += bound_iter->b[k];
                //cout << "b" << bound_iter->knot_nums[k] << "=" << slae.b[bound_iter->knot_nums[k]] << endl;
            }
        }
    }*/
}

void use_first_bounds(vector<local*>& vector_bound, vector<set<int>> L, SLAE& slae, knot*& knots) {
    for (vector <local*>::iterator iter = vector_bound.begin(); iter != vector_bound.end(); iter++) {
        local* bound = *iter;
        for (int i = 0; i < 2; i++) {
            slae.di[bound->knot_nums[i]] = 1;
            slae.b[bound->knot_nums[i]] = ub(knots, bound->knot_nums[i], bound->knot_nums[2]); //ub(knot*& knots, int i, int num) 
        }
        cout << endl;

        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < n; j++) {
                // cout << "jg" << slae.jg[j] << "  nums" << bound->knot_nums[i] << " ";
                if (slae.jg[j] == bound->knot_nums[i])
                    slae.ggu[j] = 0;
            }
            cout << endl;

            //устанавливаем начальное значение нижней границы поиска
            int ibeg = slae.ig[bound->knot_nums[i]];
            int iend = slae.ig[bound->knot_nums[i] + 1];
            for (int j = ibeg; j < iend; j++)
                slae.ggl[j] = 0;

        }
    }
}