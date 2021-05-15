#include "knots_triad_coeff.h"
#include "bounds.h"
#include <iomanip>

double determinant(triangle& triad, knot*& knots) { //detD = (x2-x1)(y3-y1)-(x3-x1)(y2-y1)
    //cout << "det";// << knots[1].x;
    return (knots[triad.knot_nums[1]].x - knots[triad.knot_nums[0]].x) *
        (knots[triad.knot_nums[2]].y - knots[triad.knot_nums[0]].y) -
        (knots[triad.knot_nums[2]].x - knots[triad.knot_nums[0]].x) *
        (knots[triad.knot_nums[1]].y - knots[triad.knot_nums[0]].y);
}

void local_G(vector<vector<double>>& G, triangle& triad, double det, double lambda, knot*& knots) {
    vector<vector<double>> a;   //a = D^-1 
    a.resize(3);
    for (int i = 0; i < 3; i++)
        a[i].resize(2);

    // в методичка 
    a[0][0] = (knots[triad.knot_nums[1]].y - knots[triad.knot_nums[2]].y) / (det);  //y2-y3 / det
    a[0][1] = (knots[triad.knot_nums[2]].x - knots[triad.knot_nums[1]].x) / (det);  //x3-x2
    a[1][0] = (knots[triad.knot_nums[2]].y - knots[triad.knot_nums[0]].y) / (det);  //y3-y1
    a[1][1] = (knots[triad.knot_nums[0]].x - knots[triad.knot_nums[2]].x) / (det);  //x1-x3
    a[2][0] = (knots[triad.knot_nums[0]].y - knots[triad.knot_nums[1]].y) / (det);  //y1-y2
    a[2][1] = (knots[triad.knot_nums[1]].x - knots[triad.knot_nums[0]].x) / (det);  //x2-x1

    cout << endl << "a" << endl;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++)
            cout << a[i][j] << " ";
        cout << endl;
    }

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            G[i][j] = lambda * abs(det) * (a[i][0] * a[j][0] + a[i][1] * a[j][1]) / 2;
            //cout << "a G" << a[i][0] << a[i][1] << "\n";
        }
    
    cout << "matrix G\n";
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++)
            cout << G[i][j] << " ";
        cout << "\n";
    }
    cout << "\n";
}

void local_M(vector<vector<double>>& M, double det, double gamma) {
    M[0][0] = M[1][1] = M[2][2] = gamma * abs(det) / 12.;
    M[0][1] = M[1][0] = M[0][2] = M[2][0] = M[2][1] = M[1][2] = gamma * abs(det) / 24.;

    /*
    cout << "matrix M\n";
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++)
        {
            cout << M[i][j] << " ";
        }
        cout << "\n";
    }*/
}

void local_A(vector<local*>& local_A_list, vector<triangle*>& triangle_list, vector<int> lambda_vector, knot* knots, vector<int> gamma_v, vector<int> f) {
    local_A_list.resize(triangle_num);
    int iA = 0;
    for (vector <triangle*>::iterator iter = triangle_list.begin(); iter != triangle_list.end(); iter++, iA++) {
        triangle* triad = *iter;
        local* local_A = new local;

        cout << "knot numbers \n";
        for (int j = 0; j < 4; j++)
            cout << triad->knot_nums[j] << " ";
        cout << "\n";

        cout << lambda_vector[triad->knot_nums[3]];
        double lambda = func_lambda(lambda_vector[triad->knot_nums[3]], knots);
        cout << "lambda " << lambda << " ";
        double det = determinant(*triad, knots);
        cout << "determinant = " << det << "\n";

        //vector<vector<double>>

        vector<vector<double>> M, G;
        M.resize(3);
        G.resize(3);
        for (int i = 0; i < 3; i++) {
            M[i].resize(3);
            G[i].resize(3);
        }

        double gamma = func_gamma(gamma_v[triad->knot_nums[3]], knots);
        cout << "gamma " << gamma;
        
        local_M(M, det, gamma);
        local_G(G, *triad, det, lambda, knots);

        local_A->knot_nums.resize(3);
        for (int i = 0; i < 3; i++)
            //local_A->knot_nums[i].insert(triad->knot_nums[i]);
            local_A->knot_nums[i] = triad->knot_nums[i];

        //сортируем для дальнейшего использования в построение глобальной матрицы вместо L(i)
        //sort(local_A->knot_nums.begin(), local_A->knot_nums.end());

        local_A->A.resize(3);
        for (int i = 0; i < 3; i++) {
            local_A->A[i].resize(3);
            for (int j = 0; j < 3; j++)
                local_A->A[i][j] = G[i][j] + M[i][j];
            //local_A->A[i][j] = iA + 1; //проверка для глобальной матрицы
        }

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++)
                cout << local_A->A[i][j] << " ";
            cout << endl;
        }


        local_A->b.resize(3);
        //локальный вектор b = f * C
        double f1 = func_f(f[triad->knot_nums[3]], triad->knot_nums[0], knots);
        double f2 = func_f(f[triad->knot_nums[3]], triad->knot_nums[1], knots);
        double f3 = func_f(f[triad->knot_nums[3]], triad->knot_nums[2], knots);

        cout << triad->knot_nums[0] << "   " << triad->knot_nums[1] << "   " << triad->knot_nums[2] << "   " << f[triad->knot_nums[3]] << endl;
        cout << f1 << "   " << f2 << "   " << f3 << "   " << f[triad->knot_nums[3]] << endl;
        local_A->b[0] = abs(det) * (2. * f1 + f2 + f3) / 24.;
        local_A->b[1] = abs(det) * (f1 + 2. * f2 + f3) / 24.;
        local_A->b[2] = abs(det) * (f1 + f2 + 2. * f3) / 24.;

        // //проверка для глобальной матрицы
        /*local_A->b[0] = iA + 1;
        local_A->b[1] = iA + 1;
        local_A->b[2] = iA + 1;*/

        cout << "local b " << endl;
        cout << local_A->b[0] << " " << local_A->b[1] << " " << local_A->b[2] << " " << endl;

        local_A_list[iA] = local_A;
    }
}

void global_A(vector<local*>& Local_A, vector<set<int>>& L, SLAE& slae) {
    for (int k = 0; k < num_knots; k++) {
        slae.di[k] = 0.;
        slae.b[k] = 0.;
    }
    cout << endl;

    //алгоритм из методички
    for (vector <local*>::iterator iter = Local_A.begin(); iter != Local_A.end(); iter++) {
        local* A_iter = *iter;
        //local* A_iter = Local_A[0];

        //заносим высе диагональные элементы
        for (int k = 0; k < 3; k++) {
            slae.di[A_iter->knot_nums[k]] += A_iter->A[k][k];
            cout << "di" << A_iter->knot_nums[k] << "=" << slae.di[A_iter->knot_nums[k]] << endl;
            slae.b[A_iter->knot_nums[k]] += A_iter->b[k];
            cout << "b" << A_iter->knot_nums[k] << "=" << slae.b[A_iter->knot_nums[k]] << endl;

        }

        //начинаем цикл по строкам нижнего
        for (int i = 0; i < 3; i++) {
            //устанавливаем начальное значение нижней границы поиска
            int ibeg = slae.ig[A_iter->knot_nums[i]];
            cout << "knot_nums: " << A_iter->knot_nums[i] << endl;
            cout << "ibeg=" << ibeg << endl;

            for (int j = 0; j < i; j++) { // do j=1,i-1
                int iend = slae.ig[A_iter->knot_nums[i] + 1] - 1;
                cout << "iend=" << iend << endl;

                //ploblems
                //while на  for
                while (slae.jg[ibeg] != A_iter->knot_nums[j]) {
                    int ind = (ibeg + iend) / 2;
                    if (slae.jg[ind] < A_iter->knot_nums[j])
                        ibeg = ind + 1;
                    else
                        iend = ind;
                }
                slae.ggu[ibeg] += A_iter->A[j][i];
                slae.ggl[ibeg] += A_iter->A[i][j];
                ibeg++;
            }
        }
    }
}

void create_L(vector<set<int>> L, vector<triangle*>& triangle_list, SLAE*& slae) {
    int a[3];
    for (vector <triangle*>::iterator iter = triangle_list.begin(); iter != triangle_list.end(); iter++) {
        triangle* triad = *iter;
        a[0] = triad->knot_nums[0];
        a[1] = triad->knot_nums[1];
        a[2] = triad->knot_nums[2];
        //cout << a[0] << a[1] << a[2] << "   ";
        vector<int> abc(3);
        abc.insert(abc.begin(), a, a + 3);
        //cout << abc[0] << abc[1] << abc[2] << "   ";
        L[abc[2]].insert(abc[1]);
        L[abc[2]].insert(abc[0]);

        L[abc[1]].insert(abc[0]);
    }

    n = 0;
    for (int i = 0; i < num_knots; i++)
        n += L[i].size();
    //cout << "nnn" << n;

    slae->jg.resize(n);
    slae->ggu.resize(n);
    slae->ggl.resize(n);
    for (int i = 0; i < n; i++) {
        slae->ggu[i] = 0.;
        slae->ggl[i] = 0.;
    }

    int i = 0;
    for (vector<set<int>>::iterator it = L.begin(); it != L.end(); it++)
        for (set<int>::const_iterator cit = it->begin(); cit != it->end(); cit++, i++)
            slae->jg[i] = (*cit);

    for (int i = 0; i < slae->jg.size(); i++)
        cout << slae->jg[i] << " ";
    cout << endl;

    slae->ig[0] = 0;
    for (int i = 1; i < num_knots + 1; i++)
        slae->ig[i] = slae->ig[i - 1] + L[i - 1].size();

    for (int i = 0; i < num_knots + 1; i++)
        cout << slae->ig[i] << " ";
}

void A_mult(SLAE& slae, vector<double> f, vector<double>& res) {
    //эх, пв пв
    for (int i = 0; i < num_knots; i++)
        res[i] = slae.di[i] * f[i];

    for (int i = 0; i < num_knots; i++) {
        for (int k = slae.ig[i]; k < slae.ig[i + 1]; k++) {
            res[i] += slae.ggl[k] * f[slae.jg[k]];
            res[slae.jg[k]] += slae.ggu[k] * f[i];
        }
    }

    /*for (int i = 0; i < num_knots; i++)
        cout << res[i] << "  ";*/
}

void Conjugate_Gradient_Method_LOS(SLAE& slae, knot* knots) {
    double alpha, betta, residual, scalar_p, sqrt_scalar_b, eps = 1E-15;
    vector<double> z, r, x, p, temp;
    temp.resize(num_knots);
    x.resize(num_knots);
    p.resize(num_knots);
    z.resize(num_knots);
    r.resize(num_knots);
    for (int i = 0; i < num_knots; i++) {
        x[i] = 0.;
        z[i] = 0.;
        r[i] = 0.;
        p[i] = 0.;
        temp[i] = 0.;
    }
    //r0 =f - A * x0
    //z0 = r0
    //p0 = A * z0
    A_mult(slae, x, temp);
    r = slae.b + (-1) * temp;
    z = r;
    A_mult(slae, z, temp);
    p = temp;
    sqrt_scalar_b = sqrt(scalar(slae.b, slae.b));
    residual = sqrt(scalar(r, r)) / sqrt_scalar_b;
    for (int k = 0; k < 100000 && residual > eps; k++) {
        scalar_p = scalar(p, p);
        alpha = scalar(p, r) / scalar_p;
        x += alpha * z;
        r += -alpha * p;
        A_mult(slae, r, temp);
        betta = -scalar(p, temp) / scalar_p;
        z = r + betta * z;
        p = temp + betta * p;
        residual = sqrt(scalar(r, r)) / sqrt_scalar_b;
        for (int i = 0; i < num_knots; i++) 
            cout << setprecision(15) << x[i] << " ";
        cout << endl;
    }

    vector<double> u(num_knots), ururur(num_knots);
    ofstream file("result.txt");
    for (int i = 0; i < num_knots; i++) {
        cout << x[i] << " ";
        file << setprecision(15) << x[i] << endl;
    }
    file << endl;

    for (int i = 0; i < num_knots; i++) {
        u[i] = pow(knots[i].x, 4);//= knots[i].x;
        file << setprecision(15) << u[i] << endl;
    }
    file << endl;

    for (int i = 0; i < num_knots; i++) {
        ururur[i] = x[i] - u[i];
        file << setprecision(15) << abs(x[i] - u[i]) << endl;
    }

    double t = 0;
    for (int i = 0; i < num_knots; i++) {
       // cout << ururur[i] << "  ";
        t += ururur[i] * ururur[i];
       // cout << t << "    ";
    }
       
    cout << sqrt(t);
    
    cout << endl << knots[4].x << " " << knots[5].x << " " << knots[7].x << " " << endl;
    cout << endl << x[4] << " " << x[5] << " " << x[7] << " " << endl;
    cout << "x/y=3.4142...." << endl;
    double asasasa = x[0] * 0.5 + x[1] * 0.25 + x[3] * 0.25;
    cout << setprecision(15) << asasasa;
    //cout << knots[4].x * 0.5 + knots[5].x * (1 - sqrt(2) / 2) + knots[7].x * (1 - sqrt(2) / 2);
}



int main(void) {
    double betta = 1;

    SLAE* slae;
    slae = new SLAE;
    knot* knots = read_knots();

    slae->ig.resize(num_knots + 1);
    slae->di.resize(num_knots);
    slae->b.resize(num_knots);

    //вектор элементов(треуголников)
    vector<triangle*> triangles;
    create_triangles(triangles);
    /*for (int i = 0; i < triangle_num; i++) {
        triangle* triad = triangles[i];
        for (int j = 0; j < 3; j++) {
            cout << triad->knot_nums[j] << " ";
        }
        cout << "\n";
    }*/

    //вектора, где размероность = кол-во подобластей, а эл - номер нужной функции
    vector<int> lambda_vector;
    vector<int> f_vector;
    vector<int> gamma;

    create_lambda_f_gamma(lambda_vector, f_vector, gamma);

    //вектор связностей
    vector<set<int>> L(num_knots, set<int>()); //вектор сетов размерности = количество узлов
    create_L(L, triangles, slae);

    /*for (int i = 0; i < num_knots; i++)
        cout << knots[i].x << " " << knots[i].y << "\n";


    for (int i = 0; i < triangle_num; i++) {
        triangle* triad = triangles[i];
        for (int j = 0; j < 3; j++) {
            cout << triad->knot_nums[j] << " ";
        }
        cout << "\n";
    }/**/

    //вектор локальных матриц
    vector<local*> Local_A;
    local_A(Local_A, triangles, lambda_vector, knots, gamma, f_vector);

    //сборка глобальной матрицы
    global_A(Local_A, L, *slae);

    cout << "ggu" << endl;
    for (int i = 0; i < n; i++)
        cout << setprecision(15) << slae->ggu[i] << " ";
    cout << endl;
    cout << "ggl" << endl;
    for (int i = 0; i < n; i++)
        cout << setprecision(15) << slae->ggl[i] << " ";
    cout << endl;
    cout << "di" << endl;
    for (int i = 0; i < num_knots; i++)
        cout << setprecision(15) << slae->di[i] << " ";
    cout << endl;
    cout << "b" << endl;
    for (int i = 0; i < num_knots; i++)
        cout << setprecision(15) << slae->b[i] << " ";
    cout << endl;

    //читаем краевые и вносим в глобальную   
    vector<local*> bounds_1, bounds_2, bounds_3;
    ifstream input_1("boundary_conditions_1.txt");
    ifstream input_2("boundary_conditions_2.txt");
    ifstream input_3("boundary_conditions_3.txt");
    
    read_bounds_2_3(bounds_2, input_2, num_bounds_2);
    read_bounds_2_3(bounds_3, input_3, num_bounds_3);
    read_bounds_1(bounds_1, input_1, num_bounds_1);

    build_bound(bounds_2, 2, knots, betta);
    build_bound(bounds_3, 3, knots, betta);
    
    //use_bounds(bounds_2, L, *slae, 2, knots);
    //use_bounds(bounds_3, L, *slae, 3, knots);

    //for (int i = 0; i < 2; i++) 
      //  cout << bounds_2[i]->knot_nums[0] << " " << bounds_2[i]->knot_nums[1] << " ";


    //решаем слау и се
    cout << endl;
    /*vector<double> x1, x2;
    x1.resize(num_knots);
    x2.resize(num_knots);
    for (int i = 0; i < num_knots; i++) {
        x1[i] = 1;
        x2[i] = 0;
    }

    slae->b[0] = 12;
    slae->b[1] = 21;
    slae->b[2] = 6;
    slae->b[3] = 9;
    slae->b[4] = 24;
    slae->b[5] = 18;*/

    
    ///*cout << "ggu" << endl;
    //for (int i = 0; i < 9; i++)
    //    cout << slae->ggu[i] << " ";
    //cout << endl;
    //cout << "ggl" << endl;
    //for (int i = 0; i < 9; i++)
    //    cout << slae->ggl[i] << " ";
    //cout << endl;
    //cout << "di" << endl;
    //for (int i = 0; i < 6; i++)
    //    cout << slae->di[i] << " ";
    //cout << endl;
    //cout << "b" << endl;
    //for (int i = 0; i < 6; i++)
    //    cout << slae->b[i] << " ";
    //cout << endl;*/

    use_first_bounds(bounds_1, L, *slae, knots);
    
    cout << "hoi" << endl;
    cout << "ggu" << endl;
    for (int i = 0; i < n; i++)
        cout << setprecision(15) << slae->ggu[i] << " ";
    cout << endl;
    cout << "ggl" << endl;
    for (int i = 0; i < n; i++)
        cout << setprecision(15) << slae->ggl[i] << " ";
    cout << endl;
    cout << "di" << endl;
    for (int i = 0; i < num_knots; i++)
        cout << setprecision(15) << slae->di[i] << " ";
    cout << endl;
    cout << "b" << endl;
    for (int i = 0; i < num_knots; i++)
        cout << setprecision(15) << slae->b[i] << " ";
    cout << endl;
    

    Conjugate_Gradient_Method_LOS(*slae, knots);
    cout << endl;

    cout << "\nHello World!\n";
}