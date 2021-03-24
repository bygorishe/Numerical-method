#pragma once
#define _CRT_SECURE_NO_WARNINGS

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <conio.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

class matrix {
private:
	int n, m, n1, n2, n3, maxk, bSize;
	double* ld, * l1, * l2, * l3, * ud, * u1, * u2, * u3, * di, * f, * x, *t, * xx, eps, w;
public:
	void read();
	double mult(int i, double* x, double* xx);
	void iter();
};

void matrix::read() {
	ifstream fin;
	fin.open("in.txt");
	fin >> n >> m >> eps >> maxk >> w;
	n1 = m + 2;
	n2 = m + 3;
	n3 = m + 4;
	x = new double[n];
	xx = new double[n];
	t = new double[n];
	f = new double [n];
	di = new double[n];
	ld = new double[n - 1];
	l1 = new double[n - n1];
	l2 = new double[n - n2];
	l3 = new double[n - n3];
	ud = new double[n - 1];
	u1 = new double[n - n1];
	u2 = new double[n - n2];
	u3 = new double[n - n3];

	/*for (int i = 0; i < n; i++) {
		x[i] = 0;
		xx[i] = 0;
	}*/
	for (int i = 0; i < n; i++)
		fin >> x[i];
	for (int i = 0; i < n; i++)
		fin >> f[i];
	for (int i = 0; i < n; i++)
		fin >> di[i];

	for (int i = 0; i < n-1; i++)
		fin >> ld[i];
	for (int i = 0; i < n - n1; i++)
		fin >> l1[i];
	for (int i = 0; i < n - n2; i++)
		fin >> l2[i];
	for (int i = 0; i < n - n3; i++)
		fin >> l3[i];

	for (int i = 0; i < n - 1; i++)
		fin >> ud[i];
	for (int i = 0; i < n - n1; i++)
		fin >> u1[i];
	for (int i = 0; i < n - n2; i++)
		fin >> u2[i];
	for (int i = 0; i < n - n3; i++)
		fin >> u3[i];
}

double norm(double *x, int n) {
	double sum = 0;
	for (int i = 0; i < n; i++)
		sum += x[i] * x[i];
	return sqrt(sum);
}

double norm2(double* x, double* v, int n) {
	double sum = 0;
	for (int i = 0; i < n; i++)
		sum += (v[i] - x[i]) * (v[i] - x[i]);
	return sqrt(sum);
}

double matrix::mult(int i, double* x, double* xx) {
	double sum = 0, a = 0;
	for (int j = 0; j < n; j++)
	{
		if (i == j) a = di[i]; // главная диагональ
		else
		{
			int r = j - i;
			if (r == 1) a = ud[i]; //верхние диагонали
			else if (r == n1) a = u1[i];
			else if (r == n2) a = u2[i];
			else if (r == n3) a = u3[i];
			else if (r == -1) a = ld[j]; //нижние диагонали
			else if (r == -n1) a = l1[j];
			else if (r == -n2) a = l2[j];
			else if (r == -n3) a = l3[j];
			else a = 0;
		}
		sum += a * x[j];
	}
	xx[i] = x[i] + (f[i] - sum) * w / di[i];
	return sum;
}

void matrix::iter() {
	double* A = new double[n];
	int k = 0,flag;
	ofstream out;
	out.open("out.txt");
	cout << "1 2\n";
	cin >> flag;
	if (flag == 1) {
		while (k < maxk && (norm2(f, A, n) / norm(f, n)) > eps) {
			for (int i = 0; i < n; i++)
				A[i] = mult(i, x, xx);
			/*for (int i = 0; i < n; i++)
				cout << " " << xx[i] << " ";*/
			swap(xx, x);
			k++;
		}
		double* mas = new double[n];
		for (int i = 1; i <= n; i++)
			mas[i - 1] = i;
		cout << "yacobi:\n";
		double ob = (norm2(mas, x, n) / norm(mas, n)) / (norm2(A, f, n) / norm(f, n));
		cout << ob << "\n" << k << "\n";
		out << ob << "\n" << k << "\n";
		for (int i = 0; i < n; i++) {
			cout << setprecision(15) << x[i] << "  ";
			out << setprecision(15) << x[i] << "\n";
		}
		cout << "\n";
	}
	else {
		while (k < maxk && norm2(f, A, n) / norm(f, n)> eps) {
			for (int i = 0; i < n; i++)
				A[i] = mult(i, x, x);
			/*for (int i = 0; i < n; i++) 
				cout << " " << x[i] << " ";*/
			k++;
		}
		double* mas = new double[n];
		for (int i = 1; i <= n; i++)
			mas[i-1] = i;
		cout << "gausszeidel:\n";
		double ob = (norm2(mas, x, n) / norm(mas, n)) / (norm2(A, f, n) / norm(f, n));
		cout << ob << "\n" << k << "\n";
		out << ob << "\n" << k << "\n";
		for (int i = 0; i < n; i++) {
			cout << setprecision(15) << x[i] << "  ";
			out << setprecision(15) << x[i] << "\n";
		}
		cout << "\n";
	}
	/*for (int i = 0; i < n; i++)
		cout << xx[i] << "   ";*/
}
