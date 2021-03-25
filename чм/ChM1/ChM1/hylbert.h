#pragma once
#define _CRT_SECURE_NO_WARNINGS

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <conio.h>
#include <iostream>

using namespace std;
//typedef double real;
//
//#define inform "%lf"
//#define form "%-20.15f"

typedef float real;

#define inform "%f"
#define form "%-10f"

void Hylb(FILE* IN, int& n, int& m, real**& matrix_L, real*& vector_di, real*& vector_b, real*& vector_x) {
	printf("n?");
	scanf("%d", &n);
	m = n - 1;

	vector_b = new real[n];
	vector_x = new real[n];
	vector_di = new real[n];
	matrix_L = new real * [n];
	for (int i = 0; i < n; i++)
		matrix_L[i] = new real[m];

	for (int i = 1; i <= n; i++) {
		vector_x[i - 1] = i;
		vector_di[i - 1] = 1 / (2*real(i)-1);
		int k = 1;
		for (int j = 0; j < m; j++) {
			if (j >= m - i + 1) {
				matrix_L[i - 1][j] = 1 / (real(i)+real(k)-1) ;
				k++;
			}
			else
				matrix_L[i - 1][j] = 0;
		}
	}


	for (int i = 0; i < n; i++) {
		real sum = 0;
		int g = m;
		sum = vector_di[i] * vector_x[i];
		for (int j = i + 1; j < n && m-g>=0; j++) {
			sum += matrix_L[j][--g] * vector_x[j];
		}
		g = m;
		for (int k = i-1; k >= 0; k--) {
			sum += matrix_L[i][--g] * vector_x[k];
		}
		vector_b[i] = sum;
	}


	fprintf(IN, "%d\n" "%d\n", n, m);
	for (int i = 0; i < n; i++)
		fprintf(IN, form, vector_b[i]);
	fprintf(IN,"\n");
	for (int i = 0; i < n; i++)
		fprintf(IN, form, vector_di[i]);
	fprintf(IN, "\n");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++)
			fprintf(IN, form, matrix_L[i][j]);
		fprintf(IN, "\n");
	}
}