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

void gauss(FILE *IN, int& n, real**& A, real*& vector_a, real*& vector_b) {
	fscanf(IN, "%d" , &n);
	vector_b = new real[n];
	vector_a = new real[n];
	A = new real *[n];
	for (int i = 0; i < n; i++)
		A[i] = new real[n+1];

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n+1; j++)
			fscanf(IN, inform, &A[i][j]);

	/*for (int i = 0; i < n; i++)
	{
		printf("\n");
		for (int j = 0; j < n + 1; j++)
			printf(form, A[i][j]);
	}*/
	printf("\n");
	real max, sum, t;
	int imax, jmax;
	int *index = new int[n];

	for (int i = 0; i < n; i++)
		index[i] = i;

	for (int k = 0; k < n-1; k++) {
		max = 0;
		for (int i = k; i < n; i++) 
			for (int j = k; j < n; j++) 
				if (max < abs(A[i][j])) {
					max = abs(A[i][j]);
					imax = i;
					jmax = j;
				}

		for (int i = 0; i < n; i++)
			swap(A[i][k], A[i][jmax]);
		swap(index[k], index[jmax]);
		for (int j = 0; j < n + 1; j++)
			swap(A[k][j], A[imax][j]);

		t = A[k][k];
		if (t != 0) {
			for (int j = k; j < n + 1; j++)
				A[k][j] /= t;
		}
		for (int i = k + 1; i < n; i++) {
			t = A[i][k];
			for (int j = k; j < n + 1; j++)
				A[i][j] -= t * A[k][j];
		}

		//for (int i = 0; i < n; i++)
		//{
		//	printf("\n");
		//	for (int j = 0; j < n + 1; j++)
		//		printf(form, A[i][j]);
		//}

		//printf("\n");
	
	}
	A[n-1][n] /= A[n - 1][n - 1];
	A[n - 1][n - 1] = 1;

	vector_a[n - 1] = A[n - 1][n];
	for (int i = n - 2; i >= 0; i--) {
		sum = 0;
		for (int j = i + 1; j < n; j++)
			sum += A[i][j] * vector_a[j];
		vector_a[i] = A[i][n] - sum;
	}

	for (int i = 0; i < n; i++)
		vector_b[index[i]] = vector_a[i];

	for (int i = 0; i < n; i++)
		printf(form, vector_b[i], "\n");
	printf("\n");
	for (int i = 0; i < n; i++)
		printf(form, real(i)+1-vector_b[i], "\n");

	FILE* OUT = fopen("outg.txt", "w");         //вывод в файл для удобства заполнения таблицы
	for (int i = 0; i < n; i++)
		fprintf(OUT, form"\n", vector_b[i]);
	fprintf(OUT, "\n");
	for (int i = 0; i < n; i++)
		fprintf(OUT, form"\n", real(i) + 1 - vector_b[i]);
}