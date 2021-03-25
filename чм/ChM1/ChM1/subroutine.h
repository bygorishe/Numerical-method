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

void read(FILE *IN, int &n, int& m, real** &matrix_L, real*& vector_di, real* &vector_b) {
	fscanf(IN, "%d" "%d", &n, &m);
	vector_b = new real[n];
	vector_di = new real[n];
	matrix_L = new real*[n];
	for (int i = 0; i < n; i++)
		matrix_L[i] = new real [m];
	for (int i = 0; i < n; i++)
		fscanf(IN, inform, &vector_b[i]);
	for (int i = 0; i < n; i++)
		fscanf(IN, inform, &vector_di[i]);
	for (int i = 0; i < n; i++) 
		for (int j = 0; j < m; j++)
			fscanf(IN, inform, &matrix_L[i][j]);	
}

void decomposition(int& n, int& m, real**& matrix_L, real*& vector_di, real*& vector_b) {
	for (int i = 0; i < n; i++) {
		int f = 0;                 //отвечает за смещение от первого индекса если ширина ленты < n-1
		int t = 0;
		int j =m- i;
		if (j < 0) j = 0;
		if (m - i < 0) f = i - m;      //смещение от нач индекса строки
		while ( j <= m-1) {
			real sum2 = 0;
			int tt = m;
			for (int k = j-1; k >= 0; k--) 
				sum2 += matrix_L[i][k] * matrix_L[t+f][--tt];  
			matrix_L[i][j] = (matrix_L[i][j] - sum2) / vector_di[f+t];
			j++;
			t++;
		}
		real sum1 = 0;
		for (int k = 1; (k <= i) && (m-k>=0); k++) {
			sum1 += matrix_L[i][m-k] * matrix_L[i][m-k];
		}
		vector_di[i] = sqrt(vector_di[i] - sum1);
	}

	//for (int i = 0; i < n; i++)
	//	printf(form, vector_di[i], "\n");
	//for (int i = 0; i < n; i++) {                                         
	//	printf("\n");
	//	for (int j = 0; j < m; j++)
	//		printf(form, matrix_L[i][j], "\n");
	//}
	//printf("\n");
	//printf("\n");
}

void vector_y(int& n, int& m, real**& matrix_L, real*& vector_di, real*& vector_b) {  
	for (int i = 0; i < n; i++) {
		int g = i;
		for (int k = m - 1; (k > m - 1 - i) && (k >= 0); k--) {
			//printf(form, vector_b[i], "   ");
			vector_b[i] -= /*(double)*/matrix_L[i][k] * /*(double)*/vector_b[--g];                        //(double) для скалярного произведения в дв точности                       
			//printf(form, vector_b[i] , "   ");
			//printf("\n");
		}
		vector_b[i] /= vector_di[i];
		printf("\n");
	}
	for (int i = 0; i < n; i++) 
		printf(form, vector_b[i], "\n");
	printf("\n");
}

void vector_x(int& n, int& m, real**& matrix_L, real*& vector_di, real*& vector_b) {                   
	for (int i = n-1 ; i >= 0; i--) {
		int j = m;
		vector_b[i] /= vector_di[i];
		for (int k = i - 1; (k >= 0) && (j >= 0); k--) 
			vector_b[k] -= /*(double)*/matrix_L[i][--j] * /*(double)*/vector_b[i];
	}

	printf("\n");
	for (int i = 0; i < n; i++)
		printf(form, vector_b[i], "\n");
	printf("\n");
	printf("\n");
	for (int i = 0; i < n; i++)
		printf(form, real(i)+1-vector_b[i], "\n");


	FILE* OUT = fopen("out.txt", "w");         //вывод в файл для удобства заполнения таблицы
	for (int i = 0; i < n; i++)
		fprintf(OUT,form"\n", vector_b[i] );
	fprintf(OUT,"\n");
	for (int i = 0; i < n; i++)
		fprintf(OUT, form"\n", real(i) + 1 - vector_b[i]);
}