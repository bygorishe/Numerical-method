#define _CRT_SECURE_NO_WARNINGS

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <conio.h>
#include <iostream>
#include "subroutine.h"
#include "gauss.h"
#include "hylbert.h"
#include "Ak.h"

using namespace std;
//typedef double real;
typedef float real;


int main()
{
    real **matrix_L, *vector_di, * vector_b, * vector_a;
    int n, m, f;

    printf("<1> LLt\n");
    printf("<2> Gauss\n");
    printf("<3> HYlbertM_generation\n");
    printf("<4> AkM_generation\n");
    scanf("%d", &f);
    switch(f){
    case(1): 
        FILE* IN;
        IN = fopen("h.txt", "r");

        read(IN, n, m, matrix_L, vector_di, vector_b);
        decomposition(n, m, matrix_L, vector_di, vector_b);
        vector_y(n, m, matrix_L, vector_di, vector_b);
        vector_x(n, m, matrix_L, vector_di, vector_b);
        break;
    case(2): 
        FILE* INg;
        INg = fopen("ingauss.txt", "r");

        gauss(INg, n, matrix_L, vector_di, vector_b);
        break;
    case(3):
        FILE* INh;
        INh = fopen("h.txt", "w");

        Hylb(INh, n, m, matrix_L, vector_di, vector_b, vector_a);
        break;
    case(4):
        FILE* INk;
        INk = fopen("ak.txt", "w");

        Ak(INk, n, m, matrix_L, vector_di, vector_b, vector_a);
        break;
    }
}

