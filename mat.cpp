/*
 * mat.cpp
 * Copyright (C) 2017 weihao <weihao@weihao-PC>
 *
 * Distributed under terms of the MIT license.
 */

#include <stdio.h>
#include <math.h>
#include "mat.h"
using namespace std;

Mat::Mat(){
    data = nullptr;
    row = 0;
    col = 0;
}
Mat::Mat(double* data, int row, int col){
    this->data = data;
    this->row = row;
    this->col = col;
}
Mat::Mat(int row, int col){
    this->row = row;
    this->col = col;
    double* data = (double*)malloc(row * col * sizeof(double));
    this->data = data;
}
// y, x start from 1
double Mat::GetElement(int y, int x){
    return data[col * (y - 1) + (x - 1)];
}
void Mat::SetElements(double* data, int row, int col){
    this->data = data;
    this->row = row;
    this->col = col;
}
void Mat::SetElement(double var, int y, int x){
    data[col * (y - 1) + (x - 1)] = var;
}
void Mat::ShowElements(){
    for(int i = 0; i < this->row; i++){
        for(int j = 0; j < this->col; j++){
            printf("%f ", this->GetElement(i + 1, j + 1));
        }
        printf("\n");
    }
}

Mat* Mat::GetCol(int col){
    if(col > this->col){
        cout << "FATAL ERROR!" << endl;
        exit(-1);
    }
    Mat* result = new Mat(this->row, 1);
    for(int i = 0; i < this->row; i++){
        result->SetElement(this->GetElement(i + 1, col), i + 1, 1);
    }
    return result;
}

Mat* Mat::GetRow(int row){
    if(row > this->row){
        cout << "FATAL ERROR!" << endl;
        exit(-1);
    }
    Mat* result = new Mat(1, this->col);
    for(int i = 0; i < this->col; i++){
        result->SetElement(this->GetElement(row, i + 1), 1, i + 1);
    }
    return result;
}

// return A.B
Mat* Mul(Mat* A, Mat* B){
    if(A->col != B->row){
        cout << "FATAL ERROR!" << endl;
        exit(-1);
    }

    /*
    Mat* result = (Mat*)malloc(sizeof(Mat));
    double* data = (double*)malloc(A->row * B->col * sizeof(double));
    result->SetElements(data, A->row, B->col);
    */
    Mat* result = new Mat(A->row, B->col);

    for(int i = 0; i < A->row; i++){
        for(int j = 0; j < B->col; j++){
            double tmp = 0;
            for(int k = 0; k < A->col; k++){
                tmp += A->GetElement(i + 1, k + 1) * B->GetElement(k + 1, j + 1);
            }
            result->SetElement(tmp, i + 1, j + 1);
        }
    }
    return result;
}

// return A - B
Mat* Sub(Mat* A, Mat* B){
    if(A->row != B->row || A->col != B->col){
        cout << "FATAL ERROR!" << endl;
        exit(-1);
    }
    Mat* result = new Mat(A->row, A->col);
    for(int i = 0; i < A->row * A->col; i++){
        result->data[i] = A->data[i] - B->data[i];
    }
    return result;
}

// return A^T
Mat* T(Mat* A){
    /*
    Mat* result = (Mat*)malloc(sizeof(Mat));
    double* data = (double*)malloc(A->row * A->col * sizeof(double));
    result->SetElements(data, A->col, A->row);
    */
    Mat* result = new Mat(A->col, A->row);
    for(int i = 0; i < A->row; i++){
        for(int j = 0; j < A->col; j++){
            double tmp = A->GetElement(i + 1, j + 1);
            result->SetElement(tmp, j + 1, i + 1);
        }
    }
    return result;
}

// return Euclidean distance between two vectors
double EDistance(Mat* A, Mat* B){
    if(A->row != B->row || A->col != B->col || A->col != 1){
        cout << "Fatal Error!!" << endl;
        exit(-1);
    }
    double tmp = 0;
    for(int i = 0; i < A-> row; i++){
        tmp += pow((A->data[i] - B->data[i]), 2);
    }
    return sqrt(tmp);
}
