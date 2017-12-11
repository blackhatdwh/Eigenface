/*
 * mat.h
 * Copyright (C) 2017 weihao <weihao@weihao-PC>
 *
 * Distributed under terms of the MIT license.
 */

#ifndef MAT_H
#define MAT_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

class Mat {
public:
    double* data;
    int row;
    int col;
    Mat();
    Mat(double* data, int row, int col);
    Mat(int row, int col);
    // y, x start from 1
    double GetElement(int y, int x);
    void SetElements(double* data, int row, int col);
    void SetElement(double var, int y, int x);
    void ShowElements();

    Mat* GetCol(int i);
    Mat* GetRow(int i);
};

Mat* Mul(Mat* A, Mat* B);
Mat* Sub(Mat* A, Mat* B);
Mat* T(Mat* A);
double EDistance(Mat* A, Mat* B);

#endif /* !MAT_H */
