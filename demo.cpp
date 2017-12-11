/*
 * a.cpp
 * Copyright (C) 2017 weihao <weihao@weihao-PC>
 *
 * Distributed under terms of the MIT license.
 */

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#define PICROW 92
#define PICCOL 112
#define PICMATSIZE 10304
#define PICSETSIZE 10

using namespace Eigen;
using namespace std;

class Mat {
public:
    double* data;
    int row;
    int col;
    Mat(){
        data = nullptr;
        row = 0;
        col = 0;
    }
    Mat(double* data, int row, int col){
        this->data = data;
        this->row = row;
        this->col = col;
    }
    // y, x start from 1
    double GetElement(int y, int x){
        return data[col * (y - 1) + (x - 1)];
    }
    void SetElements(double* data, int row, int col){
        this->data = data;
        this->row = row;
        this->col = col;
    }
    void SetElement(double var, int y, int x){
        data[col * (y - 1) + (x - 1)] = var;
    }
    void ShowElements(){
        for(int i = 0; i < this->row; i++){
            for(int j = 0; j < this->col; j++){
                printf("%f ", this->GetElement(i + 1, j + 1));
            }
            printf("\n");
        }
    }
};

Mat* Mul(Mat* A, Mat* B){
    if(A->col != B->row){
        cout << "FATAL ERROR!" << endl;
        exit(0);
    }

    Mat* result = (Mat*)malloc(sizeof(Mat));
    double* data = (double*)malloc(A->row * B->col * sizeof(double));
    result->SetElements(data, A->row, B->col);

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

Mat* T(Mat* A){
    Mat* result = (Mat*)malloc(sizeof(Mat));
    double* data = (double*)malloc(A->row * A->col * sizeof(double));
    result->SetElements(data, A->col, A->row);
    for(int i = 0; i < A->row; i++){
        for(int j = 0; j < A->col; j++){
            double tmp = A->GetElement(i + 1, j + 1);
            result->SetElement(tmp, j + 1, i + 1);
        }
    }
    return result;
}

void ReadOnePic(ifstream* pic, char* result){
    string skip;
    for(int i = 0; i < 3; i++){
        getline(*pic, skip);
    }
    (*pic).read(result, PICMATSIZE);

}

int main(){
    unsigned char face_set[PICSETSIZE][PICMATSIZE];
    for(int i = 0; i < PICSETSIZE; i++){
        char dir[20];
        sprintf(dir, "s/%d.pgm", i);
        ifstream fp (dir, ios::in | ios::binary);
        ReadOnePic(&fp, (char*)face_set[i]);
        fp.close();
    }

    int int_average_face[PICMATSIZE];
    for(int i = 0; i < PICMATSIZE; i++){
        int_average_face[i] = 0;
    }
    for(int i = 0; i < PICMATSIZE; i++){
        for(int j = 0; j < PICSETSIZE; j++){
            int_average_face[i] += face_set[j][i];
        }
        int_average_face[i] /= PICSETSIZE;
    }

    // face set matrix PICMATSIZE * PICSETSIZE
    double* face_set_mat_data = (double*)malloc(PICMATSIZE * PICSETSIZE * sizeof(double));
    Mat face_set_mat (face_set_mat_data, PICMATSIZE, PICSETSIZE);
    for(int i = 0; i < PICSETSIZE; i++){
        for(int j = 0; j < PICMATSIZE; j++){
            face_set_mat.SetElement(face_set[i][j], j + 1, i + 1);
        }
    }

    // diff
    double* diff_mat_data = (double*)malloc(PICMATSIZE * PICSETSIZE * sizeof(double));
    Mat diff_mat (diff_mat_data, PICMATSIZE, PICSETSIZE);
    for(int i = 0; i < PICSETSIZE; i++){
        for(int j = 0; j < PICMATSIZE; j++){
            diff_mat.SetElement(face_set[i][j] - int_average_face[j], j + 1, i + 1);
        }
    }
    // diff^T
    Mat* t_diff_mat = T(&diff_mat);
    // AT * A
    Mat* AtA = Mul(t_diff_mat, &diff_mat);

    //diff_mat.ShowElements();
    AtA->ShowElements();

    // using Eigen to solve eigenvectors
    MatrixXd AtA_eigen = Map<Matrix<double, PICSETSIZE, PICSETSIZE, RowMajor> >(AtA->data);   // load AtA into Eigen matrix
    EigenSolver<MatrixXd> es(AtA_eigen);
    //cout << "eigenvalues: " << endl << es.eigenvalues() << endl;
    //cout << "eigenvectors: " << endl << es.eigenvectors() << endl;
    // end
    
    // calcuate eigenface
    double eigenvector_data[PICSETSIZE];
    Mat eigenvector_mat (eigenvector_data, PICSETSIZE, 1);
    for(int i = 0; i < PICSETSIZE; i++){
        eigenvector_mat.SetElement(es.eigenvectors().col(0)(i).real(), i + 1, 1);
    }
    Mat* eigenface = Mul(&face_set_mat, &eigenvector_mat);
    //face_set_mat.ShowElements();
    //eigenvector_mat.ShowElements();
    //eigenface->ShowElements();


    // scale
    double max_ = 0;
    double min_ = 100;

    for(int i = 0; i < PICMATSIZE; i++){
        if(eigenface->data[i] > max_){
            max_ = eigenface->data[i];
        }
        if(eigenface->data[i] < min_){
            min_ = eigenface->data[i];
        }
    }
    double range = max_ - min_;
    unsigned char eigenface_pgm[PICMATSIZE];
    for(int i = 0; i < PICMATSIZE; i++){
        int tmp = (eigenface->data[i] - min_) / range * 255.0;
        eigenface_pgm[i] = tmp;
    }

    // write to image
    ofstream wfp ("data.pgm", ios::out | ios::binary);
    wfp << "P5\n92 112\n255\n";
    wfp.write((char*)eigenface_pgm, PICMATSIZE);
    wfp.close();
}
