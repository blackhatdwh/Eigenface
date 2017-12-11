/*
 * a.cpp
 * Copyright (C) 2017 weihao <weihao@weihao-PC>
 *
 * Distributed under terms of the MIT license.
 */

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <Eigen/Dense>
#define PICROW 195      // how may rows does a picture has
#define PICCOL 231      // how may cols does a picture has
#define PICMATSIZE 45045        // rows * cols
#define PICSETSIZE 10       // how many pictures are there in the training set
#define BIGGEST_PERCENTAGE 0.7      // retain how many eigenvectors

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
    Mat(int row, int col){
        this->row = row;
        this->col = col;
        double* data = (double*)malloc(row * col * sizeof(double));
        this->data = data;
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
    for(int i = 0; i < 4; i++){
        getline(*pic, skip);
    }
    (*pic).read(result, PICMATSIZE);

}

int main(){
    unsigned char face_set[PICSETSIZE][PICMATSIZE];
    for(int i = 0; i < PICSETSIZE; i++){
        char dir[20];
        sprintf(dir, "centered/%d.pgm", i);
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
    Mat face_set_mat (PICMATSIZE, PICSETSIZE);
    for(int i = 0; i < PICSETSIZE; i++){
        for(int j = 0; j < PICMATSIZE; j++){
            face_set_mat.SetElement(face_set[i][j], j + 1, i + 1);
        }
    }

    // diff
    Mat diff_mat (PICMATSIZE, PICSETSIZE);
    for(int i = 0; i < PICSETSIZE; i++){
        for(int j = 0; j < PICMATSIZE; j++){
            diff_mat.SetElement(face_set[i][j] - int_average_face[j], j + 1, i + 1);
        }
    }
    // diff^T
    Mat* t_diff_mat = T(&diff_mat);
    // AT * A
    Mat* AtA = Mul(t_diff_mat, &diff_mat);


    // using Eigen to solve eigenvectors
    MatrixXd AtA_eigen = Map<Matrix<double, PICSETSIZE, PICSETSIZE, RowMajor> >(AtA->data);   // load AtA into Eigen matrix
    EigenSolver<MatrixXd> es(AtA_eigen);
    //cout << "eigenvalues: " << endl << es.eigenvalues() << endl;
    //cout << "eigenvectors: " << endl << es.eigenvectors() << endl;
    // end


    // sort
    double eigenvalues[PICSETSIZE];
    Mat* raw_eigenvector_vec[PICSETSIZE];       // eigenvectors of AtA
    Mat* eigenvector_vec[PICSETSIZE];           // eigenvectors of AAt
    map<double, Mat*> value_corresponding_vector;
    for(int i = 0; i < PICSETSIZE; i++){
        raw_eigenvector_vec[i] = new Mat (PICSETSIZE, 1);       // use a new matrix(vector) to store one of the raw eigenvectors
        for(int j = 0; j < PICSETSIZE; j++){
            double tmp = es.eigenvectors().col(i)[j].real();
            raw_eigenvector_vec[i]->SetElement(tmp, j + 1, 1);
        }
        eigenvector_vec[i] = Mul(&face_set_mat, raw_eigenvector_vec[i]);        // calculate the eigenvectors of AAt
        eigenvalues[i] = es.eigenvalues()(i).real();        // store the corresponding eigenvalues into an array
        value_corresponding_vector[eigenvalues[i]] = eigenvector_vec[i];        // setup the link between eigenvalue and eigenvector
        delete raw_eigenvector_vec[i];
    }
    vector<double> sorted_eigenvalues (eigenvalues, eigenvalues + PICSETSIZE);
    sort(sorted_eigenvalues.begin(), sorted_eigenvalues.begin() + PICSETSIZE);
    reverse(sorted_eigenvalues.begin(), sorted_eigenvalues.begin() + PICSETSIZE);       // sort eigenvalues in descending order
    // remove those less important eigenvalues
    for(int i = 0; i < PICSETSIZE * (1 - BIGGEST_PERCENTAGE); i++){
        sorted_eigenvalues.pop_back();
    }


    // normalize
    for(int i = 0; i < sorted_eigenvalues.size(); i++){
        double tmp = 0;
        for(int j = 0; j < PICMATSIZE; j++){
            tmp += pow(value_corresponding_vector[sorted_eigenvalues[i]]->GetElement(j + 1, 1), 2);
        }
        tmp = sqrt(tmp);
        for(int j = 0; j < PICMATSIZE; j++){
            double original_value = value_corresponding_vector[sorted_eigenvalues[i]]->GetElement(j + 1, 1);
            value_corresponding_vector[sorted_eigenvalues[i]]->SetElement(original_value / tmp, j + 1, 1);
        }
    }

    // final step, setup U, the matrix constructed by all retained eigenvectors
    // PICMATSIZE * (PICSETSIZE * BIGGEST_PERCENTAGE)
    Mat U (PICMATSIZE, PICSETSIZE * BIGGEST_PERCENTAGE);
    for(int i = 0; i < PICSETSIZE * BIGGEST_PERCENTAGE; i++){
        for(int j = 0; j < PICMATSIZE; j++){
            double tmp = value_corresponding_vector[sorted_eigenvalues[i]]->GetElement(j + 1, 1);
            U.SetElement(tmp, j + 1, i);
        }
    }

    // store U on disk
    ofstream wfp ("U.mat", ios::out);
    string tmp_string = to_string(U.row) + " " + to_string(U.col) + "\n";
    wfp << tmp_string;
    for(int i = 0; i < U.row * U.col; i++){
        wfp << U.data[i] << endl;
    }
    wfp.close();
    
    /*
    // calcuate one eigenface
    Mat eigenvector_mat (PICSETSIZE, 1);
    for(int i = 0; i < PICSETSIZE; i++){
        eigenvector_mat.SetElement(es.eigenvectors().col(1)(i).real(), i + 1, 1);
    }
    Mat* eigenface_mat = Mul(&face_set_mat, &eigenvector_mat);


    // scale
    double max_ = 0;
    double min_ = 100;

    for(int i = 0; i < PICMATSIZE; i++){
        if(eigenface_mat->data[i] > max_){
            max_ = eigenface_mat->data[i];
        }
        if(eigenface_mat->data[i] < min_){
            min_ = eigenface_mat->data[i];
        }
    }
    double range = max_ - min_;
    unsigned char eigenface_pgm[PICMATSIZE];
    for(int i = 0; i < PICMATSIZE; i++){
        int tmp = (eigenface_mat->data[i] - min_) / range * 255.0;
        eigenface_pgm[i] = tmp;
    }

    // write to image
    ofstream wfp ("data.pgm", ios::out | ios::binary);
    wfp << "P5\n195 231\n255\n";
    wfp.write((char*)eigenface_pgm, PICMATSIZE);
    wfp.close();
    */
}
