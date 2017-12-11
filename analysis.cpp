/*
 * a.cpp
 * Copyright (C) 2017 weihao <weihao@weihao-PC>
 *
 * Distributed under terms of the MIT license.
 */

#include <stdio.h>
#include <math.h>
#include "mat.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <Eigen/Dense>
#define PICROW 195      // how may rows does a picture has
#define PICCOL 231      // how may cols does a picture has
#define PICMATSIZE 45045        // rows * cols
#define PICSETSIZE 15       // how many pictures are there in the training set
#define BIGGEST_PERCENTAGE 0.7      // retain how many eigenvectors

using namespace Eigen;
using namespace std;

void ReadOnePic(ifstream* pic, char* result){
    string skip;
    for(int i = 0; i < 4; i++){
        getline(*pic, skip);
    }
    (*pic).read(result, PICMATSIZE);

}

int main(){
    // load PICSETSIZE faces into a 2d array face set
    unsigned char face_set[PICSETSIZE][PICMATSIZE];
    for(int i = 0; i < PICSETSIZE; i++){
        char dir[20];
        sprintf(dir, "testset/%d.pgm", i);
        ifstream fp (dir, ios::in | ios::binary);
        ReadOnePic(&fp, (char*)face_set[i]);
        fp.close();
    }
    
    // convert face set into a matrix
    // face set matrix PICMATSIZE * PICSETSIZE
    Mat face_set_mat (PICMATSIZE, PICSETSIZE);
    for(int i = 0; i < PICSETSIZE; i++){
        for(int j = 0; j < PICMATSIZE; j++){
            face_set_mat.SetElement(face_set[i][j], j + 1, i + 1);
        }
    }

    // calculate average face
    int average_face[PICMATSIZE];
    Mat average_face_vec (PICMATSIZE, 1);
    for(int i = 0; i < PICMATSIZE; i++){        // all set to 0
        average_face[i] = 0;
    }
    for(int i = 0; i < PICMATSIZE; i++){
        for(int j = 0; j < PICSETSIZE; j++){
            average_face[i] += face_set[j][i];
        }
        average_face[i] /= PICSETSIZE;
        average_face_vec.SetElement(average_face[i], i + 1, 1);
    }


    // diff(A)
    Mat* diff_mat = new Mat(PICMATSIZE, PICSETSIZE);
    for(int i = 0; i < PICSETSIZE; i++){
        for(int j = 0; j < PICMATSIZE; j++){
            diff_mat->SetElement(face_set_mat.GetElement(j + 1, i + 1) - average_face_vec.GetElement(j + 1, 1), j + 1, i + 1);
        }
    }
    // diff^T
    Mat* t_diff_mat = T(diff_mat);
    // AT * A
    Mat* AtA = Mul(t_diff_mat, diff_mat);


    // using Eigen
    MatrixXd AtA_eigen = Map<Matrix<double, PICSETSIZE, PICSETSIZE, RowMajor> >(AtA->data);   // load AtA into Eigen matrix
    EigenSolver<MatrixXd> es(AtA_eigen);
    //cout << "eigenvalues: " << endl << es.eigenvalues() << endl;
    //cout << "eigenvectors: " << endl << es.eigenvectors() << endl;
    // end


    // calculate eigenvalues and eigenvectors
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
        //eigenvector_vec[i] = Mul(&face_set_mat, raw_eigenvector_vec[i]);        // calculate the eigenvectors of AAt
        eigenvector_vec[i] = Mul(diff_mat, raw_eigenvector_vec[i]);        // calculate the eigenvectors of AAt
        eigenvalues[i] = es.eigenvalues()(i).real();        // store the corresponding eigenvalues into an array
        value_corresponding_vector[eigenvalues[i]] = eigenvector_vec[i];        // setup the link between eigenvalue and eigenvector
        delete raw_eigenvector_vec[i];
    }

    // sort
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
    Mat U (PICMATSIZE, sorted_eigenvalues.size());
    for(int i = 0; i < sorted_eigenvalues.size(); i++){
        for(int j = 0; j < PICMATSIZE; j++){
            double tmp = value_corresponding_vector[sorted_eigenvalues[i]]->GetElement(j + 1, 1);
            U.SetElement(tmp, j + 1, i + 1);
        }
    }
    U.ShowElements();


    
    // store average face on disk
    ofstream wfp ("average_face.mat", ios::out);
    string tmp_string = to_string(average_face_vec.row) + " " + to_string(average_face_vec.col) + "\n";
    wfp << tmp_string;
    for(int i = 0; i < average_face_vec.row * average_face_vec.col; i++){
        wfp << average_face_vec.data[i] << endl;
    }
    wfp.close();

    // store U on disk
    ofstream wfp2 ("U.mat", ios::out);
    tmp_string = to_string(U.row) + " " + to_string(U.col) + "\n";
    wfp2 << tmp_string;
    for(int i = 0; i < U.row * U.col; i++){
        wfp2 << U.data[i] << endl;
    }
    wfp2.close();

    // store face_set_mat on disk
    ofstream wfp3 ("face_set.mat", ios::out);
    tmp_string = to_string(face_set_mat.row) + " " + to_string(face_set_mat.col) + "\n";
    wfp3 << tmp_string;
    for(int i = 0; i < face_set_mat.row * face_set_mat.col; i++){
        wfp3 << face_set_mat.data[i] << endl;
    }
    wfp3.close();

    




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
