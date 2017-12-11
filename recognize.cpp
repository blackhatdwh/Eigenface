/*
 * recognize.cpp
 * Copyright (C) 2017 weihao <weihao@weihao-PC>
 *
 * Distributed under terms of the MIT license.
 */

#include <iostream>
#include <fstream>
#include "mat.h"
#define PICROW 195      // how may rows does a picture has
#define PICCOL 231      // how may cols does a picture has
#define PICMATSIZE 45045        // rows * cols

using namespace std;
int main(){
    // load mat U
    ifstream fp ("U.mat", ios::in);
    int row, col;
    fp >> row >> col;

    Mat* U = new Mat(row, col);

    int index = 0;
    string tmp_str;
    while(getline(fp, tmp_str)){
        U->data[index++] = atof(tmp_str.c_str());
    }
    fp.close();

    // load mat average face
    ifstream fp2 ("average_face.mat", ios::in);
    fp2 >> row;

    Mat* average_face_mat = new Mat(row, 1);

    index = 0;
    while(getline(fp2, tmp_str)){
        average_face_mat->data[index++] = atof(tmp_str.c_str());
    }
    fp2.close();

    // load face_set_mat
    ifstream fp3 ("face_set.mat", ios::in);
    fp3 >> row >> col;

    Mat* face_set_mat = new Mat(row, col);

    index = 0;
    while(getline(fp3, tmp_str)){
        face_set_mat->data[index++] = atof(tmp_str.c_str());
    }
    fp3.close();

    // load a new face
    ifstream fp4 ("centered/18.pgm", ios::in | ios::binary);
    unsigned char new_face[PICMATSIZE];
    for(int i = 0; i < 4; i++){
        getline(fp4, tmp_str);
    }
    fp4.read((char*)new_face, PICMATSIZE);
    fp4.close();
    Mat* new_face_mat = new Mat(PICMATSIZE, 1);
    for(int i = 0; i < PICMATSIZE; i++){
        new_face_mat->data[i] = new_face[i];
    }

    // omiga = Ut(new_face - average_face)
    Mat* Ut = T(U);
    Mat* new_avg_diff = Sub(new_face_mat, average_face_mat);
    Mat* omiga_new_face = Mul(Ut, new_avg_diff);

    // face_set_mat->col represents the number of faces in the training set
    for(int i = 0; i < face_set_mat->col; i++){
        Mat* tmp_mat = face_set_mat->GetCol(i);
        Mat* tmp_mat_2 = Sub(tmp_mat, average_face_mat);
        Mat* tmp_omiga = Mul(Ut, tmp_mat_2);
        cout << EDistance(tmp_omiga, omiga_new_face) << endl;
        delete tmp_mat;
        delete tmp_mat_2;
        delete tmp_omiga;
    }


}

