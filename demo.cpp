/*
 * a.cpp
 * Copyright (C) 2017 weihao <weihao@weihao-PC>
 *
 * Distributed under terms of the MIT license.
 */

#include <stdio.h>
#include <iostream>
#include <fstream>
#define PICROW 92
#define PICCOL 112
#define PICMATSIZE 10304
#define PICSETSIZE 10

using namespace std;

class Mat {
    int* data;
public:
    int row;
    int col;
    Mat(){
        data = nullptr;
        row = 0;
        col = 0;
    }
    Mat(int* data, int row, int col){
        this->data = data;
        this->row = row;
        this->col = col;
    }
    // y, x start from 1
    char GetElement(int y, int x){
        return data[col * (y - 1) + (x - 1)];
    }
    void SetElements(int* data, int row, int col){
        this->data = data;
        this->row = row;
        this->col = col;
    }
    void SetElement(int var, int y, int x){
        data[col * (y - 1) + (x - 1)] = var;
    }
};

class Mat* Mul(class Mat* A, class Mat* B){
    if(A->col != B->row){
        cout << "FATAL ERROR!" << endl;
        exit(0);
    }

    class Mat* result = (class Mat*)malloc(sizeof(class Mat));
    int* data = (int*)malloc(A->row * B->col * sizeof(int));
    result->SetElements(data, A->row, B->col);

    for(int i = 0; i < A->row; i++){
        for(int j = 0; j < B->col; j++){
            int tmp = 0;
            for(int k = 0; k < A->col; k++){
                tmp += A->GetElement(i + 1, k + 1) * B->GetElement(k + 1, j + 1);
            }
            result->SetElement(tmp, i + 1, j + 1);
        }
    }
    return result;
}

class Mat* T(class Mat* A){
    class Mat* result = (class Mat*)malloc(sizeof(class Mat));
    int* data = (int*)malloc(A->row * A->col * sizeof(int));
    result->SetElements(data, A->col, A->row);
    for(int i = 0; i < A->row; i++){
        for(int j = 0; j < A->col; j++){
            int tmp = A->GetElement(i + 1, j + 1);
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

    // diff
    int* diff_mat_data = (int*)malloc(PICMATSIZE * PICSETSIZE * sizeof(int));
    Mat diff_mat (diff_mat_data, PICMATSIZE, PICSETSIZE);
    // diff^T
    int* t_diff_mat_data = (int*)malloc(PICMATSIZE * PICSETSIZE * sizeof(int));
    Mat t_diff_mat (t_diff_mat_data, PICSETSIZE, PICMATSIZE);
    // AT * A
    Mat* AtA = Mul(&t_diff_mat, &diff_mat);


    /*
    // write to image
    char char_average_face[PICMATSIZE];
    for(int i = 0; i < PICMATSIZE; i++){
        char_average_face[i] = int_average_face[i];
    }
    ofstream wfp ("data.pgm", ios::out | ios::binary);
    wfp << "P5\n92 112\n255\n";
    wfp.write(char_average_face, PICMATSIZE);
    wfp.close();
    */
}
