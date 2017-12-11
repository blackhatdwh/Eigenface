/*
 * mat_mul.cpp
 * Copyright (C) 2017 weihao <weihao@weihao-PC>
 *
 * Distributed under terms of the MIT license.
 */

#include <iostream>
using namespace std;
class Mat {
    char* data;
public:
    int row;
    int col;
    Mat(){
        data = nullptr;
        row = 0;
        col = 0;
    }
    Mat(char* data, int row, int col){
        this->data = data;
        this->row = row;
        this->col = col;
    }
    // y, x start from 1
    char GetElement(int y, int x){
        return data[col * (y - 1) + (x - 1)];
    }
    void SetElements(char* data, int row, int col){
        this->data = data;
        this->row = row;
        this->col = col;
    }
    void SetElement(char var, int y, int x){
        data[col * (y - 1) + (x - 1)] = var;
    }
};

class Mat* Mul(class Mat* A, class Mat* B){
    if(A->col != B->row){
        cout << "FATAL ERROR!" << endl;
        exit(0);
    }

    class Mat* result = (class Mat*)malloc(sizeof(class Mat));
    char* data = (char*)malloc(A->row * B->col * sizeof(char));
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
    char* data = (char*)malloc(A->row * A->col * sizeof(char));
    result->SetElements(data, A->col, A->row);
    for(int i = 0; i < A->row; i++){
        for(int j = 0; j < A->col; j++){
            int tmp = A->GetElement(i + 1, j + 1);
            result->SetElement(tmp, j + 1, i + 1);
        }
    }
    return result;
}

int main(){
    char a_data[] = {
        1,2,3,4,
        5,6,7,8,
        1,2,3,4,
    };
    char b_data[] = {
        1,2,3,
        1,2,3,
        1,2,3,
        1,2,3,
    };
    Mat a (a_data, 3, 4);
    Mat b (b_data, 4, 3);

    //Mat* c = Mul(&a, &b);
    Mat* c = T(&b);
    printf("%d\n", c->GetElement(1, 4));
}

