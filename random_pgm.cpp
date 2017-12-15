/*
 * random_pgm.cpp
 * Copyright (C) 2017 weihao <weihao@weihao-PC>
 *
 * Distributed under terms of the MIT license.
 */

#include <iostream>
#include <fstream>
#include <cstdlib>
using namespace std;
int main(){
    unsigned char a[195*231];
    for(int i = 0; i < 195*231; i++){
        a[i] = rand() % 256;
        //a[i] = 255;
    }
    ofstream wfp ("rand.pgm", ios::out | ios::binary);
    wfp << "P5\n195 231\n255\n";
    wfp.write((char*)a, 195*231);
    wfp.close();
}

