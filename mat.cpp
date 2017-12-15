/*
 * mat.cpp
 * Copyright (C) 2017 weihao <weihao@weihao-PC>
 *
 * Distributed under terms of the MIT license.
 */

#include <stdio.h>
#include <cstring>
#include <math.h>
#include "mat.h"
using namespace std;

Mat::Mat(){
    data = NULL;
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


/**
* @brief 求实对称矩阵的特征值及特征向量的雅克比法 
* 利用雅格比(Jacobi)方法求实对称矩阵的全部特征值及特征向量 
* @param pMatrix				长度为n*n的数组，存放实对称矩阵
* @param nDim					矩阵的阶数 
* @param pdblVects				长度为n*n的数组，返回特征向量(按列存储) 
* @param dbEps					精度要求 
* @param nJt					整型变量，控制最大迭代次数 
* @param pdbEigenValues			特征值数组
* @return 
*/
bool JacbiCor(double * pMatrix,int nDim, double *pdblVects, double *pdbEigenValues, double dbEps,int nJt)
{
	for(int i = 0; i < nDim; i ++) 
	{   
		pdblVects[i*nDim+i] = 1.0f; 
		for(int j = 0; j < nDim; j ++) 
		{ 
			if(i != j)   
				pdblVects[i*nDim+j]=0.0f; 
		} 
	} 

	int nCount = 0;		//迭代次数
	while(1)
	{
		//在pMatrix的非对角线上找到最大元素
		double dbMax = pMatrix[1];
		int nRow = 0;
		int nCol = 1;
		for (int i = 0; i < nDim; i ++)			//行
		{
			for (int j = 0; j < nDim; j ++)		//列
			{
				double d = fabs(pMatrix[i*nDim+j]); 

				if((i!=j) && (d> dbMax)) 
				{ 
					dbMax = d;   
					nRow = i;   
					nCol = j; 
				} 
			}
		}

		if(dbMax < dbEps)     //精度符合要求 
			break;  

		if(nCount > nJt)       //迭代次数超过限制
			break;

		nCount++;

		double dbApp = pMatrix[nRow*nDim+nRow];
		double dbApq = pMatrix[nRow*nDim+nCol];
		double dbAqq = pMatrix[nCol*nDim+nCol];

		//计算旋转角度
		double dbAngle = 0.5*atan2(-2*dbApq,dbAqq-dbApp);
		double dbSinTheta = sin(dbAngle);
		double dbCosTheta = cos(dbAngle);
		double dbSin2Theta = sin(2*dbAngle);
		double dbCos2Theta = cos(2*dbAngle);

		pMatrix[nRow*nDim+nRow] = dbApp*dbCosTheta*dbCosTheta + 
			dbAqq*dbSinTheta*dbSinTheta + 2*dbApq*dbCosTheta*dbSinTheta;
		pMatrix[nCol*nDim+nCol] = dbApp*dbSinTheta*dbSinTheta + 
			dbAqq*dbCosTheta*dbCosTheta - 2*dbApq*dbCosTheta*dbSinTheta;
		pMatrix[nRow*nDim+nCol] = 0.5*(dbAqq-dbApp)*dbSin2Theta + dbApq*dbCos2Theta;
		pMatrix[nCol*nDim+nRow] = pMatrix[nRow*nDim+nCol];

		for(int i = 0; i < nDim; i ++) 
		{ 
			if((i!=nCol) && (i!=nRow)) 
			{ 
				int u = i*nDim + nRow;	//p  
				int w = i*nDim + nCol;	//q
				dbMax = pMatrix[u]; 
				pMatrix[u]= pMatrix[w]*dbSinTheta + dbMax*dbCosTheta; 
				pMatrix[w]= pMatrix[w]*dbCosTheta - dbMax*dbSinTheta; 
			} 
		} 

		for (int j = 0; j < nDim; j ++)
		{
			if((j!=nCol) && (j!=nRow)) 
			{ 
				int u = nRow*nDim + j;	//p
				int w = nCol*nDim + j;	//q
				dbMax = pMatrix[u]; 
				pMatrix[u]= pMatrix[w]*dbSinTheta + dbMax*dbCosTheta; 
				pMatrix[w]= pMatrix[w]*dbCosTheta - dbMax*dbSinTheta; 
			} 
		}

		//计算特征向量
		for(int i = 0; i < nDim; i ++) 
		{ 
			int u = i*nDim + nRow;		//p   
			int w = i*nDim + nCol;		//q
			dbMax = pdblVects[u]; 
			pdblVects[u] = pdblVects[w]*dbSinTheta + dbMax*dbCosTheta; 
			pdblVects[w] = pdblVects[w]*dbCosTheta - dbMax*dbSinTheta; 
		} 

	}

	//对特征值进行排序以及重新排列特征向量,特征值即pMatrix主对角线上的元素
	std::map<double,int> mapEigen;
	for(int i = 0; i < nDim; i ++) 
	{   
		pdbEigenValues[i] = pMatrix[i*nDim+i];

		mapEigen.insert(make_pair( pdbEigenValues[i],i ) );
	} 

	double *pdbTmpVec = new double[nDim*nDim];
	std::map<double,int>::reverse_iterator iter = mapEigen.rbegin();
	for (int j = 0; iter != mapEigen.rend(),j < nDim; ++iter,++j)
	{
		for (int i = 0; i < nDim; i ++)
		{
			pdbTmpVec[i*nDim+j] = pdblVects[i*nDim + iter->second];
		}

		//特征值重新排列
		pdbEigenValues[j] = iter->first;
	}

	//设定正负号
	for(int i = 0; i < nDim; i ++) 
	{
		double dSumVec = 0;
		for(int j = 0; j < nDim; j ++)
			dSumVec += pdbTmpVec[j * nDim + i];
		if(dSumVec<0)
		{
			for(int j = 0;j < nDim; j ++)
				pdbTmpVec[j * nDim + i] *= -1;
		}
	}

	memcpy(pdblVects,pdbTmpVec,sizeof(double)*nDim*nDim);
	delete []pdbTmpVec;

	return 1;
}

MyEigen GetEigen(Mat* A)
{
	MyEigen ME;
	double dbEps=0.0000001;
	int nJt=1000;
	ME.EValue=(double *)malloc((A->row)*sizeof(double));
	ME.EVector=(double *)malloc((A->row)*(A->row)*sizeof(double));
	JacbiCor(A->data,A->row,ME.EVector,ME.EValue,dbEps,nJt);
	return ME;
}



