#include "setFunctions.h"
#include <math.h>
#pragma once

void copy(double a[nr][nw], double b[nr][nw]){
	for (int i = 0; i < nr; i++){
		for (int j = 0; j < nw; j++){
			a[i][j] = b[i][j];
		}
	}
}

void newKt(){
	//From Rmax  r integration
	for (int j = 0; j < nw; j++){
		for (int i = nr - 1; i>0; i--){
			kt[i - 1][j] = kt[i][j]-dr*(coeff(2, 2, i, j)*kt[i][j] + coeff(2, 3, i, j));
		}
	}
}

void newKr(){
	//From w=pi/2 midpoint  w integration
	int midW = (nw - 1) / 2;
	for (int i = 0; i>0; i--){
		for (int j = midW; j > 0; j--){
			kr[i][j-1] = kr[i][j]+dw*(coeff(3, 2, i, j)*kr[i][j] + coeff(3, 3, i, j))/coeff(3,1,i,j);
		}
		for (int j = midW; j < nw; j++){
			kr[i][j + 1] = kr[i][j] - dw*(coeff(3, 2, i, j)*kr[i][j] + coeff(3, 3, i, j)) / coeff(3, 1, i, j);
		}
	}
}

double F[nw][7];
double ne[nw];
double error[nr][nw];

void setError(int _i){
	double _gr, dgr, ddgr;
	error[_i][0] = 0;
	for (int j = 1; j < nw - 1; j++){
		_gr = gr[_i ][j];
		dgr = (gr[_i ][j + 1] - gr[_i ][j - 1]) / (2.*dw);
		ddgr = (gr[_i ][j + 1] - 2.*_gr + gr[_i ][j - 1]) / (dw*dw);

		error[_i][j] = F[j][1] * (gr[_i][j] - gr[_i - 1][j]) / dr +
			F[j][2] * _gr*ddgr + F[j][3] * _gr*dgr + F[j][4] * dgr*dgr + F[j][5] * _gr*_gr + F[j][6] * _gr;
		
		if (error[_i][0] < error[_i][j] * error[_i][j]){
			error[_i][0] = error[_i][j] * error[_i][j];
		}
	}
}

void setInitSlice(int _i,int opt){
	if (opt == 0){
		for (int j = 0; j < nw; j++){
			ne[j] = gr[_i][j];
		}
	}
	if (opt == 1){
		double _gr,dgr, ddgr;
		ne[0] = gr[_i][0];
		ne[nw-1] = gr[_i][nw-1];
		for (int j = 1; j < nw-1; j++){
			_gr = gr[_i - 1][j];
			dgr = (gr[_i - 1][j + 1] - gr[_i - 1][j - 1]) / (2.*dw);
			ddgr = (gr[_i - 1][j + 1] - 2.*_gr + gr[_i - 1][j - 1]) / (dw*dw);

			gr[_i][j] = gr[_i - 1][j] -
				(dr / F[j][1])*(
				F[j][2] * _gr*ddgr + F[j][3] * _gr*dgr + F[j][4] * dgr*dgr + F[j][5] * _gr*_gr + F[j][6]*_gr
				);
			ne[j] = gr[_i][j];
		}
	}
	setError(_i);
}

void improveSlice(int i){
	double eq4;
	for (int j = 1; j < nw - 1; j++){
		eq4 = F[j][1] / dr + F[j][2] * (ne[j + 1] - 4.*ne[j] + ne[j - 1]) / (dw*dw) 
			+ F[j][3]*(ne[j + 1] - ne[j - 1]) / (2.*dw) + 2.*F[j][5]*ne[j]+F[j][6];
		gr[i][j] = ne[j] - error[i][j] / eq4;
	}
	for (int j = 1; j < nw - 1; j++){
		ne[j] = gr[i][j];
	}
	setError(i);
}

void newgr(){
	double a = 0;
	for (int i = 1; i < nr; i++){
		//setCoeff
		for (int j = 1; j < nw - 1; j++){
			for (int k = 1; k < 7; k++){
				F[j][k] = coeff(1, k, i, j);
			}
		}
		//setInit
		setInitSlice(i, 0);
		a = 0;
		while (error[i][0] > 0.1){
			std::cout << "{" << i << "," << a << "}" << std::endl;
			a++;
			std::cout << error[i][0] << std::endl;
			improveSlice(i);
		}
	}
}