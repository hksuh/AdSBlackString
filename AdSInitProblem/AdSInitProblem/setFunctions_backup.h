#include "variables.h"
#include <math.h>
#pragma once

void setUnpertSol(){
	double r, w;
	for (int i = 0; i < nr; i++){
		for (int j = 0; j < nw; j++){
			//r=rmin+i*dr
			//w=(j+1)*dw
			//gr[i][j]=
			//kr[i][j]=
			//kt[i][j]=
		}
	}
}

double unpertCoeff(int eqN, int coeffN, int i, int j){
	double r = rmin + i*dr;
	double w = (j+1)*dw;
	if (eqN == 1){
		if (coeffN == 1){
			//F1 = 4 / r
			return 4 / r;
		}
		if (coeffN == 2){
			//F2 = -2
			return -2;
		}
		if (coeffN == 3){
			//F3 = 6 Cot[z]
			return 6. / tan(w);
		}
		if (coeffN == 4){
			//F4 = 1
			return 1.;
		}
		if (coeffN == 5){
			//F5 = 4 (3 + 1 / r ^ 2 + (K\[Theta]\[Theta][r, z] ^ 2 Sin[z] ^ 2) / (l5 ^ 2 r ^ 4))
			double a = (kt[i][j] * sin(w)) / (l5 * r *r);
			return 4.*(3. + 1. / (r*r) + a*a);
		}
		if (coeffN == 6){
			//F6 = -((4 (l5 ^ 2 - 2 Krr[r, z] K\[Theta]\[Theta][r, z] Sin[z] ^ 2)) / (l5 ^ 2 r ^ 2))
			return -1.*((4. *(l5 *l5 - 2.* kr[i][j] * kt[i][j] * sin(w)* sin(w))) / (l5 *l5 *r * r));
		}
	}
	else if (eqN == 2){
		if (coeffN == 1){
			//G1 = -1
			return -1;
		}
		if (coeffN == 2){
			//G2 = 1 / r
			return 1. / r;
		}
		if (coeffN == 3){
			//G3 = (r Krr[r, z]) / \[Gamma]rr[r, z]
			return (r*kt[i][j]) / gr[i][j];
		}
	}
	else if (eqN == 3){
		if (coeffN == 1){
			return 2.;
		}
		if (coeffN == 2){
			return 2. / tan(w) - (gr[i][j + 1] - gr[i][j - 1]) / (2.*dw*gr[i][j]);
		}
		if (coeffN == 3){
			return (4. *gr[i][j] * (kt[i][j] / tan(w) + (kt[i][j + 1] - kt[i][j - 1]) / (2.*dw))) / (r *r);
		}
	}
	return 0;
}

double coeff(int eqN, int coeffN, int i, int j){
	double r = rmin + i*dr;
	double w = (j+1)*dw;
	if (eqN == 1){
		if (coeffN == 1){
			//return F1	
		}
		if (coeffN == 2){
			//return F2	
		}
		if (coeffN == 3){
			//return F3	
		}
		if (coeffN == 4){
			//return F4	
		}
		if (coeffN == 5){
			//return F5	
		}
		if (coeffN == 6){
			//return F6	
		}
	}
	else if (eqN == 2){
		if (coeffN == 1){
			//return G1	
		}
		if (coeffN == 2){
			//return G2	
		}
		if (coeffN == 3){
			//return G3	
		}
	}
	else if (eqN == 3){
		if (coeffN == 1){
			//return H1	
		}
		if (coeffN == 2){
			//return H2	
		}
		if (coeffN == 3){
			//return H3	
		}
	}
	return 0;
}

double constraint(int eqN, int i, int j){
	if (eqN == 1){
		//return eq1
	}
	else if (eqN == 2){
		//return eq2
	}
	else if (eqN == 3){
		//eq3 = H1 Krr ^ (1, 0)[r, z] + H2 Krr[r, z] + H3
		return unpertCoeff(3, 1, i, j)*(kr[i + 1][j] - kr[i - 1][j]) / (2.*dr)
			+ unpertCoeff(3, 2, i, j)*kr[i][j]
			+ unpertCoeff(3, 3, i, j);
	}
	return -1e14;
}