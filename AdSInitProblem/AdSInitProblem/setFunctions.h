#include "variables.h"
#include <math.h>
#pragma once

void setUnpertSol(){
	double r, w;
	for (int i = 0; i < nr; i++){
		for (int j = 0; j < nw; j++){
			r = rmin + i*dr;
			w = (j + 1)*dw;
			gr[i][j] = (r*(2. - r + r*r*r + r0)) / 
				((1. + r*r*r)*(1. + r*r*r));
			kr[i][j] = ((-2. + 2.* pow(r, 4.) - 3.* r0 + r*r0 - r0*r0 - 4.* pow(r, 3.)*(1 + r0))*
				sqrt(pow(l5*(1. + pow(r, 3.))*(1./ sin(w)),2) / (r*(2. - r + pow(r, 3.) + r0))))
				/ (2.* pow((1. + pow(r, 3.)), 3.));
			kt[i][j] = -1.*((r*(-1. + r - r0)*sqrt((pow(l5*(1. + pow(r, 3.))/sin(w), 2)) / (r*(2. - r + pow(r, 3.) + r0))))
				/ (1. + pow(r, 3.)));
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
			return -1.;
		}
		if (coeffN == 2){
			//G2 = 1 / r
			return 1. / r;
		}
		if (coeffN == 3){
			//G3 = (r Krr[r, z]) / \[Gamma]rr[r, z]
			return (r*kr[i][j]) / gr[i][j];
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
			double F1 = -1.*((4.* (10.*exp(pow((5. - 2.*r), 2.)) + (1. + 2.* (5. - 2.* r)* r)*sin(2.*w))) / r);
				return F1;
		}
		if (coeffN == 2){
			double F2 = 2.*(10.*exp(pow((5. - 2.*r), 2.)) + sin(2.*w));
				return F2;
		}
		if (coeffN == 3){
			double F3 = -2.*(3. + cos(2.* w) + 30.* exp(pow((5. - 2.*r), 2.)) / tan(w));
				return F3;
		}
		if (coeffN == 4){
			double F4 = -10.*exp(pow((5. - 2.*r), 2.)) - sin(2.*w);
				return F4;
		}
		if (coeffN == 5){
			double F5 = -1.*(1. / (l5*l5*pow(r,4.)* (10.*exp(pow((5. - 2.*r), 2.)) + sin(2.*w))))*4.*(1./sin(w))*
				(-1.*10.*exp(pow((5. - 2.*r), 2.))*l5*l5*r*r*cos(w)*(-1. - 10.* r*r + (1. + 4.*r*r)* cos(2.*w))
				+ l5*l5*pow(r,4.)*(7. + 6.* cos	(2.*w))*sin(w)+ 100.*exp(2.*pow((5. - 2.*r),2.))*sin(w)*(l5*l5*r*r*(1. + 3.*r*r)
				+ kt[i][j]* kt[i][j]*sin(w)*sin(w)));
				return F5;
		}
		if (coeffN == 6){
			double F6 = (1. / (l5*l5*pow(r, 2.)* (10.*exp(pow((5. - 2.*r), 2.)) + sin(2.*w))))*4.*
				(100. * exp(2.*pow((5. - 2.*r),2.))*(l5*l5 - 2.*kr[i][j] * kt[i][j] * sin(w)*sin(w)) + 20.*exp(pow((5. - 2.*r),2.))*
				(l5*l5*(1. + 2.* r*(l5 + 4.* r*(23. + 4.*(-5. + r)* r))) - kr[i][j] * kt[i][j] * sin(w)*sin(w))*sin(2.*w)
					+ l5*l5*(1. + 4.* (-3. + r)* r*(-5. + 12.* (-2. + r)* r))*sin(2.*w)*sin(2.*w));
				return F6;
		}
	}
	else if (eqN == 2){
		if (coeffN == 1){
			double G1 = -1.;
			return G1;
		}
		if (coeffN == 2){
			double G2 = (10.* exp(pow((5. - 2.*r), 2.)) + (1. + 2.* (5. - 2.* r)* r)*sin(2.*w)) / (r*(10.*exp(pow((5. - 2.*r), 2.)) + sin(2.*w)));
			return G2;
		}
		if (coeffN == 3){
			double G3 = (exp(-1.*pow((5. - 2.*r), 2.))* r*kr[i][j] * (10.* exp(pow((5. - 2.*r), 2.)) + (1.* +10.* r - 4.*r*r)*sin(2.*w))) / (10.*gr[i][j]);
			return G3;
		}
	}
	else if (eqN == 3){
		if (coeffN == 1){
			double H1 = 2.;
			return H1;
		}
		if (coeffN == 2){
			double H2 = 2.*(1. / tan(w)) - (gr[i][j + 1] - gr[i][j - 1]) / (2.*dw) / gr[i][j];
			return H2;
		}
		if (coeffN == 3){
			double H3 = (40.*exp(pow((5. - 2.*r), 2.))*gr[i][j] * ((1. + 10 * exp(pow((5. - 2.*r), 2.))*(1 / tan(w)))*kt[i][j] +
				(10.*exp(pow((5. - 2.*r), 2.)) + sin(2.*w))*(kt[i][j + 1] - kt[i][j - 1]) / (2.*dw)))
				/ (r*r*pow((10.*exp(pow((5. - 2.*r), 2.)) + sin(2.*w)), 2.));
			return H3;
		}
	}
	return 0;
}

double constraint(int eqN, int i, int j){
	if (eqN == 1){
		return 	  coeff(1, 1, i, j)*(gr[i + 1][j] - gr[i - 1][j]) / (2.*dr)
				+ coeff(1, 2, i, j)*(gr[i][j])*(gr[i][j + 1]-2.*gr[i][j] + gr[i][j - 1]) /(dw*dw)
				+ coeff(1, 3, i, j)*gr[i][j] * (gr[i][j + 1] - gr[i][j - 1]) / (2.*dw)
				+ coeff(1, 4, i, j)*(pow((gr[i][j + 1] - gr[i][j - 1]) / (2.*dw), 2.))
				+ coeff(1, 5, i, j)*pow(gr[i][j], 2.)
				+ coeff(1, 6, i, j)*gr[i][j];
		
	}
	else if (eqN == 2){
		return  coeff(2, 1, i, j)*(kt[i + 1][j] - kt[i - 1][j]) / (2.*dr)
			  + coeff(2, 2, i, j)*(kt[i][j]) 
			  + coeff(2, 3, i, j);

	}
	else if (eqN == 3){
		return coeff(3, 1, i, j)*(kr[i][j+1] - kr[i][j-1]) / (2.*dw)
			+ coeff(3, 2, i, j)*kr[i][j]
			+ coeff(3, 3, i, j);
	}
	return -1e14;
}

double unpertConstraint(int eqN, int i, int j){
	if (eqN == 1){
		return 	  unpertCoeff(1, 1, i, j)*(gr[i + 1][j] - gr[i - 1][j]) / (2.*dr)
			+ unpertCoeff(1, 2, i, j)*(gr[i][j])*(gr[i][j + 1] - 2.*gr[i][j] + gr[i][j - 1]) / (dw*dw)
			+ unpertCoeff(1, 3, i, j)*gr[i][j] * (gr[i][j + 1] - gr[i][j - 1]) / (2.*dw)
			+ unpertCoeff(1, 4, i, j)*(pow((gr[i][j + 1] - gr[i][j - 1]) / (2.*dw), 2.))
			+ unpertCoeff(1, 5, i, j)*pow(gr[i][j], 2.)
			+ unpertCoeff(1, 6, i, j)*gr[i][j];

	}
	else if (eqN == 2){
		return  unpertCoeff(2, 1, i, j)*(kt[i + 1][j] - kt[i - 1][j]) / (2.*dr)
			+ unpertCoeff(2, 2, i, j)*(kt[i][j])
			+ unpertCoeff(2, 3, i, j);

	}
	else if (eqN == 3){
		return unpertCoeff(3, 1, i, j)*(kr[i][j + 1] - kr[i][j - 1]) / (2.*dw)
			+ unpertCoeff(3, 2, i, j)*kr[i][j]
			+ unpertCoeff(3, 3, i, j);
	}
	return -1e14;
}