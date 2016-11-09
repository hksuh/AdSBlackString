#include <iostream>
#include "variables.h"
#include "setFunctions.h"
#include "improveFunctions.h"


using namespace std;

int testInitCond(){
	setUnpertSol();
	cout << "eq1(unpert) : "<< unpertConstraint(1, 120, 300) << endl;
	cout << "eq2(unpert) : " << unpertConstraint(2, 120, 16) << endl;
	cout << "eq3(unpert) : " << unpertConstraint(3, 150, 230) << endl;

	cout << "eq1 : "<< constraint(1, 120, 300) << endl;
	cout << "eq2 : " << constraint(2, 120, 16) << endl;
	cout << "eq3 : " << constraint(3, 150, 230) << endl;
	int alter = 0;
	cin >> alter;
	return alter;
}

int calculateInitialData(){
	setUnpertSol();
	double e = 0;
	newKt();
	newKr();
	for (int i = 0; i < 10; i++){
		newgr();
		newKt();
		newKr();
		e = error[1][0];
		for (int i2 = 2; i2 < nr; i2++){
			if (e < error[i2][0]){ e = error[i2][0]; }
		}
		cout << e << endl;
	}
	double alter = 0;
	cin >> alter;
	return alter;
}

int main(){
	//return calculateInitialData();
	return testInitCond();
}