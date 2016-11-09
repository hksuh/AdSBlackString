#pragma once

double dr = 0.01, dw = 3.1415926538/401;
double rmin = 0.5;
double r0 = 1.;
double l5 = 1.;

const unsigned int nr = 400;
const unsigned int nw = 401;

double gr[nr][nw];
double kr[nr][nw];
double kt[nr][nw];