#include "../ComplexNum.h"
#include "../LA.h"
#pragma once

class kernelFunctions {
 static ComplexNum Linear(Matrix x, Matrix y);
 static ComplexNum Polynomial(Matrix x, Matrix y, ComplexNum constant, double degree);
 static ComplexNum Gaussian(Matrix x, Matrix y, ComplexNum gamma);
};
