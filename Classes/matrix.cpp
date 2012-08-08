#include "matrix.h"

matrix::matrix() {

}

matrix::matrix(unsigned int d) { m_.reserve(d*d); }

matrix::~matrix() {};
