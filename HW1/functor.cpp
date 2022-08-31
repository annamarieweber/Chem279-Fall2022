#include <iostream>
#include <cmath>

class quadratic {
private:
	double a;
	double b;
	double c;
public:
	quadratic(double A, double B, double C) {
		a = A;
		b = B;
		c = C;
	}
	double operator() (double x) {
		return a*std::pow(x,2) + b*x + c;	
	}
};

int main() {
	quadratic func(1,-2,1);
	for (int i = 0; i < 5; i++) {
		std::cout << "func(" << i << ") = " << func(i) << std::endl;
	}
	return 0;
}
