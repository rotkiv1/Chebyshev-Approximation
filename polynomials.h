struct Polynomial {
	int deg;
	double* coeff;
};

struct Polynomial* substractPolynomials(struct Polynomial* a, struct Polynomial* b);

struct Polynomial* multiplyPolynomials(struct Polynomial* a, struct Polynomial* b);

double horner(struct Polynomial* p, double x);

struct Polynomial* ddx(struct Polynomial* p, int n);

//void freePolynomial(struct Polynomial* p);
