#include "makespl.h"
#include "piv_ge_solver.h"
#include "polynomials.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>


void
make_spl(points_t * pts, spline_t * spl)
{

	matrix_t       *eqs= NULL;
	double         *x = pts->x;
	double         *y = pts->y;
	double		a = x[0];
	double		b = x[pts->n - 1];
	int		i, j, k;
	int		nb = pts->n - 3 > 10 ? 10 : pts->n - 3;
  char *nbEnv= getenv( "APPROX_BASE_SIZE" );

	if( nbEnv != NULL && atoi( nbEnv ) > 0 )
		nb = atoi( nbEnv );

	eqs = make_matrix(nb, nb + 1);

	// wielomian postaci 2x
	struct Polynomial* p = malloc(sizeof(*p));
	p->deg = 1;
        p->coeff=malloc((p->deg + 1) * sizeof(*p->coeff));
	p->coeff[0] = 0;
	p->coeff[1] = 2;
	
	// tablica zawierajaca nb wielomianow Czebyszewa
	struct Polynomial* c[nb];
	
	// wielomian postaci 1
	c[0] = malloc(sizeof(*c[0]));
	c[0]->deg = 0;
	c[0]->coeff = malloc(sizeof(*c[0]->coeff));
	c[0]->coeff[0] = 1;

	// wielomian postaci x
	c[1]=malloc(sizeof(*c[1]));
	c[1]->deg = 1;
	c[1]->coeff = malloc((c[1]->deg + 1) * sizeof(*c[i]->coeff));
	c[1]->coeff[0] = 0;
	c[1]->coeff[1] = 1;

	// generowanie nb kolejnych wielomianow Czebyszewa
	for (int i = 2; i < nb; i++) {
		c[i] = malloc(sizeof(*c[i]));
		c[i] = multiplyPolynomials(p, c[i - 1]);
		c[i] = substractPolynomials(c[i], c[i - 2]);
	}

	for (j = 0; j < nb; j++) {
		for (i = 0; i < nb; i++)
			for (k = 0; k < pts->n; k++)
		             add_to_entry_matrix(eqs, j, i, horner(c[i], x[k]) * horner(c[j], x[k]));

		for (k = 0; k < pts->n; k++)
                             add_to_entry_matrix(eqs, j, nb, y[k] * horner(c[j], x[k]));
	}

#ifdef DEBUG
	write_matrix(eqs, stdout);
#endif

	if (piv_ge_solver(eqs)) {
		spl->n = 0;
		return;
	}
#ifdef DEBUG
	write_matrix(eqs, stdout);
#endif

	if (alloc_spl(spl, nb) == 0) {
		for (int z = 0; z < spl->n; z++) {
			double xx = spl->x[z] = a + z*(b-a)/(spl->n-1);
			xx+= 10.0*DBL_EPSILON;  // zabezpieczenie przed ulokowaniem punktu w poprzednim przedziale
		
			spl->f[z] = 0;
			spl->f1[z] = 0;
			spl->f2[z] = 0;
			spl->f3[z] = 0;
			spl->f4[z] = 0;
			for (int i = 0; i < nb; i++) {
				double ck = get_entry_matrix(eqs, i, nb);
				for (int j = 0; j < 5; j++) {
					struct Polynomial* temp = malloc(sizeof(*temp));
					temp->deg = c[i]->deg;
					temp->coeff=malloc((c[i]->deg + 1) * sizeof(*temp->coeff));
					for (int j = 0; j < temp->deg + 1; j++) 
						temp->coeff[j] = c[i]->coeff[j];					                				       if (j == 0) 
						spl->f[z] += ck * horner(ddx(c[i], j), xx);					
					if (j == 1) 
						spl->f1[z] += ck * horner(ddx(c[i], j), xx);   
					if (j == 2) 
						spl->f2[z] += ck * horner(ddx(c[i], j), xx);   
					if (j == 3) 
						spl->f3[z] += ck * horner(ddx(c[i], j), xx);   
					if (j == 4) 
						spl->f4[z] += ck * horner(ddx(c[i], j), xx);
					free(c[i]);
					c[i] = malloc(sizeof(*c[i]));
					c[i]->deg = temp->deg;
					c[i]->coeff = malloc((c[i]->deg + 1) * sizeof(*c[i]->coeff));
				        for (int k = 0; k < c[i]->deg+1; k++) 
			                        c[i]->coeff[k] = temp->coeff[k];
		           	 }				
			}
		}
	}

	//zwalnianie pamieci wielomianow
//	for (int i = 0; i < nb; i++) 
//		freePolynomial(c[i]);

	//zwalnianie pamieci macierzy
//	freeMatrix(eqs);


}
