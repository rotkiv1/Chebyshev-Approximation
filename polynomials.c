#include "polynomials.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct Polynomial* substractPolynomials(struct Polynomial* a, struct Polynomial* b) {
    int newDeg = a->deg > b->deg ? a->deg : b->deg;
    struct Polynomial* n = malloc(sizeof(*n));
    n->coeff = malloc((newDeg + 1) * sizeof(*n->coeff));
    
    for (int i = newDeg - abs(a->deg - b->deg); i >= 0; i--)
        n->coeff[i] = a->coeff[i] - b->coeff[i];
    
    if (a->deg > b->deg)
        for (int i = a->deg; i > a->deg - abs(a->deg - b->deg); i--)
            n->coeff[i] = a->coeff[i];
    else
        for (int i = b->deg; i > b->deg - abs(a->deg - b->deg); i--)
            n->coeff[i] = -b->coeff[i];
    
    for (int i = newDeg; i >= 0; i--)
        if (n->coeff[i] != 0.0) {
            n->deg = i;
            return n;
        }
        
    n->deg = 0;
    return n;
}

struct Polynomial* multiplyPolynomials(struct Polynomial* a, struct Polynomial* b) {
    struct Polynomial* n = malloc(sizeof(*n));
    n->deg = a->deg + b->deg;
    n->coeff = calloc((a->deg + b->deg + 1), sizeof(*n->coeff));
    
    for (int i = 0; i < a->deg + 1; i++)
        for (int j = 0; j < b->deg + 1; j++)
            n->coeff[i + j] += a->coeff[i] * b->coeff[j];
            
    return n;
}


double horner(struct Polynomial* p, double x) {
    double v = p->coeff[p->deg];
    
    for (int i = p->deg - 1; i >= 0; i--)
        v = v * x + p->coeff[i];
    
    return v;
}


struct Polynomial* ddx(struct Polynomial* p, int n) {
    if (n == 0) 
        return p;
    
    struct Polynomial* temp = malloc(sizeof(*temp));
    
    if (n > p->deg) {
        temp->deg = 0;
        temp->coeff = malloc(sizeof(*temp->coeff));
        temp->coeff[0] = 0;
        return temp;
    }
    
    p->deg -= 1;
    temp->deg = p->deg;
    temp->coeff = malloc((p->deg + 1) * sizeof(*p->coeff));
    
    for (int i = 0; i < temp->deg + 1; i++)
        temp->coeff[i] = (i + 1) * p->coeff[i + 1];
    
    ddx(temp, n - 1);
}
