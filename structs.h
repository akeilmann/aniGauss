#ifndef STRUCTS_H
#define STRUCTS_H

struct coefficients{
    double b[5];
    double M[9];
};

struct filterLines{
    double* line_n;
    double* line_n1;
    double* line_n2;
    double* line_n3;
};


#endif // STRUCTS_H
