#ifndef LEVMARMIRROR_H
#define LEVMARMIRROR_H
#include "levmar.h"
#include <cstddef>
#include <fstream>
#include <vector>
#include "CDCommon.h"
class LevMarMirror
{
public:
    LevMarMirror();
   int LevMarMirroring(int *features);
   static void candide1(double *p, double *x, int m, int n, void *data);
   static void candide2(double *p, double *x, int m, int n, void *data);
public:
     static int m1;
     static int m2;
     static int n;
     static int AsmPoints[];
     static int CandideVertices[];
    static  int CandideParams1[];
    static int CandideParams2[];
    static double p1[];
    static double p2[];
};

#endif // LEVMARMIRROR_H
