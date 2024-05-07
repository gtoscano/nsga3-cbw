/* Test problem definitions */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>

#include "global.h"
#include "rand.h"
#include <misc.hpp>

/* # define sch1 */
/* # define sch2 */
/* # define fon */
/* # define kur */
/* # define pol */
/* # define vnt */
/* # define zdt1 */
/* # define zdt2 */
/* # define zdt3 */
/* # define zdt4 */
/* # define zdt5 */
/* # define zdt6 */
/* # define bnh */
/* # define osy */
/* # define srn */
/* # define tnk */
/* # define ctp1 */
/* # define ctp2 */
/* # define ctp3 */
/* # define ctp4 */
/* # define ctp5 */
/* # define ctp6 */
/* # define ctp7 */
/* #define dtlz1 */
#define cbw2

/************************************************************ MANY OBJECTIVE PROBLEMS ************************************************************/

#ifdef cbw2
int test_problem2(double *xreal, [[maybe_unused]] double *xbin, int **gene, double *obj, double *constr, std::string emo_uuid, std::string exec_uuid)
{
    std::vector<double> x(nreal, 0.0);
    std::vector<double> fx(nobj, 0.0);
    std::vector<double> g(ncon, 0.0);
    std::vector<double> g2(ncon, 0.0);

    std::copy(xreal, xreal+nreal, x.begin());
    int nrows = eval_f( x, fx, g, emo_uuid, exec_uuid);
    std::copy(fx.begin(), fx.end(), obj);
    //std::copy(g2.begin(), g2.end(), constr);
    std::copy(g.begin(), g.end(), constr);
    std::copy(x.begin(), x.end(), xreal);
    //constr[0]=g[0];
    return nrows;
}
#endif




/*********************** Test Problems ***************************/

/* Unconstrained */

#ifdef dtlz1
void test_problem(double *xreal, __attribute__((unused)) double *xbin, __attribute__((unused)) int **gene, double *obj, __attribute__((unused)) double *constr)
{
    double g, prod;
    int j, k;

    for (j = nobj - 1, g = 0.0; j < nreal; j++)
    {
        g += 100.0 * (1.0 + (xreal[j] - 0.5) * (xreal[j] - 0.5) - cos(20.0 * PI * (xreal[j] - 0.5)));
    }

    for (j = (nobj - 1); j >= 0; j--)
    {
        prod = 1.0;
        for (k = 0; k < nobj - j - 1; k++)
        {
            prod *= (k >= 0) ? xreal[k] : 1.0;
        }
        obj[j] = 0.5 * (1.0 + g) * prod * ((j > 0) ? (1.0 - xreal[nobj - j - 1]) : 1.0);
    }
}
#endif

#ifdef dtlz1_inv
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    double g, prod;
    int j, k;

    for (j = nobj - 1, g = 0.0; j < nreal; j++)
    {
        g += 100.0 * (1.0 + (xreal[j] - 0.5) * (xreal[j] - 0.5) - cos(20.0 * PI * (xreal[j] - 0.5)));
    }

    for (j = (nobj - 1); j >= 0; j--)
    {
        prod = 1.0;
        for (k = 0; k < nobj - j - 1; k++)
        {
            prod *= (k >= 0) ? xreal[k] : 1.0;
        }
        obj[j] = 0.5 * (1.0 + g) * prod * ((j > 0) ? (1.0 - xreal[nobj - j - 1]) : 1.0);
        obj[j] = 0.5 * (1 + g) - obj[j];
    }
}
#endif

#ifdef dtlz2_constr
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{

    double g, prod;
    int j, k;

    for (j = nobj - 1, g = 0.0; j < nreal; j++)
    {
        g += (xreal[j] - 0.5) * (xreal[j] - 0.5);
    }
    for (j = (nobj - 1); j >= 0; j--)
    {
        prod = 1.0;
        for (k = 0; k < nobj - j - 1; k++)
        {
            prod *= (k >= 0) ? cos(xreal[k] * PI / 2.0) : 1.0;
        }
        if (j == 2)
        {
            obj[j] = (1.0 + g) *
                     prod * ((j > 0) ? sin(xreal[nobj - j - 1] * PI / 2.0) : 1.0);
        }
        else
        {
            obj[j] = (1.0 + g) *
                     prod * ((j > 0) ? sin(xreal[nobj - j - 1] * PI / 2.0) : 1.0);
        }
    }
    constr[0] = (0.15 - obj[nobj - 1]) * (0.85 - obj[nobj - 1]);
}
#endif

#ifdef dtlz2
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{

    double g, prod;
    int j, k;

    for (j = nobj - 1, g = 0.0; j < nreal; j++)
    {
        g += (xreal[j] - 0.5) * (xreal[j] - 0.5);
    }
    for (j = (nobj - 1); j >= 0; j--)
    {
        prod = 1.0;
        for (k = 0; k < nobj - j - 1; k++)
        {
            prod *= (k >= 0) ? cos(xreal[k] * PI / 2.0) : 1.0;
        }
        if (j == 2)
        {
            obj[j] = (1.0 + g) *
                     prod * ((j > 0) ? sin(xreal[nobj - j - 1] * PI / 2.0) : 1.0);
        }
        else
        {
            obj[j] = (1.0 + g) *
                     prod * ((j > 0) ? sin(xreal[nobj - j - 1] * PI / 2.0) : 1.0);
        }
    }
}
#endif

#ifdef dtlz2_inv
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{

    double g, prod;
    int j, k;

    for (j = nobj - 1, g = 0.0; j < nreal; j++)
    {
        g += (xreal[j] - 0.5) * (xreal[j] - 0.5);
    }
    for (j = (nobj - 1); j >= 0; j--)
    {
        prod = 1.0;
        for (k = 0; k < nobj - j - 1; k++)
        {
            prod *= (k >= 0) ? cos(xreal[k] * PI / 2.0) : 1.0;
        }
        if (j == 2)
        {
            obj[j] = (1.0 + g) *
                     prod * ((j > 0) ? sin(xreal[nobj - j - 1] * PI / 2.0) : 1.0);
        }
        else
        {
            obj[j] = (1.0 + g) *
                     prod * ((j > 0) ? sin(xreal[nobj - j - 1] * PI / 2.0) : 1.0);
        }
    }
    for (j = 0; j < nobj - 1; j++)
        obj[j] = pow((1.0 + g), 4.0) - pow(obj[j], 4.0);
    obj[nobj - 1] = pow((1.0 + g), 4.0) - pow(obj[nobj - 1], 2.0);
}
#endif

#ifdef dtlz2_convex
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{

    double g, prod;
    int j, k;

    for (j = nobj - 1, g = 0.0; j < nreal; j++)
    {
        g += (xreal[j] - 0.5) * (xreal[j] - 0.5);
    }
    for (j = (nobj - 1); j >= 0; j--)
    {
        prod = 1.0;
        for (k = 0; k < nobj - j - 1; k++)
        {
            prod *= (k >= 0) ? cos(xreal[k] * PI / 2.0) : 1.0;
        }
        if (j == 2)
        {
            obj[j] = (1.0 + g) *
                     prod * ((j > 0) ? sin(xreal[nobj - j - 1] * PI / 2.0) : 1.0);
        }
        else
        {
            obj[j] = (1.0 + g) *
                     prod * ((j > 0) ? sin(xreal[nobj - j - 1] * PI / 2.0) : 1.0);
        }
    }
    for (j = 0; j < nobj - 1; j++)
    {
        obj[j] = pow(obj[j], 4.0);
    }
    obj[nobj - 1] = obj[nobj - 1] * obj[nobj - 1];
}
#endif

#ifdef dtlz3
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{

    double g, prod;
    int j, k;

    for (j = nobj - 1, g = 0.0; j < nreal; j++)
    {
        g += 100.0 * (1.0 + (xreal[j] - 0.5) * (xreal[j] - 0.5) - cos(20.0 * PI * (xreal[j] - 0.5)));
    }
    for (j = (nobj - 1); j >= 0; j--)
    {
        prod = 1.0;
        for (k = 0; k < nobj - j - 1; k++)
        {
            prod *= (k >= 0) ? cos(xreal[k] * PI / 2.0) : 1.0;
        }
        obj[j] = (1.0 + g) *
                 prod * ((j > 0) ? sin(xreal[nobj - j - 1] * PI / 2.0) : 1.0);
    }
}
#endif

#ifdef dtlz4
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    double prod, g;
    int j, k;
    g = 0.0;
    for (j = nobj - 1, g = 0.0; j < nreal; j++)
    {
        g += (xreal[j] - 0.5) * (xreal[j] - 0.5);
    }

    for (j = (nobj - 1); j >= 0; j--)
    {
        prod = 1.0;
        for (k = 0; k < nobj - j - 1; k++)
        {
            prod *= (k >= 0) ? cos(0.5 * (1.0 + pow((xreal[k] - 0.5) / 0.5, 99.0)) * PI / 2.0) : 1.0;
        }
        obj[j] = (1.0 + g) *
                 prod * ((j > 0) ? sin(0.5 * (1.0 + pow((xreal[nobj - j - 1] - 0.5) / 0.5, 99.0)) * PI / 2.0) : 1.0);
    }
}
#endif

#ifdef dtlz5
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    double prod, g;
    int j, k;
    g = 0.0;
    for (j = nobj - 1, g = 0.0; j < nreal; j++)
    {
        g += pow(xreal[j], 0.3);
    }

    for (j = (nobj - 1); j >= 0; j--)
    {
        prod = 1.0;
        for (k = 0; k < nobj - j - 1; k++)
        {
            prod *= (k >= 0) ? ((k == 0) ? cos(xreal[k] * PI / 2.0) : cos((2.0 * g * xreal[k] + 1.0) * PI / (4.0 * (1.0 + g))))
                             : 1.0;
            /* prod *= (k>=0) ? ((k==0) ? cos(xreal[k]*PI/2.0) :
                                cos(xreal[k]*acos(cos(PI/50)/(1+g))))
                  : 1.0;*/
        }

        obj[j] = (1.0 + g) *
                 prod * ((j > 0) ? ((j == nobj - 1) ? sin(xreal[nobj - j - 1] * PI / 2.0) : sin((2.0 * g * xreal[k] + 1.0) * PI / (4.0 * (1.0 + g)))) : 1.0);
        /*obj[j] = (1.0+g) *
              prod * ((j>0) ? ((j==nobj-1) ?
                               sin(xreal[nobj-j-1]*PI/2.0) :
                               sin((1-xreal[k])*sin(PI/50)/(1+g)))
                      : 1.0);*/
    }
}
#endif

#ifdef dtlz6
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    double sm, g;
    int j;
    g = 0.0;
    for (j = nobj - 1, g = 0.0; j < nreal; j++)
    {
        g += xreal[j];
    }
    g = 1.0 + 9.0 * g / (nreal - nobj + 1.0);

    for (j = 0; j < nobj - 1; j++)
    {
        obj[j] = xreal[j];
    }

    for (j = 0, sm = 0.0; j < nobj - 1; j++)
    {
        sm += obj[j] * (1.0 + sin(3.0 * PI * obj[j])) / (1 + g);
    }
    obj[nobj - 1] = (1 + g) * (nobj - sm);
}
#endif

/* Constrained */
#ifdef c1_dtlz1
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    double g, prod;
    int j, k;

    for (j = nobj - 1, g = 0.0; j < nreal; j++)
    {
        g += 100.0 * (1.0 + (xreal[j] - 0.5) * (xreal[j] - 0.5) - cos(20.0 * PI * (xreal[j] - 0.5)));
    }

    for (j = (nobj - 1); j >= 0; j--)
    {
        prod = 1.0;
        for (k = 0; k < nobj - j - 1; k++)
        {
            prod *= (k >= 0) ? xreal[k] : 1.0;
        }
        obj[j] = 0.5 * (1.0 + g) * prod * ((j > 0) ? (1.0 - xreal[nobj - j - 1]) : 1.0);
    }
    constr[0] = 0.0;
    for (j = 0; j < nobj - 1; j++)
        constr[0] += obj[j] / 0.5;
    constr[0] += obj[nobj - 1] / 0.6;
    constr[0] = 1 - constr[0];
}
#endif

#ifdef c1_dtlz3
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{

    double g, prod, r;
    int j, k;

    if (nobj < 5)
        r = 9.0;
    else if (nobj <= 12 && nobj >= 5)
        r = 12.5;
    else
        r = 15.0;

    for (j = nobj - 1, g = 0.0; j < nreal; j++)
    {
        g += 10.0 * (1.0 + (xreal[j] - 0.5) * (xreal[j] - 0.5) - cos(20.0 * PI * (xreal[j] - 0.5)));
    }
    for (j = (nobj - 1); j >= 0; j--)
    {
        prod = 1.0;
        for (k = 0; k < nobj - j - 1; k++)
        {
            prod *= (k >= 0) ? cos(xreal[k] * PI / 2.0) : 1.0;
        }
        obj[j] = (1.0 + g) *
                 prod * ((j > 0) ? sin(xreal[nobj - j - 1] * PI / 2.0) : 1.0);
    }
    prod = 0.0;
    for (j = 0; j < nobj; j++)
    {
        prod += obj[j] * obj[j];
    }
    constr[0] = (prod - (4.0 * 4.0)) * (prod - (r * r));
}
#endif

#ifdef c2_dtlz2
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{

    double g, prod, r;
    int j, k;

    if (nobj < 5)
        r = 0.5;
    else
        r = 0.5;

    for (j = nobj - 1, g = 0.0; j < nreal; j++)
    {
        g += (xreal[j] - 0.5) * (xreal[j] - 0.5);
    }
    for (j = (nobj - 1); j >= 0; j--)
    {
        prod = 1.0;
        for (k = 0; k < nobj - j - 1; k++)
        {
            prod *= (k >= 0) ? cos(xreal[k] * PI / 2.0) : 1.0;
        }
        if (j == 2)
        {
            obj[j] = (1.0 + g) *
                     prod * ((j > 0) ? sin(xreal[nobj - j - 1] * PI / 2.0) : 1.0);
        }
        else
        {
            obj[j] = (1.0 + g) *
                     prod * ((j > 0) ? sin(xreal[nobj - j - 1] * PI / 2.0) : 1.0);
        }
    }
    constr[0] = -1 * INF;
    for (j = 0; j < nobj; j++)
    {
        g = 0.0;
        for (k = 0; k < nobj; k++)
        {
            if (k == j)
                g += (obj[k] - 1.0) * (obj[k] - 1.0);
            else
                g += obj[k] * obj[k];
        }
        g = r * r - g;
        if (g > constr[0])
            constr[0] = g;
    }
    /*g=0.0;
	for(k=0;k<nobj;k++)	
		g += (obj[k]-1.0/sqrt((double)nobj))*(obj[k]-1.0/sqrt((double)nobj));
	
	g = r*r - g;
	if(g>constr[0])
		constr[0]=g;*/
}
#endif

#ifdef dtlz1_hole
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    double g, prod, lambda, r;
    int j, k;

    if (nobj < 4)
        r = 0.2;
    else if (nobj > 3 && nobj < 6)
        r = 0.2;
    else
        r = 0.225;

    for (j = nobj - 1, g = 0.0; j < nreal; j++)
    {
        g += 100.0 * (1.0 + (xreal[j] - 0.5) * (xreal[j] - 0.5) - cos(20.0 * PI * (xreal[j] - 0.5)));
    }

    for (j = (nobj - 1); j >= 0; j--)
    {
        prod = 1.0;
        for (k = 0; k < nobj - j - 1; k++)
        {
            prod *= (k >= 0) ? xreal[k] : 1.0;
        }
        obj[j] = 0.5 * (1.0 + g) * prod * ((j > 0) ? (1.0 - xreal[nobj - j - 1]) : 1.0);
    }

    lambda = 0.0;
    for (j = 0; j < nobj; j++)
    {
        lambda = lambda + obj[j];
    }
    lambda = lambda / ((double)nobj);
    constr[0] = 0.0;
    for (j = 0; j < nobj; j++)
        constr[0] = constr[0] + pow((obj[j] - lambda), 2.0);

    constr[0] = r * r - constr[0];
}
#endif

#ifdef c3_dtlz1
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    double g, prod;
    int i, j, k;

    for (j = nobj - 1, g = 0.0; j < nreal; j++)
    {
        g += 100.0 * (1.0 + (xreal[j] - 0.5) * (xreal[j] - 0.5) - cos(20.0 * PI * (xreal[j] - 0.5)));
    }

    for (j = (nobj - 1); j >= 0; j--)
    {
        prod = 1.0;
        for (k = 0; k < nobj - j - 1; k++)
        {
            prod *= (k >= 0) ? xreal[k] : 1.0;
        }
        obj[j] = 0.5 * (1.0 + g) * prod * ((j > 0) ? (1.0 - xreal[nobj - j - 1]) : 1.0);
    }
    for (i = 0; i < nobj; i++)
    {
        constr[i] = 0.0;
        for (j = 0; j < nobj; j++)
        {
            if (j == i)
                constr[i] += obj[j];
            else
                constr[i] += obj[j] / 0.5;
        }
        constr[i] = constr[i] - 1.0;
    }
}
#endif

#ifdef c3_dtlz4
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    double prod, g;
    int j, k;
    g = 0.0;
    for (j = nobj - 1, g = 0.0; j < nreal; j++)
    {
        g += 10 * (xreal[j] - 0.5) * (xreal[j] - 0.5);
    }

    for (j = (nobj - 1); j >= 0; j--)
    {
        prod = 1.0;
        for (k = 0; k < nobj - j - 1; k++)
        {
            prod *= (k >= 0) ? cos(0.5 * (1.0 + pow((xreal[k] - 0.5) / 0.5, 99.0)) * PI / 2.0) : 1.0;
        }
        obj[j] = (1.0 + g) *
                 prod * ((j > 0) ? sin(0.5 * (1.0 + pow((xreal[nobj - j - 1] - 0.5) / 0.5, 99.0)) * PI / 2.0) : 1.0);
    }
    for (j = 0; j < nobj; j++)
    {
        constr[j] = 0.0;
        for (k = 0; k < nobj; k++)
        {
            if (k == j)
                constr[j] += obj[k] * obj[k] / 4.0;
            else
                constr[j] += obj[k] * obj[k];
        }
        constr[j] = constr[j] - 1.0;
    }
}
#endif

/*********************** Real World Problems ***************************/

/* Unconstrained */

#ifdef crash
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    obj[0] = 1640.2823 + 2.3573285 * xreal[0] + 2.3220035 * xreal[1] +
             4.5688768 * xreal[2] + 7.7213633 * xreal[3] + 4.4559504 * xreal[4];
    obj[1] = 6.5856 + 1.15 * xreal[0] - 1.0427 * xreal[1] + 0.9738 * xreal[2] +
             +0.8364 * xreal[3] - 0.3695 * xreal[0] * xreal[3] + 0.0861 * xreal[0] * xreal[4] +
             0.3628 * xreal[1] * xreal[3] - 0.1106 * xreal[0] * xreal[0] - 0.3437 * xreal[2] * xreal[2] + 0.1764 * xreal[3] * xreal[3];
    obj[2] = -0.0551 + 0.0181 * xreal[0] + 0.1024 * xreal[1] + 0.0421 * xreal[2] - 0.0073 * xreal[0] * xreal[1] + 0.0204 * xreal[1] * xreal[2] - 0.0118 * xreal[1] * xreal[3] - 0.0204 * xreal[2] * xreal[3] - 0.008 * xreal[2] * xreal[4] - 0.0241 * xreal[1] * xreal[1] + 0.0109 * xreal[3] * xreal[3];
}
#endif

/* Constrained */

#ifdef welded_beam
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    double t, t1, t2, s, p;

    /*
	h=x0
	l=x1
	t=x2
	b=x3
	*/
    obj[0] = 1.10471 * xreal[0] * xreal[0] * xreal[1] + 0.04811 * xreal[2] * xreal[3] * (14.0 + xreal[1]);
    obj[1] = 2.1952 / (xreal[2] * xreal[2] * xreal[2] * xreal[3]);
    t1 = 6000.0 / (sqrt(2) * xreal[0] * xreal[1]);
    t2 = (6000.0 * (14.0 + 0.5 * xreal[1]) * sqrt(xreal[1] * xreal[1] + (xreal[0] + xreal[2]) * (xreal[0] + xreal[2]))) /
         (2.0 * 0.707 * xreal[0] * xreal[1] * (xreal[1] * xreal[1] / 12.0 + 0.25 * (xreal[0] + xreal[2]) * (xreal[0] + xreal[2])));
    t = sqrt(t1 * t1 + t2 * t2 + (xreal[1] * t1 * t2) / (sqrt(xreal[1] * xreal[1] + (xreal[0] + xreal[2]) * (xreal[0] + xreal[2]))));
    s = 504000.0 / (xreal[2] * xreal[2] * xreal[3]);
    p = 64746.022 * (1 - 0.0282346 * xreal[2]) * xreal[2] * xreal[3] * xreal[3] * xreal[3];
    constr[0] = 13600.0 - t;
    constr[0] = constr[0] / 13600.0;
    constr[1] = 30000.0 - s;
    constr[1] = constr[1] / 30000.0;
    constr[2] = xreal[3] - xreal[0];
    constr[3] = (p - 6000.0) / 6000;
}
#endif

#ifdef car
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{

    double yy[2];
    yy[0] = 0.345;
    yy[1] = 0.192;
    constr[0] = 1.16 - 0.3717 * xreal[1] * xreal[3] - 0.484 * xreal[2] * yy[1];
    constr[1] = 0.261 - 0.0159 * xreal[0] * xreal[1] - 0.188 * xreal[0] * yy[0] - 0.019 * xreal[1] * xreal[6] + 0.0144 * xreal[2] * xreal[4] + 0.08045 * xreal[5] * yy[1];
    constr[2] = 0.214 + 0.00817 * xreal[4] - 0.131 * xreal[0] * yy[0] - 0.0704 * xreal[0] * yy[1] + 0.03099 * xreal[1] * xreal[5] - 0.018 * xreal[1] * xreal[6] +
                0.0208 * xreal[2] * yy[0] + 0.121 * xreal[2] * yy[1] - 0.00364 * xreal[4] * xreal[5] - 0.018 * xreal[1] * xreal[1];
    constr[3] = 0.74 - 0.61 * xreal[1] - 0.163 * xreal[2] * yy[0] - 0.166 * xreal[6] * yy[1] + 0.227 * xreal[1] * xreal[1];
    constr[4] = 28.98 + 3.818 * xreal[2] - 4.2 * xreal[0] * xreal[1] + 6.63 * xreal[5] * yy[1] - 7.77 * xreal[6] * yy[0];
    constr[5] = 33.86 + 2.95 * xreal[2] - 5.057 * xreal[0] * xreal[1] - 11 * xreal[1] * yy[0] - 9.98 * xreal[6] * yy[0] + 22 * yy[0] * yy[1];
    constr[6] = 46.36 - 9.9 * xreal[1] - 12.9 * xreal[0] * yy[0];
    constr[7] = 4.72 - 0.5 * xreal[3] - 0.19 * xreal[1] * xreal[2];
    constr[8] = 10.58 - 0.674 * xreal[0] * xreal[1] - 1.95 * xreal[1] * yy[0];
    constr[9] = 16.45 - 0.489 * xreal[2] * xreal[6] - 0.843 * xreal[4] * xreal[5];
    obj[0] = 1.98 + 4.9 * xreal[0] + 6.67 * xreal[1] + 6.98 * xreal[2] + 4.01 * xreal[3] + 1.78 * xreal[4] + 0.00001 * xreal[5] + 2.73 * xreal[6];
    obj[1] = constr[7];
    obj[2] = (constr[8] + constr[9]) / 2.0;
    constr[0] = 1 - constr[0] / 1.0;
    constr[1] = 1 - constr[1] / 0.32;
    constr[2] = 1 - constr[2] / 0.32;
    constr[3] = 1 - constr[3] / 0.32;
    constr[4] = 1 - constr[4] / 32.0;
    constr[5] = 1 - constr[5] / 32.0;
    constr[6] = 1 - constr[6] / 32.0;
    constr[7] = 1 - constr[7] / 4.0;
    constr[8] = 1 - constr[8] / 9.9;
    constr[9] = 1 - constr[9] / 15.7;
    return;
}
#endif
#ifdef water
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    obj[0] = 106780.37 * (xreal[1] + xreal[2]) + 61704.67;
    obj[0] = obj[0] / (8.0 * 10000.0);
    obj[1] = 3000.0 * xreal[0];
    obj[1] = obj[1] / 1500.0;
    obj[2] = 305700.0 * 2289.0 * xreal[1] / pow(0.06 * 2289.0, 0.65);
    obj[2] = obj[2] / (3.0 * 1000000.0);
    obj[3] = 250.0 * 2289.0 * exp(-39.75 * xreal[1] + 9.9 * xreal[2] + 2.74);
    obj[3] = obj[3] / (6.0 * 1000000.0);
    obj[4] = 25.0 * ((1.39 / (xreal[0] * xreal[1])) + 4940.0 * xreal[2] - 80.0);
    obj[4] = obj[4] / 8000.0;
    constr[0] = 0.00139 / (xreal[0] * xreal[1]) + 4.94 * xreal[2] - 0.08;
    constr[0] = 1.0 - constr[0];
    constr[1] = 0.000306 / (xreal[0] * xreal[1]) + 1.082 * xreal[2] - 0.0986;
    constr[1] = 1.0 - constr[1];
    constr[2] = 12.307 / (xreal[0] * xreal[1]) + 49408.24 * xreal[2] + 4051.02;
    constr[2] = 50000.0 - constr[2];
    constr[2] = constr[2] / 50000.0;
    constr[3] = 2.098 / (xreal[0] * xreal[1]) + 8046.33 * xreal[2] - 696.71;
    constr[3] = 16000.0 - constr[3];
    constr[3] = constr[3] / 16000.0;
    constr[4] = 2.138 / (xreal[0] * xreal[1]) + 7883.39 * xreal[2] - 705.04;
    constr[4] = 10000.0 - constr[4];
    constr[4] = constr[4] / 10000.0;
    constr[5] = 0.417 * (xreal[0] * xreal[1]) + 1721.26 * xreal[2] - 136.54;
    constr[5] = 2000.0 - constr[5];
    constr[5] = constr[5] / 2000.0;
    constr[6] = 0.164 / (xreal[0] * xreal[1]) + 631.13 * xreal[2] - 54.48;
    constr[6] = 550.0 - constr[6];
    constr[6] = constr[6] / 550.0;
}
#endif

#ifdef machining
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    obj[0] = 7.49 - 0.44 * xreal[0] + 1.16 * xreal[1] - 0.61 * xreal[2];
    obj[1] = 4.13 - 0.92 * xreal[0] + 0.16 * xreal[1] - 0.43 * xreal[2];
    obj[2] = -21.90 + 1.94 * xreal[0] + 0.30 * xreal[1] + 1.04 * xreal[2];
    obj[3] = 11.331 - xreal[0] - xreal[1] - xreal[2];
    constr[0] = 0.44 * xreal[0] - 1.16 * xreal[1] + 0.61 * xreal[2] - 3.1725;
    constr[1] = 0.92 * xreal[0] - 0.16 * xreal[1] + 0.43 * xreal[2] - 8.0420;
    constr[2] = -1.94 * xreal[0] + 0.30 * xreal[1] + 1.04 * xreal[2] + 18.4988;
}
#endif

#ifdef wiper
void test_problem(double *xreal, double *xbin, int **gene, double *obj,
                  double *constr)
{
    double g1leq, g2leq, g3leq, g4leq, g5leq, v, f, a, t1, n1, t2, t3, t4;
    double f1, f2, f3;

    v = xreal[0];
    f = xreal[1];
    a = xreal[2];
    t1 = ((1.25 * 3.1419) / (v * f));
    n1 = pow(v, 1.672184) * pow(f, 0.036654) * pow(a, 0.072133);

    t2 = pow(993.402604 * pow(v, -0.200600) * pow(f, 0.623620) * pow(a, 0.660314) * pow(1800.0, 0.158012), 2.0);
    t3 = pow(320.402695 * pow(v, -0.213734) * pow(f, 0.293166) * pow(a, 0.373255) * pow(1800.0, 0.240112), 2.0);
    t4 = pow(814.912603 * pow(v, -0.388726) * pow(f, 0.386365) * pow(a, 1.502916) * pow(1800.0, 0.285235), 2.0);
    f1 = t1 * (1 + ((n1) / (85219.833))) + 1.6;
    f2 = t1 * (0.08 + ((n1) / (1693.304))) + 0.14;
    f3 = sqrt(t1 + t2 + t3);
    g1leq = 0.634065 - 0.004599 * v - 1.221571 * f - 1.125925 * a + 0.000010 * v * v + 2.265811 * f * f + 0.007205 * v * f + 0.007658 * v * a + 0.0000001 * v * 1800.0 + 6.020526 * f * a - 0.039339 * v * f * a - 0.8;
    g2leq = 0.2 - (0.634065 - 0.004599 * v - 1.221571 * f - 1.125925 * a + 0.000010 * v * v + 2.265811 * f * f + 0.007205 * v * f + 0.007658 * v * a + 0.0000001 * v * 1800.0 + 6.020526 * f * a - 0.039339 * v * f * a);
    g3leq = 2.678195 - 0.015906 * v - 3.961551 * f - 4.844458 * a + 0.000027 * v * v + 3.816645 * f * f + 0.046635 * v * f + 0.033541 * v * a + 28.458806 * f * a + 0.001495 * a * 1800.0 - 0.243340 * v * f * a - 0.008445 * f * a * 1800.0 + 0.000049 * v * f * a * 1800.0 - 4.0;
    g4leq = 1.0 - (2.678195 - 0.015906 * v - 3.961551 * f - 4.844458 * a + 0.000027 * v * v + 3.816645 * f * f + 0.046635 * v * f + 0.033541 * v * a + 28.458806 * f * a + 0.001495 * a * 1800.0 - 0.243340 * v * f * a - 0.008445 * f * a * 1800.0 + 0.000049 * v * f * a * 1800.0);
    g5leq = 2.775193 - 0.016948 * v - 2.123522 * f + 0.000042 * v * v + 5.430604 * f * f + 0.028112 * v * f - 0.0588322 * v * f * a + 0.000008 * v * a * 1800.0 - 6.3;

    obj[0] = f1;
    obj[1] = f2;
    obj[2] = f3;
    constr[0] = -1.0 * g1leq;
    constr[1] = -1.0 * g2leq;
    constr[2] = -1.0 * g3leq;
    constr[3] = -1.0 * g4leq;
    constr[4] = -1.0 * g5leq;

    return;
}
#endif

/************************************************************ TWO OBJECTIVE PROBLEMS ************************************************************/

/*  Test problem SCH1
    # of real variables = 1
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 0
    */

#ifdef sch1
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    obj[0] = pow(xreal[0], 2.0);
    obj[1] = pow((xreal[0] - 2.0), 2.0);
    return;
}
#endif

/*  Test problem SCH2
    # of real variables = 1
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 0
    */

#ifdef sch2
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    if (xreal[0] <= 1.0)
    {
        obj[0] = -xreal[0];
        obj[1] = pow((xreal[0] - 5.0), 2.0);
        return;
    }
    if (xreal[0] <= 3.0)
    {
        obj[0] = xreal[0] - 2.0;
        obj[1] = pow((xreal[0] - 5.0), 2.0);
        return;
    }
    if (xreal[0] <= 4.0)
    {
        obj[0] = 4.0 - xreal[0];
        obj[1] = pow((xreal[0] - 5.0), 2.0);
        return;
    }
    obj[0] = xreal[0] - 4.0;
    obj[1] = pow((xreal[0] - 5.0), 2.0);
    return;
}
#endif

/*  Test problem FON
    # of real variables = n
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 0
    */

#ifdef fon
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    double s1, s2;
    int i;
    s1 = s2 = 0.0;
    for (i = 0; i < nreal; i++)
    {
        s1 += pow((xreal[i] - (1.0 / sqrt((double)nreal))), 2.0);
        s2 += pow((xreal[i] + (1.0 / sqrt((double)nreal))), 2.0);
    }
    obj[0] = 1.0 - exp(-s1);
    obj[1] = 1.0 - exp(-s2);
    return;
}
#endif

/*  Test problem KUR
    # of real variables = 3
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 0
    */

#ifdef kur
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    int i;
    double res1, res2;
    res1 = -0.2 * sqrt((xreal[0] * xreal[0]) + (xreal[1] * xreal[1]));
    res2 = -0.2 * sqrt((xreal[1] * xreal[1]) + (xreal[2] * xreal[2]));
    obj[0] = -10.0 * (exp(res1) + exp(res2));
    obj[1] = 0.0;
    for (i = 0; i < 3; i++)
    {
        obj[1] += pow(fabs(xreal[i]), 0.8) + 5.0 * sin(pow(xreal[i], 3.0));
    }
    return;
}
#endif

/*  Test problem POL
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 0
    */

#ifdef pol
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    double a1, a2, b1, b2;
    a1 = 0.5 * sin(1.0) - 2.0 * cos(1.0) + sin(2.0) - 1.5 * cos(2.0);
    a2 = 1.5 * sin(1.0) - cos(1.0) + 2.0 * sin(2.0) - 0.5 * cos(2.0);
    b1 = 0.5 * sin(xreal[0]) - 2.0 * cos(xreal[0]) + sin(xreal[1]) - 1.5 * cos(xreal[1]);
    b2 = 1.5 * sin(xreal[0]) - cos(xreal[0]) + 2.0 * sin(xreal[1]) - 0.5 * cos(xreal[1]);
    obj[0] = 1.0 + pow((a1 - b1), 2.0) + pow((a2 - b2), 2.0);
    obj[1] = pow((xreal[0] + 3.0), 2.0) + pow((xreal[1] + 1.0), 2.0);
    return;
}
#endif

/*  Test problem VNT
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 3
    # of constraints = 0
    */

#ifdef vnt
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    obj[0] = 0.5 * (xreal[0] * xreal[0] + xreal[1] * xreal[1]) + sin(xreal[0] * xreal[0] + xreal[1] * xreal[1]);
    obj[1] = (pow((3.0 * xreal[0] - 2.0 * xreal[1] + 4.0), 2.0)) / 8.0 + (pow((xreal[0] - xreal[1] + 1.0), 2.0)) / 27.0 + 15.0;
    obj[2] = 1.0 / (xreal[0] * xreal[0] + xreal[1] * xreal[1] + 1.0) - 1.1 * exp(-(xreal[0] * xreal[0] + xreal[1] * xreal[1]));
    return;
}
#endif

/*  Test problem ZDT1
    # of real variables = 30
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 0
    */

#ifdef zdt1
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    double f1, f2, g, h;
    int i;
    f1 = xreal[0];
    g = 0.0;
    for (i = 1; i < 30; i++)
    {
        g += xreal[i];
    }
    g = 9.0 * g / 29.0;
    g += 1.0;
    h = 1.0 - sqrt(f1 / g);
    f2 = g * h;
    obj[0] = f1;
    obj[1] = f2;
    return;
}
#endif

/*  Test problem ZDT2
    # of real variables = 30
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 0
    */

#ifdef zdt2
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    double f1, f2, g, h;
    int i;
    f1 = xreal[0];
    g = 0.0;
    for (i = 1; i < 30; i++)
    {
        g += xreal[i];
    }
    g = 9.0 * g / 29.0;
    g += 1.0;
    h = 1.0 - pow((f1 / g), 2.0);
    f2 = g * h;
    obj[0] = f1;
    obj[1] = f2;
    return;
}
#endif

/*  Test problem ZDT3
    # of real variables = 30
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 0
    */

#ifdef zdt3
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    double f1, f2, g, h;
    int i;
    f1 = xreal[0];
    g = 0.0;
    for (i = 1; i < 30; i++)
    {
        g += xreal[i];
    }
    g = 9.0 * g / 29.0;
    g += 1.0;
    h = 1.0 - sqrt(f1 / g) - (f1 / g) * sin(10.0 * PI * f1);
    f2 = g * h;
    obj[0] = f1;
    obj[1] = f2;
    return;
}
#endif

/*  Test problem ZDT4
    # of real variables = 10
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 0
    */

#ifdef zdt4
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    double f1, f2, g, h;
    int i;
    f1 = xreal[0];
    g = 0.0;
    for (i = 1; i < 10; i++)
    {
        g += xreal[i] * xreal[i] - 10.0 * cos(4.0 * PI * xreal[i]);
    }
    g += 91.0;
    h = 1.0 - sqrt(f1 / g);
    f2 = g * h;
    obj[0] = f1;
    obj[1] = f2;
    return;
}
#endif

/*  Test problem ZDT5
    # of real variables = 0
    # of bin variables = 11
    # of bits for binvar1 = 30
    # of bits for binvar2-11 = 5
    # of objectives = 2
    # of constraints = 0
    */

#ifdef zdt5
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    int i, j;
    int u[11];
    int v[11];
    double f1, f2, g, h;
    for (i = 0; i < 11; i++)
    {
        u[i] = 0;
    }
    for (j = 0; j < 30; j++)
    {
        if (gene[0][j] == 1)
        {
            u[0]++;
        }
    }
    for (i = 1; i < 11; i++)
    {
        for (j = 0; j < 4; j++)
        {
            if (gene[i][j] == 1)
            {
                u[i]++;
            }
        }
    }
    f1 = 1.0 + u[0];
    for (i = 1; i < 11; i++)
    {
        if (u[i] < 5)
        {
            v[i] = 2 + u[i];
        }
        else
        {
            v[i] = 1;
        }
    }
    g = 0;
    for (i = 1; i < 11; i++)
    {
        g += v[i];
    }
    h = 1.0 / f1;
    f2 = g * h;
    obj[0] = f1;
    obj[1] = f2;
    return;
}
#endif

/*  Test problem ZDT6
    # of real variables = 10
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 0
    */

#ifdef zdt6
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    double f1, f2, g, h;
    int i;
    f1 = 1.0 - (exp(-4.0 * xreal[0])) * pow((sin(4.0 * PI * xreal[0])), 6.0);
    g = 0.0;
    for (i = 1; i < 10; i++)
    {
        g += xreal[i];
    }
    g = g / 9.0;
    g = pow(g, 0.25);
    g = 1.0 + 9.0 * g;
    h = 1.0 - pow((f1 / g), 2.0);
    f2 = g * h;
    obj[0] = f1;
    obj[1] = f2;
    return;
}
#endif

/*  Test problem BNH
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 2
    */

#ifdef bnh
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    obj[0] = 4.0 * (xreal[0] * xreal[0] + xreal[1] * xreal[1]);
    obj[1] = pow((xreal[0] - 5.0), 2.0) + pow((xreal[1] - 5.0), 2.0);
    constr[0] = 1.0 - (pow((xreal[0] - 5.0), 2.0) + xreal[1] * xreal[1]) / 25.0;
    constr[1] = (pow((xreal[0] - 8.0), 2.0) + pow((xreal[1] + 3.0), 2.0)) / 7.7 - 1.0;
    return;
}
#endif

/*  Test problem OSY
    # of real variables = 6
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 6
    */

#ifdef osy
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    obj[0] = -(25.0 * pow((xreal[0] - 2.0), 2.0) + pow((xreal[1] - 2.0), 2.0) + pow((xreal[2] - 1.0), 2.0) + pow((xreal[3] - 4.0), 2.0) + pow((xreal[4] - 1.0), 2.0));
    obj[1] = xreal[0] * xreal[0] + xreal[1] * xreal[1] + xreal[2] * xreal[2] + xreal[3] * xreal[3] + xreal[4] * xreal[4] + xreal[5] * xreal[5];
    constr[0] = (xreal[0] + xreal[1]) / 2.0 - 1.0;
    constr[1] = 1.0 - (xreal[0] + xreal[1]) / 6.0;
    constr[2] = 1.0 - xreal[1] / 2.0 + xreal[0] / 2.0;
    constr[3] = 1.0 - xreal[0] / 2.0 + 3.0 * xreal[1] / 2.0;
    constr[4] = 1.0 - (pow((xreal[2] - 3.0), 2.0)) / 4.0 - xreal[3] / 4.0;
    constr[5] = (pow((xreal[4] - 3.0), 2.0)) / 4.0 + xreal[5] / 4.0 - 1.0;
    return;
}
#endif

/*  Test problem SRN
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 2
    */

#ifdef srn
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    obj[0] = 2.0 + pow((xreal[0] - 2.0), 2.0) + pow((xreal[1] - 1.0), 2.0);
    obj[1] = 9.0 * xreal[0] - pow((xreal[1] - 1.0), 2.0);
    constr[0] = 1.0 - (pow(xreal[0], 2.0) + pow(xreal[1], 2.0)) / 225.0;
    constr[1] = 3.0 * xreal[1] / 10.0 - xreal[0] / 10.0 - 1.0;
    return;
}
#endif

/*  Test problem TNK
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 2
    */

#ifdef tnk
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    obj[0] = xreal[0];
    obj[1] = xreal[1];
    if (xreal[1] == 0.0)
    {
        constr[0] = -1.0;
    }
    else
    {
        constr[0] = xreal[0] * xreal[0] + xreal[1] * xreal[1] - 0.1 * cos(16.0 * atan(xreal[0] / xreal[1])) - 1.0;
    }
    constr[1] = 1.0 - 2.0 * pow((xreal[0] - 0.5), 2.0) + 2.0 * pow((xreal[1] - 0.5), 2.0);
    return;
}
#endif

/*  Test problem CTP1
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 2
    */

#ifdef ctp1
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    double g;
    g = 1.0 + xreal[1];
    obj[0] = xreal[0];
    obj[1] = g * exp(-obj[0] / g);
    constr[0] = obj[1] / (0.858 * exp(-0.541 * obj[0])) - 1.0;
    constr[1] = obj[1] / (0.728 * exp(-0.295 * obj[0])) - 1.0;
    return;
}
#endif

/*  Test problem CTP2
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 1
    */

#ifdef ctp2
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    double g;
    double theta, a, b, c, d, e;
    double exp1, exp2;
    theta = -0.2 * PI;
    a = 0.2;
    b = 10.0;
    c = 1.0;
    d = 6.0;
    e = 1.0;
    g = 1.0 + xreal[1];
    obj[0] = xreal[0];
    obj[1] = g * (1.0 - sqrt(obj[0] / g));
    exp1 = (obj[1] - e) * cos(theta) - obj[0] * sin(theta);
    exp2 = (obj[1] - e) * sin(theta) + obj[0] * cos(theta);
    exp2 = b * PI * pow(exp2, c);
    exp2 = fabs(sin(exp2));
    exp2 = a * pow(exp2, d);
    constr[0] = exp1 / exp2 - 1.0;
    return;
}
#endif

/*  Test problem CTP3
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 1
    */

#ifdef ctp3
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    double g;
    double theta, a, b, c, d, e;
    double exp1, exp2;
    theta = -0.2 * PI;
    a = 0.1;
    b = 10.0;
    c = 1.0;
    d = 0.5;
    e = 1.0;
    g = 1.0 + xreal[1];
    obj[0] = xreal[0];
    obj[1] = g * (1.0 - sqrt(obj[0] / g));
    exp1 = (obj[1] - e) * cos(theta) - obj[0] * sin(theta);
    exp2 = (obj[1] - e) * sin(theta) + obj[0] * cos(theta);
    exp2 = b * PI * pow(exp2, c);
    exp2 = fabs(sin(exp2));
    exp2 = a * pow(exp2, d);
    constr[0] = exp1 / exp2 - 1.0;
    return;
}
#endif

/*  Test problem CTP4
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 1
    */

#ifdef ctp4
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    double g;
    double theta, a, b, c, d, e;
    double exp1, exp2;
    theta = -0.2 * PI;
    a = 0.75;
    b = 10.0;
    c = 1.0;
    d = 0.5;
    e = 1.0;
    g = 1.0 + xreal[1];
    obj[0] = xreal[0];
    obj[1] = g * (1.0 - sqrt(obj[0] / g));
    exp1 = (obj[1] - e) * cos(theta) - obj[0] * sin(theta);
    exp2 = (obj[1] - e) * sin(theta) + obj[0] * cos(theta);
    exp2 = b * PI * pow(exp2, c);
    exp2 = fabs(sin(exp2));
    exp2 = a * pow(exp2, d);
    constr[0] = exp1 / exp2 - 1.0;
    return;
}
#endif

/*  Test problem CTP5
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 1
    */

#ifdef ctp5
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    double g;
    double theta, a, b, c, d, e;
    double exp1, exp2;
    theta = -0.2 * PI;
    a = 0.1;
    b = 10.0;
    c = 2.0;
    d = 0.5;
    e = 1.0;
    g = 1.0 + xreal[1];
    obj[0] = xreal[0];
    obj[1] = g * (1.0 - sqrt(obj[0] / g));
    exp1 = (obj[1] - e) * cos(theta) - obj[0] * sin(theta);
    exp2 = (obj[1] - e) * sin(theta) + obj[0] * cos(theta);
    exp2 = b * PI * pow(exp2, c);
    exp2 = fabs(sin(exp2));
    exp2 = a * pow(exp2, d);
    constr[0] = exp1 / exp2 - 1.0;
    return;
}
#endif

/*  Test problem CTP6
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 1
    */

#ifdef ctp6
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    double g;
    double theta, a, b, c, d, e;
    double exp1, exp2;
    theta = 0.1 * PI;
    a = 40.0;
    b = 0.5;
    c = 1.0;
    d = 2.0;
    e = -2.0;
    g = 1.0 + xreal[1];
    obj[0] = xreal[0];
    obj[1] = g * (1.0 - sqrt(obj[0] / g));
    exp1 = (obj[1] - e) * cos(theta) - obj[0] * sin(theta);
    exp2 = (obj[1] - e) * sin(theta) + obj[0] * cos(theta);
    exp2 = b * PI * pow(exp2, c);
    exp2 = fabs(sin(exp2));
    exp2 = a * pow(exp2, d);
    constr[0] = exp1 / exp2 - 1.0;
    return;
}
#endif

/*  Test problem CTP7
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 1
    */

#ifdef ctp7
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    double g;
    double theta, a, b, c, d, e;
    double exp1, exp2;
    theta = -0.05 * PI;
    a = 40.0;
    b = 5.0;
    c = 1.0;
    d = 6.0;
    e = 0.0;
    g = 1.0 + xreal[1];
    obj[0] = xreal[0];
    obj[1] = g * (1.0 - sqrt(obj[0] / g));
    exp1 = (obj[1] - e) * cos(theta) - obj[0] * sin(theta);
    exp2 = (obj[1] - e) * sin(theta) + obj[0] * cos(theta);
    exp2 = b * PI * pow(exp2, c);
    exp2 = fabs(sin(exp2));
    exp2 = a * pow(exp2, d);
    constr[0] = exp1 / exp2 - 1.0;
    return;
}
#endif

/*  Test problem CTP8
    # of real variables = 2
    # of bin variables = 0
    # of objectives = 2
    # of constraints = 2
    */

#ifdef ctp8
void test_problem(double *xreal, double *xbin, int **gene, double *obj, double *constr)
{
    double g;
    double theta, a, b, c, d, e;
    double exp1, exp2;
    g = 1.0 + xreal[1];
    obj[0] = xreal[0];
    obj[1] = g * (1.0 - sqrt(obj[0] / g));
    theta = 0.1 * PI;
    a = 40.0;
    b = 0.5;
    c = 1.0;
    d = 2.0;
    e = -2.0;
    exp1 = (obj[1] - e) * cos(theta) - obj[0] * sin(theta);
    exp2 = (obj[1] - e) * sin(theta) + obj[0] * cos(theta);
    exp2 = b * PI * pow(exp2, c);
    exp2 = fabs(sin(exp2));
    exp2 = a * pow(exp2, d);
    constr[0] = exp1 / exp2 - 1.0;
    theta = -0.05 * PI;
    a = 40.0;
    b = 2.0;
    c = 1.0;
    d = 6.0;
    e = 0.0;
    exp1 = (obj[1] - e) * cos(theta) - obj[0] * sin(theta);
    exp2 = (obj[1] - e) * sin(theta) + obj[0] * cos(theta);
    exp2 = b * PI * pow(exp2, c);
    exp2 = fabs(sin(exp2));
    exp2 = a * pow(exp2, d);
    constr[1] = exp1 / exp2 - 1.0;
    return;
}
#endif
