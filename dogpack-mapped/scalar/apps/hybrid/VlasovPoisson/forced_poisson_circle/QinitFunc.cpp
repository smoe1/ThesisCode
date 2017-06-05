#include <cmath>
#include "dogdefs.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
// Exact solution: f = A(t) * (1-x^2-y^2) * exp( -vx^2 ) / sqrt(pi)
//
// Where, A(t) = ( 0.75 + 0.25*cos(2*pi*t) )
//
// Note, A(0) = 1, so it doesn't show up here.
//
void QinitFunc(
    const dTensor2* vel_vec,
    const dTensor2& xpts, dTensor2& qvals)
{

    const int numpts  = qvals.getsize(1);
    const int meqn    = qvals.getsize(2);
    const double one_over_sqpi = 1.0/sqrt(pi);
    for (int i=1; i<=numpts; i++)
    for (int m=1; m<=meqn; m++)
    {

        const double x  = xpts.get(i,1);
        const double y  = xpts.get(i,2);

        const double vx = vel_vec->get(m,1);
        const double vy = vel_vec->get(m,2);

        double r2 = pow(x,2) + pow(y,2);
        qvals.set(i,m, (1.0 - r2)*exp(-vx*vx ) * one_over_sqpi );

    }
}

void ExactQ(
    const double t, 
    const dTensor2* vel_vec,
    const dTensor2& xpts, const dTensor2& qvals, 
    const dTensor2& auxvals, dTensor2& qex)
{

    const int numpts  = qex.getsize(1);
    const int meqn    = qex.getsize(2);
    const double one_over_sqpi = 1.0/sqrt(pi);

    const double A  = (0.75 + 0.25 * cos(2.0*pi*t) );
    for (int i=1; i<=numpts; i++)
    for (int m=1; m<=meqn; m++)
    {

        const double x  = xpts.get(i,1);
        const double y  = xpts.get(i,2);

        const double vx = vel_vec->get(m,1);
        const double vy = vel_vec->get(m,2);

        double r2 = pow(x,2) + pow(y,2);
        qex.set(i,m, A*(1.0 - r2)*exp(-vx*vx ) * one_over_sqpi );

    }
}

