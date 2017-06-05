#include "dogdefs.h"
#include "DogSolverCart4.h"

// Function that is called after a full time step
// (i.e., after all stages are complete)
void AfterFullTimeStep(DogSolverCart4& solver)
{

    const double t  = solver.get_state().get_time();
    dTensorBC6& aux = solver.fetch_state().fetch_aux();
    dTensorBC6& q   = solver.fetch_state().fetch_q();

    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int mz   = q.getsize(3);
    const int mw   = q.getsize(4);
    const int meqn = q.getsize(5);
    const int kmax = q.getsize(6);
    const int mbc  = q.getmbc();
    const int maux = aux.getsize(5);
}
