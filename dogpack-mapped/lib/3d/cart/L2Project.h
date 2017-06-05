#ifndef L2Project_h
#define L2Project_h
#include "tensors.h"

void L2ProjectGradAdd(const int istart, 
		      const int iend, 
		      const int jstart, 
		      const int jend,
		      const int kstart,
		      const int kend,
                      const dTensorBC5& qin,
                      const dTensorBC5& auxin, 
		      dTensorBC5& Fout,
                      void (*FluxFunc)(const dTensor2& xpts, const dTensor2& Q,
				       const dTensor2& Aux, dTensor3& flux));

void L2ProjectGradAdd(const int istart, 
                      const int iend, 
                      const int jstart, 
                      const int jend,
		      const int kstart,
		      const int kend,
                      const int QuadOrder, 
                      const int BasisOrder_qin,
                      const int BasisOrder_auxin,
                      const int BasisOrder_fout,
                      const dTensorBC5* qin,
                      const dTensorBC5* auxin, 
                      dTensorBC5* fout,
                      void Func(const dTensor2&,const dTensor2&,
                                const dTensor2&,dTensor3&));


void L2Project(const int istart, 
	       const int iend, 
	       const int jstart, 
	       const int jend,
	       const int kstart, 
	       const int kend,
	       const dTensorBC5& q,
	       const dTensorBC5& aux, 
	       dTensorBC5& Fout,
	       void (*Func)(const dTensor2& xpts,
			    const dTensor2& qvals,
			    const dTensor2& auxvals,
			    dTensor2& source));

void L2Project(const int istart,
	       const int iend,
	       const int jstart,
	       const int jend,
	       const int kstart,
	       const int kend,
	       const int QuadOrder,
	       const int BasisOrder_qin,
	       const int BasisOrder_auxin,
	       const int BasisOrder_fout,    
	       const dTensorBC5* qin,
	       const dTensorBC5* auxin,
	       dTensorBC5* fout,
	       void (*Func)(const dTensor2&,const dTensor2&,
			    const dTensor2&,dTensor2&));

void L2ProjectInitialCond(const int istart,
			  const int iend,
			  const int jstart,
			  const int jend,
			  const int kstart,
			  const int kend,
			  const int QuadOrder,
			  const int BasisOrder_qin,
			  const int BasisOrder_auxin,
			  const int BasisOrder_fout,    
			  const dTensorBC5* qin,
			  const dTensorBC5* auxin,
			  dTensorBC5* fout,
			  void (*Func)(const dTensor2&,const dTensor2&,
				       const dTensor2&,dTensor2&));


void L2Project_extra(const int istart,
		     const int iend, 
		     const int jstart, 
		     const int jend,
		     const int kstart,
		     const int kend,
		     const dTensorBC5& q,
		     const dTensorBC5& aux, 
		     dTensorBC5& Fout,
		     void (*Func)(const dTensor2&,const dTensor2&,
				  const dTensor2&,dTensor2&, void* data),
		     void* data);

void L2Project_extra(const int istart,
		     const int iend,
		     const int jstart,
		     const int jend,
		     const int kstart,
		     const int kend,
		     const int QuadOrder,
		     const int BasisOrder_qin,
		     const int BasisOrder_auxin,
		     const int BasisOrder_fout,    
		     const dTensorBC5* qin,
		     const dTensorBC5* auxin,
		     dTensorBC5* fout,
		     void (*Func)(const dTensor2&,const dTensor2&,
				  const dTensor2&,dTensor2&, void* data),
		     void* data);

#endif // L2Project_h