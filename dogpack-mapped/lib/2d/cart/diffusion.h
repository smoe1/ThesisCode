#ifndef DIFFUSION_H
#define DIFFUSION_H
#include "constants.h"

const double AU[6][6] = {{1.0, sq3, 0.0, 0.0, sq5, 0.0}, 
			 {-sq3, 3, 0.0, 0.0, sq3*sq5, 0.0}, 
			 {0.0, 0.0, 1.0, sq3, 0.0, 0.0}, 
			 {0.0, 0.0, -sq3, 3, 0.0, 0.0}, 
			 {sq5, -sq3*sq5, 0.0, 0.0, 5, 0.0}, 
			 {0.0, 0.0, 0.0, 0.0, 0.0, 1.0}};
  
const double BU[6][6] = {{1.0, sq3, 0.0, 0.0, sq5, 0.0},
			 {-sq3, -3, 0.0, 0.0, -sq3*sq5, 0.0},
			 {0.0, 0.0, 1.0, sq3, 0.0, 0.0},
			 {0.0, 0.0, -sq3, -3, 0.0, 0.0},
			 {sq5, sq3*sq5, 0.0, 0.0, 5.0, 0.0},
			 {0.0, 0.0, 0.0, 0.0, 0.0, 1.0}};
  
const double CU[6][6] = {{1.0, 0.0, sq3, 0.0, 0.0, sq5},
			 {0.0, 1.0, 0.0, sq3, 0.0, 0.0},
			 {-sq3, 0.0, 3.0, 0.0, 0.0, sq3*sq5},
			 {0.0, -sq3, 0.0, 3.0, 0.0, 0.0},
			 {0.0, 0.0, 0.0, 0.0, 1.0, 0.0},
			 {sq5, 0.0, -sq3*sq5, 0.0, 0.0, 5.0}};
  
const double DU[6][6] = {{1.0, 0.0, sq3, 0.0, 0.0, sq5},
			 {0.0, 1.0, 0.0, sq3, 0.0, 0.0},
			 {-sq3, 0.0, -3, 0.0, 0.0, -sq3*sq5},
			 {0.0, -sq3, 0.0, -3, 0.0, 0.0},
			 {0.0, 0.0, 0.0, 0.0, 1.0, 0.0},
			 {sq5, 0.0, sq3*sq5, 0.0, 0.0, 5.0}};
  
const double AV[6][6] = {{1.0, -sq3, 0.0, 0.0, sq5, 0.0},
			 {sq3, -3, 0.0, 0.0, sq3*sq5, 0.0},
			 {0.0, 0.0, 1.0, -sq3, 0.0, 0.0},
			 {0.0, 0.0, sq3, -3, 0.0, 0.0},
			 {sq5, -sq3*sq5, 0.0, 0.0, 5.0, 0.0},
			 {0.0, 0.0, 0.0, 0.0, 0.0, 1.0}};
  
const double BV[6][6] = {{1.0, -sq3, 0.0, 0.0, sq5, 0.0},
			 {sq3, 3.0, 0.0, 0.0, -sq3*sq5, 0.0},
			 {0.0, 0.0, 1.0, -sq3, 0.0, 0.0},
			 {0.0, 0.0, sq3, 3.0, 0.0, 0.0},
			 {sq5, sq3*sq5, 0.0, 0.0, 5.0, 0.0},
			 {0.0, 0.0, 0.0, 0.0, 0.0, 1.0}};
  
const double AW[6][6] = {{1.0, 0.0, -sq3, 0.0, 0.0, sq5},
			 {0.0, 1.0, 0.0, -sq3, 0.0, 0.0},
			 {sq3, 0.0, -3, 0.0, 0.0, sq3*sq5},
			 {0.0, sq3, 0.0, -3, 0.0, 0.0},
			 {0.0, 0.0, 0.0, 0.0, 1.0, 0.0},
			 {sq5, 0.0, -sq3*sq5, 0.0, 0.0, 5.0}};
  
const double BW[6][6] = {{1.0, 0.0, -sq3, 0.0, 0.0, sq5},
			 {0.0, 1.0, 0.0, -sq3, 0.0, 0.0},
			 {sq3, 0.0, 3.0, 0.0, 0.0, -sq3*sq5},
			 {0.0, sq3, 0.0, 3.0, 0.0, 0.0},
			 {0.0, 0.0, 0.0, 0.0, 1.0, 0.0},
			 {sq5, 0.0, sq3*sq5, 0.0, 0.0, 5.0}};

#endif