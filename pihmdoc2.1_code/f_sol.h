#ifndef F_SOL_H
#define F_SOL_H

#include <stdio.h>

#include "sundials_types.h"


int f_sol(FILE **, realtype t, N_Vector CV_Y, N_Vector CV_Ydot, void *DS);


#endif
