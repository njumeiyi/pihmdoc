#include <stdio.h>

#include "sundials_types.h"
#include "pihm.h"

realtype Interpolation(TSD *Data, realtype t);

realtype returnVal(realtype rArea, realtype rPerem, realtype eqWid,realtype ap_Bool);

realtype CS_AreaOrPerem(int rivOrder, realtype rivDepth, realtype rivCoeff, realtype a_pBool);

OverlandFlow(realtype **flux, int loci, int locj, realtype avg_y, realtype grad_y, realtype avg_sf, realtype crossA, realtype avg_rough);

OLFeleToriv(realtype eleYtot,realtype EleZ,realtype cwr,realtype rivZmax,realtype rivYtot,realtype **fluxriv,int loci,int locj,realtype length);

realtype avgY(realtype diff, realtype yi, realtype yinabr);

realtype effKV(realtype ksatFunc,realtype gradY,realtype macKV,realtype KV,realtype areaF);

realtype effKH(int mp,realtype tmpY, realtype aqDepth, realtype MacD, realtype MacKsatH, realtype areaF, realtype ksatH);



