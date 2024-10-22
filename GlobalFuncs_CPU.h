#pragma once

#include "MacroVaribles.h"
#include "StructComplex.h"
#include "StructSrc.h"
#include "Operators.h"
#include "CoeffGenerator.h"
#include "RootsFinder.h"
// #include "../device/device_macros.h"

static const double RelErrAllowed = 1e-10;

template<isFloating f_T>
void MarginBatch(int blockIdx_x, int blockDim, int const * srclist, src_ext_t<f_T>* srcs);



template<isFloating f_T>
void SolveBatch(int blockIdx_x, int blockDim, int const * srclist,src_ext_t<f_T>*srcs,const src_params_t<f_T>* src_params, int batchidx, bool muted = false, bool check=true);

// size: <<<(NWALKERS*NSRCS,1, 1),BATCH_SIZE>>>
template<isFloating f_T>
void FindNextBatch3(int blockIdx_x, int blockDim, int* srclist, src_ext_t<f_T>* srcs, int batch_idx, bool muted = false);

template <isFloating f_T>
void AreaErrBatchSrc4(int blockIdx_x, int blockDim, int* srclist, src_ext_t<f_T>* srcs, int batch_idx, bool muted);

// size <<<[(NSRCS-1) // 64 +1,1,1],64>>>
template<isFloating f_T>
void AreaErrBatchSrcCross(int blockIdx_x, int* srclist, src_ext_t<f_T>* srcs,const src_params_t<f_T>* src_params, int batch_idx, int prev_Ncal, bool muted = false);

template <isFloating f_T>
void SumArea0(int blockIdx_x, int blockDim, int* srclist, src_ext_t<f_T>* srcs, int batchidx, bool muted = false, f_T RELTOL=1e-4);

template <isFloating f_T>
void SumArea3(int blockIdx_x, int blockDim, int* srclist, src_ext_t<f_T>* srcs, int batchidx, bool muted = false, f_T RELTOL=1e-4);


static const float C_Q = 6;
static const float C_G = 2;
static const float C_P = 2;

template<isFloating f_T>
void QuadruTst(int blockIdx_x, src_ext_t<f_T>* srcs, const src_params_t<f_T>* src_params, const int Nsrcs=NSRCS, f_T RELTOL=1e-4);

template <isFloating f_T>
void SetPath(int blockIdx_z, src_ext_t<f_T>*srcs, const f_T* time, const src_params_t<f_T>* move_params, int Nsrcs=NSRCS);

// size:<<<NaxisY,NaxisX>>>; 
template <isFloating f_T>
void SetArray(int blockIdx_x, int threadIdx_x, src_ext_t<f_T>*srcs, const src_params_t<f_T>* move_params, const int NaxisX, const int NaxisY, const f_T xmin,const f_T xmax,const f_T ymin,const f_T ymax, int Nsrcs=NSRCS);

// size:<<<NaxisY,NaxisX>>>; 
template <isFloating f_T>
void SetRhoLoc(int blockIdx_z, src_ext_t<f_T>*srcs, const f_T* rhos, const complex_t<f_T>* zetas, int Nsrcs=NSRCS);

// size:<<<1,1>>>; 
template <isFloating f_T>
void SetPoint(src_ext_t<f_T>*srcs, const src_params_t<f_T>* move_params, const f_T x, const f_T y);

template <isFloating f_T>
void AdaptiveLocLoop(int blockIdx_x, int blockDim, int* srclist, src_ext_t<f_T>*srcs, int batchidx, bool muted = false, f_T RELTOL=1e-4);

// size <<<[(NWALKERS*NSRCS-1) // 64 +1,1,1],64>>>
template<isFloating f_T>
void PhysicalTest(int blockIdx_x, int* srclist, src_ext_t<f_T>*srcs, int batch_idx, int prev_Ncal, bool muted = false);


// size: <<<(Ncal-1) / 64+1 , min(64,Ncal)>>>
template<isFloating f_T>
void SuccessCheck(int blockIdx_x, int* Ncal, int* srclist, src_ext_t<f_T>* srcs, int prev_Ncal, bool muted = false);
