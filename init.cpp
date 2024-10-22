#include"init.h"

// Size: <<< (Nsrcs-1)/64+1 , min(64,Nsrcs)>>>
template<isFloating f_T>
_Global void SrcsInit(int blockIdx_x, src_ext_t<f_T>* srcs)
{
    // int srcidx = blockIdx.x * blockDim.x + threadIdx.x;
    int srcidx = blockIdx_x;
    if(srcidx<NSRCS)
    {
        srcs[srcidx].SolveSucceed = false;
        srcs[srcidx].Break = false;
        srcs[srcidx].Area = 0;
        srcs[srcidx].Err_A = 0;        
    }

    // changed 赋值在核函数内进行
    // for(int i=0;i<NPOINTS;i++)
    // {
    //     srcs[srcidx].margin_pts[i].changed = 0;
    // }
}
template _Global void SrcsInit(int blockIdx_x, src_ext_t<float>* srcs);
template _Global void SrcsInit(int blockIdx_x, src_ext_t<double>* srcs);

_Global void ResTest()
{

}



template <isFloating f_T>
void Twinkle<f_T>::malloc(int device_num)
{
    // // dev.init(device_num);
    // // stream = dev.yield_streams();
    // h_srcs.init(dev,Nsrcs,on_dev_t::host);
    // d_srcs.init(dev,Nsrcs,on_dev_t:: dev);
    // h_CC.init(dev,NCRIT,on_dev_t::host);
    // d_CC.init(dev,NCRIT,on_dev_t:: dev);        

    // h_srcidx.init(dev,Nsrcs,on_dev_t::host);
    // d_srcidx.init(dev,Nsrcs,on_dev_t:: dev);
    // h_params.init(dev,1,on_dev_t::host);
    // d_params.init(dev,1,on_dev_t:: dev);
    // h_Nfail.init(dev,1,on_dev_t::host);
    // d_Nfail.init(dev,1,on_dev_t:: dev);
    // h_Ncal.init(dev,1,on_dev_t::host);
    // d_Ncal.init(dev,1,on_dev_t:: dev);

    // h_Mag.init(dev,Nsrcs,on_dev_t::host);
    // d_Mag.init(dev,Nsrcs,on_dev_t:: dev);
    // h_Err.init(dev,Nsrcs,on_dev_t::host);
    // d_Err.init(dev,Nsrcs,on_dev_t:: dev);

    // h_time.init(dev,Nsrcs,on_dev_t::host);
    // d_time.init(dev,Nsrcs,on_dev_t:: dev);
    // // printf("malloc finished\n");

    h_srcs = new src_ext_t<f_T>[Nsrcs];
    // h_CC = new CC_t<f_T>[NCRIT];

    h_srcidx = new int[Nsrcs];
    h_params = new src_params_t<f_T>[1];
    h_Nfail = new int[1];
    h_Ncal = new int[1];

    h_Mag = new f_T[Nsrcs];
    h_Err = new f_T[Nsrcs];

    h_time = new f_T[Nsrcs];


}
template void Twinkle<float> ::malloc(int device_num);
template void Twinkle<double>::malloc(int device_num);


template <isFloating f_T>
void Twinkle<f_T>::free()
{
    // dev.sync_stream(stream);

    // h_srcs.free(dev);
    // d_srcs.free(dev);

    // h_CC.free(dev);
    // d_CC.free(dev);        


    // h_srcidx.free(dev);
    // d_srcidx.free(dev);
    // h_params.free(dev);
    // d_params.free(dev);
    // h_Nfail.free(dev);
    // d_Nfail.free(dev);
    // h_Ncal.free(dev);
    // d_Ncal.free(dev);

    // h_Mag.free(dev);
    // d_Mag.free(dev);
    // h_Err.free(dev);
    // d_Err.free(dev);

    // h_time.free(dev);
    // d_time.free(dev);

    delete[] h_srcs;
    // delete[] h_CC;
    delete[] h_srcidx;
    delete[] h_params;
    delete[] h_Nfail;
    delete[] h_Ncal;
    delete[] h_Mag;
    delete[] h_Err;
    delete[] h_time;

    // 将指针设置为 nullptr，防止悬空指针
    h_srcs = nullptr;
    // h_CC = nullptr;
    h_srcidx = nullptr;
    h_params = nullptr;
    h_Nfail = nullptr;
    h_Ncal = nullptr;
    h_Mag = nullptr;
    h_Err = nullptr;
    h_time = nullptr;

}
template void Twinkle<float> ::free();
template void Twinkle<double>::free();

// initialize for each solving step. no memory operation.
template <isFloating f_T>
void Twinkle<f_T>::init()
{
    // dev.sync_stream(stream);
    // Ncal = Nsrcs;
    // *h_Ncal = Nsrcs;
    // *h_Nfail = 0;
    // *(this->h_params) = *params;
    // dev.cp_a(d_Ncal,h_Ncal,sizeof(int),stream);
    // dev.cp_a(d_Nfail,h_Nfail,sizeof(int),stream);
    // dev.cp_a(d_params,h_params,sizeof(src_params_t<f_T>),stream);

    Ncal = Nsrcs;
    *h_Ncal = Nsrcs;
    *h_Nfail = 0;
    for(int i=0;i<Nsrcs;i++)
    {
        h_srcidx[i] = i;
    }


    // h_Ncal.cp_to(dev,d_Ncal,stream);
    // h_Nfail.cp_to(dev,d_Nfail,stream);
    // h_srcidx.cp_to(dev,d_srcidx,stream);




    // dev.sync_stream(stream);
    // // dev.sync_stream(stream);
    // // printf("initialize finished\n");
    // auto lpar = std::make_tuple
    //         < dim3, dim3, int > ( (Nsrcs-1)/BATCH_SIZE + 1 , min(BATCH_SIZE,Nsrcs), 0 );
    // dev.launch( SrcsInit<f_T>, lpar, stream, d_srcs.data);

    // // pt mode
    // // SetPoint<<<1,1>>>(d_srcs,d_params,-0.0000272125,1.4e-6);
    // SetPoint<<<1,1>>>(d_srcs,d_params,0.2,0.2);
    // // pt mode over


    for(int blockIdx=0;blockIdx<Nsrcs;blockIdx++)
    {
        SrcsInit<f_T>(blockIdx, h_srcs);
    }

}
template void Twinkle<float> ::init();
template void Twinkle<double>::init();

template <isFloating f_T>
void Twinkle<f_T>::set_time()
{
    // h_time.cp_to(dev,d_time,stream);
}
template void Twinkle<float> ::set_time();
template void Twinkle<double>::set_time();


template <isFloating f_T>
void Twinkle<f_T>::set_time(f_T* time)
{
    for(int i=0;i<Nsrcs;i++)
    {
        h_time[i] = time[i];
    }
    // h_time.cp_to(dev,d_time,stream);
}
template void Twinkle<float> ::set_time(float* time);
template void Twinkle<double>::set_time(double* time);


template <isFloating f_T>
void Twinkle<f_T>::set_params()
{
    // h_params.cp_to(dev,d_params,stream);
}
template void Twinkle<float> ::set_params();
template void Twinkle<double>::set_params();

template <isFloating f_T>
void Twinkle<f_T>::set_params(src_params_t<f_T>* params)
{
    *h_params = *params;
    // h_params.cp_to(dev,d_params,stream);
}
template void Twinkle<float> ::set_params(src_params_t<float>* params);
template void Twinkle<double>::set_params(src_params_t<double>* params);

template <isFloating f_T>
void Twinkle<f_T>::set_path()
{

    // dev.sync_stream(stream);

    // auto lpar = std::make_tuple
    //             < dim3, dim3, int > ( dim3(1,1,((Nsrcs-1)/BATCH_SIZE)+1), min(Nsrcs,BATCH_SIZE), 0 );

    // // printf("gridDimz: %d, blockDim.x: %d\n", ((Nsrcs-1)/BATCH_SIZE)+1, min(Nsrcs,BATCH_SIZE) );

    // dev.launch( SetPath<f_T>, lpar, stream, d_srcs.data, d_time.data, d_params.data, this->Nsrcs);

    // dev.sync_stream(stream);


    for(int blockIdx=0;blockIdx<Nsrcs;blockIdx++)
    {
        SetPath<f_T>(blockIdx, h_srcs, h_time, h_params, this->Nsrcs);
    }

}
template void Twinkle<float>::set_path();
template void Twinkle<double>::set_path();


template <isFloating f_T>
void Twinkle<f_T>::set_path(f_T* time)
{

    this->set_time(time);
    this->set_path();

    // dev.sync_stream(stream);

    // auto lpar = std::make_tuple
    //             < dim3, dim3, int > ( dim3(1,1,((Nsrcs-1)/BATCH_SIZE)+1), min(Nsrcs,BATCH_SIZE), 0 );

    // // printf("gridDimz: %d, blockDim.x: %d\n", ((Nsrcs-1)/BATCH_SIZE)+1, min(Nsrcs,BATCH_SIZE) );

    // dev.launch( SetPath<f_T>, lpar, stream, d_srcs.data, d_time.data, d_params.data, this->Nsrcs);

    // dev.sync_stream(stream);

}
template void Twinkle<float>::set_path(float* time);
template void Twinkle<double>::set_path(double* time);

template <isFloating f_T>
void Twinkle<f_T>::set_path(src_params_t<f_T>* params)
{
    // *h_params.data = *params;
    // h_params.cp_to(dev,d_params,stream);
    
    this->set_params(params);
    this->set_path();

    // dev.sync_stream(stream);

    // auto lpar = std::make_tuple
    //             < dim3, dim3, int > ( dim3(1,1,((Nsrcs-1)/BATCH_SIZE)+1), min(Nsrcs,BATCH_SIZE), 0 );

    // // printf("gridDimz: %d, blockDim.x: %d\n", ((Nsrcs-1)/BATCH_SIZE)+1, min(Nsrcs,BATCH_SIZE) );

    // dev.launch( SetPath<f_T>, lpar, stream, d_srcs.data, d_time.data, d_params.data, this->Nsrcs);

    // dev.sync_stream(stream);

}
template void Twinkle<float>::set_path(src_params_t<float>* params);
template void Twinkle<double>::set_path(src_params_t<double>* params);


template <isFloating f_T>
void Twinkle<f_T>::set_path(f_T* time, src_params_t<f_T>* params)
{
    // *h_params.data = *params;
    // h_params.cp_to(dev,d_params,stream);
    
    this->set_time(time);
    this->set_params(params);
    this->set_path();

    // dev.sync_stream(stream);

    // auto lpar = std::make_tuple
    //             < dim3, dim3, int > ( dim3(1,1,((Nsrcs-1)/BATCH_SIZE)+1), min(Nsrcs,BATCH_SIZE), 0 );

    // // printf("gridDimz: %d, blockDim.x: %d\n", ((Nsrcs-1)/BATCH_SIZE)+1, min(Nsrcs,BATCH_SIZE) );

    // dev.launch( SetPath<f_T>, lpar, stream, d_srcs.data, d_time.data, d_params.data, this->Nsrcs);

    // dev.sync_stream(stream);

}
template void Twinkle<float>::set_path(float* time, src_params_t<float>* params);
template void Twinkle<double>::set_path(double* time, src_params_t<double>* params);


// template <isFloating f_T>
// void Twinkle<f_T>::set_array(int NaxisX, int NaxisY, f_T xmin, f_T xmax, f_T ymin, f_T ymax)
// {
//     auto lpar = std::make_tuple
//             < dim3, dim3, int > ( NaxisX, NaxisY, 0 );
//     dev.launch( SetArray<f_T>, lpar, stream, d_srcs.data,d_params.data, NaxisX, NaxisY, xmin,xmax,ymin,ymax,this->Nsrcs);
//     // dev.launch( SetArray<f_T>, lpar, stream, d_srcs.data,d_params.data, -2e-2,4e-2,0e-3,1e-2);

// }

// template <isFloating f_T>
// void Twinkle<f_T>::set_rho_loc(f_T* rhos, complex_t<f_T>* zetas)
// {
//     array_t<f_T> h_rhos, d_rhos;
//     array_t<complex_t<f_T>> h_zetas, d_zetas;

//     // h_srcs.init(dev,Nsrcs,on_dev_t::host);
//     h_rhos.init (dev,Nsrcs,on_dev_t::host);
//     d_rhos.init (dev,Nsrcs,on_dev_t:: dev);
//     h_zetas.init(dev,Nsrcs,on_dev_t::host);
//     d_zetas.init(dev,Nsrcs,on_dev_t:: dev);
//     for(int i=0;i<Nsrcs;i++)
//     {
//         h_rhos[i] = rhos[i];
//         h_zetas[i] = zetas[i];
//     }
//     h_rhos.cp_to(dev,d_rhos,stream);
//     h_zetas.cp_to(dev,d_zetas,stream);

//     auto lpar = std::make_tuple
//                 < dim3, dim3, int > ( dim3(1,1,((Nsrcs-1)/BATCH_SIZE)+1), min(Nsrcs,BATCH_SIZE), 0 );

//     // printf("gridDimz: %d, blockDim.x: %d\n", ((Nsrcs-1)/BATCH_SIZE)+1, min(Nsrcs,BATCH_SIZE) );

//     dev.launch( SetRhoLoc<f_T>, lpar, stream, d_srcs.data, d_rhos.data, d_zetas.data, this->Nsrcs);

//     h_rhos.free(dev);
//     d_rhos.free(dev);
//     h_zetas.free(dev);
//     d_zetas.free(dev);
// }

template <isFloating f_T>
void Twinkle<f_T>::solve()
{

    init();

    // auto lpar = std::make_tuple
    //         < dim3, dim3, int > ( Nsrcs, BATCH_SIZE, 0 );


    // lpar = std::make_tuple < dim3, dim3, int > ( (Nsrcs-1) / BATCH_SIZE+1,min(BATCH_SIZE,Nsrcs), 0 );
    // dev.launch(QuadruTst<f_T>, lpar, stream, d_srcs.data,d_params.data, this->Nsrcs);

    for(int blockIdx=0;blockIdx<Nsrcs;blockIdx++)
    {
        QuadruTst<f_T>(blockIdx,h_srcs,h_params,Nsrcs,RELTOL);
    }

    for(int Batchidx =0;Batchidx<(NPOINTS/BATCH_SIZE);Batchidx++)
    // for(int Batchidx =0;Batchidx<14;Batchidx++)
    {
        // printf("batchidx: %d\n",Batchidx);
        // *h_Ncal.data = 0;
        // h_Ncal.cp_to(dev,d_Ncal,stream);
        // // printf("check1: Batchidx: %d, Ncal: %d\n",Batchidx,Ncal);
        // lpar = std::make_tuple < dim3, dim3, int > ( (Ncal-1) / BATCH_SIZE+1 , min(BATCH_SIZE,Ncal), 0 );

        // dev.launch(SuccessCheck<f_T>, lpar, stream, d_Ncal.data, d_srcidx.data, d_srcs.data, Ncal, muted);
        // d_Ncal.cp_to(dev,h_Ncal,stream);
        
        // dev.sync_stream( stream );
        // Ncal = *h_Ncal.data;


        *h_Ncal = 0;
        for(int blockIdx=0;blockIdx<Nsrcs;blockIdx++)
        {
            SuccessCheck<f_T>(blockIdx, h_Ncal, h_srcidx, h_srcs, Ncal, muted);
        }
        Ncal = *h_Ncal;


        // printf("check2: Batchidx: %d, Ncal: %d\n",Batchidx,Ncal);

        if(Ncal==0){break;}

        if(Batchidx==0)
        {
            // lpar = std::make_tuple < dim3, dim3, int > ( Ncal,BATCH_SIZE, 0 );
            // dev.launch(MarginBatch<f_T>, lpar, stream, d_srcidx.data, d_srcs.data);

            for(int blockIdx=0;blockIdx<Ncal;blockIdx++)
            {
                MarginBatch<f_T>(blockIdx, BATCH_SIZE, h_srcidx, h_srcs);
            }
        }
        else
        {
            // lpar = std::make_tuple < dim3, dim3, int > ( Ncal,min((BATCH_SIZE*Batchidx),512),(BATCH_SIZE)*(Batchidx)*sizeof(f_T) + (BATCH_SIZE)*sizeof(int) );
            // dev.launch(AdaptiveLocLoop<f_T>, lpar, stream, d_srcidx.data,d_srcs.data,Batchidx,muted);
            for(int blockIdx=0;blockIdx<Ncal;blockIdx++)
            {
                AdaptiveLocLoop<f_T>(blockIdx, min((BATCH_SIZE*Batchidx),512), h_srcidx, h_srcs, Batchidx, muted, RELTOL);
            }            
        }

        // if(Batchidx==1){break;}

        
        // lpar = std::make_tuple < dim3, dim3, int > ( Ncal,BATCH_SIZE, 0 );
        // dev.launch(SolveBatch<f_T>, lpar, stream, d_srcidx.data,d_srcs.data,d_params.data,Batchidx,muted,true);
        for(int blockIdx=0;blockIdx<Ncal;blockIdx++)
        {
            SolveBatch<f_T>(blockIdx, BATCH_SIZE, h_srcidx, h_srcs, h_params, Batchidx, muted, true);
        }

        // dev.launch(FindNextBatch3<f_T>, lpar, stream, d_srcidx.data,d_srcs.data,Batchidx,muted);
        for(int blockIdx=0;blockIdx<Ncal;blockIdx++)
        {
            FindNextBatch3<f_T>(blockIdx, BATCH_SIZE, h_srcidx, h_srcs, Batchidx, muted);
        }

        // lpar = std::make_tuple < dim3, dim3, int > ( (Ncal-1)/BATCH_SIZE + 1,min(BATCH_SIZE,Ncal), 0 );
        // dev.launch(PhysicalTest<f_T>, lpar, stream, d_srcidx.data,d_srcs.data,Batchidx,Ncal,muted);
        for(int blockIdx=0;blockIdx<Ncal;blockIdx++)
        {
            PhysicalTest<f_T>(blockIdx, h_srcidx, h_srcs, Batchidx, Ncal, muted);
        }        

        // lpar = std::make_tuple < dim3, dim3, int > ( Ncal,(BATCH_SIZE), 0 );
        // dev.launch(AreaErrBatchSrc4<f_T>, lpar, stream, d_srcidx.data,d_srcs.data,Batchidx,muted);
        for(int blockIdx=0;blockIdx<Ncal;blockIdx++)
        {
            AreaErrBatchSrc4<f_T>(blockIdx, BATCH_SIZE, h_srcidx, h_srcs, Batchidx, muted);
        }     
        
        // lpar = std::make_tuple < dim3, dim3, int > ( (Ncal-1)/BATCH_SIZE + 1,min(BATCH_SIZE,Ncal), 0 );
        // dev.launch(AreaErrBatchSrcCross<f_T>, lpar, stream, d_srcidx.data,d_srcs.data,d_params.data,Batchidx,Ncal,muted);
        for(int blockIdx=0;blockIdx<Ncal;blockIdx++)
        {
            AreaErrBatchSrcCross<f_T>(blockIdx, h_srcidx, h_srcs, h_params, Batchidx, Ncal, muted);
        }         
        
        if(Batchidx>4)
        // if(false)
        {
            // lpar = std::make_tuple < dim3, dim3, int > ( Ncal,BATCH_SIZE,4*sizeof(f_T)*BATCH_SIZE );
            // dev.launch(SumArea3<f_T>, lpar, stream, d_srcidx.data,d_srcs.data,Batchidx,muted); 
            for(int blockIdx=0;blockIdx<Ncal;blockIdx++)
            {
                SumArea3<f_T>(blockIdx, BATCH_SIZE, h_srcidx, h_srcs, Batchidx, muted, RELTOL);
            }             
        }
        else
        {
            // lpar = std::make_tuple < dim3, dim3, int > ( Ncal,min((BATCH_SIZE*(Batchidx+1)),512),2*sizeof(f_T)*min((BATCH_SIZE*(Batchidx+1)),512) );
            // dev.launch(SumArea0<f_T>, lpar, stream, d_srcidx.data,d_srcs.data,Batchidx,muted);
            for(int blockIdx=0;blockIdx<Ncal;blockIdx++)
            {
                SumArea0<f_T>(blockIdx, min((BATCH_SIZE*(Batchidx+1)),512), h_srcidx, h_srcs, Batchidx, muted, RELTOL);
            }
        }    

        // lpar = std::make_tuple < dim3, dim3, int > ( Ncal,min((BATCH_SIZE*(Batchidx+1)),512),0 );
        // dev.launch(ChangedCheck<f_T>, lpar, stream, d_srcidx.data,d_srcs.data,Batchidx,muted);        

        // printf("%d,%d\n",Batchidx,Ncal);
    }

    // // printf("Ncal here: %d\n",Ncal);
    // int Batchidx = 14;
    // *h_Ncal.data = 0;
    // h_Ncal.cp_to(dev,d_Ncal,stream);
    // lpar = std::make_tuple < dim3, dim3, int > ( (Ncal-1) / BATCH_SIZE+1 , min(BATCH_SIZE,Ncal), 0 );
    // dev.launch(SuccessCheck<f_T>, lpar, stream, d_Ncal.data, d_srcidx.data, d_srcs.data, Ncal, muted);
    // d_Ncal.cp_to(dev,h_Ncal,stream);
    
    // dev.sync_stream( stream );
    // Ncal = *h_Ncal.data;
    // // printf("Ncal here: %d\n",Ncal);


    // // printf("batchidx: %d\n",Batchidx);
    // lpar = std::make_tuple < dim3, dim3, int > ( Ncal,min((BATCH_SIZE*Batchidx),512),(BATCH_SIZE)*(Batchidx)*sizeof(f_T) + (BATCH_SIZE)*sizeof(int) );
    // dev.launch(AdaptiveLocLoop<f_T>, lpar, stream, d_srcidx.data,d_srcs.data,Batchidx,muted);

        
    // lpar = std::make_tuple < dim3, dim3, int > ( Ncal,BATCH_SIZE, 0 );
    // dev.launch(SolveBatch<f_T>, lpar, stream, d_srcidx.data,d_srcs.data,d_params.data,Batchidx,muted,true);
    
    // dev.launch(FindNextBatch3<f_T>, lpar, stream, d_srcidx.data,d_srcs.data,Batchidx,muted);
    
    // lpar = std::make_tuple < dim3, dim3, int > ( (Ncal-1)/BATCH_SIZE + 1,min(BATCH_SIZE,Ncal), 0 );
    // dev.launch(PhysicalTest_debug<f_T>, lpar, stream, d_srcidx.data,d_srcs.data,Batchidx,Ncal,muted);

    // lpar = std::make_tuple < dim3, dim3, int > ( Ncal,(BATCH_SIZE), 0 );
    // dev.launch(AreaErrBatchSrc4<f_T>, lpar, stream, d_srcidx.data,d_srcs.data,Batchidx,muted);
    
    // lpar = std::make_tuple < dim3, dim3, int > ( (Ncal-1)/BATCH_SIZE + 1,min(BATCH_SIZE,Ncal), 0 );
    // dev.launch(AreaErrBatchSrcCross<f_T>, lpar, stream, d_srcidx.data,d_srcs.data,d_params.data,Batchidx,Ncal,muted);
    
    // if(Batchidx>4)
    // // if(true)
    // {
    //     lpar = std::make_tuple < dim3, dim3, int > ( Ncal,BATCH_SIZE,4*sizeof(f_T)*BATCH_SIZE );
    //     dev.launch(SumArea3<f_T>, lpar, stream, d_srcidx.data,d_srcs.data,Batchidx,muted);
        
    //     // lpar = std::make_tuple < dim3, dim3, int > ( Ncal,min((BATCH_SIZE*(Batchidx+1)),512),2*sizeof(f_T)*min((BATCH_SIZE*(Batchidx+1)),512) );
    //     // dev.launch(SumArea0<f_T>, lpar, stream, d_srcidx.data,d_srcs.data,Batchidx,muted);            // attention: sum0 can't bigger than 1024!!!
            
    // }
    // else
    // {
    //     lpar = std::make_tuple < dim3, dim3, int > ( Ncal,min((BATCH_SIZE*(Batchidx+1)),1024),2*sizeof(f_T)*min((BATCH_SIZE*(Batchidx+1)),1024) );
    //     dev.launch(SumArea0<f_T>, lpar, stream, d_srcidx.data,d_srcs.data,Batchidx,muted);            
    // } 

 
}
template void Twinkle<float> ::solve();
template void Twinkle<double>::solve();

template<isFloating f_T>
void Twinkle<f_T>::cp_back()
{
    // d_Nfail.cp_to(dev,h_Nfail,stream);
    // // d_params_back.cp_to(dev,h_params_back,stream);
    // d_srcs.cp_to(dev,h_srcs,stream);
    // // d_crit.cp_to(dev,h_crit,stream);
    // // d_caus.cp_to(dev,h_caus,stream);
    // // d_connect.cp_to(dev,h_connect,stream);
    // // if(critcaus)
    // // {
    //     d_CC.cp_to(dev,h_CC,stream);
    // // }

    // dev.sync_stream( stream );   
    if(!muted)
    {
        printf("total failed points: %d\n",*h_Nfail);
    }
}
template void Twinkle<float> ::cp_back();
template void Twinkle<double>::cp_back();

template<isFloating f_T>
void Twinkle<f_T>::writeto(std::string path, std::string suffix)
{
    std::ofstream ofs;     // 创建流对象
    ofs.precision(16);

    std::string name = "src_ext";
    std::string txt = ".txt";
    name = path + name + suffix + txt;

    ofs.open(name, std::ios::out);      // 以 write 模式打开 test.txt
    for(int iii=0;iii<Nsrcs;iii++)
    {

        ofs<<h_srcs[iii].shape.loc_centre.re<<",";
        ofs<<h_srcs[iii].shape.loc_centre.im<<",";
        ofs<<h_srcs[iii].shape.rho<<",";
        ofs<<h_srcs[iii].Mag<<",";
        // ofs<<h_srcs[i].Mag_p<<",";
        ofs<<h_srcs[iii].Err<<",";
        ofs<<h_srcs[iii].Area<<",";
        ofs<<h_srcs[iii].Err_A<<",";
        ofs<<h_srcs[iii].SolveSucceed<<",";
        #if DetailedRecord
            ofs<<h_srcs[iii].points_used<<",";
        #endif
        ofs<<h_srcs[iii].Mag_P<<",";
        ofs<<h_srcs[iii].Err_Mag_P<<",";
        ofs<<h_srcs[iii].Ncross<<",";
        ofs<<std::endl;
        
    }
    ofs.close();
    // dev.sync_stream( stream );
    if(!muted)
    {
        std::cout<<"print file: "<<name<<std::endl;        
    }

    // if(critcaus){
    //     name = "critcaus";
    //     name = path + name + suffix + txt;
    //     // cout<<h_params[0].q<<endl;
    //     // cout<<h_params[0].s<<endl;
    //     ofs.open(name, std::ios::out);      // 以 write 模式打开 test.txt
    //     for(int jjj=0;jjj<NBODY*NBODY;jjj++)
    //     {    
    //         for(int iii=0;iii<NCRIT;iii++)
    //         {
    //         // ofs<<h_cc_bool[iii]<<",";
    //         // ofs<<h_crit[iii].re<<",";
    //         // ofs<<h_crit[iii].im<<",";
    //         // ofs<<h_caus[iii].re<<",";
    //         // ofs<<h_caus[iii].im<<",";    
    //         // ofs<<h_connect[iii];
    //             ofs<<iii<<",";
    //             ofs<<h_CC[iii].psi<<",";
    //             ofs<<h_CC[iii].crit[jjj].re<<",",
    //             ofs<<h_CC[iii].crit[jjj].im<<",",
    //             ofs<<h_CC[iii].caus[jjj].re<<",",
    //             ofs<<h_CC[iii].caus[jjj].im<<",",
    //             ofs<<h_CC[iii].next_idx<<",";
    //             ofs<<h_CC[iii].next_j[jjj]<<",";
    //             // ofs<<h_CC[iii].dist2next<<",";
    //             ofs<<h_CC[iii].skip;
    //             ofs<<std::endl;            
    //         }

    //     }
    //     ofs.close();
    //     if(!muted)
    //     {
    //         std::cout<<"print file: "<<name<<std::endl;
    //     }
    // }

    for(int iii=0;iii<Nsrcs;iii++)
    {
        // if(iii == 663 || iii == 701 || iii==576)
        // if(iii>=505 && iii<= 515)
        // if(iii == 0)
        if(false)
        {    

            name= "img_pt" + std::to_string(iii);
            name = path + name + txt;

            ofs.open(name, std::ios::out);      // 以 write 模式打开 test.txt


                // 第i个src，第j个root，实部，虚部
                for (int j = 0; j < LENGTH-1; j++) {
                    for (int i = 0; i < NPOINTS; i++) {
                        ofs << i << ",";
                        ofs << j << ",";
                        ofs<<h_srcs[iii].margin_pts[i].images[j].position.re<<",";
                        ofs<<h_srcs[iii].margin_pts[i].images[j].position.im<<",";
                        ofs<<h_srcs[iii].margin_pts[i].images[j].physical << ",";
                        ofs<<h_srcs[iii].margin_pts[i].images[j].parity << "," ;
                        ofs<<h_srcs[iii].margin_pts[i].deltaS[j]<<",";
                        ofs<<h_srcs[iii].margin_pts[i].next_idx[j]<<",";
                        ofs<<h_srcs[iii].margin_pts[i].next_j[j]<<",";
                        ofs<<h_srcs[iii].margin_pts[i].deltaS_Err[j]<<",";
                        ofs<<h_srcs[iii].margin_pts[i].deltaS_new[j]<<",";
                        ofs<<h_srcs[iii].margin_pts[i].Err_new[j]<<",";
                        // if constexpr (DetailRecord)
                        #if DetailedRecord
                        // {
                            ofs<<h_srcs[iii].margin_pts[i].jacobi[j]<<",";
                            ofs<<h_srcs[iii].margin_pts[i].E1[j]<<",";
                            ofs<<h_srcs[iii].margin_pts[i].E2[j]<<",";
                            ofs<<h_srcs[iii].margin_pts[i].E3[j]<<",";
                            ofs<<h_srcs[iii].margin_pts[i].deltaS_p[j]<<",";
                            ofs<<h_srcs[iii].margin_pts[i].phys_err[j]<<",";
                            
                        // }
                        #endif
                        // ofs<<h_srcs[iii].margin_pts[i].images[j].deltaz.re<<",";
                        // ofs<<h_srcs[iii].margin_pts[i].images[j].deltaz.im<<",";
                        ofs << std::endl;
                    }        
                }
            ofs.close();
            std::cout<<"print file: "<<name<<std::endl;


            name= "src_pt" + std::to_string(iii);
            name = path + name + txt;
            ofs.open(name, std::ios::out);      // 以 write 模式打开 test.txt
            
            for(int i = 0;i<NPOINTS;i++)
            {
                ofs<<h_srcs[iii].margin_pts[i].position.re<<",";
                ofs<<h_srcs[iii].margin_pts[i].position.im<<",";
                ofs<<h_srcs[iii].margin_pts[i].next_src_idx<<",";
                ofs<<h_srcs[iii].margin_pts[i].Q<<",";
                // ofs<<h_srcs[iii].margin_pts[i].quantile;
                ofs<<h_srcs[iii].margin_pts[i].prev_src_idx<<",";
                ofs<<h_srcs[iii].margin_pts[i].skip<<",";
                ofs<<h_srcs[iii].margin_pts[i].Nphys<<",";
                if(critcaus)
                {
                    ofs<<h_srcs[iii].margin_pts[i].NcrossCaus<<",";
                }
                else{
                    ofs<<0<<",";
                }
                ofs<<h_srcs[iii].margin_pts[i].changed<<",";
                ofs<<h_srcs[iii].margin_pts[i].error_interval<<",";
                // if(i<BATCH_SIZE)
                // {
                //     cout<<h_srcs[iii].margin_pts[i].Nphys<<endl;
                // }
                


                ofs<<std::endl;
            }
            ofs.close();
            std::cout<<"print file: "<<name<<std::endl;

        }

    }


}
template void Twinkle<float> ::writeto(std::string path, std::string suffix);
template void Twinkle<double>::writeto(std::string path, std::string suffix);



template<isFloating f_T>
_Global void DataZip(int srcidx, src_ext_t<f_T>* srcs, f_T* Mags, f_T* Errs)
{
    // int srcidx = threadIdx.x + blockIdx.x*blockDim.x;
    if(srcidx<NSRCS)
    {
        Mags[srcidx] = srcs[srcidx].Mag;
        Errs[srcidx] = srcs[srcidx].Err;        
    }
}
_Global void DataZip(src_ext_t<float>* srcs, float* Mags, float* Errs);
_Global void DataZip(src_ext_t<double>* srcs, double* Mags, double* Errs);

template<isFloating f_T>
void Twinkle<f_T>::cp_back_ME()
{
    // auto lpar = std::make_tuple < dim3, dim3, int > ( (Nsrcs-1)/BATCH_SIZE + 1,min(BATCH_SIZE,Nsrcs), 0 );
    // dev.launch(DataZip<f_T>, lpar, stream, d_srcs.data, d_Mag.data, d_Err.data);

    // d_Mag.cp_to(dev,h_Mag,stream);
    // d_Err.cp_to(dev,h_Err,stream);    

    // dev.sync_stream( stream );   
}

// Magnification and Error
template<isFloating f_T>
void Twinkle<f_T>::writeto_ME(std::string path, std::string suffix)
{
    std::ofstream ofs;     // 创建流对象
    ofs.precision(16);

    std::string name = "src_ext";
    std::string txt = ".txt";
    name = path + name + suffix + txt;

    ofs.open(name, std::ios::out);      // 以 write 模式打开 test.txt
    for(int iii=0;iii<Nsrcs;iii++)
    {

        ofs<<h_Mag[iii]<<",";
        // ofs<<h_srcs[i].Mag_p<<",";
        ofs<<h_Err[iii];
        ofs<<std::endl;
        
    }
    ofs.close();
    // dev.sync_stream( stream );
    if(!muted)
    {
        std::cout<<"print file: "<<name<<std::endl;        
    }
}
template void Twinkle<float> ::cp_back_ME();
template void Twinkle<double>::cp_back_ME();


// template<isFloating f_T>
// void Twinkle<f_T>::writeto_caus(std::string path, std::string suffix)
// {
//     d_CC.cp_to(dev,h_CC,stream);
//     dev.sync_stream( stream );   


//     std::ofstream ofs;     // 创建流对象
//     ofs.precision(16);

//     std::string name = "critcaus";
//     std::string txt = ".txt";
//     name = path + name + suffix + txt;

//     // cout<<h_params[0].q<<endl;
//     // cout<<h_params[0].s<<endl;
//     ofs.open(name, std::ios::out);      // 以 write 模式打开 test.txt
//     for(int jjj=0;jjj<NBODY*NBODY;jjj++)
//     {    
//         for(int iii=0;iii<NCRIT;iii++)
//         {
//         // ofs<<h_cc_bool[iii]<<",";
//         // ofs<<h_crit[iii].re<<",";
//         // ofs<<h_crit[iii].im<<",";
//         // ofs<<h_caus[iii].re<<",";
//         // ofs<<h_caus[iii].im<<",";    
//         // ofs<<h_connect[iii];
//             ofs<<iii<<",";
//             ofs<<h_CC[iii].psi<<",";
//             ofs<<h_CC[iii].crit[jjj].re<<",",
//             ofs<<h_CC[iii].crit[jjj].im<<",",
//             ofs<<h_CC[iii].caus[jjj].re<<",",
//             ofs<<h_CC[iii].caus[jjj].im<<",",
//             ofs<<h_CC[iii].next_idx<<",";
//             ofs<<h_CC[iii].next_j[jjj]<<",";
//             // ofs<<h_CC[iii].dist2next<<",";
//             ofs<<h_CC[iii].skip;
//             ofs<<std::endl;            
//         }

//     }
//     ofs.close();
//     if(!muted)
//     {
//         std::cout<<"print file: "<<name<<std::endl;
//     }
// }



// template <isFloating f_T>
// void Twinkle<f_T>::solve_caus()
// {
//     // 这一步仅用于定义 lpar
//     auto lpar = std::make_tuple
//             < dim3, dim3, int > ( Nsrcs, BATCH_SIZE, 0 );

//     for(int batch_c=0;batch_c<NCRIT/BATCH_SIZE;batch_c++)
//     // for(int batch_c=0;batch_c<1;batch_c++)
//     {
//         // printf("%d\n",batch_c);
//         lpar = std::make_tuple < dim3, dim3, int > ( dim3(1,1,1), max(batch_c,1)*BATCH_SIZE, (BATCH_SIZE)*(batch_c)*sizeof(f_T) + (BATCH_SIZE)*sizeof(int) );
//         dev.launch(psiLoc<f_T>, lpar, stream, d_CC.data,batch_c);

//         lpar = std::make_tuple < dim3, dim3, int > ( dim3(1,1,1), BATCH_SIZE, 0 );
//         dev.launch(SolveCritCaus<f_T>, lpar, stream, d_params.data, d_CC.data, batch_c);

//         lpar = std::make_tuple < dim3, dim3, int > ( dim3(1,1,(batch_c+1)), BATCH_SIZE, 0 );
//         dev.launch(CNext2<f_T>, lpar, stream, d_CC.data,0);            // 所有点放进去连接
//     }
// }
// template void Twinkle<float> ::solve_caus();
// template void Twinkle<double>::solve_caus();



// template<isFloating f_T>
// void Twinkle<f_T>::cp_back_caus()
// {
//     d_CC.cp_to(dev,h_CC,stream);
//     dev.sync_stream( stream );   
// }

template<isFloating f_T>
f_T Twinkle<f_T>::Mag(f_T s, f_T q, f_T y1, f_T y2, f_T Rs)
{
    this->Nsrcs = 1;
    h_srcs[0].shape.rho = Rs;
    h_srcs[0].shape.loc_centre.re = y1;
    h_srcs[0].shape.loc_centre.im = y2;
    h_srcs[0].src_area = Area_src(h_srcs[0].shape);
    h_params[0].s = s;
    h_params[0].q = q;

    this->solve();

    // printf("Mag = %.16f\n",h_srcs[0].Mag);
    return h_srcs[0].Mag;
}
template float Twinkle<float>::Mag(float s, float q, float y1, float y2, float Rs);
template double Twinkle<double>::Mag(double s, double q, double y1, double y2, double Rs);