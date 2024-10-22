// #include<iostream>

// #include"calculation/init.h"

// #define cout std::cout
// #define endl std::endl

// #include <fstream>


// int main()
// {
//     Twinkle<double> LittleStar;
//     LittleStar.malloc(0);
    
//     LittleStar.RELTOL = 1e-3;

//     ////// INPUT 
//     // input observation time
//     double time[NSRCS];
//     for(int iii=0;iii<NSRCS;iii++)
//     {
//         time[iii] = double(iii) / double(NSRCS) * 2 - 1;
//     }


//     // input parameters
//     src_params_t<double> params;
//     params.t_0 = 0;
//     params.t_E = 0.218999842;
//     params.u_0 = 0.01;
//     params.alpha = 0.53;
//     params.shape.rho = 1.939366667;
//     params.q = 0.125;
//     params.s = 0.437613333;
//     // LittleStar.set_params(&params);

//     LittleStar.set_path(time,&params);


//     ////// SOLVE
//     LittleStar.solve();


//     ////// COPY BACK ALL DATA
//     // LittleStar.cp_back();
//     LittleStar.writeto("./");


//     // // single point
//     // double Magnification;
//     // Magnification = LittleStar.Mag(0.437613333, 0.125,-0.2411797387902984,-0.1529013886069155, 1.939366667);

    
    

//     // FREE
//     LittleStar.free();   


//     return 0;
// }