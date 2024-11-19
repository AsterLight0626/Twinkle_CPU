#include<iostream>

#include"calculation/init.h"


int main()
{
    Twinkle<double> LittleStar;
    LittleStar.malloc(0);

    LittleStar.RELTOL = 1e-5;


    // single point
    double Magnification;
    // LittleStar.Mag(s, q, y1, y2, Rs);
    Magnification = LittleStar.Mag(0.5, 1e-6,-1.499934,0.0035, 1e-4);
    printf("Magnification = %.12f\n",Magnification);



    LittleStar.writeto("./");

    // FREE
    LittleStar.free();   


    return 0;
}
