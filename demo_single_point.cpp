#include<iostream>

#include"calculation/init.h"


int main()
{
    Twinkle<double> LittleStar;
    LittleStar.malloc(0);

    LittleStar.RELTOL = 1e-2;


    // single point
    double Magnification;
    // LittleStar.Mag(s, q, y1, y2, Rs);
    Magnification = LittleStar.Mag(0.437613333, 0.125,-0.2411797387902984,-0.1529013886069155, 1.939366667e-2);
    printf("Magnification = %.12f\n",Magnification);

    // FREE
    LittleStar.free();   


    return 0;
}