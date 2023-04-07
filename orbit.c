//GSL
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_errno.h>

//C
#include <stdio.h>
#include <math.h>

#define Sin gsl_sf_sin
#define Cos gsl_sf_cos

/*
y0 = theta
y1 = phi
y2 = p_theta
y3 = p_phi
*/

int Func(double t, const double y[], double dydt[], void* params)
{
    (void)t;
    double* temp = (double*)params;
    double M = temp[0], L = temp[1], m = temp[2], l = temp[3], g = temp[4];
    double del = y[0] - y[1];
    
    dydt[0] = (l*y[2]-L*y[3]*Cos(del))/(l*L*L*(M+m*Sin(del)*Sin(del)));
    dydt[1] = ((m+M)*L*y[3] - m*l*y[2]*Cos(del))/(m*L*l*l*(M+m*Sin(del)*Sin(del)));
    
    double T1 = (y[2]*y[3]*Sin(del))/(l*L*(M+m*Sin(del)*Sin(del)));
    double T2 = ((m*l*l*y[2]*y[2]+(M+m)*L*L*y[3]*y[3]-2*m*l*L*y[2]*y[3]*Cos(del))*Sin(2*del))/
                (2*l*l*L*L*(M+m*Sin(del)*Sin(del))*(M+m*Sin(del)*Sin(del)));
    //The error was right here. sin of theta, not of theta-phi.
    dydt[2] = T2 - T1 - (M+m)*g*L*Sin(y[0]);
    dydt[3] = T1 - T2 - m*g*l*Sin(y[1]);
    
    return GSL_SUCCESS;
}

int main(void)
{
    //M,L,m,l,g
    double param[5] = {1.0, 1.0, 0.1, 1.0, 1.0};
    double L = param[1], l = param[3];
    gsl_odeiv2_system system = {Func, NULL, 4, param};
    gsl_odeiv2_driver* driver = 
            gsl_odeiv2_driver_alloc_y_new(&system, gsl_odeiv2_step_rkf45, 1e-8, 1e-8, 0.0);
    double t = 0.0, t1 = 1.0;
    
    double y[4] = {3.14159/4, 0.0, 0.0, 0.5};
    
    for(int i = 1; i <= 1000000; ++i) //180 seconds (3 minutes).
    {
        y[0] = fmod(y[0], 2*M_PI); //Theta.
        y[1] = fmod(y[1], 2*M_PI); //Phi.
    
        double ti = i*t1/100.0;
        int status = gsl_odeiv2_driver_apply(driver, &t, ti, y);
        if(status != GSL_SUCCESS)
        {
            fprintf(stderr, "Error with return: %d\n", status);
            break;
        }
        
        double x1 =  L*Sin(y[0]);
        double y1 = -L*Cos(y[0]);
        double x2 =  L*Sin(y[0]) + l*Sin(y[1]);
        double y2 = -L*Cos(y[0]) - l*Cos(y[1]);
        
        /*Here we have time, angle of the first rod, momentum of the first rod
          angle of the second rod and momentum of the second rod.*/
        printf("%.5e %.5e %.5e %.5e\n", x1, y1, x2, y2);
    }
    
    gsl_odeiv2_driver_free(driver);
    return 0;
}