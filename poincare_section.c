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
    double M = param[0], L = param[1], m = param[2], l = param[3], g = param[4];
    gsl_odeiv2_system system = {Func, NULL, 4, param};
    gsl_odeiv2_driver* driver = 
            gsl_odeiv2_driver_alloc_y_new(&system, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);
    double t = 0.0, t1 = 1.0;
    
    double y[4] = {3.14159/6, 3.14159/2, 0.0, 0.0};
    
    double prev_y[4];
    double section_width = 0.1; //For phi = 0.
    
    for(int i = 1; i <= 1000000; ++i)
    {
        for(int j = 0; j < 4; ++j)
        {
            y[j] = fmod(y[j], 2*M_PI);
            prev_y[j] = y[j];
        }
        
        double ti = i*t1/100.0;
        int status = gsl_odeiv2_driver_apply(driver, &t, ti, y);
        if(status != GSL_SUCCESS)
        {
            fprintf(stderr, "Error with return: %d\n", status);
            break;
        }
        
        double del = y[0] - y[1];
        double theta_dot = (l*y[2]-L*y[3]*Cos(del))/(l*L*L*(M+m*Sin(del)*Sin(del)));
        
        //Check if crossed the section.
        if((y[1] - section_width)*(prev_y[1] - section_width) < 0 && theta_dot > 0.1)
        {
            double data[4];
            for(int j = 0; j < 4; ++j)
                data[j] = prev_y[j] + (y[j]-prev_y[j])*(section_width - prev_y[1])/(y[1]-prev_y[1]);
            printf("%.5e %.5e\n", data[0], data[2]);
        }
    }
    
    gsl_odeiv2_driver_free(driver);
    return 0;
}