#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_errno.h>
//#include <gsl/gsl_matrix.h>
#include <stdio.h>

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
    double* temp = (double*)params;
    double M = temp[0], L = temp[1], m = temp[2], l = temp[3], g = temp[4];
    double del = y[0] - y[1];
    
    dydt[0] = (l*y[2]-L*y[3]*Cos(del))/(l*L*L*(M+m*Sin(del)*Sin(del)));
    dydt[1] = ((m+M)*L*y[3] - m*l*y[2]*Cos(del))/(m*L*l*l*(M+m*Sin(del)*Sin(del)));
    
    double T1 = (y[2]*y[3]*Sin(del))/(l*L*(M+m*Sin(del)*Sin(del)));
    double T2 = ((m*l*l*y[2]*y[2]+(M+m)*L*L*y[3]*y[3]-2*m*l*L*y[2]*y[3]*Cos(del))*Sin(2*del))/
                (2*l*l*L*L*(M+m*Sin(del)*Sin(del))*(M+m*Sin(del)*Sin(del)));
    
    dydt[2] = T2 - T1 - (M+m)*g*L*Sin(del);
    dydt[3] = T1 - T2 - m*g*l*Sin(y[1]);
    
    return GSL_SUCCESS;
}

int main(void)
{
    //M,L,m,l,g
    double param[5] = {1.0, 1.0, 0.01, 0.01, 1.0};
    gsl_odeiv2_system system = {Func, NULL, 4, param};
    gsl_odeiv2_driver* driver = 
            gsl_odeiv2_driver_alloc_y_new(&system, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);
    double t = 0.0, t1 = 50.0;
    
    double y[4] = {3.14159/6, 0.0, 0.0, 0.0};
    
    for(int i = 1; i <= 360; ++i) //180 seconds (3 minutes).
    {
        double ti = i*t1/100.0;
        int status = gsl_odeiv2_driver_apply(driver, &t, ti, y);
        if(status != GSL_SUCCESS)
        {
            fprintf(stderr, "Error with return: %d\n", status);
            break;
        }
        
        /*Here we have time, angle of the first rod, momentum of the first rod
          angle of the second rod and momentum of the second rod.*/
        printf("%.5e %.5e %.5e %.5e %.5e\n", t, y[0], y[1], y[2], y[3]);
    }
    
    gsl_odeiv2_driver_free(driver);
    return 0;
}