#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <cstdio>

#define Sin gsl_sf_sin
#define Cos gsl_sf_cos

//We have the following parameters in this exact order:
//Mass and length of the first part and mass and length of the second part.
//and as the last parameter we have gravity.
int Func(double t, const double y[], double dydt[], void* params)
{
    double* temp = static_cast<double*>(params);
    double M = temp[0], L = temp[1], m = temp[2], l = temp[3], g = temp[4];
    double del = y[0] - y[2];
    dydt[0] = y[1];
    double denom1 = L*(m+M) - m*L*Cos(del)*Cos(del);
    dydt[1] = (-m*Cos(del)*L*y[1]*y[1]*Sin(del)+m*Cos(del)*g*Sin(y[2])
              -m*l*y[3]*y[3]*Sin(del)-(m+M)*Sin(y[0]))/denom1;
    dydt[2] = y[3];
    double denom2 = M*l+m*l*Sin(del)*Sin(del);
    dydt[3] = (m+M)*(L*y[1]*y[1]*Sin(del)+(m*l/(m+M))*y[3]*y[3]*Sin(del)*Cos(del)
              -g*Sin(y[2])+g*Sin(y[0])*Cos(del))/denom2;
    
    return GSL_SUCCESS;
}

int main()
{
    //M, L, m, l, g
    double param[5] = {2.0, 2.0, 0.05, 0.01, 9.8};
    gsl_odeiv2_system system{Func, nullptr, 4, param};
    gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new(&system, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);
    
    double t = 0.0, t1 = 50.0;
    //double y[4] = {0.0, 0.0, 0.0, 0.0};
    //double y[4] = {0.0, 0.0, 3.14159/6, 0.0};
    //double y[4] = {3.14159/4, 0.0, 3.14159/4, 0.0};
    //double y[4] = {0.0, 0.0, 3.14159/2, 0.0};
    //double y[4] = {3.14159, 0.0, 3.14159, 0.0};
    double y[4] = {3.14159/2, 0.0, 0.0, 0.0};
    
    for(int i = 1; i <= 360; ++i) //Run from t=0 to t=6 minutes(360 seconds).
    {
        double ti = i*t1/100.0;
        int status = gsl_odeiv2_driver_apply(driver, &t, ti, y);
        if(status != GSL_SUCCESS)
        {
            std::fprintf(stderr, "Error with return: %d\n", status);
            break;
        }
        /*Here we have time, angle of the first rod, angular velocity of the first rod
          angle of the second rod and angular velocity of the second rod.*/
        //if(i>=280)
        std::printf("%.5e %.5e %.5e %.5e %.5e\n", t, y[0], y[1], y[2], y[3]);
    }
    
    gsl_odeiv2_driver_free(driver);
}