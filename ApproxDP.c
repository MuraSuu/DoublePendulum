//GSL
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>

//C
#include <stdio.h>
#include <math.h>

int Func(double t, const double y[], double dydt[], void* params)
{
    (void)t;
    double* temp = (double*)params;
    double M = temp[0], L = temp[1], m = temp[2], l = temp[3], g = temp[4];
    //double del = y[0] - y[1];
    
    dydt[0] = -(-((1.0/2.0)/M + m/pow(M, 2))*pow(-y[1] + y[0], 2) + 1.0/M)*y[3]/(L*l) + (-m*(-y[1] + y[0])/M + 1.0/M)*y[2]/pow(L, 2);
    dydt[1] = (M + m)*(-m*(-y[1] + y[0])/M + 1.0/M)*y[3]/(pow(l, 2)*m) - (-((1.0/2.0)/M + m/pow(M, 2))*pow(-y[1] + y[0], 2) + 1.0/M)*y[2]/(L*l);
    dydt[2] = -(L*g*(-M - m)*((1.0/6.0)*pow(y[0], 3) - y[0]) + ((1.0/2.0)/M + m/pow(M, 2))*(-2*y[1] + 2*y[0])*y[3]*y[2]/(L*l) - 1.0/2.0*(pow(L, 2)*(M + m)*pow(y[3], 2) + pow(l, 2)*m*pow(y[2], 2))/(pow(L, 2)*M*pow(l, 2)));
    dydt[3] = -(-g*l*m*((1.0/6.0)*pow(y[1], 3) - y[1]) + ((1.0/2.0)/M + m/pow(M, 2))*(2*y[1] - 2*y[0])*y[3]*y[2]/(L*l) + (1.0/2.0)*(pow(L, 2)*(M + m)*pow(y[3], 2) + pow(l, 2)*m*pow(y[2], 2))/(pow(L, 2)*M*pow(l, 2)));
    
    return GSL_SUCCESS;
}

int main(void)
{
    double param[5] = {1.0, 1.0, 0.5, 1.0, 2.0};
    //double M = param[0], L = param[1], m = param[2], l = param[3], g = param[4];
    
    gsl_odeiv2_system system = {Func, NULL, 4, param};
    gsl_odeiv2_driver* driver = 
            gsl_odeiv2_driver_alloc_y_new(&system, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);
    double t = 0.0;
    const double t_step_size = 0.01;
    
    double y[4] = {M_PI/10, M_PI/100, 0.0, 0.0};
    
    for(int i = 1; i <= 900; ++i)
    {
        double ti = i*t_step_size; //Current time.
        int status = gsl_odeiv2_driver_apply(driver, &t, ti, y);
        if(status != GSL_SUCCESS)
        {
            fprintf(stderr, "Error with return: %d\n", status);
            break;
        }
        
        printf("%.5e %.5e %.5e %.5e %.5e\n", t, y[0], y[1], y[2], y[3]);
    }
    
    gsl_odeiv2_driver_free(driver);
    return 0;
}