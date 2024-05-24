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
    
    double denom = m*Sin(del)*Sin(del)+M;
    dydt[0] = (l*y[2]-y[3]*L*Cos(del))/(l*L*L*denom);
    dydt[1] = ((M+m)*l*y[3]-l*m*y[2]*Cos(del))/(l*l*L*m*denom);
    
    double T1 = (y[2]*y[3]*Sin(del))/(L*l*denom);
    double T2 = ((L*L*y[3]*y[3]*(M+m)+m*l*l*y[2]*y[2]-2*m*l*L*y[2]*y[3]*Cos(del))*Sin(2*del))/(2*l*l*L*L*denom*denom);
    
    dydt[2] = T2 - T1 - L*g*(M+m)*Sin(y[0]);
    dydt[3] = T1 - T2 - g*l*m*Sin(y[1]);
    
    return GSL_SUCCESS;
}

int main(void)
{
    //M,L,m,l,g
    double param[5] = {1.0, 1.0, 1.0, 1.0, 1.0};
    double M = param[0], L = param[1], m = param[2], l = param[3], g = param[4];
    
    gsl_odeiv2_system system = {Func, NULL, 4, param};
    gsl_odeiv2_driver* driver = 
            gsl_odeiv2_driver_alloc_y_new(&system, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);
    double t = 0.0, prev_t = 0.0;
    const double t_step_size = 0.01;
    
    //Initial conditions.
    double y[4] = {M_PI/6, 0.0, 0.0, 0.0}, prev_y[4], del = y[0] - y[1];
    const double E = -2.0;
    
    //Calculate new pφ based on choice of E.
    y[3] = (m*l*y[2]*Cos(y[0])+l*sqrt(m*m*y[2]*y[2]*Cos(y[0])*Cos(y[0])-m*(M+m)*
    (y[2]*y[2]-2*L*L*(M+m*Sin(y[0])*Sin(y[0]))*(E+(M+m)*g*L*Cos(y[0])+m*g*l))))/((M+m)*L);
    
    //printf("%.5e \n", y[3]);
    
    for(int i = 1; i <= 100000; ++i)
    {
        double ti = i*t_step_size; //Current time.
        int status = gsl_odeiv2_driver_apply(driver, &t, ti, y);
        if(status != GSL_SUCCESS)
        {
            fprintf(stderr, "Error with return: %d\n", status);
            break;
        }
        
        double Et = (m*l*l*y[2]*y[2]+(m+M)*L*L*y[3]*y[3] - 2*m*l*L*y[2]*y[3]*Cos(y[0]-y[1]))/(2*m*l*l*L*L*(M+m*Sin(y[0]-y[1])*Sin(y[0]-y[1]))) 
                -(m+M)*g*L*Cos(y[0]) - m*g*l*Cos(y[1]);
        printf("%.5e\n", Et);
        
        /*Here we have time, angle of the first rod, momentum of the first rod
          angle of the second rod and momentum of the second rod.
        printf("%.5e %.5e %.5e %.5e %.5e\n", t, y[0], y[1], y[2], y[3]); */
    }
    
    gsl_odeiv2_driver_free(driver);
    return 0;
}