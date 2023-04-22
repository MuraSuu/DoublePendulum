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
    dydt[2] = T2 - T1 - (M+m)*g*L*Sin(y[0]); //The error was right here. sin(θ), not of θ-φ.
    dydt[3] = T1 - T2 - m*g*l*Sin(y[1]);
    
    return GSL_SUCCESS;
}

//Helper function.
void Swap(double* a, double* b)
{
    if(*a > *b)
    {
        double c = *a;
        *a = *b;
        *b = c;
    }
}

int main(void)
{
    //M,L,m,l,g
    double param[5] = {1.0, 1.0, 1.0, 1.0, 1.0};
    double M = param[0], L = param[1], m = param[2], l = param[3], g = param[4];
    
    gsl_odeiv2_system system = {Func, NULL, 4, param};
    gsl_odeiv2_driver* driver = 
            gsl_odeiv2_driver_alloc_y_new(&system, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);
    double t = 0.0, t1 = 1.0, prev_t = 0.0;
    
    //Initial conditions.
    double y[4] = {M_PI/6.0, 0.0, 0.0, 0.0};
    double prev_y[4];
    double E = 5.0, del = y[0] - y[1];
    
    //Calculate new pφ based on choice of E.
    y[3] = (m*l*y[2]*Cos(y[0])+l*sqrt(m*m*y[2]*y[2]*Cos(y[0])*Cos(y[0])-m*(M+m)*
    (y[2]*y[2]-2*L*L*(M+m*Sin(y[0])*Sin(y[0]))*(E+(M+m)*g*L*Cos(y[0])+m*g*l))))/((M+m)*L);
    
    //Solve eq.
    for(int i = 1; i <= 1000000; ++i)
    {
        for(int i = 0; i < 4; ++i) prev_y[i] = y[i];
        prev_t = t;
    
        double ti = i*t1/100.0; //Current time.
        int status = gsl_odeiv2_driver_apply(driver, &t, ti, y);
        if(status != GSL_SUCCESS)
        {
            fprintf(stderr, "Error with return: %d\n", status);
            break;
        }
        
        //x (mod 2pi).
        y[0] = fmod(y[0], 2*M_PI); //θ.
        y[1] = fmod(y[1], 2*M_PI); //φ
        
        double theta_dot = (l*y[2]-L*y[3]*Cos(del))/(l*L*L*(M+m*Sin(del)*Sin(del)));
        
        //Check if crossed the section.
        if(y[1]*prev_y[1] < 0 && theta_dot > 0)
        {
            double delta_phi = y[1] - prev_y[1], delta_t = t - prev_t;
            double delH_delp = delta_phi/delta_t; // This will be very close to the exact time.
            
            double t_new = ti - delta_t + prev_y[1]/delH_delp;
            Swap(&prev_t, &t_new); //Necessary in order for the limits of integration to make sense.
            
            status = gsl_odeiv2_driver_apply(driver, &prev_t, t_new, prev_y);
            if(status != GSL_SUCCESS)
            {
                fprintf(stderr, "Error with return: %d\n", status);
                break;
            }else{
                printf("%.5e %.5e\n", y[0], y[2]);
            }
        }
        
        if(status != GSL_SUCCESS) break;
    }
    printf("End\n");
    
    gsl_odeiv2_driver_free(driver);
    return 0;
}