//GSL
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_errno.h>

//C
#include <stdio.h>
#include <string.h>
#include <math.h>

#define Sin gsl_sf_sin
#define Cos gsl_sf_cos

/*
y0 = θ
y1 = φ
y2 = p_θ
y3 = p_φ
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

int main(void)
{
    //M,L,m,l,g
    double param[5] = {3.0, 2.0, 1.0, 1.0, 1.0};
    double M = param[0], L = param[1], m = param[2], l = param[3], g = param[4];
    
    gsl_odeiv2_system system = {Func, NULL, 4, param};
    gsl_odeiv2_driver* driver = 
            gsl_odeiv2_driver_alloc_y_new(&system, gsl_odeiv2_step_rkf45, 1e-8, 1e-8, 0.0);
    double t = 0.0, prev_t = 0.0;
    const double t_step_size = 0.01;
    
    //Initial conditions.
    double y[4] = {M_PI/6.0, 0.0, 0.0, 0.0}, prev_y[4], del = y[0] - y[1];
    const double E = -8.85;
    
    //Calculate new pφ based on choice of E.
    y[3] = (m*l*y[2]*Cos(y[0])+l*sqrt(m*m*y[2]*y[2]*Cos(y[0])*Cos(y[0])-m*(M+m)*
    (y[2]*y[2]-2*L*L*(M+m*Sin(y[0])*Sin(y[0]))*(E+(M+m)*g*L*Cos(y[0])+m*g*l))))/((M+m)*L);
    
    //Solve eq.
    for(int i = 1; i <= 500000; ++i)
    {
        for(int k = 0; k < 4; ++k){ prev_y[k] = y[k]; }
        prev_t = t;
        
        double ti = i*t_step_size; //Current time.
        int status = gsl_odeiv2_driver_apply(driver, &t, ti, y);
        if(status != GSL_SUCCESS)
        {
            fprintf(stderr, "Error with return: %d\n", status);
            break;
        }
        
        double theta_dot = (l*y[2]-L*y[3]*Cos(del))/(l*L*L*(M+m*Sin(del)*Sin(del)));
        
        //Check if crossed the section.
        if(y[1]*prev_y[1] < 0 && theta_dot > 0)
        {
            double ya[4], yb[4], ta = prev_t, tb = t;
            memcpy(ya, prev_y, sizeof(prev_y));
            memcpy(yb, y, sizeof(y));
            
            double temp_y[4]; memcpy(temp_y, ya, sizeof(ya));
            double temp_t = ta;
            for(int i = 1; i <= 3; ++i)
            {
                double phi_dot = ((m+M)*L*temp_y[3] - m*l*temp_y[2]*Cos(temp_y[0]-temp_y[1]))/(m*L*l*l*(M+m*Sin(temp_y[0]-temp_y[1])*Sin(temp_y[0]-temp_y[1])));
                
                double t_next = temp_t - temp_y[1]/phi_dot;
                
                //If a swap is necessary it will also be necessary to swap the data.
                //For example, if y_temp(yb) is far ahead of the t_next point, then I will need to use the data
                //from ya to integrate.
                if(t_next < temp_t)
                {
                    temp_t = ta;
                    memcpy(temp_y, ya, sizeof(ya));
                }
                
                int status2 = gsl_odeiv2_driver_apply(driver, &temp_t, t_next, temp_y);
                if(status2 != GSL_SUCCESS)
                {
                    fprintf(stderr, "Error with return: %d\n", status2);
                    break;
                }
                
                if(temp_y[1]*yb[1] > 0)
                {
                    tb = temp_t;
                    memcpy(yb, temp_y, sizeof(temp_y));
                }else{
                    ta = temp_t;
                    memcpy(ya, temp_y, sizeof(temp_y));
                }
            }
            
            printf("%.5e %.5e\n", temp_y[0], temp_y[2]);
            
        }
    }
    
    gsl_odeiv2_driver_free(driver);
    return 0;
}