#ifndef AUXFUNCS_H
#define AUXFUNCS_H

#include "parameters.h"

/*----------------------------------------
Auxiliary functions
-----------------------------------------*/
// Initialize variables
void initialize_variables(int N, double **V, double **X_r1, double **X_r2, double **X_s, double **m, double **h, double **j, double **d, double **f, double **f2, double **fCass, double **s, double **r, double **Ca_i, double **Ca_SR, double **Ca_SS, double **R_prime, double **Na_i, double **K_i)
{
    for (int i = 0; i < N; i++)
    {
        for (int k = 0; k < N; k++)
        {
            V[i][k] = V_init;
            X_r1[i][k] = X_r1_init;
            X_r2[i][k] = X_r2_init;
            X_s[i][k] = X_s_init;
            m[i][k] = m_init;
            h[i][k] = h_init;
            j[i][k] = j_init;
            d[i][k] = d_init;
            f[i][k] = f_init;
            f2[i][k] = f2_init;
            fCass[i][k] = fCass_init;
            s[i][k] = s_init;
            r[i][k] = r_init;
            Ca_i[i][k] = Ca_i_init;
            Ca_SR[i][k] = Ca_SR_init;
            Ca_SS[i][k] = Ca_SS_init;
            R_prime[i][k] = R_prime_init;
            Na_i[i][k] = Na_i_init;
            K_i[i][k] = K_i_init;
        }
    }
}

// Adapted for 2nd order approximation
void thomas_algorithm_2nd(double *d, double *solution, unsigned long N, double phi, double *c_, double *d_)
{   
    // Coefficients
    double a = -phi;    // subdiagonal
    double b = 1 + 2 * phi; // diagonal (1st and last row)
    double c = - 2 * phi;    // superdiagonal
    
    // 1st: update auxiliary arrays
    c_[0] = c / b;
    d_[0] = d[0] / b;

    c = -phi;
    
    for (int i = 1; i <= N - 2; i++)
    {
        c_[i] = c / (b - a * c_[i - 1]);
        d_[i] = (d[i] - a * d_[i - 1]) / (b - a * c_[i - 1]);
    }
    
    a = - 2 * phi;
    d_[N - 1] = (d[N - 1] - a * d_[N - 2]) / (b - a * c_[N - 2]);

    a = -phi;

    // 2nd: update solution
    solution[N - 1] = d_[N - 1];
    
    for (int i = N - 2; i >= 0; i--)
    {
        solution[i] = d_[i] - c_[i] * solution[i + 1];
    }
}

// Adapted for 2nd order approximation
double diffusion_i_2nd(int i, int j, int N, double **v)
{
    double result = 0.0;
    if (i == 0)
    {
        result = - 2.0*v[i][j] + 2.0*v[i + 1][j]; 
    }
    else if (i == N - 1)
    {
        result = 2.0*v[i - 1][j] - 2.0*v[i][j]; 
    }
    else
    {
        result = v[i - 1][j] - 2.0*v[i][j] + v[i + 1][j];
    }

    return result;
}

// Adapted for 2nd order approximation
double diffusion_j_2nd(int i, int j, int N, double **v)
{
    double result = 0.0;
    if (j == 0)
    {
        result = - 2.0*v[i][j] + 2.0*v[i][j + 1]; 
    }
    else if (j == N - 1)
    {
        result = 2.0*v[i][j - 1] - 2.0*v[i][j]; 
    }
    else
    {
        result = v[i][j - 1] - 2.0*v[i][j] + v[i][j + 1];
    }

    return result;
}

// Normalize voltage
void normalize_V_2d(int N, double **V)
{
    double V_max = V[0][0];
    double V_min = V[0][0];

    for (int i = 0; i < N; i++)
    {
        for (int k = 0; k < N; k++)
        {
            if (V[i][k] > V_max)
            {
                V_max = V[i][k];
            }
            if (V[i][k] < V_min)
            {
                V_min = V[i][k];
            }
        }
    }

    double V_range = V_max - V_min;

    for (int i = 0; i < N; i++)
    {
        for (int k = 0; k < N; k++)
        {
            V[i][k] = (V[i][k] - V_min) / V_range;
        }
    }
}

void normalize_V_1d(int N, double *V)
{
    double V_max = V[0];
    double V_min = V[0];

    for (int i = 0; i < N; i++)
    {
        if (V[i] > V_max)
        {
            V_max = V[i];
        }
        if (V[i] < V_min)
        {
            V_min = V[i];
        }
    }

    double V_range = V_max - V_min;

    for (int i = 0; i < N; i++)
    {
        V[i] = (V[i] - V_min) / V_range;
    }
}

#endif // AUXFUNCS_H