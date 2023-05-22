/*-----------------------------------------------------
Cable equation with the ten Tusscher model 2006
Author: Guilherme Couto
FISIOCOMP - UFJF
------------------------------------------------------*/

#include "./include/parameters.h"
#include "./include/functions.h"

/*-----------------------------------------------------
Model parameters
-----------------------------------------------------*/
double sigma_long = 1.334;          // Conductivity -> mS/cm
double sigma_trans = 0.176;         // Conductivity -> mS/cm
double sigma = 1.171;               // Diffusion coefficient -> cm²/s


/*-----------------------------------------------------
Simulation parameters
-----------------------------------------------------*/
int L = 2.0;            // Length of each side -> cm
double dx = 0.02;       // Spatial step -> cm
double sim_time = 500.0;       // Simulation time -> ms


/*-----------------------------------------------------
Stimulation parameters
-----------------------------------------------------*/
double stim_strength = -38;          // Stimulation strength -> uA/cm^2 

double t_s1_begin = 0.0;            // Stimulation start time -> ms
double stim_duration = 2.0;         // Stimulation duration -> ms
double s1_x_limit = 0.04;            // Stimulation x limit -> cm


/*----------------------------------------
Main function
-----------------------------------------*/
int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        printf("Usage: %s <num_threads> <delta_t (ms)>\n", argv[0]);
        exit(1);
    }

    int num_threads = atoi(argv[1]);
    double dt = atof(argv[2]);

    if (num_threads <= 0)
    {
        fprintf(stderr, "Number of threads must greater than 0\n");
        exit(1);
    }

    // Number of steps
    int N = (int)(L / dx);  // Number of spatial steps (square tissue)
    int M = (int)(T / dt);  // Number of time steps

    // Variables
    double *V, *X_r1, *X_r2, *X_s, *m, *h, *j, *d, *f, *f2, *fCass, *s, *r, *Ca_i, *Ca_SR, *Ca_SS, *R_prime, *Na_i, *K_i;

    // Allocate memory for variables
    V = (double *)malloc(N * sizeof(double));
    X_r1 = (double *)malloc(N * sizeof(double ));
    X_r2 = (double *)malloc(N * sizeof(double ));
    X_s = (double *)malloc(N * sizeof(double ));
    m = (double *)malloc(N * sizeof(double ));
    h = (double *)malloc(N * sizeof(double ));
    j = (double *)malloc(N * sizeof(double ));
    d = (double *)malloc(N * sizeof(double ));
    f = (double *)malloc(N * sizeof(double ));
    f2 = (double *)malloc(N * sizeof(double ));
    fCass = (double *)malloc(N * sizeof(double ));
    s = (double *)malloc(N * sizeof(double ));
    r = (double *)malloc(N * sizeof(double ));
    Ca_i = (double *)malloc(N * sizeof(double ));
    Ca_SR = (double *)malloc(N * sizeof(double ));
    Ca_SS = (double *)malloc(N * sizeof(double ));
    R_prime = (double *)malloc(N * sizeof(double ));
    Na_i = (double *)malloc(N * sizeof(double ));
    K_i = (double *)malloc(N * sizeof(double ));

    double dR_prime_dt, dCa_SR_dt, dCa_SS_dt, dCa_i_dt, dNa_i_dt, dK_i_dt, diff_term = 0.0;
    double V_step, X_r1_step, X_r2_step, X_s_step, m_step, h_step, j_step, d_step, f_step, f2_step, fCass_step, s_step, r_step, Ca_i_step, Ca_SR_step, Ca_SS_step, R_prime_step, Na_i_step, K_i_step;

    int step = 0;
    double tstep = 0.0;
    double *time = (double *)malloc(M * sizeof(double));
    int n;  // Time index
    for (n = 0; n < M; n++)
    {
        time[n] = n * dt;
    }

    // Stim variables
    double I_stim = 0.0;
    double I_total = 0.0;
    int x_lim = s1_x_limit / dx;

    // Phi
    // double phi_long = sigma_long * dt / (chi * dx * dx);
    // double phi_trans = sigma_trans * dt / (chi * dx * dx);
    double phi = sigma * dt / (chi * dx * dx);      // Isotropic diffusion

    // Initial conditions
    int i;
    for (i = 0; i < N; i++)
    {
        V[i] = V_init;
        X_r1[i] = X_r1_init;
        X_r2[i] = X_r2_init;
        X_s[i] = X_s_init;
        m[i] = m_init;
        h[i] = h_init;
        j[i] = j_init;
        d[i] = d_init;
        f[i] = f_init;
        f2[i] = f2_init;
        fCass[i] = fCass_init;
        s[i] = s_init;
        r[i] = r_init;
        Ca_i[i] = Ca_i_init;
        Ca_SR[i] = Ca_SR_init;
        Ca_SS[i] = Ca_SS_init;
        R_prime[i] = R_prime_init;
        Na_i[i] = Na_i_init;
        K_i[i] = K_i_init;
    }

    // Prepare files to save data
    system("mkdir -p simulation-files");

    // Open the file to write for complete gif
    char fname_complete[100];
    sprintf(fname_complete, "./simulation-files/tnnp-cable-eq-%.03f.txt", dt);
    FILE *fp_all = NULL;
    fp_all = fopen(fname_complete, "w");

    int save_rate = ceil(M / 150.0);

    // Open the file to write for times
    char fname_times[100];
    sprintf(fname_times, "./simulation-files/times-cable-eq-%.03f.txt", dt);
    FILE *fp_times = NULL;
    fp_times = fopen(fname_times, "w");

    // For velocity
    bool tag = true;
    double velocity = 0.0;
    
    // Timer
    double start, finish, elapsed = 0.0;

    start = omp_get_wtime();

    // Forward Euler
    #pragma omp parallel num_threads(num_threads) default(none) \
    private(i, I_stim, dR_prime_dt, dCa_SR_dt, dCa_SS_dt, dCa_i_dt, dNa_i_dt, dK_i_dt, diff_term, I_total, \
    V_step, X_r1_step, X_r2_step, X_s_step, m_step, h_step, j_step, d_step, f_step, f2_step, fCass_step, s_step, r_step, Ca_i_step, Ca_SR_step, Ca_SS_step, R_prime_step, Na_i_step, K_i_step) \
    shared(N, M, V, X_r1, X_r2, X_s, m, h, j, d, f, f2, fCass, s, r, Ca_i, Ca_SR, Ca_SS, R_prime, Na_i, K_i, \
    dt, L, stim_strength, s1_x_limit, t_s1_begin, stim_duration, x_lim, time, phi, sim_time, tstep, step, fp_all, velocity, save_rate, fp_times, tag)
    {
        while (step < M)
        {
            // Get time step
            tstep = time[step];

            #pragma omp for
            for (i = 1; i < N - 1; i++)
            {
                // Stimulus 1
                if (tstep >= t_s1_begin && tstep <= t_s1_begin + stim_duration && i <= x_lim)
                {
                    I_stim = stim_strength;
                }
                else 
                {
                    I_stim = 0.0;
                }

                // Get values at current time step and space
                V_step = V[i];
                X_r1_step = X_r1[i];
                X_r2_step = X_r2[i];
                X_s_step = X_s[i];
                m_step = m[i];
                h_step = h[i];
                j_step = j[i];
                d_step = d[i];
                f_step = f[i];
                f2_step = f2[i];
                fCass_step = fCass[i];
                s_step = s[i];
                r_step = r[i];
                Ca_i_step = Ca_i[i];
                Ca_SR_step = Ca_SR[i];
                Ca_SS_step = Ca_SS[i];
                R_prime_step = R_prime[i];
                Na_i_step = Na_i[i];
                K_i_step = K_i[i];

                // Update total current
                I_total = Itotal(I_stim, V_step, m_step, h_step, j_step, Na_i_step, K_i_step, r_step, s_step, X_r1_step, X_r2_step, X_s_step, d_step, f_step, f2_step, fCass_step, Ca_SS_step, Ca_i_step);

                // Update voltage
                V[i] = V_step + ((-I_total) * dt) + (phi * (V[i - 1] - 2.0 * V[i] + V[i + 1]));

                // Update concentrations
                dR_prime_dt = dRprimedt(Ca_SS_step, R_prime_step);
                dCa_i_dt = dCaidt(Ca_i_step, Ca_SR_step, Ca_SS_step, V_step, Na_i_step); 
                dCa_SR_dt = dCaSRdt(Ca_SR_step, Ca_i_step, Ca_SS_step, R_prime_step);
                dCa_SS_dt = dCaSSdt(Ca_SS_step, V_step, d_step, f_step, f2_step, fCass_step, Ca_SR_step, R_prime_step, Ca_i_step);
                dNa_i_dt = dNaidt(V_step, m_step, h_step, j_step, Na_i_step, Ca_i_step);
                dK_i_dt = dKidt(I_stim, V_step, K_i_step, r_step, s_step, X_r1_step, X_r2_step, X_s_step, Na_i_step);

                R_prime[i] = R_prime_step + dR_prime_dt * dt;
                Ca_SR[i] = Ca_SR_step + dCa_SR_dt * dt;
                Ca_SS[i] = Ca_SS_step + dCa_SS_dt * dt;
                Ca_i[i] = Ca_i_step + dCa_i_dt * dt;
                Na_i[i] = Na_i_step + dNa_i_dt * dt;
                K_i[i] = K_i_step + dK_i_dt * dt;

                // Update variables - Rush Larsen
                X_r1[i] = updateXr1(X_r1_step, V_step, dt);
                X_r2[i] = updateXr2(X_r2_step, V_step, dt);
                X_s[i] = updateXs(X_s_step, V_step, dt);
                r[i] = updater(r_step, V_step, dt);
                s[i] = updates(s_step, V_step, dt);
                m[i] = updatem(m_step, V_step, dt);
                h[i] = updateh(h_step, V_step, dt);
                j[i] = updatej(j_step, V_step, dt);
                d[i] = updated(d_step, V_step, dt);
                f[i] = updatef(f_step, V_step, dt);
                f2[i] = updatef2(f2_step, V_step, dt);
                fCass[i] = updatefCass(fCass_step, V_step, dt);
            }

            // Boundary conditions
            #pragma omp for nowait
            for (i = 0; i < N; i++)
            {
                V[0] = V[1];
                V[N-1] = V[N-2];
            }

            // Save data to file
            #pragma omp master
            {
                // Write to file
                if (step % save_rate == 0)
                {
                    for (int i = 0; i < N; i++)
                    {
                        fprintf(fp_all, "%lf\n", V[i]);
                    }
                    fprintf(fp_times, "%lf\n", time[step]);
                }

                // Check S1 velocity
                if (V[N-1] > 20.0 && tag)
                {
                    velocity = ((10*(L - s1_x_limit)) / (time[step]));
                    printf("S1 velocity: %lf\n", velocity);
                    tag = false;
                }
            } 
            
            // Update step
            #pragma omp master
            {
                step++;
            }
            #pragma omp barrier
        }
    }

    // Check time
    finish = omp_get_wtime();
    elapsed = finish - start;

    printf("\nElapsed time = %e seconds\n", elapsed);

    // Close files
    fclose(fp_all);
    fclose(fp_times);
    
    // Free alocated memory
    free(time);
    free(V);
    free(X_r1);
    free(X_r2);
    free(X_s);
    free(m);
    free(h);
    free(j);
    free(d);
    free(f);
    free(f2);
    free(fCass);
    free(s);
    free(r);
    free(Ca_i);
    free(Ca_SR);
    free(Ca_SS);
    free(R_prime);
    free(Na_i);
    free(K_i);

    return 0;
}
