/*-----------------------------------------------------
Monodomain with the ten Tusscher model 2006
Author: Guilherme Couto
FISIOCOMP - UFJF
------------------------------------------------------*/

/*------------------------------------------------------
Electrophysiological parameters for ten Tusscher model 2006 (https://journals.physiology.org/doi/full/10.1152/ajpheart.00109.2006)
from https://tbb.bio.uu.nl/khwjtuss/SourceCodes/HVM2/Source/Main.cc - ten Tusscher code
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3263775/ - Benchmark
and https://github.com/rsachetto/MonoAlg3D_C/blob/master/src/models_library/ten_tusscher/ten_tusscher_2006_RS_CPU.c - Sachetto MonoAlg3D
--------------------------------------------------------*/

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
double dy = 0.02;       // Spatial step -> cm
double sim_time = 400.0;       // Simulation time -> ms


/*-----------------------------------------------------
Stimulation parameters
-----------------------------------------------------*/
double stim_strength = -38;          // Stimulation strength -> uA/cm^2 

double t_s1_begin = 0.0;            // Stimulation start time -> ms
double stim_duration = 2.0;         // Stimulation duration -> ms
double s1_x_limit = 0.04;            // Stimulation x limit -> cm

double t_s2_begin = 310.0;          // Stimulation start time -> ms
double stim2_duration = 2.0;        // Stimulation duration -> ms
double s2_x_max = 1.0;              // Stimulation x max -> cm
double s2_y_max = 1.0;              // Stimulation y limit -> cm
double s2_x_min = 0.0;              // Stimulation x min -> cm
double s2_y_min = 0.0;              // Stimulation y min -> cm


/*----------------------------------------
Main function
-----------------------------------------*/
int main(int argc, char *argv[])
{
    if (argc != 4)
    {
        printf("Usage: %s <num_threads> <delta_t (ms)> <method>\n", argv[0]);
        exit(1);
    }

    int num_threads = atoi(argv[1]);
    double dt = atof(argv[2]);
    char *method = argv[3];

    if (num_threads <= 0)
    {
        fprintf(stderr, "Number of threads must greater than 0\n");
        exit(1);
    }
    if (strcmp(method, "ADI2") != 0 && strcmp(method, "FE") != 0)
    {
        fprintf(stderr, "Method must be ADI2 (second order) or FE\n");
        exit(1);
    }

    // Number of steps
    int N = (int)(L / dx);  // Number of spatial steps (square tissue)
    int M = (int)(sim_time / dt);  // Number of time steps

    // Variables
    double **V, **X_r1, **X_r2, **X_s, **m, **h, **j, **d, **f, **f2, **fCass, **s, **r, **Ca_i, **Ca_SR, **Ca_SS, **R_prime, **Na_i, **K_i;
    double **V_aux, **X_r1_aux, **X_r2_aux, **X_s_aux, **m_aux, **h_aux, **j_aux, **d_aux, **f_aux, **f2_aux, **fCass_aux, **s_aux, **r_aux, **Ca_i_aux, **Ca_SR_aux, **Ca_SS_aux, **R_prime_aux, **Na_i_aux, **K_i_aux;
    
    // Auxiliar arrays
    double **reaction_V, **rightside, **solution, **c_, **d_;

    // Allocate memory for variables
    V = (double **)malloc(N * sizeof(double *));
    X_r1 = (double **)malloc(N * sizeof(double *));
    X_r2 = (double **)malloc(N * sizeof(double *));
    X_s = (double **)malloc(N * sizeof(double *));
    m = (double **)malloc(N * sizeof(double *));
    h = (double **)malloc(N * sizeof(double *));
    j = (double **)malloc(N * sizeof(double *));
    d = (double **)malloc(N * sizeof(double *));
    f = (double **)malloc(N * sizeof(double *));
    f2 = (double **)malloc(N * sizeof(double *));
    fCass = (double **)malloc(N * sizeof(double *));
    s = (double **)malloc(N * sizeof(double *));
    r = (double **)malloc(N * sizeof(double *));
    Ca_i = (double **)malloc(N * sizeof(double *));
    Ca_SR = (double **)malloc(N * sizeof(double *));
    Ca_SS = (double **)malloc(N * sizeof(double *));
    R_prime = (double **)malloc(N * sizeof(double *));
    Na_i = (double **)malloc(N * sizeof(double *));
    K_i = (double **)malloc(N * sizeof(double *));

    V_aux = (double **)malloc(N * sizeof(double *));
    X_r1_aux = (double **)malloc(N * sizeof(double *));
    X_r2_aux = (double **)malloc(N * sizeof(double *));
    X_s_aux = (double **)malloc(N * sizeof(double *));
    m_aux = (double **)malloc(N * sizeof(double *));
    h_aux = (double **)malloc(N * sizeof(double *));
    j_aux = (double **)malloc(N * sizeof(double *));
    d_aux = (double **)malloc(N * sizeof(double *));
    f_aux = (double **)malloc(N * sizeof(double *));
    f2_aux = (double **)malloc(N * sizeof(double *));
    fCass_aux = (double **)malloc(N * sizeof(double *));
    s_aux = (double **)malloc(N * sizeof(double *));
    r_aux = (double **)malloc(N * sizeof(double *));
    Ca_i_aux = (double **)malloc(N * sizeof(double *));
    Ca_SR_aux = (double **)malloc(N * sizeof(double *));
    Ca_SS_aux = (double **)malloc(N * sizeof(double *));
    R_prime_aux = (double **)malloc(N * sizeof(double *));
    Na_i_aux = (double **)malloc(N * sizeof(double *));
    K_i_aux = (double **)malloc(N * sizeof(double *));

    reaction_V = (double **)malloc(N * sizeof(double *));
    rightside = (double **)malloc(N * sizeof(double *));
    solution = (double **)malloc(N * sizeof(double *));
    c_ = (double **)malloc(N * sizeof(double *));
    d_ = (double **)malloc(N * sizeof(double *));

    for (int i = 0; i < N; i++)
    {
        V[i] = (double *)malloc(N * sizeof(double));
        X_r1[i] = (double *)malloc(N * sizeof(double));
        X_r2[i] = (double *)malloc(N * sizeof(double));
        X_s[i] = (double *)malloc(N * sizeof(double));
        m[i] = (double *)malloc(N * sizeof(double));
        h[i] = (double *)malloc(N * sizeof(double));
        j[i] = (double *)malloc(N * sizeof(double));
        d[i] = (double *)malloc(N * sizeof(double));
        f[i] = (double *)malloc(N * sizeof(double));
        f2[i] = (double *)malloc(N * sizeof(double));
        fCass[i] = (double *)malloc(N * sizeof(double));
        s[i] = (double *)malloc(N * sizeof(double));
        r[i] = (double *)malloc(N * sizeof(double));
        Ca_i[i] = (double *)malloc(N * sizeof(double));
        Ca_SR[i] = (double *)malloc(N * sizeof(double));
        Ca_SS[i] = (double *)malloc(N * sizeof(double));
        R_prime[i] = (double *)malloc(N * sizeof(double));
        Na_i[i] = (double *)malloc(N * sizeof(double));
        K_i[i] = (double *)malloc(N * sizeof(double));

        V_aux[i] = (double *)malloc(N * sizeof(double));
        X_r1_aux[i] = (double *)malloc(N * sizeof(double));
        X_r2_aux[i] = (double *)malloc(N * sizeof(double));
        X_s_aux[i] = (double *)malloc(N * sizeof(double));
        m_aux[i] = (double *)malloc(N * sizeof(double));
        h_aux[i] = (double *)malloc(N * sizeof(double));
        j_aux[i] = (double *)malloc(N * sizeof(double));
        d_aux[i] = (double *)malloc(N * sizeof(double));
        f_aux[i] = (double *)malloc(N * sizeof(double));
        f2_aux[i] = (double *)malloc(N * sizeof(double));
        fCass_aux[i] = (double *)malloc(N * sizeof(double));
        s_aux[i] = (double *)malloc(N * sizeof(double));
        r_aux[i] = (double *)malloc(N * sizeof(double));
        Ca_i_aux[i] = (double *)malloc(N * sizeof(double));
        Ca_SR_aux[i] = (double *)malloc(N * sizeof(double));
        Ca_SS_aux[i] = (double *)malloc(N * sizeof(double));
        R_prime_aux[i] = (double *)malloc(N * sizeof(double));
        Na_i_aux[i] = (double *)malloc(N * sizeof(double));
        K_i_aux[i] = (double *)malloc(N * sizeof(double));

        reaction_V[i] = (double *)malloc(N * sizeof(double));
        rightside[i] = (double *)malloc(N * sizeof(double));
        solution[i] = (double *)malloc(N * sizeof(double));
        c_[i] = (double *)malloc(N * sizeof(double));
        d_[i] = (double *)malloc(N * sizeof(double));
    }

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
    int x_max = s2_x_max / dx;
    int x_min = s2_x_min / dx;
    int y_max = N;
    int y_min = N - s2_y_max / dy;

    // Phi
    // double phi_long = sigma_long * dt / (chi * dx * dx);
    // double phi_trans = sigma_trans * dt / (chi * dx * dx);
    double phi = sigma * dt / (chi * dx * dx);      // Isotropic diffusion

    // Initial conditions
    int i, k;   // i for y-axis and k for x-axis
    initialize_variables(N, V, X_r1, X_r2, X_s, m, h, j, d, f, f2, fCass, s, r, Ca_i, Ca_SR, Ca_SS, R_prime, Na_i, K_i);

    // Prepare files to save data
    // Convert dt to string
    char s_dt[10];
    sprintf(s_dt, "%.03f", dt);

    system("mkdir -p simulation-files");

    // Open the file to write for complete gif
    char fname_complete[100] = "./simulation-files/tnnp-";
    strcat(fname_complete, method);
    strcat(fname_complete, "-");
    strcat(fname_complete, s_dt);
    strcat(fname_complete, ".txt");
    FILE *fp_all = NULL;
    fp_all = fopen(fname_complete, "w");
    int save_rate = ceil(M / 100.0);

    // Open the file to write for times
    char fname_times[100] = "./simulation-files/sim-times-";
    strcat(fname_times, method);
    strcat(fname_times, "-");
    strcat(fname_times, s_dt);
    strcat(fname_times, ".txt");
    FILE *fp_times = NULL;
    fp_times = fopen(fname_times, "w");

    // For velocity
    bool tag = true;
    double velocity = 0.0;
    
    // Timer
    double start, finish, elapsed = 0.0;
    double start_ode, finish_ode, elapsed_ode = 0.0;
    double start_pde, finish_pde, elapsed_pde = 0.0;

    start = omp_get_wtime();

    // Forward Euler
    if (strcmp(method, "FE") == 0)
    {
        // Forward Euler
        #pragma omp parallel num_threads(num_threads) default(none) \
        private(i, k, I_stim, dR_prime_dt, dCa_SR_dt, dCa_SS_dt, dCa_i_dt, dNa_i_dt, dK_i_dt, diff_term, I_total, \
        V_step, X_r1_step, X_r2_step, X_s_step, m_step, h_step, j_step, d_step, f_step, f2_step, fCass_step, s_step, r_step, Ca_i_step, Ca_SR_step, Ca_SS_step, R_prime_step, Na_i_step, K_i_step) \
        shared(N, M, V, X_r1, X_r2, X_s, m, h, j, d, f, f2, fCass, s, r, Ca_i, Ca_SR, Ca_SS, R_prime, Na_i, K_i, \
        dt, L, stim_strength, s1_x_limit, t_s1_begin, stim_duration, x_lim, t_s2_begin, stim2_duration, \
        x_max, y_max, x_min, y_min, time,  c_, d_, V_aux, X_r1_aux, X_r2_aux, X_s_aux, m_aux, h_aux, j_aux, d_aux, f_aux, f2_aux, fCass_aux, s_aux, r_aux, Ca_i_aux, Ca_SR_aux, Ca_SS_aux, R_prime_aux, Na_i_aux, K_i_aux, \
        phi, sim_time, tstep, step, reaction_V, rightside, solution, start_ode, finish_ode, elapsed_ode, start_pde, finish_pde, elapsed_pde, \
        fp_all, velocity, save_rate, fp_times, tag)
        {
            while (step < M)
            {   
                // Get time step
                tstep = time[step];

                // Start ode timer
                #pragma omp master
                {
                    start_ode = omp_get_wtime();
                }
                
                // ODEs - Reaction
                #pragma omp for collapse(2)
                for (i = 1; i < N - 1; i++)
                {   
                    for (k = 1; k < N - 1; k++)
                    {   
                        // Stimulus 1
                        if (tstep >= t_s1_begin && tstep <= t_s1_begin + stim_duration && k <= x_lim)
                        {
                            I_stim = stim_strength;
                        }
                        // Stimulus 2
                        else if (tstep >= t_s2_begin && tstep <= t_s2_begin + stim2_duration && k >= x_min && k <= x_max && i >= y_min && i <= y_max)
                        {
                            I_stim = stim_strength;
                        }
                        else 
                        {
                            I_stim = 0.0;
                        }

                        // Get values at current time step and space
                        V_step = V[i][k];
                        X_r1_step = X_r1[i][k];
                        X_r2_step = X_r2[i][k];
                        X_s_step = X_s[i][k];
                        m_step = m[i][k];
                        h_step = h[i][k];
                        j_step = j[i][k];
                        d_step = d[i][k];
                        f_step = f[i][k];
                        f2_step = f2[i][k];
                        fCass_step = fCass[i][k];
                        s_step = s[i][k];
                        r_step = r[i][k];
                        Ca_i_step = Ca_i[i][k];
                        Ca_SR_step = Ca_SR[i][k];
                        Ca_SS_step = Ca_SS[i][k];
                        R_prime_step = R_prime[i][k];
                        Na_i_step = Na_i[i][k];
                        K_i_step = K_i[i][k];

                        // Update total current
                        I_total = Itotal(I_stim, V_step, m_step, h_step, j_step, Na_i_step, K_i_step, r_step, s_step, X_r1_step, X_r2_step, X_s_step, d_step, f_step, f2_step, fCass_step, Ca_SS_step, Ca_i_step);

                        // Update voltage
                        V[i][k] = V_step + (-I_total) * dt;
                        V_aux[i][k] = V[i][k];

                        // Update concentrations
                        dR_prime_dt = dRprimedt(Ca_SS_step, R_prime_step);
                        dCa_i_dt = dCaidt(Ca_i_step, Ca_SR_step, Ca_SS_step, V_step, Na_i_step); 
                        dCa_SR_dt = dCaSRdt(Ca_SR_step, Ca_i_step, Ca_SS_step, R_prime_step);
                        dCa_SS_dt = dCaSSdt(Ca_SS_step, V_step, d_step, f_step, f2_step, fCass_step, Ca_SR_step, R_prime_step, Ca_i_step);
                        dNa_i_dt = dNaidt(V_step, m_step, h_step, j_step, Na_i_step, Ca_i_step);
                        dK_i_dt = dKidt(I_stim, V_step, K_i_step, r_step, s_step, X_r1_step, X_r2_step, X_s_step, Na_i_step);

                        R_prime[i][k] = R_prime_step + dR_prime_dt * dt;
                        Ca_SR[i][k] = Ca_SR_step + dCa_SR_dt * dt;
                        Ca_SS[i][k] = Ca_SS_step + dCa_SS_dt * dt;
                        Ca_i[i][k] = Ca_i_step + dCa_i_dt * dt;
                        Na_i[i][k] = Na_i_step + dNa_i_dt * dt;
                        K_i[i][k] = K_i_step + dK_i_dt * dt;

                        // Update variables - Rush Larsen
                        X_r1[i][k] = updateXr1(X_r1_step, V_step, dt);
                        X_r2[i][k] = updateXr2(X_r2_step, V_step, dt);
                        X_s[i][k] = updateXs(X_s_step, V_step, dt);
                        r[i][k] = updater(r_step, V_step, dt);
                        s[i][k] = updates(s_step, V_step, dt);
                        m[i][k] = updatem(m_step, V_step, dt);
                        h[i][k] = updateh(h_step, V_step, dt);
                        j[i][k] = updatej(j_step, V_step, dt);
                        d[i][k] = updated(d_step, V_step, dt);
                        f[i][k] = updatef(f_step, V_step, dt);
                        f2[i][k] = updatef2(f2_step, V_step, dt);
                        fCass[i][k] = updatefCass(fCass_step, V_step, dt);
                    }
                }

                // Finish ode timer
                #pragma omp master
                {
                    finish_ode = omp_get_wtime();
                    elapsed_ode += finish_ode - start_ode;
                }

                // Boundary conditions
                #pragma omp for
                for (i = 0; i < N; i++)
                {
                    V_aux[i][0] = V_aux[i][1];
                    V_aux[i][N-1] = V_aux[i][N-2];
                    V_aux[0][i] = V_aux[1][i];
                    V_aux[N-1][i] = V_aux[N-2][i];
                }

                // Start pde timer
                #pragma omp master
                {
                    start_pde = omp_get_wtime();
                }

                // Diffusion
                #pragma omp barrier
                #pragma omp for collapse(2)
                for (i = 1; i < N - 1; i++)
                {
                    for (k = 1; k < N - 1; k++)
                    {   
                        V[i][k] = V_aux[i][k] + (phi * ((V_aux[i - 1][k] - 2.0 * V_aux[i][k] + V_aux[i + 1][k])));
                        V[i][k] += (phi * ((V_aux[i][k - 1] - 2.0 * V_aux[i][k] + V_aux[i][k + 1])));
                    }
                }

                // Finish pde timer
                #pragma omp master
                {
                    finish_pde = omp_get_wtime();
                    elapsed_pde += finish_pde - start_pde;
                }

                // Boundary conditions
                #pragma omp for nowait
                for (i = 0; i < N; i++)
                {
                    V[i][0] = V[i][1];
                    V[i][N-1] = V[i][N-2];
                    V[0][i] = V[1][i];
                    V[N-1][i] = V[N-2][i];
                }

                // Save data to file
                #pragma omp master
                {
                    // Write to file
                    if (step % save_rate == 0)
                    {
                        for (int i = 0; i < N; i++)
                        {
                            for (int j = 0; j < N; j++)
                            {
                                fprintf(fp_all, "%lf\n", V[i][j]);
                            }
                        }
                        fprintf(fp_times, "%lf\n", time[step]);
                    }

                    // Check S1 velocity
                    if (V[0][N-1] > 20.0 && tag)
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
    }

    

    // Check time
    finish = omp_get_wtime();
    elapsed = finish - start;

    printf("\nElapsed time = %e seconds\n", elapsed);

    // Comparison file
    FILE *fp_comp = NULL;
    fp_comp = fopen("comparison.txt", "a");
    fprintf(fp_comp, "%s  \t|\t%d threads\t|\t%.3f ms\t|\t%lf m/s\t|\t%e seconds\n", method, num_threads, dt, velocity, elapsed);

    // ODE/PDE times file
    FILE *fp_all_times = NULL;
    fp_all_times = fopen("times.txt", "a");
    fprintf(fp_all_times, "%s  \t|\t%d threads\t|\t%.3f ms\t|\t%e seconds\t|\t%e seconds\t|\t%e seconds\n", method, num_threads, dt, elapsed_ode, elapsed_pde, elapsed);

    // Close files
    fclose(fp_all);
    fclose(fp_times);
    fclose(fp_comp);
    fclose(fp_all_times);

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

    free(V_aux);
    free(X_r1_aux);
    free(X_r2_aux);
    free(X_s_aux);
    free(m_aux);
    free(h_aux);
    free(j_aux);
    free(d_aux);
    free(f_aux);
    free(f2_aux);
    free(fCass_aux);
    free(s_aux);
    free(r_aux);
    free(Ca_i_aux);
    free(Ca_SR_aux);
    free(Ca_SS_aux);
    free(R_prime_aux);
    free(Na_i_aux);
    free(K_i_aux);
    
    free(reaction_V);
    free(rightside);
    free(solution);
    free(c_);
    free(d_);

    return 0;
}
