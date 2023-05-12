/*-----------------------------------------------------
Single cell with the ten Tusscher model 2006
Author: Guilherme Couto
FISIOCOMP - UFJF
------------------------------------------------------*/

#include "./include/parameters.h"
#include "./include/functions.h"

/*-----------------------------------------------------
Simulation parameters
-----------------------------------------------------*/
double sim_time = 600.0;       // Simulation time -> ms


/*-----------------------------------------------------
Stimulation parameters
-----------------------------------------------------*/
double stim_strength = -38;          // Stimulation strength -> uA/cm^2 

double t_s1_begin = 50.0;            // Stimulation start time -> ms
double stim_duration = 1.0;         // Stimulation duration -> ms


/*----------------------------------------
Main function
-----------------------------------------*/
int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        printf("Usage: %s <delta_t (ms)>\n", argv[0]);
        exit(1);
    }

    double dt = atof(argv[1]);

    if (dt <= 0)
    {
        fprintf(stderr, "Delta_t must greater than 0\n");
        exit(1);
    }

    // Number of steps
    int M = (int)(sim_time / dt);  // Number of time steps

    // Variables
    double *V_array = (double *)malloc(M * sizeof(double));
    double V, X_r1, X_r2, X_s, m, h, j, d, f, f2, fCass, s, r, Ca_i, Ca_SR, Ca_SS, R_prime, Na_i, K_i;

    double dR_prime_dt, dCa_SR_dt, dCa_SS_dt, dCa_i_dt, dNa_i_dt, dK_i_dt;
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

    // Initial conditions
    V = V_init;
    X_r1 = X_r1_init;
    X_r2 = X_r2_init;
    X_s = X_s_init;
    m = m_init;
    h = h_init;
    j = j_init;
    d = d_init;
    f = f_init;
    f2 = f2_init;
    fCass = fCass_init;
    s = s_init;
    r = r_init;
    Ca_i = Ca_i_init;
    Ca_SR = Ca_SR_init;
    Ca_SS = Ca_SS_init;
    R_prime = R_prime_init;
    Na_i = Na_i_init;
    K_i = K_i_init;

    // Prepare files to save data
    // Convert dt to string
    char s_dt[10];
    sprintf(s_dt, "%.03f", dt);

    system("mkdir -p simulation-files");

    // Open the file to write for complete gif
    char fname_complete[100] = "./simulation-files/tnnp-";
    strcat(fname_complete, "cell-norm");
    strcat(fname_complete, "-");
    strcat(fname_complete, s_dt);
    strcat(fname_complete, ".txt");
    FILE *fp_all = NULL;
    fp_all = fopen(fname_complete, "w");
    
    // Open the file to write for times
    char fname_times[100] = "./simulation-files/sim-times-";
    strcat(fname_times, "cell-norm");
    strcat(fname_times, "-");
    strcat(fname_times, s_dt);
    strcat(fname_times, ".txt");
    FILE *fp_times = NULL;
    fp_times = fopen(fname_times, "w");
    
    // Timer
    double start, finish, elapsed = 0.0;

    start = omp_get_wtime();

    // Forward Euler
    while (step < M)
    {
        // Get time step
        tstep = time[step];

        // Stimulus 1
        if (tstep >= t_s1_begin && tstep <= t_s1_begin + stim_duration)
        {
            I_stim = stim_strength;
        }
        else 
        {
            I_stim = 0.0;
        }

        // Get values at current time step and space
        V_step = V;
        X_r1_step = X_r1;
        X_r2_step = X_r2;
        X_s_step = X_s;
        m_step = m;
        h_step = h;
        j_step = j;
        d_step = d;
        f_step = f;
        f2_step = f2;
        fCass_step = fCass;
        s_step = s;
        r_step = r;
        Ca_i_step = Ca_i;
        Ca_SR_step = Ca_SR;
        Ca_SS_step = Ca_SS;
        R_prime_step = R_prime;
        Na_i_step = Na_i;
        K_i_step = K_i;

        // Update total current
        I_total = Itotal(I_stim, V_step, m_step, h_step, j_step, Na_i_step, K_i_step, r_step, s_step, X_r1_step, X_r2_step, X_s_step, d_step, f_step, f2_step, fCass_step, Ca_SS_step, Ca_i_step);

        // Update voltage
        V = V_step + (-I_total) * dt;
        V_array[step] = V;

        // Update concentrations
        dR_prime_dt = dRprimedt(Ca_SS_step, R_prime_step);
        dCa_i_dt = dCaidt(Ca_i_step, Ca_SR_step, Ca_SS_step, V_step, Na_i_step); 
        dCa_SR_dt = dCaSRdt(Ca_SR_step, Ca_i_step, Ca_SS_step, R_prime_step);
        dCa_SS_dt = dCaSSdt(Ca_SS_step, V_step, d_step, f_step, f2_step, fCass_step, Ca_SR_step, R_prime_step, Ca_i_step);
        dNa_i_dt = dNaidt(V_step, m_step, h_step, j_step, Na_i_step, Ca_i_step);
        dK_i_dt = dKidt(I_stim, V_step, K_i_step, r_step, s_step, X_r1_step, X_r2_step, X_s_step, Na_i_step);

        R_prime = R_prime_step + dR_prime_dt * dt;
        Ca_SR = Ca_SR_step + dCa_SR_dt * dt;
        Ca_SS = Ca_SS_step + dCa_SS_dt * dt;
        Ca_i = Ca_i_step + dCa_i_dt * dt;
        Na_i = Na_i_step + dNa_i_dt * dt;
        K_i = K_i_step + dK_i_dt * dt;

        // Update variables - Rush Larsen
        X_r1 = updateXr1(X_r1_step, V_step, dt);
        X_r2 = updateXr2(X_r2_step, V_step, dt);
        X_s = updateXs(X_s_step, V_step, dt);
        r = updater(r_step, V_step, dt);
        s = updates(s_step, V_step, dt);
        m = updatem(m_step, V_step, dt);
        h = updateh(h_step, V_step, dt);
        j = updatej(j_step, V_step, dt);
        d = updated(d_step, V_step, dt);
        f = updatef(f_step, V_step, dt);
        f2 = updatef2(f2_step, V_step, dt);
        fCass = updatefCass(fCass_step, V_step, dt);

        // Save data to file
        // Write to file
        fprintf(fp_times, "%lf\n", time[step]);
        
        // Update step
        step++;
    }

    // Save data to file
    // Write to file
    normalize_V_1d(M, V_array);
    for (int i = 0; i < M; i++)
    {
        fprintf(fp_all, "%lf\n", V_array[i]);
    }

    // Check time
    finish = omp_get_wtime();
    elapsed = finish - start;

    printf("\nElapsed time = %e seconds\n", elapsed);

    // Close files
    fclose(fp_all);
    fclose(fp_times);
    
    // Free alocated memory
    free(V_array);
    free(time);

    return 0;
}