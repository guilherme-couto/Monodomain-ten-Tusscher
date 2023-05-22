/*-----------------------------------------------------
APD90 x BCL with the ten Tusscher model 2006
Author: Guilherme Couto
FISIOCOMP - UFJF
------------------------------------------------------*/

#include "./include/parameters.h"
#include "./include/functions.h"

/*----------------------------------------
Main function
-----------------------------------------*/
int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        printf("Usage: %s <BCL (ms)>\n", argv[0]);
        exit(1);
    }

    double dt = 0.05;                           // Same dt used with the Genetic Algorithm
    double BCL = atof(argv[1]);                 // Basic Cycle Length -> ms
    int num_stim = 20 + 1;                      // Number of stimulations
    double sim_time = BCL * num_stim + 500.0;   // Simulation time -> ms (500 ms for guarantee)

    if (BCL <= 0)
    {
        fprintf(stderr, "BCL must greater than 0\n");
        exit(1);
    }

    // Stimulation parameters
    double stim_strength = -38;          // Stimulation strength -> uA/cm^2 
    double t_stim_begin = 50.0;          // Stimulation start time -> ms
    double stim_duration = 1.0;          // Stimulation duration -> ms
    int stim_count = 0;                  // Stimulation counter

    // Number of steps
    int M = (int)(sim_time / dt);  // Number of time steps

    // Variables
    double V, X_r1, X_r2, X_s, m, h, j, d, f, f2, fCass, s, r, Ca_i, Ca_SR, Ca_SS, R_prime, Na_i, K_i;
    double V_max = -100.0, V_rest;

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
    system("mkdir -p simulation-files");

    // Open the file to write for complete gif
    char fname_complete[100];
    sprintf(fname_complete, "./simulation-files/tnnp-%s-%.02f-%.01f.txt", "pacing", dt, BCL);
    FILE *fp_all = NULL;
    fp_all = fopen(fname_complete, "w");
    
    // Open the file to write for times
    char fname_times[100];
    sprintf(fname_times, "./simulation-files/times-%s-%.02f-%.01f.txt", "pacing", dt, BCL);
    FILE *fp_times = NULL;
    fp_times = fopen(fname_times, "w");

    // APD controlers
    double APD90 = 0.0;
    bool v_max_found = false;
    bool compute_apd = false;
    
    // Timer
    double start, finish, elapsed = 0.0;

    start = omp_get_wtime();

    // Forward Euler
    while (step < M)
    {
        // Get time step
        tstep = time[step];

        // Stimulus
        if (tstep >= t_stim_begin && tstep <= t_stim_begin + stim_duration)
        {
            I_stim = stim_strength;

            // Update stimulus counter
            if (tstep == t_stim_begin)
            {
                V_rest = V;
                stim_count++;
            }
            // Update new stimulus time
            else if (tstep == t_stim_begin + stim_duration && stim_count < num_stim)
            {
                t_stim_begin += BCL;
            }
        }
        else 
        {
            I_stim = 0.0;
        }

        // For the last stimulus, calculate the amplitude of the action potential for APD calculation
        if (stim_count == num_stim && !v_max_found)
        {
            if (V > V_max)
            {
                V_max = V;
            }
            else
            {
                v_max_found = true;
                compute_apd = true;
            }
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
        fprintf(fp_all, "%lf\n", V);
        fprintf(fp_times, "%lf\n", time[step]);

        // Calculate APD90
        if (stim_count == num_stim && tstep > t_stim_begin + stim_duration && V < V_rest + (V_max - V_rest) * 0.1 && compute_apd)
        {
            APD90 = tstep - (t_stim_begin + stim_duration);

            // Print APD90
            printf("APD90 = %.03f\n", APD90);

            // Open file to write APD90
            FILE *fp_apd = NULL;
            fp_apd = fopen("./simulation-files/apd90.txt", "a");

            // Write APD90 to file
            fprintf(fp_apd, "%.01f | ", BCL);
            fprintf(fp_apd, "%.03lf\n", APD90);

            // APD90 calculated
            compute_apd = false;

            // Close file
            fclose(fp_apd);
        }
        
        // Update step
        step++;
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

    return 0;
}