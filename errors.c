#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

const int N = 200;

void print_errors(FILE *fp, double sum, double sum_abs, double norm_01)
{
    double rmse = sqrt(sum / (N*N));    // Root Mean Square Error
    double ed = sqrt(sum);    // Euclidean distance

    fprintf(fp, "Euclidean distance = %lf\n", ed);
    fprintf(fp, "RMSE = %lf\n", rmse);
    fprintf(fp, "Mean absolute error = %lf\n", sum_abs / (N*N));
    fprintf(fp, "ED / Norm_ref = %lf\n", (ed / norm_01));
    fprintf(fp, "Relative %% = %lf\n", (ed / norm_01) * 100);
    fprintf(fp, "\n\n");
}

int main()
{
    FILE *fp;
    fp = fopen("tnnp-errors.txt", "w");
    FILE *fp_1;
    char *ptr;


    // Compare EXP with dt = 0.02 and EXP with dt_edo = dt_pde = 0.02
    fp_1 = fopen("last-1.0-0.020-0.020-final.txt", "r");
    FILE *fp_adi;
    fp_adi = fopen("last-1.0-0.020-0.020-final.txt", "r");
    fprintf(fp, "Error (EXP with dt = 0.02 and EXP with dt_edo = dt_pde = 0.02):\n");
    double sum = 0;
    double sum_norm_01 = 0;
    double sum_abs = 0;
    for(int i = 0; i < N*N; i ++)
    {
        char line_1[20];
        char line_adi[20];
        fgets(line_1, 20, fp_1);
        fgets(line_adi, 20, fp_adi);
        double num_1 = strtod(line_1, &ptr);
        double num_exp = strtod(line_adi, &ptr);

        double p = pow(num_exp - num_1, 2); 
        sum += p;
        sum_abs += fabs(num_exp - num_1);
        sum_norm_01 += pow(num_1, 2);
    }
    double norm_01 = sqrt(sum_norm_01);    // Normalization
    print_errors(fp, sum, sum_abs, norm_01);
    fclose(fp_adi);
    fclose(fp_1);


    // Compare EXP with dt = 0.02 and EXP dt = 0.04
    fp_1 = fopen("last-1.0-0.020-0.020-final.txt", "r");
    FILE *fp_exp_04;
    fp_exp_04 = fopen("last-1.0-0.040-0.040-final.txt", "r");
    fprintf(fp, "Error (EXP with dt = 0.02 and EXP dt = 0.04):\n");
    sum = 0;
    sum_abs = 0;
    for(int i = 0; i < N*N; i ++)
    {
        char line_1[20];
        char line_exp_04[20];
        fgets(line_1, 20, fp_1);
        fgets(line_exp_04, 20, fp_exp_04);
        double num_1 = strtod(line_1, &ptr);
        double num_exp_04 = strtod(line_exp_04, &ptr);

        double p = pow(num_exp_04 - num_1, 2); 
        sum += p;
        sum_abs += fabs(num_exp_04 - num_1);
    }
    print_errors(fp, sum, sum_abs, norm_01);
    fclose(fp_exp_04);
    fclose(fp_1);

    
    // Compare EXP with dt = 0.02 and EXP dt = 0.06
    fp_1 = fopen("last-1.0-0.020-0.020-final.txt", "r");
    FILE *fp_exp_05;
    fp_exp_05 = fopen("last-1.0-0.060-0.060-final.txt", "r");
    fprintf(fp, "Error (EXP with dt = 0.02 and EXP dt = 0.06):\n");
    sum = 0;
    sum_abs = 0;
    for(int i = 0; i < N*N; i ++)
    {
        char line_1[20];
        char line_exp_05[20];
        fgets(line_1, 20, fp_1);
        fgets(line_exp_05, 20, fp_exp_05);
        double num_1 = strtod(line_1, &ptr);
        double num_exp_05 = strtod(line_exp_05, &ptr);

        double p = pow(num_exp_05 - num_1, 2); 
        sum += p;
        sum_abs += fabs(num_exp_05 - num_1);
    }
    print_errors(fp, sum, sum_abs, norm_01);
    fclose(fp_exp_05);
    fclose(fp_1);


    // Compare EXP with dt = 0.02 and EXP dt = 0.08
    fp_1 = fopen("last-1.0-0.020-0.020-final.txt", "r");
    fp_exp_05 = fopen("last-1.0-0.080-0.080-final.txt", "r");
    fprintf(fp, "Error (EXP with dt = 0.02 and EXP dt = 0.08):\n");
    sum = 0;
    sum_abs = 0;
    for(int i = 0; i < N*N; i ++)
    {
        char line_1[20];
        char line_exp_05[20];
        fgets(line_1, 20, fp_1);
        fgets(line_exp_05, 20, fp_exp_05);
        double num_1 = strtod(line_1, &ptr);
        double num_exp_05 = strtod(line_exp_05, &ptr);

        double p = pow(num_exp_05 - num_1, 2); 
        sum += p;
        sum_abs += fabs(num_exp_05 - num_1);
    }
    print_errors(fp, sum, sum_abs, norm_01);
    fclose(fp_exp_05);
    fclose(fp_1);

    // ----------------------------------------------------------------------------------------------------------
    fprintf(fp, "\n\n============================================================================\n\n");
    // ----------------------------------------------------------------------------------------------------------

    // Compare ADI with dt_edo = dt_pde = 0.02 and ADI with dt_edo = dt_pde = 0.02
    fp_1 = fopen("last-2.0-0.020-0.020-final.txt", "r");
    fp_adi = fopen("last-2.0-0.020-0.020-final.txt", "r");
    fprintf(fp, "Error (ADI with dt_edo = dt_pde = 0.02 and ADI with dt_edo = dt_pde = 0.02):\n");
    sum = 0;
    double sum_norm_02 = 0;
    sum_abs = 0;
    for(int i = 0; i < N*N; i ++)
    {
        char line_1[20];
        char line_adi[20];
        fgets(line_1, 20, fp_1);
        fgets(line_adi, 20, fp_adi);
        double num_1 = strtod(line_1, &ptr);
        double num_exp = strtod(line_adi, &ptr);

        double p = pow(num_exp - num_1, 2); 
        sum += p;
        sum_abs += fabs(num_exp - num_1);
        sum_norm_02 += pow(num_1, 2);
    }
    double norm_02 = sqrt(sum_norm_02);    // Normalization
    print_errors(fp, sum, sum_abs, norm_02);
    fclose(fp_adi);
    fclose(fp_1);


    // Compare ADI with dt_edo = dt_pde = 0.02 and ADI dt_edo = dt_pde = 0.04
    fp_1 = fopen("last-2.0-0.020-0.020-final.txt", "r");
    fp_adi = fopen("last-2.0-0.040-0.040-final.txt", "r");
    fprintf(fp, "Error (ADI with dt_edo = dt_pde = 0.02 and ADI dt_edo = dt_pde = 0.04):\n");
    sum = 0;
    sum_abs = 0;
    for(int i = 0; i < N*N; i ++)
    {
        char line_1[20];
        char line_adi_03[20];
        fgets(line_1, 20, fp_1);
        fgets(line_adi_03, 20, fp_adi);
        double num_1 = strtod(line_1, &ptr);
        double num_adi_03 = strtod(line_adi_03, &ptr);

        double p = pow(num_adi_03 - num_1, 2); 
        sum += p;
        sum_abs += fabs(num_adi_03 - num_1);
    }
    print_errors(fp, sum, sum_abs, norm_02);
    fclose(fp_adi);
    fclose(fp_1);


    // Compare ADI with dt_edo = dt_pde = 0.02 and ADI dt_edo = dt_pde = 0.06
    fp_1 = fopen("last-2.0-0.020-0.020-final.txt", "r");
    fp_adi = fopen("last-2.0-0.060-0.060-final.txt", "r");
    fprintf(fp, "Error (ADI with dt_edo = dt_pde = 0.02 and ADI dt_edo = dt_pde = 0.06):\n");
    sum = 0;
    sum_abs = 0;
    for(int i = 0; i < N*N; i ++)
    {
        char line_1[20];
        char line_adi_03[20];
        fgets(line_1, 20, fp_1);
        fgets(line_adi_03, 20, fp_adi);
        double num_1 = strtod(line_1, &ptr);
        double num_adi_03 = strtod(line_adi_03, &ptr);

        double p = pow(num_adi_03 - num_1, 2); 
        sum += p;
        sum_abs += fabs(num_adi_03 - num_1);
    }
    print_errors(fp, sum, sum_abs, norm_02);
    fclose(fp_adi);
    fclose(fp_1);

    // Compare ADI with dt_edo = dt_pde = 0.02 and ADI dt_edo = dt_pde = 0.08
    fp_1 = fopen("last-2.0-0.020-0.020-final.txt", "r");
    fp_adi = fopen("last-2.0-0.080-0.080-final.txt", "r");
    fprintf(fp, "Error (ADI with dt_edo = dt_pde = 0.02 and ADI dt_edo = dt_pde = 0.08):\n");
    sum = 0;
    sum_abs = 0;
    for(int i = 0; i < N*N; i ++)
    {
        char line_1[20];
        char line_adi_03[20];
        fgets(line_1, 20, fp_1);
        fgets(line_adi_03, 20, fp_adi);
        double num_1 = strtod(line_1, &ptr);
        double num_adi_03 = strtod(line_adi_03, &ptr);

        double p = pow(num_adi_03 - num_1, 2); 
        sum += p;
        sum_abs += fabs(num_adi_03 - num_1);
    }
    print_errors(fp, sum, sum_abs, norm_02);
    fclose(fp_adi);
    fclose(fp_1);


    fclose(fp);
    return 0;
}