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


    // Compare EXP with dt = 0.02 and ADI with dt = 0.02
    fp_1 = fopen("last-exp-0.02.txt", "r");
    FILE *fp_adi;
    fp_adi = fopen("last-adi-0.02.txt", "r");
    fprintf(fp, "Error (EXP with dt = 0.02 and ADI with dt = 0.02):\n");
    double sum = 0;
    double sum_norm_01 = 0;
    double sum_abs = 0;
    for(int i = 0; i < N*N; i ++)
    {
        char line_1[10];
        char line_adi[10];
        fgets(line_1, 10, fp_1);
        fgets(line_adi, 10, fp_adi);
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


    // Compare EXP with dt = 0.02 and EXP dt = 0.03
    fp_1 = fopen("last-exp-0.02.txt", "r");
    FILE *fp_exp_03;
    fp_exp_03 = fopen("last-exp-0.03.txt", "r");
    fprintf(fp, "Error (EXP with dt = 0.02 and EXP dt = 0.03):\n");
    sum = 0;
    sum_abs = 0;
    for(int i = 0; i < N*N; i ++)
    {
        char line_1[10];
        char line_exp_03[10];
        fgets(line_1, 10, fp_1);
        fgets(line_exp_03, 10, fp_exp_03);
        double num_1 = strtod(line_1, &ptr);
        double num_exp_03 = strtod(line_exp_03, &ptr);

        double p = pow(num_exp_03 - num_1, 2); 
        sum += p;
        sum_abs += fabs(num_exp_03 - num_1);
    }
    print_errors(fp, sum, sum_abs, norm_01);
    fclose(fp_exp_03);
    fclose(fp_1);


    // Compare EXP with dt = 0.02 and EXP dt = 0.04
    fp_1 = fopen("last-exp-0.02.txt", "r");
    FILE *fp_exp_04;
    fp_exp_04 = fopen("last-exp-0.04.txt", "r");
    fprintf(fp, "Error (EXP with dt = 0.02 and EXP dt = 0.04):\n");
    sum = 0;
    sum_abs = 0;
    for(int i = 0; i < N*N; i ++)
    {
        char line_1[10];
        char line_exp_04[10];
        fgets(line_1, 10, fp_1);
        fgets(line_exp_04, 10, fp_exp_04);
        double num_1 = strtod(line_1, &ptr);
        double num_exp_04 = strtod(line_exp_04, &ptr);

        double p = pow(num_exp_04 - num_1, 2); 
        sum += p;
        sum_abs += fabs(num_exp_04 - num_1);
    }
    print_errors(fp, sum, sum_abs, norm_01);
    fclose(fp_exp_04);
    fclose(fp_1);

    
    // Compare EXP with dt = 0.02 and EXP dt = 0.05
    fp_1 = fopen("last-exp-0.02.txt", "r");
    FILE *fp_exp_05;
    fp_exp_05 = fopen("last-exp-0.05.txt", "r");
    fprintf(fp, "Error (EXP with dt = 0.02 and EXP dt = 0.05):\n");
    sum = 0;
    sum_abs = 0;
    for(int i = 0; i < N*N; i ++)
    {
        char line_1[10];
        char line_exp_05[10];
        fgets(line_1, 10, fp_1);
        fgets(line_exp_05, 10, fp_exp_05);
        double num_1 = strtod(line_1, &ptr);
        double num_exp_05 = strtod(line_exp_05, &ptr);

        double p = pow(num_exp_05 - num_1, 2); 
        sum += p;
        sum_abs += fabs(num_exp_05 - num_1);
    }
    print_errors(fp, sum, sum_abs, norm_01);
    fclose(fp_exp_05);
    fclose(fp_1);


    fclose(fp);
    return 0;
}