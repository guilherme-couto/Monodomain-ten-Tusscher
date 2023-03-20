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
    fp_1 = fopen("tnnp-exp-0.02.txt", "r");
    FILE *fp_adi;
    fp_adi = fopen("tnnp-adi-0.02.txt", "r");
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


    // Compare EXP with dt = 0.02 and EXP dt = 0.05
    fp_1 = fopen("tnnp-exp-0.02.txt", "r");
    FILE *fp_exp_05;
    fp_exp_05 = fopen("tnnp-exp-0.05.txt", "r");
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


    // Compare EXP with dt = 0.02 and EXP dt = 0.1
    fp_1 = fopen("tnnp-exp-0.02.txt", "r");
    FILE *fp_exp_1;
    fp_exp_1 = fopen("tnnp-exp-0.1.txt", "r");
    fprintf(fp, "Error (EXP with dt = 0.02 and EXP dt = 0.1):\n");
    sum = 0;
    sum_abs = 0;
    for(int i = 0; i < N*N; i ++)
    {
        char line_1[10];
        char line_exp_1[10];
        fgets(line_1, 10, fp_1);
        fgets(line_exp_1, 10, fp_exp_1);
        double num_1 = strtod(line_1, &ptr);
        double num_exp_1 = strtod(line_exp_1, &ptr);

        double p = pow(num_exp_1 - num_1, 2); 
        sum += p;
        sum_abs += fabs(num_exp_1 - num_1);
    }
    print_errors(fp, sum, sum_abs, norm_01);
    fclose(fp_exp_1);
    fclose(fp_1);


    // Compare EXP with dt = 0.02 and ADI dt = 0.03
    fp_1 = fopen("tnnp-exp-0.02.txt", "r");
    FILE *fp_3;
    fp_3 = fopen("tnnp-adi-0.03.txt", "r");
    fprintf(fp, "Error (EXP with dt = 0.02 and ADI dt = 0.03):\n");
    sum = 0;
    sum_abs = 0;
    for(int i = 0; i < N*N; i ++)
    {
        char line_1[10];
        char line_3[10];
        fgets(line_1, 10, fp_1);
        fgets(line_3, 10, fp_3);
        double num_1 = strtod(line_1, &ptr);
        double num_3 = strtod(line_3, &ptr);

        double p = pow(num_3 - num_1, 2); 
        sum += p;
        sum_abs += fabs(num_3 - num_1);
    }
    print_errors(fp, sum, sum_abs, norm_01);
    fclose(fp_3);
    fclose(fp_1);
    

    // Compare EXP with dt = 0.02 and ADI dt = 0.04
    fp_1 = fopen("tnnp-exp-0.02.txt", "r");
    FILE *fp_4;
    fp_4 = fopen("tnnp-adi-0.04.txt", "r");
    fprintf(fp, "Error (EXP with dt = 0.02 and ADI dt = 0.04):\n");
    sum = 0;
    sum_abs = 0;
    for(int i = 0; i < N*N; i ++)
    {
        char line_1[10];
        char line_4[10];
        fgets(line_1, 10, fp_1);
        fgets(line_4, 10, fp_4);
        double num_1 = strtod(line_1, &ptr);
        double num_4 = strtod(line_4, &ptr);

        double p = pow(num_4 - num_1, 2); 
        sum += p;
        sum_abs += fabs(num_4 - num_1);
    }
    print_errors(fp, sum, sum_abs, norm_01);
    fclose(fp_4);
    fclose(fp_1);


    // Compare EXP with dt = 0.02 and ADI dt = 0.05
    fp_1 = fopen("tnnp-exp-0.02.txt", "r");
    FILE *fp_5;
    fp_5 = fopen("tnnp-adi-0.05.txt", "r");
    fprintf(fp, "Error (EXP with dt = 0.02 and ADI dt = 0.05):\n");
    sum = 0;
    sum_abs = 0;
    for(int i = 0; i < N*N; i ++)
    {
        char line_1[10];
        char line_5[10];
        fgets(line_1, 10, fp_1);
        fgets(line_5, 10, fp_5);
        double num_1 = strtod(line_1, &ptr);
        double num_5 = strtod(line_5, &ptr);

        double p = pow(num_5 - num_1, 2); 
        sum += p;
        sum_abs += fabs(num_5 - num_1);
    }
    print_errors(fp, sum, sum_abs, norm_01);
    fclose(fp_5);
    fclose(fp_1);


    // Compare EXP with dt = 0.02 and ADI dt = 0.1
    fp_1 = fopen("tnnp-exp-0.02.txt", "r");
    FILE *fp_10;
    fp_10 = fopen("tnnp-adi-0.1.txt", "r");
    fprintf(fp, "Error (EXP with dt = 0.02 and ADI dt = 0.1):\n");
    sum = 0;
    sum_abs = 0;
    for(int i = 0; i < N*N; i ++)
    {
        char line_1[10];
        char line_10[10];
        fgets(line_1, 10, fp_1);
        fgets(line_10, 10, fp_10);
        double num_1 = strtod(line_1, &ptr);
        double num_10 = strtod(line_10, &ptr);

        double p = pow(num_10 - num_1, 2); 
        sum += p;
        sum_abs += fabs(num_10 - num_1);
    }
    print_errors(fp, sum, sum_abs, norm_01);
    fclose(fp_10);
    fclose(fp_1);


    // Compare EXP with dt = 0.02 and ADI dt = 0.2
    fp_1 = fopen("tnnp-exp-0.02.txt", "r");
    FILE *fp_20;
    fp_20 = fopen("tnnp-adi-0.2.txt", "r");
    fprintf(fp, "Error (EXP with dt = 0.02 and ADI dt = 0.2):\n");
    sum = 0;
    sum_abs = 0;
    for(int i = 0; i < N*N; i ++)
    {
        char line_1[10];
        char line_20[10];
        fgets(line_1, 10, fp_1);
        fgets(line_20, 10, fp_20);
        double num_1 = strtod(line_1, &ptr);
        double num_20 = strtod(line_20, &ptr);

        double p = pow(num_20 - num_1, 2); 
        sum += p;
        sum_abs += fabs(num_20 - num_1);
    }
    print_errors(fp, sum, sum_abs, norm_01);
    fclose(fp_20);
    fclose(fp_1);


    // Compare EXP with dt = 0.02 and EXP dt = 0.005
    fp_1 = fopen("tnnp-exp-0.02.txt", "r");
    FILE *fp_05;
    fp_05 = fopen("tnnp-exp-0.005.txt", "r");
    fprintf(fp, "Error (EXP with dt = 0.02 and EXP dt = 0.005):\n");
    sum = 0;
    sum_abs = 0;
    for(int i = 0; i < N*N; i ++)
    {
        char line_1[10];
        char line_05[10];
        fgets(line_1, 10, fp_1);
        fgets(line_05, 10, fp_05);
        double num_1 = strtod(line_1, &ptr);
        double num_05 = strtod(line_05, &ptr);

        double p = pow(num_05 - num_1, 2); 
        sum += p;
        sum_abs += fabs(num_05 - num_1);
    }
    print_errors(fp, sum, sum_abs, norm_01);
    fclose(fp_05);
    fclose(fp_1);




    fclose(fp);
    return 0;
}