#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define SMALL_DATA_SIZE 32
#define MAX_STARS 50  // Maximum number of stars for histogram scaling

double variance(int large_data[], double mean_large, int large_data_size);
int generate_poisson(double lambda_hat);
void create_frequency_table(int data[], int size, int freq[], int *max_value);
void print_star_histogram(int freq[], int max_value, const char *title);
void print_gaussian_histogram(double gaussian_freq[], int max_value, const char *title);
double gaussian_pmf(int k, double lambda);

int main(){
    
    srand(time(NULL));
    long large_data_size;
    int small_data[SMALL_DATA_SIZE];
    double sum = 0;
    
    // User inputs
    char description[100];
    char output_filename[100];
    
    // Get user inputs
    printf("Enter description of the small dataset: ");
    fgets(description, sizeof(description), stdin);
    description[strcspn(description, "\n")] = 0;  // Remove newline
    
    do{
        printf("Enter an N for a large dataset version (10000 - 20000): ");
        scanf("%ld", &large_data_size);
    } while (large_data_size < 10000 || large_data_size > 20000);
    
    printf("Enter output filename (e.g., simulated.txt): ");
    scanf("%s", output_filename);
    
    // Clear input buffer after scanf
    int c;
    while ((c = getchar()) != '\n' && c != EOF);

    FILE *file, *file1, *file2;
    file = fopen("number_of_likes_per_min.csv", "r");

    if (file == NULL) {
        printf("Error opening file.\n");
        return 1;
    }
    
    printf("\n=== SMALL DATASET ===\n");
    printf("Description: %s\n", description);
    printf("Size: %d\n", SMALL_DATA_SIZE);
    printf("These are the number of likes per minute of a celebrity post\n");
    
    int i = 0;
    while (i < SMALL_DATA_SIZE && fscanf(file, "%d", &small_data[i]) == 1){
        printf("%d\n", small_data[i]);
        sum += small_data[i];
        i++;
    }
    fclose(file);
    
    /* CALCULATING THE LAMBDA */
    float lambda = sum / SMALL_DATA_SIZE;
    printf("\n=== SUMMARY STATISTICS ===\n");
    printf("Estimated parameter lambda = %.3f\n", lambda);
    
    /* OPENING ANOTHER FILE FOR WRITE FOR THE POISSON DATA OUTPUT */
    file1 = fopen(output_filename, "w");

    if (file1 == NULL) {
        printf("Error opening file.\n");
        return 1;
    }
    
    printf("\nGenerating %ld Poisson(%.3f) samples...\n", large_data_size, lambda);
    
    /* ACTUAL GENERATION OF LARGE DATA */
    for (int n = 0; n < large_data_size; n++){
        int poisson_value = generate_poisson(lambda);
        fprintf(file1, "%d\n", poisson_value);
    }
    fclose(file1);
    
    /* READING THE LARGE DATA SET THAT WAS CREATED BY KNUTH'S ALGORITHM */
    int large_data[large_data_size];
    long sum_large = 0;
    int o = 0;
    file2 = fopen(output_filename, "r");
    if (file2 == NULL) {
        printf("Error opening file.\n");
        return 1;
    }
    
    while (o < large_data_size && fscanf(file2, "%d", &large_data[o]) == 1){
        sum_large += large_data[o];
        o++;
    }
    fclose(file2);
    
    /* STATISTICAL ANALYSIS */
    double mean_large = sum_large / (float)large_data_size;
    double variance_large = variance(large_data, mean_large, large_data_size);
    
    printf("Sample mean of large dataset: %.3lf\n", mean_large);
    printf("Sample variance of large dataset: %.3lf\n", variance_large);
    
    /* FREQUENCY TABLE AND HISTOGRAM */
    printf("\n=== FREQUENCY TABLE ===\n");
    
    /* MAX VALUE */
    int max_value = 0;
    for (int j = 0; j < large_data_size; j++) {
        if (large_data[j] > max_value) {
            max_value = large_data[j];
        }
    }
    
    /* Create frequency table (size is max_value + 1 to include 0) */
    int *freq = (int *)calloc(max_value + 1, sizeof(int));
    if (freq == NULL) {
        printf("Memory allocation failed!\n");
        return 1;
    }
    
    for (int j = 0; j < large_data_size; j++) {
        freq[large_data[j]]++;
    }
    
    printf("\nFrequency Table:\n");
    printf("Value | Frequency\n");
    printf("------|----------\n");
    for (int k = 0; k <= max_value; k++) {
        if (freq[k] > 0) {
            printf("%5d | %9d\n", k, freq[k]);
        }
    }
    
    /* GAUSSIAN APPROXIMATION */
    printf("\n=== GAUSSIAN APPROXIMATION ===\n");
    
    // Create Gaussian expected frequencies
    double *gaussian_freq = (double *)calloc(max_value + 1, sizeof(double));
    if (gaussian_freq == NULL) {
        printf("Memory allocation failed!\n");
        free(freq);
        return 1;
    }
    
    // Calculate Gaussian PMF for each k and scale to dataset size
    double total_gaussian = 0.0;
    for (int k = 0; k <= max_value; k++) {
        gaussian_freq[k] = gaussian_pmf(k, lambda) * large_data_size;
        total_gaussian += gaussian_freq[k];
    }
    
    // Normalize to match total count
    double scale_factor = large_data_size / total_gaussian;
    for (int k = 0; k <= max_value; k++) {
        gaussian_freq[k] *= scale_factor;
    }
    
    /* PRINT BOTH HISTOGRAMS */
    printf("\n=== POISSON HISTOGRAM (from simulation) ===\n");
    print_star_histogram(freq, max_value, "Poisson");
    
    printf("\n=== GAUSSIAN HISTOGRAM (theoretical approximation) ===\n");
    print_gaussian_histogram(gaussian_freq, max_value, "Gaussian");
    
    // Free allocated memory
    free(freq);
    free(gaussian_freq);
    
    printf("\n=== PROGRAM COMPLETE ===\n");
    printf("Data saved to: %s\n", output_filename);
    
    return 0;
}

/* KNUTH'S ALGORITHM */
int generate_poisson(double lambda){
    double L = exp(-lambda);
    int k = 0;
    double p = 1.0;
    
    do {
        k++;
        p *= (double)rand() / RAND_MAX;
    } while (p > L);
    
    return k - 1;
}

/* CALCULATING THE VARIANCE */
double variance(int large_data[], double mean_large, int large_data_size){
    double sum_of_squares = 0.0;
    
    for (int p = 0; p < large_data_size; p++){
        double deviation = large_data[p] - mean_large;
        sum_of_squares += deviation * deviation; 
    }
    
    return sum_of_squares / large_data_size;
}

/* BUILD FREQUENCY TABLE */
void create_frequency_table(int data[], int size, int freq[], int *max_value) {
    // Initialize frequency array to zero
    for (int i = 0; i <= *max_value; i++) {
        freq[i] = 0;
    }
    
    // Count frequencies
    for (int i = 0; i < size; i++) {
        freq[data[i]]++;
    }
}

/* PRINT STAR HISTOGRAM WITH SCALING */
void print_star_histogram(int freq[], int max_value, const char *title) {
    // Find maximum frequency for scaling
    int max_freq = 0;
    for (int k = 0; k <= max_value; k++) {
        if (freq[k] > max_freq) {
            max_freq = freq[k];
        }
    }
    
    // Scale factor: divide by max_freq and multiply by MAX_STARS
    double scale_factor = (max_freq > MAX_STARS) ? (double)MAX_STARS / max_freq : 1.0;
    
    printf("\nValue | Histogram (%s)\n", title);
    printf("------|");
    for (int i = 0; i < MAX_STARS + 5; i++) printf("-");
    printf("\n");
    
    for (int k = 0; k <= max_value; k++) {
        if (freq[k] > 0) {
            printf("%5d | ", k);
            
            // Calculate number of stars (with scaling)
            int num_stars = (int)(freq[k] * scale_factor + 0.5); // +0.5 for rounding
            
            // Ensure at least 1 star for non-zero frequencies if scaling would make it 0
            if (num_stars == 0 && freq[k] > 0) {
                num_stars = 1;
            }
            
            // Print stars
            for (int s = 0; s < num_stars; s++) {
                printf("*");
            }
            
            // Show actual frequency in parentheses
            printf(" (%d)", freq[k]);
            printf("\n");
        }
    }
    
    // Print scaling information
    printf("\nScaling: Each star represents approximately %.1f counts\n", 
           (double)max_freq / MAX_STARS);
}

/* PRINT GAUSSIAN HISTOGRAM */
void print_gaussian_histogram(double gaussian_freq[], int max_value, const char *title) {
    // Find maximum frequency for scaling
    double max_freq = 0;
    for (int k = 0; k <= max_value; k++) {
        if (gaussian_freq[k] > max_freq) {
            max_freq = gaussian_freq[k];
        }
    }
    
    // Scale factor: divide by max_freq and multiply by MAX_STARS
    double scale_factor = (max_freq > MAX_STARS) ? (double)MAX_STARS / max_freq : 1.0;
    
    printf("\nValue | Histogram (%s)\n", title);
    printf("------|");
    for (int i = 0; i < MAX_STARS + 5; i++) printf("-");
    printf("\n");
    
    for (int k = 0; k <= max_value; k++) {
        if (gaussian_freq[k] > 0.1) {  // Only show values with significant probability
            printf("%5d | ", k);
            
            // Calculate number of stars (with scaling)
            int num_stars = (int)(gaussian_freq[k] * scale_factor + 0.5); // +0.5 for rounding
            
            // Ensure at least 1 star for non-zero frequencies if scaling would make it 0
            if (num_stars == 0 && gaussian_freq[k] > 0.5) {
                num_stars = 1;
            }
            
            // Print stars
            for (int s = 0; s < num_stars; s++) {
                printf("*");
            }
            
            // Show actual frequency in parentheses
            printf(" (%.1f)", gaussian_freq[k]);
            printf("\n");
        }
    }
    
    // Print scaling information
    printf("\nScaling: Each star represents approximately %.1f counts\n", 
           max_freq / MAX_STARS);
}

/* GAUSSIAN PMF APPROXIMATION FOR POISSON */
double gaussian_pmf(int k, double lambda) {
    if (lambda <= 0) return 0.0;
    
    double sigma = sqrt(lambda);
    double z1 = (k - 0.5 - lambda) / sigma;
    double z2 = (k + 0.5 - lambda) / sigma;
    
    // Using the error function (erf) to calculate standard normal CDF
    return 0.5 * (erf(z2 / sqrt(2.0)) - erf(z1 / sqrt(2.0)));
}