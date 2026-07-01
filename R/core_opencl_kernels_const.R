compute_permutation_diff_kernel <- c("
      __kernel void compute_row_means_diff(
        __global float* mean_differences, // Output vector for the differences of means
        const int n_permutations, // Number of n_permutations in the independent_var_permutated
        const int n_samples, // Number of columns in the independent_var_permutated
        const int n_areas, // Number of areas
        __global const float* burden, // Input independent_var_permutated
        __global const int* independent_var_permutated // burden vector
    ) {
        int n_permutation = get_global_id(0); // Global ID for the current n_permutation

        // Check that the n_permutation ID is within bounds
        if (n_permutation < n_permutations) {

          for(int n_area = 0; n_area < n_areas; n_area++) {

              // Initialize sums and counts for the two burden_of_area categories
              float sum0 = 0.0f;
              float sum1 = 0.0f;
              int count0 = 0;
              int count1 = 0;
              float burden_of_area = 0.0f;

              // Calculate sums for the two burden_of_area categories
              for (int n_sample = 0; n_sample < n_samples; n_sample++) {

                  burden_of_area = burden[(n_area + n_samples) + n_sample];
                  int index = n_permutation + (n_permutations * n_sample); // Index in the independent_var_permutated;
                  if ( independent_var_permutated[index] == 0) {
                      sum0 += burden_of_area;
                      count0++;
                  } else if (independent_var_permutated[index] == 1) {
                      sum1 += burden_of_area;
                      count1++;
                  }
              }

              // Calculate means for the two categories
              float mean0 = count0 > 0 ? sum0 / count0 : 0.0f;
              float mean1 = count1 > 0 ? sum1 / count1 : 0.0f;

              // Calculate the difference of the means and store it in the output vector
              mean_differences[n_permutation * n_area] = mean1 - mean0;
          }
        }
    }
  ")
