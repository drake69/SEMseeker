# rm(list=ls())
# library(OpenCL)
# # create random phenotype two categories
# n_samples <- 5
# n_permutations <- 10
#
# phenotype <- sample(c(0,1),n_samples,replace=TRUE)
# burden <- rnorm(n_samples, mean = 0, sd = 1)
#
# # create a burden_permutated of n_permutations x n_samples where each n_permutation is a permutation of burden
# burden_permutated <- matrix(rep(burden, n_permutations), nrow = n_permutations, byrow = TRUE)
# burden_permutated <- t(apply(burden_permutated, 1, function(x) sample(x, length(x))))
#
# # check first n_permutation has all the values of burden
# # print(all(burden_permutated[2,] %in% burden))
#
# context <- OpenCL::oclContext(precision = "single")
#
# burden_buf <- OpenCL::as.clBuffer(as.single(burden), context, mode = "single")
# phenotype_buf <- OpenCL::as.clBuffer(phenotype, context, mode = "single")
# output_diffs_buf <- OpenCL::as.clBuffer(as.single(rep(0,n_permutations)), context, mode = "single")
#
# compute_permutation_diff_kernel_ok <- c("
#   __kernel void compute_row_means_diff(
#     __global float* mean_differences, // Output vector for the differences of means
#     const int n_permutations, // Number of n_permutations in the burden_permutated
#     const int n_samples, // Number of columns in the burden_permutated
#     __global const float* burden_permutated, // Input burden_permutated
#     __global const int* phenotype // Phenotype vector
# ) {
#     int n_permutation = get_global_id(0); // Global ID for the current n_permutation
#
#     // Check that the n_permutation ID is within bounds
#     if (n_permutation < n_permutations) {
#         float sum0 = 0.0f;
#         float sum1 = 0.0f;
#         int count0 = 0;
#         int count1 = 0;
#
#         // Calculate sums for the two phenotype categories
#         for (int n_sample = 0; n_sample < n_samples; n_sample++) {
#             int index = n_permutation + (n_permutations * n_sample); // Index in the burden_permutated
#             if (phenotype[n_sample] == 0) {
#                 sum0 += burden_permutated[index];
#                 count0++;
#             } else if (phenotype[n_sample] == 1) {
#                 sum1 += burden_permutated[index];
#                 count1++;
#             }
#         }
#
#         // Calculate means for the two categories
#         float mean0 = count0 > 0 ? sum0 / count0 : 0.0f;
#         float mean1 = count1 > 0 ? sum1 / count1 : 0.0f;
#
#         // Calculate the difference of the means and store it in the output vector
#         mean_differences[n_permutation] = mean0 - mean1;
#     }
# }
# ")
#
# compute_permutation_diff_kernel_obj_ok <- OpenCL::oclSimpleKernel(context, "compute_row_means_diff", compute_permutation_diff_kernel_ok, "single")
#
# f_ok <- function(compute_permutation_diff_kernel_obj, n_permutations,burden_buf, phenotype_buf,n_samples)
#   as.numeric(
#     OpenCL::oclRun(
#       kernel = compute_permutation_diff_kernel_obj,
#       dim = c(n_permutations), # First argument: number of n_permutations (vector)
#       as.integer(n_permutations), # Second argument: number of n_permutations (scalar)
#       as.integer(n_samples), # Fifth argument: number of columns (scalar)
#       OpenCL::as.clBuffer(as.single(burden), context, mode = "single"), # Third argument: buffer for the burden_permutated
#       OpenCL::as.clBuffer(as.integer(phenotype_buf), context, mode = "int") # Fourth argument: buffer for the phenotype vector
#     )
#   )
#
# # get execution time
# system.time(res <- f_ok(compute_permutation_diff_kernel_obj_ok, n_permutations,burden_buf, phenotype_buf,n_samples))
# res[1]
# system.time(res_ok <- apply(burden_permutated, 1, function(x) mean(x[phenotype == 0]) - mean(x[phenotype == 1])))
# res_ok[1]
# all.equal(res, res_ok)
#
#
