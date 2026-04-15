# VI
variation_of_information <- function(original, quantized) {

  original <- as.vector(original)
  quantized <- as.vector(quantized)

  entropy <- function(p) {
    p <- p[p > 0]  # Remove zeros to avoid log(0)
    -sum(p * log2(p))
  }

  joint_table <- table(original, quantized)
  joint_prob <- joint_table / sum(joint_table)
  original_prob <- margin.table(joint_prob, 1)
  quantized_prob <- margin.table(joint_prob, 2)

  H_original <- entropy(original_prob)
  H_quantized <- entropy(quantized_prob)
  H_joint <- entropy(joint_prob)

  # VI = H(X|Y) + H(Y|X) = 2*H(X,Y) - H(X) - H(Y)  (always >= 0)
  VI <- 2 * H_joint - H_original - H_quantized


  # # Function to calculate entropy
  # entropy <- function(p) {
  #   p <- p[p > 0]  # Remove zeros to avoid log(0)
  #   -sum(p * log2(p))
  # }

  # # Create a joint frequency table for the original and quantized distributions
  # joint_table <- table(original, quantized)
  #
  # # Convert joint frequency table to joint probability table
  # joint_prob <- joint_table / sum(joint_table)
  #
  # # Calculate marginal probabilities for the original and quantized distributions
  # original_prob <- margin.table(joint_prob, 1)
  # quantized_prob <- margin.table(joint_prob, 2)
  #
  # # Calculate entropy for original and quantized distributions
  # H_original <- entropy(original_prob)
  # H_quantized <- entropy(quantized_prob)
  #
  # # Calculate joint entropy by flattening the joint probability table
  # H_joint <- entropy(as.vector(joint_prob))
  #
  # # Calculate Variation of Information (VI)
  # VI <- H_original + H_quantized - 2 * H_joint

  # if(VI<0)
  #   

  # Return the calculated Variation of Information

  return(VI)
}
