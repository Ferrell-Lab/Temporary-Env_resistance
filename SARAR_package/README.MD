# SARAR: Statistical Analysis of Perturbation Response in R

The `SARAR` package provides an implementation of the Statistical Analysis of Perturbation Response (SARA), a method designed to quantify the response of single cells to stimulation by comparing the distribution of phosphoprotein intensities between basal and stimulated conditions. This is an R implementation of the method used by Levine, et al. in the article ["Data-Driven Phenotypic Dissection of AML Reveals Progenitor-like Cells that Correlate with Prognosis, Cell 2015"](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4508757/).

## Installation

To install the package from source, follow these steps:

```r
# Install devtools if you don't have it
install.packages("devtools")

# Install the SARA package
devtools::install("path_to_SARA_package")

# Example basal and stimulated values
basal <- rnorm(100)
stimulated <- rnorm(100, mean = 0.5)

# Calculate the SARA score
sara_score <- sara(basal, stimulated)
print(sara_score)

```

## Detailed Steps of the SARA Method

The`SARAR` package implements an updated version of the original SARA method, focused on comparing distributions between two conditions (e.g., unstimulated and stimulated) using normalized histograms and Earth Mover’s Distance (EMD). This method provides a robust approach for detecting and quantifying differences between these conditions.

### Key Adjustments:
- **Normalized Histograms**: Instead of using empirical cumulative distribution functions (ECDFs), we compute and compare normalized histograms of the two conditions. This ensures that the area under each distribution sums to 1, providing a consistent basis for comparison.
- **Earth Mover’s Distance (EMD)**: We calculate the EMD using the emdR() function, which measures how much distribution mass must be moved to transform one distribution into the other. This method is particularly well-suited for comparing distributions of different shapes.
  
SARAR applies the following steps:  

### 1. Create Normalized Histograms

- For both conditions (e.g., unstimulated and stimulated), histograms are created with the same binning strategy.
- The data are normalized such that the total density across all bins sums to 1. This normalization step ensures that the histograms can be meaningfully compared, even if the two conditions have different numbers of observations.
  
### 2. Compute Earth Mover's Distance (EMD)

- The emdR() function is used to calculate the EMD between the two normalized histograms.
- EMD measures the minimum amount of "work" required to move the mass from one histogram (unstimulated) to the other (stimulated), providing a quantitative measure of distributional differences.

### 3. Permutation Test for Statistical Significance

- A permutation test is performed to assess whether the observed EMD is statistically significant.
- The data from both conditions are randomly shuffled, and the EMD is recalculated for each permutation. This process is repeated many times (e.g., 1000 permutations).
- The p-value is calculated as the proportion of permuted EMD values that are greater than or equal to the observed EMD. This p-value indicates the likelihood that the observed difference between the distributions is due to random chance.

### 4. Calculate the Final SARA Score

- Compute the difference between the means of the two distributions (basal and stimulated).
- The final score is given by the formula:

  $` \text{score} = \text{EMD} \times \text{sign}(\mathbb{E}[\phi_s] - \mathbb{E}[\phi_b]) \times (1 - p) `$

  Where:
  - $` \mathbb{E}[\phi_s] `$ is the mean of the stimulated condition.
  - $` \mathbb{E}[\phi_b] `$ is the mean of the basal condition.
  - $` p `$ is the permutation p-value.

### Interpreting the SARA Score

The `sara_score` provides a quantitative measure of the difference between the distributions of the basal and stimulated conditions, with respect to both magnitude and statistical significance. Here's how to interpret the score:

#### Components of the SARA Score

1. **Earth Mover's Distance (EMD):**:
   - **Magnitude of difference**: EMD quantifies how different the two distributions are. A higher EMD value means the distributions are more different, while a lower EMD value means they are more similar.

2. **Sign of the Difference** (`sign(mean_diff)`):
   - The sign component indicates whether the stimulated condition has a higher or lower mean compared to the basal condition:
     - **Positive `sara_score`**: The stimulated distribution tends to have higher values than the basal distribution.
     - **Negative `sara_score`**: The stimulated distribution tends to have lower values than the basal distribution.

3. **Statistical Significance (`1 - p_value`)**:
   - The p-value from the permutation test gives a measure of the likelihood that the observed difference is due to random chance. A smaller p-value (near 0) means the difference is statistically significant, while a larger p-value (near 1) indicates that the difference may not be significant.
   - The final `sara_score` is scaled by the significance of the observed difference: differences that are not statistically significant will result in a smaller score (closer to 0).

#### Interpreting the SARA Score

- **Large positive `sara_score`**: Indicates a strong, statistically significant increase in the stimulated condition compared to the basal condition.
- **Large negative `sara_score`**: Indicates a strong, statistically significant decrease in the stimulated condition compared to the basal condition.
- **Near-zero `sara_score`**: Either there is no significant difference between the basal and stimulated conditions, or any observed differences are not statistically significant.
  
To sum, the sign of the score tells you the direction of the response, while the magnitude reflects both the size of the difference and its statistical significance.

#### Example Interpretations

- A **`sara_score` of +2.5**: The stimulated condition has a higher distribution compared to the basal condition, and this difference is statistically significant.
- A **`sara_score` of -1.0**: The stimulated condition has a lower distribution compared to the basal condition, but the difference is likely not very statistically significant.
- A **`sara_score` of 0.1**: The difference between basal and stimulated conditions is small and statistically insignificant.

#### Earth Mover's Distance `emdR()` function:
- The `emdR()` uses the Hungarian algorithm-inspired method where $` d_0 = 0 `$ and iteratively apply $` d_{i+1} = f_i + d_i - m_i `$. 
- The Hungarian algorithm approach calculates the difference between the cumulative sums (via $` d_i `$) of the two distributions. The EMD is then the sum of the absolute values of these differences.
