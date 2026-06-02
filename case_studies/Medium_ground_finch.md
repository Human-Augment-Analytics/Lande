## Title: "Evolutionary Selection Analysis: Medium Ground Finch (Geospiza fortis)"

### 1. Introduction

Natural selection acts on phenotypic traits and determines survival under environmental stress. The medium ground finch (Geospiza fortis) population at El Garrapatero on Santa Cruz Island exhibits a well-documented bimodal distribution of beak size, providing an excellent system for studying disruptive selection and its role in adaptive radiation.

In this study, we analyze selection on beak size (PC1, derived from beak length, width, and depth) across multiple years (2003-2011) following the methods of Beausoleil et al. (2019). Our analysis focuses on:

- **Beak size (PC1)**: The primary axis of variation in beak morphology
- **Body size (PC.body1)**: A secondary trait representing overall body size

Univariate analyses were conducted for beak size across years to assess temporal variation in disruptive selection. Multivariate analyses incorporating body size were used to examine how selection operates on trait combinations.

---

### 2. Setup and Script
`R/scripts/test_birds.R`

---

### 3. Results

#### 3.1 Trait Correlation 

Before examining selection patterns, we assessed the correlation between beak size (Beak_PC1) and body size (PC.body1).

| Trait Pair | Correlation | Interpretation |
|------------|-------------|----------------|
| Beak_PC1 × PC.body1 | **0.816** | Strong positive correlation |

The strong positive correlation (r = 0.816) indicates the biological relationship between the traits: individuals with larger beaks also tend to have larger body size. This correlation is important for interpreting selection results, as it can lead to indirect selection where selection on one trait influences the apparent selection on the other.

---

#### 3.2 Summary of Disruptive Selection on Beak Size

| Year | N   | Survival_Rate | Gamma   | SE       | P_Value | Significant |
|------|-----|---------------|---------|----------|---------|-------------|
| 2003 | 49  | 0.408         | -0.0731 | 0.1271   | 0.4125  | FALSE       |
| 2004 | 110 | 0.345         | 0.0267  | 0.1028   | 0.8062  | FALSE       |
| 2005 | 185 | 0.276         | -0.0648 | 0.0635   | 0.2935  | FALSE       |
| 2006 | 233 | 0.172         | 0.0070  | 0.0403   | 0.8610  | FALSE       |
| 2007 | 61  | 0.328         | -0.0792 | 0.1381   | 0.5536  | FALSE       |
| 2008 | 127 | 0.291         | -0.0006 | 0.0658   | 0.8647  | FALSE       |
| 2009 | 196 | 0.173         | 0.1071  | 0.0390   | 0.0319  | TRUE        |
| 2010 | 175 | 0.120         | 0.0027  | 0.0508   | 0.7143  | FALSE       |
| 2011 | 211 | 1.000         | 0.0000  | 0.0000   | 1.0000  | FALSE       |

A total of 1,347 individuals were captured and recaptured across the 9 year study period, with annual sample sizes ranging from 49 to 233 individuals. Survival rates varied substantially among years, ranging from 0.120 (2010) to 0.408 (2003), reflecting the environmental variability characteristic of the Galapagos archipelago.

Significant disruptive selection on beak size was detected **only in 2009** (γ = 0.1071, P = 0.0319), indicating that individuals with extreme beak sizes (both small and large) had higher survival than those with intermediate beak sizes in that year. No other years showed significant quadratic selection (all P > 0.05), suggesting that disruptive selection on beak size is temporally variable in this population.

---

#### 3.3 Selection Coefficients

| Year | Trait | Linear (β) | P (Linear) | Quadratic (γ) | P (Quadratic) | Correlational (γ_ij) | P (Correlational) |
|------|-------|-----------|------------|---------------|---------------|---------------------|-------------------|
| 2009 | Beak_PC1 | 1.143 | 0.002 | -1.281 | 0.320 | 1.945 | 0.117 |
| 2009 | PC.body1 | -0.780 | 0.036 | -2.009 | 0.127 | - | - |
| 2006 | Beak_PC1 × PC.body1 | - | - | - | - | 1.146 | 0.099 |
| All others | All terms | - | > 0.05 | - | > 0.05 | - | > 0.05 |

**2009**: Significant directional selection on both traits - larger beaks and smaller body size favored
**No nonlinear selection** detected in any year (all quadratic P > 0.05)
**No correlational selection** detected in any year (all P > 0.05, with 2006 showing marginal non-significant trend)

*Complete coefficient tables for all years are provided in results folder.*

---


#### 3.4 Univariate Fitness on Beak Size

<table>
<tr>
<td><b>(A) 2003</b><br><img src="../results/medium_grouped_finch_results/figures/univariate_2003.png" width="100%"></td>
<td><b>(B) 2004</b><br><img src="../results/medium_grouped_finch_results/figures/univariate_2004.png" width="100%"></td>
<td><b>(C) 2005</b><br><img src="../results/medium_grouped_finch_results/figures/univariate_2005.png" width="100%"></td>
</tr>

<tr>
<td><b>(D) 2006</b><br><img src="../results/medium_grouped_finch_results/figures/univariate_2006.png" width="100%"></td>
<td><b>(E) 2007</b><br><img src="../results/medium_grouped_finch_results/figures/univariate_2007.png" width="100%"></td>
<td><b>(F) 2008</b><br><img src="../results/medium_grouped_finch_results/figures/univariate_2008.png" width="100%"></td>
</tr>

<tr>
<td><b>(G) 2009</b><br><img src="../results/medium_grouped_finch_results/figures/univariate_2009.png" width="100%"></td>
<td><b>(H) 2010</b><br><img src="../results/medium_grouped_finch_results/figures/univariate_2010.png" width="100%"></td>
<td><b>(I) 2011</b><br><img src="../results/medium_grouped_finch_results/figures/univariate_2011.png" width="100%"></td>
</tr>
</table>

The fitness function for 2009 shows statistically significant disruptive selection, indicating a non-linear relationship where extreme phenotypes (either small or large beaks) had a relative fitness advantage over intermediate beak sizes. This contrasts with other years where no significant non-linear selection was detected.

---

#### 3.5 Correlations Fitness 

<table>
<tr>
<td><b>(A) 2003</b><br><img src="../results/medium_grouped_finch_results/figures/cfs_2003.png" width="100%"></td>
<td><b>(B) 2004</b><br><img src="../results/medium_grouped_finch_results/figures/cfs_2004.png" width="100%"></td>

</tr>

<tr>
<td><b>(D) 2005</b><br><img src="../results/medium_grouped_finch_results/figures/cfs_2005.png" width="100%"></td>
<td><b>(E) 2006</b><br><img src="../results/medium_grouped_finch_results/figures/cfs_2006.png" width="100%"></td>
</tr>

<tr>
<td><b>(D) 2007</b><br><img src="../results/medium_grouped_finch_results/figures/cfs_2007.png" width="100%"></td>
<td><b>(E) 2008</b><br><img src="../results/medium_grouped_finch_results/figures/cfs_2008.png" width="100%"></td>

</tr>

<tr>
<td><b>(G) 2009</b><br><img src="../results/medium_grouped_finch_results/figures/cfs_2009.png" width="100%"></td>
<td><b>(H) 2010</b><br><img src="../results/medium_grouped_finch_results/figures/cfs_2010.png" width="100%"></td>

</tr>
</table>

---

#### 3.6 Adaptive Landscape 

<table>
<tr>
<td width="50%">
<b>(A) 2D Adaptive Landscape</b><br>
<img src="../results/medium_grouped_finch_results/figures/adaptive_landscape_2009.png" width="100%">
</td>

<td width="50%">
<b>(B) 3D Adaptive Landscape</b><br>
<img src="../results/medium_grouped_finch_results/figures/adaptive_landscape_3d_2009.png" width="100%">
</td>
</tr>
</table>


#### 3.7 Correlated Fitness Surface (Individual) vs. Adaptive Landscape (Population)

| Year | Trait 1   | Trait 2    | Correlation | n_points |
|------|-----------|------------|-------------|----------|
| 2009 | Beak_PC1  | PC.body1   | 0.9827      | 2500     |

*Note: The correlation of 0.9827 represents the high degree of similarity between the individual-level statistical surface (Correlated Fitness) and the population-level surface (Adaptive Landscape). This is distinct from the biological trait correlation (r = 0.816) reported in section 3.1.*


<table>
<tr>
<td width="40%" align="center">
<b>(A) Overlay</b><br>
<img src="../results/medium_grouped_finch_results/figures/comparison_2009_overlay.png" width="95%">
</td>

<td width="60%" align="center">
<b>(B) Side-by-Side</b><br>
<img src="../results/medium_grouped_finch_results/figures/comparison_2009_side_by_side.png" width="95%">
</td>
</tr>
</table>

---

### 4. Conclusion

**Disruptive selection was rare**: Only 2009 showed significant quadratic selection (γ = 0.1071, P = 0.0319)
**Directional selection dominated**: In 2009, larger beaks (β = 1.14) and smaller bodies (β = -0.78) were favored
**Strong trait correlation**: Beak size and body size were highly correlated (r = 0.82), complicating interpretation

**Comparison with Beausoleil et al. (2019):**
Our results are consistent with previous work in showing temporally variable selection. However, we detected weaker and less frequent disruptive selection, with only one significant year，this difference is more likely due to environmental variation or differences in the dataset rather than methodology.