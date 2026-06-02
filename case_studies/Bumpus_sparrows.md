
## Title: "Evolutionary Selection Analysis: Bumpus Sparrows"

### 1. Introduction

Natural selection acts on phenotypic traits and determines survival under environmental stress.  
The Bumpus sparrow dataset provides a classical example of viability selection following a severe winter storm.

In this study, we analyze selection on four morphological traits:

- Total length (TL)
- Humerus length (HL)
- Body weight (WT)
- Keel length (KL)

Total length (TL) was used as a primary trait for comparison with classical analyses, 
as previous studies have demonstrated strong directional selection on this trait.

Univariate analyses were conducted for each trait (TL, HL, WT, KL) 
to assess their independent effects on survival.

Multivariate analyses were then performed by progressively adding traits 
to examine how selection operates on combinations of traits.

In addition, Type III analyses were used to evaluate the contribution of each trait 
while controlling for the effects of other traits in the model.

---

### 2. Setup and Script
`R/tests/test_bumpus_sparrows.R`
 
---

### 3. Results

#### 3.1 Selection Differentials 

| Trait | S (Selection Differential) | Direction | Interpretation |
|------|----------------------------|----------|----------------|
| TL   | -0.235 | Negative | Larger individuals have lower survival probability |
| HL   | 0.163  | Positive | Longer humerus length increases survival |
| WT   | -0.202 | Negative | Heavier individuals have lower survival |
| KL   | 0.128  | Positive | Larger keel length is associated with higher survival |

Selection differentials indicate that survival favors smaller body size (TL, WT) and larger structural traits (HL, KL).

Unlike the logistic regression results presented in the lecture, which use raw trait values for interpretability, 
selection differentials in this study are calculated using standardized traits. This allows for direct comparison of selection strength across traits.

---

#### 3.2 Selection Coefficients 

| Term | Type | Beta_Coefficient | Standard_Error | Variance | P_Value |
|------|------|------------------|----------------|----------|--------|
| TL   | Linear | -1.142 | 0.309 | 0.095 | 0.000 |
| HL   | Linear | 1.177  | 0.311 | 0.097 | 0.000 |
| WT   | Linear | -1.023 | 0.307 | 0.094 | 0.001 |
| KL   | Linear | 0.869  | 0.295 | 0.087 | 0.003 |
| TL²  | Quadratic | -0.101 | 0.736 | 0.542 | 0.891 |
| HL²  | Quadratic | -0.337 | 0.674 | 0.454 | 0.617 |
| WT²  | Quadratic | -0.532 | 0.691 | 0.478 | 0.441 |
| KL²  | Quadratic | 0.582  | 0.739 | 0.546 | 0.430 |
| TL × HL | Correlational | 0.368 | 0.492 | 0.242 | 0.454 |
| TL × WT | Correlational | -0.061 | 0.557 | 0.310 | 0.913 |
| TL × KL | Correlational | -1.018 | 0.598 | 0.358 | 0.089 |
| HL × WT | Correlational | 0.348 | 0.574 | 0.330 | 0.545 |
| HL × KL | Correlational | -0.073 | 0.435 | 0.190 | 0.867 |
| WT × KL | Correlational | 0.351 | 0.525 | 0.275 | 0.503 |

**Directional selection** was strong and consistent: body size traits (TL, WT) showed negative selection (β ≈ -1.1), while structural traits (HL, KL) showed positive selection (β ≈ 1.2 and 0.9, respectively). All linear coefficients were significant (P < 0.01), indicating robust directional selection.

**Nonlinear selection (quadratic terms)** was weak and not statistically significant (all P > 0.4), suggesting the fitness surface is primarily linear with no evidence of stabilizing or disruptive selection.

**Correlational selection** among traits was generally weak (all P > 0.09), with only marginal evidence for an interaction between total length and keel length (P = 0.089). This suggests trait combinations do not strongly influence survival in this dataset.

---

#### 3.3 Disruptive (Nonlinear) Selection 

| Trait | Linear (β) | Quadratic (γ) | p (Linear) | p (Quadratic) | Interpretation |
|------|-----------|---------------|------------|---------------|----------------|
| TL   | -0.609 | -0.829 | 0.002 | 0.020 | Stabilizing selection |
| HL   | 0.297  | -0.374 | 0.116 | 0.203 | No nonlinear selection |
| WT   | -0.416 | -0.318 | 0.036 | 0.322 | No nonlinear selection |
| KL   | 0.277  | -0.117 | 0.118 | 0.685 | No nonlinear selection |

Univariate analysis of quadratic selection gradients revealed significant stabilizing selection for total length (TL: γ = -0.83, P = 0.02), indicating that intermediate-sized individuals had the highest survival probability. This finding aligns with classic studies of Bumpus sparrows, which documented higher mortality among extreme-sized individuals.

No other traits showed significant quadratic selection (all P > 0.20), suggesting that nonlinear selection is generally weak or absent for the remaining morphological traits in this dataset. Directional selection patterns (linear coefficients) were consistent across traits: body size traits (TL, WT) showed negative directional selection, while structural traits (HL, KL) showed positive directional selection.

---

#### 3.4 Univariate Fitness 

<table>
<tr>
<td><b>(A) TL</b><br><img src="../results/bumpus_sparrows_results/figures/univariate_TL.png" width="100%"></td>
<td><b>(B) HL</b><br><img src="../results/bumpus_sparrows_results/figures/univariate_HL.png" width="100%"></td>
</tr>

<tr>
<td><b>(C) WT</b><br><img src="../results/bumpus_sparrows_results/figures/univariate_WT.png" width="100%"></td>
<td><b>(D) KL</b><br><img src="../results/bumpus_sparrows_results/figures/univariate_KL.png" width="100%"></td>
</tr>
</table>

---

#### 3.5 Correlations Fitness 

<table>
<tr>
<td><b>(A) TL × HL</b><br><img src="../results/bumpus_sparrows_results/figures/cfs_TL_HL_enhanced.png" width="100%"></td>
<td><b>(B) TL × WT</b><br><img src="../results/bumpus_sparrows_results/figures/cfs_TL_WT_enhanced.png" width="100%"></td>
<td><b>(C) TL × KL</b><br><img src="../results/bumpus_sparrows_results/figures/cfs_TL_KL_enhanced.png" width="100%"></td>
</tr>

<tr>
<td><b>(D) HL × WT</b><br><img src="../results/bumpus_sparrows_results/figures/cfs_HL_WT_enhanced.png" width="100%"></td>
<td><b>(E) HL × KL</b><br><img src="../results/bumpus_sparrows_results/figures/cfs_HL_KL_enhanced.png" width="100%"></td>
<td><b>(F) WT × KL</b><br><img src="../results/bumpus_sparrows_results/figures/cfs_WT_KL_enhanced.png" width="100%"></td>
</tr>

<tr>
<td colspan="3" align="center"><i>All surfaces show largely linear gradients with minimal curvature.</i></td>
</tr>
</table>

Fitness increases with HL and decreases with TL, showing opposing directional selection with no strong interaction.

Fitness is highest at intermediate TL and lower WT, suggesting weak curvature but primarily directional selection.

Fitness increases with KL and decreases with TL, consistent with directional selection and weak interaction effects.

Fitness increases with HL and decreases with WT, with a smooth gradient indicating additive effects of both traits.

Fitness increases with both HL and KL, forming a planar surface consistent with strong positive selection on structural traits.

Fitness increases with KL and decreases with WT, showing a clear directional gradient with minimal interaction.

---


#### 3.6 Adaptive Landscape 

| Metric | Value |
|--------|------|
| Optimal TL | -2.961 |
| Optimal KL | 3.169 |
| Optimal Fitness | 0.995 |

<table>
<tr>
<td width="50%">
<b>(A) 2D Adaptive Landscape</b><br>
<img src="../results/bumpus_sparrows_results/figures/adaptive_landscape_2d.png" width="100%">
</td>

<td width="50%">
<b>(B) 3D Adaptive Landscape</b><br>
<img src="../results/bumpus_sparrows_results/figures/adaptive_landscape_3d.png" width="100%">
</td>
</tr>
</table>

The adaptive landscape was visualized using the two most functionally distinct traits:
total length (TL) representing overall body size, and keel length (KL) representing skeletal/flying apparatus structure. These traits showed the strongest directional selection in the multivariate analysis (β_TL = -1.14, β_TL = 0.87) and represent different morphological axes.

---

#### 3.7 Correlated Fitness Surface (Individual) vs. Adaptive Landscape (Population)

| Trait Pair | Correlation | N (Points) |
|-----------|------------|------------|
| TL × KL | 0.969 | 2500 |

For the Bumpus sparrow dataset, comparing TL and HL surfaces revealed 
a correlation of 0.969, indicating that populations with 
mean phenotypes favored by individual selection have the highest 
mean fitness. This strong agreement supports the use of individual 
selection gradients to predict evolutionary change。


<table>
<tr>
<td width="40%" align="center">
<b>(A) Overlay</b><br>
<img src="../results/bumpus_sparrows_results/figures/comparison_TL_KL_overlay.png" width="95%">
</td>

<td width="60%" align="center">
<b>(B) Side-by-Side</b><br>
<img src="../results/bumpus_sparrows_results/figures/comparison_TL_KL_side_by_side.png" width="95%">
</td>
</tr>
</table>


The overlay (A) shows the alignment between individual-level and population-level fitness patterns, while the side-by-side view (B) highlights their structural similarity. 