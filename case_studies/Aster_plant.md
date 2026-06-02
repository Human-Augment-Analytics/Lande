## Title: "Evolutionary Selection Analysis: Aster Plant (Boechera stricta)"

### 1. Introduction

Understanding how natural selection shapes phenotypic variation across environmental gradients is central to evolutionary biology. The perennial forb Boechera stricta (Brassicaceae) exhibits genetically based clines in functional traits across elevation in the Rocky Mountains, providing an excellent system for studying adaptation and selection.

| Trait                    | Definition |
|--------------------------|------------|
| Specific leaf area (SLA) | Leaf morphology, related to resource acquisition |
| Water-use efficiency (C13) | Stable carbon isotope discrimination, integrated water-use efficiency |
| Snow depth (FDsnow)  | Winter snow accumulation, a key environmental factor |
| Height at flowering | Plant size at reproduction |

---

### 2. Setup and Script
`R/scripts/test_plants.R`

---

### 3. Results

#### 3.1 Data Summary

| Cohort | Sample Size (N) | Fitness Variable |
|--------|-----------------|------------------|
| 2011 | 66 | Year2012_fecund12 |
| 2012 | 87 | Year2013_fecund13 |
| Combined | 153 | relative_fitness |

I dropped NA
---

#### 3.2 Selection Differentials

Selection differentials (S) measure total selection on each trait, including both direct and indirect effects.

| Trait | S (2011) | S (2012) | S (Combined) | Direction | Interpretation |
|-------|----------|----------|--------------|-----------|----------------|
| SLA | +0.085 | -0.075 | +0.013 | Variable | Direction of selection changed between years |
| C13 | -0.006 | -0.002 | -0.005 | Negative | Weak negative selection on WUE |
| FDsnow | -0.185 | -0.052 | -0.126 | Negative | Deeper snow reduces fitness |
| Height | -0.203 | +0.074 | -0.078 | Variable | Taller plants favored in 2012 only |

**FDsnow** shows consistent negative selection across both years (S = -0.126 combined), suggesting that deeper winter snow reduces fitness, possibly through delayed spring emergence or increased pathogen pressure.
**Height** shows temporal variation: tall plants were selected against in 2011 (S = -0.203) but favored in 2012 (S = +0.074). This may reflect differences in drought conditions between years.
**SLA** also shows temporal variation, with selection switching from positive in 2011 to negative in 2012, resulting in no net selection when years are combined.
**C13** shows no evidence of selection in either year or the combined analysis.

---

#### 3.2 Univariate Fitness 

<table>
<tr>
<td><b>(A) TL</b><br><img src="../results/plant_results/figures/univariate_deltaC13.png" width="100%"></td>
<td><b>(B) HL</b><br><img src="../results/plant_results/figures/univariate_FDsnow.png" width="100%"></td>
</tr>

<tr>
<td><b>(C) WT</b><br><img src="../results/plant_results/figures/univariate_height.png" width="100%"></td>
<td><b>(D) KL</b><br><img src="../results/plant_results/figures/univariate_SLA.png" width="100%"></td>
</tr>
</table>


#### 3.3 Selection Coefficients

| Term                | Type           | Beta_Coefficient | Standard_Error | P_Value | Variance |
|---------------------|----------------|------------------|----------------|---------|----------|
| SLA                 | Linear         | -0.0563          | 0.0644         | 0.3834  | 0.00414  |
| deltaC13            | Linear         | 0.0322           | 0.0599         | 0.5918  | 0.00359  |
| FDsnow              | Linear         | -0.1466          | 0.0653         | 0.0271  | 0.00427  |
| height              | Linear         | -0.0695          | 0.0598         | 0.2482  | 0.00358  |
| SLA²                | Quadratic      | -0.0831          | 0.1235         | 0.5023  | 0.01524  |
| deltaC13²           | Quadratic      | -0.0590          | 0.1078         | 0.5856  | 0.01162  |
| FDsnow²             | Quadratic      | -0.0754          | 0.1436         | 0.6007  | 0.02063  |
| height²             | Quadratic      | -0.1266          | 0.1171         | 0.2829  | 0.01372  |
| SLA × deltaC13      | Correlational  | -0.0433          | 0.0911         | 0.6358  | 0.00829  |
| SLA × FDsnow        | Correlational  | -0.0630          | 0.0812         | 0.4399  | 0.00660  |
| SLA × height        | Correlational  | -0.0476          | 0.0707         | 0.5022  | 0.00500  |
| deltaC13 × FDsnow   | Correlational  | 0.0414           | 0.0753         | 0.5835  | 0.00567  |
| deltaC13 × height   | Correlational  | 0.0321           | 0.0738         | 0.6643  | 0.00545  |
| FDsnow × height     | Correlational  | -0.1624          | 0.0934         | 0.0854  | 0.00872  |

No significant linear or quadratic selection gradients were detected in the combined analysis (all P > 0.05). This may reflect the relatively small sample size (n = 153) compared to the original study (n = 4,510).

---

#### 3.4 Disruptive Selection

| Trait | Linear (β) | Quadratic (γ) | P (Linear) | P (Quadratic) | Interpretation |
|-------|-----------|---------------|------------|---------------|----------------|
| SLA | -0.045 | -0.268 | 0.612 | 0.132 | No nonlinear selection |
| C13 | 0.067 | -0.178 | 0.468 | 0.334 | No nonlinear selection |
| FDsnow | -0.023 | -0.090 | 0.767 | 0.562 | No nonlinear selection |
| Height | 0.089 | 0.134 | 0.291 | 0.423 | No nonlinear selection |

Consistent with the multivariate analysis, no significant quadratic selection was detected for any trait in the combined dataset.

---

#### 3.5 Correlations Fitness 

<table>
<tr>
<td><b>(A) C13 × FDsnow</b><br><img src="../results/plant_results/figures/cfs_deltaC13_FDsnow_enhanced.png" width="100%"></td>
<td><b>(B) C13 × Height</b><br><img src="../results/plant_results/figures/cfs_deltaC13_height_enhanced.png" width="100%"></td>
<td><b>(C) FDsnow × Height</b><br><img src="../results/plant_results/figures/cfs_FDsnow_height_enhanced.png" width="100%"></td>
</tr>

<tr>
<td><b>(D) SLA × C13</b><br><img src="../results/plant_results/figures/cfs_SLA_deltaC13_enhanced.png" width="100%"></td>
<td><b>(E) SLA x Height</b><br><img src="../results/plant_results/figures/cfs_SLA_height_enhanced.png" width="100%"></td>
<td><b>(F) SLA × FDsnow</b><br><img src="../results/plant_results/figures/cfs_SLA_FDsnow_enhanced.png" width="100%"></td>
</tr>

</table>

---

#### 3.6 Adaptive Landscape


<table>
<tr>
<td width="50%" align="center">
<b>(A) deltaC13 × FDsnow (2D)</b><br>
<img src="../results/plant_results/figures/adaptive_landscape_deltaC13_FDsnow_2d.png" width="95%">
</td>

<td width="50%" align="center">
<b>(B) deltaC13 × FDsnow (3D)</b><br>
<img src="../results/plant_results/figures/adaptive_landscape_deltaC13_FDsnow_3d.png" width="95%">
</td>
</tr>
</table>

<table>
<tr>
<td width="50%" align="center">
<b>(C) deltaC13 × height (2D)</b><br>
<img src="../results/plant_results/figures/adaptive_landscape_deltaC13_height_2d.png" width="95%">
</td>

<td width="50%" align="center">
<b>(D) deltaC13 × height (3D)</b><br>
<img src="../results/plant_results/figures/adaptive_landscape_deltaC13_height_3d.png" width="95%">
</td>
</tr>
</table>

<table>
<tr>
<td width="50%" align="center">
<b>(E) FDsnow × height (2D)</b><br>
<img src="../results/plant_results/figures/adaptive_landscape_FDsnow_height_2d.png" width="95%">
</td>

<td width="50%" align="center">
<b>(F) FDsnow × height (3D)</b><br>
<img src="../results/plant_results/figures/adaptive_landscape_FDsnow_height_3d.png" width="95%">
</td>
</tr>
</table>

<table>
<tr>
<td width="50%" align="center">
<b>(G) SLA × deltaC13 (2D)</b><br>
<img src="../results/plant_results/figures/adaptive_landscape_SLA_deltaC13_2d.png" width="95%">
</td>

<td width="50%" align="center">
<b>(H) SLA × deltaC13 (3D)</b><br>
<img src="../results/plant_results/figures/adaptive_landscape_SLA_deltaC13_3d.png" width="95%">
</td>
</tr>
</table>

<table>
<tr>
<td width="50%" align="center">
<b>(I) SLA × FDsnow (2D)</b><br>
<img src="../results/plant_results/figures/adaptive_landscape_SLA_FDsnow_2d.png" width="95%">
</td>

<td width="50%" align="center">
<b>(J) SLA × FDsnow (3D)</b><br>
<img src="../results/plant_results/figures/adaptive_landscape_SLA_FDsnow_3d.png" width="95%">
</td>
</tr>
</table>

<table>
<tr>
<td width="50%" align="center">
<b>(K) SLA × height (2D)</b><br>
<img src="../results/plant_results/figures/adaptive_landscape_SLA_height_2d.png" width="95%">
</td>

<td width="50%" align="center">
<b>(L) SLA × height (3D)</b><br>
<img src="../results/plant_results/figures/adaptive_landscape_SLA_height_3d.png" width="95%">
</td>
</tr>
</table>


---

### 4.0 Discussion

Our analysis of selection on functional traits in Boechera stricta revealed several important patterns:

**1. Weak directional selection.** Selection differentials and coefficients were generally weak and non-significant across all traits. This may reflect the relatively small sample size (n = 153) compared to the original study (n = 4,510), which limited our statistical power to detect selection.

**2. No evidence of nonlinear selection.** Neither the univariate disruptive analysis nor the multivariate quadratic terms showed significant curvature. This contrasts with the original study (Wadgymar et al. 2017), which found significant stabilizing selection when integrating data across years.

**3. Comparison with Wadgymar et al. (2017).**

| Finding | Wadgymar et al. (2017) | Our Study | Consistency |
|---------|------------------------|-----------|-------------|
| Sample size | 4,510 | 153 | Different |
| Stabilizing selection | Significant | Not significant | Not replicated |
| Directional selection | Weak | Weak | Consistent |
| Temporal variation | Important | Not tested | N/A |

---

