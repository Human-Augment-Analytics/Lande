## Title: "Evolutionary Selection Analysis: Caribbean Pupfish"

### 1. Introduction

The adaptive radiation of pupfishes (*Cyprinodon* spp.) on San Salvador Island, Bahamas, represents one of the most rapid and spectacular examples of trophic diversification in vertebrates. Within the past 10,000 years, three sympatric species have evolved from a generalist ancestor: the generalist (*C. variegatus*), the scale-eater (*C. desquamator*), and the molluscivore (*C. brontotheroides*). These species exhibit extreme morphological divergence, particularly in jaw length, nasal protrusion, and body depth.

In this study, we analyze selection on six functional traits using a large F2 hybrid population (n = 1,347 individuals) from a previous field experiment (Martin 2016). We measure two fitness components:
- **Survival** (binary): Whether individuals survived the 3-month field enclosure period
- **Growth** (continuous): Log-transformed growth rate of surviving individuals


| Trait | Definition |
|-------|------------|
| jaw   | Length or size of the jaw, related to feeding ability and bite reach |
| eye   | Size of the eye, associated with vision and prey detection |
| body  | Body depth or shape, influencing swimming performance and stability |
| nasal | Degree of nasal protrusion, linked to specialized feeding strategies |
| mouth | Size or shape of the mouth opening, affecting food intake |
| SL    | Standard length (overall body size from snout to tail base) |

---

### 2. Setup and Script
`R/scripts/test_fish.R`

---

### 3. Results

#### 3.1 Selection Differentials

| Trait | Binary (Survival) | Continuous (Growth) | Interpretation |
|-------|------------------|---------------------|----------------|
| jaw | +0.057 | +0.325 | Jaw length positively affects growth |
| eye | +0.004 | -0.208 | Eye diameter negatively affects growth |
| body | -0.054 | -0.107 | Body depth negatively affects both |
| nasal | -0.046 | -0.211 | Nasal protrusion negatively affects growth |
| mouth | -0.010 | -0.055 | Weak negative effect |
| SL | -0.071 | -0.106 | Body size negatively affects both |

**Jaw length** shows strong positive selection on growth rate (S = 0.325), while body size (SL) shows consistent negative selection on both survival and growth.

**Body size (SL)** shows consistent negative selection on both survival and growth, suggesting that smaller individuals have a fitness advantage in the high-density enclosure environment.

**Eye diameter** shows opposing patterns: weakly positive on survival but negative on growth, potentially reflecting a trade-off between visual acuity and metabolic costs.

---


#### 3.2 Selection Coefficients (Binary)

| Term        | Type           | Beta_Coefficient | Standard_Error | P_Value | Variance |
|-------------|----------------|------------------|----------------|---------|----------|
| jaw         | Linear         | 0.2790           | 0.1709         | 0.1025  | 0.0292   |
| eye         | Linear         | 0.1927           | 0.1141         | 0.0912  | 0.0130   |
| body        | Linear         | 0.0154           | 0.1339         | 0.9087  | 0.0179   |
| nasal       | Linear         | -0.3393          | 0.0990         | 0.0006  | 0.0098   |
| mouth       | Linear         | 0.0363           | 0.0992         | 0.7144  | 0.0098   |
| SL          | Linear         | -0.3674          | 0.1333         | 0.0059  | 0.0178   |
| jaw²        | Quadratic      | -0.6784          | 0.5076         | 0.1814  | 0.2577   |
| eye²        | Quadratic      | -0.3411          | 0.2329         | 0.1431  | 0.0542   |
| body²       | Quadratic      | -0.1813          | 0.3091         | 0.5575  | 0.0956   |
| nasal²      | Quadratic      | 0.0637           | 0.1693         | 0.7068  | 0.0287   |
| mouth²      | Quadratic      | -0.2349          | 0.1704         | 0.1681  | 0.0290   |
| SL²         | Quadratic      | 0.1963           | 0.3182         | 0.5373  | 0.1012   |
| jaw × eye   | Correlational  | -0.1539          | 0.2824         | 0.5857  | 0.0797   |
| jaw × body  | Correlational  | -0.3722          | 0.3265         | 0.2543  | 0.1066   |
| jaw × nasal | Correlational  | -0.2299          | 0.2304         | 0.3184  | 0.0531   |
| jaw × mouth | Correlational  | 0.3915           | 0.2459         | 0.1114  | 0.0605   |
| jaw × SL    | Correlational  | -0.3998          | 0.3060         | 0.1913  | 0.0936   |
| eye × body  | Correlational  | -0.3862          | 0.2094         | 0.0651  | 0.0438   |
| eye × nasal | Correlational  | 0.1372           | 0.1535         | 0.3714  | 0.0236   |
| eye × mouth | Correlational  | -0.0126          | 0.1573         | 0.9359  | 0.0248   |
| eye × SL    | Correlational  | -0.2554          | 0.1970         | 0.1949  | 0.0388   |
| body × nasal| Correlational  | -0.2843          | 0.1696         | 0.0937  | 0.0288   |
| body × mouth| Correlational  | 0.3237           | 0.1749         | 0.0642  | 0.0306   |
| body × SL   | Correlational  | -0.2508          | 0.2278         | 0.2708  | 0.0519   |
| nasal × mouth| Correlational | 0.2327           | 0.1380         | 0.0918  | 0.0190   |
| nasal × SL  | Correlational  | 0.2068           | 0.1675         | 0.2170  | 0.0281   |
| mouth × SL  | Correlational  | 0.1534           | 0.1675         | 0.3597  | 0.0280   |

#### 3.3 Selection Coefficients (Continuouse)

| Term        | Type           | Beta_Coefficient | Standard_Error | P_Value | Variance |
|-------------|----------------|------------------|----------------|---------|----------|
| jaw         | Linear         | 0.4915           | 0.1767         | 0.0055  | 0.0312   |
| eye         | Linear         | -0.0303          | 0.1218         | 0.8039  | 0.0148   |
| body        | Linear         | 0.1230           | 0.1484         | 0.4072  | 0.0220   |
| nasal       | Linear         | -0.2154          | 0.1106         | 0.0518  | 0.0122   |
| mouth       | Linear         | -0.0021          | 0.1094         | 0.9844  | 0.0120   |
| SL          | Linear         | 0.1306           | 0.1488         | 0.3804  | 0.0221   |
| jaw²        | Quadratic      | 0.4971           | 0.4392         | 0.2581  | 0.1929   |
| eye²        | Quadratic      | -0.1070          | 0.2075         | 0.6064  | 0.0431   |
| body²       | Quadratic      | -0.2846          | 0.3218         | 0.3766  | 0.1035   |
| nasal²      | Quadratic      | 0.1907           | 0.1704         | 0.2633  | 0.0290   |
| mouth²      | Quadratic      | 0.4548           | 0.1651         | 0.0060  | 0.0273   |
| SL²         | Quadratic      | 0.1278           | 0.3397         | 0.7070  | 0.1154   |
| jaw × eye   | Correlational  | 0.1860           | 0.2489         | 0.4550  | 0.0619   |
| jaw × body  | Correlational  | -0.0586          | 0.3074         | 0.8487  | 0.0945   |
| jaw × nasal | Correlational  | -0.0171          | 0.2144         | 0.9365  | 0.0460   |
| jaw × mouth | Correlational  | 0.0656           | 0.2223         | 0.7679  | 0.0494   |
| jaw × SL    | Correlational  | 0.2401           | 0.2919         | 0.4110  | 0.0852   |
| eye × body  | Correlational  | -0.2547          | 0.1951         | 0.1922  | 0.0381   |
| eye × nasal | Correlational  | 0.1695           | 0.1511         | 0.2624  | 0.0228   |
| eye × mouth | Correlational  | -0.0737          | 0.1464         | 0.6149  | 0.0214   |
| eye × SL    | Correlational  | 0.2071           | 0.1880         | 0.2710  | 0.0353   |
| body × nasal| Correlational  | -0.1829          | 0.1681         | 0.2771  | 0.0283   |
| body × mouth| Correlational  | 0.0285           | 0.1642         | 0.8621  | 0.0270   |
| body × SL   | Correlational  | 0.2360           | 0.2401         | 0.3258  | 0.0576   |
| nasal × mouth| Correlational | -0.1835          | 0.1360         | 0.1775  | 0.0185   |
| nasal × SL  | Correlational  | 0.0475           | 0.1668         | 0.7761  | 0.0278   |
| mouth × SL  | Correlational  | 0.3349           | 0.1664         | 0.0445  | 0.0277   |


#### 3.4 Compare with paper
| Trait | Martin (2016) - edf | Results (Binary) | Results (Continuous) | Consistency |
|-------|---------------------|----------------------|--------------------------|-------------|
| jaw | 2.1 (moderate) | Non-significant | Significant positive linear | Partial |
| eye | 1.8 (low) | Non-significant | Non-significant | Consistent |
| body | 8.8 (high) | Non-significant | Non-significant | Discrepant |
| nasal | 1.8 (low) | Significant negative linear | Marginally significant | Consistent |
| mouth | 2.0 (moderate) | Non-significant | Significant nonlinear (γ > 0) | Consistent |
| SL | 1.5 (low) | Significant negative linear | Non-significant | Partial |

*edf = effective degrees of freedom; higher values indicate greater nonlinearity*

Our results are generally consistent with Martin (2016), who described the adaptive fitness landscape in pupfishes. Both studies show positive selection on jaw length and negative selection on nasal protrusion, reflecting differences in feeding strategies. We also found nonlinear selection on mouth width, with a significant positive quadratic effect (γ = 0.455, P = 0.006), suggesting that extreme values may be favored.

However, there are some differences. Martin (2016) found strong nonlinear selection on body depth, but we did not detect significant effects for this trait. This may be due to differences in data, environmental conditions, or the inclusion of growth as an additional fitness component in our analysis.

Importantly, our results show that growth reveals stronger selection patterns than survival, including strong positive selection on jaw length and evidence of trait interactions. Overall, our findings support the main conclusions of Martin (2016) but also show that selection can vary depending on the fitness measure used.

---

#### 3.5 Univariate Fitness 

<table>
<tr>
<td><b>(A) Binary_body</b><br><img src="../results/pupfish_results/figures/univariate_binary_body.png" width="100%"></td>
<td><b>(B) Binary_eye</b><br><img src="../results/pupfish_results/figures/univariate_binary_eye.png" width="100%"></td>
<td><b>(C) Binary_jaw</b><br><img src="../results/pupfish_results/figures/univariate_binary_jaw.png" width="100%"></td>
</tr>

<tr>
<td><b>(D) Binary_mouth</b><br><img src="../results/pupfish_results/figures/univariate_binary_mouth.png" width="100%"></td>
<td><b>(E) Binary_nasal</b><br><img src="../results/pupfish_results/figures/univariate_binary_nasal.png" width="100%"></td>
<td><b>(F) Binary_SL</b><br><img src="../results/pupfish_results/figures/univariate_binary_SL.png" width="100%"></td>
</tr>

<tr>
<td><b>(G) Contiuous_body</b><br><img src="../results/pupfish_results/figures/univariate_continuous_body.png" width="100%"></td>
<td><b>(H) Contiuous_eye</b><br><img src="../results/pupfish_results/figures/univariate_continuous_eye.png" width="100%"></td>
<td><b>(I) Contiuous_jaw</b><br><img src="../results/pupfish_results/figures/univariate_continuous_jaw.png" width="100%"></td>
</tr>

<tr>
<td><b>(G) Contiuous_mouth</b><br><img src="../results/pupfish_results/figures/univariate_continuous_mouth.png" width="100%"></td>
<td><b>(H) Contiuous_nasal</b><br><img src="../results/pupfish_results/figures/univariate_continuous_nasal.png" width="100%"></td>
<td><b>(I) Contiuous_SL</b><br><img src="../results/pupfish_results/figures/univariate_continuous_SL.png" width="100%"></td>
</tr>
</table>

---

#### 3.6 Correlations Fitness 

<table>

<tr>
<td><b>(A) Binary_body_mouth</b><br><img src="../results/pupfish_results/figures/cfs_binary_body_mouth_enhanced.png" width="100%"></td>
<td><b>(B) Binary_body_nasal</b><br><img src="../results/pupfish_results/figures/cfs_binary_body_nasal_enhanced.png" width="100%"></td>
<td><b>(C) Binary_body_SL</b><br><img src="../results/pupfish_results/figures/cfs_binary_body_SL_enhanced.png" width="100%"></td>
<td><b>(D) Binary_eye_body</b><br><img src="../results/pupfish_results/figures/cfs_binary_eye_body_enhanced.png" width="100%"></td>
<td><b>(E) Binary_eye_mouth</b><br><img src="../results/pupfish_results/figures/cfs_binary_eye_mouth_enhanced.png" width="100%"></td>
</tr>

<tr>
<td><b>(F) Binary_eye_nasal</b><br><img src="../results/pupfish_results/figures/cfs_binary_eye_nasal_enhanced.png" width="100%"></td>
<td><b>(G) Binary_eye_SL</b><br><img src="../results/pupfish_results/figures/cfs_binary_eye_SL_enhanced.png" width="100%"></td>
<td><b>(H) Binary_jaw_body</b><br><img src="../results/pupfish_results/figures/cfs_binary_jaw_body_enhanced.png" width="100%"></td>
<td><b>(I) Binary_jaw_eye</b><br><img src="../results/pupfish_results/figures/cfs_binary_jaw_eye_enhanced.png" width="100%"></td>
<td><b>(J) Binary_jaw_mouth</b><br><img src="../results/pupfish_results/figures/cfs_binary_jaw_mouth_enhanced.png" width="100%"></td>
</tr>

<tr>
<td><b>(K) Binary_jaw_nasal</b><br><img src="../results/pupfish_results/figures/cfs_binary_jaw_nasal_enhanced.png" width="100%"></td>
<td><b>(L) Binary_jaw_SL</b><br><img src="../results/pupfish_results/figures/cfs_binary_jaw_SL_enhanced.png" width="100%"></td>
<td><b>(M) Binary_mouth_SL</b><br><img src="../results/pupfish_results/figures/cfs_binary_mouth_SL_enhanced.png" width="100%"></td>
<td><b>(N) Binary_nasal_mouth</b><br><img src="../results/pupfish_results/figures/cfs_binary_nasal_mouth_enhanced.png" width="100%"></td>
<td><b>(O) Binary_nasal_SL</b><br><img src="../results/pupfish_results/figures/cfs_binary_nasal_SL_enhanced.png" width="100%"></td>
</tr>

<tr>
<td><b>(P) Continuous_body_mouth</b><br><img src="../results/pupfish_results/figures/cfs_continuous_body_mouth_enhanced.png" width="100%"></td>
<td><b>(Q) Continuous_body_nasal</b><br><img src="../results/pupfish_results/figures/cfs_continuous_body_nasal_enhanced.png" width="100%"></td>
<td><b>(R) Continuous_body_SL</b><br><img src="../results/pupfish_results/figures/cfs_continuous_body_SL_enhanced.png" width="100%"></td>
<td><b>(S) Continuous_eye_body</b><br><img src="../results/pupfish_results/figures/cfs_continuous_eye_body_enhanced.png" width="100%"></td>
<td><b>(T) Continuous_eye_mouth</b><br><img src="../results/pupfish_results/figures/cfs_continuous_eye_mouth_enhanced.png" width="100%"></td>
</tr>

<tr>
<td><b>(U) Continuous_eye_nasal</b><br><img src="../results/pupfish_results/figures/cfs_continuous_eye_nasal_enhanced.png" width="100%"></td>
<td><b>(V) Continuous_eye_SL</b><br><img src="../results/pupfish_results/figures/cfs_continuous_eye_SL_enhanced.png" width="100%"></td>
<td><b>(W) Continuous_jaw_body</b><br><img src="../results/pupfish_results/figures/cfs_continuous_jaw_body_enhanced.png" width="100%"></td>
<td><b>(X) Continuous_jaw_eye</b><br><img src="../results/pupfish_results/figures/cfs_continuous_jaw_eye_enhanced.png" width="100%"></td>
<td><b>(Y) Continuous_jaw_mouth</b><br><img src="../results/pupfish_results/figures/cfs_continuous_jaw_mouth_enhanced.png" width="100%"></td>
</tr>

<tr>
<td><b>(Z) Continuous_jaw_nasal</b><br><img src="../results/pupfish_results/figures/cfs_continuous_jaw_nasal_enhanced.png" width="100%"></td>
<td><b>(AA) Continuous_jaw_SL</b><br><img src="../results/pupfish_results/figures/cfs_continuous_jaw_SL_enhanced.png" width="100%"></td>
<td><b>(AB) Continuous_mouth_SL</b><br><img src="../results/pupfish_results/figures/cfs_continuous_mouth_SL_enhanced.png" width="100%"></td>
<td><b>(AC) Continuous_nasal_mouth</b><br><img src="../results/pupfish_results/figures/cfs_continuous_nasal_mouth_enhanced.png" width="100%"></td>
<td><b>(AD) Continuous_nasal_SL</b><br><img src="../results/pupfish_results/figures/cfs_continuous_nasal_SL_enhanced.png" width="100%"></td>
</tr>

</table>

##### Binary Fitness (Survival)

| Trait Pair | Key Finding |
|------------|-------------|
| nasal × SL | Both traits exhibit negative directional selection |
| jaw × SL | Contrasting directional selection: jaw (+) vs SL (-) |
| body × mouth | Weak or marginal nonlinear effects (p = 0.06) |

##### Continuous Fitness (Growth)

| Trait Pair | Key Finding |
|------------|-------------|
| jaw × mouth | Strong positive selection on jaw length and nonlinear curvature in mouth width |
| mouth × SL | Significant correlational selection (γ = 0.335, P = 0.045) |
| jaw × SL | Strong directional selection on jaw with a weaker trend for SL |

---

#### 3.5 Disruptive (Nonlinear) Selection 

| Trait | Fitness     | Beta (Linear) | Gamma (Quadratic) | P_Linear | P_Quadratic | Interpretation |
|-------|-------------|---------------|-------------------|----------|-------------|----------------|
| jaw   | Binary      | 0.4206        | 0.1001            | <0.001   | 0.5182      | Strong directional selection |
| jaw   | Continuous  | 0.3433        | 0.0958            | <0.001   | 0.5536      | Strong directional selection |
| eye   | Binary      | 0.0723        | -0.4118           | 0.4731   | 0.0109      | Significant stabilizing selection |
| eye   | Continuous  | -0.1796       | -0.2865           | 0.0696   | 0.0388      | Stabilizing selection |
| body  | Binary      | -0.3307       | 0.0909            | <0.001   | 0.4589      | Directional (negative) |
| body  | Continuous  | -0.0906       | 0.1084            | 0.3671   | 0.4407      | No strong selection |
| nasal | Binary      | -0.1610       | 0.2920            | 0.0732   | 0.0028      | Significant nonlinear (potential disruptive) |
| nasal | Continuous  | -0.1736       | 0.1269            | 0.0951   | 0.2686      | Weak nonlinear |
| mouth | Binary      | -0.0414       | 0.1472            | 0.6247   | 0.1693      | No strong selection |
| mouth | Continuous  | -0.0180       | 0.3428            | 0.8560   | 0.0096      | Significant nonlinear (potential disruptive) |
| SL    | Binary      | -0.4002       | 0.1639            | <0.001   | 0.1946      | Directional (negative) |
| SL    | Continuous  | -0.0757       | 0.1368            | 0.4637   | 0.3430      | No strong selection |


Quadratic selection revealed both stabilizing and nonlinear patterns. **Eye diameter** showed consistent stabilizing selection (γ = -0.41 to -0.29, P < 0.05), while **nasal protrusion** (binary: γ = 0.29, P = 0.003) and **mouth width** (continuous: γ = 0.34, P = 0.010) showed significant positive quadratic effects, indicating disruptive selection.

Directional selection was stronger overall, with **positive selection on jaw length** (β ≈ 0.35–0.42, P < 0.001) and **negative selection on body size (SL)** (β = -0.40, P < 0.001).

Overall, selection is predominantly directional, with nonlinear effects confined to specific traits (eye, nasal, mouth).


#### 3.6 Adaptive Landscape 

<table>

<tr>
<td><b>(A) Binary_body_mouth (2D)</b><br>
<img src="../results/pupfish_results/figures/adaptive_landscape_binary_body_mouth_2d.png" width="100%"></td>

<td><b>(B) Binary_body_mouth (3D)</b><br>
<img src="../results/pupfish_results/figures/adaptive_landscape_binary_body_mouth_3d.png" width="100%"></td>
</tr>

<tr>
<td><b>(C) Binary_jaw_SL (2D)</b><br>
<img src="../results/pupfish_results/figures/adaptive_landscape_binary_jaw_SL_2d.png" width="100%"></td>

<td><b>(D) Binary_jaw_SL (3D)</b><br>
<img src="../results/pupfish_results/figures/adaptive_landscape_binary_jaw_SL_3d.png" width="100%"></td>
</tr>

<tr>
<td><b>(E) Binary_nasal_SL (2D)</b><br>
<img src="../results/pupfish_results/figures/adaptive_landscape_binary_nasal_SL_2d.png" width="100%"></td>

<td><b>(F) Binary_nasal_SL (3D)</b><br>
<img src="../results/pupfish_results/figures/adaptive_landscape_binary_nasal_SL_3d.png" width="100%"></td>
</tr>

<tr>
<td><b>(G) Continuous_jaw_mouth (2D)</b><br>
<img src="../results/pupfish_results/figures/adaptive_landscape_continuous_jaw_mouth_2d.png" width="100%"></td>

<td><b>(H) Continuous_jaw_mouth (3D)</b><br>
<img src="../results/pupfish_results/figures/adaptive_landscape_continuous_jaw_mouth_3d.png" width="100%"></td>
</tr>

<tr>
<td><b>(I) Continuous_jaw_SL (2D)</b><br>
<img src="../results/pupfish_results/figures/adaptive_landscape_continuous_jaw_SL_2d.png" width="100%"></td>

<td><b>(J) Continuous_jaw_SL (3D)</b><br>
<img src="../results/pupfish_results/figures/adaptive_landscape_continuous_jaw_SL_3d.png" width="100%"></td>
</tr>

<tr>
<td><b>(K) Continuous_mouth_SL (2D)</b><br>
<img src="../results/pupfish_results/figures/adaptive_landscape_continuous_mouth_SL_2d.png" width="100%"></td>

<td><b>(L) Continuous_mouth_SL (3D)</b><br>
<img src="../results/pupfish_results/figures/adaptive_landscape_continuous_mouth_SL_3d.png" width="100%"></td>
</tr>

</table>

---

#### 3.7 Correlated Fitness Surface (Individual) vs. Adaptive Landscape (Population)

| Fitness     | Trait1 | Trait2 | Correlation      | N    |
|-------------|--------|--------|------------------|------|
| Binary      | jaw    | SL     | 0.9405           | 2500 |
| Binary      | nasal  | SL     | 0.9935           | 2500 |
| Continuous  | jaw    | mouth  | 0.9565           | 2500 |
| Continuous  | mouth  | SL     | 0.8806           | 2500 |


<table>
<tr>
<td width="40%" align="center">
<b>(A) Binary_jaw_SL (Overlay)</b><br>
<img src="../results/pupfish_results/figures/comparison_binary_jaw_SL_overlay.png" width="95%">
</td>

<td width="60%" align="center">
<b>(B) Binary_jaw_SL (Side-by-Side)</b><br>
<img src="../results/pupfish_results/figures/comparison_binary_jaw_SL_side_by_side.png" width="95%">
</td>
</tr>
</table>

<table>
<tr>
<td width="40%" align="center">
<b>(C) Binary_nasal_SL (Overlay)</b><br>
<img src="../results/pupfish_results/figures/comparison_binary_nasal_SL_overlay.png" width="95%">
</td>

<td width="60%" align="center">
<b>(D) Binary_nasal_SL (Side-by-Side)</b><br>
<img src="../results/pupfish_results/figures/comparison_binary_nasal_SL_side_by_side.png" width="95%">
</td>
</tr>
</table>

<table>
<tr>
<td width="40%" align="center">
<b>(E) Continuous_jaw_mouth (Overlay)</b><br>
<img src="../results/pupfish_results/figures/comparison_continuous_jaw_mouth_overlay.png" width="95%">
</td>

<td width="60%" align="center">
<b>(F) Continuous_jaw_mouth (Side-by-Side)</b><br>
<img src="../results/pupfish_results/figures/comparison_continuous_jaw_mouth_side_by_side.png" width="95%">
</td>
</tr>
</table>

<table>
<tr>
<td width="40%" align="center">
<b>(G) Continuous_mouth_SL (Overlay)</b><br>
<img src="../results/pupfish_results/figures/comparison_continuous_mouth_SL_overlay.png" width="95%">
</td>

<td width="60%" align="center">
<b>(H) Continuous_mouth_SL (Side-by-Side)</b><br>
<img src="../results/pupfish_results/figures/comparison_continuous_mouth_SL_side_by_side.png" width="95%">
</td>
</tr>
</table>


### 4. Comparison with Martin (2016)

Our results are broadly consistent with the key findings of Martin (2016), who first characterized the adaptive landscape in this pupfish radiation. Both studies identify:

- **Positive selection on jaw length** (β = 0.49, P = 0.006)
- **Negative selection on nasal protrusion** (β = -0.34, P = 0.001)
- **Negative selection on body size (SL)** (β = -0.37, P = 0.006)
- **Stabilizing selection on eye diameter** (γ = -0.41, P = 0.011)
- **Nonlinear selection on mouth width** (γ = 0.34, P = 0.010)

However, several differences emerge:

- **Body depth nonlinearity**: Martin (2016) reported strong nonlinear selection on body depth (edf = 8.8), whereas we detected no significant quadratic effects (γ = 0.09, P = 0.46).

- **Correlational selection**: We detected significant correlational selection between mouth width and standard length (γ = 0.335, P = 0.045), a pattern not emphasized in the original study.

- **Growth as a fitness component**: By including growth alongside survival, we reveal stronger and more complex selection patterns, indicating that different fitness components capture distinct aspects of selection. Notably, directional selection on jaw length was stronger for growth (β = 0.49) than for survival (β = 0.28).

Overall, our results support the main conclusions of Martin (2016), while suggesting that selection is primarily directional with additional nonlinear effects confined to specific traits (eye, nasal, mouth). The inclusion of growth as a fitness component reveals novel correlational selection and highlights the multidimensional nature of the adaptive landscape in this system.
