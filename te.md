# Methodological Details

Authors: Maximilian Sprengholz and Maik Hamjediers

**Table of contents**
1. [Ñopo's (2008) matching decomposition](#ñopos-2008-matching-decomposition)
   1. [On the Choice of the Matching Procedure](#on-the-choice-of-the-matching-procedure)
   2. [Relationship to Regression-Based Decompositions](#relationship-to-regression-based-decompositions)
   3. [Relationship to Estimation of Treatment Effects via Matching](#relationship-to-estimation-of-treatment-effects-via-matching)
2. [Post-Estimation Statistics](#post-estimation-statistics)
   1. [Descriptive statistics by matching status (`summarize`)](#descriptive-statistics-by-matching-status-summarize)
   2. [Decomposition-components across the distribution of $Y$ (`gapoverdist`)](#decomposition-components-across-the-distribution-of-y-gapoverdist)
   3. [Contribution of characteristics to the components $D_A$ and $D_B$ (`dadb`)](#contribution-of-characteristics-to-the-components-d_a-and-d_b-dadb)
3. [References](#references)


## Ñopo's (2008) matching decomposition

We are interested in the decomposition of the raw gap $D$ in outcome $Y$ between group A and group B

$$
D = \overline{Y}_B - \overline{Y}_A \quad ,
$$

into an "explained" component that is based on differences between groups in characteristics that predict $Y$ and a remaining "unexplained" component. However, we do not want to compare apples and oranges, so that the decomposition also needs to calculate which part of $D$ is due to units in both groups which have combinations of predictive characteristics out of common support (characteristics that do not occur in both groups).

Ñopo's matching achieves these goals as follows. In a one-to-many exact matching, each individual from group $B$ is matched to all individuals from group $A$ with the same combination of characteristics (each unique combination of characteristics represents one stratum). The matching flags for all observations of groups $A$ and $B$ if their characteristics are in common support (matched $m$) or out of common support (unmatched $u$). Over common support, the ratio of $B$ to $A$ units in each stratum can be used to create a reweighted group $A^B$ which has the exact same distribution across all strata as group $B$. The outcome of this counterfactual group $A^B$ can be interpreted in two ways: 
1. as the average outcome of group $A$ if it had the same characteristics as group $B$ (rendering group $B$ as the reference group in terms of characteristics, which we denote throughout the documentation and output as `xref');
2. as the average outcome of group $B$ if it had the same returns to characteristics as group $A$ (rendering group $A$ as the reference group in terms of returns, which we denote throughout the documentation and output as `bref').

The overall gap can then be additively decomposed into four parts:

```math
\begin{equation} \tag{1}
\begin{array}{rcccccc}
D &=& D_0  &+& D_X &+& D_A + D_B \\
&=& \overbrace{\overline{Y}_{B,m} - \overline{Y}_{A^B,m}} &+& \overbrace{\overline{Y}_{A^B,m} - \overline{Y}_{A,m}} &+& \underbrace{D_A + D_B}_{\mathclap{\substack{\text{out of} \\ \text{support}}}}
\end{array}
\end{equation}
```

$D_X$ is the average gap between the matched units of re-weighted group $A^B$ and the matched units of group $A$, which is explained by the fact that groups $A$ and $B$ are differently distributed across matched strata (some sets of characteristics being more likely in one group than the other). 
By contrast, $D_0$ is the average gap between matched units of group $B$ and the re-weighted group $A^B$. Since $B$ and $A^B$ are equally distributed across matched strata, $D_0$ captures how much of the raw gap remains unexplained by differences in the considered characteristics. When compositional differences between groups limit common support, the effect that unmatched individuals in both groups have on $D$ is captured by the components $D_A$ and $D_B$:

```math
\begin{equation} \tag{2}
	D_A =  \underbrace{(\overline{Y}_{A,m} - \overline{Y}_{A,u})}_{\mathclap{\substack{\text{gap between matched} \\ \text{and unmatched }A}}} \cdot \underbrace{(N_{A,u}/N_A)}_{\mathclap{\substack{\text{share of} \\ \text{unmatched }A}}}  \quad\quad
	D_B =  \underbrace{(\overline{Y}_{B,u} - \overline{Y}_{B,m})}_{\mathclap{\substack{\text{gap between unmatched} \\ \text{and matched }B}}} \cdot \underbrace{(N_{B,u}/N_B)}_{\mathclap{\substack{\text{share of} \\ \text{unmatched }B}}} 
\end{equation}
```

$D_A$ is the gap between the averages of the outcome $\overline{Y}$ for the unmatched $u$ and matched $m$ units within group $A$, weighted by the frequency of unmatched $A$ units ($N_{A,u}$) in relation to the overall size of group $A$ ($N_{A}$), so that $D_A$ approaches zero with fewer observations out of common support. $D_A$ denotes how much of the raw gap is due to unmatched $A$ units having higher or lower values in the outcome than matched $A$ units, where $D_A < 0$ if the outcome is lower among the matched and $D_A > 0$ if the outcome is lower among the unmatched (reversed for $D_B$).

Naturally, one could also change the direction of the matching and create a counterfactual group of $B^A$.

### On the Choice of the Matching Procedure

As any kind of matching on characteristics, the matching-based decomposition faces a trade-off of, on the one hand, balancing the counterfactual distribution of $X_{A^B,m}$ and its target distribution of $X_{A,m}$, and on the other hand, maintaining sufficient common support (Iacus et al. 2012). The latter is also known as the curse of dimensionality, which describes that increasing the number of characteristics and/or their distinct values, increases the number of strata for the matching, which in return decreases the likelihood of observing exact matches of both groups in finite samples. To address this, several alternate matching procedures have been proposed, including coarsened exact matching (Iacus et al. 2012) as well as matching on propensity scores or multivariate (Mahalanobis) distances via nearest neighbor or kernel-bandwidth approaches. Naturally, the results of the matching-based decompositions hinge on the specifics of the underlying matching procedure (e.g., the extent of coarsening of continuous variables before exact matching, or the bandwidth selection for determining matches in terms of the propensity score or multivariate distance).

While in contrast to (coarsened) exact matching, these alternative matching procedures tackle the curse of dimensionality to some extent, they also increase the likelihood that substantial and theoretically meaningful lacks of common support being not attributed to $D_A$ and $D_B$ (e.g., due to a glass-ceiling one group is never observed in some strata or institutional barriers prompting one group to work in low wage occupations despite high educational attainment and this overqualification does not occur among another group). Thus, the trade-off between balancing and the curse of dimensionality cannot be avoided, and careful exploration of instances of unmatched units as well as the discrepancy between the distributions of $X_{A^B,m}$ and $X_{A,m}$ (see [`nopo summarize`](te.md#descriptive-statistics-by-machting-status-summarize)) is needed. 

Note that only (coarsened) exact matching ensures that the interpretations of the unexplained $D_0$ and explained $D_X$ components directly refer to the characteristics $X$ (e.g., $D_0$ being the remaining gap if both groups had the same characteristics as group $B$). This interpretation changes slightly in the case of propensity-score or multivariate distance matching, as $D_0$ then refers to both groups having an equal likelihood to be the group specified in `xref` based on the characteristics $X$. 

### Relationship to Regression-Based Decompositions

$D_0$ and $D_X$ from the matching-based decomposition are analogous to the components of a two-fold regression-based composition as proposed by Oaxaca (1973), Blinder (1973), and Kitagawa (1955), but they only pertain to matched units. 

These decompositions build on (1) group-specific vectors of the mean values $\overline{X}_A$ and $\overline{X}_B$ for the specified predictors; and (2) group-specific vectors of coefficients $\hat{\beta}_A$ and $\hat{\beta}_B$ obtained from a regression of the outcome on a set of characteristics for each group. The mean-differences in characteristics represent compositional differences that make up the explained component $D_X$, whereas differences in the associated regression coefficients represent differences in returns that make up the unexplained component $D_0$:

```math
\begin{equation}
\begin{array}{rcccc}
    D &=&  \overline{X}'_B \underbrace{\big(\hat{\beta}_B - \hat{\beta}_A\big)}_{\mathclap{\substack{\text{difference} \\ \text{in returns}}}} &+& \underbrace{\big(\overline{X}_B - \overline{X}_A\big)'}_{\mathclap{\substack{\text{compositional} \\ \text{difference}}}}\hat{\beta}_A\\
    &=& \qquad D_0  &+& \mkern-18mu D_X
\end{array}
\end{equation}
```

Therein, the mean-differences in characteristics are multiplied with a coefficient vector that serves as the reference vector for the invoked counterfactual scenario.  As long as it is taken from one of the groups underlying the decomposition (e.g., here as $\hat{\beta}_A$ from the group $A$), it provides an analogous interpretation for the unexplained component $D_0$ as in the matching-based decomposition (with `bref` denoting the group whose coefficient-vector is applied). Alternatively, the reference vector could be defined in any other way (e.g., from a pooled regression of both groups), which does not ensure a direct analogy between the matching- and regression-based decomposition components anymore.

The matching-based decomposition has some advantages and disadvantages compared to a regression-based decomposition. On the one hand, matching allows for non-parametric estimation of all decomposition components that avoid biases due to mis-specifications of the functional form between the characteristics $X$ (including their interactions) and the outcome. It also ensures that the components $D_0$ and $D_X$ only apply to matched units, which avoids model-based extrapolations that could cover meaningful lacks of common support. On the other hand, matching suffers from the curse of dimensionality and runs at the risk of attributing too much of the gap to $D_A$ and/or $D_B$ if lacking common support is only due to limited sample size. It also does not allow to disentangle how each characteristic contributes to $D_X$ and $D_0$. Overall, these arguments are very similar to the literature on regression- vs. matching-based adjustment when estimating (local) treatment effects. Based on a simulation study that highlights and discusses these issues (Hamjediers & Sprengholz, 2023), we would recommend assessing the severity of lacking common support and functional form assumptions by first using a matching-based decomposition (and potentially restricting the analytic sample to observations with common support).

### Relationship to Estimation of Treatment Effects via Matching 

Depending on the direction of the matching underlying the decomposition, $D_0$ corresponds to either $ATT$ (Average Treatment effect on the Treated) or $ATC$ (Average Treatment effect on the Unreated/Controls) from the treatment effects literature. 

Let $PO$ denote the (potential) outcome $Y$, the observed treatment status is denoted by $T$, with $T=1$ being treated and $T=0$ being untreated, and a potential treatment status that may or may not be observed is denoted by $t$. Using the potential outcome framework (Rubin, 2008), the $ATT$ is then defined as $ATT = PO_{t=1}^{T=1} - PO_{t=0}^{T=1}$, thus, as the average difference between the observed outcome for the treated as if they were treated ($PO_{t=1}^{T=1}$) and the unobservable, counterfactual outcome for the treated as if they were untreated ($PO_{t=0}^{T=1}$). 

To estimate the $ATT$, we can use the observed outcome for all treated units to estimate $PO_{t=1}^{T=1}$, as the treatment status is observed coherently. By contrast, $PO_{t=0}^{T=1}$ can be estimated as the outcome of the untreated, control group $T=0$ as if it would have the characteristics of the treated group under the assumption of conditional independence. This counterfactual potential outcome is given by matching the treated and untreated groups and producing weights that allow to weigh untreated group to match the strata of the treated group (producing a counterfactual group of $t=0^{t=1}$:

```math
\begin{equation}
\begin{array}{r c c c c }
	ATT &=& PO_{t=1}^{T=1} &-& PO_{t=0}^{T=1} \\
		&=& \overline{Y}_{t=1,m} &-& \overline{Y}_{t=0^{t=1},m} \\
\end{array}
\end{equation}
```

If group $B$ from our decomposition notation is assigned as the treated group $t=1$ (and group $A$ is untreated $t=0$), we see that $ATT$ equals the unexplained component of $D_0$.

Alternatively, we could also switch the matching direction to estimate the average treatment effect on the controls $ATC$, which gives us:

```math
\begin{equation}
\begin{array}{r c c c c c c }
	ATC &=& PO_{t=1}^{T=0} &-& PO_{t=0}^{T=0} && \\
		&=& \overline{Y}_{t=1^{t=0},m} &-& \overline{Y}_{t=0,m} && \\
		&=& \overline{Y}_{B^A,m} &-& \overline{Y}_{A,m} &=&  D_0 \\
\end{array}
\end{equation}
```

Thus, the $ATC$ equals again the unexplained decomposition-component $D_0$, however, with the reversed matching direction. To indicate which group's characteristics are used to estimate the counterfactual potential outcome, refer to which group is assigned the `xref`-label. 

Note that the estimation of any treatment effects ($ATT$ as well as $ATC$) via matching hinges on the assumption of conditional independence. This is also the reason why matching-based treatment effect estimation only provides estimates of $ATT$ or $ATC$ (or even $ATE$) and omits other components ($D_X$, $D_A$, or $D_B$), which are not of direct interest. While the decomposition framework also operates the same counterfactual as the treatment effect estimation (e.g., what would group's $A$ outcome be, if it had the characteristics of group $B$), it does not necessarily aim at a causal interpretation. If the decomposition components are of descriptive interest, the matching-based decomposition does not invoke the assumption of conditional independence. 


## Post-Estimation Statistics

### Descriptive statistics by matching status (`summarize`)

It provides descriptive statistics (e.g., means and standard deviations) for each group and by matching status. Additionally, the descriptive statistics for the invoked counterfactual group of $A^B$ is reported, which allows, on the one hand, to assess the balancing between the counterfactual weighted group $A^B$ and the target distribution of the characteristics-wise reference group $B$ (`xref`), and on the other hand, the extent of differences in the characteristics between $A^B$ and $A$, which underlay the explained component $D_X$. 

### Decomposition-components across the distribution of $Y$ (`gapoverdist`)

Plots decomposition-components over the distribution of $Y$.  

Therefore, we compare the means of $Y$ at each quantile for the respective groups (e.g. the matched and unmatched of group $B$ for $D_B$). This means for each quantile that the single component values *do not* add up to $D$. But the quantile values of each component sum to the overall decomposition component values.


### Contribution of characteristics to the components $D_A$ and $D_B$ (`dadb`)

Creates a plot showing how the different levels of one single characteristic $x$ contribute to the components $D_A$ and $D_B$. This contribution emerges because these levels are associated with either many unmatched units and/or large differences in the outcome $Y$ by matching status within groups $A$ and $B$.

To acknowledge the latter, the difference in the means between the matched and unmatched units of each group is calculated separetely for each level and depicted. However, the components $D_A$ and $D_B$ not only depend on the difference between matched and unmatched units, but also on the share of unmatched units (see Equation 2). Accordingly, the mean differences for each level are additionally weighted by the share of unmatched units for each level to obtain the level's contribution. In other words, we apply the equations that define $D_A$ and $D_B$ separetely for each level of the respective characteristic. 

Note that this contribution is *not the same* as a detailed decomposition in regression-based approaches (which is generally not possible with matching). The contribution to $D_A$ and $D_B$ pertains only to the comparison between matched and unmatched units among group $A$ and $B$ and is interdepent with the matching across all other characteristics of the matching set.


## References

Blinder A.S. (1973). Wage Discrimination: Reduced Form and Structural Estimates. Journal of Human Resources, 8(4), 436–55. https://doi.org/10.2307/144855

Hamjediers, M. & Sprengholz, M. (2023). Comparing the Incomparable? Issues of Lacking  Common Support, Functional Form Mis-Specification, and Insufficient Sample Size in  Decompositions. Sociological Methodology, 53(2), 344-365. https://doi.org/10.1177/00811750231169729

Iacus, S.M., King, G. & Porro, G. (2012). Causal Inference without Balance Checking: Coarsened Exact Matching. Political Analysis 20, 1–24. https://doi.org/10.1093/pan/mpr013

Kitagawa E.M. (1955). Components of a Difference between Two Rates. Journal of the American Statistical Association, 50(272), 1168–94. https://doi.org/10.1080/01621459.1955.10501299

Ñopo, H. (2008). Matching as a Tool to Decompose Wage Gaps. The Review of Economics and Statistics, 90(2), 290–299. https://doi.org/10/b6tqwq

Oaxaca R. (1973). Male-Female Wage Differentials in Urban Labor Markets. International Economic Review, 14(3), 693–709. https://doi.org/10.2307/2525981
