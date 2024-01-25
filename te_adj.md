# Methodological Details

Authors: Maximilian Sprengholz and Maik Hamjediers

## Ñopo's (2008) matching decomposition:

We are interested in the decomposition of the raw gap $D$ in outcome $Y$ between group A and group B

$$
D = \overline{Y}_B - \overline{Y}_A \quad ,
$$

into an “explained” component that is based on differences between groups in characteristics that predict $Y$ and a remaining “unexplained” component. However, we do not want to compare apples and oranges, so that the decomposition also needs to calculate which part of $D$ is due to units in both groups which have combinations of predictive characteristics out of common support (characteristics which do not occur in both groups).

Ñopo's matching achieves these goals as follows. In a one-to-many exact matching, each individual from group $B$ is matched to all individuals from group $A$ with the same combination of characteristics (each unique combination of characteristics represents one stratum). The matching flags for all observations of groups $A$ and $B$ if their characteristics are in common support (matched $m$) or out of common support (unmatched $u$). Over common support, the ratio of $B$ to $A$ units in each stratum can be used to create a reweighted group $A^B$ which has the exact same distribution across all strata as group $B$. The outcome of this counterfactual group $A^B$ can be interpreted in two ways: (1) as the average outcome of group $A$ if it had the same characteristics as group $B$ and (2) as the average outcome of group $B$ if it had the same returns to characteristics as group $A$.

The overall gap can then be additively decomposed into four parts:

```math
\begin{equation}
\begin{array}{rcccccc}
D &=& D_0  &+& D_X &+& D_A + D_B \\
&=& \overbrace{\overline{Y}_{B,m} - \overline{Y}_{A^B,m}} &+& \overbrace{\overline{Y}_{A^B,m} - \overline{Y}_{A,m}} &+& \underbrace{D_A + D_B}_{\mathclap{\substack{\text{out of} \\ \text{support}}}}
\end{array}
\end{equation}
```

$D_X$ is the average gap between the matched units of re-weighted group $A^B$ and the matched units of group $A$, which is explained by the fact that groups $A$ and $B$ are differently distributed across matched strata (some sets of characteristics being more likely in one group than the other). 
By contrast, $D_0$ is the average gap between matched units of group $B$ and the re-weighted group $A^B$. Since $B$ and $A^B$ are equally distributed across matched strata, $D_0$ captures how much of the raw gap remains unexplained by differences in the considered characteristics. $D_0$ and $D_X$ are analogous to the components of a twofold (Kitagawa-)Blinder-Oaxaca decomposition, but they only pertain to matched units. When compositional differences between groups limit common support, the effect that unmatched individuals in both groups have on $D$ is captured by the components $D_A$ and $D_B$:

```math
\begin{equation}
	D_A =  \underbrace{(\overline{Y}_{A,m} - \overline{Y}_{A,u})}_{\mathclap{\substack{\text{gap between matched} \\ \text{and unmatched }A}}} \cdot \underbrace{(N_{A,u}/N_A)}_{\mathclap{\substack{\text{share of} \\ \text{unmatched }A}}}  \quad\quad
	D_B =  \underbrace{(\overline{Y}_{B,u} - \overline{Y}_{B,m})}_{\mathclap{\substack{\text{gap between unmatched} \\ \text{and matched }B}}} \cdot \underbrace{(N_{B,u}/N_B)}_{\mathclap{\substack{\text{share of} \\ \text{unmatched }B}}} 
\end{equation}
```

$D_A$ is the gap between the averages of the outcome $\overline{Y}$ for the unmatched $u$ and matched $m$ units within group $A$, weighted by the frequency of unmatched $A$ units ($N_{A,u}$) in relation to the overall size of group $A$ ($N_{A}$), so that $D_A$ approaches zero with fewer observations out of common support. $D_A$ denotes how much of the raw gap is due to unmatched $A$ units having higher or lower values in the outcome than matched $A$ units, where $D_A < 0$ if the outcome is lower among the matched and $D_A > 0$ if the outcome is lower among the unmatched (reversed for $D_B$).

Note that the unexplained component $D_0$ has two possible interpretations. As it is not attributable to compositional differences, it is either

- the gap that remains if one group had the same characteristics as the other, reference group; or 
- the gap that remains if one group had the same returns as the other, reference group.

Throughout our description and in the output, we denote the reference group for characteristics-based interpretation as the `xref()`-group and for the return-based interpretation as the `bref()`-group.


### On the Choice of the Matching Procedure

Trade-off between reaching balance on predictors between $A^B,m$ and $A^,m$ vs. curse of dimensionality inducing lack of common support (Iacus et al. 2012). 

Accordingly, the results of the decomposition may hinge on the specifics of either matching procedure (e.g., the extent of coarsening of continuous variables before exact matching, or the bandwidth selection for determining matches in terms of the propensity-score or multivariate-distance).

Note that only exact matching ensures that the interpretation of the unexplained $D_0$ and explained $D_X$ components directly refer to the characteristics $X$ (e.g., $D_0$ being the remaining gap if both groups had the same characteristics as group $B$). This interpretation changes slightly in the case of multivariate-distance and propensity-score matching, as $D_0$ then refers to both groups having an equal likelihood to be the group specified in `xref()` based on the characteristics $X$.


### Relationship to Regression-Based Decompositions

Generally, very similar to a two-fold regression-based composition as proposed by Oaxaca (1973), Blinder (1973), and Kitagawa (19XX). These build on (1) group-specific vectors of the mean values $\overline{X}_A$ and $\overline{X}_B$ for the specified predictors; and (2) group-specific vectors of coefficients $\hat{\beta}_A$ and $\hat{\beta}_B$ obtained from a regression of the outcome on a set of predictors for each group. The mean-differences in predictors represent compositional differences that make up the explained component $D_X$ (e.g., wage differences due to differences in labor market experience), whereas differences in the associated regression coefficients represent differences in returns that make up the unexplained component $D_0$ (e.g., the same educational attainment might have different wage returns for each group):

```math
\begin{equation}
\begin{array}{rcccc}
    D &=&  \overline{X}'_B \underbrace{\big(\hat{\beta}_B - \hat{\beta}_A\big)}_{\mathclap{\substack{\text{difference} \\ \text{in returns}}}} &+& \underbrace{\big(\overline{X}_B - \overline{X}_A\big)'}_{\mathclap{\substack{\text{compositional} \\ \text{difference}}}}\hat{\beta}_A\\
    &=& \qquad D_0  &+& \mkern-18mu D_X
\end{array}
\end{equation}
```math

Therein, $\hat{\beta}_A$ serves as the coefficient-vector, which is denoted as `bref()` and allows for the return-based interpretation of the unexplained $D_0$ component in matching based decompositions detailled above.


In Hamjediers & Sprengholz (2023), we use simulations to highlight the advantages and disadvantages across both approaches, which we summarize here:

Advantages of matching-based decompositions:
+ Non-parametric estimation $\rightarrow$ no assumptions about functional form
+ $D_0$ \& $D_X$ apply only to matched units $\rightarrow$ no model-based extrapolation


Disdvantages of matching-based decompositions:
- Suffers from curse of dimensionality $\rightarrow$ risk of attributing too much to $D_A$ \& $D_B$
- Does not allow to disentangle explained component across predictors

Overall, these arguments are very similar to the literature on regression- vs. matching-based adjustment for confounders in estimating (local) treatment effects. We would recomend to assess severity of lacking common support and functional form-assumptions by first using a matching based decomposition (and potentially restricting the analytic sample to observations with common support).


### Relationship to Estimation of Treatment Effects via Matching 

Depending on the direction of the matching underlying the decomposition, $D_0$ corresponds to either $ATT$ (Average Treatment effect on the Treated) or $ATC$ (Average Treatment effect on the Unreated/Controls) from the treatment effects literature. 

Let $PO$ denote the (potential) outcome $Y$, $T$ denotes the observed treatment status of $T=1$ being treated and $T=0$ being untreated, and $t$ denotes a potential treatment status that may or may not be observed. Using the potential outcome framework (Rubin, 2008), the $ATT$ is then definied as $ATT = PO_{t=1}^{T=1} - PO_{t=0}^{T=1}$, thus, as the average difference between the observed outcome for the treated as if they were treated ($PO_{t=1}^{T=1}$) and the unobservable, counterfactural outcome for the treated as if they were untreated ($PO_{t=0}^{T=1}$). 

To estimate the $ATT$, we can use observed outcome for all treated units to estimate PO_{t=1}^{T=1}, as the treatment status is observed coherently. By contrast, $PO_{t=0}^{T=1}$ can be estimated as outcome of the untreated, controll group $T=0$ as if it would have the characteristics of the treated group under the assumption of no unobservable counfounding (a.k.a. exogeneity and independence assumption). This counterfactual potential outcome is given by matching the treated and untreated group and producing weights that allow to weight the outcome of the untreated group accordingly:

```math
\begin{equation}
\begin{array}{r c c c c }
	ATT &=& PO_{t=1}^{T=1} &-& PO_{t=0}^{T=1} \\
		&=& \overline{Y}_{t=1,m} &-& \overline{Y}_{t=0^{t=1},m} \\
\end{array}
\end{equation}
```

If group $B$ from our decomposition notation is assinged as the treated group $t=1$ (and group $A$ is untreated $t=0$), we see that $ATT$ equals the unexplained component of $ D_0.

Alternatively, we could also switch the matching direction to estimate the average treatement effect on the controls $ATC$, which gives us:

```math
\begin{equation}
\begin{array}{r c c c c c c }
	ATC &=& PO_{t=1}^{T=0} &-& PO_{t=0}^{T=0} && \\
		&=& \overline{Y}_{t=1^{t=0},m} &-& \overline{Y}_{t=0,m} && \\
		&=& \overline{Y}_{B^A,m} &-& \overline{Y}_{A,m} &=&  D_0 \\
\end{array}
\end{equation}
```

Thus, the $ATC$ equals again the unexplained decomposition-component $D_0$, however, with the reversed matching direction. To indicate which group's characteristics are used to estimate the counterfactual potential outcome, refer to which group is assigned the `xref()`-label; to switch the matching direction, you can use this option to adjust it accordingly. 


Note that the estimation of any treatment effects ($ATT$ as well as $ATC$) via mtaching hinges on the assumption of no unobservable confounding (a.k.a. exogeneity and independence assumption). This is also the reason why matching based treatment effect estimation only provides estimates of $ATT$ or $ATC$ (or $ATEs$) and omits
other components ($D_X$, $D_A$, or $D_B$), which are not of direct interest. While the decomposition framework also operates the same counterfactual as the treatment effect estimation (e.g., what would group's $A$ outcome if it had the characteristics of group $B$), it does not necessarily aim at a causal interpretation (see XXX for a discussion). If the decomposition components are of descriptive interest, the machting-based decomposition does not hinge on the assumption of no unobservable confounding (a.k.a. exogeneity and independence assumption). 




## Post-Estimation Statistics

### Descriptive statistics by machting status (`summarize`)

### Decomposition-components across the distribution of $Y$$ (`gapoverdist`)

### Contribution of characteristics to the components $D_A$ and $D_B$ (`dada`)




## References:
Ñopo, H. (2008). Matching as a Tool to Decompose Wage Gaps. The Review of Economics and Statistics, 90(2), 290–299. https://doi.org/10/b6tqwq

Hamjediers, M. & Sprengholz, M. (2023). Comparing the Incomparable? Issues of Lacking  Common Support, Functional Form Mis-Specification, and Insufficient Sample Size in  Decompositions. Sociological Methodology, 53(2), 344-365. https://doi.org/10.1177/00811750231169729

Iacus, S.M., King, G. & Porro, G. (2012). Causal Inference without Balance Checking: Coarsened Exact Matching. Political Analysis 20, 1–24. https://doi.org/10.1093/pan/mpr013.
