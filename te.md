# Methodological Details

## Ñopo's (2008) matching decomposition:

We are interested in the decomposition of the raw gap $D$ in outcome $Y$ between group A and group B

$$
D = \overline{Y}_B - \overline{Y}_A \quad ,
$$

into an “explained” component that is based on differences between groups in characteristics that predict $Y$ and a remaining “unexplained” component. However, we do not want to compare apples and oranges, so that the decomposition also needs to calculate which part of $D$ is due to units in both groups which have combinations of predictive characteristics out of common support (characteristics which do not occur in both groups).

Ñopo's matching achieves these goals as follows. In a one-to-many exact matching, each individual from group $B$ is matched to all individuals from group $A$ with the same combination of characteristics (each unique combination of characteristics represents one stratum). The matching flags for all observations of groups $A$ and $B$ if their characteristics are in common support (matched $m$) or out of common support (unmatched $u$). Over common support, the ratio of $B$ to $A$ units in each stratum can be used to create a reweighted group $A^B$ which has the exact same distribution across all strata as group $B$. The outcome of this counterfactual group $A^B$ can be interpreted in two ways: (1) as the average outcome of group $A$ if it had the same characteristics as group $B$ and (2) as the average outcome of group $B$ if it had the same returns to characteristics as group $A$.

The overall gap can then be additively decomposed into four parts:

$$\begin{array}{rcccccc}
D &=& D_0  &+& \rlap{$\overbrace{\phantom{\qquad D_X + \quad D_A + D_B}}^{\text{compositional difference}}$} D_X &+&  D_A +D_B \\[5pt]
&=& \overbrace{\overline{Y}_{B,m} - \rlap{$\underbrace{\phantom{\overline{Y}_{A^B,m} + \quad \overline{Y}_{A^B,m}}}_{\mathclap{\substack{\text{splitting difference} \\ \text{among matched by} \\ \text{reweighted group A}}}}$} \overline{Y}_{A^B,m}} &+& \overbrace{\overline{Y}_{A^B,m} - \overline{Y}_{A,m}} &+& \underbrace{D_A + D_B}_{\mathclap{\substack{\text{out of} \\ \text{support}}}}
\end{array}$$

$D_X$ is the average gap between the matched units of re-weighted group $A^B$ and the matched units of group $A$, which is explained by the fact that groups $A$ and $B$ are differently distributed across matched strata (some sets of characteristics being more likely in one group than the other). 
By contrast, $D_0$ is the average gap between matched units of group $B$ and the re-weighted group $A^B$. Since $B$ and $A^B$ are equally distributed across matched strata, $D_0$ captures how much of the raw gap remains unexplained by differences in the considered characteristics. $D_0$ and $D_X$ are analogous to the components of a twofold KBO decomposition, but they only pertain to matched units. When compositional differences between groups limit common support, the effect that unmatched individuals in both groups have on the outcome gap is captured by the components $D_A$ and $D_B$:

$$\begin{equation}
	D_A =  \underbrace{(\overline{Y}_{A,m} - \overline{Y}_{A,u})}_{\mathclap{\substack{\text{gap between matched} \\ \text{and unmatched }A}}} \cdot \underbrace{(N_{A,u}/N_A)}_{\mathclap{\substack{\text{share of} \\ \text{unmatched }A}}}  \quad\quad
	D_B =  \underbrace{(\overline{Y}_{B,u} - \overline{Y}_{B,m})}_{\mathclap{\substack{\text{gap between unmatched} \\ \text{and matched }B}}} \cdot \underbrace{(N_{B,u}/N_B)}_{\mathclap{\substack{\text{share of} \\ \text{unmatched }B}}} 
\end{equation}$$

$D_A$ is the gap between the averages of the outcome $\overline{Y}$ for the unmatched $u$ and matched $m$ units within group $A$, weighted by the frequency of unmatched $A$ units ($N_{A,u}$) in relation to the overall size of group $A$ ($N_{A}$), so that $D_A$ approaches zero with fewer observations out of common support. $D_A$ denotes how much of the raw gap is due to unmatched $A$ units having higher or lower values in the outcome than matched $A$ units, where $D_A < 0$ if the outcome is lower among the matched and $D_A > 0$ if the outcome is lower among the unmatched (reversed for $D_B$).

## Terminology: The relationship between $D_0$, $ATT$ and $ATC$

$D_0$ corresponds to either $ATT$ (Average Treatment effect on the Treated) or $ATC$ (Average Treatment effect on the Unreated/Control) from the treatment effects literature. Given the following setup

- Treatment $T[0;1]$
- $Group A = T == 0$ (untreated; control)
- $Group B = T == 1$ (treated) $\quad ,$

$ATT$ and $ATC$ are given by

- $ATT = D_0 = \overline{Y}_{B,m} - \overline{Y}_{A^B,m}$
- $ATC = D_0 = \overline{Y}_{B^A,m} - \overline{Y}_{A,m} \quad .$

In both cases, $D_0$ is due to differences returns, which remain after _matching and weighting:_

- $ATT$: $D_0$ remains if $A$ had characteristics like $B$ ($\rightarrow A^B$)
- $ATC$: $D_0$ remains if $B$ had characteristics like $A$ ($\rightarrow B^A$)

In the Potential-Outcome-Framework:

- ATT:
  - $ATT = PO_{t=1}^{T=1} - PO_{t=0}^{T=1} = D_0$
  - For $PO_{t=1}^{T=1} = B_m$ to be untreated, we assume it to have $B$ characteristics but $A$ returns ($\rightarrow A^B$)
- ATC:
  - $ATC = PO_{t=1}^{T=0} - PO_{t=0}^{T=0} = D_0$
  - For $PO_{t=0}^{T=0} = A_m$ to be treated, we assume it to have $A$ characteristics but $B$ returns ($\rightarrow B^A$)

Thus, $D_0$ is the part of $D$ which would remain in the counterfactual setting of $A$ and $B$ having the characteristics of the group specified in option `xref()`, which can either be $A$ or $B$. Please note that the reference in terms of _returns to characteristics_ would always be the respective other group (in the (Kitagawa-)Blinder-Oaxaca parlance, this would be the coefficient-vector).

Ñopo, H. (2008). Matching as a Tool to Decompose Wage Gaps. The Review of Economics and Statistics, 90(2), 290–299. https://doi.org/10/b6tqwq
