Setup:

- Treatment $T = [0;1]$
- $Group A = T == 0$ (untreated; control)
- $Group B = T == 1$ (treated)
- $D = \overline{Y}_B - \overline{Y}_A$

Treatment Effects:

- $ATT = D_0 = \overline{Y}_{B,m} - \overline{Y}_{A^B,m}$
- $ATC = D_0 = \overline{Y}_{B^A,m} - \overline{Y}_{A,m}$

$D_0$ is due to differences returns, which remain after _matching and weighting_:

- $ATT$: $D_0$ remains if $A$ had characteritics like $B$ ($\rightarrow A^B$)
- $ATC$: $D_0$ remains if $B$ had characteritics like $A$ ($\rightarrow B^A$)

In the Potential-Outcome-Framework:

- ATT:
  - $ATT = PO_{t=1}^{T=1} - PO_{t=0}^{T=1} = D_0$
  - For $PO_{t=1}^{T=1} = B_m$ to be untreated, we assume it to have $B$ characteristics but $A$ returns ($\rightarrow A^B$)
- ATC:
  - $ATC = PO_{t=1}^{T=0} - PO_{t=0}^{T=0} = D_0$
  - For $PO_{t=0}^{T=0} = A_m$ to be untreated, we assume it to have $A$ characteristics but $B$ returns ($\rightarrow B^A$)


> $D_0$ is the part of $D$ which would remain if $A$ and $B$ had the characteristics of the group specified in `xref()`, which can either be $A$ or $B$ (only one group has to be reweighted so that it's characteristics correspond to the respective other group).