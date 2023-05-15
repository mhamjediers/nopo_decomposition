//cd "Z:\Projekte\stata_nopo_decomp\nopo_decomposition"

run "postestimation.ado"
run "nopodecomp.ado"

clear
set seed 1234
set obs 10000

gen age = _n - 5*floor((_n-1)/5)
gen edu = _n - 4*floor((_n-1)/4)

// decomposition comparison
gen t = -0.25*edu - 0.1*age + runiform()
qui sum t, d
replace t = t > 0.7 * `r(mean)'
replace t = 1 if age==1 & edu==1

// wage
gen wage = 5*t + 0.5*age + edu + age*edu - (age*edu*t/5) + runiform()
bys t: sum wage, d

// norm. / log.
sum wage if t==0, meanonly

// dummies/strata
egen strata = group(age edu)
qui foreach v in age edu strata {
	tab `v', gen(`v'_)
}

// timer clear
//
// // mi check: DA/DB AND N UNMATCHED NOT CORRECT!
// // see comment in nopodecomp; same with nopomatch, see below
// // compare to nopoplot2 output below (which is correct I believe)
// // needs `touse' handling
// //replace age = . if edu == 3 & t == 1
//
// // nopo
// timer on 1
//nopomatch age edu, outcome(wage) by(t) replace abs
// nopomatch age edu if !mi(age), outcome(wage) by(t) replace abs
// timer off 1
//
// // new command
// timer on 2
// nopodecomp wage age edu, by(t)  //normalize replace
// timer off 2
//
// timer list
//

// labels are appropriately captured in matching table
recode t (0 = 0 "Immigrant women") (1 = 4 "Native men"), gen(groups)
lab var groups "Groups"
lab var edu "Edu"
lab def edu 1 "Edu 1" 2 "Edu 2" 3 "Edu 3" 4 "Edu 4"
lab val edu edu


//nopomatch age edu, outcome(wage) by(groups) replace abs
kmatch em groups age edu (wage), ate atc att wgenerate generate tval(4)
qui estimates store kmatch
nopopost decomp, att
qui estimates restore kmatch
nopopost decomp, atc
qui kmatch em groups age edu (wage), ate atc att wgenerate generate tval(0) replace
qui estimates store kmatch
nopopost decomp, att
qui estimates restore kmatch
nopopost decomp, atc
tempfile test
nopopost dadb edu, save(`test')
nopopost gapoverdist
nopopost summarize
mat list r(npsum)
ereturn list
