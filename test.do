cd "Z:\Projekte\stata_nopo_decomp\nopo_decomposition"

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
recode t (0 = 0 "Immigrant women") (1 = 1 "Native men"), gen(groups)
lab var groups "Groups"
lab var edu "Edu"
lab def edu 1 "Edu 1" 2 "Edu 2" 3 "Edu 3" 4 "Edu 4"
lab val edu edu


*Standalone
nopopost decomp wage age edu, by(groups) 

/* Standalone to dos: 

*Change the atc/att naming
*The default should be reference in gap estimation is also reference in vector (now atc)
*Swap option flips group-indicator (and thereby both references)
*Reference option allows groups == # to indicate which reference vector 

*for standalone option kmatch(em|ps|md) (default is em)
* and potentially some of the further options kmatch_options(...)
*and ereturn internally called kmatch line for potential adjustment (with atc and att); then everybody can work with this is they want to
*noisily option to display kmatch-output with all its specifications/bandwith  --> add the notable and nose options, as these should not be part of it
*scalar passthrough

*This should be somehow the final result:
nopopost decomp wage age edu, by(groups) swap ref(groups == 1) kmatch(md) kmatch_options(bw(0.2)) noisly passthrough(bw)

*/



/* Output to dos:

*Output tables should be shown as in previous version (indicating group A/B and which one is reference)
*Show both tables also when used as post-estimation
*indicate also the matching algorithm in output
*indicate number of strata in exact matching or some other scalars(bw?) when other matching algorithm

*output for summarize via matlist (and option with labels for rownames)
*/

nopomatch age edu, outcome(wage) by(groups) replace abs sd
kmatch em groups age edu (wage), ate atc att wgenerate generate tval(1)
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


*bootstrap: nopopost decomp wage age edu, by(groups)
