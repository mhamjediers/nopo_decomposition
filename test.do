//cd "Z:\Projekte\stata_nopo_decomp\nopo_decomposition"

run "nopo.ado"
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
//replace t = 1 if age==1 & edu==1

// wage
gen wage = 5*t + 0.5*age + edu + age*edu - (age*edu*t/5) + runiform()
bys t: sum wage, d
lab var wage "Hourly wage"

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

// nopomatch comparison
nopomatch age edu, outcome(wage) by(groups) replace abs sd

//
// Standalone
//

// see swap and bref
nopo decomp wage age edu, by(groups) norm
nopo decomp wage age edu, by(groups) swap norm
nopo decomp wage age edu, by(groups) swap bref(groups == 1) norm
nopo decomp wage age edu, by(groups) bref(groups == 0) norm
stop

// see noisily and passthru options
nopo decomp wage age edu, by(groups) kmpassthru(_N generate k_omit bwidth) kmnoisily
ereturn list

// see subcmd changes
nopo decomp wage age edu, by(groups) kmatch(ps)

//
// Post kmatch
//

kmatch em groups age edu (wage), tval(1) atc att generate wgenerate replace
nopo decomp // defaults to ATT if present, uses ATC if only ATC estimates present

kmatch ps groups age edu (wage), tval(1) atc att bwidth(0.5) generate wgenerate replace
nopo decomp

//
// Postestimation (see check estimates, these can be omitted after testing)
//

nopo decomp wage age edu, by(groups)
tempfile test
nopo dadb edu, save(`test')
nopo gapoverdist
nopo summarize, label // columns of unmactched omitted if not present (better than all mi, I think)

/* Standalone to dos: ALL DONE

*Change the atc/att naming
*The default should be reference in gap estimation is also reference in vector (now atc)
*Swap option flips group-indicator (and thereby both references)
*Reference option allows groups == # to indicate which reference vector 

*for standalone option kmatch(em|ps|md) (default is em)
* and potentially some of the further options kmatch_options(...)

*and ereturn internally called kmatch line for potential adjustment (with atc and att); then everybody can work with this is they want to
*scalar passthrough

*This should be somehow the final result:
nopo decomp wage age edu, by(groups) swap ref(groups == 1) kmatch(md) kmatch_options(bw(0.2)) noisily passthrough(bw)

*noisily option to display kmatch-output with all its specifications/bandwith  --> add the notable and nose options, as these should not be part of it

*/


/* Output to dos: ALL DONE

*output for summarize via matlist (and option with labels for rownames)
*Output tables should be shown as in previous version (indicating group A/B and which one is reference)
*Show both tables also when used as post-estimation
*indicate also the matching algorithm in output
*indicate number of strata in exact matching or some other scalars(bw?) when other matching algorithm

*/

*bootstrap: nopo decomp wage age edu, by(groups)
