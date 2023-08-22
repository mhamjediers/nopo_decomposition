cap cd "Z:\Projekte\stata_nopo_decomp\nopo_decomposition"

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

replace wage = round(wage) if t == 1

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
recode t (0 = 0 "Native Men") (1 = 1 "Immigrant Women"), gen(groups)
lab var groups "Groups"
lab var age "Current age in years"
lab def age 1 "18-27" 2 "28-37" 3 "38-47" 4 "48-57" 5 "58-67"
lab val age age
lab var edu "Education"
lab def edu 1 "Edu 1" 2 "Edu 2" 3 "Edu 3" 4 "Edu 4"
lab val edu edu

// IA
egen grp = group(age edu)
tab grp, gen(grp_)

// nopomatch comparison
nopomatch age edu, outcome(wage) by(groups) replace abs sd
eststo nm
// according to lit, w(1) would be standard given that for oaxaca D = A - B (gaps expressed positively)
// oaxaca wage grp_* if _supp == 2, by(groups) w(1) nodetail relax // weight: groups == 1
// eststo o_d_w0
// oaxaca wage grp_* if _supp == 2, by(groups) w(0) nodetail relax // weight: groups == 1
// eststo o_d_w1
// oaxaca wage grp_* if _supp == 2, by(groups) swap w(1) nodetail relax // weight: groups == 1
// eststo o_swap_w0
// oaxaca wage grp_* if _supp == 2, by(groups) swap w(0) nodetail relax // weight: groups == 1
// eststo o_swap_w1
//oaxaca wage grp_* if _supp == 2, by(groups) w(0) nodetail relax // weight: groups == 1
//kmatch em groups age edu (wage), nate ate atc att po generate wgenerate replace
nopo decomp wage i.age i.edu, by(groups)
nopo decomp wage i.age i.edu, by(groups) swap
stop
nopo decomp wage i.age i.edu, by(groups) xref(0) norm
nopo decomp wage i.age i.edu, by(groups) swap norm
nopo decomp wage i.age i.edu, by(groups) swap xref(1) norm
stop

//
// nopo summarize implementation with factors and labels
//

// if md or ps: use factor notation in nopo call
nopo decomp wage i.age i.edu, by(groups) kmatch(ps)
nopo summarize
nopo summarize, label

// or use factor notation directly in nopo summarize
// i guess the only usecase would be dummies, otherwise misspecified model for ps?
// but ok for md?
nopo decomp wage age edu, by(groups) kmatch(ps)
nopo summarize i.edu, label

//
// Plot options
//

// only xsize and ysize directly adjustable (should allow to correct most scale issues)
// the other options per plot are returned in r()
// either
// - you adjust these returns and feed them back to the plot command or
// - build your own plot via `save' and these returns
nopo decomp wage i.age i.edu, by(groups)
nopo gapoverdist
return list
nopo dadb edu
return list

stop

//
// Standalone
//

// see swap and bref
nopo decomp wage i.age edu, by(groups)
nopo decomp wage age edu, by(groups) swap
nopo decomp wage age edu, by(groups) swap bref(groups == 1)
nopo decomp wage age edu, by(groups) bref(groups == 0)

// see normalize
nopo decomp wage age edu, by(groups) norm

// see noisily and passthru options
nopo decomp wage age edu, by(groups) kmpassthru(_N generate k_omit bwidth) kmnoisily
ereturn list

// see subcmd changes
nopo decomp wage age edu, by(groups) kmatch(ps)
nopo decomp wage age edu, by(groups) kmatch(md)
nopo decomp wage age edu, by(groups) kmatch(md) kmopts(nn(4))
nopo decomp wage age edu, by(groups) kmatch(md) kmopts(ridge)

// see kmkeepgen
nopo decomp wage age edu, by(groups) kmatch(ps) kmkeepgen

//
// Post kmatch
//

kmatch md groups age edu (wage), tval(1) atc att generate wgenerate replace nn(1) wor
nopo decomp // defaults to ATT if present, uses ATC if only ATC estimates present

kmatch ps groups age edu (wage), tval(1) atc att bwidth(0.5) generate wgenerate replace
nopo decomp

// see error if generated variables missing from kmatch
// there is no way around this bc. of user intervention
// theoretically it is also possible that the generated vars are there but
// do not correspond to the estimates which are restored
kmatch md t age edu (wage), atc generate wgenerate replace
estimates store km
nopo decomp
estimates restore km
nopo decomp

//
// Postestimation (see check estimates, these can be omitted after testing)
//

tempfile test
nopo dadb edu, save(`test')
nopo gapoverdist, save(`test')
nopo summarize, label // columns of unmactched omitted if not present (better than all mi, I think)

*bootstrap: nopo decomp wage age edu, by(groups)


