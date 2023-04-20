//cd "Z:\Projekte\stata_nopo_decomp\nopo_decomposition"

run "postestimation.ado"
run "nopodecomp.ado"

clear
set seed 1234
set sortseed 1234
set obs 10000


gen age = _n - 5*floor((_n-1)/5)
gen edu = _n - 4*floor((_n-1)/4)

gen x = runiform(1, 10000)

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

timer clear

// nopo
timer on 1
nopomatch age edu, outcome(wage) by(t) replace abs
timer off 1

// new command
timer on 2
nopodecomp wage age edu, by(t) replace //normalize replace
timer off 2

ereturn list

timer list

// bootstrapping
nopodecomp wage age edu, by(t) replace bootstrap bsopts(seed(123))
ereturn list


// labels are appropriately captured in matching table
recode t (0 = 0 "immigrant women") (1 = 4 "native men"), gen(groups)
nopodecomp wage age edu, by(groups) swap prefix(new) // normalize replace swap prefix

ereturn list

// Post-Estimation commands
nopoplot
