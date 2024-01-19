cap cd "Z:\Projekte\stata_nopo_decomp\nopo_decomposition"

run "nopomatch.ado"
run "nopo.ado"

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
// egen strata = group(age edu)
// qui foreach v in age edu strata {
// 	tab `v', gen(`v'_)
// }

// labels are appropriately captured in matching table
recode t (0 = 0 "M") (1 = 1 "F"), gen(groups)
lab var groups "Groups"
lab var age "Current age in years"
lab def age 1 "18-27" 2 "28-37" 3 "38-47" 4 "48-57" 5 "58-67"
lab val age age
lab var edu "Education"
lab def edu 1 "Edu 1" 2 "Edu 2" 3 "Edu 3" 4 "Edu 4"
lab val edu edu

gen cluster = floor(_n/100) + t * 1000

bys groups: gen n = _n
replace age = 6 if n == 1 & groups == 1
nopo decomp wage i.age i.edu, by(groups) xref(0)
nopo summarize age edu t, label fvdummies
stop

bootstrap, noisily reps(2): nopo decomp wage i.age i.edu, by(groups) xref(0)

ereturn list
nopo decomp wage i.age i.edu, by(groups) xref(1)
ereturn list
nopo decomp wage i.age i.edu, by(groups) bref(0)
ereturn list
nopo decomp wage i.age i.edu, by(groups) bref(1)
ereturn list
stop
kmatch em groups age edu (wage), att replace
kmatch em groups age edu (wage), atc replace
nopomatch age edu, outcome(wage) by(groups) abs replace sd

stop


// nopomatch comparison: SEs somewhere in between suest and bootstrap
eststo clear
nopomatch age edu, outcome(wage) by(groups) abs replace sd
eststo nopomatch
kmatch em groups age edu (wage), atc replace
// eststo kmatch
nopo decomp wage i.age i.edu, by(groups) atc kmkeepgen
eststo nopo_diff
stop
bootstrap, reps(100): nopo decomp wage i.age i.edu, by(groups)
eststo nopo_bs
esttab, compress se rename(DF DA DM DB ATC D0) keep(D0 DX) nodepvars


// nopo_ex 1
//
// kmatch em mbsmoke mage_c fage_c prenatal1 mmarried fbaby foreign alcohol deadkids (bweight), atc
// nopo decomp bweight mage_c fage_c prenatal1 mmarried fbaby foreign alcohol deadkids, by(mbsmoke)
