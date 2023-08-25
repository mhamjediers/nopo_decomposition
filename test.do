cap cd "Z:\Projekte\stata_nopo_decomp\nopo_decomposition"

run "nopo.ado"
run "nopomatch.ado"
run "nopodecomp.ado"

clear
set seed 1234
set obs 1000

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

// labels are appropriately captured in matching table
recode t (0 = 0 "F") (1 = 1 "M"), gen(groups)
lab var groups "Groups"
lab var age "Current age in years"
lab def age 1 "18-27" 2 "28-37" 3 "38-47" 4 "48-57" 5 "58-67"
lab val age age
lab var edu "Education"
lab def edu 1 "Edu 1" 2 "Edu 2" 3 "Edu 3" 4 "Edu 4"
lab val edu edu

nopo decomp wage i.age i.edu, by(groups) dtable
stop


// IA
egen grp = group(age edu)
tab grp, gen(grp_)

// program to get counterfactual outcome and matching set values after 
// teffects nnmatch
cap program drop getcf
program define getcf

    // syntax
    syntax varname [if] [in] [aweight] ///
		, id(varlist max=1) matchid(str) [strata(varlist max=1)] treat(varlist max=1)

    // sample
    tempvar touse
    mark `touse' `if' `in'
	
	// weight
	local weightvar = subinstr("`exp'", "=", "", 1)
	
	dis "Fetching counterfactual outcome"
	
	//quietly {
		
		//preserve
			
			// esample
			keep if `touse'
            
            // save input
            keep `id' `matchid'* `treat' `strata' `varlist' `weightvar'
            reshape long `matchid', i(`id') j(n)
 
            collapse ///
                (mean) `varlist'_cfm = `varlist' ///
                (sd) `varlist'_cfv = `varlist' ///
                [`weight'`exp'] if `treat' == 1, by(`strata')
            replace `varlist'_cfv = `varlist'_cfv^2
            sum `varlist'_cfm
            sum `varlist'_cfv
            stop
            
            tempfile input
            save `input', replace

            // file for matches
            keep if `treat' == 1
            keep `id' `strata' `varlist'
            rename `id' `matchid'
            rename `varlist' `varlist'_cf
            tempfile matches
            save `matches', replace
            
            // gen long dataset for obs and their matches
            use `input', clear
            keep if `treat' == 0
            unab _m : `matchid'*
            egen nonmi = rownonmiss(`_m')
            keep if nonmi > 0
			reshape long `matchid', i(`id') j(n)
            drop if mi(`matchid')
            save `input', replace
            
            // save outcome per id
            use `input', clear
            merge m:1 `matchid' using `matches' // merge cf outcome
            keep if _merge == 3
            drop _merge
            
            // collapse unit
            if ("`strata'" != "") local group = "`strata'"
                else local group = "`id'"

			// collapse & rename counterfactuals
			collapse ///
                (mean) `varlist'_cfm = `varlist'_cf ///
                (sd) `varlist'_cfv = `varlist'_cf ///
                [`weight'`exp'], by(`group')
			replace `varlist'_cfv = `varlist'_cfv^2
            //drop if mi(`varlist'_cfv)
            //replace `varlist'_cfv = 0 if mi(`varlist'_cfv) & !mi(`varlist'_cfm)
            tempfile target
			save `target', replace
            
            noisily sum `varlist'_cfm
            noisily sum `varlist'_cfv
		
		//restore

		if ("`strata'" != "") merge m:1 `strata' using `target', nogen
            else merge 1:1 `id' using `target', nogen
		//drop `matchid'*
	
	//}
	
	dis "Counterfactuals merged by `id' for: `varlist'."

end

// // kmatch
// nopomatch age edu, outcome(wage) by(groups) abs replace sd
// nopo decomp wage i.age i.edu, by(groups) atc kmkeepgen
// //qui kmatch em groups age edu (wage), ate atc att po replace idgen(mid_) wgenerate generate dygen
// //mat list e(b)\
// gen w = 1
// preserve
//     collapse (sum) _wf = w if _nopo_matched == 1 & t==0, by(strata)
//     replace _wf = _wf/470
//     tempfile wf
//     save `wf'
// restore
// preserve
//     keep if _nopo_matched == 1 & t == 1
//     collapse (mean) mwage = wage (sd) vwage = wage, by(strata) 
//     replace vwage = vwage^2
//     merge 1:1 strata using `wf'
//     noisily sum _wf
// 	noisily sum mwage
// 	noisily sum vwage
//     // test for part1
//     quietly gen _part1=(_wf*(1-_wf)*(mwage^2)) / ((1.236842105263158)^2) + vwage*(_wf^2)
//     quietly summ _part1
//     local _total1=(r(sum))/380
//     dis "TOTAL1: `_total1'"
// restore
// stop
// kmatch em groups age edu (wage), atc replace wgenerate generate dygen
// bys t: sum _KM_mw
// bys t: sum _KM_nm _KM_nc
// gen nm = _KM_nm
// gen nc = _KM_nc
// kmatch em groups age edu (wage), att replace wgenerate generate dygen
// bys t: sum _KM_nm _KM_nc
// replace nm = _KM_nm if inlist(nm, ., 1)
// replace nc = _KM_nc if inlist(nm, ., 0)
// bys t: sum nm nc



// kmatch
nopomatch age edu, outcome(wage) by(groups) abs replace sd
nopo decomp wage i.age i.edu, by(groups) atc kmkeepgen tttt
stop
gen id = _n
qui kmatch em groups age edu (wage), atc replace idgen(mid_) wgenerate generate dygen
sort id
mat list e(b)
quietly count if t==1 & _nopo_matched == 1
local _nm = r(N)
quietly count if t==0 & _nopo_matched == 1
local _nf = r(N)
local _alpha = `_nf' / `_nm'
sum wage if t == 0 & _nopo_matched == 1
local _total0 = r(Var) / `_nf'
noisily dis "TOTAL0: `_total0'"

noisily getcf wage if _nopo_matched == 1, id(id) matchid(mid_) strata(strata) treat(t)
//browse id t strata wage wage_cfm wage_cfv cfobs
// sum wage_cfm
// sum wage_cfv
preserve
    gen w = 1
    collapse (sum) _wf = w if _nopo_matched == 1 & t==0, by(strata)
    replace _wf = _wf/470
    tempfile wf
    save `wf'
restore

keep if _nopo_matched == 1 & t == 0
bys strata: gen n = _n
keep if n==1
merge 1:1 strata using `wf'

// sum cfobs if t == 1 & _nopo_matched == 1
// gen _wf = (cfobs / r(mean)) / `_nf'
// keep if t == 1 & _nopo_matched == 1
// gen _wf = cfobs / `_nf'
// sum _wf if t == 1 & _nopo_matched == 1
rename wage_cfm _ym
rename wage_cfv _varym
// // why is the same as the mean of _ym? -> unweighted mean over strata; interesting!
// gen wage_ate = wage if t == 1
// replace wage_ate = wage_cfm if t == 0
// sum wage_ate
// test for part1
sum _wf
sum _ym
sum _varym
noisily table strata, stat(mean _varym)
stop
quietly gen _part1=(_wf*(1-_wf)*(_ym^2)) / ((`_alpha')^2) + _varym*(_wf^2)
quietly summ _part1 if t == 1 & _nopo_matched == 1
local _total1=(r(sum))/`_nm'
dis "TOTAL1: `_total1'"

count if t == 1 & _nopo_matched == 1
local _K = _result(1) /*numero de celdas*/
local _total2 = 0
local _j = 1
gen _wfym = _wf * _ym if t == 1 & _nopo_matched == 1
egen _sumwfym1 = sum(_wfym)
gen _sumwfym2 = _sumwfym1 - _wfym
gen _sumwfym3 = _sumwfym2 * _wfym

sum _sumwfym3
local _total2 = r(sum)/2

quietly drop _wfym _sumwfym1 _sumwfym2 _sumwfym3

local _total2 = 2*(`_total2')/((`_nm')*((`_alpha')^2))
local _dev = sqrt(`_total0'+`_total1'-`_total2')
display in yellow "Std.error DO = " `_dev'

stop

// nopomatch comparison: SEs somewhere in between suest and bootstrap
eststo clear
nopomatch age edu, outcome(wage) by(groups) abs replace sd
eststo nopomatch
// kmatch em groups age edu (wage), atc replace
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
