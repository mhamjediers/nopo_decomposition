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
qui bys t: sum wage, d
lab var wage "Hourly wage"

// dummies/strata
egen strata = group(age edu)
qui foreach v in age edu strata {
	tab `v', gen(`v'_)
}

lab var t "T"
lab def t 0 "0 [F]" 1 "1 [M]"
lab val t t

// program to get counterfactual outcome mean and variance per id / stratum
cap program drop getcf
program define getcf

    // syntax
    syntax varname [if] [in] [aweight] ///
		, id(varlist max=1) matchid(str) treat(varlist max=1)
    
    // dis
    dis "Fetching counterfactual outcome"
    
    // sample
    tempvar touse
    mark `touse' `if' `in'
	
	// weight
    tempvar w
    if ("`exp'" != "") gen `w' `exp'
        else gen `w' = 1
	
	quietly {
		
		preserve
			
			// esample
			keep if `touse'

            // save input
            keep `id' `matchid'* `treat' `varlist' `w'
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
            
            // save outcome per id
            merge m:1 `matchid' using `matches' // merge cf outcome
            keep if _merge == 3
            drop _merge
            
            // get outcome mean and var per id
            // further collapse by stratum if present
            collapse ///
                (mean) `varlist' `w' `varlist'_cfm = `varlist'_cf ///
                (sd) `varlist'_cfv = `varlist'_cf ///
                [aw=`w'], by(`id' `strata')
            replace `varlist'_cfv = `varlist'_cfv^2
            
            noisily sum `varlist'_cfm
            noisily sum `varlist'_cfv
            
            tempvar strata
            egen `strata' = group(`varlist'_cfm `varlist'_cfv)
            levelsof `strata'
            noisily dis r(levels)
            //if ("`strata'" != "") {
                collapse (mean) `varlist'_cfm `varlist'_cfv [aw=`w'], by(`strata')
			//}
            replace `varlist'_cfv = 0 if mi(`varlist'_cfv) & !mi(`varlist'_cfm)
            tempfile target
			save `target', replace
            
            noisily sum `varlist'_cfm
            noisily sum `varlist'_cfv
		
		restore

		if ("`strata'" != "") merge m:1 `strata' using `target', nogen
            else merge 1:1 `id' using `target', nogen
		drop `matchid'*
	
	}
	
    // dis
	dis "Counterfactuals merged by `id' for: `varlist'."

end

// kmatch
nopomatch age edu, outcome(wage) by(t) abs replace sd
// qui nopo decomp wage i.age i.edu, by(t) atc kmkeepgen
gen id = _n
//qui kmatch em t age edu (wage), atc replace idgen(mid_) wgenerate generate
kmatch em t age edu (wage), ate atc att po replace dygenerate wgenerate generate

bys t _supp: sum _DY_wage

table strata t, stat(mean _KM_nc _KM_nm)



mat list e(b)
mat list e(V)
quietly count if t==1 & _nopo_matched == 1
local _nm = r(N)
quietly count if t==0 & _nopo_matched == 1
local _nf = r(N)

// alpha
local _alpha = `_nf' / `_nm'
sum wage if t == 0 & _nopo_matched == 1

// total0
local _total0 = r(Var) / `_nf'
noisily dis "TOTAL0: `_total0'"

// getcf
noisily getcf wage if _nopo_matched == 1, id(id) matchid(mid_) treat(t) // strata(strata) 

// wf
// preserve    
//     keep if t == 0 & _nopo_matched == 1
//     gen w = 1
//     collapse (sum) _wf = w, by(strata)
//     replace _wf = _wf/`_nf'
//     tempfile wf
//     save `wf'
// restore

// cf mean var
preserve
    keep if t == 0 & _nopo_matched == 1
    //collapse _ym = wage_cfm _varym = wage_cfv, by(strata)
    //merge 1:1 strata using `wf'
    rename wage_cfm _ym
    rename wage_cfv _varym
    //gen _wf = 1/`_nf'
    sum _KM_nc
    gen _wf = (_KM_nc / r(sum)) / `_nf'
    sum _wf
    sum _ym
    sum _varym

    gen _part1=(_wf*(1-_wf)*(_ym^2))/((`_alpha')^2)+_varym*(_wf^2)
    sum _part1
    local _total1=(r(sum))/`_nm'
    dis "TOTAL1: `_total1'"

    gen _wfym = _wf * _ym
    egen _sumwfym1 = sum(_wfym)
    gen _sumwfym2 = _sumwfym1 - _wfym
    gen _sumwfym3 = _sumwfym2 * _wfym

    sum _sumwfym3
    local _total2 = r(sum)/2

    quietly drop _wfym _sumwfym1 _sumwfym2 _sumwfym3

    local _total2 = 2*(`_total2')/((`_nm')*((`_alpha')^2))
    local _dev = sqrt(`_total0'+`_total1'-`_total2')
    display in yellow "Std.error DO = " `_dev'

restore


//     Variable |        Obs        Mean    Std. dev.       Min        Max
// -------------+---------------------------------------------------------
//          _wf |         17    .0588235    .0287622   .0021277   .0957447
//
//     Variable |        Obs        Mean    Std. dev.       Min        Max
// -------------+---------------------------------------------------------
//          _ym |         17    13.93009    4.083944   7.824759   22.98011
//
//     Variable |        Obs        Mean    Std. dev.       Min        Max
// -------------+---------------------------------------------------------
//       _varym |         17    .0843229    .0301279   .0586971   .1890679
//
// Std.error DO = .25941368
