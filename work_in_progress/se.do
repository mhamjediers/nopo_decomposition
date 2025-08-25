cap cd "Z:\Projekte\stata_nopo_decomp\nopo_decomposition"

// run "work_in_progress/nopomatch.ado"

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

// labels
lab var t "T"
lab def t 0 "0 [F]" 1 "1 [M]"
lab val t t

// program to get counterfactual outcome mean and variance per id / stratum
/*
 Since we do not have strata to collapse both group by, we need to merge
 the matches by stratum first, then compute mean and variance of outcome.
 Pretty slow!
 */
cap program drop getcf
program define getcf

    // syntax
    syntax varname [if] [in] [aweight] ///
		, id(varlist max=1) matchid(str) treat(varlist max=1) strata(varlist)
    
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
            keep `id' `matchid'* `treat' `varlist' `w' `strata'
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
            
            // get outcome mean and variance per stratum
            // which can be unique for ps/md
            collapse ///
                (mean) `varlist' `w' `varlist'_cfm = `varlist'_cf ///
                (sd) `varlist'_cfv = `varlist'_cf ///
                [aw=`w'], by(`strata')
            replace `varlist'_cfv = `varlist'_cfv^2 // variance
            replace `varlist'_cfv = 0 if mi(`varlist'_cfv) & !mi(`varlist'_cfm)
            tempfile target
			save `target', replace
		
		restore

		merge m:1 `strata' using `target', nogen
		drop `matchid'*
	
	}
	
    // dis
	dis "Counterfactuals merged by `id' for: `varlist'."

end

// continuous + MD
gen agec = round((age + runiform()) * 10)
egen strata = group(agec edu)

// kmatch
gen id = _n
qui kmatch md t agec edu (wage), atc replace idgen(mid_) wgenerate generate bwidth(0.5)
ereturn display

qui {
    // weight 1
    tempvar _fexp
    quietly gen `_fexp' = 1

    // normalize outcome to be comparable to original ado
    sum wage if t == 0 [iw =`_fexp']
    gen _rwage = wage / r(mean)

    // keep support
    keep if (_KM_nm > 0 & !mi(_KM_nm)) | (_KM_nc > 0 & !mi(_KM_nc))

    // group n
    quietly count if t == 1
    local _nm = _result(1)
    quietly count if t == 0
    local _nf = _result(1)

    // alpha: group w sum
    quietly summ `_fexp' if t == 1
    local _Nm = r(sum)
    quietly summ `_fexp' if t == 0
    local _Nf = r(sum)
    local _alpha = `_Nf'/`_Nm'

    // total0: Variance of the second term of the right hand side of equation 9
    sum _rwage if t == 0 
    local _total0 = r(Var) / `_nf'
    noisily dis "TOTAL0: `_total0'"

    // DX should be analogous, just the equation reversed and total0 calculated for males
    quietly summ _rwage if t == 1 [iw=`_fexp']
    local _total0_m = _result(4)/`_nf'

    /*Now I construct the first term*/
    /*1. The sample proportion of females that exhibit the set of caracteristics*/

    // wf
    preserve
        quietly keep if t == 0
        collapse (sum) `_fexp', by(strata)
        quietly gen _wf = `_fexp'/`_Nf'
        tempfile wf
        quietly save `wf', replace
    restore

    // getcf
    noisily getcf _rwage, id(id) matchid(mid_) treat(t) strata(strata) 

    // cf mean var
    preserve
        keep if t == 1
        collapse _ym = _rwage_cfm _varym = _rwage_cfv, by(strata)
        merge 1:1 strata using `wf'
        tab _merge
        keep if _merge==3
        drop _merge

        // part 1
        quietly gen _part1 = (_wf*(1-_wf)*(_ym^2))/((`_alpha')^2) + _varym*(_wf^2)
        quietly summ _part1
        local _total1=(r(sum))/`_nm'

        // try to compute total2 by yourself: sum up covariances
        qui count
        local k = r(N)
        scalar _t2 = 0
        forvalues i = 1/`k' {
            local d = `i' + 1
            forvalues j = `d' / `k' {
                scalar _t2 = _t2 + (_wf[`i'] * _wf[`j'] * _ym[`i'] * _ym[`j']) / (`_alpha'^2)
            }
        }
        scalar _t2 = 2 * _t2 / `_nm'
        noi dis "total2 by hand:" _t2

        // original computation
        quietly gen _wfym=_wf*_ym
        quietly egen _sumwfym1 = sum(_wfym)
        quietly gen _sumwfym2 = _sumwfym1 - _wfym
        quietly gen _sumwfym3 = _sumwfym2*_wfym

        quietly sum _sumwfym3
        local _total2 = r(sum)/2

        quietly drop _wfym _sumwfym1 _sumwfym2 _sumwfym3

        local _total2 = 2*(`_total2')/((`_nm')*((`_alpha')^2))
        noi dis "total2 original:" `_total2'

        local _dev = sqrt(`_total0'+`_total1'-`_total2')
        noi dis "Std.error DO = " `_dev'
        noi dis "Std.error DX = " sqrt(`_total0_m'+`_total1'-`_total2')

    restore

}
