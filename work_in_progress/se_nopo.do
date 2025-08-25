cap cd "Z:\Projekte\stata_nopo_decomp\nopo_decomposition"

//run "work_in_progress/nopomatch.ado" // code below already implemented here

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

/*
 Intution behind computation is:

 SE = sqrt(   Var(PO_f)/N_m   + Var(Y_f)/N_f )
    = sqrt( (total1 - total2) + total0 )

 Computation via delta method of first part (if I understood correctly) based on the product
 of m wage x f weighting factor:

 https://stats.stackexchange.com/questions/62916/confidence-interval-for-the-product-of-two-parameters

 Var =   (P1^2) * VarP2
       + (P2^2) * VarP1
       + 2 * P1 * P2 * CoV 
       (last term: covariance is OK, but P1*P2 is not: 
        looks like variance estimation for weighted sum of variables, but that does not fit
        with the delta method for a product...what to do?!)


 And by Nopos definitions:
 ------------------------
 MATCHED ONLY; PER STRATUM

 P1 = _ym (male mean wage)
 P2 = _wf (weighting factor; share of women in stratum)

 VarP1 = _varym (variance of male wages in stratum)
 CovP1 across strata = 0
 
 VarP2 = _wf(1 - _wf) / alpha^2
 CovP2 across strata (i,j) = -(_wf_i x _wf_j) / alpha^2
 
 => across strata K:

 Var(P0_f) =   SUM(i=1 -> K) [ (_wf(1 - _wf) / alpha^2) * _ym^2 ]
             + SUM(i=1 -> K) [ _varym * _wf^2 ]
             + 2 * SUM(i=1 -> K) [ 
                SUM (j=1 -> K; j != i) [ (-(_wf_i x _wf_j) / alpha^2) * _ym_i * _ym_j ] 
                ]
 

*/

// get matching indicator from nopo; keep support
nopomatch age edu, outcome(wage) by(t) replace sd
egen strata = group(age edu)

qui {
    
    // weight 1
    tempvar _fexp
    quietly gen `_fexp' = 1
    
    // normalize outcome to be comparable to original ado
    sum wage if t == 0 [iw =`_fexp']
    gen _rwage = wage / r(mean)
    
    // keep support
    keep if _supp == 2

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

    /*2. The sample average of earnings for males that exhibit the set of characteristics*/
    /*the population variance of male wages that exhibit the set of characteristics*/

    preserve
        keep if t == 1
        collapse (mean) _ym = _rwage (sd) _varym = _rwage [iw=`_fexp'], by (strata)
        replace _varym = (_varym)^2
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
