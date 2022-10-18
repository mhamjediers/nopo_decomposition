cap program drop nopodecomp
program define nopodecomp, eclass
	version 17
	syntax varlist [if] , BY(varname) 
		
	quietly {

	gettoken y varlist:varlist
	
	tempfile __start
	save `__start', replace
	
	if "`if'" != "" {
		keep `if'
	}
	
	count if !mi(`y') // adjust for whole varlist
	local N = `r(N)'
	
	ttest  `y', by(`by')
	scalar gap = `r(mu_2)' - `r(mu_1)'

	preserve 
		keep if `by'
		gen y_1 = `y' 
		collapse (count) N_1 = y_1 (mean) y_1, by(`varlist')
		tempfile `by'
		save ``by'', replace
	restore

	
	keep if !`by'
	gen y_0 =  `y'

	bys `varlist': egen N_0 = count(y_0)
	bys `varlist': gen i = _n
	keep `varlist' y_0 N_0 i 

	merge m:1 `varlist' using ``by'' 

	***
	count if _merge == 3 //machted from ref_`by'oup
	local m_a = `r(N)'
	count if _merge == 1 //unmachted from ref_`by'oup
	local um_a = `r(N)'
	qui sum  N_1 if _merge == 2 //unmachted from match_`by'oup
	local um_b = `r(sum)'
	qui sum N_1 if _merge == 3 & i == 1 //machted from match_`by'oup
	local m_b = `r(sum)'
	
	scalar N_a = `m_a' + `um_a'
	scalar r_a = `m_a' / N_a
	scalar N_b = `m_b' + `um_b'
	scalar r_b = `m_b' / N_b

	*Weighting is already implied by collapsing
	ttest y_0 == y_1 // unexplained gap

	scalar d_0 = `r(mu_2)' - `r(mu_1)' 

	if `um_a' != 0 {
		ttest y_0 if inlist(_merge,1,3), by(_merge)
		scalar d_a = (`r(mu_2)' - `r(mu_1)') * (`r(N_1)' / (`r(N_2)' + `r(N_1)'))
	}
	else {
		scalar d_a = 0
	}
	
	if `um_b' != 0 {
		bys `varlist': gen h = _n
		reg y_1  ib2._merge [fw = N_1] if h == 1 & inlist(_merge,2,3)
		scalar d_b = -(_b[3._merge]) * (`um_b' / (`um_b' + `m_b')) 
		}
	else {
		scalar d_b = 0
	}

	scalar d_1 = gap - d_0 - d_a - d_b
	
	
	matrix b = [gap, d_0, d_1, d_a, d_b]
	matrix colnames b = "Gap" "Unexplained" "Explained" "Unmatched A" "Unmatched B"
	
	matrix N = [N_a, r_a, N_b, r_b]
	matrix colnames N = "N_a" "%m_a" "N_b" "%m_b"
	
	ereturn post  b, depname("`y'") obs(`N')  properties(b)
	
	use `__start', clear 
}
ereturn display
end
