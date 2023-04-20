
*********
*
*********
cap program drop nopodecomp
program define nopodecomp , eclass
	syntax varlist , BY(varname) [PREFix(string) REPLACE SWITCH NORMalize]
	

	marksample touse
	gettoken outcome match_set:varlist

	qui {
		
		*handle by-variable (indicator `gr', get labels, assert 2-levels, allow reversed)
		tempvar gr
		gen `gr' = .
		local lbe : value label `by'
		levelsof `by', local(levels)
		assert `r(r)' == 2
		local i = 0
		if "`switch'" != "" { // reverse order for reference group
			local levels = strreverse("`levels'")
		}
		foreach l of local levels {
			replace `gr' = `i' if `by' == `l'
			local ref`i' = "`by'==`l'"
			if "`lbe'" !=  "" {
				local gr`i++' : label `lbe' `l'
			}
			else {
				local gr`i++' = "`by'==`l'"
			}
		}
		
	
		tempvar y 
		gen `y' = `outcome' if !mi(`gr')
	
	
		*Option replace
		if "`replace'" != "" {
			cap drop `prefix'_strata
			cap drop `prefix'_matched
			cap drop `prefix'_weights
		}
		
		egen `prefix'_strata = group(`match_set')
		lab var `prefix'_strata "strata of nopo decomp."
		gen `prefix'_matched = .
		lab var `prefix'_matched "matching ind. of nopo decomp."
		gen `prefix'_weights = .
		lab var `prefix'_weights "weights of nopo decomp."
	}
	
	*Matching
	qui nopo_match `y' `gr' `prefix'_strata  `prefix'_matched `prefix'_weights
	
	*Result table of matching
	di as text _newline
	di as text "Exact matching for Nopo decomposition:"	_col(60) "N(strata) = " 		_col(85) nstrata
	di as text " " 										_col(60) "N(matched strata) = " _col(85) mstrata
	di as text _newline

	di as text " Group " _col(20) "{c |}" /*
		*/ _col(25) "Matched" _col(41) "{c |}" /*
		*/ _col(46) "Unmatched"  _col(62) "{c |}" /*
		*/ _col(67) "Total"  _col(73) "{c |}" /*
		*/ _col(76) "Mean(" abbrev("`outcome'",9) ")"
	di as text "{hline 19}{c +}{hline 20}{c +}{hline 20}{c +}{hline 10}{c +}{hline 15}"
	di as text abbrev("`gr0'",13) " (ref)"  _col(20) "{c |}" /*
		*/ as result _col(21) %8.0g `=mtab[1,1]' _col(30) "("%4.1f `=mtab[1,2]*100' "%)" _col(41) "{c |}" /*
		*/ _col(43) %8.0g `=mtab[1,3]' _col(52) "("%4.1f `=mtab[1,4]*100' "%)"  _col(62) "{c |}" /*
		*/ _col(64) %8.0g `=mtab[1,5]'  _col(73) "{c |}" /*
		*/ _col(78) %05.3g `=mtab[1,6]' 
	di as text  abbrev("`gr1'",13)   _col(20) "{c |}" /*
		*/ as result _col(21) %8.0g `=mtab[2,1]' _col(30) "("%4.1f `=mtab[2,2]*100' "%)" _col(41) "{c |}" /*
		*/ _col(43) %8.0g `=mtab[2,3]' _col(52) "("%4.1f `=mtab[2,4]*100' "%)"  _col(62) "{c |}" /*
		*/ _col(64) %8.0g `=mtab[2,5]'  _col(73) "{c |}" /*
		*/ _col(78) %05.3g `=mtab[2,6]' 
	di as text "{hline 19}{c +}{hline 20}{c +}{hline 20}{c +}{hline 10}{c +}{hline 15}"
	di as text _newline

	
	
	*Estimating gaps
	if "`normalize'" != "" {
		qui sum `outcome' if `gr' == 0 
		qui replace `y' = `outcome' / `r(mean)' if !mi(`gr')
	}
	qui nopo_gaps `y' `gr' `prefix'_matched `prefix'_weights
	matrix colnames b = "D (raw gap)" "D0 (unexpl.)" "DX (expl.)" "DA (unmatch. A)" "DB (unmatch. B)" 

	
	* Returns
	ereturn post b, depname(`outcome') esamp(`touse') 
	ereturn matrix match_table = mtab
	ereturn local ref = "`ref0'"
	ereturn local match_set = strltrim("`match_set'")
	ereturn local prefix = "`prefix'"
	ereturn scalar n_strata = nstrata
	ereturn scalar match_strata = mstrata
	ereturn local by = "`by'"
	*Result table of gap-estimation
	di as text "Estimates for Nopo decomposition:"
	ereturn display, noomitted
	
end

/****/
/* Matching and Strata definition */
/****/
cap program drop nopo_match
program define nopo_match
	args y by _strata _matched _weights
	
	/* generating matching indicator */
	tempvar _min _max
	bys `_strata': egen `_min' = min(`by')
	bys `_strata': egen `_max' = max(`by')  
	replace `_matched' = `_min' == 0 & `_max' == 1
	
	levelsof `_strata'
	scalar nstrata = `r(r)'
	levelsof `_strata' if `_matched' == 1
	scalar mstrata = `r(r)'
	
	/* Matching Table */	
	tab `_strata' `by', matcell(strtable)
	mata: strtable = st_matrix("strtable")
	mata: mtab = J(2,6,.)
	mata: mtab[,5] = colsum(strtable)'
	mata: mtab[,1] = colsum(select(strtable,rowmin(strtable) :>0))'
	mata: mtab[,2] = mtab[,1] :/ mtab[,5]
	mata: mtab[,3] = colsum(select(strtable,rowmin(strtable) :==0))'
	mata: mtab[,4] = mtab[,3] :/ mtab[,5]
	mata: st_matrix("mtab", mtab)
	*Raw means of outcome for table
	qui sum `y' if `by' == 0 
	matrix mtab[1,6] = `r(mean)'
	qui sum `y' if `by' == 1
	matrix  mtab[2,6] = `r(mean)'
	
	/* Estimating weights */
	replace `_weights' = `_matched' 
	tempvar _celln _rown
	bys `_strata' `by': gen `_celln' = _N 
	bys `_strata': gen `_rown' = _N 
	replace `_weights' = (`_rown' - `_celln') / `_celln' if `_matched' == 1 & `by' == 1
	

end

		
/***/
/* Estimate components */
/***/

cap program drop nopo_gaps
program define nopo_gaps
	args y by _matched _weights
	
	matrix b = J(1,5,.)
	
	*D
	matrix b[1,1] = mtab[2,6] - mtab[1,6]
	
	*D0
	qui su `y' [iw=`_weights'] if `by'== 0 & `_matched' == 1
	local y0 = r(mean) 
	qui su `y' [iw=`_weights'] if `by'== 1 & `_matched' == 1
	matrix b[1,2] = `r(mean)' - `y0'

	*DX
	qui su `y' [iw=`_weights'] if `by'== 1 & `_matched' == 1
	local y0 = r(mean)
	qui su `y'  if `by'== 1 & `_matched'== 1
	matrix b[1,3] = `r(mean)' - `y0'
	
	*DA
	qui su `y'  if `by'== 0 & `_matched' == 0
	local y0=r(mean)
	qui su `y'  if `by'== 0 & `_matched' == 1
	matrix b[1,4] = cond((`r(mean)'-`y0')* mtab[1,3] / mtab[1,5] != ., (`r(mean)'-`y0')* mtab[1,3] / mtab[1,5], 0)

	*DB
	qui su `y' if `by'== 1 & `_matched' == 1
	local y0=r(mean)
	qui su `y'  if `by'== 1 & `_matched' == 0
	matrix b[1,5] =  cond((`r(mean)'-`y0')* mtab[2,3] / mtab[2,5] != ., (`r(mean)'-`y0')* mtab[2,3] / mtab[2,5], 0)
	

end
