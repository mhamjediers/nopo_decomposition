
*********
*
*********
cap program drop nopodecomp
program define nopodecomp , eclass
	syntax varlist , BY(varname) [PREFix(string) REPLACE swap NORMalize BOOTSTRAP BSOPTS(string)]
	

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
		if "`swap'" != "" { // reverse order for reference group
			local levels = strreverse("`levels'")
			local swap = 1
		}
		else {
			local swap = 0
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
	qui match_table `y' `gr' `prefix'_strata `prefix'_matched
	di as text _newline
	di as text "Exact matching for Nopo decomposition:"	_col(60) "N(strata) = " 		_col(85) nstrata
	di as text " " 										_col(60) "N(matched strata) = " _col(85) mstrata
	di as text _newline

	di as text " Group " _col(22) "{c |}" /*
		*/ _col(29) "Matched" _col(40) "{c |}" /*
		*/ _col(45) "Unmatched"  _col(58) "{c |}" /*
		*/ _col(62) "Total"  _col(71) "{c |}" /*
		*/ _col(73) "Mean(" abbrev("`outcome'",9) ")"
	di as text "{hline 21}{c +}{hline 17}{c +}{hline 17}{c +}{hline 12}{c +}{hline 15}"
	di as text "A: " abbrev("`gr0'",13) " (ref)"  _col(22) "{c |}" /*
		*/ as result _col(22) %8.0g `=mtab[1,1]' _col(32) "("%4.1f `=mtab[1,2]*100' "%)" _col(40) "{c |}" /*
		*/ _col(41) %8.0g `=mtab[1,3]' _col(50) "("%4.1f `=mtab[1,4]*100' "%)"  _col(58) "{c |}" /*
		*/ _col(59) %8.0g `=mtab[1,5]'  _col(71) "{c |}" /*
		*/ _col(75) %05.3g `=mtab[1,6]' 
	di as text "B: "  abbrev("`gr1'",13)   _col(22) "{c |}" /*
		*/ as result _col(22) %8.0g `=mtab[2,1]' _col(32) "("%4.1f `=mtab[2,2]*100' "%)" _col(40) "{c |}" /*
		*/ _col(41) %8.0g `=mtab[2,3]' _col(50) "("%4.1f `=mtab[2,4]*100' "%)"  _col(58) "{c |}" /*
		*/ _col(59) %8.0g `=mtab[2,5]'  _col(71) "{c |}" /*
		*/ _col(75) %05.3g `=mtab[2,6]' 
	di as text "{hline 21}{c +}{hline 17}{c +}{hline 17}{c +}{hline 12}{c +}{hline 15}"

	
	*Estimating gaps
	if "`normalize'" != "" {
		qui sum `outcome' if `gr' == 0 
		qui replace `y' = `outcome' / `r(mean)' if !mi(`gr')
	}
	
	if "`bootstrap'" == "" {
		nopo_gaps  `y' `gr' `prefix'_strata  `prefix'_matched `prefix'_weights
		matrix b = e(b)
		matrix colnames b = "raw gap (D)" "unexpl. (D0)" "explain. (DX)" "unmatch. A (DA)" "unmatch. B (DB)" 
		ereturn repost b = b, rename
	}	
	else {
		di as text _newline
		*with repeatedly generating new weights or not?? (if yes, just add the repeat_weights at the end)
		*with bootstrap-specific weights, the SE get much smaller (I expected the opposite)
		bootstrap, noheader nolegend nowarn notable `bsopts': ///
			nopo_gaps `y' `gr' `prefix'_strata  `prefix'_matched `prefix'_weights 
		matrix b = e(b)
		matrix V = e(V)
		matrix colnames b = "raw gap (D)" "unexpl. (D0)" "explain. (DX)" "unmatch. A (DA)" "unmatch. B (DB)" 
		matrix colnames V = "raw gap (D)" "unexpl. (D0)" "explain. (DX)" "unmatch. A (DA)" "unmatch. B (DB)" 
		matrix rownames V = "raw gap (D)" "unexpl. (D0)" "explain. (DX)" "unmatch. A (DA)" "unmatch. B (DB)" 
		ereturn repost b = b V = V, rename
	}
	
	* Returns
	ereturn matrix match_table = mtab
	ereturn local ref = "`ref0'"
	ereturn local match_set = strltrim("`match_set'")
	ereturn local prefix = "`prefix'"
	ereturn scalar n_strata = nstrata
	ereturn scalar match_strata = mstrata
	ereturn local by = "`by'"
	ereturn scalar swap = `swap'
	ereturn local aweight "XXX" // for passthru; only relevant if we implement (distplot already prepped)
	ereturn local cmd "nopodecomp" // for postestimation checking
	*Result table of gap-estimation
	di as text _newline
	di as text "Estimates for Nopo decomposition:"
	ereturn display, noomitted
	
end

/****/
/* Matching and weights */
/****/
cap program drop nopo_match
program define nopo_match
	args y by _strata _matched _weights
	
	/* generating matching indicator */
	tempvar _min _max
	bys `_strata': egen `_min' = min(`by')
	bys `_strata': egen `_max' = max(`by')  
	replace `_matched' = `_min' == 0 & `_max' == 1
					
	/* Estimating weights */
	replace `_weights' = `_matched' 
	tempvar _celln _rown
	bys `_strata' `by': gen `_celln' = _N 
	bys `_strata': gen `_rown' = _N 
	replace `_weights' = (`_rown' - `_celln') / `_celln' if `_matched' == 1 & `by' == 1


end


/***/
/* Matching table */
/***/
cap program drop match_table
program define match_table 
	args y by _strata _matched 

		mata: st_numscalar("nstrata", colmax(st_data(., "`_strata'")))
		mata: st_numscalar("mstrata", length(uniqrows(st_data(., "`_strata'","`_matched'"))))

		*does not work with a high number of strata
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
end
		
		
/***/
/* Estimate components */
/***/

cap program drop nopo_gaps
program define nopo_gaps, eclass
	args y by _strata _matched _weights repeat_weights
	
	
	if "`repeat_weights'" != "" {
		nopo_match `y' `by' `_strata' `_matched' `_weights'
	}

	matrix b = J(1,5,.)
	
	*D
	qui sum `y' if `by' == 0 
	local y0 = `r(mean)'
	qui sum `y' if `by' == 1
	matrix b[1,1] = `r(mean)' - `y0'

	
	*D0
	qui su `y' [iw=`_weights'] if `by'== 0 & `_matched' == 1
	local y0 = `r(mean)'
	qui su `y' [iw=`_weights'] if `by'== 1 & `_matched' == 1
	matrix b[1,2] = `r(mean)' - `y0'

	*DX
	qui su `y' [iw=`_weights'] if `by'== 1 & `_matched' == 1
	local y0 = `r(mean)'
	qui su `y'  if `by'== 1 & `_matched'== 1
	matrix b[1,3] = `r(mean)' - `y0'
	
	*DA
	qui su `y'  if `by'== 0 & `_matched' == 0
	local y0 = `r(mean)'
	qui su `y'  if `by'== 0 & `_matched' == 1
	matrix b[1,4] = cond((`r(mean)'-`y0')* mtab[1,3] / mtab[1,5] != ., (`r(mean)'-`y0')* mtab[1,3] / mtab[1,5], 0)

	*DB
	qui su `y' if `by'== 1 & `_matched' == 1
	local y0 = `r(mean)'
	qui su `y'  if `by'== 1 & `_matched' == 0
	matrix b[1,5] =  cond((`r(mean)'-`y0')* mtab[2,3] / mtab[2,5] != ., (`r(mean)'-`y0')* mtab[2,3] / mtab[2,5], 0)
	
	ereturn post b, esamp(`touse')
	end
