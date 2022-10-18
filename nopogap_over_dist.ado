// Program to produce plots for gaps across wage distribution
cap program drop nopogap_over_dist
program define nopogap_over_dist
syntax varname [if] [aweight], ///
	BY(varlist) ///
	[Quantiles(integer 100)] ///
	[QMIN(integer 1)] ///
	[QMAX(integer 100)] ///
	[Stat(string)] ///
	[scaleumdiff(string)] ///
	[revsign] /// reverse sign
	[TWtype(string)] ///
	[RELative(integer 1)] ///
	[SAVE(string asis)] ///
	[USE(string asis)] *

quietly {
	preserve
		if "`if'" != "" {
			keep `if'
		}
		levelsof `by', local(levels) // get levels
		if `r(r)' != 2 {
			dis in red "by-variable has more than two values"
			assert `r(r)' == 2
		}
		// allow for quantile mean (or any other stat) estimation'
		/*
		 Correction: xtile aggregates quantiles if they contain constant values.
		 Fill up by hand to avoid empty cells.
		*/
		gen p = .
		gen r = runiform()
		foreach l in `levels' {
		    xtile p_`l' = `varlist' if `by'==`l' [`weight'`exp'], nquantiles(`quantiles')
		    replace p = p_`l' if `by'==`l'
		}
		// collapse, use sum of weights (passed via `exp') as N
		collapse (`stat') `stat'_wage=`varlist' (sum) n=`=(subinstr("`exp'","=","",.))' [`weight'`exp'] if !mi(`varlist'), by(`by' p)
		reshape wide `stat'_wage n, i(p) j(`by')
		foreach l in `levels' {
			replace `stat'_wage`l' = `stat'_wage`l'[_n-1] if mi(`stat'_wage`l')
			// check for a maximum of 3 merged ranks
			replace n`l' = n`l'/2 if mi(n`l'[_n+1]) & !mi(p[_n+1]) & !mi(n`l'[_n+2])
			replace n`l' = n`l'/3 if mi(n`l'[_n+1]) & mi(n`l'[_n+2]) & !mi(p[_n+1]) & !mi(p[_n+2])
			replace n`l' = n`l'[_n-1] if mi(n`l') & !mi(n`l'[_n-1])
		}
		// scale to percentiles
		replace p = round((p/`quantiles')*100,1)
		rename p rank
		lab var rank "Percentile of group-specific wage distribution"

		// gen diff (scale if requested for DT/DC)
		tokenize `levels'
		if (`relative'	!= 1) {
			if ("`scaleumdiff'"=="by2") gen diff = (`stat'_wage`2' - `stat'_wage`1') * (n`2'/(n`1'+n`2'))
				else if ("`scaleumdiff'"=="by1") gen diff = (`stat'_wage`2' - `stat'_wage`1') * (n`1'/(n`1'+n`2'))
				else gen diff = `stat'_wage`2' - `stat'_wage`1'
			lab var diff "Diff. in `stat' values per quantile"
		}
		else {
			if ("`scaleumdiff'"=="by2") gen diff = (1 - (`stat'_wage`1' / `stat'_wage`2')) * (n`2'/(n`1'+n`2'))
				else if ("`scaleumdiff'"=="by1") gen diff = (1 - (`stat'_wage`1' / `stat'_wage`2')) * (n`1'/(n`1'+n`2'))
				else gen diff = (1 - (`stat'_wage`1' / `stat'_wage`2'))
			replace diff = diff * 100 // percent
			lab var diff "Relative diff. in `stat' values per quantile"
		}

		if ("`revsign'" != "") replace diff = -diff

		if "`twtype'" == "" {
			local twtype line
		}

		if `"`use'"' == "" {
			twoway `twtype' diff rank, `options'
		}
		else {
			describe using `use'
			if "`r(datalabel)'" != "generated by tw_gap_over_dist" {
				dis in red `"Dataset in "use" was not generated by tw_gap_over_dist"'
				assert "`r(datalabel)'" == "generated by tw_gap_over_dist"
			}
			append using `use', gen(prev)
			twoway `twtype' diff rank if prev == 1 & inrange(rank,`qmin',`qmax') ///
				|| `twtype' diff rank if prev == 0 & inrange(rank,`qmin',`qmax'), `options'
		}

		if `"`save'"' != "" {
			drop if !inrange(rank,`qmin',`qmax')
			label data "generated by tw_gap_over_dist"
			save `save'
		}

	restore
}

end
