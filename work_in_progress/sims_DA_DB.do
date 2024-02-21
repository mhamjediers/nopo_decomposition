

clear

cap program drop sims
program define sims
	set obs 100
	gen double xn = .
	gen double xu = .
	gen double mu = .
	gen double di = .
	gen double sh = .
	gen double D = .
	gen double V = .
	gen double IF_xn = .
	gen double IF_xu = .
	gen double IF_di = .
	gen double IF_sh = .
	gen double IF_D  = .
	
	count 
	qui forval i = 1/`r(N)' {
		preserve 
			clear
			set obs 1000
			gen u = runiform() < `1'
			gen y = rnormal(1,1) + `2'*u
			sum u,
			local share = `r(mean)'
			ttest y, by(u)
			
			// MATA IF - ESTIMATION
				mata: N = st_nobs()
				mata: D = st_data(.,"u")
				mata: Y = st_data(.,"y")
				
				// compute mean of matched / unmatched
				mata: mA = mean(Y, D)
				mata: IF_m = N/sum(D) * D :* (Y :- mA)
				mata: uA = mean(Y, !D)
				mata: IF_u = N/sum(!D) * !D :* (Y :- uA)
				// compute mean difference
				mata: gap = mA - uA
				mata: IF_gap = IF_m - IF_u 
				
				// compute share of matched
				mata: share = sum(!D) / N
				mata: IF_share = (!D :- share)
				
				// compute weighted mean-diff (this IF does not work yet)
				// (but we also do not have to figure it out anymore, when treating the share of matched as exogeneous)
				mata: DA = mA * share - uA * share
				mata: IF_DA = N/sum(D) * D :* (Y :- mA) :* (!D :- share) - N/sum(!D) * !D :* (Y :- uA) :* (!D :- share)
				mata: (mA, uA, gap, share, DA)', sqrt(diagonal(variance((IF_m, IF_u, IF_gap, IF_share, IF_DA)) / N))
				
				mata: IF_m = diagonal(variance((IF_m, IF_u, IF_gap, IF_share, IF_DA)) / N)
				mata: st_matrix("IF_m", IF_m)
				
		restore 
		replace xn = `r(mu_1)' in `i'
		replace xu = `r(mu_2)' in `i'
		replace mu = xn*(1-`share') + xu*`share'
		replace di = xn - xu in `i'
		replace sh = `share' in `i'
		replace D = di * `share' in `i'
		
		replace V = di^2*(`share'*(1-`share') / 999) ///
			+ `share'^2 * `r(se)'^2 + ///
			(`share'*(1-`share') / 999)*`r(se)'^2 ///
			in `i'
			
		replace IF_xu = IF_m[1,1] in `i'
		replace IF_xn = IF_m[2,1] in `i'
		replace IF_di = IF_m[3,1] in `i'
		replace IF_sh = IF_m[4,1] in `i'
		replace IF_D  = IF_m[5,1] in `i'
		
	}
	corr di sh, cov
	gen C = `r(cov_12)'
	for any xn xu mu di sh D: gen double se_X = X
	collapse (mean) s d xn xu mu di sh D V C IF_* (sd) se_*
end


set obs 1 
gen s = .
gen d = .


tempfile sims
save `sims', replace

qui foreach s of num  0.01 0.05 0.1 0.2 0.3 { // 
	foreach d of num  0.1 0.5 3 5 7 {  // 
		dis in red  `s' " and " `d'
		clear
		set obs 1
		gen s = `s'
		gen d = `d'
		sims `s' `d'
		append using `sims'
		save `sims', replace
	}
}

replace V = sqrt(V)
foreach v of var IF_* {
	replace `v' = sqrt(`v')
}



gen t = sqrt(V^2 - (2*C))
corr t se_D V
*---> no correction for correlation of sh and di needed 
