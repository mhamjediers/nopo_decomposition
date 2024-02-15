

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
		
	}
	corr di sh, cov
	gen C = `r(cov_12)'
	for any xn xu mu di sh D: gen double se_X = X
	collapse (mean) s d xn xu mu di sh D V C (sd) se_*
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


gen t = sqrt(V^2 - (2*C))
corr t se_D V
*---> no correction for correlation of sh and di needed 
