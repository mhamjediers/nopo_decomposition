
clear


cap program drop sims
program define sims
	set obs 100
	gen double x = .
	gen double w = .
	gen double xw = .
	gen double di = .
	count 
	qui forval i = 1/`r(N)' {
		preserve 
			clear
			set obs 100
			gen w = runiform(0,5)
			gen x = rnormal(1,1)
			sum x 
			local x = `r(mean)'
			sum x [iw = w]
			local xw = `r(mean)'
			sum w
			local w = `r(mean)'
		restore 
		replace x = `x' in `i'
		replace w = `w' in `i'
		replace xw = `xw' in `i'
		replace di = `x' - `xw' in `i'
	}
	noisily: corr x xw, cov
	scalar cov = `r(cov_12)'
	for any x w xw di: gen double se_X = X
	collapse (mean) x w xw di (sd) se_*
end

clear
set obs 100
sims

gen V = sqrt(se_x^2 + se_xw^2 - 2*cov)

***--> covariance wieder über die Strata berechnen lassen, wie auch cov von w und y


