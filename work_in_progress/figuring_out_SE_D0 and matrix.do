

clear 
set obs 100 // number of repitions

gen N = .
gen double true1 = .
gen double true2 = .
gen double direct1 = .
gen double direct2 = .
gen double matrix1 = .
gen double matrix2 = .
gen double matrix3 = .
gen double matrix4 = .
gen double nopo = .

count
foreach s of num 1 (1) `r(N)' {

	preserve 
		clear 
		set obs 20
		gen w = runiform()
		gen x = rnormal(3,1.5)  // + 1*w
		*reg x w // enforce no covariance --> if no covariance, all estimations bring out the same result
		*predict res, res
		*replace x = res + _b[_cons]

		// matrix equation
		mata: x = st_data(.,"x",.)
		mata: w = st_data(.,"w",.)
		mata: matrix1 = (sum(x'*variance(w)*x + w'*variance(x)*w) ///
			+ trace(((w :-mean(w))*(w :-mean(w))') * ((x :-mean(x))*(x :-mean(x))'))) :/ rows(x)^2
		mata: matrix2 = (sum(x*variance(w)*x' + w*variance(x)*w') ///
			+ 2 * trace(((w :-mean(w))*(w :-mean(w))') * ((x :-mean(x))*(x :-mean(x))'))) :/ rows(x)^2
		mata: matrix3 = (sum(x*variance(w)*x' + w*variance(x)*w') ///
			+ 2 * mean(w) * mean(x) * trace(((w :-mean(w))*(w :-mean(w))') * ((x :-mean(x))*(x :-mean(x))'))) :/ rows(x)^2
	
		mata: matrix4 = (sum(x*variance(w)*x' + w*variance(x)*w') ///		
			+ trace(variance(w) * variance(x))) :/ rows(x)^2 				// differs only due to degrees of freedom (var automatically substracts 1)
		
		*mata: matrix4 = (sum(x'*variance(w*w')*x + w'*variance(x*x')*w)  /// completely different result
		*	+ trace(variance(w*w') * variance(x*x'))) :/ rows(x)^2
	
		mata: st_numscalar("matrix1", matrix1)
		mata: st_numscalar("matrix2", matrix2) // differs from 1 by switching transpose in first part and adding 2*trace  			--> works for cases without variance
		mata: st_numscalar("matrix3", matrix3) // differs from 1 by switching transpose in first part and adding 2*mean*mean*trace  --> works for cases without variance
		mata: st_numscalar("matrix4", matrix4) // differs from 1 changing the variance function to w*w'

		*Reference to footnote 5 in Jann (2008) and why the trace is vanashing: both are of order (n^-1) and (n^-2); yet in our case 
		*they have the order of the strata (or strata per n?), which means they are not only not vanashing, but also depending on N in strata
		
		*General idea of Nopo: We have an empirical distribution of the counterfactuals across strata; 
		*	and estimate the variance-covariance matrix of y and w across strata; yet this holds only assymptotically and then see above again ;)
		
		*yet, it also holds when the distributions are not highly unequal with respect to y, because then no covariance between y and w; 
		*	--> unlikely in the case of decomposition, if we are actually interested in it
		
		local matrix1 = `=matrix1' 
		local matrix2 = `=matrix2' 
		local matrix3 = `=matrix3'
		local matrix4 = `=matrix4'
		
		// direct implementation of delta method
		qui {
			sum x
			local xm = `r(mean)'
			local xv = `r(Var)'
			sum w
			local wm = `r(mean)'
			local wv = `r(Var)'
		}
		cap noisily corr x w
		cap noisily corr x w, cov 
		
		if _rc == 0 {
			local direct1 = (`xm'^2 *`wv' + `wm'^2 *`xv' + `r(cov_12)') 
			local direct2 = (`xm'^2 *`wv' + `wm'^2 *`xv' + 2 * `xm' * `wm' * `r(cov_12)')
		}
		// nopo ( is exactly equal to matrix2-equation)
		mata: x = st_data(.,"x",.)
		mata: w = st_data(.,"w",.)
		mata: comp1 = (sum(x*variance(w)*x' + w*variance(x)*w')) :/ rows(x)^2
		mata: comp2 = (sum(((w :-mean(w))*(w :-mean(w))') * ((x :-mean(x))*(x :-mean(x))')) ///
			- trace(((w :-mean(w))*(w :-mean(w))') * ((x :-mean(x))*(x :-mean(x))'))) :/ rows(x)^2
		mata: nopo = comp1 - 2*comp2  // comp2 gets 0 if no covariance occurs
		mata: st_numscalar("nopo", nopo)
		local nopo = `=nopo' 
		
		// empirical SE
		local true1 = `xm' * `wm' // difference to true1 is zero if no covariance (obviously)
		gen xw = x * w 
		sum xw
		local true2 = `r(mean)'
		local N = `r(N)'
		
		
	restore 

	replace N =    `N'  in `s'
	replace true1 =    `true1'  in `s'
	replace true2 =    `true2'  in `s'
	if ("`direct1'" != "") replace direct1 = `direct1' in `s'
	if ("`direct1'" != "") replace direct2 = `direct2' in `s'
	replace matrix1 = `matrix1' in `s'
	replace matrix2 = `matrix2' in `s'
	replace matrix3 = `matrix3' in `s'
	replace matrix4 = `matrix4' in `s'
	replace nopo =    `nopo'    in `s'

}

sum true1, d
sum true2, d
br



collapse direct* matrix* nopo N (sd) true1 true2
replace true1 = true1^2 * N 
replace true2 = true2^2 * N

***---> the problem in nopos equation is that too often the variance of y in a stratum is zero (because only one observation)
* thereby we substract too much of the covariance from the underestimated variance in y


/*
Var =   (P1^2) * VarP2
				+ (P2^2) * VarP1
				+ 2 * P1 * P2 * CoV 
				(last term: covariance is OK, but P1*P2 is not: 
				looks like variance estimation for weighted sum of variables, but that does not fit
				with the delta method for a product...what to do?!)