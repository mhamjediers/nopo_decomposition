
// Preparing example data for help-file examples
cap program drop nopo_ex
program define nopo_ex
args ex

if ("`ex'" == "1") {
  
  // use data from oaxaca-help-file
  use "http://fmwww.bc.edu/RePEc/bocode/o/oaxaca.dta", clear



  // categorize and label
	for any exper tenure: gen c_X = round(X,5)
	gen c_educ = round(educ,1)

	lab var c_educ "years of educational attainment (rounded)"
	lab var c_exper "years of work experience (5-year intervalls)"
	lab var c_tenure "years of job tenure (5-year intervalls)"
	
	lab def female 0 "Men" 1 "Women"
    lab val female female
}

end