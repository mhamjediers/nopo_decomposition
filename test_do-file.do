
adopath ++ "C:\Users\Maik\Dropbox\Scripts\do-files\nopo-decomposition-ado"

use "http://fmwww.bc.edu/RePEc/bocode/o/oaxaca.dta", clear

nopodecomp lnwage isco married kids6, by(female)