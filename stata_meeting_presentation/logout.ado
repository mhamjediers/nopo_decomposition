cap program drop logout
program define logout
qui {	
	gettoken left right2 : 0, parse(":")
	if `"`left'"' == ":" {
		local right `"`right2'"'
	}
	else {
		gettoken right3 right : right2, parse(":")
	}
	local 0 : copy local left
	syntax using, [CMD] [ONLYCMD]
	
	
	log `using', replace text
	local filename = "`r(filename)'"
	if "`cmd'" == "cmd" | "`onlycmd'" == "onlycmd"{
		noisily: display ".`right'"
	}
	if "`onlycmd'" != "onlycmd" {
		noisily: `right'
	}
	log close 
	if "`onlycmd'" == "onlycmd" {
		noisily: `right'
	}
}
	dis ""
	dis as text "Output saved as `filename'"
end 

