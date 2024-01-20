
// Preparing example data for help-file examples
cap program drop nopo_ex
program define nopo_ex
args ex

if ("`ex'" == "1") {
  
  // use cattaneo data
  webuse cattaneo2, clear

  // categorize and label
  foreach v in mage fage {
    recode `v' (min/18 = 1 "-18") (19/28 = 2 "19-28") (29/38 = 3 "29-38") (39/max = 4 "39-") ///
      , gen(`v'_c)
  }
  lab var mage_c "Mother's age"
  lab var fage_c "Father's age"

}

end