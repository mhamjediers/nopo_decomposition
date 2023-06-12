program noheadlog
        tempname in out
	capture erase "noheadlog_file.log"
        file open `in' using "`2'", read
        quietly: file open `out' using "noheadlog_file.log", write replace
        file read `in' line
		file read `in' line
        while !r(eof) {
                local w1 : word 1 of `line'
                local head = match("`w1'", "*name:*")
                if "`head'"=="1" {
                        file read `in' line
                        file read `in' line
                        file read `in' line
                        file read `in' line
				 		file read `in' line
                }
                local close = match("`line'", "*log close*")
                if "`close'"=="1" {
                        file read `in' line2
                        local w1 : word 1 of `line2'
                        local right_arrow = match("`w1'", "*>*")
                        if "`right_arrow'"=="0" {
                                local line "`line2'"
                        }
                }
                else {
                        file write `out' `"`macval(line)'"' _n
                        file read `in' line
                }
        }
file close `in' 
file close `out'
quietly: copy "noheadlog_file.log" `"`2'"', replace
capture erase "noheadlog_file.log"
end
