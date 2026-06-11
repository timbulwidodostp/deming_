*! version 1.0.2  06mar2008
program deming, eclass
	version 9.0
        if replay() {
                if "`e(cmd)'" != "deming" {
                        error 301
                }
		syntax [, level(cilevel) ]
		if `"`e(duplicates)'"' != "" {
			DiDupNotes `"`e(expr1)'"' `"`e(expr2)'"' first
			local varnames _depvar _indepvar
	
		}
		else {
			local varnames `e(varnames)'
		}
                DiHeader e(N) e(meany) e(meanx) e(sdy) e(sdx) e(ratio) ///
			 e(rmse) `varnames' e(ErrVary) e(ErrVarx)
                eret display, level(`level')
		if `"`e(duplicates)'"' != "" {
			DiDupNotes
		}
		else if e(ErrCVy)!=. | e(ErrCVx)!=. {
			DiCVNote `e(varnames)' e(ErrCVy) e(ErrCVx)
		}
		exit
        }
	syntax varlist(min=2 max=2 numeric) [if] [in] [, ///
			 cv1(string) cv2(string) LAMBda(string) ///
			 var1(string) var2(string) 		///
			 DUPlicates(varlist numeric min=2 max=2) * ]
	marksample touse
	qui count if `touse'
	local N = r(N)
	tokenize `varlist'
	local depvar `1'
	local indepvar `2'
	local errvarx = .
	local errvary = .
	local errcvy = .
	local errcvx = .
	local parms `lambda'`cv1'`cv2'`var1'`var2'
	if `"`duplicates'"' != "" {
		if `"`parms'"' != "" {
di as err "duplicates() may not be combined with other options"
exit 198
		}
		tempvar diff2y diff2x
		tokenize `duplicates'
		local expr1 _depvar=(`depvar'+`1')/2
		local expr2 _indepvar=(`indepvar'+`2')/2
		DiDupNotes `"`expr1'"' `"`expr2'"' first
		// obtain error variances from duplicate
		qui gen double `diff2y' = (`depvar' - `1')^2 if `touse'
		qui  summ `diff2y' if `touse', meanonly
		local errvary = r(sum)/(2*`N')
		qui gen double `diff2x' = (`indepvar' - `2')^2 if `touse'
		qui  summ `diff2x' if `touse', meanonly
		local errvarx = r(sum)/(2*`N')
		local ratio = `errvarx'/`errvary'
		// create means of duplicates
		qui gen double _depvar   = (`depvar'+`1')/2 if `touse'
		qui gen double _indepvar = (`indepvar'+`2')/2 if `touse'
		local depvar _depvar
		local indepvar _indepvar
	}
	qui summ `depvar' if `touse'
	local meany = r(mean)
	local sdy = r(sd)
	qui summ `indepvar' if `touse'
	local meanx = r(mean)
	local sdx = r(sd)
	if `"`duplicates'"' != "" {
		local errcvy = 100*sqrt(`errvary')/`meany'
		local errcvx = 100*sqrt(`errvarx')/`meanx'
	}
	local parms1 `cv1' `var1'
	local parms2 `cv2' `var2'
	if `"`lambda'"' != "" {
		if `"`parms1'"' != "" | `"`parms2'"' != "" {
di as err "lambda() may not be combined with cv1(), cv2(), var1(), or var2()"
exit 198
		}
		cap confirm number `lambda'
		if _rc {
			di as err "lambda() must be positive number"
			exit 198
		}
		cap assert `lambda'>0
		if _rc {
			di as err "lambda() must be positive number"
			exit 198
		}
		local ratio = `lambda'
	}
	else {
		if `"`parms1'`parms2'"' == "" & `"`duplicates'"' == "" {
			local lambda = 1
			local ratio = 1
		}
	}
	if `"`parms1'`parms2'"' != "" {
		local nword: word count `parms1'
		if `nword' != 1 {
			if `nword' == 0 {
di as err "cv1() or var1() must also be specified"
exit 198
			}
			if `nword' == 2 {
di as err "cv1() and var1() may not be combined"
exit 198
			}
		}
		local nword: word count `parms2'
		if `nword' != 1 {
			if `nword' == 0 {
di as err "cv2() or var2() must also be specified"
exit 198
			}
			if `nword' == 2 {
di as err "cv2() and var2() may not be combined"
exit 198
			}
		}
		if `"`cv1'"' != "" {
			cap confirm number `cv1'
			if _rc {
di as err "cv1() must be expressed as percentages between 0 and 100"
exit 198
			}
			if `cv1'<= 0 | `cv1'>=100 {
di as err "cv1() must be expressed as percentages between 0 and 100"
exit 198
			}
			local errcvy = `cv1'
			local errvary = (`cv1'*`meany'/100)^2
		}
		if `"`cv2'"' != "" {
			cap confirm number `cv2'
			if _rc {
di as err "cv2() must be expressed as percentages between 0 and 100"
exit 198
			}
			if `cv2'<= 0 | `cv2'>=100 {
di as err "cv2() must be expressed as percentages between 0 and 100"
exit 198
			}
			local errcvx = `cv2'
			local errvarx = (`cv2'*`meanx'/100)^2
		}
		if `"`var1'"' != "" {
			cap confirm number `var1'
			if _rc {
				di as err "var1() must be positive number"
				exit 198
			}
			if `var1'<= 0 {
				di as err "var1() must be positive number"
				exit 198
			}
			local errvary = `var1'
		}
		if `"`var2'"' != "" {
			cap confirm number `var2'
			if _rc {
				di as err "var2() must be positive number"
				exit 198
			}
			if `var2'<= 0 {
				di as err "var2() must be positive number"
				exit 198
			}
			local errvarx = `var2'
		}
		local ratio = `errvarx'/`errvary'
	}
	local cmd _deming `depvar' `indepvar' `if' `in', ///
				lambda(`lambda') cv1(`cv1') cv2(`cv2') ///
				var1(`var1') var2(`var2') ///
				diffvars(`diff2y' `diff2x')
	local options rclass noh nodots title(Deming regression) `options'
	// run deming
	qui `cmd'
	local rmse = r(Syx)
	// display header
	DiHeader `N' `meany' `meanx' `sdy' `sdx' `ratio' `rmse' ///
		 `depvar' `indepvar' `errvary' `errvarx'
	// obtain standard errors
	tempname b V
	jackknife `indepvar' = r(slope) _cons = r(intercept), `options' : `cmd'
	mat `b' = e(b)
	mat `V' = e(V)
	// display note
	if `"`cv1'`cv2'"' != "" {
		DiCVNote `depvar' `indepvar' `errcvy' `errcvx'
	}
	if `"`duplicates'"' != "" {
		DiDupNotes
	}
	// repost results
	ereturn post `b' `V', obs(`e(N)') esample(`touse') dof(`e(df_r)')
	eret scalar df_r = e(df_r)
	eret scalar ratio = `ratio'
	eret scalar rmse = `rmse'
	eret scalar meany = `meany'
	eret scalar meanx = `meanx'
	eret scalar sdy = `sdy'
	eret scalar sdx = `sdx'
	eret scalar ErrVary = `errvary'
	eret scalar ErrVarx = `errvarx'
	eret scalar ErrCVy = `errcvy'
	eret scalar ErrCVx = `errcvx'
	eret local vce "jackknife"
	eret local cmd "deming"
	eret local varnames "`varlist'"
	eret local duplicates `"`duplicates'"'
	eret local expr1 `"`expr1'"'
	eret local expr2 `"`expr2'"'
	cap drop _depvar
	cap drop _indepvar
end

program DiDupNotes
	args expr1 expr2 first
	if `"`first'"' != "" {
		di as txt "{p 0 6 `=max(0,c(linesize)-78)'}Note: " ///
		   "average of the duplicates is " ///
		   `"used in the computation: `expr1' and `expr2'.{p_end}"'
	}
	else {
		di as txt "Note: measurement error variances are " ///
			"computed using duplicates."
	}
end

program DiCVNote
	args meth1 meth2 cv1 cv2
	di as txt "{p 0 6 `=max(0,c(linesize)-78)'}Note:" _c
	if `cv1' != . {
		di as txt " CV of " as res %3.2f `cv1' ///
		   as txt "% is used to compute error variance about " ///
			  "mean of `meth1'" _c
		local punct ;	
	}
	if `cv2' != . {
		di as txt "`punct' CV of " as res %3.2f `cv2' ///
		   as txt "% is used to compute error variance about " ///
			  "mean of `meth2'.{p_end}"
	}
	else {
		di as txt ".{p_end}"
	}
	
end

program DiHeader
	args N m1 m2 sd1 sd2 lambda rmse meth1 meth2 verr1 verr2
	if `verr1' != . {
		local dierrti _col(37) "Err. Var."
		local len = 23+8
		local rghtcol = 53
	}
	else {
		local dierr qui
		local len = 19
		local rghtcol = 49
	}
	di
	di as txt "Deming regression" _c
	di _col(`rghtcol') as txt "Number of obs" _col(68) "=" _col(70) ///
	   as res %9.0g `N'
	di
	di as text _col(14)"{c |}" _col(19) "Mean" ///
		   _col(26) "Std. Dev." `dierrti'
	di as text "{hline 13}{c +}{hline `len'}"
	di as text %12s abbrev("`meth1'",12) " {c |}" ///
		as result " " %7.0g `m1' "    " %7.0g `sd1' _c 
	`dierr' di as res "   " %9.0g `verr1' _c
	di _col(`rghtcol') as txt "Variance ratio" _col(68) "=" ///
	   _col(70) as res %9.0g `lambda'
	di as text %12s abbrev("`meth2'",12) " {c |}" ///
		as result " " %7.0g `m2' "    " %7.0g `sd2' _c
	`dierr' di as res "   " %9.0g `verr2' _c
	di _col(`rghtcol') as txt "Root MSE" _col(68) "=" _col(70) as res ///
        	        %9.0g `rmse'
	di
end

