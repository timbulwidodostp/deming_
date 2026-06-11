*! version 1.0.2  06mar2008
program _deming, rclass
        syntax varlist(min=2 max=2 numeric) [if] [in] [, ///
			cv1(string) cv2(string) LAMBda(string) ///
			var1(string) var2(string) 		///
			diffvars(varlist numeric min=2 max=2) ]       
	tokenize `varlist'
        // meth2 is x
        local meth1 `1'
        local meth2 `2'
        marksample touse
        qui count if `touse'
        local N = r(N)
	local eVar1 = .
	local eVar2 = .
	if `"`diffvars'"' != "" {
		// if duplicates are used
		tokenize `diffvars'
		qui  summ `1' if `touse', meanonly
		local eVar1 = r(sum)/(2*`N')
		qui  summ `2' if `touse', meanonly
		local eVar2 = r(sum)/(2*`N')
		local lambda = `eVar2'/`eVar1'
	}
	qui summ `meth1' if `touse'
	local mean1 = r(mean)
	local S21 = r(Var)
	qui summ `meth2' if `touse'
	local mean2 = r(mean)
	local S22 = r(Var)
	qui correlate `meth1' `meth2' if `touse', cov
	local cov = r(cov_12)
	if `"`lambda'"' == "" {
		if "`cv1'" == "" {
			local eVar1 = `var1'
		}
		else {
			local eVar1 = (`cv1'*`mean1'/100)^2
		}
		if "`cv2'" == "" {
			local eVar2 = `var2'
		}
		else {
			local eVar2 = (`cv2'*`mean2'/100)^2
		}
		local lambda = `eVar2'/`eVar1'
	}
	// compute Deming's intercept and slope
	local U = (`S21'-(1/`lambda')*`S22')/(2*`cov')
        local b = `U' + sign(`cov')*sqrt((`U')^2+1/`lambda')
        local a = `mean1'-`b'*`mean2'
	local Syx = sqrt((`N'-1)*(`S21'-(`b'*`cov'))/(`N'-2))        
        ret scalar ratio = `lambda'
        ret scalar N = `N'
        ret scalar slope = `b'
        ret scalar intercept = `a'
        ret scalar Syx = `Syx'
        ret scalar meany = `mean1'
        ret scalar meanx = `mean2'
	ret scalar ErrVary = `eVar1'
        ret scalar ErrVarx = `eVar2'
end
