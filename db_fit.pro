function db_gauss1, x, p
  return, p[0]*exp(-((x-p[1])/(2.*p[2]))^2.) + p[3]
end
function db_gauss2, x, p
  return, p[0]*exp(-((x-p[1])/(2.*p[2]))^2.) + p[3]*exp(-((x-p[4])/(2.*p[5]))^2.) + p[6]
end
function plot_dbgauss, x, y, p
  p33 = plot(x, y)
  p34 = plot(x, db_gauss(x, p), '--r', over=p33)
  p341 = plot(x, db_gauss1(x, p[0:2]), '--b', over=p33)
  p342 = plot(x, db_gauss1(x, p[3:5]), '--g', over=p33)
end
; amp1, cen1, wid1, amp2, cen2, wid2, lev1, chisq1, amp0, cen0, wid0, lev0, chisq0
;   0     1     2     3     4     5     6     7      8     9     10    11     12
function db_fit, wp, specp
  si_cen = 1402.77d0
  w_th_si = si_cen/3d8*sqrt(8.*alog(2.)*1.38d-23*10d0^(4.9)/(28.0855*1.6605d-27))  ; in angstrom
  ;; ~ 0.053 angstrom (https://iris.lmsal.com/itn38/diagnostics.html --> 0.05)
  w_inst = 0.026 ; in angstrom
  ooe_fac = 1./(2.*sqrt(alog(2))) ; 1/e factor for w_{1/e}
  fwhm_fac = 1./(2.*sqrt(2.*alog(2))) ; FWHM * FWHM_fac = Gaussian sigma
  w_th_si = w_th_si*fwhm_fac  ; in Gaussian sigma
  w_inst = w_inst*fwhm_fac    ; in Gaussian sigma
  min_width = sqrt(w_inst^2. + w_th_si^2.)
  if n_elements(dn2phot) eq 0 then dn2phot = 4.  ; Please check using iris_get_response(time)
  if n_elements(init) eq 0 then init0 = [40, si_cen, 0.05, 40, si_cen, 0.05, 0.] $
  else init0 = init

  lims = {value:0., fixed:0, limited:[0, 0], limits:[0., 0.]}
  lims = replicate(lims, 7)
  lims[0].limited[0] = 1 & lims[0].limits[0] = 0d                        ; amplitude
  lims[1].limited[*] = 1 & lims[1].limits = si_cen+[-1, 1]*0.5         ; central wavelength
  lims[2].limited[*] = 1 & lims[2].limits = [min_width, 0.5]             ; width
  lims[3].limited[0] = 1 & lims[3].limits[0] = 0d                        ; amplitude
  lims[4].limited[*] = 1 & lims[4].limits = si_cen+[-1, 1]*0.5           ; central wavelength
  lims[5].limited[*] = 1 & lims[5].limits = [min_width, 0.5]             ; width
  sig = sqrt(specp) + 1.
  ;err[where(specp le 0.)] = 1d5
  err = 1.
  res1 = mpfitfun('db_gauss2', wp, specp, err, init0, $
    parinfo=lims, quiet=1, weights=1d, $/specp, $
    maxiter=400, status=st, ftol=1d-9, /nan)
  ;  if abs(res1[1]-si_cen) gt abs(res1[4]-si_cen) then res1 = [res1[3:5], res1[0:2], res1[6]]
  if res1[0] lt res1[3] then res1 = [res1[3:5], res1[0:2], res1[6]]

  chisq1 = total((specp - db_gauss2(wp, res1))^2./(sig)^2., /nan)

  res0 = mpfitfun('db_gauss1', wp, specp, err, [init0[0:2], 0], $
    parinfo=lims[3:6], quiet=1, weight=1d, $
    maxiter=400, status=st, ftol=1d-9, /nan)
  chisq0 = total((specp - db_gauss1(wp, res0))^2./(sig)^2., /nan)
  ;  stop
  res = [res1, chisq1, res0, chisq0]
  res[[1, 4, 9]] = (res[[1, 4, 9]] - si_cen)/si_cen*3d5 ; im km/s
  res[[2, 5, 10]]  = ((sqrt(res[[2, 5, 10]]^2. - w_inst^2. - w_th_si^2.)>0d0)*sqrt(2)) $ ; in 1/e nonthermal
    /si_cen*3d5
  ;  stop
  return, res

end


