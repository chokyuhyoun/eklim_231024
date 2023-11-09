function gaussian, x, p
  return, p[0]*exp(-((x-p[1])/(2.*p[2]))^2.) + p[3]
end

function gaussian_fit, wv, specp, si_cen=si_cen, w_th_si=w_th_si, w_inst=w_inst, init=init
  if n_elements(si_cen) eq 0 then si_cen = 1402.77d0
  if n_elements(w_th_si) eq 0 then w_th_si = si_cen/3d8*sqrt(8.*alog(2.)*1.38d-23*10d0^(4.9)/(28.0855*1.6605d-27))  ; in angstrom
  ;; ~ 0.053 angstrom (https://iris.lmsal.com/itn38/diagnostics.html --> 0.05)
  if n_elements(w_inst) eq 0 then w_inst = 0.026 ; in angstrom

  min_width = sqrt(w_inst^2. + w_th_si^2.)
  if n_elements(init) eq 0 then init0 = [1., si_cen, min_width, 0.] $
                           else init0 = init

  lims = {value:0., fixed:0, limited:[0, 0], limits:[0., 0.]}
  lims = replicate(lims, n_elements(init0))
  lims[0].limited[0] = 1 & lims[0].limits[0] = 0d                        ; amplitude
  lims[1].limited[*] = 1 & lims[1].limits = si_cen+[-1, 1]*0.5         ; central wavelength
  lims[2].limited[*] = 1 & lims[2].limits = [min_width, 0.5]             ; width

  ;err[where(specp le 0.)] = 1d5
;  err = 1.
  err = sqrt(specp>0.)
  res = mpfitfun('gaussian', wv, specp, err, init0, $
                  parinfo=lims, quiet=1, weights=1d, $/specp, $
                  maxiter=400, status=st, ftol=1d-9, /nan)
  return, res
end

path = '/Users/khcho/Desktop/eklim/'
cd, path
iris_path = '/irisa/data/level2/2020/07/30/20200730_155928_3600011659/'
f = file_search(iris_path+'iris*_raster*.fits')
t_ref = ['20200730_170135', $
         '20200730_170152', $
         '20200730_170209', $
         '20200730_170226']
t_ref0 = anytim2tai(t_ref)
t_ref0 = t_ref0[sort(t_ref0)]

si_cen = 1402.77d0 
;si_cen = 1393.755d0  ; select one of them 

si_cen_str = string(round(si_cen), f='(i0)')

w_th_si = si_cen/3d8*sqrt(8.*alog(2.)*1.38d-23*10d0^(4.9)/(28.0855*1.6605d-27))  ; in FWHM angstrom
w_inst = 0.026 ; in FWHM angstrom
fwhm_fac = 1./(2.*sqrt(2.*alog(2))) ; FWHM * FWHM_fac = Gaussian sigma
w_th_si0 = w_th_si*fwhm_fac  ; in Gaussian sigma
w_inst0 = w_inst*fwhm_fac    ; in Gaussian sigma

resp = iris_get_response(anytim2utc(t_ref[0], /ccsds))
getmin = min(abs(resp.lambda-si_cen/10.), imin)
pix_x_size=0.33   ; arcsec
pix_y_size0=0.16635  ; arcsec
en=1.986e-8/si_cen
;;;
sji_file = file_search(iris_path+'iris*SJI*.fits')
read_iris_l2, sji_file, sji_indices, sji_data
sji_time = anytim2tai(sji_indices.DATE_OBS)
iris_lct, 'SJI_1400', rr, gg, bb
iris1400_ct = [[rr], [gg], [bb]]
;for i=0, n_elements(f)-1 do begin
i = 3
  dd = iris_obj(f[i])
  line_id = dd->getline_id()
  si_id = (where(line_id eq 'Si IV '+ si_cen_str, /null))[0]
  sg_time = dd->ti2tai()
  exp_time = dd->getexp(iwin=si_id) ;; array
  spec_bin = (dd->binning_spectral(si_id))[0]
  spat_bin = (dd->binning_region('FUV'))[0]
  pix_y_size = pix_y_size0*spat_bin
  dn2phot_sg = resp.dn2phot_sg[0]
  pix_size=!pi/(180.)*!pi/(180.)*(pix_x_size/3600.)*(pix_y_size/3600.)
  area_sg = resp.area_sg[imin, 0]
  flux_per_dn = en*dn2phot_sg/area_sg/pix_size
  xpos0 = dd->getxpos(si_id)
  ypos = dd->getypos(si_id)
;  stop
  wv0 = dd->getlam(si_id)
  eff_wv = where(wv0 gt (si_cen-1) and wv0 lt (si_cen+1), /null) 
  spectra = dblarr(n_elements(t_ref), n_elements(eff_wv), n_elements(ypos))
  wv = wv0[eff_wv]*1d0

  fit_res = fltarr(n_elements(t_ref), n_elements(ypos), 4)

  xpos = fltarr(n_elements(t_ref))
  win_sz = [8d2, 8d2]
  gap = 30
  nx = 2.
  ny = 2.
  xmar = [100, 50]
  ymar = [80, 70]
  im_xsz = (win_sz[0] - gap - total(xmar))/nx
  im_ysz = (win_sz[1] - gap - total(ymar))/ny
  xi = xmar[0] + findgen(nx)*(im_xsz + gap)
  yi = win_sz[1] - ymar[1] - findgen(ny)*(im_ysz + gap) - im_ysz
  w01 = window(dim=win_sz)
  im01 = objarr(4)
  for j=0, n_elements(t_ref)-1 do begin
    match = (where(abs(sg_time - t_ref0[j]) eq min(abs(sg_time - t_ref0[j])), /null))[0]
    xpos[j] = xpos0[match]
    sji_ind = (where(abs(sji_time-sg_Time[match]) eq $
                     min(abs(sji_time-sg_time[match])), /null))[0]
    sji_xp = (findgen(sji_indices[sji_ind].naxis1)-sji_indices[sji_ind].crpix1) $
              *sji_indices[sji_ind].cdelt1 + sji_indices[sji_ind].crval1
    sji_yp = (findgen(sji_indices[sji_ind].naxis2)-sji_indices[sji_ind].crpix2) $
              *sji_indices[sji_ind].cdelt2 + sji_indices[sji_ind].crval2
    im01[j] = image(sji_data[*, *, sji_ind], sji_xp, sji_yp, /current, /dev, $
                    pos = [xi[j mod 2], yi[j/2], xi[j mod 2]+im_xsz, yi[j/2]+im_ysz], $ 
                    axis=2, min=5, max=3d1, xr=[-5, 5], yr=[-80, -70], $
                    xtickdir=1, ytickdir=1, $
                    xtitle=(j/2) eq 1 ? 'Solar X (arcsec)' : '', $
                    xtickformat=(j/2) eq 1 ? '(i0)' : '(a1)', $
                    ytitle=(j mod 2) eq 0 ? 'Solar Y (arcsec)' : '', $
                    ytickformat=(j mod 2) eq 0 ? '(i0)' : '(a1)', $
                    font_size=13, rgb_table=iris1400_ct)
    t011 = text(xi[j mod 2]+22, yi[j/2]+im_ysz-22, 'IRIS SJI '+t_ref[j], $
                /current, /dev, $
               font_size=15, color='black', align=0, vertical_align=1)
    t012 = text(xi[j mod 2]+20, yi[j/2]+im_ysz-20, 'IRIS SJI '+t_ref[j], $
                /current, /dev, $
               font_size=15, color='white', align=0, vertical_align=1)
               
    p01 = plot(xpos[j]+[0, 0], minmax(ypos), '--b3', transp=50, over=im01[j])

    if abs(sg_time[match] - t_ref0[j]) gt 1 then continue 
    spec0 = (dd->getvar(si_id))[*, *, match]
    spec = dd->descale_array(spec0)*1.*flux_per_dn/exp_time[match]
    spectra[j, *, *] = spec[eff_wv, *]
    ys = (size(spec))[2]
    spec_max = reform(max(spectra[j, *, *], dim=2))
    for k=0, ys-1 do begin
      res = gaussian_fit(wv, (spectra[j, *, k])/max(spectra[j, *, k]), $
                         si_cen=si_cen, w_th_si=w_th_si0, w_inst=w_inst0) ;;amp0, cen0, wid0, lev0,
                         ;gaussian, not FWHM
      fit_res[j, k, *] = res
    endfor
;    print, minmax(fit_res[j, *, 0])
    fit_res[j, *, 0] *= spec_max
    fit_res[j, *, 3] *= spec_max
  endfor
 
  amplitude = reform(fit_res[*, *, 0])
  velocity = reform((fit_res[*, *, 1] - si_cen)/si_cen*3d5) ; in km/s
  nth_width = (sqrt(reform(fit_res[*, *, 2])^2. - w_inst0^2. - w_th_si0^2.)>0d0)*sqrt(2)/si_cen*3d5 ; gaussian -> 1/e width

  save, spectra, wv, amplitude, velocity, nth_width, xpos, ypos, t_ref, fit_res, $
        filename = 'SI_IV_'+si_cen_str+'_analysis.sav'

  obj_destroy, dd
  t0 = 3
  y0 = 270
  p021 = plot(wv, spectra[t0, *, y0], hist=1, $
              pos = [0.2, 0.15, 0.9, 0.9], $
              xtitle='Wavelength ($\AA$)', ytitle='Intensity (erg s$^{-1} cm^{-2} sr^{-1}$)', $
              font_size=15, xr=minmax(wv), $
              title='t = '+t_ref[t0]+',   y = '+string(ypos[y0], f='(f5.1)')+' arcsec ('+string(y0, f='(i0)')+' pix)')
  p022 = plot(wv, gaussian(wv, fit_res[t0, y0, *]), '-r3', transp=50, over=p021)
  p023 = plot(si_cen*[1, 1], p021.yr, '-b', over=p021)
;endfor
  win_sz = [15, 5]*1d2
  gap = 150
  nx = 3.
  ny = 1.
  xmar = [80, 80]
  ymar = [80, 70]
  im_xsz = (win_sz[0] - gap*(nx-1) - total(xmar))/nx
  im_ysz = (win_sz[1] - gap*(ny-1) - total(ymar))/ny
  xi = xmar[0] + findgen(nx)*(im_xsz + gap)
  yi = ymar[0]
;
  w03 = window(dim=[win_sz])
  im31 = objarr(3)
  for ii=0, 2 do begin
    obj = ii eq 0 ? amplitude : (ii eq 1 ? velocity : nth_width)
    title = 'Si IV '+si_cen_str+' '+$
            (['Amplitude', 'Doppler Velocity', 'Nonthermal Width'])[ii]
    cb_title = (['Intensity (erg s$^{-1} cm^{-2} sr^{-1}$)', $
              'Velocity (km s$^{-2}$)', $
              'Velocity (km s$^{-2}$)'])[ii]
    ct = ii eq 0 ? iris1400_ct : (ii eq 1 ? colortable(70, /rev) : 4)   
    min = ([0, -30, 0])[ii]
    max = ([2.3d2, 30, 40])[ii]           
    im31[ii] = image(obj, findgen(4)-0.5, ypos, axis=2, /current, /dev, $
                    aspect_ratio=0, pos=[xi[ii], yi, xi[ii]+im_xsz, yi+im_ysz], $
                    rgb_table=ct, yr=[-80, -70], min=min, max=max, $
                    xtickname=strmid(t_ref, 5, 6, /rev), xtickval=findgen(4), $
                    ytitle=ii eq 0 ? 'Solar Y (arcsec)' : '', xtitle='Time', $
                    xtickdir=1, ytickdir=1, font_size=13, xminor=0, title=title)
    im31[ii].xticklen = im31[ii].yticklen                   
    cb31 = colorbar(target=im31[ii], /relative, pos=[1, 0, 1.05, 1], $
                    orient=1, textpos=1, /border, title=cb_title, font_size=10)
  endfor
;
  w01.save, 'SJI_and_slit_'+si_cen_str+'.png', resol=200
  p021.save, 'Fit_result_'+si_cen_str+'.png', resol=200
  w03.save, 'parameters_'+si_cen_str+'.png', resol=200
end