load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
; ==============================================================
begin
latS   =  7.
latN   =  52. 
lonL   =  97. 
lonR   =  142.
neof=2
yrStrt = 1979
yrLast = 2010
ymStrt = yrStrt*100 +  1
ymLast = yrLast*100 + 12


f =addfile("~/data/precip.mon.mean.nc","r")
fu=addfile("~/data/uwnd.mon.mean.nc","r")
fv=addfile("~/data/vwnd.mon.mean.nc","r")

TIME   = f->time
YYYY   = cd_calendar(TIME,-1)/100                 
iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
pre=f->precip(iYYYY,{latS:latN},{lonL:lonR})
pre=dim_standardize_n_Wrap(pre, 0, 0)
precip=month_to_season(pre, "JJA")
; a      = short2flt(f->precip(iYYYY,:,:))
; b      = short2flt(f->precip(iYYYY,:,:))
; xclm   = clmMonTLL(b)
; xAom   = calcMonAnomTLL(a,xclm)
; xAom@long_name = "precip anomalies"
; gpcp      = xAom(0::12,:,:)
delete(TIME)
delete(YYYY)
delete(iYYYY)
;delete(a)
;delete(b)
;delete(xAom)
;delete(xclm)

TIME   = fu->time
lon    = fu->lon
lat    = fu->lat  
YYYY   = cd_calendar(TIME,-1)/100                 
iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
uwnd=short2flt(fu->uwnd(iYYYY,{850},{latS:latN},{lonL:lonR}))
uwnd=dim_standardize_n_Wrap(uwnd, 0, 0)
u   =month_to_season(uwnd, "JJA")
delete(TIME)
delete(YYYY)
delete(iYYYY)

TIME   = fv->time
YYYY   = cd_calendar(TIME,-1)/100                 
iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
vwnd=short2flt(fv->vwnd(iYYYY,{850},{latS:latN},{lonL:lonR}))
vwnd=dim_standardize_n_Wrap(vwnd, 0, 0)
v   =month_to_season(vwnd,"JJA")
printVarSummary(v)
printMinMax(v,False)


dimw    = dimsizes(precip)
year = ispan(1979,1979+dimw(0)-1,1)

mtim    = dimw(0)
mlat    = dimw(1)
mlon    = dimw(2)
;>------------------------------------------------------------<
  ;;              计算各个数据的9年滑动平均
  ;>------------------------------------------------------------<
  nave=9
  precip_mean = runave_n_Wrap(precip, nave, 0, 0)

  do m = 0, (nave/2 - 1), 1
      precip_mean(m,:,:) = dim_avg_n_Wrap(precip(0:nave/2+m,:,:), 0)
      precip_mean(dimw(0)-1-m,:,:) = dim_avg_n_Wrap(precip(dimw(0)-nave/2 -1 -m:,:,:), 0)
  end do

  nave=9
  u_mean = runave_n_Wrap(u, nave, 0, 0)

  do m = 0, (nave/2 - 1), 1
      u_mean(m,:,:) = dim_avg_n_Wrap(u(0:nave/2+m,:,:), 0)
      u_mean(dimw(0)-1-m,:,:) = dim_avg_n_Wrap(u(dimw(0)-nave/2 -1 -m:,:,:), 0)
  end do

  nave=9
  v_mean = runave_n_Wrap(v, nave, 0, 0)

  do m = 0, (nave/2 - 1), 1
      v_mean(m,:,:) = dim_avg_n_Wrap(v(0:nave/2+m,:,:), 0)
      v_mean(dimw(0)-1-m,:,:) = dim_avg_n_Wrap(v(dimw(0)-nave/2 -1 -m:,:,:), 0)
  end do
;>------------------------------------------------------------<
  ;;                    分离年际尺度
  ;>------------------------------------------------------------<
  precip_year = precip - precip_mean
  u_year = u - u_mean
  v_year = v - v_mean

  copy_VarCoords(precip, precip_year)
  copy_VarCoords(precip, precip_mean)
  
  copy_VarCoords(u, u_year)
  copy_VarCoords(u, u_mean)

  copy_VarCoords(v, v_year)
  copy_VarCoords(v, v_mean)

  precip=precip_year
  u=u_year
  v=v_year

; ; =================================================================
;  ; 权重与插值
;  ;precip = linint2_Wrap(gpcp&lon,gpcp&lat,gpcp, True, lon,lat(::-1),0)
; ; =================================================================
; nlat=dimsizes(uwnd(0,:,0))
; nlon=dimsizes(uwnd(0,0,:))
; re   = 6.37122e06
; rad  = 4.0 * atan(1.0) / 180.0
; con  = re * rad
; clat = cos(lat * rad)           ; cosine of latitude   
; dlon = (lon(2) - lon(1))        ; assume dlon is constant
; dlat = (lat(2) - lat(1))        ; assume dlat is constant
; dx   = con * dlon * clat        ; dx at each latitude
; dy   = con * dlat               ; dy is constant
; dxdy = dy * dx                  ; dydx(nlat)
; wgt  = new((/nlat, nlon/), typeof(uwnd))
; ;wgt  = conform (wgt,dxdy,0)
; wgt!0   = "lat"
; wgt!1   = "lon"


; do n=0,dimsizes(u(:,0,0))-1
;    u(n,:,:)   = u(n,:,:)*wgt(:,:)
;    v(n,:,:)   = v(n,:,:)*wgt(:,:)
;    ;precip(n,:,:)   = precip(n,:,:)*wgt(:,:)
; end do
; u = u(:,::-1,:)
; v = v(:,::-1,:)

; u = dim_standardize_n_Wrap(u, 0,0)
; v = dim_standardize_n_Wrap(v, 0,0)
; precip = dim_standardize_n_Wrap(precip, 0,0)
;*******************Write into one dataset*****************************
;               Combine the normalized data into one variable                                                       ;
;**********************************************************************
dimw    = dimsizes(precip)
year = ispan(1979,1979+dimw(0)-1,1)

mtim    = dimw(0)
mlat    = dimw(1)
mlon    = dimw(2)
cdata   = new ( (/3*mlat,3*mlon,mtim/), typeof(v), getFillValue(v))
do m1=0,mlat-1
    do m2=0,mlon-1
        cdata(m1,m2,:)=(/precip(:,m1,m2)/)
        copy_VarMeta(precip,cdata(m1,m2,:) )
        cdata(m1+mlat,m2+mlon,:)=(/u (:,m1,m2)/)
        copy_VarMeta(precip, cdata(m1+mlat, m2+mlon,  :))
        cdata(m1+2*mlat,m2+2*mlon,:)=(/v (:,m1,m2)/)
        copy_VarMeta(precip, cdata(m1+2*mlat,m2+2*mlon,:))
    end do 
end do
;***************************************t*****************************
;               Compute **combined** EOF; Sign of EOF is arbitrary                                                      ;
;**********************************************************************
eof_cdata    = eofunc(cdata   , neof, False)
eof_cdata=-eof_cdata;JJA时使用
eof_ts_cdata = eofunc_ts(cdata,eof_cdata,False)
eof_ts_cdata = dim_standardize_n(eof_ts_cdata,1,1);时间序列标准化

printVarSummary(eof_cdata)
printVarSummary(eof_ts_cdata)
nvar=3                                 ;precip,u,v
ceof = new( (/nvar,neof,mlat,mlon/), typeof(cdata), getFillValue(cdata))
do n=0,neof-1
     ceof(0,n,:,:) = eof_cdata(n,0:mlat-1,0:mlon-1)           ; precip
     ceof(1,n,:,:) = eof_cdata(n,mlat:mlat*2-1,mlon:2*mlon-1) ; u
     ceof(2,n,:,:) = eof_cdata(n,2*mlat:,2*mlon:)             ; v
  end do
ceof!0   = "var"
ceof!1   = "eof"
ceof!2   = "lat"
ceof!3   = "lon"   
ceof&lat = precip&lat
ceof&lon = precip&lon
printVarSummary(ceof)
printMinMax(ceof(0,1,:,:),False)
;***************************************t*****************************
;               plot now                                              ;
;**********************************************************************
wks=gsn_open_wks("pdf","~/data/MVEOF_JJA")

plot = new(neof,graphic)
plot2= new(neof,graphic)
gsn_define_colormap(wks,"cmp_b2r")
  res                      = True         
  res@gsnDraw              = False        ; don't draw yet
  res@gsnFrame             = False        ; don't advance frame yet
; res@gsnPolar             = "SH"
  res@gsnAddCyclic         = False
  res@gsnMaximize          = True
  res@mpFillOn             = True        ; turn off map fill
  res@mpOutlineOn          = True
  res@mpMinLatF            = 10
  res@mpMaxLatF            = 50
  res@mpMaxLonF            = 140
  res@mpMinLonF            = 100

  res@cnFillOn             = True         ; turn on color fill
  ;res@cnFillPalette        = "BlueWhiteOrangeRed"  
  res@cnLinesOn            = False        ; True is default
  res@cnLineLabelsOn       = False        ; True is default
  res@lbLabelBarOn         = False        ; turn off individual lb's
  ;res@cnLevelSelectionMode = "ManualLevels"
  ;res@cnMaxLevelValF       = -0.02
  ;res@cnMinLevelValF       = 0.02
  ;res@cnLevelSpacingF      = 0.005
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels             = (/-0.02,-0.015,-0.01,-0.005,0,0.005,0.01,0.015,0.02/)
  ;res@cnFillColors         = (/5,25,38,56,0,0,80,94,108,121/)
; panel plot only resources
  
  vcres                    = True
  vcres@gsnDraw            = False
  vcres@gsnFrame           = False
  vcres@vcRefLengthF       =0.1
  vcres@vcMinDistanceF     =0.02
  vcres@vcRefAnnoOn             = True
  vcres@vcRefAnnoOrthogonalPosF = -1.0             ; move ref vector up
  vcres@vcRefMagnitudeF         = 0.5            ; define vector ref mag
  vcres@vcRefAnnoString1On      = False
  vcres@vcRefAnnoString2On      = True
  vcres@vcRefAnnoFontHeightF    = 0.020
  vcres@vcRefAnnoString2        = "0.5 m/s"
  ;res@vcGlyphStyle            = "LineArrow"    ; turn on curly vectors
  ;res@vcVectorDrawOrder       = "PreDraw"
  vcres@vcMinDistanceF          = 0.020
  vcres@vcLineArrowColor        = "black"
  ;vcres@vcGlyphStyle     = "CurlyVector"

  resP                     = False         ; modify the panel plot
  resP@gsnMaximize         = True         ; large format
  resP@gsnPanelLabelBar    = True         ; add common colorbar
  resP@gsnPaperOrientation = "portrait"   ; force portrait

  resP@txString            = ""

;*******************************************
; first plot
;*******************************************
  do n=0,neof-1
     ;res@tiMainString   = "MVEOF "+(n+1)
     res@gsnRightString = ""

     x_variance = decimalPlaces(eof_cdata@pcvar,0,True)
     res@gsnLeftString = "(a) EOF1"+(n+1)+" ("+x_variance(n)+"%"+")"
     res@gsnLeftStringFontHeightF = 0.015

     plot(n) =gsn_csm_contour_map(wks,ceof(0,n,:,:),res)
     plot2(n)=gsn_csm_vector     (wks,ceof(1,n,:,:), ceof(2,n,:,:), vcres)
     overlay(plot(n),plot2(n))
     draw(plot(n))
     frame(wks)
  end do
  gsn_panel(wks,plot,(/neof,1/),resP)     ; now draw as one plot
  ;>------------------------------------------------------------<
  ;                  绘制时间序列
  ;>------------------------------------------------------------<
  wks=gsn_open_wks("pdf","~/data/PC2_JJA")


    res1=res;时间系数曲线
    res1@vpXF = 0.45
    res1@vpYF = 0.8
    res1@vpHeightF = 0.3
    res1@vpWidthF = 0.55

    res1@trXMinF = yrStrt
    res1@trXMaxF = yrLast
    res1@trYMinF = -3.5
    res1@trYMaxF = 3.5
    res1@gsnYRefLine = 0     ;设定参考值

    res1@tmXMajorGrid = True
    res1@tmXMajorGridThicknessF = 0.3
    res1@tmXMajorGridLineDashPattern = 2
    res1@tmYMajorGrid = True
    res1@tmYMajorGridThicknessF = 0.3
    res1@tmYMajorGridLineDashPattern = 2

    res1@tmXBMajorLengthF= 0.005
    res1@tmYLMajorLengthF= 0.005
    res1@tmYLLabelDeltaF = -0.7
    res1@tmXBLabelDeltaF = -0.7
    res1@tmXBLabelFontHeightF = 0.01
    res1@tmYLLabelFontHeightF = 0.01

    res1@tmYLMode = "Manual"
    res1@tmYLTickSpacingF = 1.0
    res1@tmYLMinorPerMajor = 1

    res1@tiXAxisString = "Year"   ;设置标题
    res1@tiXAxisFontHeightF = 0.015
    res1@tmXTOn = False
    res1@tmYROn = False
;>------------------------------------------------------------<
    pre_toto = dim_avg_n_Wrap(precip, (/1,2/))
    pre_toto = dim_standardize(pre_toto, 1)
    pcpc1 = new((/1,dimsizes(eof_ts_cdata(1,:))/), "double")
    pcpc1(0,:) = (/eof_ts_cdata(1,:)/)
    ; pcpc1(1,:) = (/pc_mean(0,:)/)
    ; pcpc1(2,:) = (/pc_year(0,:)/)
    ; pcpc1(3,:) = (/pre_toto/)
    res1@pmLegendDisplayMode = "Always"
    res1@pmLegendWidthF = 0.07
    res1@pmLegendHeightF = 0.05
    res1@pmLegendOrthogonalPosF = -1.12 ;-1.22
    res1@pmLegendParallelPosF = 0.08
    res1@lgLabelFontHeightF = 0.01
    res1@lgItemOrder = (/2,1,0/)
    res1@lgBoxMinorExtentF = 0.4
    res1@xyExplicitLabels = (/" PC1"/);," Pre_total"
    res1@xyLineColors = (/"black"/);,"gray"
    res1@xyLineThicknesses = 2
    res1@xyDashPatterns = 0

    res1@gsnLeftStringOrthogonalPosF = -0.005
    res1@gsnRightStringOrthogonalPosF = -0.015
    res1@gsnLeftString = "(b) Principle Component"
    res1@gsnLeftStringFontHeightF = 0.015

    ; plres = True
    ; plres@gsLineThicknessF = 0.5
    ; plres@gsLineDashPattern = 2
    ; eof_ts_devk = onedtond(eof_ts_dev, (/dimsizes(eof_ts_dev),dimsizes(year)/))
    ; dum1 = gsn_add_polyline(wks2, plot20, year, eof_ts_devk(0,:), plres)
    ; dum2 = gsn_add_polyline(wks2, plot20, year, -eof_ts_devk(0,:), plres)

    ;wks2 = gsn_open_wks("eps","pre_anom_"+season+"_ave_pc"+nave)
    plot20 = gsn_csm_xy(wks,year,pcpc1, res1)
    draw(plot20)
    frame(wks)
end