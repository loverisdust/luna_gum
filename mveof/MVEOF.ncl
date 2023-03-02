load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
; ==============================================================
begin
optEOF = True       
optEOF@jopt = 0 

optETS = True
optETS@jopt = 0

latS   =  10.
latN   =  50. 
lonL   =  110. 
lonR   =  140.

neof   =2
yrStrt = 197901
yrLast = 202012

yrStrt1 = 197901
yrLast1 = 202012

f =addfile("~/data/precip.mon.mean.nc","r")
fu=addfile("~/data/uwnd.mon.mean.nc","r")
fv=addfile("~/data/vwnd.mon.mean.nc","r")
fp =addfile("~/data/sst.mon.mean.nc","r")

;;;;========读取 数据========================================
TIME   = f->time
YYYY   = cd_calendar(TIME,-1)                
iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
iYYYY1 = ind(YYYY.ge.yrStrt1 .and. YYYY.le.yrLast1)
; a      = short2flt(f->precip(iYYYY,:,:))
; b      = short2flt(f->precip(iYYYY1,:,:))
pre=f->precip(iYYYY,{latS:latN},{lonL:lonR})
pre=dim_standardize_n_Wrap(pre, 0, 0)
p=month_to_season(pre, "JJA")
; xclm   = clmMonTLL(b)
; xAom   = calcMonAnomTLL(a,xclm)
; xAom@long_name = "precip anomalies"
; gpcp      = xAom(0::12,:,:)
delete(TIME)
delete(YYYY)
delete(iYYYY)
delete(iYYYY1)
; delete(a)
; delete(b)
; delete(xAom)
; delete(xclm)

TIME   = fu->time
YYYY   = cd_calendar(TIME,-1)
lon    = fu->lon
lat    = fu->lat          
iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
iYYYY1 = ind(YYYY.ge.yrStrt1 .and. YYYY.le.yrLast1)
uwnd=short2flt(fu->uwnd(iYYYY,{850},{latS:latN},{lonL:lonR}))
uwnd=dim_standardize_n_Wrap(uwnd, 0, 0)
u   =month_to_season(uwnd, "JJA")
; a      = short2flt(fu->uwnd(iYYYY,{850},:,:))
; b      = short2flt(fu->uwnd(iYYYY1,{850},:,:))
; xclm   = clmMonTLL(b);长期平均值
; xAom   = calcMonAnomTLL(a,xclm);长期异常值
; xAom@long_name = "uwnd anomalies"
; uwnd      = xAom(0::12,:,:) 
delete(TIME)
delete(YYYY)
delete(iYYYY)
delete(iYYYY1)
; delete(a)
; delete(b)
; delete(xAom)
; delete(xclm)

TIME   = fv->time
YYYY   = cd_calendar(TIME,-1)                
iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
iYYYY1 = ind(YYYY.ge.yrStrt1 .and. YYYY.le.yrLast1)
vwnd=short2flt(fv->vwnd(iYYYY,{850},{latS:latN},{lonL:lonR}))
vwnd=dim_standardize_n_Wrap(vwnd, 0, 0)
v   =month_to_season(vwnd,"JJA")
; a      = short2flt(fv->vwnd(iYYYY,{850},:,:))
; b      = short2flt(fv->vwnd(iYYYY1,{850},:,:))
; xclm   = clmMonTLL(b)
; xAom   = calcMonAnomTLL(a,xclm)
; xAom@long_name = "vwnd anomalies"
; vwnd      = xAom(0::12,:,:)
delete(TIME)
delete(YYYY)
delete(iYYYY)
delete(iYYYY1)
; delete(a)
; delete(b)
; delete(xAom)
; delete(xclm)

; TIME   = fp->time
; YYYY   = cd_calendar(TIME,-1)                
; iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
; iYYYY1 = ind(YYYY.ge.yrStrt1 .and. YYYY.le.yrLast1)
; a      = short2flt(fp->sst(iYYYY,:,:))
; b      = short2flt(fp->sst(iYYYY1,:,:))
; xclm   = clmMonTLL(b)
; xAom   = calcMonAnomTLL(a,xclm)
; xAom@long_name = "sst anomalies"
; sst      = xAom(0::12,:,:) 
; delete(TIME)
; delete(YYYY)
; delete(iYYYY)
; delete(iYYYY1)
; delete(a)
; delete(b)
; delete(xAom)
; delete(xclm)

;;;======================插值======================================
;lat1   = sst&lat
;precip = linint2_Wrap(gpcp&lon,gpcp&lat,gpcp, True, lon,lat(::-1),0)
;nosst  = linint2_Wrap(sst&lon,lat1(::-1),sst, True, lon,lat(::-1),0)
; =================================================================
sqrt(cos(lat))  ; [or sqrt(gw) ] 权重
; =================================================================
nlat=dimsizes(uwnd(0,:,0))
nlon=dimsizes(uwnd(0,0,:))
re   = 6.37122e06
rad  = 4.0 * atan(1.0) / 180.0
con  = re * rad
clat = cos(lat * rad)           ; cosine of latitude   
dlon = (lon(2) - lon(1))        ; assume dlon is constant
dlat = (lat(2) - lat(1))        ; assume dlat is constant
dx   = con * dlon * clat        ; dx at each latitude
dy   = con * dlat               ; dy is constant
dxdy = dy * dx                  ; dydx(nlat)
wgt  = new((/nlat, nlon/), typeof(uwnd))
;wgt  = conform (wgt,dxdy,0)
wgt!0   = "lat"
wgt!1   = "lon"
; wgt&lat =  lat
; wgt&lon =  lon

;s = nosst

do n=0,dimsizes(u(:,0,0))-1
   u(n,:,:)   = u(n,:,:)*wgt(:,:)
   v(n,:,:)   = v(n,:,:)*wgt(:,:)
   p(n,:,:)   = p(n,:,:)*wgt(:,:)
   ;s(n,:,:)   = nosst(n,:,:)*wgt(:,:)
end do
u = u(:,::-1,:)
v = v(:,::-1,:)

u = dim_standardize_n_Wrap(u, 0,0)
v = dim_standardize_n_Wrap(v, 0,0)
p = dim_standardize_n_Wrap(p, 0,0)
;s = dim_standardize_n_Wrap(s, 0,0)
;*******************Write into one dataset*****************************
;               Combine the normalized data into one variable                                                       ;
;**********************************************************************
u1      = u(:,{latS:latN},{lonL:lonR})
v1      = v(:,{latS:latN},{lonL:lonR})
p1      = p(:,{latS:latN},{lonL:lonR})
;s1      = s(:,{latS:latN},{lonL:lonR})
dimw    = dimsizes(u1)
mtim    = dimw(0)
mlat    = dimw(1)
mlon    = dimw(2)
cdata   = new ((/4*mlat,4*mlon,mtim/), typeof(u1), getFillValue(u1))

do m1=0,mlat-1
    do m2=0,mlon-1
        cdata(m1      ,        m2,:) = (/u1(:,m1,m2)/)
;        copy_VarMeta(u1, cdata(m1       ,       m2,:))
        cdata(m1+mlat ,   m2+mlon,:) = (/v1(:,m1,m2)/)
;        copy_VarMeta(u1, cdata(m1+mlat,  m2+mlon,:))
        cdata(m1+2*mlat,m2+2*mlon,:) = (/p1(:,m1,m2)/)
;        copy_VarMeta(u1, cdata(m1+2*mlat,m2+2*mlon,:))
        cdata(m1+3*mlat,m2+3*mlon,:) = (/s1(:,m1,m2)/)
;        copy_VarMeta(u1, cdata(m1+3*mlat,m2+3*mlon,:))
    end do 
end do
time       = ispan(1979,2021,1)
cdata!2    = "time"
;cdata&time = time
cdata!0    = "lat"
cdata!1    = "lon"
printVarSummary(cdata)
;***************************************t*****************************
;               Compute **combined** EOF; Sign of EOF is arbitrary                                                      ;
;**********************************************************************
eof_cdata    = eofunc_n_Wrap(cdata, neof,optEOF,2)
eof_ts_cdata = eofunc_ts_Wrap(cdata,eof_cdata,optETS)
eof_ts_cdata = dim_standardize_n_Wrap(eof_ts_cdata, 0, 1)
;printVarSummary(eof_cdata)
;printVarSummary(eof_ts_cdata)

;;;;==========================================================
nvar=4                                 ;u,v,precip,sst
ceof = new( (/nvar,neof,mlat,mlon/), typeof(cdata), getFillValue(cdata))
do n=0,neof-1
     ceof(0,n,:,:) = eof_cdata(n,0:mlat-1,0:mlon-1)           ; SLP
     ceof(1,n,:,:) = eof_cdata(n,mlat:mlat*2-1,mlon:2*mlon-1) ; u200
     ceof(2,n,:,:) = eof_cdata(n,2*mlat:mlat*3-1,2*mlon:mlon*3-1)  ; u850
     ceof(3,n,:,:) = eof_cdata(n,3*mlat:mlat*4-1,3*mlon:mlon*4-1)  ; v200
end do
ceof!0   = "var"
ceof!1   = "eof"
ceof!2   = "lat"
ceof!3   = "lon"   
ceof&lat = u1&lat
ceof&lon = u1&lon

tval      = new( (/nvar,neof,mlat,mlon/), "float")
nxy       = new( (/nvar,neof,mlat,mlon/), "integer")
regres    = ceof
tval      = ceof
pro       = ceof
do n = 0,neof-1
   regres(0,n,:,:) = regcoef(eof_ts_cdata(n,:),u1(lat|:,lon|:,time|:), tval(0,n,:,:), nxy(0,n,:,:))
   regres(1,n,:,:) = regcoef(eof_ts_cdata(n,:),v1(lat|:,lon|:,time|:), tval(1,n,:,:), nxy(1,n,:,:))
   regres(2,n,:,:) = regcoef(eof_ts_cdata(n,:),p1(lat|:,lon|:,time|:), tval(2,n,:,:), nxy(2,n,:,:))
   regres(3,n,:,:) = regcoef(eof_ts_cdata(n,:),s1(lat|:,lon|:,time|:), tval(3,n,:,:), nxy(3,n,:,:))
end do 
pro = student_t(tval,mtim)
do i = 0,nvar-1
    do n = 0,neof-1
        copy_VarCoords(u1(0,:,:),regres(i,n,:,:))
        copy_VarCoords(u1(0,:,:),pro(i,n,:,:))
    end do 
end do

;;;;============================================
uu    = ceof(0,:,:,:)
vv    = ceof(0,:,:,:)
pruu  = ceof(0,:,:,:)
prvv  = ceof(0,:,:,:)

uu    = where(pro(0,:,:,:).le.0.05.or.pro(1,:,:,:).le.0.05,regres(0,:,:,:),u1@_FillValue)
vv    = where(pro(0,:,:,:).le.0.05.or.pro(1,:,:,:).le.0.05,regres(1,:,:,:),u1@_FillValue)

pruu  = where(pro(0,:,:,:).le.0.05.or.pro(1,:,:,:).le.0.05,u1@_FillValue,regres(0,:,:,:))
prvv  = where(pro(0,:,:,:).le.0.05.or.pro(1,:,:,:).le.0.05,u1@_FillValue,regres(1,:,:,:))
printVarSummary(uu)
;***************************************t*****************************
;               plot now                                              ;
;**********************************************************************
wks=gsn_open_wks("png","/home/luna/data/MV-EOF")
gsn_define_colormap(wks,"GMT_polar")
plot = new(4   ,graphic)
map  = new(neof,graphic)
vcplt= new(neof,graphic)
shad = new(neof,graphic)
sstt = new(neof,graphic)

shres                            = True            
shres@gsnDraw                    = False
shres@gsnFrame                   = False
; set for the plot
shres@cnFillOn            = True
shres@cnLinesOn           = False
shres@cnLineLabelsOn      = False
shres@cnInfoLabelOn       = False
shres@lbLabelBarOn        = False
shres@cnMonoFillPattern   = False
shres@cnLevelSelectionMode    = "ExplicitLevels"
shres@cnLevels                = (/0.01,0.05/)
shres@cnFillPatterns          = (/17,17,-1/)
shres@cnFillColors            = (/"black","red","white"/)
shres@cnFillDotSizeF          = 0.002
shres@gsnRightString      = ""
shres@gsnLeftString       = ""
shad(0)                   = gsn_csm_contour(wks,pro(2,0,:,:),shres)
shad(1)                   = gsn_csm_contour(wks,pro(2,1,:,:),shres)

sstt(0)                   = gsn_csm_contour(wks,pro(3,0,:,:),shres)
sstt(1)                   = gsn_csm_contour(wks,pro(3,1,:,:),shres)

  res                            = True            
  res@gsnDraw                    = False
  res@gsnFrame                   = False
;---------------------------------------------------------------
; set for the map 
;---------------------------------------------------------------
; set for the plot
  res@vcRefAnnoOn             = True
  res@vcRefAnnoOrthogonalPosF = -1.0             ; move ref vector up
  res@vcRefMagnitudeF         = 0.5            ; define vector ref mag
  res@vcRefAnnoString1On      = False
  res@vcRefAnnoString2On      = True
  res@vcRefAnnoFontHeightF    = 0.020
  res@vcRefAnnoString2        = "0.5 m/s"
  res@vcRefLengthF            = 0.030            ; define length of vec ref
  ;res@vcGlyphStyle            = "LineArrow"    ; turn on curly vectors
  ;res@vcVectorDrawOrder       = "PreDraw"
  res@vcMinDistanceF          = 0.020
  res@gsnRightString             = ""
  res@gsnLeftString              = ""
  res@gsnCenterString            = ""
  res@gsnLeftStringFontHeightF = 0.020
  res@gsnRightStringFontHeightF = 0.020
  res@vcLineArrowColor        = "black"
  res@gsnLeftString = ""
  map(0) = gsn_csm_vector(wks,uu(0,:,:),vv(0,:,:),res) 
  res@gsnLeftString = ""
  map(1) = gsn_csm_vector(wks,uu(1,:,:),vv(1,:,:),res)
  
  res@gsnLeftString = ""
  res@gsnRightString= ""
  res@vcRefAnnoOn             = False
  res@vcLineArrowColor        = "gray65"
  vcplt(0) = gsn_csm_vector(wks,pruu(0,:,:),prvv(0,:,:),res) 
  vcplt(1) = gsn_csm_vector(wks,pruu(1,:,:),prvv(1,:,:),res)


  cnres                         = True
  cnres@gsnAddCyclic            = False
  cnres@gsnDraw                 = False
  cnres@gsnFrame                = False
  cnres@cnLinesOn               = False
  cnres@cnFillOn                = True
  cnres@mpOutlineOn           = True
  cnres@mpFillOn              = False
  cnres@tmXTLabelsOn         = False              ; do not draw bottom labels
  cnres@tmXTOn               = False 
  cnres@tmYRLabelsOn         = False              ; do not draw bottom labels
  cnres@tmYROn               = False 
  cnres@mpMinLatF                  = latS
  cnres@mpMaxLatF                  = latN
  cnres@mpMinLonF                  = lonL
  cnres@mpMaxLonF                  = lonR
  cnres@mpCenterLonF               = (lonR+lonL)/2.
  cnres@mpNationalLineThicknessF   = 2
  cnres@pmTickMarkDisplayMode = "Always"
  cnres@cnInfoLabelOn           = False
  cnres@lbLabelBarOn            = False
  cnres@cnLineLabelsOn          = False
  cnres@gsnContourNegLineDashPattern = 1
  cnres@gsnContourZeroLineThicknessF = 0
  cnres@gsnLeftString           = ""
  cnres@gsnRightString          = ""
  cnres@gsnLeftStringFontHeightF = 0.025
  cnres@gsnRightStringFontHeightF = 0.025
  cnres@tmXBLabelFontHeightF  = 0.020
  cnres@tmYRLabelFontHeightF  = 0.020
  cnres@cnLevelSelectionMode  = "ExplicitLevels"
  cnres@cnLevels              = (/-0.8,-0.6,-0.4,-0.2,-0.1,0.1,0.2,0.4,0.6,0.8/)
  cnres@cnFillColors          = (/4,5,6,7,8,0,15,16,17,18,19/)
  cnres@gsnLeftString = "(a) MV-EOF1"
  plot(0) = gsn_csm_contour_map(wks,regres(2,0,:,:), cnres)
  cnres@gsnLeftString = "(b) MV-EOF2"
  cnres@gsnRightString= ""
  plot(1) = gsn_csm_contour_map(wks,regres(2,1,:,:), cnres)
  
;;;=================================================
do i =0,1
    overlay(plot(i),map(i))
    overlay(plot(i),vcplt(i))
    overlay(plot(i),shad(i))
    overlay(plot(i+2),sstt(i))
end do 


resP = True
resP@gsnPanelLabelBar = True
resP@lbOrientation    = "Vertical"
resP@lbLabelFontHeightF = 0.009
resP@lbTitleString    = ""
resP@lbTitleFontHeightF = 0.015
resP@lbLabelFontHeightF = 0.015
resP@gsnPanelXWhiteSpacePercent = 3.
resP@gsnPanelYWhiteSpacePercent = 4.
resP@gsnPanelMainString = ""
gsn_panel(wks,plot,(/2,2/),resP)

end
