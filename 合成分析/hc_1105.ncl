load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
; ==============================================================
begin
latS   =  -20.
latN   =  80. 
lonL   =  0. 
lonR   =  360.


i=0
yystr = 1979
yyend = 2020
spe_h1=2003
spe_h2=2011
season="MAM"
siglvl=0.05
sig=siglvl
wgty=1.0

f1=addfile("~/data/precip.mon.mean.nc","r")
f2=addfile("~/data/uwnd.mon.mean.nc","r")
f3=addfile("~/data/vwnd.mon.mean.nc","r")
f4=addfile("~/data/hgt.mon.mean.nc", "r") 
f5=addfile("~/data/sst.mon.mean.nc", "r") 

;>------------------------------------------------------------<
;;              数据预处理 
;>------------------------------------------------------------<

TIME   = f1->time
YYYY   = cd_calendar(TIME,-1)/100                 
iYYYY  = ind(YYYY.ge.yystr .and. YYYY.le.yyend)
iY     =ind(YYYY.eq.spe_h1 .or. YYYY.eq.2011)

pre_anom_ave=f1->precip(iYYYY,:,:)
uwnd_anom_ave = f2->uwnd(iYYYY,{200},:,:)
vwnd_anom_ave = f3->vwnd(iYYYY,{200},:,:)
hgt_anom_ave = f4->hgt(iYYYY,{200},:,:)
uwnd_anom_ave850 = f2->uwnd(iYYYY,{850},:,:)
vwnd_anom_ave850 = f3->vwnd(iYYYY,{850},:,:)
hgt_anom_ave850 = f4->hgt(iYYYY,{850},:,:)
sst_anom_ave = f5->sst(iYYYY,:,:)

pre_anom_ave_22=f1->precip(iY,:,:)
uwnd_anom_ave_22 = f2->uwnd(iY,{200},:,:)
vwnd_anom_ave_22 = f3->vwnd(iY,{200},:,:)
hgt_anom_ave_22 = f4->hgt(iY,{200},:,:)
uwnd_anom_ave850_22 = f2->uwnd(iY,{850},:,:)
vwnd_anom_ave850_22 = f3->vwnd(iY,{850},:,:)
hgt_anom_ave850_22 = f4->hgt(iY,{850},:,:)
sst_anom_ave_22 = f5->sst(iY,:,:)


pre=dim_standardize_n_Wrap(pre_anom_ave, 0, 0)
precip=month_to_season(pre, season)
precipvartime=dim_variance_n_Wrap(precip,0)

pre_22=dim_standardize_n_Wrap(pre_anom_ave_22, 0, 0)
precip_22=month_to_season(pre_22, season)
precip_2vartime=dim_variance_n_Wrap(precip_22,0)


uwnd_200=dim_standardize_n_Wrap(uwnd_anom_ave, 0, 0)
u_200=month_to_season(uwnd_200, season)
u_200vartime=dim_variance_n_Wrap(u_200,0)

uwnd_200_22=dim_standardize_n_Wrap(uwnd_anom_ave_22, 0, 0)
u_200_22=month_to_season(uwnd_200_22, season)
u_200_2vartime=dim_variance_n_Wrap(u_200_22,0)

vwnd_200=dim_standardize_n_Wrap(vwnd_anom_ave, 0, 0)
v_200=month_to_season(vwnd_200, season)
v_200vartime=dim_variance_n_Wrap(v_200,0)

vwnd_200_22=dim_standardize_n_Wrap(vwnd_anom_ave_22, 0, 0)
v_200_22=month_to_season(vwnd_200_22, season)
v_200_2vartime=dim_variance_n_Wrap(v_200_22,0)

hgt_200_1=dim_standardize_n_Wrap(hgt_anom_ave, 0, 0)
hgt_200=month_to_season(hgt_200_1, season)
hgt_200vartime=dim_variance_n_Wrap(hgt_200,0)

hgt_200_1_22=dim_standardize_n_Wrap(hgt_anom_ave_22, 0, 0)
hgt_200_22=month_to_season(hgt_200_1_22, season)
hgt_200_2vartime=dim_variance_n_Wrap(hgt_200_22,0)

uwnd_850=dim_standardize_n_Wrap(uwnd_anom_ave850, 0, 0)
u_850=month_to_season(uwnd_850, season)
u_850vartime=dim_variance_n_Wrap(u_850,0)

uwnd_850_22=dim_standardize_n_Wrap(uwnd_anom_ave850_22, 0, 0)
u_850_22=month_to_season(uwnd_850_22, season)
u_850_2vartime=dim_variance_n_Wrap(u_850_22,0)

vwnd_850=dim_standardize_n_Wrap(vwnd_anom_ave850, 0, 0)
v_850=month_to_season(vwnd_850, season)
v_850vartime=dim_variance_n_Wrap(v_850,0)

vwnd_850_22=dim_standardize_n_Wrap(vwnd_anom_ave850_22, 0, 0)
v_850_22=month_to_season(vwnd_850_22, season)
v_850_2vartime=dim_variance_n_Wrap(v_850_22,0)


hgt_850_1=dim_standardize_n_Wrap(hgt_anom_ave850, 0, 0)
hgt_850=month_to_season(hgt_850_1, season)
hgt_850vartime=dim_variance_n_Wrap(hgt_850,0)

hgt_850_1_22=dim_standardize_n_Wrap(hgt_anom_ave850_22, 0, 0)
hgt_850_22=month_to_season(hgt_850_1_22, season)
hgt_850_2vartime=dim_variance_n_Wrap(hgt_850_22,0)

sst_anom=dim_standardize_n_Wrap(sst_anom_ave, 0, 0)
sst=month_to_season(sst_anom, season)
sstvartime=dim_variance_n_Wrap(sst,0)

sst_anom_22=dim_standardize_n_Wrap(sst_anom_ave_22, 0, 0)
sst_22=month_to_season(sst_anom_22, season)
sst_2vartime=dim_variance_n_Wrap(sst_22,0)

dimw    = dimsizes(precip)
year = ispan(1979,1979+dimw(0)-1,1)

mtim    = dimw(0)
mlat    = dimw(1)
mlon    = dimw(2)

;>------------------------------------------------------------<
;;              计算各个数据的异常场
;>------------------------------------------------------------<

  lat=fspan(-88.75,88.75,72)
  lon=fspan(0.0,357.5,144)
  lon@units="degrees_east"
  lat@units="degrees_north"

  precip_mean=dim_avg_n(precip, 0)
  precip_2=(precip(24,:,:)+precip(32,:,:))/2
  precip_ano=precip_2-precip_mean
  ;printVarSummary(precip_ano)
  precip_ano!0="lat"
  precip_ano!1="lon"
  precip_ano&lat=lat
  precip_ano&lon=lon
  ; dims=dimsizes(precip_ano)
  ; nlon=dims(0)
  ; nlat=dims(1)
  ; lon@units="Degrees_east"
  ; lat@units="Degrees_north"
  printVarSummary(precip_ano)


  
  lat_u=fspan(-90,90,73)
  lat_u@units="degrees_north"

  u_200_mean=dim_avg_n(u_200, 0)
  u_200_2=(u_200(24,:,:)+u_200(32,:,:))/2
  u_200_ano=u_200_2-u_200_mean
  ;printVarSummary(u_200_ano)
  u_200_ano!0="lat"
  u_200_ano!1="lon"
  u_200_ano&lat=lat_u
  u_200_ano&lon=lon
  printVarSummary(u_200)
  printVarSummary(u_200_ano)


  v_200_mean=dim_avg_n(v_200, 0)
  v_200_2=(v_200(24,:,:)+v_200(32,:,:))/2
  v_200_ano=v_200_2-v_200_mean
  v_200_ano!0="lat"
  v_200_ano!1="lon"
  v_200_ano&lat=lat_u
  v_200_ano&lon=lon
  printVarSummary(v_200_ano)


  hgt_200_mean=dim_avg_n(hgt_200, 0)
      hgt_200_2=(hgt_200(24,:,:)+hgt_200(32,:,:))/2
  hgt_200_ano=hgt_200_2-hgt_200_mean
  hgt_200_ano!0="lat"
  hgt_200_ano!1="lon"
  hgt_200_ano&lat=lat_u
  hgt_200_ano&lon=lon
    printVarSummary(hgt_200_ano)

  u_850_mean=dim_avg_n(u_850, 0)
      u_850_2=(u_850(24,:,:)+u_850(32,:,:))/2
  u_850_ano=u_850_2-u_850_mean
  u_850_ano!0="lat"
  u_850_ano!1="lon"
  u_850_ano&lat=lat_u
  u_850_ano&lon=lon
   printVarSummary(u_850_ano)


  v_850_mean=dim_avg_n(v_850, 0)
      v_850_2=(v_850(24,:,:)+v_850(32,:,:))/2
  v_850_ano=v_850_2-v_850_mean
  v_850_ano!0="lat"
  v_850_ano!1="lon"
  v_850_ano&lat=lat_u
  v_850_ano&lon=lon
printVarSummary(v_850_ano)

  hgt_850_mean=dim_avg_n(hgt_850, 0)
      hgt_850_2=(hgt_850(24,:,:)+hgt_850(32,:,:))/2
  hgt_850_ano=hgt_850_2-hgt_850_mean
  hgt_850_ano!0="lat"
  hgt_850_ano!1="lon"
  hgt_850_ano&lat=lat_u
  hgt_850_ano&lon=lon
printVarSummary(hgt_850_ano)

lat_sst=fspan(-89.5,89.5,180)
lon_sst=fspan(0.5,359.5,360)
lon_sst@units="degrees_east"
lat_sst@units="degrees_north"


  sst_mean=dim_avg_n(sst, 0)
      sst_2=(sst(24,:,:)+sst(32,:,:))/2
  sst_ano=sst_2-sst_mean
  sst_ano!0="lat"
  sst_ano!1="lon"
  sst_ano&lat=lat_sst
  sst_ano&lon=lon_sst

printVarSummary(sst_ano)

;>------------------------------------------------------------<
;;                    t_test
;>------------------------------------------------------------<
iflag=True
sigr=0.05

preciptmp=precip({lat|:},{lon|:},{time|:})
precip_2tmp=precip_22({lat|:},{lon|:},{time|:})

u_200tmp=u_200(lat|:,lon|:,time|:)
u_200_2tmp=u_200_22(lat|:,lon|:,time|:)

v_200tmp=v_200(lat|:,lon|:,time|:)
v_200_2tmp=v_200_22(lat|:,lon|:,time|:)

hgt_200tmp=hgt_200(lat|:,lon|:,time|:)
hgt_200_2tmp=hgt_200_22(lat|:,lon|:,time|:)

u_850tmp=u_850(lat|:,lon|:,time|:)
u_850_2tmp=u_850_22(lat|:,lon|:,time|:)

v_850tmp=v_850(lat|:,lon|:,time|:)
v_850_2tmp=v_850_22(lat|:,lon|:,time|:)

hgt_850tmp=hgt_850(lat|:,lon|:,time|:)
hgt_850_2tmp=hgt_850_22(lat|:,lon|:,time|:)

ssttmp=sst(lat|:,lon|:,time|:)
sst_2tmp=sst_22(lat|:,lon|:,time|:)

precipEqv=equiv_sample_size(preciptmp, siglvl, 0)
precip_2Eqv=equiv_sample_size(precip_2tmp, siglvl, 0)
; precipN=wgt_areaave(precipEqv, wgty, 1., 0)
; precip_2N=wgt_areaave(precip_2Eqv, wgty, 1., 0)
precip_2N=dim_num(.not.ismissing(precip_2tmp))
precipN=dim_num(.not.ismissing(preciptmp))

u_200Eqv=equiv_sample_size(u_200tmp, siglvl, 0)
u_200_2Eqv=equiv_sample_size(u_200_2tmp, siglvl, 0)
; u_200N=wgt_areaave(u_200Eqv, wgty, 1., 0)
;u_200_2N=wgt_areaave(u_200_2Eqv, wgty, 1., 0)
u_200_2N=dim_num(.not.ismissing(u_200_2tmp))
u_200N=dim_num(.not.ismissing(u_200tmp))

v_200Eqv=equiv_sample_size(v_200tmp, siglvl, 0)
v_200_2Eqv=equiv_sample_size(v_200_2tmp, siglvl, 0)
; v_200N=wgt_areaave(v_200Eqv, wgty, 1., 0)
; v_200_2N=wgt_areaave(v_200_2Eqv, wgty, 1., 0)
v_200_2N=dim_num(.not.ismissing(v_200_2tmp))
v_200N=dim_num(.not.ismissing(v_200tmp))

hgt_200Eqv=equiv_sample_size(hgt_200tmp, siglvl, 0)
hgt_200_2Eqv=equiv_sample_size(hgt_200_2tmp, siglvl, 0)
; hgt_200N=wgt_areaave(hgt_200Eqv, wgty, 1., 0)
; hgt_200_2N=wgt_areaave(hgt_200_2Eqv, wgty, 1., 0)
hgt_200_2N=dim_num(.not.ismissing(hgt_200_2tmp))
hgt_200N=dim_num(.not.ismissing(hgt_200tmp))

u_850Eqv=equiv_sample_size(u_850tmp, siglvl, 0)
u_850_2Eqv=equiv_sample_size(u_850_2tmp, siglvl, 0)
; u_850N=wgt_areaave(u_850Eqv, wgty, 1., 0)
; u_850_2N=wgt_areaave(u_850_2Eqv, wgty, 1., 0)
u_850_2N=dim_num(.not.ismissing(u_850_2tmp))
u_850N=dim_num(.not.ismissing(u_850tmp))

v_850Eqv=equiv_sample_size(v_850tmp, siglvl, 0)
v_850_2Eqv=equiv_sample_size(v_850_2tmp, siglvl, 0)
; v_850N=wgt_areaave(v_850Eqv, wgty, 1., 0)
; v_850_2N=wgt_areaave(v_850_2Eqv, wgty, 1., 0)
v_850_2N=dim_num(.not.ismissing(v_850_2tmp))
v_850N=dim_num(.not.ismissing(v_850tmp))

hgt_850Eqv=equiv_sample_size(hgt_850tmp, siglvl, 0)
hgt_850_2Eqv=equiv_sample_size(hgt_850_2tmp, siglvl, 0)
; hgt_850N=wgt_areaave(hgt_850Eqv, wgty, 1., 0)
; hgt_850_2N=wgt_areaave(hgt_850_2Eqv, wgty, 1., 0)
hgt_850_2N=dim_num(.not.ismissing(hgt_850_2tmp))
hgt_850N=dim_num(.not.ismissing(hgt_850tmp))

sstEqv=equiv_sample_size(ssttmp, siglvl, 0)
sst_2Eqv=equiv_sample_size(sst_2tmp, siglvl, 0)
; sstN=wgt_areaave(sstEqv, wgty, 1., 0)
; sst_2N=wgt_areaave(sst_2Eqv, wgty, 1., 0)
sst_2N=dim_num(.not.ismissing(sst_2tmp))
sstN=dim_num(.not.ismissing(ssttmp))

prob_precip = ttest(precip_2,precip_2vartime,precip_2N, precip_mean,precipvartime,precipN, iflag, False) 
prob_u_200 = ttest(u_200_2,u_200_2vartime,u_200_2N, u_200_mean,u_200vartime,u_200N, iflag, False) 
prob_v_200 = ttest(v_200_2,v_200_2vartime,v_200_2N, v_200_mean,v_200vartime,v_200N, iflag, False) 
prob_hgt_200 = ttest(hgt_200_2,hgt_200_2vartime,hgt_200_2N, hgt_200_mean,hgt_200vartime,hgt_200N, iflag, False) 
prob_u_850 = ttest(u_850_2,u_850_2vartime,u_850_2N, u_850_mean,u_850vartime,u_850N, iflag, False) 
prob_v_850 = ttest(v_850_2,v_850_2vartime,v_850_2N, v_850_mean,v_850vartime,v_850N, iflag, False) 
prob_hgt_850 = ttest(hgt_850_2,hgt_850_2vartime,hgt_850_2N, hgt_850_mean,hgt_850vartime,hgt_850N, iflag, False) 
prob_sst = ttest(sst_2,sst_2vartime,sst_2N, sst_mean,sstvartime,sstN, iflag, False) 

  prob_precip!0="lat"
  prob_precip!1="lon"
  prob_precip&lat=lat
  prob_precip&lon=lon

  prob_sst!0="lat"
  prob_sst!1="lon"
  prob_sst&lat=lat_sst
  prob_sst&lon=lon_sst

;>------------------------------------------------------------<
;;                    绘图
;>------------------------------------------------------------<
    wks1y = gsn_open_wks("pdf","~/data/850_1106")
    res                   = True    
    
    res@gsnDraw           = False         
    res@gsnFrame          = False        
    res@gsnLeftString     = ""
    res@gsnRightString    = ""
    res@gsnCenterString   = ""  
    res@gsnTickMarksOn = True ; 关闭经度标签

    res1 = res      ;时间系数曲线属性
    res2 = res      ;时间系数曲线滑动平均属性
    resc = res      ;创建resc，以绘制合成场的等值线       
    res0 = res      ;置信度属性
    resalt = res 
    resp = res
    resuv = res

    res@gsnAddCyclic          = False
    res@cnFillOn              = True 
    ;res@lbLabelBarOn          = False   ;关闭图例
    res@cnLinesOn             = False   ;  关闭等值线显示
    res@cnInfoLabelOn         = False   ; 关闭右下角信息标签
    res@cnLineLabelsOn        = False   ; 关闭等值线上标签
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMinLevelValF = 0.
    res@cnMaxLevelValF = 0.3
    res@cnLevelSpacingF = 0.05
    res@tmXTOn = False
    res@tmYROn = False
    res@tmXBMajorLengthF= 0.005
    res@tmYLMajorLengthF= 0.005
    res@tmXBLabelDeltaF = -0.7
    res@tmYLLabelDeltaF = -0.7
    res@tmXBLabelFontHeightF = 0.01
    res@tmYLLabelFontHeightF = 0.01
    ; res@pmTickMarkDisplayMode = "Always"
    res@pmLabelBarOrthogonalPosF = 0.1
    res@tmXBTickSpacingF = 5
    res@tmYLTickSpacingF = 3

    res@pmLabelBarHeightF = 0.03
    res@gsnSpreadColorStart = 129
    res@gsnSpreadColorEnd = -2

  ;>------------------------------------------------------------<
  ;                    绘制置信度
  ;>------------------------------------------------------------<


    res0@cnFillOn              = True
    res0@lbLabelBarOn          = False   ;关闭图例
    res0@cnLinesOn             = False   ;关闭等值线显示
    res0@cnInfoLabelOn         = False   ;关闭右下角信息标签
    res0@cnLineLabelsOn        = False   ;关闭等值线上标签
    res0@cnMonoFillPattern     = False

    res0@cnLevelSelectionMode  = "ExplicitLevels"
    res0@cnLevels              = sig ; -20与20均对应着0.05的置信度，前者对应负值异常，后者对应正值异常
    res0@cnFillPatterns        = (/17,-1/)
    res0@cnFillColors          = (/1,0/)
    res0@cnFillDrawOrder = "PostDraw"

    res0@gsnLeftStringFontHeightF = 17;20
    res0@gsnLeftStringOrthogonalPosF = 0.01
    res0@mpShapeMode = "FreeAspect"
    res0@vpHeightF = 0.3
    res0@vpWidthF = 0.7
    res0@mpGeophysicalLineThicknessF = 0.4
    res0@mpLandFillColor = "white"
    res0@mpCenterLonF = 180
    res0@mpMinLatF = -25
    res0@mpMaxLatF = 85

    res0@tmYROn = False
    res0@tmXTOn = False
    res0@tmXBMajorLengthF= 0.005
    res0@tmXBMinorOn = False
    res0@tmYLMajorLengthF= 0.005
    res0@tmYLMinorOn = False
    res0@tmXBLabelDeltaF = -0.7
    res0@tmYLLabelDeltaF = -0.7
    ; res0@tmXBLabelFontHeightF = 0.015
    ; res0@tmYLLabelFontHeightF = 0.015
    res0@tmXBLabelFontHeightF = 0.01
    res0@tmYLLabelFontHeightF = 0.01
    ; res0@pmTickMarkDisplayMode = "Always"
    res0@tmXBTickSpacingF = 30
    res0@tmYLTickSpacingF = 15
  ;>------------------------------------------------------------<
  ;;                   绘制地形
  ;>------------------------------------------------------------<      
    resalt@cnFillOn = True
    resalt@gsnAddCyclic = False
    resalt@cnLinesOn = False
    resalt@cnLineLabelsOn = False
    resalt@lbLabelBarOn = False
    resalt@cnInfoLabelOn = False
    resalt@cnFillPalette = "gsltod"
    resalt@cnLevelSelectionMode = "ManualLevels"
    resalt@cnMinLevelValF = 0.
    resalt@cnMaxLevelValF = 9000.
    resalt@cnLevelSpacingF = 100.
  ;>------------------------------------------------------------<
  ;                  位势高度合成场等值线设置
  ;>------------------------------------------------------------<
    resc@cnFillOn              = False 
    resc@lbLabelBarOn          = False   ;关闭图例
    resc@cnLinesOn             = True   ;  关闭等值线显示
    resc@cnInfoLabelOn         = False   ; 关闭右下角信息标签
    resc@cnLineLabelsOn        = True   ; 关闭等值线上标签
    resc@cnLevelSelectionMode = "ExplicitLevels"
    ; resc@cnLineThicknessF = 1.0
    resc@cnLineThicknessF = 2.0
    ; resc@gsnContourZeroLineThicknessF = 0
    resc@cnLineLabelFontHeightF = 0.008
    resc@cnLineLabelFontThicknessF = 0.8
    resc@cnLineLabelBackgroundColor = "white"
    resc@cnLineLabelAngleF = 0.0
    resc@cnLineLabelDensityF = 1
    resc@cnLineLabelInterval = 1
    resc@cnLineLabelPerimOn = True
    resc@cnLineLabelPerimColor = "white"
    resc@cnLineLabelPerimThicknessF = 1.5
  ;>------------------------------------------------------------<
  ;                  sst/降水等值线设置
  ;>------------------------------------------------------------<
    resp@lbLabelBarOn          = True
    ; resp@lbBoxEndCapStyle  = "TriangleBothEnds"
    resp@lbOrientation = "Vertical"
    resp@lbBoxLinesOn = True
    resp@lbBoxLineColor = "black"
    resp@lbBoxLineThicknessF = 0.01
    resp@pmLabelBarWidthF = 0.04
    ; resp@gsnLeftStringOrthogonalPosF = -1
    resp@gsnAddCyclic = True
    resp@cnInfoLabelOn = False
    resp@cnFillOn = True
    resp@cnLinesOn = False
    resp@cnLevelSelectionMode = "ExplicitLevels"
    resp@lbLabelFontHeightF = 0.01
  ;>------------------------------------------------------------<
  ;                  绘制uv
  ;>------------------------------------------------------------<
    resuv@vcGlyphStyle = "CurlyVector"   ;设置向量的显示形状（大小方向）
    resuv@vcGlyphStyle = "LineArrow"   ;设置向量的显示形状（大小方向）
    resuv@vcLineArrowThicknessF = 0.7;1.0
    resuv@vcMinDistanceF = 0.02   ;设置箭头间最小距离

    resuv@vcRefAnnoOn = True
    resuv@vcRefLengthF = 0.045              ;标准长度箭头在单位平方中的大小
    resuv@vcRefMagnitudeF = 0.1           ;标准长度箭头所表示的大小
    resuv@vcRefAnnoBackgroundColor = "white";背景颜色
    resuv@vcRefAnnoPerimOn = True           ;关闭边框

    resuv@vcRefAnnoFontHeightF = 0.015      ;字体大小

    resuv@vcRefAnnoString1On =   False      ;设置标注箭头上方不显示注释
    resuv@vcRefAnnoString2On = True         ;设置标注箭头下方显示注释
    resuv@vcRefAnnoString2  ="0.1 m/s"

    resuv@vcRefAnnoSide            = "Top" ; 参考箭头放至图形上方
    resuv@vcRefAnnoOrthogonalPosF  = -0.12 ; 调整其位置
    resuv@vcRefAnnoParallelPosF    = 0.95 

  ;>------------------------------------------------------------<
  ;                        绘图
  ;>------------------------------------------------------------<
    plot1y = new(3, graphic)
    plot2y = new(3, graphic)
    plot1d = new(3, graphic)
    lres = True
    lres@gsLineColor = "black"
    lres@gsLineThicknessF = 1.5
    ts = True
  ;>------------------------------------------------------------<
    ts@txFontHeightF = 0.03;0.04
    ts@txFontThicknessF = 0.7
    
    resp@cnFillPalette = "BlueWhiteOrangeRed"

    resc@cnLevels := (/-0.1,10/); 
    resuv@vcRefMagnitudeF = 1.0;
    resuv@vcRefAnnoString2 = "1.0"
    resp@cnLevels := (/-0.9,-0.6,-0.3,-0.1,0,0.1,0.3,0.6,0.9/); 
    res0@gsnLeftString = "(b) 850hgt_uv_sst"

      base = gsn_csm_contour_map(wks1y, prob_sst, res0)
      plot1y(0) = gsn_csm_contour(wks1y, sst_ano, resp)
      plot1y(1) = gsn_csm_vector(wks1y, u_850_ano, v_850_ano, resuv)    ;uv
      plot1y(2) = gsn_csm_contour(wks1y,hgt_850_ano,resc)    ;hgt
      plot1y(2) = ColorNegDashZeroPosContour(plot1y(2),"blue","black","red")


      overlay(base, plot1y(0))
      overlay(base, plot1y(1))
      overlay(base, plot1y(2))
       draw(base)
    
  ;>------------------------------------------------------------<
    ; ts@txFontHeightF := 0.03;0.04
    ; ts@txFontThicknessF := 0.7


    ; ; else if (level .eq. 200) then
    ;   ; gsn_define_colormap(wks1y, "ViBlGrWhYeOrRe")
    ;   resp@cnFillPalette := "ViBlGrWhYeOrRe"

    ;   resc@cnLevels := (/-5,-1,20,30/);;
    ;   resuv@vcRefMagnitudeF := 1.5;回归到ALI MCI
    ;   resuv@vcRefAnnoString2 := "1.5"
    ;   resp@cnLevels := (/-1.2,-0.8,-0.4,-0.05,0,0.05,0.4,0.8,1.2/);
    ;   res0@gsnLeftString = "(a) Reg_200hgt_uv_pre"

    ;   resp@cnFillColors := (/100,94,86,78,51,51,42,38,34,30/)
    ;   base1 = gsn_csm_contour_map(wks1y, prob_precip, res0)
    ;   plot2y(0)  = gsn_csm_contour(wks1y, precip_ano,resp)
    ;   plot2y(1) = gsn_csm_vector(wks1y, u_200_ano,v_200_ano,resuv)    ;uv
    ;   plot2y(2) = gsn_csm_contour(wks1y,hgt_200_ano,resc)    ;hgt
    ;   plot2y(2) = ColorNegDashZeroPosContour(plot2y(2),"blue","black","red")

  

    ;   overlay(base1, plot2y(0))
    ;   overlay(base1, plot2y(1))
    ;   overlay(base1, plot2y(2))
    ;    draw(base1)
    ;>------------------------------------------------------------< 
       frame(wks1y)
end
    ; gsn_panel(wks1y, (/base1, base/), (/2,1/), False)