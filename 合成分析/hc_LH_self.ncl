load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
; ==============================================================
begin
latS   =  10.
latN   =  60. 
lonL   =  70. 
lonR   =  140.


i=0
yystr = 1979
yyend = 2019
spyear=(/1980,1981,1983,1984,1987,1989,1992,2012,2013/)
spyear_minus=(/1,2,4,5,8,10,13,33,34/)

season="JJA"
siglvl=0.05
sig=siglvl
wgty=1.0

f1=addfile("~/data/pre_mon.nc","r")
f2=addfile("~/data/uwnd.mon.mean.nc","r")
f3=addfile("~/data/vwnd.mon.mean.nc","r")
f4=addfile("~/data/hgt.mon.mean.nc", "r") 
f5=addfile("~/data/sst.mon.mean.nc", "r") 
f6=addfile("~/data/air.2m.mon.mean.nc", "r")
f7=addfile("~/data/ano2time_025.nc", "r")

;>------------------------------------------------------------<
;;              数据预处理 
;>------------------------------------------------------------<
;               非降水数据挑选
;>------------------------------------------------------------<

TIME   = f2->time

YYYY   = cd_calendar(TIME,-1)/100  
spyear_79=new(dimsizes(spyear)*12,integer)
iYYYY  = ind(YYYY.ge.yystr .and. YYYY.le.yyend)

a=0
do i=0,8
  spyear_79(a:a+11) = ind(YYYY.eq.spyear(i))
  a=a+12
end do
;>------------------------------------------------------------<

TIME_precip = f1->time

pre=f1->Pre(:,:,:)

uwnd_200 = f2->uwnd(iYYYY,{200},:,:)
vwnd_200 = f3->vwnd(iYYYY,{200},:,:)
hgt_200_1 = f4->hgt(iYYYY,{200},:,:)
uwnd_850 = f2->uwnd(iYYYY,{850},:,:)
vwnd_850 = f3->vwnd(iYYYY,{850},:,:)
hgt_850_1 = f4->hgt(iYYYY,{850},:,:)
uwnd_500= f2->uwnd(iYYYY,{500},:,:)
vwnd_500= f3->vwnd(iYYYY,{500},:,:)
hgt_500_1= f4->hgt(iYYYY,{500},:,:)
sst_anom = f5->sst(iYYYY,:,:)
air_anom = f6->air(iYYYY,:,:)

pre_22=f7->Pre(:,:,:)
uwnd_200_22 = f2->uwnd(spyear_79,{200},:,:)
vwnd_200_22 = f3->vwnd(spyear_79,{200},:,:)
hgt_200_1_22 = f4->hgt(spyear_79,{200},:,:)
uwnd_850_22 = f2->uwnd(spyear_79,{850},:,:)
vwnd_850_22 = f3->vwnd(spyear_79,{850},:,:)
hgt_850_1_22 = f4->hgt(spyear_79,{850},:,:)
uwnd_500_22 = f2->uwnd(spyear_79,{500},:,:)
vwnd_500_22 = f3->vwnd(spyear_79,{500},:,:)
hgt_500_1_22 = f4->hgt(spyear_79,{500},:,:)
sst_anom_22 = f5->sst(spyear_79,:,:)
air_anom_22 = f6->air(spyear_79,:,:)

precip=month_to_season(pre, season)
precipvartime=dim_variance_n_Wrap(precip,0)

precip_22=month_to_season(pre_22, season)
precip_2=dim_avg_n_Wrap(precip_22,0)
precip_2vartime=dim_variance_n_Wrap(precip_22,0)



u_200=month_to_season(uwnd_200, season)
u_200vartime=dim_variance_n_Wrap(u_200,0)


u_200_22=month_to_season(uwnd_200_22, season)
u_200_2=dim_avg_n_Wrap(u_200_22,0)
u_200_2vartime=dim_variance_n_Wrap(u_200_22,0)


v_200=month_to_season(vwnd_200, season)
v_200vartime=dim_variance_n_Wrap(v_200,0)


v_200_22=month_to_season(vwnd_200_22, season)
v_200_2=dim_avg_n_Wrap(v_200_22,0)
v_200_2vartime=dim_variance_n_Wrap(v_200_22,0)


hgt_200=month_to_season(hgt_200_1, season)
hgt_200vartime=dim_variance_n_Wrap(hgt_200,0)


hgt_200_22=month_to_season(hgt_200_1_22, season)
hgt_200_2=dim_avg_n_Wrap(hgt_200_22,0)
hgt_200_2vartime=dim_variance_n_Wrap(hgt_200_22,0)


u_850=month_to_season(uwnd_850, season)
u_850vartime=dim_variance_n_Wrap(u_850,0)


u_850_22=month_to_season(uwnd_850_22, season)
u_850_2=dim_avg_n_Wrap(u_850_22,0)
u_850_2vartime=dim_variance_n_Wrap(u_850_22,0)


v_850=month_to_season(vwnd_850, season)
v_850vartime=dim_variance_n_Wrap(v_850,0)


v_850_22=month_to_season(vwnd_850_22, season)
v_850_2=dim_avg_n_Wrap(v_850_22,0)
v_850_2vartime=dim_variance_n_Wrap(v_850_22,0)



hgt_850=month_to_season(hgt_850_1, season)
hgt_850vartime=dim_variance_n_Wrap(hgt_850,0)

hgt_850_22=month_to_season(hgt_850_1_22, season)
hgt_850_2=dim_avg_n_Wrap(hgt_850_22,0)
hgt_850_2vartime=dim_variance_n_Wrap(hgt_850_22,0)


u_500=month_to_season(uwnd_500, season)
u_500vartime=dim_variance_n_Wrap(u_500,0)


u_500_22=month_to_season(uwnd_500_22, season)
u_500_2=dim_avg_n_Wrap(u_500_22,0)
u_500_2vartime=dim_variance_n_Wrap(u_500_22,0)


v_500=month_to_season(vwnd_500, season)
v_500vartime=dim_variance_n_Wrap(v_500,0)


v_500_22=month_to_season(vwnd_500_22, season)
v_500_2=dim_avg_n_Wrap(v_500_22,0)
v_500_2vartime=dim_variance_n_Wrap(v_500_22,0)


hgt_500=month_to_season(hgt_500_1, season)
hgt_500vartime=dim_variance_n_Wrap(hgt_500,0)


hgt_500_22=month_to_season(hgt_500_1_22, season)
hgt_500_2=dim_avg_n_Wrap(hgt_500_22,0)
hgt_500_2vartime=dim_variance_n_Wrap(hgt_500_22,0)


sst=month_to_season(sst_anom, season)
sstvartime=dim_variance_n_Wrap(sst,0)

sst_22=month_to_season(sst_anom_22, season)
sst_2=dim_avg_n_Wrap(sst_22,0)
sst_2vartime=dim_variance_n_Wrap(sst_22,0)


air=month_to_season(air_anom, season)
airvartime=dim_variance_n_Wrap(air,0)

air_22=month_to_season(air_anom_22, season)
air_2=dim_avg_n_Wrap(air_22,0)
air_2vartime=dim_variance_n_Wrap(air_22,0)

dimw    = dimsizes(precip)
year = ispan(1979,1979+dimw(0)-1,1)

mtim    = dimw(0)
mlat    = dimw(1)
mlon    = dimw(2)

;>------------------------------------------------------------<
;;              计算各个数据的异常场
;>------------------------------------------------------------<



  precip_mean=dim_avg_n_Wrap(precip, 0)
  precip_zong_ano=dim_rmvmean_n_Wrap(precip,0)
  precip_ano=dim_avg_n_Wrap(precip_zong_ano(spyear_minus,:,:),0)
  copy_VarCoords(precip_zong_ano(0,:,:), precip_ano)
  ; printVarSummary(precip_ano)
  ; lat1=precip_ano&lat
  ; printVarSummary(lat1)
 
  u_200_mean=dim_avg_n_Wrap(u_200, 0)
  u_200_zong_ano=dim_rmvmean_n_Wrap(u_200,0)
  u_200_ano=dim_avg_n_Wrap(u_200_zong_ano(spyear_minus,:,:),0)
  copy_VarCoords(u_200_zong_ano(0,:,:), u_200_ano)
  ; printVarSummary(u_200_ano)
  ; lat2=u_200_ano&lat
  ; printVarSummary(lat2)

  v_200_mean=dim_avg_n_Wrap(v_200, 0)
  v_200_zong_ano=dim_rmvmean_n_Wrap(v_200,0)
v_200_ano=dim_avg_n_Wrap(v_200_zong_ano(spyear_minus,:,:),0)
  copy_VarCoords(v_200_zong_ano(0,:,:), v_200_ano)

  hgt_200_mean=dim_avg_n_Wrap(hgt_200, 0)
  hgt_200_zong_ano=dim_rmvmean_n_Wrap(hgt_200,0)
hgt_200_ano=dim_avg_n_Wrap(hgt_200_zong_ano(spyear_minus,:,:),0)
  copy_VarCoords(hgt_200_zong_ano(0,:,:), hgt_200_ano)

  u_850_mean=dim_avg_n_Wrap(u_850, 0)
  u_850_zong_ano=dim_rmvmean_n_Wrap(u_850,0)
u_850_ano=dim_avg_n_Wrap(u_850_zong_ano(spyear_minus,:,:),0)
  copy_VarCoords(u_850_zong_ano(0,:,:), u_850_ano)

  v_850_mean=dim_avg_n_Wrap(v_850, 0)
  v_850_zong_ano=dim_rmvmean_n_Wrap(v_850,0)
v_850_ano=dim_avg_n_Wrap(v_850_zong_ano(spyear_minus,:,:),0)
  copy_VarCoords(v_850_zong_ano(0,:,:), v_850_ano)  

  hgt_850_mean=dim_avg_n_Wrap(hgt_850, 0)
  hgt_850_zong_ano=dim_rmvmean_n_Wrap(hgt_850,0)
hgt_850_ano=dim_avg_n_Wrap(hgt_850_zong_ano(spyear_minus,:,:),0)
  copy_VarCoords(hgt_850_zong_ano(0,:,:), hgt_850_ano)

  u_500_mean=dim_avg_n_Wrap(u_500, 0)
  u_500_zong_ano=dim_rmvmean_n_Wrap(u_500,0)
u_500_ano=dim_avg_n_Wrap(u_500_zong_ano(spyear_minus,:,:),0)
  copy_VarCoords(u_500_zong_ano(0,:,:), u_500_ano)

  v_500_mean=dim_avg_n_Wrap(v_500, 0)
  v_500_zong_ano=dim_rmvmean_n_Wrap(v_500,0)
v_500_ano=dim_avg_n_Wrap(v_500_zong_ano(spyear_minus,:,:),0)
  copy_VarCoords(v_500_zong_ano(0,:,:), v_500_ano)

  hgt_500_mean=dim_avg_n_Wrap(hgt_500, 0)
  hgt_500_zong_ano=dim_rmvmean_n_Wrap(hgt_500,0)
hgt_500_ano=dim_avg_n_Wrap(hgt_500_zong_ano(spyear_minus,:,:),0)
  copy_VarCoords(hgt_500_zong_ano(0,:,:), hgt_500_ano)


  sst_mean=dim_avg_n_Wrap(sst, 0)
  sst_zong_ano=dim_rmvmean_n_Wrap(sst,0)
sst_ano=dim_avg_n_Wrap(sst_zong_ano(spyear_minus,:,:),0)
  copy_VarCoords(sst_zong_ano(0,:,:), sst_ano)

  air_mean=dim_avg_n_Wrap(air, 0)
  air_zong_ano=dim_rmvmean_n_Wrap(air,0)
air_ano=dim_avg_n_Wrap(air_zong_ano(spyear_minus,:,:),0)
  copy_VarCoords(air_zong_ano(0,:,:), air_ano)

;>------------------------------------------------------------<
;;                    t_test
;>------------------------------------------------------------<
iflag=True
sigr=0.05

preciptmp=precip({lat|:},{lon|:},{time|:})
precip_2tmp=precip_22({lat|:},{lon|:},{time|:})

ssttmp=sst(lat|:,lon|:,time|:)
sst_2tmp=sst_22(lat|:,lon|:,time|:)

airtmp=air(lat|:,lon|:,time|:)
air_2tmp=air_22(lat|:,lon|:,time|:)

precipEqv=equiv_sample_size(preciptmp, siglvl, 0)
precip_2Eqv=equiv_sample_size(precip_2tmp, siglvl, 0)
; precipN=wgt_areaave(precipEqv, wgty, 1., 0)
; precip_2N=wgt_areaave(precip_2Eqv, wgty, 1., 0)
precip_2N=dim_num(.not.ismissing(precip_2tmp))
precipN=dim_num(.not.ismissing(preciptmp))

sstEqv=equiv_sample_size(ssttmp, siglvl, 0)
sst_2Eqv=equiv_sample_size(sst_2tmp, siglvl, 0)
; sstN=wgt_areaave(sstEqv, wgty, 1., 0)
; sst_2N=wgt_areaave(sst_2Eqv, wgty, 1., 0)
sst_2N=dim_num(.not.ismissing(sst_2tmp))
sstN=dim_num(.not.ismissing(ssttmp))

airEqv=equiv_sample_size(airtmp, siglvl, 0)
air_2Eqv=equiv_sample_size(air_2tmp, siglvl, 0)
; airN=wgt_areaave(airEqv, wgty, 1., 0)
; air_2N=wgt_areaave(air_2Eqv, wgty, 1., 0)
air_2N=dim_num(.not.ismissing(air_2tmp))
airN=dim_num(.not.ismissing(airtmp))

prob_precip = ttest(precip_2,precip_2vartime,precip_2N, precip_mean,precipvartime,precipN, iflag, False) 
prob_sst = ttest(sst_2,sst_2vartime,sst_2N, sst_mean,sstvartime,sstN, iflag, False) 
prob_air = ttest(air_2,air_2vartime,air_2N, air_mean,airvartime,airN, iflag, False) 
printVarSummary(prob_precip)

  lat=pre&lat
  lon=pre&lon
  lon@units="degrees_east"
  lat@units="degrees_north"

  lat_sst=fspan(-89.5,89.5,180)
  lon_sst=fspan(0.5,359.5,360)
  lon_sst@units="degrees_east"
  lat_sst@units="degrees_north"

  lat_air=fspan(-88.54,88.54,94)
  lon_air=fspan(0.0,358.1,192)
  lon_air@units="degrees_east"
  lat_air@units="degrees_north"


  prob_precip!0="lat"
  prob_precip!1="lon"
  prob_precip&lat=lat
  prob_precip&lon=lon

  prob_sst!0="lat"
  prob_sst!1="lon"
  prob_sst&lat=lat_sst
  prob_sst&lon=lon_sst

  prob_air!0="lat"
  prob_air!1="lon"
  prob_air&lat=lat_air
  prob_air&lon=lon_air

; ;>------------------------------------------------------------<
; ;;              挑选出中国区域
; ;>------------------------------------------------------------<
; reshape(precip_ano, (/time*1,nlat,nlon/))
; printVarSummary(precip_ano)
;   shp_filename    = "bou2_4p.shp"
;   opt             = True
;   opt@return_mask = True
;  ;>------------------------------------------------------------< 
;   china_mask   = shapefile_mask_data(precip_ano(0,:,:),shp_filename,opt)
  
;   precip_ano@missing_value = -99999.99
;   precip_ano@_FillValue    = -99999.99
; ;>------------------------------------------------------------<
;   printMinMax(china_mask,0)

;   china_mask_nd = conform_dims(dimsizes(precip_ano),china_mask,(/0,1/))
;   precip_ano_mask = mask(precip_ano,china_mask_nd,1)
;   copy_VarMeta(precip_ano,precip_ano_mask)

; ;---Print information about precip_ano and precip_ano_mask
;   printVarSummary(precip_ano)            
;   printVarSummary(precip_ano_mask)       ; ditto
; ; ;>------------------------------------------------------------<  
;  precip_ano=precip_ano_mask
;>------------------------------------------------------------<
;;                    绘图
;>------------------------------------------------------------<
    wks1y = gsn_open_wks("pdf","~/data/JJA_LH_self")
    res                   = True    
    
    res@gsnDraw           = False         
    res@gsnFrame          = False        
    res@gsnLeftString     = ""
    res@gsnRightString    = ""
    res@gsnCenterString   = ""  
    res@gsnTickMarksOn = True ; 关闭经度标签
    res@vpHeightF = 0.5;设置画布大小
    res@vpWidthF = 0.6

;>------------------------------------------------------------<
;;                    中国地图与河流设置
;>------------------------------------------------------------<
  res@mpMinLatF            = 10
  res@mpMaxLatF            = 60 
  res@mpMinLonF            = 70
  res@mpMaxLonF            = 140
  res@mpDataBaseVersion = "MediumRes"
  res@mpDataSetName = "./database/Earth..4"
  res@mpOutlineBoundarySets   = "NoBoundaries"
  res@mpFillBoundarySets      = "NoBoundaries"
  res@gsnLeftStringOrthogonalPosF = -0.005
  res@gsnRightStringOrthogonalPosF = -0.015

  res@mpOutlineOn           = True         ; Turn on map outlines
  res@mpOutlineSpecifiers   =  (/"China:states","Taiwan"/)
  res@mpAreaMaskingOn = True   ;使能填充覆盖
  res@mpMaskAreaSpecifiers =  (/"China:states","Taiwan"/)   ;China:states
  ; res@mpOceanFillColor = 0     ;用白色填充海洋  0是colormap的索引值
  ; res@mpInlandWaterFillColor = 0  ;用白色填充内陆湖水
  ;res@mpLandFillColor=0
  ;res@mpGeophysicalLineColor=0
  
  res@mpNationalLineColor        = "black"
  res@mpProvincialLineColor      = "black"
  ;res@mpGeophysicalLineColor     = "black"
  res@mpNationalLineThicknessF   = 2
  res@mpProvincialLineThicknessF = 1  
  res@gsnAddCyclic = False

river=True
river@gsLineThicknessF=0.8
river@gsLineColor="blue"
;>------------------------------------------------------------<
    res1 = res      ;时间系数曲线属性
    res2 = res      ;时间系数曲线滑动平均属性
    resc = res      ;创建resc，以绘制合成场的等值线       
    res0 = res      ;置信度属性
    resalt = res 
    resp = res
    resuv = res
  ;>------------------------------------------------------------<
  ;                    绘制置信度与地图
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
    res0@vpHeightF = 0.5;设置画布大小
    res0@vpWidthF = 0.6
    res0@mpGeophysicalLineThicknessF = 0.4
    res0@mpLandFillColor = "white"
    res0@mpCenterLonF = 105
    res0@mpMinLatF = 10
    res0@mpMaxLatF = 60

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
    resc@gsnAddCyclic = False;不受分辨率的干扰
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
    resp@lbBottomMarginF=1;底部偏移量
    resp@lbBoxMinorExtentF=0.4;bar粗细
    resp@lbLeftMarginF=0.2;偏移量
     ;resp@gsnLeftStringOrthogonalPosF = -1
    resp@gsnAddCyclic = True
    resp@cnInfoLabelOn = False
    resp@cnFillOn = True
    resp@cnLinesOn = False
    resp@cnLevelSelectionMode = "ExplicitLevels"
    resp@lbLabelFontHeightF = 0.01
    resp@gsnAddCyclic = False
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
    resuv@vcRefAnnoPerimOn = False           ;关闭边框

    resuv@vcRefAnnoFontHeightF = 0.015      ;字体大小

    resuv@vcRefAnnoString1On =   False      ;设置标注箭头上方不显示注释
    resuv@vcRefAnnoString2On = True         ;设置标注箭头下方显示注释
    resuv@vcRefAnnoString2  ="5m/s"

    resuv@vcRefAnnoSide            = "Top" ; 参考箭头放至图形上方
    resuv@vcRefAnnoOrthogonalPosF  = -1.1 ; 调整其位置
    resuv@vcRefAnnoParallelPosF    = 1.12
    resuv@gsnAddCyclic = False


  ;>------------------------------------------------------------<
  ;                        绘图
  ;>------------------------------------------------------------<
    plot1y = new(4, graphic)
    plot2y = new(3, graphic)
    plot1d = new(3, graphic)
    lres = True
    lres@gsLineColor = "black"
    lres@gsLineThicknessF = 1.5
    ts = True
    ;>------------------------------------------------------------< 
    ;850
    ;>------------------------------------------------------------< 
    ts@txFontHeightF = 0.03;0.04
    ts@txFontThicknessF = 0.7
    
    resp@cnFillPalette = "BlueWhiteOrangeRed"

    resc@cnLevels := (/-0.9,-0.6,-0.3,-0.1,0,0.1,0.3,0.6,0.9/); 
    resuv@vcRefMagnitudeF = 2.0;
    resuv@vcRefAnnoString2 = "2.0m/s"
    resp@cnLevels := (/-0.9,-0.6,-0.3,-0.1,0,0.1,0.3,0.6,0.9/); 
    res0@gsnLeftString = "850hgt_uv_sst_t2m"

      base = gsn_csm_contour_map(wks1y, prob_sst, res0)
      plot1y(0) = gsn_csm_contour(wks1y, sst_ano, resp)
      plot1y(1)  = gsn_csm_contour(wks1y, air_ano,resp)
      plot1y(2) = gsn_csm_vector(wks1y, u_850_ano, v_850_ano, resuv)    ;uv
      plot1y(3) = gsn_csm_contour(wks1y,hgt_850_ano,resc)    ;hgt
      plot1y(3) = ColorNegDashZeroPosContour(plot1y(3),"blue","black","red")

      plotrv=gsn_add_shapefile_polylines(wks1y,base,"~/cnmap_NetCDF/rivers.nc",river)
      overlay(base, plot1y(0))
      overlay(base, plot1y(1))
      overlay(base, plot1y(2))
      overlay(base, plot1y(3))
       draw(base)
      frame(wks1y)
    
    ;>------------------------------------------------------------< 
    ;200
    ;>------------------------------------------------------------< 
    ts@txFontHeightF := 0.03;0.04
    ts@txFontThicknessF := 0.7

      resp@cnFillPalette := "ViBlGrWhYeOrRe"

      resc@cnLevels := (/-2,0,2/);;
      resuv@vcRefMagnitudeF := 2.0;
      resuv@vcRefAnnoString2 := "2.0"
      resp@cnLevels := (/-1.2,-0.8,-0.4,-0.2,0,0.2,0.4,0.8,1.2/);
      res0@gsnLeftString = "200hgt_uv_pre"
      resp@cnFillColors := (/100,94,86,78,51,51,42,38,34,30/)

      base1 = gsn_csm_contour_map(wks1y, prob_precip, res0)
      plot2y(0)  = gsn_csm_contour(wks1y, precip_ano,resp)
      plot2y(1) = gsn_csm_vector(wks1y, u_200_ano,v_200_ano,resuv)    ;uv
      plot2y(2) = gsn_csm_contour(wks1y,hgt_200_ano,resc)    ;hgt
      plot2y(2) = ColorNegDashZeroPosContour(plot2y(2),"blue","black","red")

  
      plotrv=gsn_add_shapefile_polylines(wks1y,base1,"~/cnmap_NetCDF/rivers.nc",river)
      overlay(base1, plot2y(0))
      overlay(base1, plot2y(1))
      overlay(base1, plot2y(2))
       draw(base1)
        frame(wks1y)
    ;>------------------------------------------------------------< 
    ;500
    ;>------------------------------------------------------------< 
    ts@txFontHeightF := 0.1;0.04
    ts@txFontThicknessF := 0.7
    resc@tmYROn = False
    resc@tmXTOn = False
    resc@tmXBMajorLengthF= 0.005
    resc@tmXBMinorOn = False
    resc@tmYLMajorLengthF= 0.005
    resc@tmYLMinorOn = False
    resc@tmXBLabelDeltaF = -0.7
    resc@tmYLLabelDeltaF = -0.7
    ; res0@tmXBLabelFontHeightF = 0.015
    ; res0@tmYLLabelFontHeightF = 0.015
    resc@tmXBLabelFontHeightF = 0.01
    resc@tmYLLabelFontHeightF = 0.01
    ; res0@pmTickMarkDisplayMode = "Always"
    resc@tmXBTickSpacingF = 30
    resc@tmYLTickSpacingF = 15
    resc@mpCenterLonF = 135
    resc@mpMinLatF = 10
    resc@mpMaxLatF = 60
    resuv@vpHeightF = 0.5;设置画布大小
    resuv@vpWidthF = 0.6

      resuv@vcRefMagnitudeF := 2.0;
      resuv@vcRefAnnoString2 := "2.0"

      resuv@gsnLeftString = "500hgt_uv"
 
      base2= gsn_csm_contour_map(wks1y,hgt_500_ano,resc) ;hgt
      plot1d(0) = gsn_csm_vector(wks1y, u_500_ano,v_500_ano,resuv)    ;uv
      base2 = ColorNegDashZeroPosContour(base2,"blue","black","red")

  
      plotrv=gsn_add_shapefile_polylines(wks1y,base2,"~/cnmap_NetCDF/rivers.nc",river)
      overlay(base2, plot1d(0))

       draw(base2)
        frame(wks1y)
;>------------------------------------------------------------< 
      
end
    ; gsn_panel(wks1y, (/base1, base/), (/2,1/), False)