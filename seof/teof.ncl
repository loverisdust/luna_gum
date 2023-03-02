load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

load "./shapefile_utils.ncl";mask
; ==============================================================
begin
;>------------------------------------------------------------<
;;              EOF范围
;>------------------------------------------------------------<
latS = 18
latN = 32
lonL = 107
lonR = 123
neof=2
i=0
spyear=0;特殊年份
yrStrt = 1979
yrLast = 2019
ymStrt = yrStrt*100 +  1
ymLast = yrLast*100 + 12


f_M =addfile("/mnt/f/data/CN05.1_pre_MAM.nc","r")
f_J =addfile("/mnt/f/data/CN05.1_pre_JJA.nc","r")

;>------------------------------------------------------------<
;;              数据预处理 precip,u,v
;>------------------------------------------------------------<
precip_M=f_M->pre(18:,{latS:latN},{lonL:lonR})
precip_J=f_J->pre(18:,{latS:latN},{lonL:lonR})

printVarSummary(precip_M)
printVarSummary(precip_J)

precip_M&time=ispan(1979, 2019, 1)
precip_J&time=ispan(1979, 2019, 1)
;>------------------------------------------------------------<
;;              计算各个数据的9年滑动平均
;>------------------------------------------------------------<
  nave=9
  precip_M_mean = runave_n_Wrap(precip_M, nave, 1, 0)
  precip_J_mean = runave_n_Wrap(precip_J, nave, 1, 0)
;>------------------------------------------------------------<
;;                    分离年际尺度
;>------------------------------------------------------------<
  
  precip_M_year = precip_M - precip_M_mean
  precip_J_year = precip_J - precip_J_mean

  copy_VarCoords(precip_M, precip_M_year)
  copy_VarCoords(precip_M, precip_M_mean)
  copy_VarCoords(precip_J, precip_J_year)
  copy_VarCoords(precip_J, precip_J_mean)


  precip_M=precip_M_year
  precip_J=precip_J_year
;>------------------------------------------------------------<
;;                    获取各维大小
;>------------------------------------------------------------<
dimw    = dimsizes(precip_M)
year = ispan(1979,2019,1)
mtim    = dimw(0)
mlat    = dimw(1)
mlon    = dimw(2)
  precip_M=dim_standardize_n_Wrap(precip_M, 0, 0)
  precip_J=dim_standardize_n_Wrap(precip_J, 0, 0)
;>------------------------------------------------------------<
;;                  EEOF分解
;>------------------------------------------------------------<
; std_M = dim_stddev_n(precip_M,0)  ; local (grid point, station)  standard deviations
; std_J = dim_stddev_n(precip_J,0)

printVarSummary(std_M)

nvar = 2
                    ; T-mode
                    ; 0  ,  1  ,  2  ,  3   ;  dimension numbers
  new_data= new((/2*mtim, mlat, mlon/), typeof(precip_M),getFillValue(precip_M))
  do m=0,mtim-1
  new_data(0,:,:,:) =  precip_M/std_M

  new_data(1,:,:,:) =  precip_J/std_J

  copy_VarCoords(precip_M,new_data(0,:,:,:))
  copy_VarCoords(precip_M,new_data(1,:,:,:))
  
  printVarSummary(new_data)
  NEW_DATA = reshape(new_data,(/2*mtim,mlat,mlon/))

  copy_VarCoords(precip_M,NEW_DATA(0:40,:,:))
  copy_VarCoords(precip_M,NEW_DATA(41:81,:,:))

  printVarSummary(NEW_DATA)


  neof = 3          ; # of desired EOFs
  ndim = 0        ; dimension number for 'time'

  eeof = eofunc_n(NEW_DATA, neof, False, ndim)

  ;copy_VarCoords(precip_M(0,:,:),eeof(0,:,:))

  printVarSummary(eeof)

; x    = wp({lat|:},{lon|:},{time|:})  ;时间维放在最右边
; neof = 2                            ; 保留前1个模态
; eof  = eofunc_Wrap(x, neof, False)    ; EOF分解
; eof = -eof;春季使用
tzh = eeof@pcvar
lam = tzh * sqrt(2. / dimw(0))       ;计算解释方差误差范围
;>------------------------------------------------------------<
;;             north检验(只取一个模态也有用！)
;>------------------------------------------------------------<
 ntim = dimw(0)               
sig_pcv = eofunc_north(eeof@pcvar, ntim, False)
;>------------------------------------------------------------<
;;                 计算时间序列并标准化
;>------------------------------------------------------------<
eof_ts = eofunc_ts_n_Wrap (NEW_DATA, eeof, False,ndim)
eof_ts_stddev = dim_stddev(eof_ts)       ;计算时间系数标准差
eof_ts = dim_standardize_n(eof_ts,1,1)
eof_ts_dev = dim_stddev_n(eof_ts, 1)
printVarSummary(eof_ts)
;>------------------------------------------------------------<
;                 将标准差与空间模态相乘
;>------------------------------------------------------------<
do i=0,2 
  eeof(i,:,:) = (/ eeof(i,:,:) /) * (/ eof_ts_stddev(i) /)  
end do
;>------------------------------------------------------------<
;;                            绘图
;>------------------------------------------------------------<
 wks=gsn_open_wks("png","/mnt/f/results/eeof")
 ;res@cnFillPalette="CBR_coldhot"

 plot = new(neof,graphic)
 plot2= new(neof,graphic)
 


  res                   = True    
    
  res@gsnDraw           = False         
  res@gsnFrame          = False

  res@gsnLeftString     = ""
  res@gsnRightString    = ""
  res@gsnCenterString   = ""  
  res@gsnMaximize= True
  res@gsnTickMarksOn = True ; 关闭经度标签

  res1=res        ;时间系数曲线属性
  
  resc = res      ;创建resc，以绘制回归场的等值线       
  res0 = res      ;置信度属性
  resalt = res 
  resp = res
  resuv = res
  resv=res  
  res@gsnAddCyclic         = False
  ;res@lbLabelBarOn         = False ;关闭图例
  res@cnLinesOn            = False ;  关闭等值线显示
  res@cnInfoLabelOn         = False   ; 关闭右下角信息标签
  res@cnLineLabelsOn        = False   ; 关闭等值线上标签
  res@cnFillOn             = True         ; turn on color fill

  res@cnFillDrawOrder            = "PreDraw"
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnFillPalette := "OceanLakeLandSnow"
  res@cnLevels             = (/-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9/)
  res@cnFillColors := (/166,150,134,122,106,253,68,52,36,20,4/)
  
  res@pmLabelBarWidthF = 0.03
  res@pmLabelBarHeightF = 0.2

  res@lbOrientation = "vertical"
  res@tmXTOn = False
  res@tmYROn = False
  res@tmXBMajorLengthF= 0.005
  res@tmYLMajorLengthF= 0.005
  res@tmXBLabelDeltaF = -0.7
  res@tmYLLabelDeltaF = -0.7
  res@tmXBLabelFontHeightF = 0.01
  res@tmYLLabelFontHeightF = 0.01

  res@pmTickMarkDisplayMode = "Always"
  res@pmLabelBarOrthogonalPosF = 0.05
  res@tmXBTickSpacingF = 5
  res@tmYLTickSpacingF = 3
  ;>------------------------------------------------------------<
  ;                  设置地图 
  ;>------------------------------------------------------------<
  res@mpOutlineOn          = True
  res@mpFillOn             = True
  ;res@mpShapeMode = "FreeAspect"
  res@vpXF = 0.1
  res@vpYF = 0.8
  res@vpHeightF = 0.3
  res@vpWidthF = 0.31
 
  res@mpMinLatF            = 18
  res@mpMaxLatF            = 32
  res@mpMinLonF            = 107
  res@mpMaxLonF            = 123
  res@mpDataBaseVersion = "MediumRes"
  res@mpDataSetName = "Earth..4"
  res@mpOutlineBoundarySets   = "NoBoundaries"
  res@mpFillBoundarySets      = "NoBoundaries"
  res@gsnLeftStringOrthogonalPosF = -0.005
  res@gsnRightStringOrthogonalPosF = -0.015

  res@mpOutlineOn           = True         ; Turn on map outlines
  res@mpOutlineSpecifiers   =  (/"China","Taiwan"/)
  res@mpAreaMaskingOn = True   ;使能填充覆盖
  res@mpMaskAreaSpecifiers =  (/"China","Taiwan"/)   ;China:states
  ; res@mpOceanFillColor = 0     ;用白色填充海洋  0是colormap的索引值
  ; res@mpInlandWaterFillColor = 0  ;用白色填充内陆湖水
  ;res@mpLandFillColor=0
  ;res@mpGeophysicalLineColor=0
  
  res@mpNationalLineColor        = "black"
  res@mpProvincialLineColor      = "black"
  res@mpGeophysicalLineColor     = "black"
  res@mpProvincialLineThicknessF = 2  
  res@mpNationalLineThicknessF   = 1
  res@gsnAddCyclic = False


;>------------------------------------------------------------<
do n=0,neof-1;取多个模态的时候用

res@gsnRightString = ""

x_variance = decimalPlaces(eeof@pcvar,0,True)
res@gsnLeftString = "(a) EOF"+(n+1)+" ("+x_variance(n)+"%"+")"
res@gsnLeftStringFontHeightF = 0.01
res@gsnLeftStringOrthogonalPosF =0.01
plot(n) =gsn_csm_contour_map(wks,eeof(n,:,:),res)
;>----------------------------------------------------------<
;                    设置中国地区区域与河流
;>------------------------------------------------------------<
river=True
river@gsLineThicknessF=0.8
river@gsLineColor="blue"

plotrv=gsn_add_shapefile_polylines(wks,plot(n),"~/cnmap_NetCDF/rivers.nc",river)

draw(plot(n))
frame(wks)
end do
;gsn_panel(wks,plot,(/1,2/),resP)     ; now draw as one plot

;   ;>------------------------------------------------------------<
;   ;                  绘制时间序列
;   ;>------------------------------------------------------------<
    res1=res;时间系数曲线
    res1@vpXF = 0.45
    res1@vpYF = 0.8
    res1@vpHeightF = 0.3
    res1@vpWidthF = 0.55

    res1@trXMinF = yrStrt
    res1@trXMaxF = yrLast
    res1@trYMinF = -3.5
    res1@trYMaxF = 3.5
    res1@gsnYRefLine =(/0.5,0,-0.5/)  ;设定参考值
    res1@gsnYRefLineDashPatterns=(/11,0,11/)
    res1@gsnYRefLineColors=(/"red","black","red"/)

    res1@tmXMajorGrid = False
    res1@tmXMajorGridThicknessF = 0.3
    res1@tmXMajorGridLineDashPattern = 2
    res1@tmYMajorGrid = False
    res1@tmYMajorGridThicknessF = 0.3
    res1@tmYMajorGridLineDashPattern = 2

    res1@tmXBMajorLengthF= 0.005
    res1@tmYLMajorLengthF= 0.005
    res1@tmYLLabelDeltaF = -0.7
    res1@tmXBLabelDeltaF = -0.7
    res1@tmXBLabelFontHeightF = 0.01
    res1@tmYLLabelFontHeightF = 0.01

    res1@tmYLMode = "Manual"
    res1@tmYLTickSpacingF = 1
    res1@tmYLMinorPerMajor = 1

    res1@tiXAxisFontHeightF = 0.015
    res1@tmXTOn = False
    res1@tmYROn = False
    res1@pmLegendDisplayMode = "NoCreate"
    res1@xyLineColors = (/"black"/);,"gray"
    res1@xyLineThicknesses = 2
    res1@xyDashPatterns = 0
    res1@gsnLeftStringOrthogonalPosF = -0.005
    res1@gsnRightStringOrthogonalPosF = -0.015
    res1@gsnLeftStringFontHeightF = 0.015

    
    pcpc1 = new((/1,dimsizes(eof_ts(0,:))/), "double")

    do n=0,neof-1
      
      pcpc1(0,:) = (/eof_ts(0,:)/);PC2
      res1@gsnLeftString = "(b) PC"+(n+1)
      plot20 = gsn_csm_xy(wks,year,pcpc1(n,:), res1)

     draw(plot20)
     frame(wks)
     end do
    ;>------------------------------------------------------------<
;   ;                  绘制north检验
;   ;>------------------------------------------------------------<
eeof_n = eofunc_n(NEW_DATA, neof, False, ndim)
;eof_n = -eof_n;春季使用
tzh_n = eeof_n@pcvar
ntim = dimw(0)
lam_n = tzh_n * sqrt(2. / dimw(0))
x_variance_all = eeof_n@pcvar

resv = True
resv@trXMaxF = 10
resv@trXMinF = -1
resv@trYMaxF = 60
resv@trYMinF = 0
resv@tmXTOn = False
resv@tmYROn = False
resv@tmYLMinorOn = False  
resv@tmXBMode = "Explicit"
resv@tmXBValues = (/0,1,2,3,4,5,6,7,8,9/)
resv@tmXBLabels = (/1,2,3,4,5,6,7,8,9,10/)
plotn = gsn_csm_y(wks, x_variance_all, resv)
polyres = True
polyres@gsMarkerIndex = 1 
polyres@gsMarkerSizeF = .02
polyres@gsLineColor = "blue"
polyres@gsLineThicknessF = 2.5
centers = new(dimsizes(x_variance_all),graphic)
l1dum = new(dimsizes(x_variance_all),graphic) 
l2dum = new(dimsizes(x_variance_all),graphic)
l3dum = new(dimsizes(x_variance_all),graphic)
time   = fspan(1,10,10)
do t=0,dimsizes(x_variance_all)-1
  centers(t)   = gsn_add_polymarker(wks,plotn,t,x_variance_all(t),polyres)
  l1dum(t) = gsn_add_polyline(wks,plotn,(/t,t/),(/x_variance_all(t)-lam_n(t),x_variance_all(t)+lam_n(t)/),polyres)
  l2dum(t) = gsn_add_polyline(wks,plotn,(/t-.15,t+.15/),(/x_variance_all(t)+lam_n(t),x_variance_all(t)+lam_n(t)/),polyres)
  l3dum(t) = gsn_add_polyline(wks,plotn,(/t-.15,t+.15/),(/x_variance_all(t)-lam_n(t),x_variance_all(t)-lam_n(t)/),polyres)
end do
draw(plotn)
frame(wks)
end