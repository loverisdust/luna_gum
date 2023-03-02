load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

load "./shapefile_utils.ncl";mask
; ==============================================================
begin
;>------------------------------------------------------------<
;;              EOF范围
;>------------------------------------------------------------<
latS = 10.0
latN = 40.0
lonL = 100.0
lonR = 130.0

i=0
spyear=0;特殊年份
yrStrt = 1979
yrLast = 2019
ymStrt = yrStrt*100 +  1
ymLast = yrLast*100 + 12


f_34 =addfile("/mnt/f/data/pree_34.nc","r")
f_56 =addfile("/mnt/f/data/pree_56.nc","r")
f_78 =addfile("/mnt/f/data/pree_78.nc","r")

;>------------------------------------------------------------<
;;              数据预处理 precip,u,v
;>------------------------------------------------------------<
precip_34=f_34->pre_34(:,{latS:latN},{lonL:lonR})
precip_56=f_56->pre_56(:,{latS:latN},{lonL:lonR})
precip_78=f_78->pre_78(:,{latS:latN},{lonL:lonR})
;>------------------------------------------------------------<
;;              数据预处理 缺测处理
;>------------------------------------------------------------<
 precip_34@_FillValue=0; 首先将0设为无效值
  replace_ieeenan (precip_34, -9999.0, 0) ;然后用-9999替换无效值0
  precip_34@_FillValue = -9999.0 ;再将无效值设置为-9999

   precip_56@_FillValue=0; 
  replace_ieeenan (precip_56, -9999.0, 0) 
  precip_56@_FillValue = -9999.0 

   precip_78@_FillValue=0; 
  replace_ieeenan (precip_78, -9999.0, 0) 
  precip_78@_FillValue = -9999.0 
;>------------------------------------------------------------<
;;              计算各个数据的9年滑动平均
;>------------------------------------------------------------<
  nave=9
  precip_34_mean = runave_n_Wrap(precip_34, nave, 1, 0)
  precip_56_mean = runave_n_Wrap(precip_56, nave, 1, 0)
  precip_78_mean = runave_n_Wrap(precip_78, nave, 1, 0)
;>------------------------------------------------------------<
;;                    分离年际尺度
;>------------------------------------------------------------<
  
  precip_34_year = precip_34 - precip_34_mean
  precip_56_year = precip_56 - precip_56_mean
  precip_78_year = precip_78 - precip_78_mean

  copy_VarCoords(precip_34, precip_34_year)
  copy_VarCoords(precip_34, precip_34_mean)

  copy_VarCoords(precip_56, precip_56_year)
  copy_VarCoords(precip_56, precip_56_mean)

  copy_VarCoords(precip_78, precip_78_year)
  copy_VarCoords(precip_78, precip_78_mean)


  precip_34=precip_34_year
  precip_56=precip_56_year
  precip_78=precip_78_year

  precip_34=dim_standardize_n_Wrap(precip_34, 0, 0)
  precip_56=dim_standardize_n_Wrap(precip_56, 0, 0)
  precip_78=dim_standardize_n_Wrap(precip_78, 0, 0)


;>------------------------------------------------------------<
;;                    获取各维大小
;>------------------------------------------------------------<
dimw    = dimsizes(precip_34)
year = ispan(1979,2019,1)
mtim    = dimw(0)
mlat    = dimw(1)
mlon    = dimw(2)
neof=3
;>------------------------------------------------------------<
;;                  MVEOF分解
;>------------------------------------------------------------<
cdata   = new ( (/3*mlat,3*mlon,mtim/), typeof(precip_34), getFillValue(precip_34))
do m1=0,mlat-1
    do m2=0,mlon-1
        cdata(m1,       m2,       :)=(/precip_34(:,m1,m2)/)
        copy_VarMeta(precip_34,cdata(m1,       m2,       :) )
        cdata(m1+mlat,  m2+mlon,  :)=(/precip_56(:,m1,m2)/)
        copy_VarMeta(precip_56, cdata(m1+mlat,  m2+mlon,  :))
        cdata(m1+2*mlat,m2+2*mlon,:)=(/precip_78(:,m1,m2)/)
        copy_VarMeta(precip_78, cdata(m1+2*mlat,m2+2*mlon,:))
    end do
end do
;>------------------------------------------------------------<

eof_cdata    = eofunc(cdata , neof, False)
eof_ts_cdata = eofunc_ts(cdata,eof_cdata,False)

printVarSummary(eof_cdata)
printVarSummary(eof_ts_cdata)

nvar=3                                 ;precip
ceof = new( (/nvar,neof,mlat,mlon/), typeof(cdata), getFillValue(cdata))

do n=0,neof-1
     ceof(0,n,:,:) = eof_cdata(n,0:mlat-1,0:mlon-1)           ; precip_34
     ceof(1,n,:,:) = eof_cdata(n,mlat:mlat*2-1,mlon:2*mlon-1) ; precip_56
     ceof(2,n,:,:) = eof_cdata(n,2*mlat:,2*mlon:)             ; precip_78
end do

ceof!0   = "var"
ceof!1   = "eof"
ceof!2   = "lat"
ceof!3   = "lon"   
ceof&lat = precip_34&lat
ceof&lon = precip_34&lon
printVarSummary(ceof)

tzh = eof_cdata@pcvar
lam = tzh * sqrt(2. / dimw(0))       ;计算解释方差误差范围
;>------------------------------------------------------------<
;;             north检验
;>------------------------------------------------------------<
 ntim = dimw(0)               
 sig_pcv = eofunc_north(eof_cdata@pcvar, ntim, False)
 print(sig_pcv)
;>------------------------------------------------------------<
;;                 计算时间序列并标准化
;>------------------------------------------------------------<
eof_ts = eofunc_ts_n_Wrap (cdata, eof_cdata, False,2)
eof_ts_stddev = dim_stddev(eof_ts)       ;计算时间系数标准差
eof_ts = dim_standardize_n(eof_ts,1,1)
eof_ts_dev = dim_stddev_n(eof_ts, 1)
printVarSummary(eof_ts_stddev)
;>------------------------------------------------------------<
;                 将标准差与空间模态相乘
;>------------------------------------------------------------<
do i=0,neof-1 
  ceof(0,i,:,:) = (/ ceof(0,i,:,:) /) * (/ eof_ts_stddev(i) /)
  ceof(1,i,:,:) = (/ ceof(1,i,:,:) /) * (/ eof_ts_stddev(i) /)
  ceof(2,i,:,:) = (/ ceof(2,i,:,:) /) * (/ eof_ts_stddev(i) /)      
end do
;>------------------------------------------------------------<
;;                            绘图
;>------------------------------------------------------------<
 wks=gsn_open_wks("png","/mnt/f/results/eeof 34-56-78-3")
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
res@gsnRightString = ""
;>------------------------------------------------------------<
;第一模态
;>------------------------------------------------------------<
; res@gsnLeftString = "34_SEOF1"
; res@gsnLeftStringFontHeightF = 0.01
; res@gsnLeftStringOrthogonalPosF =0.01
; plot(0) =gsn_csm_contour_map(wks,ceof(0,0,:,:),res)

; res@gsnLeftString = "56_SEOF1"
; res@gsnLeftStringFontHeightF = 0.01
; res@gsnLeftStringOrthogonalPosF =0.01
; plot(1) =gsn_csm_contour_map(wks,ceof(1,0,:,:),res)

; res@gsnLeftString = "78_SEOF1"
; res@gsnLeftStringFontHeightF = 0.01
; res@gsnLeftStringOrthogonalPosF =0.01
; plot(2) =gsn_csm_contour_map(wks,ceof(2,0,:,:),res)

;>------------------------------------------------------------<
;第二模态
;>------------------------------------------------------------<
; res@gsnLeftString = "34_SEOF2"
; res@gsnLeftStringFontHeightF = 0.01
; res@gsnLeftStringOrthogonalPosF =0.01
; plot(0) =gsn_csm_contour_map(wks,ceof(0,1,:,:),res)

; res@gsnLeftString = "56_SEOF2"
; res@gsnLeftStringFontHeightF = 0.01
; res@gsnLeftStringOrthogonalPosF =0.01
; plot(1) =gsn_csm_contour_map(wks,ceof(1,1,:,:),res)

; res@gsnLeftString = "78_SEOF2"
; res@gsnLeftStringFontHeightF = 0.01
; res@gsnLeftStringOrthogonalPosF =0.01
; plot(2) =gsn_csm_contour_map(wks,ceof(2,1,:,:),res)
;>------------------------------------------------------------<
;第三模态
;>------------------------------------------------------------<
res@gsnLeftString = "34_SEOF3"
res@gsnLeftStringFontHeightF = 0.01
res@gsnLeftStringOrthogonalPosF =0.01
plot(0) =gsn_csm_contour_map(wks,ceof(0,2,:,:),res)

res@gsnLeftString = "56_SEOF3"
res@gsnLeftStringFontHeightF = 0.01
res@gsnLeftStringOrthogonalPosF =0.01
plot(1) =gsn_csm_contour_map(wks,ceof(1,2,:,:),res)

res@gsnLeftString = "78_SEOF3"
res@gsnLeftStringFontHeightF = 0.01
res@gsnLeftStringOrthogonalPosF =0.01
plot(2) =gsn_csm_contour_map(wks,ceof(2,2,:,:),res)
;>----------------------------------------------------------<
;                    设置中国地区区域与河流
;>------------------------------------------------------------<
river=True
river@gsLineThicknessF=0.8
river@gsLineColor="blue"
do n=0,1
plotrv=gsn_add_shapefile_polylines(wks,plot(n),"~/cnmap_NetCDF/rivers.nc",river)
draw(plot(n))
frame(wks)
end do
gsn_panel(wks,plot,(/1,3/),False)     ; now draw as one plot

; ;   ;>------------------------------------------------------------<
; ;   ;                  绘制时间序列
; ;   ;>------------------------------------------------------------<
;     res1=res;时间系数曲线
;     res1@vpXF = 0.45
;     res1@vpYF = 0.8
;     res1@vpHeightF = 0.3
;     res1@vpWidthF = 0.55

;     res1@trXMinF = yrStrt
;     res1@trXMaxF = yrLast
;     res1@trYMinF = -3.5
;     res1@trYMaxF = 3.5
;     res1@gsnYRefLine =(/0.5,0,-0.5/)  ;设定参考值
;     res1@gsnYRefLineDashPatterns=(/11,0,11/)
;     res1@gsnYRefLineColors=(/"red","black","red"/)

;     res1@tmXMajorGrid = False
;     res1@tmXMajorGridThicknessF = 0.3
;     res1@tmXMajorGridLineDashPattern = 2
;     res1@tmYMajorGrid = False
;     res1@tmYMajorGridThicknessF = 0.3
;     res1@tmYMajorGridLineDashPattern = 2

;     res1@tmXBMajorLengthF= 0.005
;     res1@tmYLMajorLengthF= 0.005
;     res1@tmYLLabelDeltaF = -0.7
;     res1@tmXBLabelDeltaF = -0.7
;     res1@tmXBLabelFontHeightF = 0.01
;     res1@tmYLLabelFontHeightF = 0.01

;     res1@tmYLMode = "Manual"
;     res1@tmYLTickSpacingF = 1
;     res1@tmYLMinorPerMajor = 1

;     res1@tiXAxisFontHeightF = 0.015
;     res1@tmXTOn = False
;     res1@tmYROn = False
;     res1@pmLegendDisplayMode = "NoCreate"
;     res1@xyLineColors = (/"black"/);,"gray"
;     res1@xyLineThicknesses = 2
;     res1@xyDashPatterns = 0
;     res1@gsnLeftStringOrthogonalPosF = -0.005
;     res1@gsnRightStringOrthogonalPosF = -0.015
;     res1@gsnLeftStringFontHeightF = 0.015

    
;     pcpc1 = new((/1,dimsizes(eof_ts(0,:))/), "double")
;     printVarSummary(pcpc1)

;     do n=0,neof-1
      
;       pcpc1(0,:) = (/eof_ts(n,:)/);根据neof修改
;       res1@gsnLeftString = "(b) PC"+(n+1)
;       plot20 = gsn_csm_xy(wks,year,pcpc1(0,:), res1)

;      draw(plot20)
;      frame(wks)
;      end do
;     ;>------------------------------------------------------------<
; ;   ;                  绘制north检验
; ;   ;>------------------------------------------------------------<

; tzh_n = eof_cdata@pcvar
; ntim = dimw(0)
; lam_n = tzh_n * sqrt(2. / dimw(0))
; x_variance_all = eof_cdata@pcvar

; resv = True
; resv@trXMaxF = 5
; resv@trXMinF = -1
; resv@trYMaxF = 20
; resv@trYMinF = 0
; resv@tmXTOn = False
; resv@tmYROn = False
; resv@tmYLMinorOn = False  
; resv@tmXBMode = "Explicit"
; resv@tmXBValues = (/0,1,2,3,4,5,6,7,8,9/)
; resv@tmXBLabels = (/1,2,3,4,5,6,7,8,9,10/)
; plotn = gsn_csm_y(wks, x_variance_all, resv)
; polyres = True
; polyres@gsMarkerIndex = 1 
; polyres@gsMarkerSizeF = .02
; polyres@gsLineColor = "blue"
; polyres@gsLineThicknessF = 2.5
; centers = new(dimsizes(x_variance_all),graphic)
; l1dum = new(dimsizes(x_variance_all),graphic) 
; l2dum = new(dimsizes(x_variance_all),graphic)
; l3dum = new(dimsizes(x_variance_all),graphic)
; time   = fspan(1,10,10)
; do t=0,dimsizes(x_variance_all)-1
;   centers(t)   = gsn_add_polymarker(wks,plotn,t,x_variance_all(t),polyres)
;   l1dum(t) = gsn_add_polyline(wks,plotn,(/t,t/),(/x_variance_all(t)-lam_n(t),x_variance_all(t)+lam_n(t)/),polyres)
;   l2dum(t) = gsn_add_polyline(wks,plotn,(/t-.15,t+.15/),(/x_variance_all(t)+lam_n(t),x_variance_all(t)+lam_n(t)/),polyres)
;   l3dum(t) = gsn_add_polyline(wks,plotn,(/t-.15,t+.15/),(/x_variance_all(t)-lam_n(t),x_variance_all(t)-lam_n(t)/),polyres)
; end do
; draw(plotn)
; frame(wks)
end