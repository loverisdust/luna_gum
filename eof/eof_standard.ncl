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
i=0
spyear=0;特殊年份
yrStrt = 1979
yrLast = 2020
ymStrt = yrStrt*100 +  1
ymLast = yrLast*100 + 12


f =addfile("~/data/precip.mon.mean.nc","r")
fu=addfile("~/data/uwnd.mon.mean.nc","r")
fv=addfile("~/data/vwnd.mon.mean.nc","r")

;>------------------------------------------------------------<
;;              数据预处理 precip,u,v
;>------------------------------------------------------------<

TIME   = f->time
YYYY   = cd_calendar(TIME,-1)/100                 
iYYYY  = ind(YYYY.ge.yrStrt .and. YYYY.le.yrLast)
pre=f->precip(iYYYY,{latS:latN},{lonL:lonR})
pre=dim_standardize_n_Wrap(pre, 0, 0)
precip=month_to_season(pre, "JJA")
printVarSummary(precip)
delete(TIME)
delete(YYYY)
delete(iYYYY)

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
;>------------------------------------------------------------<
;;        Combine the normalized data into one variable
;>------------------------------------------------------------<
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

;sig_pcv = eofunc_north(eof_cdata@pcvar, mtim, False)
;>-----------------------------t------------------------------<
;;     Compute **combined** EOF; Sign of EOF is arbitrary
;>------------------------------------------------------------<
eof_cdata    = eofunc(cdata , neof, False)

;eof_cdata=-eof_cdata;处理使用

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
;print(ceof)
printMinMax(ceof(0,1,:,:),False)
;print(eof_ts_cdata(0,1))
;>------------------------------------------------------------<
;;                            年份提取
;>------------------------------------------------------------<
; pc=0.8;特殊年份提取点
; print("PC1符合要求的年份有")
; do i=0, 42
;   if (abs(eof_ts_cdata(0,i)).ge.pc) then
;     spyear=i+1979
;     print(spyear)
;     print(eof_ts_cdata(0,i))
;   end if
; end do 

; print("PC2符合要求的年份有")
; do i=0, 42
;   if (abs(eof_ts_cdata(1,i)).ge.pc) then
;     spyear=i+1979
;     print(spyear)
;     print(eof_ts_cdata(1,i))
;   end if
; end do 
;>------------------------------------------------------------<
;;                            绘图
;>------------------------------------------------------------<
 wks=gsn_open_wks("pdf","~/data/MVEOF_JJA")

 plot = new(neof,graphic)
 plot2= new(neof,graphic)
 gsn_define_colormap(wks,"CBR_coldhot")


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

  res@gsnAddCyclic         = False

  ;res@lbLabelBarOn         = False ;关闭图例
  res@cnLinesOn            = False ;  关闭等值线显示
  res@cnInfoLabelOn         = False   ; 关闭右下角信息标签
  res@cnLineLabelsOn        = False   ; 关闭等值线上标签
  res@cnFillOn             = True         ; turn on color fill
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels             = (/-0.05,-0.03,-0.015,-0.01,-0.005,0.005,0.01,0.015,0.03,0.05/)
  ;res@cnFillColors         = (/5,25,38,56,0,0,80,94,108,121/)


  res@tmXTOn = False
  res@tmYROn = False
  res@tmXBMajorLengthF= 0.005
  res@tmYLMajorLengthF= 0.005
  res@tmXBLabelDeltaF = -0.7
  res@tmYLLabelDeltaF = -0.7
  res@tmXBLabelFontHeightF = 0.01
  res@tmYLLabelFontHeightF = 0.01

  res@pmTickMarkDisplayMode = "Always"
  res@lbOrientation="Vertical"


  res@pmLabelBarOrthogonalPosF = 0.1
  res@tmXBTickSpacingF = 5
  res@tmYLTickSpacingF = 3

  res@pmLabelBarHeightF = 0.6

  ; res@gsnSpreadColorStart = 129
  ; res@gsnSpreadColorEnd = -2 
  ;>------------------------------------------------------------<
  ;                  设置地图 
  ;>------------------------------------------------------------<
        
  res@mpOutlineOn          = True
  res@mpFillOn             = True

  res@mpMinLatF            = 10
  res@mpMaxLatF            = 50 
  res@mpMinLonF            = 100
  res@mpMaxLonF            = 140
 

  
  res@gsnLeftStringOrthogonalPosF = -0.005
  res@gsnRightStringOrthogonalPosF = -0.015

  ;>------------------------------------------------------------<
  ;plot0 = gsn_csm_contour_map(wks, eof_cdata(0,:,:), resp1)

  ;>------------------------------------------------------------<
  ;                  绘制uv
  ;>------------------------------------------------------------<
    ; resuv@vcGlyphStyle = "CurlyVector"   ;设置向量的显示形状（大小方向）
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

  resP                     = False         ; modify the panel plot
  resP@gsnMaximize         = True         ; large format
  resP@gsnPanelLabelBar    = True         ; add common colorbar
  resP@gsnPaperOrientation = "portrait"   ; force portrait

  resP@txString            = ""

;>------------------------------------------------------------<
   do n=0,neof-1

     res@gsnRightString = ""

     x_variance = decimalPlaces(eof_cdata@pcvar,0,True)
     res@gsnLeftString = "(a) EOF"+(n+1)+" ("+x_variance(n)+"%"+")"
     res@gsnLeftStringFontHeightF = 0.015

     plot(n) =gsn_csm_contour_map(wks,ceof(0,n,:,:),res)
     plot2(n)=gsn_csm_vector     (wks,ceof(1,n,:,:), ceof(2,n,:,:), resuv)
     overlay(plot(n),plot2(n))
     draw(plot(n))
     frame(wks)
  end do

  ;gsn_panel(wks,plot,(/1,2/),resP)     ; now draw as one plot

;   ;>------------------------------------------------------------<
;   ;                  绘制时间序列
;   ;>------------------------------------------------------------<
  ;wks=gsn_open_wks("pdf","~/data/PC2_MAM")


    res1=res;时间系数曲线
    res1@vpXF = 0.45
    res1@vpYF = 0.8
    res1@vpHeightF = 0.3
    res1@vpWidthF = 0.55

    res1@trXMinF = yrStrt
    res1@trXMaxF = yrLast
    res1@trYMinF = -3.5
    res1@trYMaxF = 3.5
    res1@gsnYRefLine =(/0.8,0,-0.8/)  ;设定参考值
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

    ;res1@tiXAxisString = "Year"   ;设置标题
    res1@tiXAxisFontHeightF = 0.015
    res1@tmXTOn = False
    res1@tmYROn = False

;>------------------------------------------------------------<
    pre_toto = dim_avg_n_Wrap(precip, (/1,2/))
    pre_toto = dim_standardize(pre_toto, 1)

    ; pcpc1 = new((/1,dimsizes(eof_ts_cdata(1,:))/), "double")
    ; pcpc1(0,:) = (/eof_ts_cdata(1,:)/);PC2

    ; pcpc1 = new((/1,dimsizes(eof_ts_cdata(0,:))/), "double")
    ; pcpc1(0,:) = (/eof_ts_cdata(0,:)/);PC1
    
    res1@pmLegendDisplayMode = "NoCreate"
    res1@xyLineColors = (/"black"/);,"gray"
    res1@xyLineThicknesses = 2
    res1@xyDashPatterns = 0
 

    res1@gsnLeftStringOrthogonalPosF = -0.005
    res1@gsnRightStringOrthogonalPosF = -0.015
    res1@gsnLeftString = "(b) Principle Component"
    res1@gsnLeftStringFontHeightF = 0.015

    do n=0,neof-1

   ;   res@gsnRightString = ""

   ;   x_variance = decimalPlaces(eof_cdata@pcvar,0,True)
   ;   res@gsnLeftString = "(a) EOF"+(n+1)+" ("+x_variance(n)+"%"+")"
   ;   res@gsnLeftStringFontHeightF = 0.015

   ;   pcpc1 = new((/1,dimsizes(eof_ts_cdata(n,:))/), "double")
   ;   pcpc1(0,:) = (/eof_ts_cdata(n,:)/);PC2

   ;   plot(n) =gsn_csm_contour_map(wks,ceof(0,n,:,:),res)
   ;   plot2(n)=gsn_csm_vector     (wks,ceof(1,n,:,:), ceof(2,n,:,:), resuv)
       pcpc1 = new((/1,dimsizes(eof_ts_cdata(n,:))/), "double")
       pcpc1(0,:) = (/eof_ts_cdata(n,:)/);PC2

    ; pcpc1 = new((/1,dimsizes(eof_ts_cdata(0,:))/), "double")
    ; pcpc1(0,:) = (/eof_ts_cdata(0,:)/);PC1

      plot20 = gsn_csm_xy(wks,year,pcpc1, res1)

     ; overlay(plot(n),plot2(n))
     ; draw(plot(n))
     draw(plot20)

     ;gsn_panel(wks,(/plot,plot20/),(/1,2/),resP)
     frame(wks)
     end do
    
 end