load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
; ==============================================================
begin
latS   =  7.
latN   =  52. 
lonL   =  97. 
lonR   =  142.
i=0
yystr = 1979
yyend = 2020
season="MAM"

f =addfile("~/data/precip.mon.mean.nc","r")
f2=addfile("~/data/uwnd.mon.mean.nc","r")
f3=addfile("~/data/vwnd.mon.mean.nc","r")
f4=addfile("~/data/hgt.mon.mean.nc", "r") 
f5=addfile("~/data/sst.mon.mean.nc", "r") 

;>------------------------------------------------------------<
;;              数据预处理 
;>------------------------------------------------------------<

precip=f->precip({yystr:yyend},:,:)
uwnd_anom_ave = f2->uwnd({yystr:yyend},{200},:,:)
vwnd_anom_ave = f3->vwnd({yystr:yyend},{200},:,:)
hgt_anom_ave = f4->hgt({yystr:yyend},{200},:,:)
uwnd_anom_ave850 = f2->uwnd({yystr:yyend},{850},:,:)
vwnd_anom_ave850 = f3->vwnd({yystr:yyend},{850},:,:)
hgt_anom_ave850 = f4->hgt({yystr:yyend},{850},:,:)
sst_anom_ave = f5->sst({yystr:yyend},:,:) 

pre=dim_standardize_n_Wrap(pre_anom_ave, 0, 0)
precip=month_to_season(pre, season)
printVarSummary(precip)
precipvartime=dim_variance_n_Wrap(precip,0)


uwnd_200=dim_standardize_n_Wrap(uwnd_anom_ave, 0, 0)
u_200=month_to_season(uwnd_200, season)
u_200vartime=dim_variance_n_Wrap(u_200,0)

vwnd_200=dim_standardize_n_Wrap(vwnd_anom_ave, 0, 0)
v_200=month_to_season(vwnd_200, season)
v_200vartime=dim_variance_n_Wrap(v_200,0)

hgt_200_1=dim_standardize_n_Wrap(hgt_anom_ave, 0, 0)
hgt_200=month_to_season(hgt_200_1, season)
hgt_200vartime=dim_variance_n_Wrap(hgt_200,0)

uwnd_850=dim_standardize_n_Wrap(uwnd_anom_ave850, 0, 0)
u_850=month_to_season(uwnd_850, season)
u_850vartime=dim_variance_n_Wrap(u_850,0)

vwnd_850=dim_standardize_n_Wrap(vwnd_anom_ave850, 0, 0)
v_850=month_to_season(vwnd_850, season)
v_850vartime=dim_variance_n_Wrap(v_850,0)

hgt_850_1=dim_standardize_n_Wrap(hgt_anom_ave850, 0, 0)
hgt_850=month_to_season(hgt_850_1, season)
hgt_850vartime=dim_variance_n_Wrap(hgt_850,0)

sst_anom=dim_standardize_n_Wrap(sst_anom_ave, 0, 0)
sst=month_to_season(sst_anom, season)
sstvartime=dim_variance_n_Wrap(sst,0)



dimw    = dimsizes(precip)
year = ispan(1979,1979+dimw(0)-1,1)

mtim    = dimw(0)
mlat    = dimw(1)
mlon    = dimw(2)

;>------------------------------------------------------------<
;;              计算各个数据的异常场
;>------------------------------------------------------------<
  precip_mean=dim_avg_n(precip, 0)
  do i=0,17
    do j=0,17
      precip_2=(precip(24,i,j)+precip(32,i,j))/2
    end do 
  end do
  precip_ano=precip_2-precip_mean
  precip_2vartime=dim_variance_n_Wrap(precip_2,0)


  u_200_mean=dim_avg_n(u_200, 0)
  do i=0,17
    do j=0,17
      u_200_2=(u_200(24,i,j)+u_200(32,i,j))/2
    end do 
  end do
  u_200_ano=u_200_2-u_200_mean
  u_200_2vartime=dim_variance_n_Wrap(u_200_2,0)

  v_200_mean=dim_avg_n(v_200, 0)
  do i=0,17
    do j=0,17
      v_200_2=(v_200(24,i,j)+v_200(32,i,j))/2
    end do 
  end do
  v_200_ano=v_200_2-v_200_mean
  v_200_2vartime=dim_variance_n_Wrap(v_200_2,0)

  hgt_200_mean=dim_avg_n(hgt_200, 0)
  do i=0,17
    do j=0,17
      hgt_200_2=(hgt_200(24,i,j)+hgt_200(32,i,j))/2
    end do 
  end do
  hgt_200_ano=hgt_200_2-hgt_200_mean
  hgt_200_2vartime=dim_variance_n_Wrap(hgt_200_2,0)

  u_850_mean=dim_avg_n(u_850, 0)
  do i=0,17
    do j=0,17
      u_850_2=(u_850(24,i,j)+u_850(32,i,j))/2
    end do 
  end do
  u_850_ano=u_850_2-u_850_mean
  u_850_2vartime=dim_variance_n_Wrap(u_850_2,0)

  v_850_mean=dim_avg_n(v_850, 0)
  do i=0,17
    do j=0,17
      v_850_2=(v_850(24,i,j)+v_850(32,i,j))/2
    end do 
  end do
  v_850_ano=v_850_2-v_850_mean
  v_850_2vartime=dim_variance_n_Wrap(v_850_2,0)

  hgt_850_mean=dim_avg_n(hgt_850, 0)
  do i=0,17
    do j=0,17
      hgt_850_2=(hgt_850(24,i,j)+hgt_850(32,i,j))/2
    end do 
  end do
  hgt_850_ano=hgt_850_2-hgt_850_mean
  hgt_850_2vartime=dim_variance_n_Wrap(hgt_850_2,0)

  sst_mean=dim_avg_n(sst, 0)
  do i=0,17
    do j=0,17
      sst_2=(sst(24,i,j)+sst(32,i,j))/2
    end do 
  end do
  sst_ano=sst_2-sst_mean
  sst_2vartime=dim_variance_n_Wrap(sst_2,0)
;>------------------------------------------------------------<
;;                    t_test
;>------------------------------------------------------------<
iflag=True
sigr=0.05

preciptmp=precip(lat|:,lon|:,time|:)
precip_2tmp=precip_2(lat|:,lon|:,time|:)

u_200tmp=u_200(lat|:,lon|:,time|:)
u_200_2tmp=u_200_2(lat|:,lon|:,time|:)

v_200tmp=v_200(lat|:,lon|:,time|:)
v_200_2tmp=v_200_2(lat|:,lon|:,time|:)

hgt_200tmp=hgt_200(lat|:,lon|:,time|:)
hgt_200_2tmp=hgt_200_2(lat|:,lon|:,time|:)

u_850tmp=u_850(lat|:,lon|:,time|:)
u_850_2tmp=u_850_2(lat|:,lon|:,time|:)

v_850tmp=v_850(lat|:,lon|:,time|:)
v_850_2tmp=v_850_2(lat|:,lon|:,time|:)

hgt_850tmp=hgt_850(lat|:,lon|:,time|:)
hgt_850_2tmp=hgt_850_2(lat|:,lon|:,time|:)

ssttmp=sst(lat|:,lon|:,time|:)
sst_2tmp=sst_2(lat|:,lon|:,time|:)

precipEqv=equiv_sample_size(preciptmp, siglvl, 0)
precip_2Eqv=equiv_sample_size(precip_2tmp, siglvl, 0)
precipN=wgt_areaave(precipEqv, wgty, 1., 0)
precip_2N=wgt_areaave(precip_2Eqv, wgty, 1., 0)

u_200Eqv=equiv_sample_size(u_200tmp, siglvl, 0)
u_200_2Eqv=equiv_sample_size(u_200_2tmp, siglvl, 0)
u_200N=wgt_areaave(u_200Eqv, wgty, 1., 0)
u_200_2N=wgt_areaave(u_200_2Eqv, wgty, 1., 0)

v_200Eqv=equiv_sample_size(v_200tmp, siglvl, 0)
v_200_2Eqv=equiv_sample_size(v_200_2tmp, siglvl, 0)
v_200N=wgt_areaave(v_200Eqv, wgty, 1., 0)
v_200_2N=wgt_areaave(v_200_2Eqv, wgty, 1., 0)

hgt_200Eqv=equiv_sample_size(hgt_200tmp, siglvl, 0)
hgt_200_2Eqv=equiv_sample_size(hgt_200_2tmp, siglvl, 0)
hgt_200N=wgt_areaave(hgt_200Eqv, wgty, 1., 0)
hgt_200_2N=wgt_areaave(hgt_200_2Eqv, wgty, 1., 0)

u_850Eqv=equiv_sample_size(u_850tmp, siglvl, 0)
u_850_2Eqv=equiv_sample_size(u_850_2tmp, siglvl, 0)
u_850N=wgt_areaave(u_850Eqv, wgty, 1., 0)
u_850_2N=wgt_areaave(u_850_2Eqv, wgty, 1., 0)

v_850Eqv=equiv_sample_size(v_850tmp, siglvl, 0)
v_850_2Eqv=equiv_sample_size(v_850_2tmp, siglvl, 0)
v_850N=wgt_areaave(v_850Eqv, wgty, 1., 0)
v_850_2N=wgt_areaave(v_850_2Eqv, wgty, 1., 0)

hgt_850Eqv=equiv_sample_size(hgt_850tmp, siglvl, 0)
hgt_850_2Eqv=equiv_sample_size(hgt_850_2tmp, siglvl, 0)
hgt_850N=wgt_areaave(hgt_850Eqv, wgty, 1., 0)
hgt_850_2=wgt_areaave(hgt_850_2Eqv, wgty, 1., 0)

sstEqv=equiv_sample_size(ssttmp, siglvl, 0)
sst_2Eqv=equiv_sample_size(sst_2tmp, siglvl, 0)
sstN=wgt_areaave(sstEqv, wgty, 1., 0)
sst_2=wgt_areaave(sst_2Eqv, wgty, 1., 0)

prob_precip = ttest(precip_2,precip_2vartime,precip_2N, precip_mean,precipvartime,precipN, iflag, False) 
prob_u_200 = ttest(u_200_2,u_200_2vartime,u_200_2N, u_200_mean,u_200vartime,u_200N, iflag, False) 
prob_v_200 = ttest(v_200_2,v_200_2vartime,v_200_2N, v_200_mean,v_200vartime,v_200N, iflag, False) 
prob_hgt_200 = ttest(hgt_200_2,hgt_200_2vartime,hgt_200_2N, hgt_200_mean,hgt_200vartime,hgt_200N, iflag, False) 
prob_u_850 = ttest(u_850_2,u_850_2vartime,u_850_2N, u_850_mean,u_850vartime,u_850N, iflag, False) 
prob_v_850 = ttest(v_850_2,v_850_2vartime,v_850_2N, v_850_mean,v_850vartime,v_850N, iflag, False) 
prob_hgt_850 = ttest(hgt_850_2,hgt_850_2vartime,hgt_850_2N, hgt_850_mean,hgt_850vartime,hgt_850N, iflag, False) 
prob_sst = ttest(sst_2,sst_2vartime,sst_2N, sst_mean,sstvartime,sstN, iflag, False) 
;>------------------------------------------------------------<
;;                    绘图
;>------------------------------------------------------------<
    wks1y = gsn_open_wks("pdf","~/data/850_1105")
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
    ts@txFontHeightF = 0.03;0.04
    ts@txFontThicknessF = 0.7
    
    resp@cnFillPalette = "BlueWhiteOrangeRed"

    resc@cnLevels := (/-0.1,10/);回归到MCI 
    resuv@vcRefMagnitudeF = 1.0;回归到MCI
    resuv@vcRefAnnoString2 = "1.0"
    resp@cnLevels := (/-0.9,-0.6,-0.3,-0.1,0,0.1,0.3,0.6,0.9/);回归到MCI 
    res0@gsnLeftString = "(b) 850hgt_uv_sst"

      base = gsn_csm_contour_map(wks1y, confi_rc_pre_y4, res0)


      plot1y(0) = gsn_csm_contour(wks1y,sst_ano,resp)   ;sst
      plot1y(1) = gsn_csm_vector(wks1y, u_850_ano, v_850_ano, resuv)    ;uv
      plot1y(2) = gsn_csm_contour(wks1y,hgt_850_ano,resc)    ;hgt
      plot1y(2) = ColorNegDashZeroPosContour(plot1y(2),"blue","black","red")


      overlay(base, plot1y(0))
      overlay(base, plot1y(1))
      overlay(base, plot1y(2))
      ; draw(base)
      ; frame(wks1y)
  ;>------------------------------------------------------------<
    ts@txFontHeightF := 0.03;0.04
    ts@txFontThicknessF := 0.7


    ; else if (level .eq. 200) then
      ; gsn_define_colormap(wks1y, "ViBlGrWhYeOrRe")
      resp@cnFillPalette := "ViBlGrWhYeOrRe"

      resc@cnLevels := (/-5,-1,20,30/);;回归到MCI
      resuv@vcRefMagnitudeF := 1.5;回归到ALI MCI
      resuv@vcRefAnnoString2 := "1.5"
      resp@cnLevels := (/-1.2,-0.8,-0.4,-0.05,0,0.05,0.4,0.8,1.2/);回归到MCI
      res0@gsnLeftString = "(a) Reg_200hgt_uv_pre_to_MCI"

      resp@cnFillColors := (/100,94,86,78,51,51,42,38,34,30/)
      base1 = gsn_csm_contour_map(wks1y, confi_rc_pre_y0, res0)
      plot2y(0) = gsn_csm_contour(wks1y,precip_ano,resp)   ;降水
      plot2y(1) = gsn_csm_vector(wks1y, u_200_ano, V_200_ano, resuv)    ;uv
      plot2y(2) = gsn_csm_contour(wks1y,hgt_200_ano,resc)    ;hgt
      plot2y(2) = ColorNegDashZeroPosContour(plot2y(2),"blue","black","red")

  

      overlay(base1, plot2y(0))
      overlay(base1, plot2y(1))
      overlay(base1, plot2y(2))
       draw(base1)
       frame(wks1y)
end
    ; gsn_panel(wks1y, (/base1, base/), (/2,1/), False)