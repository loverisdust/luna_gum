load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/run_cor.ncl"
load "~/shapefile/shapefile_utils.ncl"
begin
;>============================================================<
;;;;;;;;;;;;;;;;;;;;;;;;;数据处理;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
;>------------------------------------------------------------<
  level = 200              ;选取层次
  season = "spring"
  yystr = 1961
  yyend = 2012
  nave = 9      ;nave点平滑
  sig = 0.05

  fa = addfile("~/shapefile/dixing/altitude_2.5x2.5_landonly.nc", "r")
  f1 = addfile("~/data/northeast/"+season+"/Pre_1961_2012_025x025_"+season+"_anom_northeast_new.nc", "r")  ;obs 
  f2 = addfile("~/data/northeast/"+season+"/pre.noaa.mon.anom_1961-2012_"+season+".nc", "r")     ;module
  f3 = addfile("~/data/northeast/"+season+"/uwnd.mon.anom_1961-2012_"+season+".nc", "r") 
  f4 = addfile("~/data/northeast/"+season+"/vwnd.mon.anom_1961-2012_"+season+".nc", "r") 
  f5 = addfile("~/data/northeast/"+season+"/hgt.mon.anom_1961-2012_"+season+".nc", "r") 
  f6 = addfile("~/data/northeast/"+season+"/sst.mon.anom_1961-2012_"+season+".nc", "r") 
  f7 = addfile("~/iobi.nc", "r")
  f8 = addfile("/data2/home/lurui/data/jra55/Q.nc", "r")

  ali = f7->iobi

  alt = fa->alt
  pre_anom_ave = f1->pre_anom(:,:,:)  ;读取数据   ;obs44
  pree_anom_ave = f2->pre_anom({yystr:yyend},:,:) ;          ;module44
  uwnd_anom_ave = f3->uwnd_anom({yystr:yyend},{200},:,:);44
  vwnd_anom_ave = f4->vwnd_anom({yystr:yyend},{200},:,:);44
  hgt_anom_ave = f5->hgt_anom({yystr:yyend},{200},:,:);44
  uwnd_anom_ave850 = f3->uwnd_anom({yystr:yyend},{850},:,:);44
  vwnd_anom_ave850 = f4->vwnd_anom({yystr:yyend},{850},:,:);44
  hgt_anom_ave850 = f5->hgt_anom({yystr:yyend},{850},:,:);44
  sst_anom_ave = f6->sst_anom({yystr:yyend},:,:) ; 44
  q = f8->q

  lat = f3->lat
  lon = f3->lon
  ;>------------------------------------------------------------<
  ;                  获取各维大小
  ;>------------------------------------------------------------<
  dims1 = dimsizes(pre_anom_ave)
  dims2 = dimsizes(pree_anom_ave)
  dims3 = dimsizes(uwnd_anom_ave)
  dims4 = dimsizes(sst_anom_ave)
  year = ispan(1961,1961+dims1(0)-1,1)     ; 要分析的年份
  ;>------------------------------------------------------------<
  ;;                  计算权重 
  ;>------------------------------------------------------------<
  w  = sqrt(cos(0.01745329 * pre_anom_ave&lat) )     ; 权重                
  wp = pre_anom_ave * conform(pre_anom_ave, w, 1)   ; wp
  copy_VarCoords(pre_anom_ave, wp)
  ;>------------------------------------------------------------<
  ;;                  EOF分解
  ;>------------------------------------------------------------<
  x    = wp({lat|:},{lon|:},{years|:})  ;时间维放在最右边
  neof = 10                             ; 保留前10个模态
  eof  = eofunc_Wrap(x, neof, False)    ; EOF分解
  ; eof = -eof
  tzh = eof@pcvar
  lam = tzh * sqrt(2. / dims1(0))       ;计算解释方差误差范围
  ;>------------------------------------------------------------<
  ;;                  north检验
  ;>------------------------------------------------------------<
  ntim = dims1(0)               
  sig_pcv = eofunc_north(eof@pcvar, ntim, False)
  ;>------------------------------------------------------------<
  ;;                 计算时间序列并标准化
  ;>------------------------------------------------------------<
  eof_ts = eofunc_ts_Wrap (x, eof, False)
  eof_ts_stddev = dim_stddev(eof_ts)       ;计算时间系数标准差
  eof_ts = dim_standardize_n(eof_ts,1,1)
  eof_ts_dev = dim_stddev_n(eof_ts, 1)
  ;>------------------------------------------------------------<
  ;                 将标准差与空间模态相乘
  ;>------------------------------------------------------------<
  do i = 0,9
    eof(i,:,:) = (/ eof(i,:,:) /) * (/ eof_ts_stddev(i) /)  
  end do

  ; eof_ts1 = dtrend_n(eof_ts, False, 1)
  ; copy_VarCoords(eof_ts, eof_ts1)
  ; eof_ts := eof_ts1
  ;>------------------------------------------------------------<
  ;;              计算PC的9年滑动平均
  ;>------------------------------------------------------------<
  pc_mean = runave_n_Wrap(eof_ts, nave, 0, 1)

  do m = 0, (nave/2 - 1), 1
      pc_mean(:,m) = dim_avg_n_Wrap(eof_ts(:,0:nave/2+m), 1)
      pc_mean(:,dims1(0)-1-m) = dim_avg_n_Wrap(eof_ts(:,dims1(0)-nave/2 -1 -m:), 1)
  end do
  ;>------------------------------------------------------------<
  ;;                    分离年际尺度
  ;>------------------------------------------------------------<
  pc_year = eof_ts - pc_mean
  copy_VarCoords(eof_ts, pc_year)
  copy_VarCoords(eof_ts, pc_mean)

  ; m1 = dim_variance_n(eof_ts(0,:), 0)
  ; m2 = dim_variance_n(pc_mean(0,:), 0)
  ; m3 = dim_variance_n(pc_year(0,:), 0)

  ; print(m1-m2-m3)

  ; print(m2/m1)
  ; print(m3/m1)
  ; exit()

  ; f1_out = addfile("~/pc_year.nc", "c")
  ; f1_out->pcy = pc_year
  ; f2_out = addfile("~/pc_dacade.nc", "c")
  ; f2_out->pcd = pc_mean
  ; f3_out = addfile("~/pc1.nc", "c")
  ; f3_out->pc = eof_ts
  ;>------------------------------------------------------------<
  ali = dim_standardize_n(ali,1,0)
  ali_mean = runave_n_Wrap(ali, nave, 0, 0)

  do m = 0, (nave/2 - 1), 1
      ali_mean(m) = dim_avg_n_Wrap(ali(0:nave/2+m), 0)
      ali_mean(dims1(0)-1-m) = dim_avg_n_Wrap(ali(dims1(0)-nave/2 -1 -m:), 0)
  end do


  pc_year(0,:) = ali -ali_mean
  pc_mean(0,:) = (/ali_mean/)

  ;>------------------------------------------------------------<
  ;;                       回归
  ;>------------------------------------------------------------<
  pc = pc_mean(0,{yystr:yyend})
  ; print(ind(pc .gt. dim_stddev_n_Wrap(pc, 0)))
  ; print(ind(pc .lt. -dim_stddev_n_Wrap(pc, 0)))
  ; exit()

    ; pre_recon1 = eof_ts(0,:) * avg(eof(0,:,:)) 
    ; fffff = addfile("~/pre_recon_pc1.nc", "c")
    ; fffff->pre_rec_pc1 = pre_recon1
    ; exit()

  ; nao = asciiread("~/nao_index.txt", (/624,3/), "float")
  ; nao_index = nao(:,2)
  ; nao_index!0 = "time"
  ; nao_index&time = ispan(1, 624, 1)
  ; nao_spring = dim_standardize(month_to_season(nao_index, "MAM"), 1)

  ; nao_mean = runave_n_Wrap(nao_spring, nave, 0, 0)

  ; do m = 0, (nave/2 - 1), 1
  ;     nao_mean(m) = dim_avg_n_Wrap(nao_spring(0:nave/2+m), 0)
  ;     nao_mean(dims1(0)-1-m) = dim_avg_n_Wrap(nao_spring(dims1(0)-nave/2 -1 -m:), 0)
  ; end do


  ; nao_year = nao_spring -nao_mean
  ; nao_mean = nao_mean

  ; ; print(esccr(pdo_spring, pc, 0))
  ; ; print(esccr(nat_spring, pc, 0))
  ; print(esccr(nao_year, pc, 0))
  ; print(rtest(esccr(nao_year, pc, 0), 50, 0))

  ; print(nao_index)
  ; exit()


  ; oooo = regCoef(pc, pre_anom_ave({lat|:},{lon|:},{years|:}))
  ; copy_VarCoords(pre_anom_ave(0,:,:), oooo)


  pree_y = regCoef(pc, pree_anom_ave({lat|:},{lon|:},{years|:}))
  uu_y = regCoef(pc, uwnd_anom_ave({lat|:},{lon|:},{years|:}))
  vv_y = regCoef(pc, vwnd_anom_ave({lat|:},{lon|:},{years|:}))
  hh_y = regCoef(pc, hgt_anom_ave({lat|:},{lon|:},{years|:}))  
  uu_y850 = regCoef(pc, uwnd_anom_ave850({lat|:},{lon|:},{years|:}))
  vv_y850 = regCoef(pc, vwnd_anom_ave850({lat|:},{lon|:},{years|:}))
  hh_y850 = regCoef(pc, hgt_anom_ave850({lat|:},{lon|:},{years|:}))  
  ss_y = regCoef(pc, sst_anom_ave({lat|:},{lon|:},{years|:}))


;;;;;;;;;;q
  ; q_mean = month_to_season(q, "MAM")
  ; q_anom = dim_rmvmean_n_Wrap(q_mean, 0)
  ; delete(q_mean)
  ; qq_y = regCoef(pc, q_anom({level|:},{lat|:},{lon|:},{time|:}))
  ; copy_VarMeta(q(0,:,:,:), qq_y)
  ; delete(q)

  ; fff = addfile("~/q_reg_to_pc_id.nc", "c")
  ; fff->q = qq_y
  ; exit()

  ; copy_VarCoords(uwnd_anom_ave(0,:,:), uu_y)
  ; copy_VarCoords(uwnd_anom_ave(0,:,:), vv_y)
  ; printVarSummary(uu_y)

  ; f1 = addfile("~/uv200_reg_to_pc1_ia.nc", "c")
  ; f1->u = uu_y
  ; f1->v = vv_y

  ; exit()
;;;;;;;;;;;;

  ; copy_VarCoords(hgt_anom_ave(0,:,:), hh_y)
  ; copy_VarCoords(uwnd_anom_ave(0,:,:), uu_y)
  ; copy_VarCoords(uwnd_anom_ave(0,:,:), vv_y)
  ; printVarSummary(hh_y)
  ; fout = addfile("~/hgt_uv_200_reg_to_mci_id.nc", "c")
  ; fout->hgt = hh_y
  ; fout->uwnd = uu_y
  ; fout->vwnd = vv_y
  ; exit()

  ;;;;;edof
  rpc = esacr(pc, 1)
  rpre = esacr_n(pree_anom_ave, 1, 0)
  rsst = esacr_n(sst_anom_ave, 1, 0)
  ruu = esacr_n(uwnd_anom_ave, 1, 0)
  rvv = esacr_n(vwnd_anom_ave, 1, 0)
  ru850 = esacr_n(uwnd_anom_ave850, 1, 0)
  rv850 = esacr_n(uwnd_anom_ave850, 1, 0)

  edof_pre = dims2(0)*(1 - rpc(1)*rpre(1,:,:))/(1 + rpc(1)*rpre(1,:,:))
  edof_sst = dims4(0)*(1 - rpc(1)*rsst(1,:,:))/(1 + rpc(1)*rsst(1,:,:))
  edof_uu = dims3(0)*(1 - rpc(1)*ruu(1,:,:))/(1 + rpc(1)*ruu(1,:,:))
  edof_vv = dims3(0)*(1 - rpc(1)*rvv(1,:,:))/(1 + rpc(1)*rvv(1,:,:))
  edof_u850 = dims3(0)*(1 - rpc(1)*ru850(1,:,:))/(1 + rpc(1)*ru850(1,:,:))
  edof_v850 = dims3(0)*(1 - rpc(1)*rv850(1,:,:))/(1 + rpc(1)*rv850(1,:,:))

  ;>------------------------------------------------------------<
  ;;                    t-test
  ;>------------------------------------------------------------<
  rc_pre_y0 = pree_anom_ave(0,:,:)
  rc_pre_y0 = tofloat(pree_y)
  confi_rc_pre_y0 = rc_pre_y0
  p_y0 = student_t(pree_y@tval, ndtooned(edof_pre))    
  confi_rc_pre_y0 = onedtond( tofloat(p_y0), (/dims2(1),dims2(2)/))
  ;>------------------------------------------------------------<  
  rc_pre_y1 = uwnd_anom_ave(0,:,:)
  rc_pre_y1 = tofloat(uu_y)
  confi_rc_pre_y1 = rc_pre_y1
  p_y1 = student_t(uu_y@tval, ndtooned(edof_uu))
  confi_rc_pre_y1 = onedtond( tofloat(p_y1), (/dims3(1),dims3(2)/))
  ;>------------------------------------------------------------<
  rc_pre_y2 = uwnd_anom_ave(0,:,:)
  rc_pre_y2 = tofloat(vv_y)
  confi_rc_pre_y2 = rc_pre_y2
  p_y2 = student_t(vv_y@tval, ndtooned(edof_vv))
  confi_rc_pre_y2 = onedtond( tofloat(p_y2), (/dims3(1),dims3(2)/))
  ;>------------------------------------------------------------<
  rc_pre_y3 = hgt_anom_ave(0,:,:)
  rc_pre_y3 = tofloat(hh_y)
  confi_rc_pre_y3 = rc_pre_y3
  p_y3 = student_t(hh_y@tval, hh_y@nptxy-2)
  confi_rc_pre_y3 = onedtond( tofloat(p_y3), (/dims3(1),dims3(2)/))
  ;>------------------------------------------------------------<
  ;>------------------------------------------------------------<  
  rc_pre_y1850 = uwnd_anom_ave850(0,:,:)
  rc_pre_y1850 = tofloat(uu_y850)
  confi_rc_pre_y1850 = rc_pre_y1850
  p_y1850 = student_t(uu_y850@tval, ndtooned(edof_u850))
  confi_rc_pre_y1850 = onedtond( tofloat(p_y1850), (/dims3(1),dims3(2)/))
  ;>------------------------------------------------------------<
  rc_pre_y2850 = uwnd_anom_ave850(0,:,:)
  rc_pre_y2850 = tofloat(vv_y850)
  confi_rc_pre_y2850 = rc_pre_y2850
  p_y2850 = student_t(vv_y850@tval, ndtooned(edof_v850))
  confi_rc_pre_y2850 = onedtond( tofloat(p_y2850), (/dims3(1),dims3(2)/))
  ;>------------------------------------------------------------<
  rc_pre_y3850 = hgt_anom_ave850(0,:,:)
  rc_pre_y3850 = tofloat(hh_y850)
  confi_rc_pre_y3850 = rc_pre_y3850
  p_y3850 = student_t(hh_y850@tval, hh_y850@nptxy-2)
  confi_rc_pre_y3850 = onedtond( tofloat(p_y3850), (/dims3(1),dims3(2)/))
  ;>------------------------------------------------------------<
  rc_pre_y4 = sst_anom_ave(0,:,:)
  rc_pre_y4 = tofloat(ss_y)
  confi_rc_pre_y4 = rc_pre_y4
  ; p_y4 = student_t(ss_y@tval, where((ss_y@nptxy-2)/2.lt.0, 0, (ss_y@nptxy-2)/2))
  p_y4 = student_t(ss_y@tval, ndtooned(edof_sst))
  confi_rc_pre_y4 = onedtond( tofloat(p_y4), (/dims4(1),dims4(2)/))
  ;>------------------------------------------------------------<

  ; rc_pre_y1 = where(confi_rc_pre_y1 .gt. sig, rc_pre_y1@_FillValue ,rc_pre_y1)
  ; rc_pre_y2 = where(confi_rc_pre_y2 .gt. sig, rc_pre_y2@_FillValue ,rc_pre_y2)
  ; rc_pre_y1850 = where(confi_rc_pre_y1850 .gt. sig, rc_pre_y1850@_FillValue ,rc_pre_y1850)
  ; rc_pre_y2850 = where(confi_rc_pre_y2850 .gt. sig, rc_pre_y2850@_FillValue ,rc_pre_y2850)


;>============================================================<
 ;;;;;;;;;;;;;;;;;;;;;;;;绘图;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;>============================================================<
    ;wks = gsn_open_wks("eps","pre_anom_"+season+"_ave_eof")
    wks = gsn_open_wks("eps","pre_eof&pc")

    res                   = True    
    
    res@gsnDraw           = False         
    res@gsnFrame          = False        
    res@gsnLeftString     = ""
    res@gsnRightString    = ""
    res@gsnCenterString   = ""  
    res@gsnTickMarksOn = True ; 关闭经度标签

    res1 = res      ;时间系数曲线属性
    res2 = res      ;时间系数曲线滑动平均属性
    resc = res      ;创建resc，以绘制回归场的等值线       
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
  ;                  设置地图 
  ;>------------------------------------------------------------<
    res@mpShapeMode = "FreeAspect"
    res@vpXF = 0.1
    res@vpYF = 0.8
    res@vpHeightF = 0.3
    res@vpWidthF = 0.31

    res@mpFillOn = False
    res@mpMinLatF               = 38.                         
    res@mpMaxLatF               = 54.
    res@mpMinLonF               = 115.
    res@mpMaxLonF               = 135.5
    res@mpDataBaseVersion = "MediumRes"
    res@mpDataSetName = "Earth..4"
    res@mpOutlineBoundarySets   = "NoBoundaries"
    res@gsnLeftStringOrthogonalPosF = -0.005
    res@gsnRightStringOrthogonalPosF = -0.015
  ;>------------------------------------------------------------<
  ;                 在左标题处添加方差贡献
  ;>------------------------------------------------------------<
    x_variance = decimalPlaces(eof@pcvar,0,True)   

    resp1 = res
    ; resp2 = res
    ; resp3 = res
    ; resp4 = res
    resp1@gsnLeftString = "(a) EOF1"+" ("+x_variance(0)+"%"+")"
    resp1@gsnLeftStringFontHeightF = 0.015
    ; resp2@gsnLeftString = "EOF2"+"("+x_variance(1)+"%"+")"
    ; resp3@gsnLeftString = "EOF3"+"("+x_variance(2)+"%"+")"
    ; resp4@gsnLeftString = "EOF4"+"("+x_variance(3)+"%"+")"
  ;>------------------------------------------------------------<
    shpfn="~/shapefile/dongsansh_new/northeast_new.shp" 
    river ="~/shapefile/cnmap/river2_northeast.shp"
    plot0 = gsn_csm_contour_map(wks, eof(0,:,:), resp1)

    lres = True
    lres@gsLineThicknessF = 1.0
    line_mask1 = gsn_add_shapefile_polylines(wks, plot0, shpfn, lres)
    lres@gsLineThicknessF = 1.5
    lres@gsLineColor = "blue"
    line_mask2 = gsn_add_shapefile_polylines(wks, plot0, river, lres)

    ; draw(plot0)
    ; frame(wks)
  ;>------------------------------------------------------------<
  ;                    绘制置信度
  ;>------------------------------------------------------------<
    wks1y = gsn_open_wks("eps","reg_"+level+"hgt&uv&sst_ia")

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
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; res0@mpGridAndLimbOn = True
    ; res0@mpGridLatSpacingF = 5
    ; res0@mpGridLonSpacingF = 5
    ; res0@mpGridLineDashPattern = 2
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;>------------------------------------------------------------<
  ;                  位势高度回归场等值线设置
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
    ; resuv@vcGlyphStyle = "CurlyVector"   ;设置向量的显示形状（大小方向）
    resuv@vcGlyphStyle = "LineArrow"   ;设置向量的显示形状（大小方向）
    resuv@vcLineArrowThicknessF = 0.7;1.0
    resuv@vcMinDistanceF = 0.02   ;设置箭头间最小距离
    ; resuv@vcMinMagnitudeF = 0.5
    ; resuv@vcMaxMagnitudeF = 4.5
    resuv@vcRefAnnoOn = True
    resuv@vcRefLengthF = 0.025
    resuv@vcRefAnnoBackgroundColor = "white"
    resuv@vcRefAnnoPerimOn = True   ;关闭边框
    resuv@vcRefAnnoFontHeightF = 0.008
    resuv@vcRefAnnoString1On =   False ;设置标注箭头上方不显示注释
    resuv@vcRefAnnoString2On = True   ;设置标注箭头下方显示注释
    ; resuv@vcRefAnnoOrthogonalPosF = -0.18  ;垂直移动位置
    resuv@vcRefAnnoOrthogonalPosF = -0.12  ;垂直移动位置
  ;>------------------------------------------------------------<
  ;                  绘制box
  ;>------------------------------------------------------------<
    gonres = True
    gonres@gsLineThicknessF = 3.0
    gonres@gsLineColor = "red"
    gonres@gsLineDashPattern = 2

    ; lat0s = 25   ;ALI_2
    ; lat0e = 45
    ; lon0s = 135
    ; lon0e = 165

    lat1s = 0   ;ALI_2  1
    lat1e = 24
    lon1s = 285
    lon1e = 350

    ; lat0s = 45   ;ALI_2
    ; lat0e = 60
    ; lon0s = 290
    ; lon0e = 340

    ; lat1s = 5   ;ALI_2
    ; lat1e = 20
    ; lon1s = 290
    ; lon1e = 340

    ; lat2s = 25   ;ALI_2
    ; lat2e = 40
    ; lon2s = 290
    ; lon2e = 340

    lat3s = -15   ;IOBI
    lat3e = 20
    lon3s = 90;120
    lon3e = 160    

    ; lat3s = 20
    ; lat3e = 45
    ; lon3s = 170
    ; lon3e = 210

    ; x0 = (/lon0s,lon0e,lon0e,lon0s,lon0s/)
    ; y0 = (/lat0e,lat0e,lat0s,lat0s,lat0e/)
    x1 = (/lon1s,lon1e,lon1e,lon1s,lon1s/)
    y1 = (/lat1e,lat1e,lat1s,lat1s,lat1e/)
    ; x2 = (/lon2s,lon2e,lon2e,lon2s,lon2s/)
    ; y2 = (/lat2e,lat2e,lat2s,lat2s,lat2e/)
    x3 = (/lon3s,lon3e,lon3e,lon3s,lon3s/)
    y3 = (/lat3e,lat3e,lat3s,lat3s,lat3e/)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;run-cor;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; print(esccr(pc, x1x2_mean, 0))
; print(esccr(pc, x3_mean, 0))
; pcccc1 = esacr(pc, 1)
; pcccc2 = esacr(x3_mean, 1)
; eee = dims2(0)*(1 - pcccc1(1)*pcccc2(1))/(1 + pcccc1(1)*pcccc2(1))
; print(eee)

;     year@units = "years since 1961-1"
;     atl = run_cor(pc, x3_mean, year, 11)
;     aaa = run_cor(pc, x1x2_mean, year, 11)
;     pccc = new((/2,dimsizes(atl)/), "double")
;     pccc(0,:) = (/atl/)
;     pccc(1,:) = (/aaa/)

;     resr = True
;     resr@gsnDraw = False
;     resr@gsnFrame = False
;     resr@vpHeightF = 0.35
;     resr@vpWidthF = 0.6

;     resr@trXMinF = yystr
;     resr@trXMaxF = yyend
;     resr@trYMinF = -1
;     resr@trYMaxF = 1

;     resr@gsnYRefLine = 0.413
;     resr@gsnYRefLineDashPattern = 12
;     resr@gsnYRefLineThicknessF = 0.5
;     resr@gsnXRefLineColor = "red"

;     resr@tmXMajorGrid = True
;     resr@tmXMajorGridThicknessF = 0.3
;     resr@tmXMajorGridLineDashPattern = 2
;     resr@tmYMajorGrid = True
;     resr@tmYMajorGridThicknessF = 0.3
;     resr@tmYMajorGridLineDashPattern = 2

;     resr@tiXAxisString = "Year"   ;设置标题
;     resr@tmXTOn = False
;     resr@tmYROn = False

;     resr@pmLegendDisplayMode = "Always"
;     resr@pmLegendWidthF = 0.1
;     resr@pmLegendHeightF = 0.05
;     resr@pmLegendOrthogonalPosF = -1.22
;     resr@pmLegendParallelPosF = 0.1
;     resr@lgLabelFontHeightF = 0.01
;     ;resr@lgItemOrder = (/2,1,0/)
;     resr@lgBoxMinorExtentF = 0.4
;     resr@xyExplicitLabels = (/"AT","WNP"/)
;     resr@xyLineColors = (/"blue","red"/)
;     resr@xyLineThicknesses = (/3,3/)
;     resr@xyDashPatterns = (/0,0/)

;     wkss = gsn_open_wks("eps", "runcor")
;     plottt = gsn_csm_xy(wkss, year(5:46), pccc, resr)
;     draw(plottt)
;     frame(wkss)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;>------------------------------------------------------------<
  ;                        绘图
  ;>------------------------------------------------------------<
    plot1y = new(3, graphic)
    plot2y = new(3, graphic)
    plot1d = new(3, graphic)
    shpfn="~/shapefile/dongsansh_new/northeast_new_whole.shp" 
    lres = True
    lres@gsLineColor = "black"
    lres@gsLineThicknessF = 1.5
    ts = True
    ts@txFontHeightF = 0.03;0.04
    ts@txFontThicknessF = 0.7

    ; if (level .eq. 850) then
      resp@cnFillPalette = "BlueWhiteOrangeRed"

      ; resc@cnLevels := (/-4,-1,1,4/);回归到IA 
      ; resuv@vcRefMagnitudeF = 0.3;回归到IA
      ; resuv@vcRefAnnoString2 = "0.3"
      ; ; resp@cnLevels := (/-0.3,-0.2,-0.1,-0.05,0,0.05,0.1,0.2,0.3/);回归到IA 
      ; resp@cnLevels := (/-0.4,-0.2,-0.1,-0.05,0,0.05,0.1,0.2,0.4/)
      ; res0@gsnLeftString = "(b) Reg_850hgt_uv_sst_to_IA"

      ; resc@cnLevels := (/-20,-8,-5,-0.1,5,15/);回归到ALI 
      ; resuv@vcRefMagnitudeF = 0.7;回归到ALI
      ; resuv@vcRefAnnoString2 = "0.7"
      ; ; resp@cnLevels := (/-0.7,-0.4,-0.2,-0.1,0,0.1,0.2,0.4,0.7/);回归到ALI 
      ; resp@cnLevels := (/-0.4,-0.2,-0.1,-0.05,0,0.05,0.1,0.2,0.4/)
      ; res0@gsnLeftString = "(c) Reg_850hgt_uv_sst_to_ALI"

      resc@cnLevels := (/-0.1,10/);回归到MCI 
      resuv@vcRefMagnitudeF = 1.0;回归到MCI
      resuv@vcRefAnnoString2 = "1.0"
      resp@cnLevels := (/-0.9,-0.6,-0.3,-0.1,0,0.1,0.3,0.6,0.9/);回归到MCI 
      res0@gsnLeftString = "(b) Reg_850hgt_uv_sst_to_MCI"

      ; resc@cnLevels := (/-5,-0.1,10,15/);回归到ID      
      ; resuv@vcRefMagnitudeF = 1.5;回归到ID
      ; resuv@vcRefAnnoString2 = "1.5"
      ; resuv@vcMinMagnitudeF :=  0.5 ;0.1
      ; ; resp@cnLevels := (/-1.6,-0.8,-0.4,-0.2,0,0.2,0.4,0.8,1.6/);回归到ID
      ; resp@cnLevels := (/-0.9,-0.6,-0.3,-0.1,0,0.1,0.3,0.6,0.9/)
      ; res0@gsnLeftString = "(b) Reg_850hgt_uv_sst_to_ID"

      base = gsn_csm_contour_map(wks1y, confi_rc_pre_y4, res0)
      plot_altp = gsn_csm_contour(wks1y, alt, resalt)
      line_mask1y = gsn_add_shapefile_polylines(wks1y, base, shpfn, lres)
      plot1y(0) = gsn_csm_contour(wks1y,rc_pre_y4,resp)   ;sst
      plot1y(1) = gsn_csm_vector(wks1y, rc_pre_y1850, rc_pre_y2850, resuv)    ;uv
      plot1y(2) = gsn_csm_contour(wks1y,rc_pre_y3850,resc)    ;hgt
      plot1y(2) = ColorNegDashZeroPosContour(plot1y(2),"blue","black","red")
      ; dum0y = gsn_add_polyline(wks1y, plot1y(2), x0, y0, gonres)
      ; dum1y = gsn_add_polyline(wks1y, plot1y(2), x1, y1, gonres)
      ; dum2y = gsn_add_polyline(wks1y, plot1y(2), x2, y2, gonres)
      dum3y = gsn_add_polyline(wks1y, plot1y(2), x3, y3, gonres)

      ts@txFontColor = "blue"
      ; text1 = gsn_add_text(wks1y, plot1y(2), "A", 325, 60, ts);回归到IA
      ; text2 = gsn_add_text(wks1y, plot1y(2), "A", 74, 52, ts)
      ; text6 = gsn_add_text(wks1y, plot1y(2), "A", 15, 36, ts)

      ; text1 = gsn_add_text(wks1y, plot1y(2), "A", 323, 60, ts);回归到ALI
      ; text2 = gsn_add_text(wks1y, plot1y(2), "A", 69, 52, ts)
      ; text7 = gsn_add_text(wks1y, plot1y(2), "A", 5, 30, ts)

      text3 = gsn_add_text(wks1y, plot1y(2), "A", 135, 14, ts);回归到ID_new
      text2 = gsn_add_text(wks1y, plot1y(2), "A", 85, 42, ts)

      ts@txFontColor = "red"
      ; text3 = gsn_add_text(wks1y, plot1y(2), "C", 32, 66, ts);回归到IA
      ; text4 = gsn_add_text(wks1y, plot1y(2), "C", 124, 45, ts)
      ; text5 = gsn_add_text(wks1y, plot1y(2), "C", 318, 36, ts)
      ; text7 = gsn_add_text(wks1y, plot1y(2), "C", 40, 20, ts)

      ; text3 = gsn_add_text(wks1y, plot1y(2), "C", 30, 50, ts);回归到ALI
      ; text4 = gsn_add_text(wks1y, plot1y(2), "C", 124, 45, ts)
      ; text5 = gsn_add_text(wks1y, plot1y(2), "C", 320, 37, ts)
      ; text4 = gsn_add_text(wks1y, plot1y(2), "C", 40, 22, ts)

      text4 = gsn_add_text(wks1y, plot1y(2), "C", 145, 35, ts);回归到ID_new


      ts@txFontColor = "black"
      ts@txFontHeightF = 0.04;0.05
      ; text9 = gsn_add_text(wks1y, plot1y(2), "S1", 318, 12, ts)
      text10 = gsn_add_text(wks1y, plot1y(2), "S2", 125, 3, ts)

      overlay(base, plot_altp)
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

      ; resc@cnLevels := (/-15,-7.5,-0.1,7.5,10/);;回归到IA
      ; resuv@vcRefMagnitudeF := 1.0;回归到IA
      ; resuv@vcRefAnnoString2 := "1.0"
      ; ; resp@cnLevels := (/-0.6,-0.3,-0.1,-0.02,0,0.02,0.1,0.3,0.6/);回归到IA
      ; resp@cnLevels := (/-0.8,-0.3,-0.1,-0.02,0,0.02,0.1,0.3,0.8/)
      ; res0@gsnLeftString = "(a) Reg_200hgt_uv_pre_to_IA"

      ; resc@cnLevels := (/-10,-1,10,25/);;回归到ALI
      ; resuv@vcRefMagnitudeF := 1.5;回归到ALI MCI
      ; resuv@vcRefAnnoString2 := "1.5"
      ; ; resp@cnLevels := (/-1.6,-1.0,-0.5,-0.05,0,0.05,0.5,1.0,1.6/);回归到ALI
      ; resp@cnLevels := (/-0.8,-0.3,-0.1,-0.02,0,0.02,0.1,0.3,0.8/)
      ; res0@gsnLeftString = "(b) Reg_200hgt_uv_pre_to_ALI"

      resc@cnLevels := (/-5,-1,20,30/);;回归到MCI
      resuv@vcRefMagnitudeF := 1.5;回归到ALI MCI
      resuv@vcRefAnnoString2 := "1.5"
      resp@cnLevels := (/-1.2,-0.8,-0.4,-0.05,0,0.05,0.4,0.8,1.2/);回归到MCI
      res0@gsnLeftString = "(a) Reg_200hgt_uv_pre_to_MCI"

      ; resc@cnLevels := (/-1,30,42/);;回归到ID
      ; resuv@vcMinMagnitudeF = 0.5
      ; resuv@vcMaxMagnitudeF = 2.0
      ; resuv@vcRefMagnitudeF = 3.0;回归到ID
      ; resuv@vcRefAnnoString2 = "3.0"
      ; resuv@vcMinMagnitudeF =  1.5;0.2
      ; ; resp@cnLevels := (/-2.4,-1.6,-0.8,-0.1,0,0.1,0.8,1.6,2.4/);回归到ID
      ; resp@cnLevels := (/-1.2,-0.8,-0.4,-0.05,0,0.05,0.4,0.8,1.2/)
      ; res0@gsnLeftString = "(a) Reg_200hgt_uv_pre_to_ID"

      resp@cnFillColors := (/100,94,86,78,51,51,42,38,34,30/)
      base1 = gsn_csm_contour_map(wks1y, confi_rc_pre_y0, res0)
      line_mask2y = gsn_add_shapefile_polylines(wks1y, base1, shpfn, lres)
      plot2y(0) = gsn_csm_contour(wks1y,rc_pre_y0,resp)   ;降水
      plot2y(1) = gsn_csm_vector(wks1y, rc_pre_y1, rc_pre_y2, resuv)    ;uv
      plot2y(2) = gsn_csm_contour(wks1y,rc_pre_y3,resc)    ;hgt
      plot2y(2) = ColorNegDashZeroPosContour(plot2y(2),"blue","black","red")

      ts@txFontColor = "blue"
      ; text1 = gsn_add_text(wks1y, plot2y(2), "A", 317, 59, ts);回归到IA
      ; text2 = gsn_add_text(wks1y, plot2y(2), "A", 74, 51, ts)
      ; text6 = gsn_add_text(wks1y, plot2y(2), "A", 15, 36, ts)
      ; text8 = gsn_add_text(wks1y, plot2y(2), "A", 305, 20, ts)
      ; text1 := gsn_add_text(wks1y, plot2y(2), "A", 315, 61, ts);回归到ALI
      ; text2 := gsn_add_text(wks1y, plot2y(2), "A", 68, 44, ts)
      ; text6 := gsn_add_text(wks1y, plot2y(2), "A", -5, 30, ts)
      ; text8 := gsn_add_text(wks1y, plot2y(2), "A", 305, 20, ts)

      text3 = gsn_add_text(wks1y, plot2y(2), "A", 130, 15, ts);回归到ID_new
      text2 = gsn_add_text(wks1y, plot2y(2), "A", 67, 30, ts)


      ts@txFontColor = "red"
      ; text3 = gsn_add_text(wks1y, plot2y(2), "C", 15, 73, ts);回归到IA
      ; text4 = gsn_add_text(wks1y, plot2y(2), "C", 114, 44, ts)
      ; text5 = gsn_add_text(wks1y, plot2y(2), "C", 315, 38, ts)
      ; text7 = gsn_add_text(wks1y, plot2y(2), "C", 36, 23, ts)

      ; text3 := gsn_add_text(wks1y, plot2y(2), "C", 15, 60, ts);回归到ALI
      ; text4 := gsn_add_text(wks1y, plot2y(2), "C", 114, 44, ts)
      ; text5 := gsn_add_text(wks1y, plot2y(2), "C", 316, 39, ts)
      ; text7 := gsn_add_text(wks1y, plot2y(2), "C", 32, 21, ts)

      text7 = gsn_add_text(wks1y, plot2y(2), "C", 122, 43, ts);回归到ID_new

      overlay(base1, plot2y(0))
      overlay(base1, plot2y(1))
      overlay(base1, plot2y(2))
      ; draw(base1)
      ; frame(wks1y)
    ; end if
    ; end if
    ; gsn_panel(wks1y, (/base1, base/), (/2,1/), False)

  ;>------------------------------------------------------------<
  ;                  绘制时间序列
  ;>------------------------------------------------------------<
    res1@vpXF = 0.45
    res1@vpYF = 0.8
    res1@vpHeightF = 0.3
    res1@vpWidthF = 0.55

    res1@trXMinF = yystr
    res1@trXMaxF = yyend
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

    ; pcpc1 = new((/3,dimsizes(eof_ts(0,:))/), "double")
    ; pcpc1(0,:) = (/eof_ts(0,:)/)
    ; pcpc1(1,:) = (/pc_mean(0,:)/)
    ; pcpc1(2,:) = (/pc_year(0,:)/)
    ; res1@pmLegendDisplayMode = "Always"
    ; res1@pmLegendWidthF = 0.07
    ; res1@pmLegendHeightF = 0.05
    ; res1@pmLegendOrthogonalPosF = -1.12 ;-1.22
    ; res1@pmLegendParallelPosF = 0.08
    ; res1@lgLabelFontHeightF = 0.01
    ; res1@lgItemOrder = (/2,1,0/)
    ; res1@lgBoxMinorExtentF = 0.4
    ; res1@xyExplicitLabels = (/" PC1","  ID","  IA"/)
    ; res1@xyLineColors = (/"black","blue","red"/)
    ; res1@xyLineThicknesses = (/2,3,3/)
    ; res1@xyDashPatterns = (/0,2,2/)

    ; res1@gsnLeftStringOrthogonalPosF = -0.005
    ; res1@gsnRightStringOrthogonalPosF = -0.015
    ; res1@gsnLeftString = "(b) Principle Component"
    ; res1@gsnLeftStringFontHeightF = 0.015    

    pre_toto = dim_avg_n_Wrap(pre_anom_ave, (/1,2/))
    pre_toto = dim_standardize(pre_toto, 1)
    pcpc1 = new((/3,dimsizes(eof_ts(0,:))/), "double")
    pcpc1(0,:) = (/eof_ts(0,:)/)
    pcpc1(1,:) = (/pc_mean(0,:)/)
    pcpc1(2,:) = (/pc_year(0,:)/)
    ; pcpc1(3,:) = (/pre_toto/)
    res1@pmLegendDisplayMode = "Always"
    res1@pmLegendWidthF = 0.07
    res1@pmLegendHeightF = 0.05
    res1@pmLegendOrthogonalPosF = -1.12 ;-1.22
    res1@pmLegendParallelPosF = 0.08
    res1@lgLabelFontHeightF = 0.01
    res1@lgItemOrder = (/2,1,0/)
    res1@lgBoxMinorExtentF = 0.4
    res1@xyExplicitLabels = (/" PC1","  ID","  IA"/);," Pre_total"
    res1@xyLineColors = (/"black","blue","red"/);,"gray"
    res1@xyLineThicknesses = (/2,3,3/);,2
    res1@xyDashPatterns = (/0,2,2/);,0

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
    draw(plot0)
    draw(plot20)
    frame(wks)
  ;>------------------------------------------------------------<
  ;                 绘制box与pc的时间序列
  ;>------------------------------------------------------------<
    ; pcpc = new((/2,dimsizes(eof_ts(0,:))/), "double")
    ; pcpc(0,:) = (/pc_year(0,:)/)
    ; pcpc(1,:) = (/pppp1/)
    ; ; pcpc2 = new((/2,dimsizes(eof_ts(0,:))/), "double")
    ; ; pcpc2(0,:) = (/pc_mean(0,:)/)
    ; ; pcpc2(1,:) = (/x1x2_mean/)
    ; ; pcpc3 = new((/2,dimsizes(eof_ts(0,:))/), "double")
    ; ; pcpc3(0,:) = (/pc_mean(0,:)/)
    ; ; pcpc3(1,:) = (/x3_mean/)

    ; wks21 = gsn_open_wks("eps", "ALI_1")
    ; ; wks22 = gsn_open_wks("eps", "ALI_2") 
    ; ; wks23 = gsn_open_wks("eps", "IOBI")   

    ; res2@vpHeightF = 0.35
    ; res2@vpWidthF = 0.6
    ; res2@trXMinF = yystr
    ; res2@trXMaxF = yyend
    ; res2@trYMinF = -3.0
    ; res2@trYMaxF = 3.0
    ; res2@gsnYRefLine = 0     ;设定参考值

    ; res2@tmXMajorGrid = True
    ; res2@tmXMajorGridThicknessF = 0.3
    ; res2@tmXMajorGridLineDashPattern = 2
    ; res2@tmYMajorGrid = True
    ; res2@tmYMajorGridThicknessF = 0.3
    ; res2@tmYMajorGridLineDashPattern = 2

    ; res2@tmXBMajorLengthF= 0.005
    ; res2@tmYLMajorLengthF= 0.005
    ; res2@tmYLLabelDeltaF = -0.7
    ; res2@tmXBLabelDeltaF = -0.7
    ; res2@tmXBLabelFontHeightF = 0.015
    ; res2@tmYLLabelFontHeightF = 0.015

    ; res2@tmYLMode = "Manual"
    ; res2@tmYLTickSpacingF = 1.0
    ; res2@tmYLMinorPerMajor = 1

    ; res2@tiXAxisString = "Year"   ;设置标题
    ; res2@tiXAxisFontHeightF = 0.015
    ; res2@tmXTOn = False
    ; res2@tmYROn = False

    ; res2@lgBoxMinorExtentF = 0.4
    ; res2@xyLineThicknesses := (/2,3/)
    ; res2@xyDashPatterns := (/0,2/)
    ; res2@lgItemOrder := (/1,0/)
    ; res2@xyLineColors := (/"black","red"/)
    ; res2@pmLegendDisplayMode = "Always"
    ; res2@pmLegendWidthF = 0.1
    ; res2@pmLegendHeightF = 0.05
    ; res2@pmLegendOrthogonalPosF = -1.15 ;-1.25
    ; res2@pmLegendParallelPosF = 0.09
    ; res2@lgLabelFontHeightF = 0.01

    ; res2@gsnLeftStringFontHeightF = 15;17
    ; res2@gsnLeftStringOrthogonalPosF = 0

    ; ; res2@gsnLeftString = "TCC="+ttt1+"(99%)"
    ; res2@gsnLeftString = "(a) Time series of PC1_IA and ALI"
    ; res2@xyExplicitLabels := (/" PC1_IA"," ALI"/)
    ; plot21 = gsn_csm_xy(wks21,year,pcpc, res2)
    ; ; res2@gsnLeftString = "TCC=0.69(99%)";"+ttt2+"
    ; ; res2@xyExplicitLabels := (/"ID","ALI_2"/)
    ; ; plot22 = gsn_csm_xy(wks22,year,pcpc2, res2)
    ; ; res2@gsnLeftString = "TCC="+ttt3+"(99%)"
    ; ; res2@xyExplicitLabels := (/"PC1_ID","IOBI"/)
    ; ; plot23 = gsn_csm_xy(wks23,year,pcpc3, res2)

    ; draw(plot21)
    ; frame(wks21)
    ; ; draw(plot22)
    ; ; frame(wks22)
    ; draw(plot23)
    ; frame(wks23)

end