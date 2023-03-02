 begin
 a = addfile("/disk1/lvcy/cru_ts4.04.1901.2019.tmp.dat.nc","r")
  tm= short2flt(a->tmp(684:1415,:,:));1958-2018 {0:90},{0:180}
   time1 = tm&time
   tmp=dtrend_msg_n(time1,tm, False, True, 0)
   copy_VarMeta(tm, tmp)
; --------------------处理sat数据-------------------------
   t1=cd_calendar(tmp&time,0)
   month=toint(t1(:,1))

   latS   =  -90
   latN   =  90.
   lonL   =  -180.
   lonR   =  180.
   yrStrt = 1958
   yrLast = 2018
   e1=ind(month.eq.6)
   e2=ind(month.eq.7)
   e3=ind(month.eq.8)
   tmp_june=tmp(e1,:,:)
   tmp_july=tmp(e2,:,:)
   tmp_aug=tmp(e3,:,:)
   tmp_60=(tmp_june+tmp_july+tmp_aug)/3
   copy_VarMeta(tmp_june, tmp_60 )
   rad=3.1415926/180
   lat1=tmp_60&lat
   clat=cos(lat1*rad)

   ;---------------------------去除九年滑动平均--------------------
tmp_june_mean=runave_n_Wrap(tmp_june,9,0,0)
tmp_july_mean=runave_n_Wrap(tmp_july,9,0,0)
tmp_aug_mean=runave_n_Wrap(tmp_aug,9,0,0)
tmp6=tmp_june(:,:,:)-tmp_june_mean(:,:,:)
tmp7=tmp_july(:,:,:)-tmp_july_mean(:,:,:)
tmp8=tmp_aug(:,:,:)-tmp_aug_mean(:,:,:)
tmp_60_summer=(tmp6+tmp7+tmp8)/3
; tmp_60_summer=(tmp_june_mean+tmp_june_mean+tmp_aug_mean)/3

copy_VarMeta(tmp_june_mean, tmp_60_summer )

tmp_jp   =dim_rmvmean_n_Wrap(tmp_60_summer, 0)
copy_VarMeta(tmp_june_mean, tmp_jp )
      ; --------------求区域平均值的时间序列-------------
   t_678=wgt_areaave(tmp_60_summer(:,{45:70},{30:75}), clat(90:139), 1.0, 0)
   tmp_678=dim_rmvmean_n_Wrap(t_678, 0)
   data       = new((/4,61/), "float")
   data(0,:) = tmp_678(:)
   printVarSummary(tmp_60_summer)
  ; --------------------------计算ENSO指数------------------------------
; ===================================================================

  b = addfile("/disk1/lvcy/data/sst.mnmean.ersstv5.nc","r")
; sst= short2flt(b->sst(1248:1979,:,:));1958-2018年
  sstt= short2flt(b->sst(:,:,:));1958-2018年
  time2 =sstt&time
  sst=dtrend_msg_n(time2,sstt, False, True, 0)
  copy_VarMeta(sstt, sst)
  printVarSummary(sst)
  istart=195801
  iend=201812
  dym=cd_calendar(sst&time,-1)
  ist=ind(istart.eq.dym)
  ied=ind(iend.eq.dym)
  n_day = dimsizes(sst&time)
  n_lat = dimsizes(sst&lat)
  n_lon = dimsizes(sst&lon)
  sst_djf=month_to_season(sst(ist:ied,:,:),"JJA")
  printVarSummary(sst_djf)
  nino3_4 =dim_avg_n_Wrap(sst_djf(:,{-5:5},{190:240}),(/1,2/))
  printVarSummary(nino3_4)
  nino3_4_pic=dim_standardize_n_Wrap(nino3_4,0,0)
  printVarSummary(nino3_4_pic)
  ind_lnn = ind( nino3_4_pic .lt. -0.5);+1958
  ind_enn = ind( nino3_4_pic .gt. 0.5);+1958
  ; print(ind_enn)
  ; print(ind_lnn)
   latS   =  -90
   latN   =  90.
   lonL   =  0
   lonR   =  360.
  nino3 =dim_avg_n_Wrap(sst_djf(:,{-5:5},{210:270}),(/1,2/))
  nino4 =dim_avg_n_Wrap(sst_djf(:,{-5:5},{160:210}),(/1,2/))
  nino3_pic=dim_standardize_n_Wrap(nino3,0,0)
  nino4_pic=dim_standardize_n_Wrap(nino4,0,0)
  copy_VarMeta(nino3_4_pic, tmp_678)
; --------------------------------回归------------------------------
    rc =new((/360,720/),float,"No_FillValue")
    rc(:,:) =regCoef(nino3_4_pic, tmp_60_summer({lat|:},{lon|:},{time|:}));回归!!!!!!!!!!

    ; rc_u =new((/2,37,73/),float,"No_FillValue")
    ; rc_u(:,:,:) =regCoef(eof_ts(:,0:58), u_jp1({lat|:},{lon|:},{time|:}));回归!!!!!!!!!!
    ; rc_v =new((/2,37,73/),float,"No_FillValue")
    ; rc_v(:,:,:) =regCoef(eof_ts(:,0:58), v_jp1({lat|:},{lon|:},{time|:}));回归!!!!!!!!!!

    rc!0 = "lat"
    rc!1 = "lon"
    ; env = (/1,2/)
    ; lat = fspan(2.5,0,90)
    ; lon = fspan(2.5,0,180)
    lat = tmp_60_summer&lat
    lon = tmp_60_summer&lon
    rc&lat = lat
    rc&lon = lon
    lat@units = "degrees_north"
    lon@units = "degrees_east"
    rc&lat@units = "degrees_north"
    rc&lon@units = "degrees_east"  
; ;-----------------------------去除ENSO影响--------------------------
; ----------------------------减去线性回归的温度场-------------------
noenso_sat =new((/61,360,720/),float,"No_FillValue")
do n=0, 60
noenso_sat(n,:,:)= tmp_60_summer(n,:,:) - rc(:,:)
end do
copy_VarMeta(tmp_60_summer, noenso_sat)
printVarSummary(noenso_sat)
printMinMax(noenso_sat, -1)
; exit
; ------------------------将去除ENSO影响的SAT写成nc文件----------------
fout = addfile("no_ENSO_SAT2.nc", "c")
filedimdef(fout, "time", -1, True)
fout->tmp = noenso_sat
end
