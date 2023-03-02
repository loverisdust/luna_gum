begin

diri = "./cloud/"
fili = systemfunc("cd "+diri+" ;ls *.txt")
nfili= dimsizes(fili)   ;;;nfili = 684
cloudy = new((/nfili,2072,6/), "float")
do i = 0,nfili-1,1
cloudy(i,:,:) = asciiread(diri+fili(i), (/2072,6/) , "float")
end do

;;;挑选东部站点
lon = cloudy(0,:,1)
lat = cloudy(0,:,2)
sta = ind(lon .ge. 105)
lon := lon(sta)
lat := lat(sta)
; printMinMax(lon, False)    ;;;min = 105,max = 132.93
; printMinMax(lat, False)    ;;;min = 16.83,max = 53.48
cloud_Jan = cloudy(0:672:12,:,5)
cloud_Jan := cloud_Jan(:,sta)
cloud_Feb = cloudy(1:673:12,:,5)
cloud_Feb := cloud_Feb(:,sta)
cloud_Dec = cloudy(11:(nfili-1):12,:,5)
cloud_Dec := cloud_Dec(:,sta)
nsta= dimsizes(sta)  ;;;nsta = 1624

;;;异常数据处理
cloud_Jan = where(cloud_Jan .gt. 31, 32766, cloud_Jan)
cloud_Feb = where(cloud_Feb .gt. 28, 32766, cloud_Feb)
cloud_Dec = where(cloud_Dec .gt. 31, 32766, cloud_Dec)


;;;冬季数据
cloud_win = cloud_Jan+cloud_Feb+cloud_Dec/3
cloud_win@_FillValue = 32766
cloud_win@long_name="cloudy_days"
cloud_win@units = "days"
cloud_win!0 = "time"
time = ispan(1961,2017,1)
time@units = "year"
cloud_win&time = time
cloud_win!1 = "var"
var= ispan(1,nsta,1)
var@units = "days"
cloud_win&var = var
cloud_win = where(cloud_win .gt. 30, 32766, cloud_win)
cloud_win@_FillValue = 32766
;printVarSummary(cloud_win)

;;;气候态
cloud_cli = new(nsta, "float")
do i = 0,nsta-1,1
cloud_cli(i) = avg( cloud_win(:,i) )
end do

; cloud_win = dtrend(cloud_win, False)
;;eof、时间序列及north检验
neval = 4
optEOF = True
optEOF@jopt = 0
eof = eofunc_Wrap(cloud_win(var|:,time|:), neval, optEOF)
;printVarSummary(eof)

eof_ts=eofunc_ts_Wrap(cloud_win(var|:,time|:), eof, False)
eof_ts = dim_standardize_Wrap(eof_ts,1)
;printVarSummary(eof_ts)
eof_ts@long_name = "Amplitude"

eval = eof@eval
e = eval*sqrt(2.0/nsta)
nor = new(neval, "float")
north = new(neval, "float")
do i = 0,neval-2
nor(i) = eval(i)-eval(i+1)
if (nor(i) .ge. e(i)) then 
north(i) = 1
else
north(i) = 0
end if
end do

option = True
option@row = False
option@tspace = 0
option@fout = "./tu/hhh_winter.txt"
option@title = "eig"
hhh = new((/3,neval/),"float")
hhh(0,:) = eval
hhh(1,:) = eof@pcvar
hhh(2,:) = north
write_matrix(hhh, "4f20.2", option)

;;;插值
olon = new(29,"float");
olat = new(39,"float");
do i=0,28
olon(i) =105+i
end do
do l=0,38
olat(l) = 16+l
end do
olon!0          = "lon" 
olon@long_name  = "lon"
olon@units      = "degrees-east"
olon&lon        = olon
olat!0          = "lat"
olat@long_name  = "lat"
olat@units      = "degrees_north"
olat&lat        = olat
rscan = (/5,2.5,1/)

EOF = new((/neval,39,29/), "float")
do i = 0,neval-1
EOF(i,:,:)  = obj_anal_ic_Wrap(lon,lat,eof(i,:),olon,olat,rscan, False)
end do
;printVarSummary(EOF)

;;;模态极值
option = True
option@row = False
option@tspace = 0
option@fout = "./tu/jizhi_winter.txt"
option@title = "模态极值"
k = new((/4,2/),"float")
do i = 0,3,1
k(i,0) = max(EOF(i,:,:))
k(i,1) = min(EOF(i,:,:))
end do
write_matrix(k, "4f10.5", option)

;;;气候态出图
CLI = obj_anal_ic_Wrap(lon,lat,cloud_cli,olon,olat,rscan, False)
wks = gsn_open_wks("png", "./tu/climate_winter")
gsn_define_colormap(wks, "BlRe")
res1= True
res1@gsnLeftString  = "cloudy"
res1@gsnRightString = "days"
res1@pmTickMarkDisplayMode="Always"
res1@gsnAddCyclic  = False      ;由于我们的数据不是循环地球一周的，因此必须把这个置否
res1@mpDataSetName         = "Earth..4"   ; This new database contains
;                                        ; divisions for other countries.
res1@mpDataBaseVersion     = "MediumRes"  ; Medium resolution database
res1@mpOutlineOn           = True         ; Turn on map outlines
res1@mpOutlineSpecifiers   = (/"China:states","Taiwan"/)
res1@cnLevelSelectionMode = "ExplicitLevels"
res1@cnLevels    = ispan(3, 27, 3)
res1@cnFillColors = (/85,75,65,55,48,45,40,35,25,20/) 
res1@mpMinLatF=min(olat)
res1@mpMaxLatF=max(olat)
res1@mpMinLonF=min(olon)
res1@mpMaxLonF=max(olon)
res1@lbLabelBarOn = True       ;LabelBar显示
res1@lbOrientation="Vertical"
res1@pmLabelBarOrthogonalPosF = .20
res1@mpAreaMaskingOn = True   ;使能填充覆盖
res1@mpMaskAreaSpecifiers = (/"China:states","Taiwan"/)   ;China:states
res1@mpOceanFillColor = 0     ;用白色填充海洋  0是colormap的索引值
res1@mpInlandWaterFillColor = 0  ;用白色填充内陆湖水
res1@mpLandFillColor=0
res1@mpGeophysicalLineColor=0
res1@cnFillOn      = True
res1@cnLinesOn     = False

res1@cnFillDrawOrder = "PreDraw"         ; draw contours first 
plot2 = gsn_csm_contour_map(wks,CLI,res1)

;;;模态出图
n = 2
nscan = 2000             ; usually *much* fewer
eps   = 0.001            ; variable depended
gtype = False            ; "gendat" does not generate cyclic fields
guess = 0                ; use zonal means
relc  = 0.6              ; standard relaxation coef
opt   = 0
poisson_grid_fill(EOF(n,:,:), gtype, guess, nscan, eps, relc, opt)  ;;;高斯插值
printMinMax(EOF(n,:,:),False)

wks1 = gsn_open_wks("png", "./tu/motai"+(n+1)+"_winter")
gsn_define_colormap(wks1, "hotcolr_19lev")
res= True
res@pmTickMarkDisplayMode="Always"
res@gsnAddCyclic  = False      ;由于我们的数据不是循环地球一周的，因此必须把这个置否
res@mpDataSetName         = "Earth..4"   ; This new database contains
;                                        ; divisions for other countries.
res@mpDataBaseVersion     = "MediumRes"  ; Medium resolution database
res@mpOutlineOn           = True         ; Turn on map outlines
res@mpOutlineSpecifiers   = (/"China:states","Taiwan"/)
res@mpMinLatF=min(olat)
res@mpMaxLatF=max(olat)
res@mpMinLonF=min(olon)
res@mpMaxLonF=max(olon)

res@mpAreaMaskingOn = True   ;使能填充覆盖
res@mpMaskAreaSpecifiers = (/"China:states","Taiwan"/)   ;China:states
res@mpOceanFillColor = 0     ;用白色填充海洋  0是colormap的索引值
res@mpInlandWaterFillColor = 0  ;用白色填充内陆湖水
res@mpLandFillColor=0
res@mpGeophysicalLineColor=0

res@cnFillOn      = True
res@cnLinesOn     = False
res@cnFillDrawOrder = "PreDraw"         ; draw contours first 
res@cnLevelSelectionMode = "ExplicitLevels"       ; set explicit contour levels
; res@cnLevels    = (/-0.003,0,0.02,0.04,0.006/) 
; res@cnFillColors = (/6,8,13,15,17,20/) 
; res@cnLevels    = fspan(-0.04, 0.06, 11)
; res@cnFillColors = (/3,4,7,8,9,13,14,17,18,19,20,21/) 
res@cnLevels    = fspan(-0.04, 0.06, 11)
res@cnFillColors = (/3,4,6,7,8,12,13,14,17,18,19,20/) 
; res@cnLevels    = (/-0.08,-0.06,-0.04,-0.02,-0.01,0,0.01,0.02,0.03,0.04/) 
; res@cnFillColors = (/3,4,5,6,7,8,13,17,18,19,20/) 
res@lbLabelBarOn = True       ;LabelBar显示
res@lbOrientation="Vertical"
res@pmLabelBarOrthogonalPosF = .20

; do n=0,neval-1
; res@gsnLeftString  = "EOF "+(n+1)
; res@gsnRightString = sprintf("%5.1f", eof@pcvar(n)) +"%"
; plot(n) = gsn_csm_contour_map(wks+i,EOF(n,:,:),res)
; end do
res@gsnLeftString  = "EOF"+(n+1)
res@gsnRightString = sprintf("%5.1f", eof@pcvar(n)) +"%"
plot = gsn_csm_contour_map(wks1,EOF(n,:,:),res)

m = 2
printMinMax(eof_ts(m,:),False)
wks2 = gsn_open_wks("png", "./tu/motai_ts"+(m+1)+"_winter")
rts = True
rts@vpHeightF = 0.40        ; Changes the aspect ratio
rts@vpWidthF  = 0.85
rts@vpXF      = 0.10        ; change start locations
rts@vpYF      = 0.75        ; the plot  
rts@trYMaxF = 3
rts@trYMinF = -3
rts@gsnYRefLine           = 0.              ; reference line   
rts@gsnXYBarChart         = True            ; create bar chart 
rts@gsnAboveYRefLineColor = "red"           ; above ref line fill red
rts@gsnBelowYRefLineColor = "blue"          ; below ref line fill blue
rts@gsnLeftString  = "PC"+(m+1)
rts@gsnRightString = sprintf("%5.1f", eof@pcvar(m)) +"%"
plot1 = gsn_csm_xy(wks2,time,eof_ts(m,:),rts)

end