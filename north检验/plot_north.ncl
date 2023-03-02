load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"

begin
;
; Create workstation.
;
  wks = gsn_open_wks("png","add_primitive")

;
; Create a data object.
;
  npts = 500
  x    = fspan(0,npts-1,npts)
  y    = 500.+ 0.9 * x * sin(0.031415926535898*x)

  xyres = True
  xyres@gsnMaximize = True
  xy = gsn_y(wks,y,xyres)

;
; Create data for primitives.
;
  plx = x
  ply = 500. + 0.5 * x * sin(0.031415926535898*x)
  pgx = (/ 100., 200., 200., 100., 100. /)
  pgy = (/ 200., 200., 300., 300., 200. /)
  pmx = (/ 200., 150., 200., 250. /)
  pmy = (/ 900., 800., 700., 800. /)

;
; Set up three separate resource lists, although we could have 
; used the same one here.
;
  pmres                  = True
  plres                  = True
  pgres                  = True
  plres@gsLineColor      = "orange"
  plres@gsLineThicknessF = 2.0
  pgres@gsFillColor      = "Navy"
  pmres@gsMarkerIndex    = 12
  pmres@gsMarkerSizeF    = 0.02
  pmres@gsMarkerColor    = "yellow"

;
; Make sure each variable name is unique.
;
  dum1 = gsn_add_polyline  (wks, xy, plx, ply, plres)
  dum2 = gsn_add_polygon   (wks, xy, pgx, pgy, pgres)
  dum3 = gsn_add_polymarker(wks, xy, pmx, pmy, pmres)

  draw(xy)
  frame(wks)

;
; Resize the plot and draw in the middle of the frame, and
; see how the primitives automatically get adjusted as well.
; 
  setvalues xy
    "vpXF"        : 0.3
    "vpYF"        : 0.7
    "vpWidthF"    : 0.4
    "vpHeightF"   : 0.4
  end setvalues

  draw(xy)
  frame(wks)
end