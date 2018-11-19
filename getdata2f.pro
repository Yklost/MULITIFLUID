function div, bx,by,x,y,n_x,n_y
divb=dblarr(n_x,n_y)
for i=0,n_y-1 do divb(*,i)=deriv(bx(*,i))
for j=0,n_x-1 do divb(j,*)=divb(j,*)+deriv(by(j,*))
return,divb
end


gamma=5./3.

path = './'

file = file_search(path+'data_*', COUNT=nFiles)


n_x = 256
n_y = 256
n_var = 12

x=dblarr(n_x,n_y)
y=dblarr(n_x,n_y)
rho=dblarr(n_x,n_y)
rho_i_t=dblarr(nFiles)
rho_n_t=dblarr(nFiles)
rho_i_bound_left=dblarr(nFiles)
rho_i_bound_right=dblarr(nFiles)

for i = 0, nFiles-1 do begin
  close, 1
  print,file[i]

      
  openr,1, file[i] 
  a = assoc(1, dblarr(n_x, n_y))
  x   = a(0)
  y   = a(1)
  rho_i = a(2)
  m_x_i = a(3)
  m_y_i = a(4)
  m_z_i = a(5)
  ene_i = a(6)
  e_x = a(7)
  e_y = a(8)
  e_z = a(9)
  rho_n = a(10)
  m_x_n = a(11)
  m_y_n = a(12)
  m_z_n = a(13)
  ene_n = a(14)
  ene_e = a(15)
  v_x_i = a(16)
  v_y_i = a(17)
  v_z_i = a(18)
  pre_i = a(19)
  v_x_n = a(20)
  v_y_n = a(21)

  v_z_n = a(22)
  pre_n = a(23)
  pre_e = a(24)
  k_ion = a(25)
  k_re = a(26)
  a_ni = a(27)
  a_ne = a(28)
  qqq = a(29)
  q_n = a(30)
  q_x = a(31)
  q_y = a(32)
  q_z = a(33) 
  tem_i = a(34)
  tem_n = a(35)
  tem_e = a(36)
  vidte = a(37)

  rho_i_t[i]=total(rho_i,/double)
  rho_n_t[i]=total(rho_n,/double)

  rho_i_bound_left[i]=rho_i(22,0)
  ;rho_i_bound_right[i]=rho_i(22,255)


  close,1





e_k_i = rho_i * (v_x_i*v_x_i + v_y_i*v_y_i + v_z_i*v_z_i) / 2.d0
e_k_n = rho_n * (v_x_n*v_x_n + v_y_n*v_y_n + v_z_n*v_z_n) / 2.d0



e_t_i = ene_i - e_k_i
e_t_n = ene_n - e_k_n    


cs=2.0

!p.multi=[0,2,2,0,1]

;PLOT,y(22,*),rho_i(22,*),charsize=1
;w = WINDOW( $
;  WINDOW_TITLE='Mass density and temperature profile')
;  p1=PLOT(y(22,*),rho_i(22,*),LAYOUT=[2,2,1])
;  p2=PLOT(y(22,*),rho_n(22,*),LAYOUT=[2,2,2])prot_velocities_128_256_512_1024_domain_size.pngprot_velocities_128_256_512_1024_domain_size.png

;  p3=PLOT(y(22,*),tem_i(22,*),LAYOUT=[2,2,1])
;  p4=PLOT(y(22,*),tem_n(22,*),LAYOUT=[2,2,1])

;p1=PLOT(y(22,*),ALOG10(rho_i(22,*)),/OVERPLOT,YTITLE='Proton mass density with fixed boun cond, log',XTITLE='Height,cm')
;p1=PLOT(y(22,*),rho_n(22,*),/OVERPLOT,XTITLE='b')
;p2=PLOT(y(22,*),tem_i(22,*),/OVERPLOT,XTITLE='c')
;p3=PLOT(y(22,*),tem_n(22,*),/OVERPLOT,XTITLE='d')

;p3=PLOT(y(22,*),rho_n(22,*),/OVERPLOT)

;PLOT,rho_i_t,/ys,CHARSIZe=2
;PLOT,rho_n_t,/ys,CHARSIZe=2

;PLOT,rho_i_bound_left,/ys,charsize=2
;p3=PLOT(rho_i_bound_right,XTITLE='Timesteps*100',YTITLE='Proton mass density top boundary')

;p3=PLOT(y(22,*),tem_i(22,*),/OVERPLOT,XTITLE='Height,cm',YTITLE='Proton temperature, eV')




;p1=PLOT(y(22,*),rho_i(22,*),/OVERPLOT,XTITLE='Height,cm',YTITLE='Proton mass density')
;p2=PLOT(y(22,*),v_y_i(22,*),/OVERPLOT,XTITLE='Height,cm',YTITLE='Proton vertical velocity',COLOR='Green')
;p3=PLOT(y(22,*),tem_i(22,*),/OVERPLOT, XTITLE='Height,cm',YTITLE='Proton temperature')
;p4=PLOT(y(22,*),ene_i(22,*),/OVERPLOT, XTITLE='Height,cm',YTITLE='Proton energy')
;p9=PLOT(y(22,*),pre_i(22,*),/OVERPLOT, XTITLE='Height,cm',YTITLE='Proton pressure')

;p5=PLOT(y(22,*),rho_n(22,*),/OVERPLOT, XTITLE='Height,cm',YTITLE='Hydrogen mass density')
;p6=PLOT(y(22,*),v_y_n(22,*),/OVERPLOT, XTITLE='Height,cm',YTITLE='Hydrogen vertical velocity')
;p7=PLOT(y(22,*),tem_n(22,*),/OVERPLOT, XTITLE='Height,cm',YTITLE='Hydrogen temperature')
;p8=PLOT(y(22,*),ene_n(22,*),/OVERPLOT, XTITLE='Height,cm',YTITLE='Hydrogen energy')
;p1=PLOT(y(22,*),pre_n(22,*),/OVERPLOT, XTITLE='Height,cm',YTITLE='Hydrogen pressure')

;p6=SCATTERPLOT(y(22,*),v_y_n(22,*),/OVERPLOT, SYMBOL='dot',XTITLE='Height,cm',YTITLE='Hydrogen vertical velocity')



;tvframe,rho_i,/ba,/sa,btitle='rho_i',charsize=cs,/as
;tvframe,v_x_i,/ba,/sa,btitle='v_x_i',charsize=cs,/as
;tvframe,v_y_i,/ba,/sa,btitle='v_y_i',charsize=cs,/as
;tvframe,ene_i,/ba,/sa,btitle='ene_i',charsize=cs,/as
;tvframe,tem_i,/ba,/sa,btitle='tem_i',charsize=cs,/as


;tvframe,rho_n,/ba,/sa,btitle='rho_n',charsize=cs,/as
;tvframe,v_x_n,/ba,/sa,btitle='v_x_n',charsize=cs,/as
;tvframe,v_y_n,/ba,/sa,btitle='v_y_n',charsize=cs,/as
;tvframe,rho_i/rho_n,/ba,/sa,btitle='rho_i/rho_n',charsize=cs,/as
;tvframe,tem_n,/ba,/sa,btitle='tem_n',charsize=cs,/as

;print,size(rho_i)

;tvframe,ene_e,/ba,/sa,btitle='ene_e',charsize=cs,/as
;tvframe,vidte,/ba,/sa,btitle='vidte',charsize=cs,/as

as=''
read,as



endfor

close, 1

;openw,lun,'data.dat',/get_lun
;printf,lun,256
;printf,lun,format='(3e20.10)',80000000,150000000,273437.5
;for j=0,255 do begin
;printf,lun,format='(15e20.10)',y(22,j),rho_i(22,j)/2000,rho_i(22,j),rho_n(22,j),v_y_i(22,j),v_y_n(22,j),pre_e(22,j),pre_i(22,j),pre_n(22,j)
;endfor
;Free_lun,lun
;close,lun

end
