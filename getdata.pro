function div, bx,by,x,y,n_x,n_y
divb=dblarr(n_x,n_y)
for i=0,n_y-1 do divb(*,i)=deriv(bx(*,i))
for j=0,n_x-1 do divb(j,*)=divb(j,*)+deriv(by(j,*))
return,divb
end


gamma=5./3.

path = './'

file = file_search(path+'data_*', COUNT=nFiles)


n_x = 512
n_y = 512
n_var = 12

x=dblarr(n_x,n_y)
y=dblarr(n_x,n_y)
rho=dblarr(n_x,n_y)


for i = 0, nFiles-1 do begin
  close, 1
  print,file[i]

      
  openr,1, file[i]  
  a = assoc(1, dblarr(n_x, n_y))
  x   = a(0)
  y   = a(1)
  rho = a(2)
  m_x = a(3)
  m_y = a(4)
  m_z = a(5)
  ene = a(6)
  b_x = a(7)
  b_y = a(8)
  b_z = a(9)
  close,1

v_x=m_x/rho
v_y=m_y/rho
v_z=m_z/rho

absv=sqrt(v_x^2.0+v_y^2.0+v_z^2.0)

p_gas = (ene - (b_x*b_x + b_y*b_y + b_z*b_z)*0.5 - (m_x*m_x + m_y*m_y + m_z*m_z)/(2.0*rho))*(gamma-1.0)
p_mag = (b_x*b_x + b_y*b_y + b_z*b_z)*0.5
c_sou = sqrt(gamma*p_gas/rho)
mach = absv/c_sou

cs=2.0

;!p.multi=[0,4,2]
;!p.multi=0
!p.multi=[0,2,1]

divb = div(b_x,b_y,x,y,n_x,n_y)
beta = p_mag/p_gas

tvframe,rho,/ba,/sa,btitle='rho',charsize=cs,/as
;tvframe,m_x,/ba,/sa,btitle='m_x',charsize=cs
;tvframe,m_y,/ba,/sa,btitle='m_y',charsize=cs
tvframe,divb,/ba,/sa,btitle='divb',charsize=cs
;tvframe,ene,/ba,/sa,btitle='ene',charsize=cs
;tvframe,b_x,/ba,/sa,btitle='b_x',charsize=cs
;tvframe,b_y,/ba,/sa,btitle='b_y',charsize=cs
;tvframe,beta,/ba,/sa,btitle='beta',charsize=cs

as=''
read,as



endfor


close, 1

end
