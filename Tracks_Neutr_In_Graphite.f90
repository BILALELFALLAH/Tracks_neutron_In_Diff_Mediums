!   realized by BILAL and IBRAHIM
program tracks_neutron_graphite
implicit none
real::diff
integer::i,n
real::Eint,sigmat,sigmad,sigmaa,phi,teta,r,L,x0,y0,z0,gsigmat,cm,moy,somL,u
n=150 
moy=0.
u=0.
WRITE(*,*)'  E_int     nbr_scattering      path_neutron'
do i=1,n
Eint=2.E6
diff=0.
somL=0.
x0=0;y0=0;z0=0
do
call g(Eint,sigmat,sigmad,sigmaa)
gsigmat=(((0.6023*2.2)/12))*sigmat ! gsigmat symbolyse sigmma macrosco total
L=-log(rand())/gsigmat
somL=somL+L
r=rand()
teta=acos(2*r-1)
phi=2.*4*atan(1.)*rand()
x0=x0+L*sin(teta)*cos(phi)            
y0=y0+L*sin(teta)*sin(phi)
z0=z0+L*cos(teta)
     if (r<=(sigmaa/sigmat)) then
     exit
         else
    
         cm=2*rand()-1      !cm=cos(tetacm)
         Eint=(Eint*(1+12*12+26*cm))/(13*13)
         diff=diff+1
        
              if(Eint<=1.) then
                 exit
              end if
             
      endif

enddo
u=u+somL
print*,Eint,diff,L
moy=moy+diff
enddo
print*,'le nombre de diffusion moy est :',moy/n
end program tracks_neutron_graphite
!---------------------------interpolation and read cross sections of C-12---------------------------
subroutine g(Eint,sigmat,sigmad,sigmaa)
implicit none
real,dimension(1055)::sigmatotal,Tn,sigmadiff,sigmaabso
real::sigmat,Eint,sigmad,sigmaa
integer::i
open(unit=5,file='carbone.txt')
do i=1,1021
     read(5,*)Tn(i),sigmatotal(i),sigmadiff(i),sigmaabso(i)
enddo
close(5)
do i=1,1021
    if(Eint>=Tn(i) .and. Eint<=Tn(i+1)) then
       sigmat=((sigmatotal(i)-sigmatotal(i+1))/(Tn(i)-Tn(i+1)))*(Eint-Tn(i))+sigmatotal(i)
       sigmad=((sigmadiff(i)-sigmadiff(i+1))/(Tn(i)-Tn(i+1)))*(Eint-Tn(i))+sigmadiff(i)
        sigmaa=((sigmaabso(i)-sigmaabso(i+1))/(Tn(i)-Tn(i+1)))*(Eint-Tn(i))+sigmaabso(i)
    endif
enddo
end subroutine
