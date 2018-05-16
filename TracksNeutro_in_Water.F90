!-------realized by Bilal Elfallah and ibrahim amrazguiou
program track_neutron_in_water
implicit none
real::diff
integer::i,n,mh2o,fuit
real::Eint,sigmath,sigmadh,sigmaah,phi,teta,r,L,x0,y0,z0,gsigmath,cm,moy,somL,u,r1
real::sigmato,sigmado,sigmaao,gsigmat,gsigmato,sigmaa,sigmat
n=200 
  r1=20.     ! Rayon de la sphere 
moy=0.  
x0=0;y0=0;z0=0
u=0.
fuit=0
print*,' E_int       nbr_scattering    '
do i=1,n
Eint=2.E6
diff=0.
somL=0.

do
call g(Eint,sigmath,sigmadh,sigmaah)
call f(Eint,sigmato,sigmado,sigmaao)
r=rand()
gsigmath=2*((0.6023)*sigmath/18)!section effecace total macroscopique de H
gsigmato=(((0.6023)*sigmato)/18) !section effecace total macroscopique de O-16
gsigmat=gsigmato+gsigmath
teta=acos(1-2*r)
phi=2*4*atan(1.)*r
L=-log(rand())/gsigmat
somL=somL+L
x0=x0+L*sin(teta)*cos(phi)            
y0=y0+L*sin(teta)*sin(phi)
z0=z0+L*cos(teta)
if(r<(gsigmath/gsigmat)) then !si le noyau heurté celui de Hy
   mh2o=1
   sigmat=sigmath
   sigmaa=sigmaah
   else                     !si l'oxygene qui est heurté
   mh2o=16
   sigmat=sigmato
   sigmaa=sigmaao
        
end if
   

       if (r<=(sigmaa/sigmat)) then
       exit
       else
       cm=2*rand()-1         !cm simbolyse cos(teta) dans le SCM
       Eint=(((1+mh2o*mh2o+2*mh2o*cm))/((1+mh2o)*(1+mh2o)))*Eint
       diff=diff+1
          if(Eint<1.) then
          exit
          end if             
endif        
enddo
print*,Eint,diff
u=u+somL
moy=moy+diff
if(sqrt(x0*x0+y0*y0+z0*z0)>r1) then !Test sur une géometrie spherique
fuit=fuit+1
endif
enddo
print*,'la valeur moyenne du nombre de diffusions est:',moy/n
write(*,*)'le nombre de fuite est :',fuit
end  program track_neutron_in_water
!------------------Lecture et interpolation des sections effecases pour H ---------------------------
subroutine g(Eint,sigmath,sigmadh,sigmaah)
implicit none
real,dimension(307)::sigmatotalh,Tn,sigmadiffh,sigmaabsoh
real::sigmath,Eint,sigmadh,sigmaah
integer::i
     open(unit=2,file='1001.txt')  ! read cross section of hydrogene
  do i=1,307
     read(2,*)Tn(i),sigmatotalh(i),sigmadiffh(i),sigmaabsoh(i)
  enddo
     close(2)
    do i=1,307
      if(Eint>=Tn(i) .and. Eint<=Tn(i+1)) then
       sigmath=((sigmatotalh(i)-sigmatotalh(i+1))/(Tn(i)-Tn(i+1)))*(Eint-Tn(i))+sigmatotalh(i)
       sigmadh=((sigmadiffh(i)-sigmadiffh(i+1))/(Tn(i)-Tn(i+1)))*(Eint-Tn(i))+sigmadiffh(i)
       sigmaah=((sigmaabsoh(i)-sigmaabsoh(i+1))/(Tn(i)-Tn(i+1)))*(Eint-Tn(i))+sigmaabsoh(i)
       endif
   enddo
end subroutine 
!---------------------------Lecture et interpolation des sections effecases pour H ---------------------------
subroutine f(Eint,sigmato,sigmado,sigmaao)
implicit none
real,dimension(2509)::sigmatotalo,Tn,sigmadiffo,sigmaabsoo
real::sigmato,Eint,sigmado,sigmaao
integer::i
   open(unit=2,file='8016.txt')  ! read cross section of Oxygene
  do i=1,2503
     read(2,*)Tn(i),sigmatotalo(i),sigmadiffo(i),sigmaabsoo(i)
 enddo
    close(2)
          do i=1,2503
    if(Eint>=Tn(i) .and. Eint<=Tn(i+1)) then
sigmato=((sigmatotalo(i)-sigmatotalo(i+1))/(Tn(i)-Tn(i+1)))*(Eint-Tn(i))+sigmatotalo(i)
sigmado=((sigmadiffo(i)-sigmadiffo(i+1))/(Tn(i)-Tn(i+1)))*(Eint-Tn(i))+sigmadiffo(i)
sigmaao=((sigmaabsoo(i)-sigmaabsoo(i+1))/(Tn(i)-Tn(i+1)))*(Eint-Tn(i))+sigmaabsoo(i)
    endif
         enddo
end subroutine 
 
