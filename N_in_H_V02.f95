! realized by bilal el fallah
program N_chocs_in_Hy
implicit none 
real ::siigmma_scattering,siigmma_absorption,siigmma_total
real ::nbr_scattering,nbr_scattering_moy,E_int,grand_sigmat
real :: COS_THETA_CM ,R ,x0,y0,z0,L,theta,phi,mean_path
integer ::i
x0=0 ; y0=0 ; z0=0 
print*,'   E_out         Nbr_scatering      Path' 
nbr_scattering_moy=0.
 mean_path=0.
do i=1,500
         E_int=2.E6 
         nbr_scattering=0.
         do
 CALL g(E_int,siigmma_absorption,siigmma_scattering,siigmma_total)
               
               grand_sigmat=((6.023E23)*siigmma_total)/9 
               CALL RANDOM_NUMBER(R)
               
                theta=acos(2*R-1)
                 phi=2*4*atan(1.)*R
                 L=-log(R)/grand_sigmat
                 x0=x0+L*sin(theta)*cos(phi)            
                 y0=y0+L*sin(theta)*sin(phi)
                 z0=z0+L*cos(theta)
               
       if(R <= siigmma_absorption/siigmma_total) then 
       exit
       else        
               COS_THETA_CM=-1+2*R
               E_int=E_int*(1+COS_THETA_CM)/2
               nbr_scattering=nbr_scattering+1
               if (E_int <= 1.) then 
               exit
               end if
               endif
         enddo
        
    print*,E_int,nbr_scattering,L   
   nbr_scattering_moy = nbr_scattering_moy + nbr_scattering
   mean_path=mean_path+L
 print*,
end do
          print*,'the mean number of scatering is :',nbr_scattering_moy/i
           print*,'the mean path of neutron is :',mean_path/i                  
end program N_chocs_in_Hy
subroutine  g(E_int,siigmma_absorption,siigmma_scattering,siigmma_total)
implicit none 
real,dimension(307):: sigmma_absorption,sigmma_scattering,sigmma_total,Tn
real  :: siigmma_absorption,siigmma_scattering,siigmma_total,E_int
integer :: i
open(unit=4,file='1001.csv')
do i=1,307
read(4,*)Tn(i),sigmma_total(i),sigmma_scattering(i),sigmma_absorption(i)
end do 
close(4)
     do i=1,307
            if (E_int>=Tn(i) .and. E_int <= Tn(i+1))  then 
siigmma_total=((sigmma_total(i)-sigmma_total(i+1))/(Tn(i)-Tn(i+1)))*(E_int-Tn(i))+sigmma_total(i)
siigmma_scattering=((sigmma_scattering(i)-sigmma_scattering(i+1))/(Tn(i)-Tn(i+1)))*(E_int-Tn(i))+sigmma_scattering(i)
siigmma_absorption=((sigmma_absorption(i)-sigmma_absorption(i+1))/(Tn(i)-Tn(i+1)))*(E_int-Tn(i))+sigmma_absorption(i)              
           endif
           end do
end subroutine           


          
