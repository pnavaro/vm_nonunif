module diagno

use zone

implicit none

integer, private :: i, j, k

contains

subroutine plot_part( ele, time, iplot )

type (particle) :: ele
integer :: iplot
real(kind=prec) :: time, umod

if (nomcas == 'viry__') then
   open(40, file="viry.dat", position="append")
   if( iplot == 1 ) rewind(40)
   
   umod = sqrt(ele%vit(1,1)**2+ele%vit(1,2)**2)
   
   write(40,"(G15.3,1X,5F12.7)") sngl(time)	&
        ,sngl(ele%pos(1,1)), sngl(ele%pos(1,2))    &
        ,sngl(ele%vit(1,1)), sngl(ele%vit(1,2))    &
        ,umod
   close(40)
end if

end subroutine plot_part

subroutine diag_champ_part( ele, time, iplot )

type(particle) :: ele
real(kind=prec) :: time
integer :: iplot

if ( nbpart > 0 ) then

   open(17, file="chpart", position="append"  )
   if( iplot == 1 ) rewind(17)
   
   write(17,*) sngl(time), sngl(ele%epx(1)), &
        sngl(ele%epy(1)), sngl(ele%bpz(1))

   close(17)

endif

end subroutine diag_champ_part

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine plot_champ( tm, iplot, time )

type(mesh_fields) :: tm
integer :: iplot
real(kind=prec) :: time

integer :: kk0, kk1, kk2, kk3, kk4
character(len=4) :: fin
character(len=1) :: aa,bb,cc,dd

kk0 = iplot
kk1 = kk0/1000
aa  = char(kk1 + 48)
kk2 = (kk0 - kk1*1000)/100
bb  = char(kk2 + 48)
kk3 = (kk0 - (kk1*1000) - (kk2*100))/10
cc  = char(kk3 + 48)
kk4 = (kk0 - (kk1*1000) - (kk2*100) - (kk3*10))/1
dd  = char(kk4 + 48)
fin = aa//bb//cc//dd
   
open(  11, file = 'ex.gnu', position="append" )
if ( iplot .eq. 1 ) rewind(11)
write(11,"(A18,G10.3,A1)")"set title 'Time = ",time,"'"

open(  12, file = 'ey.gnu', position="append" )
if ( iplot .eq. 1 ) rewind(12)
write(12,"(A18,G10.3,A1)")"set title 'Time = ",time,"'"

open(  13, file = 'bz.gnu', position="append" )
if ( iplot .eq. 1 ) rewind(13)
write(13,"(A18,G10.3,A1)")"set title 'Time = ",time,"'"

open(  14, file = 'jx.gnu', position="append" )
if ( iplot .eq. 1 ) rewind(14)
write(14,"(A18,G10.3,A1)")"set title 'Time = ",time,"'"

open(  15, file = 'jy.gnu', position="append" )
if ( iplot .eq. 1 ) rewind(15)
write(15,"(A18,G10.3,A1)")"set title 'Time = ",time,"'"

write(11,*)"splot  'ex_"//fin//"' w l"!, 'ex_"//fin//"' u 1:2:4 w l"
write(11,*)"pause 1"
write(12,*)"splot  'ey_"//fin//"' w l"!, 'ey_"//fin//"' u 1:2:4 w l"
write(12,*)"pause 1"
write(13,*)"splot  'bz_"//fin//"' w l"!, 'bz_"//fin//"' u 1:2:4 w l"
write(13,*)"pause 1"
write(14,*)"splot  'jx_"//fin//"' w l"!, 'jx_"//fin//"' u 1:2:4 w l"
write(14,*)"pause 1"
write(15,*)"splot  'jy_"//fin//"' w l"!, 'jy_"//fin//"' u 1:2:4 w l"
write(15,*)"pause 1"

open( 16, file = 'ex_'//fin )
open( 17, file = 'ey_'//fin )
open( 18, file = 'bz_'//fin )
open( 19, file = 'jx_'//fin )
open( 20, file = 'jy_'//fin )

do i=0,nx-1
   do j=0,ny
      write(16,*) 0.5*(x(i)+x(i+1)), y(j), tm%ex(i,j), c*cos(2*pi*(y(j)+c*time))
      write(19,*) 0.5*(x(i)+x(i+1)), y(j), tm%jx(i,j)
   end do
   write(16,*)
   write(19,*)
end do

do i=0,nx
   do j=0,ny-1
      write(17,*) x(i), 0.5*(y(j)+y(j+1)), tm%ey(i,j) 
      write(20,*) x(i), 0.5*(y(j)+y(j+1)), tm%jy(i,j) 
   end do
   write(17,*)
   write(20,*)
end do

do i=0,nx-1
   do j=0,ny-1
      write(18,*) 0.5*(x(i)+x(i+1)), 0.5*(y(j)+y(j+1)), tm%bz(i,j),  cos(2*pi*(0.5*(y(j)+y(j+1))+c*time))
   end do
   write(18,*)
end do

do k =11,20 
  close(k)
end do

end subroutine plot_champ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine plot_phases( ele, iplot, time )

type(particle) :: ele
integer :: iplot, ipart
real(kind=prec) :: time
real(kind=prec) :: gama, aux, speed
integer :: kk0, kk1, kk2, kk3, kk4, k, l
character(len=4) :: fin
character(len=1) :: aa,bb,cc,dd

kk0 = iplot
kk1 = kk0/1000
aa  = char(kk1 + 48)
kk2 = (kk0 - kk1*1000)/100
bb  = char(kk2 + 48)
kk3 = (kk0 - (kk1*1000) - (kk2*100))/10
cc  = char(kk3 + 48)
kk4 = (kk0 - (kk1*1000) - (kk2*100) - (kk3*10))/1
dd  = char(kk4 + 48)
fin = aa//bb//cc//dd

open( 11, file = 'part_'//nomcas//jname//'.gnu', position="append" )
open( 12, file = 'xvx_'//nomcas//jname//'.gnu', position="append" )
open( 13, file = 'yvy_'//nomcas//jname//'.gnu', position="append" )

if ( iplot .eq. 1 ) then
   rewind(11)
   write(11,*)"set xr[0:",sngl(dimx),"]"
   write(11,*)"set yr[0:",sngl(dimy),"]"
   rewind(12)
   write(12,*)"set xr[0:",sngl(dimx),"]"
   write(12,*)"set yr[-1:1]"
   rewind(13)
   write(13,*)"set xr[0:",sngl(dimx),"]"
   write(13,*)"set yr[-1:1]"
end if

write(11,"(A18,G10.3,A1)")"set title 'Time = ",time,"'"
write(11,*)"plot 'part_"//nomcas//jname//fin//"' w d "
close(11)
write(12,"(A18,G10.3,A1)")"set title 'Time = ",time,"'"
write(12,*)"plot 'part_"//nomcas//jname//fin//"' u 1:3 w d "
close(12)
write(13,"(A18,G10.3,A1)")"set title 'Time = ",time,"'"
write(13,*)"plot 'part_"//nomcas//jname//fin//"' u 2:4 w d "
close(13)

open( 14, file = 'part_'//nomcas//jname//fin )
do ipart=1,nbpart
   speed = sqrt( ele%vit(ipart,1)*ele%vit(ipart,1) + &
        &        ele%vit(ipart,2)*ele%vit(ipart,2) )
   write(14,*) sngl(ele%pos(ipart,1)),sngl(ele%pos(ipart,2))    &
              ,sngl(ele%vit(ipart,1)),sngl(ele%vit(ipart,2)), sngl(speed)
end do
close(14)

end subroutine plot_phases

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Calcul de la fct de distribution en fonction de (vx,vy)

subroutine distribution_v(ele, iplot, time)  

integer :: i, ipart, iplot
real(kind=prec), dimension(:,:), allocatable :: densite
type(particle) :: ele
real(kind=prec) :: speed, time, vx, vy, aux, vth=1.
real(kind=prec) :: pas_v, vmin, vmax
integer :: kk0, kk1, kk2, kk3, kk4, nv=100
character(len=4) :: fin
character(len=1) :: aa,bb,cc,dd

kk0 = iplot
kk1 = kk0/1000
aa  = char(kk1 + 48)
kk2 = (kk0 - kk1*1000)/100
bb  = char(kk2 + 48)
kk3 = (kk0 - (kk1*1000) - (kk2*100))/10
cc  = char(kk3 + 48)
kk4 = (kk0 - (kk1*1000) - (kk2*100) - (kk3*10))/1
dd  = char(kk4 + 48)
fin = aa//bb//cc//dd

allocate(densite(nv,nv))

vmin = -6.d0
vmax = 6.d0
pas_v = (vmax-vmin)/nv

densite = 0.d0
do ipart=1,nbpart
   do i=1,nv
      do j=1,nv
         if (vmin+(i-1)*pas_v <= ele%vit(ipart,1) .and. &
              & ele%vit(ipart,1) < vmin+i*pas_v  .and. & 
              & vmin+(j-1)*pas_v <= ele%vit(ipart,2) .and. &
              & ele%vit(ipart,2) < vmin+j*pas_v) then
            densite(i,j) = densite(i,j) + ele%p(ipart)
         endif
      enddo
   enddo
enddo

open( 27, file = 'densite_v.gnu', position="append" )
if ( iplot .eq. 1 ) rewind(27)
write(27,"(A18,G10.3,A1)")"set title 'Time = ",time,"'"

open( 28, file = 'densite_v_'//nomcas//jname//fin )

write(27,*)"splot  'densite_v_"//nomcas//jname//fin//"' w l, 'densitetheo_"//nomcas//jname//"' w l"
write(27,*)"pause 1"

do i=1,nv
   do j=1,nv
      vx = vmin+(i-0.5)*pas_v
      vy = vmin+(j-0.5)*pas_v
      write(28,*) vx,vy,densite(i,j)/(pas_v*pas_v)
   end do
   write(28,*)
enddo

close(27)
close(28)

if ( iplot .eq. 1 ) then
   open( 37, file = 'densitetheo.gnu' )
   rewind(37)
   write(37,"(A18,G10.3,A1)")"set title 'Time = ",time,"'"

   open( 38, file = 'densitetheo_'//nomcas//jname )

   write(37,*)"splot  'densitetheo_"//nomcas//jname//"' w l"
   write(37,*)"pause 1"

   do i=1,nv
      do j=1,nv
         vx = vmin+(i-0.5)*pas_v
         vy = vmin+(j-0.5)*pas_v
         aux = exp(-(vx*vx+vy*vy)/(2*vth*vth))/(2*pi*vth*vth)
         write(38,*) vx, vy, aux*dimx*dimy 
      end do
      write(38,*)
   enddo
   
   close(37)
   close(38)

endif

end subroutine distribution_v

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Calcul de la fct de distribution en fonction de (x,y)

subroutine distribution_x(ele, iplot, time)  

integer :: i, ipart, iplot
real(kind=prec), dimension(100,100) :: densite
type(particle) :: ele
real(kind=prec) :: time, x, y, pas_x, pas_y
character(len=4) :: fin
character(len=1) :: aa,bb,cc,dd
integer :: kk0, kk1, kk2, kk3, kk4

kk0 = iplot
kk1 = kk0/1000
aa  = char(kk1 + 48)
kk2 = (kk0 - kk1*1000)/100
bb  = char(kk2 + 48)
kk3 = (kk0 - (kk1*1000) - (kk2*100))/10
cc  = char(kk3 + 48)
kk4 = (kk0 - (kk1*1000) - (kk2*100) - (kk3*10))/1
dd  = char(kk4 + 48)
fin = aa//bb//cc//dd

densite = 0.d0
pas_x = dimx/100
pas_y = dimy/100
do ipart=1,nbpart
   do i=1,100
      do j=1,100
         if ((i-1)*pas_x <= ele%pos(ipart,1) .and. &
              & ele%pos(ipart,1) < i*pas_x  .and. & 
              & (j-1)*pas_y <= ele%pos(ipart,2) .and. &
              & ele%pos(ipart,2) < j*pas_y) then
            densite(i,j) = densite(i,j) + ele%p(ipart)
         endif
      enddo
   enddo
enddo

open( 27, file = 'densite_x.gnu', position="append" )
if ( iplot .eq. 1 ) rewind(27)
write(27,"(A18,G10.3,A1)")"set title 'Time = ",time,"'"

open( 28, file = 'densite_x_'//nomcas//jname//fin )

write(27,*)"splot  'densite_x_"//nomcas//jname//fin//"' w l"
write(27,*)"pause 1"

do i=1,100  
   do j=1,100 
      x = (i-0.5)*pas_x
      y = (j-0.5)*pas_y
      write(28,*) x,y,densite(i,j)  
   end do
   write(28,*)
enddo

close(27)
close(28)

end subroutine distribution_x

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine modeE( tm, iplot, time )

type(mesh_fields) :: tm, sol
real(kind=prec) :: time, aux
integer :: iplot

aux =0.d0
do i=0,nx-1
   do j=0,ny-1
      aux = aux + tm%ex(i,j)*tm%ex(i,j)*hx(i)*hhy(j)
   end do
end do
aux = 0.5*log(aux)

open(34,file='modeE.dat',position="append")
if (iplot==1) rewind(34)
write(34,*) time, aux
close(34)

end subroutine modeE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module diagno
