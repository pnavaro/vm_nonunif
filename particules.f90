module particules

use zone
use quietstart

implicit none

integer, private :: ipart 
integer, private :: i, j, k
integer, private :: ivarx, ivary
real(kind=prec), private :: varx, vary

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine creapa( ele, time )

type  (particle) :: ele
real(kind=prec) :: time

call plasma( ele, time )

end subroutine creapa

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine interpol_eb( tm1, ele )

type(particle) :: ele
type(mesh_fields) :: tm1
real(kind = prec) :: a1, a2, a3, a4
real(kind = prec) :: xp, yp, dum
!   ______________
!  |     |        |
!  | a2  |  a1    |
!  |_____|________|
!  |     |        |
!  | a3  |  a4    |
!  |     |        |
!  |_____|________|

do ipart=1,nbpart
   i = ele%case(ipart,1)
   j = ele%case(ipart,2)
   xp = ele%pos(ipart,1)
   yp = ele%pos(ipart,2)

   dum = 1./(hx(i)*hy(j))
   a1 = (x(i+1)-xp) * (y(j+1)-yp) * dum
   a2 = (xp-x(i)) * (y(j+1)-yp) * dum
   a3 = (xp-x(i)) * (yp-y(j)) * dum
   a4 = (x(i+1)-xp) * (yp-y(j)) * dum

   ele%epx(ipart) = a1 * tm1%ex(i,j) + a2 * tm1%ex(i+1,j) &
        & + a3 * tm1%ex(i+1,j+1) + a4 * tm1%ex(i,j+1) 
   ele%epy(ipart) = a1 * tm1%ey(i,j) + a2 * tm1%ey(i+1,j) &
        & + a3 * tm1%ey(i+1,j+1) + a4 * tm1%ey(i,j+1) 
   ele%bpz(ipart) =  a1 * tm1%bz(i,j) + a2 * tm1%bz(i+1,j) &
        & + a3 * tm1%bz(i+1,j+1) + a4 * tm1%bz(i,j+1) 
end do

end subroutine interpol_eb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine avancee_vitesse( ele )

type (particle) :: ele
real(kind=prec) :: dum, u2
real(kind=prec) :: tantheta, sintheta
real(kind=prec) :: xpp, ypp
real(kind=prec) :: gamma

do ipart = 1, nbpart

   !*** Changement de variable u = gamma*vit

   if( relativ ) then

      u2    =   ele%vit(ipart,1)*ele%vit(ipart,1)	&
           + ele%vit(ipart,2)*ele%vit(ipart,2)
      if ( u2 >= csq ) then 
         print*,'Erreur : u2 >= c2 dans le calcul de la vitesse'
         print*,'ipart = ',ipart,' vx = ',ele%vit(ipart,1),' vy = ',ele%vit(ipart,2)
         stop
      else
         gamma = 1./sqrt( 1. - u2/csq )
      endif

      ele%vit(ipart,1) = gamma*ele%vit(ipart,1)
      ele%vit(ipart,2) = gamma*ele%vit(ipart,2)

   else

      gamma=1.

   end if


   !*** Separation des effets electriques et magnetiques

   !*** On ajoute la moitie de l'effet champ electrique E

   dum = 0.5 * dt * q_sur_m
   ele%vit(ipart,1) = ele%vit(ipart,1) + dum*(ele%epx(ipart)+exext)
   ele%vit(ipart,2) = ele%vit(ipart,2) + dum*(ele%epy(ipart)+eyext)

   !*** Algorithme de Buneman pour les effets magnetiques
 
   tantheta = dum * (ele%bpz(ipart)+bzext) / gamma 
   sintheta = 2.0 * tantheta / ( 1. + tantheta*tantheta)

   ele%vit(ipart,1) = ele%vit(ipart,1) + ele%vit(ipart,2)*tantheta
   ele%vit(ipart,2) = ele%vit(ipart,2) - ele%vit(ipart,1)*sintheta
   ele%vit(ipart,1) = ele%vit(ipart,1) + ele%vit(ipart,2)*tantheta

   !*** Autre moitie de l'effet du champ electrique E

   ele%vit(ipart,1) = ele%vit(ipart,1) + dum*(ele%epx(ipart)+exext)
   ele%vit(ipart,2) = ele%vit(ipart,2) + dum*(ele%epy(ipart)+eyext)

   !*** On repasse a la vitesse (changement de variable inverse)

   if( relativ ) then

      u2 =   ele%vit(ipart,1)*ele%vit(ipart,1)	&
           + ele%vit(ipart,2)*ele%vit(ipart,2)

      gamma = sqrt( 1. + u2/csq )
 
      ele%vit(ipart,1) = ele%vit(ipart,1) / gamma
      ele%vit(ipart,2) = ele%vit(ipart,2) / gamma

   end if

end do

end subroutine avancee_vitesse

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine avancee_part( ele, coef )  !Avancee de coef * dt

type(particle) :: ele
real(kind=prec) :: coef

do ipart=1,nbpart     
   ele%pos(ipart,1) = ele%pos(ipart,1) + ele%vit(ipart,1)*dt*coef
   ele%pos(ipart,2) = ele%pos(ipart,2) + ele%vit(ipart,2)*dt*coef
enddo

!*** Mise a jour des "cases"

do ipart=1,nbpart
   i = 0
   do while (ele%pos(ipart,1) >= x(i) .and. ele%pos(ipart,1)<dimx) 
      i=i+1
   enddo
   if ( ele%pos(ipart,1) >= dimx ) i=nx+1
   ele%case(ipart,1) = i-1
   j = 0
   do while (ele%pos(ipart,2) >= y(j) .and. ele%pos(ipart,2)<dimy) 
      j=j+1 
   enddo
   if ( ele%pos(ipart,2) >= dimy ) j=ny+1
   ele%case(ipart,2) = j-1
end do

end subroutine avancee_part

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sortie_part( ele )

type(particle) :: ele

!*** Traitement de la sortie des particules

if (bcname == 'period') then

   do ipart=1,nbpart
      if( ele%pos(ipart,1) >= dimx ) ele%pos(ipart,1) = ele%pos(ipart,1) - dimx
      if( ele%pos(ipart,2) >= dimy ) ele%pos(ipart,2) = ele%pos(ipart,2) - dimy
      if( ele%pos(ipart,1) < 0.d0 )  ele%pos(ipart,1) = ele%pos(ipart,1) + dimx
      if( ele%pos(ipart,2) < 0.d0 )  ele%pos(ipart,2) = ele%pos(ipart,2) + dimy
   end do   

else
   ipart = 1
   do while ( ipart <= nbpart )  
      if (      ele%pos(ipart,1) < 0.0	        &
           .or. ele%pos(ipart,1) >= dimx     	&	
           .or. ele%pos(ipart,2) < 0.         	&
           .or. ele%pos(ipart,2) >= dimy ) then	
         !*** Recuperation du trou laisse par la particule sortie
         ele%pos(ipart,1) = ele%pos(nbpart,1)
         ele%pos(ipart,2) = ele%pos(nbpart,2)
         ele%vit(ipart,1) = ele%vit(nbpart,1)
         ele%vit(ipart,2) = ele%vit(nbpart,2)
         ele%p(ipart)     = ele%p(nbpart)
	 nbpart = nbpart - 1
         if (nbpart == 0) stop 'plus de particule'
      else
         ipart = ipart + 1
      end if
   end do
end if


!*** Mise a jour des "cases"

do ipart=1,nbpart
   i = 0
   do while (ele%pos(ipart,1) >= x(i) .and. ele%pos(ipart,1)<=dimx) 
      i=i+1
   enddo
   ele%case(ipart,1) = i-1
   j = 0
   do while (ele%pos(ipart,2) >= y(j) .and. ele%pos(ipart,2)<=dimy) 
      j=j+1 
   enddo
   ele%case(ipart,2) = j-1
end do

end subroutine sortie_part

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calcul_rho( ele, tm )

type(particle) :: ele
type(mesh_fields) :: tm
real(kind=prec) :: a1, a2, a3, a4, dum, xp, yp
real(kind=prec) :: rho_total

tm%r0 = tm%r1   
tm%r1 = 0.d0    
                
!   ______________
!  |     |        |
!  | a2  |  a1    |
!  |_____|________|
!  |     |        |
!  | a3  |  a4    |
!  |     |        |
!  |_____|________|

do ipart=1,nbpart
   i = ele%case(ipart,1)
   j = ele%case(ipart,2)
   xp = ele%pos(ipart,1)
   yp = ele%pos(ipart,2)
   dum = ele%p(ipart) / (hx(i)*hy(j))
   a1 = (x(i+1)-xp) * (y(j+1)-yp) * dum
   a2 = (xp-x(i)) * (y(j+1)-yp) * dum
   a3 = (xp-x(i)) * (yp-y(j)) * dum
   a4 = (x(i+1)-xp) * (yp-y(j)) * dum
   tm%r1(i,j)     = tm%r1(i,j)     + a1/(hhx(i)*hhy(j)) !charge unite = 1
   tm%r1(i+1,j)   = tm%r1(i+1,j)   + a2/(hhx(i+1)*hhy(j)) 
   tm%r1(i+1,j+1) = tm%r1(i+1,j+1) + a3/(hhx(i+1)*hhy(j+1))
   tm%r1(i,j+1)   = tm%r1(i,j+1)   + a4/(hhx(i)*hhy(j+1))
end do

if (bcname == 'period') then
   do i=0,nx
      tm%r1(i,0)  = tm%r1(i,0) + tm%r1(i,ny)
      tm%r1(i,ny) = tm%r1(i,0)
   end do
   do j=0,ny
      tm%r1(0,j)  = tm%r1(0,j) + tm%r1(nx,j)
      tm%r1(nx,j) = tm%r1(0,j)
   end do
end if

if (nomcas == 'gaussv' .or. nomcas == 'plasma' ) then
   rho_total = 0.d0
   do i=0,nx-1
      do j=0,ny-1
         rho_total = rho_total + tm%r1(i,j)*hhx(i)*hhy(j)
      enddo
   enddo
   print*,'rho total',rho_total 
   ! Neutralisation du milieu
   tm%r1 = tm%r1 - rho_total/dimx/dimy
endif

end subroutine calcul_rho

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calcul_j_cic( ele, tm, tm1 )

type(particle) :: ele
type(mesh_fields) :: tm, tm1
real(kind=prec) :: a1, a2, a3, a4, dum, xp, yp

tm1%jx = 0.d0
tm1%jy = 0.d0

do ipart=1,nbpart
   i = ele%case(ipart,1)
   j = ele%case(ipart,2)
   xp = ele%pos(ipart,1)
   yp = ele%pos(ipart,2)
   dum = ele%p(ipart) / (hx(i)*hy(j))
   a1 = (x(i+1)-xp) * (y(j+1)-yp) * dum
   a2 = (xp-x(i)) * (y(j+1)-yp) * dum
   a3 = (xp-x(i)) * (yp-y(j)) * dum
   a4 = (x(i+1)-xp) * (yp-y(j)) * dum
   dum = ele%vit(ipart,1) / (hx(i)*hy(j)) !charge unite = 1
   tm1%jx(i,j)     = tm1%jx(i,j)     + a1*dum  
   tm1%jx(i+1,j)   = tm1%jx(i+1,j)   + a2*dum 
   tm1%jx(i+1,j+1) = tm1%jx(i+1,j+1) + a3*dum 
   tm1%jx(i,j+1)   = tm1%jx(i,j+1)   + a4*dum 
   dum = ele%vit(ipart,2) / (hx(i)*hy(j)) 
   tm1%jy(i,j)     = tm1%jy(i,j)     + a1*dum  
   tm1%jy(i+1,j)   = tm1%jy(i+1,j)   + a2*dum 
   tm1%jy(i+1,j+1) = tm1%jy(i+1,j+1) + a3*dum 
   tm1%jy(i,j+1)   = tm1%jy(i,j+1)   + a4*dum 
end do

if (bcname == 'period') then
   do i=0,nx
      tm1%jx(i,0)  = tm1%jx(i,0) + tm1%jx(i,ny)
      tm1%jx(i,ny) = tm1%jx(i,0)
      tm1%jy(i,0)  = tm1%jy(i,0) + tm1%jy(i,ny)
      tm1%jy(i,ny) = tm1%jy(i,0)
   end do
   do j=0,ny
      tm1%jx(0,j)  = tm1%jx(0,j) + tm1%jx(nx,j)
      tm1%jx(nx,j) = tm1%jx(0,j)
      tm1%jy(0,j)  = tm1%jy(0,j) + tm1%jy(nx,j)
      tm1%jy(nx,j) = tm1%jy(0,j)
   end do
end if


do i=0,nx-1
do j=0,ny
   tm%jx(i,j) = 0.5 * (tm1%jx(i,j)+tm1%jx(i+1,j))
end do
end do

do i=0,nx
do j=0,ny-1
   tm%jy(i,j) = 0.5 * (tm1%jy(i,j)+tm1%jy(i,j+1))
end do
end do

end subroutine calcul_j_cic

subroutine plasma( ele, time )

type (particle) :: ele
real(kind=prec) :: time, speed, theta, vth, n, aux
real(kind=prec) :: a, b, eps, R
integer :: k

if( time == 0 ) then

   eps = 1.d-12

   vth =  1.
   nbpart = 100*(nx)*(ny)
   n = 1.d0/nbpart

   allocate(ele%pos(nbpart,2))
   allocate(ele%case(nbpart,2))
   allocate(ele%vit(nbpart,2))
   allocate(ele%epx(nbpart))
   allocate(ele%epy(nbpart))
   allocate(ele%bpz(nbpart))
   allocate(ele%p(nbpart))

   do k=0,nbpart-1

      speed = vth * sqrt(-2 * log( (k+0.5)*n ))

      theta = trinary_reversing( k ) * 2 * pi

      a = 0; b = dimx ! 2*pi/kx 
      R = bit_reversing( k )
      call dichotomie_x(a,b,R,eps) 
      ele%pos(k+1,1) = a
      ele%pos(k+1,2) = dimy * penta_reversing( k ) 

      i = 0
      do while (ele%pos(k+1,1) >= x(i)) 
         i=i+1
      enddo
      ele%case(k+1,1) = i-1
      
      j = 0
      do while (ele%pos(k+1,2) >= y(j)) 
         j=j+1 
      enddo
      ele%case(k+1,2) = j-1
      
      ele%vit(k+1,1) = speed * cos(theta)  !
      ele%vit(k+1,2) = speed * sin(theta)  !

      ele%p(k+1) = poids * n

   enddo


end if

end subroutine plasma

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module particules
