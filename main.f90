program VM_non_unif_2D

use zone
use particules
use maxwell

implicit none

type(mesh_fields) :: f0, f1
type(particle) :: p

real(kind=prec) :: time
integer :: istep
integer :: iargc, n, i, j
character(len=72) :: argv

n = iargc()
if (n == 0) stop 'Usage: ./VM.exe data file'
do i = 1, n
   call getarg( i, argv); write(*,'(i2, 1x, a)') i, argv
end do

call readin( trim(argv) )

allocate(f0%ex(0:nx-1,0:ny))
allocate(f0%ey(0:nx,0:ny-1))
allocate(f0%bz(0:nx-1,0:ny-1))
allocate(f0%jx(0:nx-1,0:ny))
allocate(f0%jy(0:nx,0:ny-1))
allocate(f1%ex(0:nx,0:ny)) !decales sur maillage de Maxwell
allocate(f1%ey(0:nx,0:ny))
allocate(f1%bz(0:nx,0:ny))
allocate(f1%jx(0:nx,0:ny))
allocate(f1%jy(0:nx,0:ny))

time  = 0.d0

if( nstep > nstepmax ) nstep = nstepmax

!********************************************************************

istep = 1

f0%ex = 0.d0; f0%ey = 0.d0; f0%bz = 0.d0
do i=0,nx-1
   do j=0,ny
      f0%ex(i,j) = alpha/kx * sin(kx*x(i))
   enddo
enddo

call plasma( p ) !creation des particules

do istep = 1, nstep

   if (istep > 1) call faraday( f0 )

   call decalage( f0, f1 )
   call interpol_eb( f1, p )

   call avancee_vitesse( p )

   call avancee_part( p, 0.5d0 )  ! x(n) --> x(n+1/2)
   call calcul_j_cic( p, f0, f1 )
   call avancee_part( p, 0.5d0 )  ! x(n+1/2) -- x(n+1)
        
   call faraday( f0 )   !Calcul de B(n) --> B(n+1/2)
   call ampere( f0 )    !Calcul de E(n) --> E(n+1)

   time = time + dt
   print*,'time = ',time, ' nbpart = ', nbpart

   call modeE( f0, istep, time )

end do


contains

subroutine modeE( tm, iplot, time )

type(mesh_fields) :: tm
real(kind=prec) :: time, aux
integer :: iplot, i, j

aux =0.d0
do i=0,nx-1
   do j=0,ny-1
      aux = aux + tm%ex(i,j)*tm%ex(i,j)*dx*dy
   end do
end do
aux = 0.5*log(aux)

open(34,file='modeE.dat',position="append")
if (iplot==1) rewind(34)
write(34,*) time, aux
close(34)

end subroutine modeE

end program VM_non_unif_2D
