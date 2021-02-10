program VM_non_unif_2D

use zone
use particules
use pstd

implicit none

type(mesh_fields) :: f1
type(particle) :: p
type(maxwell_pstd) :: solver

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

allocate(f1%ex(0:nx,0:ny)) !decales sur maillage de Maxwell
allocate(f1%ey(0:nx,0:ny))
allocate(f1%bz(0:nx,0:ny))
allocate(f1%jx(0:nx,0:ny))
allocate(f1%jy(0:nx,0:ny))

time  = 0.d0

if( nstep > nstepmax ) nstep = nstepmax

!********************************************************************

istep = 1

f1%ex = 0.d0; f1%ey = 0.d0; f1%bz = 0.d0
do i=0,nx
   do j=0,ny
      f1%ex(i,j) = alpha/kx * sin(kx*x(i))
   enddo
enddo

call init_pstd( solver, dimx, nx, dimy, ny)

call plasma( p ) !creation des particules
call calcul_j_cic( p, f1 )
dt = 0.001

do istep = 1, nstep

   if (istep > 1) call faraday_pstd( solver, f1%ex, f1%ey, f1%bz, 0.5*dt )

   call interpol_eb( f1, p )
   call avancee_vitesse( p )

   call avancee_part( p, 0.5d0 )  ! x(n) --> x(n+1/2)
   call calcul_j_cic( p, f1 )
   call avancee_part( p, 0.5d0 )  ! x(n+1/2) -- x(n+1)
        
   call faraday_pstd( solver, f1%ex, f1%ey, f1%bz, 0.5*dt )
   call ampere_pstd( solver, f1%ex, f1%ey, f1%bz, dt, f1%jx, f1%jy )

   time = time + dt
   print*,'time = ',time, ' nbpart = ', nbpart

   call modeE( f1, istep, time )

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
