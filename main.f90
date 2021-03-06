program VM_non_unif_2D

use zone
use particules
use initialisation
use villasenor
use maxwell
use diagno

implicit none

type(mesh_fields) :: f0, f1
type(particle) :: p

real(kind=prec) :: time
integer :: istep, icall, iplot
integer :: iargc, n, i
character(len=72) :: argv

n = iargc()
if (n == 0) stop 'Usage: ./VM.exe data file'
do i = 1, n
   call getarg( i, argv); write(*,'(i2, 1x, a)') i, argv
end do

call readin( f0, trim(argv) )

allocate(f0%ex(0:nx-1,0:ny))
allocate(f0%ey(0:nx,0:ny-1))
allocate(f0%bz(0:nx-1,0:ny-1))
allocate(f0%jx(0:nx-1,0:ny))
allocate(f0%jy(0:nx,0:ny-1))
allocate(f0%r0(0:nx,0:ny))  !rho au temps n
allocate(f0%r1(0:nx,0:ny))  !rho au temps n+1
allocate(f1%ex(0:nx,0:ny)) !decales sur maillage de Maxwell
allocate(f1%ey(0:nx,0:ny))
allocate(f1%bz(0:nx,0:ny))
allocate(f1%jx(0:nx,0:ny))
allocate(f1%jy(0:nx,0:ny))

time  = 0.d0
iplot = 0

if( nstep > nstepmax ) nstep = nstepmax

!********************************************************************

istep = 1

call init( f0 )                 !initialisation des champs et densites

do istep = 1, nstep

   if ((nomcas == "faisce") .or. istep == 1) then
      call creapa( p, time ) !creation des particules
   endif
   if (istep > 1) then
      call faraday( f0 ) 	!Calcul de B(n-1/2) --> B(n)			
   end if

   call decalage( f0, f1 )
   call interpol_eb( f1, p )

   write(*,*) p%pos(1,1:2), p%vit(1,1:2)
   call avancee_vitesse( p )

   if (jname == 'jcico1') then
      call avancee_part( p, 0.5d0 )  ! x(n) --> x(n+1/2)
      call sortie_part( p )
      call calcul_j_cic( p, f0, f1 )
      call avancee_part( p, 0.5d0 )  ! x(n+1/2) -- x(n+1)
      call sortie_part( p )
   else if (jname == 'jcoco1') then
      call avancee_part( p, 1.d0 )
      call calcul_j_villa( p, f0 )
      call sortie_part( p )
   else
      call avancee_part( p, 1.d0 )
      call sortie_part( p )
   end if
        
   !call calcul_rho( p, f0 )

   call faraday( f0 )   !Calcul de B(n) --> B(n+1/2)
   call ampere( f0 )    !Calcul de E(n) --> E(n+1)
   call conditions_limites( f0, time )

   time = time + dt
   print*,'time = ',time, ' nbpart = ', nbpart

   if ( istep==1 .or. mod(istep,idiag) == 0 .or. istep==nstep ) then
      iplot = iplot + 1
     ! call diag_coc( f0, p, time, iplot )
     ! call diag_champ_part( p, time, iplot )
     ! call plot_champ( f0, iplot, time )
     ! call plot_phases( p, iplot, time )
     ! call distribution_v( p, iplot, time )  
     ! call distribution_x( p, iplot, time )
      call modeE( f0, iplot, time )
   endif

end do

end program VM_non_unif_2D
