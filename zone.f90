module zone

integer, parameter :: prec=8

type mesh_fields
   real(kind=prec), dimension(:,:), pointer :: ex, ey
   real(kind=prec), dimension(:,:), pointer :: bz
   real(kind=prec), dimension(:,:), pointer :: jx, jy
end type mesh_fields

type particle
   real(kind=prec)   , pointer :: pos(:,:)
   integer           , pointer :: case(:,:)
   real(kind=prec)   , pointer :: vit(:,:)
   real(kind=prec)   , pointer :: epx(:)
   real(kind=prec)   , pointer :: epy(:)
   real(kind=prec)   , pointer :: bpz(:)
end type particle

logical :: relativ = .false.

real(kind=prec) :: pi 

character(len=6) :: nomcas
character(len=6) :: bcname 
character(len=6) :: jname


integer :: nx, ny
integer :: nstep, nstepmax
integer :: icrea, idiag
integer :: nbpart

integer, private :: i, j

real(kind=prec) :: dt, alpha, kx, ky, c, csq, e0
real(kind=prec) :: dx, dy, dx1, dy1, dx2, dy2
real(kind=prec), dimension(:), pointer :: x, y
real(kind=prec) :: dimx, dimy, dimx1, dimy1, dimx2, dimy2
real(kind=prec) :: cfl
real(kind=prec) :: tfinal

real(kind=prec) :: exext, eyext, bzext

real(kind=prec) :: charge, masse, poids
real(kind=prec) :: q_sur_m 


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine readin( filename )

implicit none

character(len=*) :: filename

namelist/donnees/ dimx,  &      !dimensions du domaine
                  dimy,  & 
                    nx,  &      !nbre de pas
                    ny,  &
                   cfl,  &      !nbre de Courant
                tfinal,  &      !duree maxi
         nstepmax,  &   !nbre d'iterations maxi
                nomcas,  &      !nom du cas ce calcul
                 jname,  &      !calcul de j   
       icrea,  &   !frequence d'emission des particules
       idiag,  &   !frequence des diagnostics
       bcname, &    !type de conditions limites
       exext,    &   !champ electrique exterieur
       eyext,    &   !champ electrique exterieur
       bzext,  &   !champ magnetique exterieur
                 charge, &   !charge d'une macroparticule
                 masse,  &      !masse d'une macroparticule
                     c,  &      !vitesse de la lumiere
                    e0,  &      !permittivite du vide
                 relativ        !calcul relativiste de la vitesse

!***Initialisation  des valeurs pas default

pi = 4. * atan(1.)

open(93,file=filename,status='old')
read(93,donnees) 
close(93)

csq = c*c
q_sur_m = charge / masse
poids = charge

alpha = 0.1
kx = 0.5
ky = 0.
dimx = 2*pi/kx
poids = dimx * dimy ! car int(f0) = dimx*dimy

!Creation du maillage

allocate(x(-1:nx+1))  !0:nx))
allocate(y(-1:ny+1))  !0:ny))

dx = dimx / nx
dy = dimy / ny

x(0) = 0.
y(0) = 0.

do i=1,nx
   x(i) = i*dx 
enddo
do j=1,ny
   y(j) = j*dy
enddo

x(-1)   = x(0) - dx
x(nx+1) = x(nx) + dx
y(-1)   = y(0) - dy
y(ny+1) = y(ny) + dy

dt    = cfl  / sqrt (1./(dx*dx)+1./(dy*dy)) / c

nstep = floor(tfinal/dt)

write(*,*) " cfl = ", cfl
write(*,*) " dx = ", dx, " dy = ", dy, " dt = ", dt
if( nstep > nstepmax ) nstep = nstepmax
write(*,*) " Nombre d'iteration nstep = ", nstep
write(*,*)

end subroutine readin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module zone
