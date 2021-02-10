module pstd

use, intrinsic :: iso_c_binding

implicit none

include 'fftw3.f03'

type :: maxwell_pstd

   private
   integer          :: nc_eta1      !< x cells number
   integer          :: nc_eta2      !< y cells number
   real(8)          :: eta1_min     !< left side
   real(8)          :: eta1_max     !< right side
   real(8)          :: delta_eta1   !< step size
   real(8)          :: eta2_min     !< bottom side
   real(8)          :: eta2_max     !< top side
   real(8)          :: delta_eta2   !< step size
   real(8), pointer :: d_dx(:)      !< field x derivative
   real(8), pointer :: d_dy(:)      !< field y derivative
   real(8), pointer :: kx(:)        !< x wave number
   real(8), pointer :: ky(:)        !< y wave number
   integer(8)       :: fwx          !< forward fft plan along x
   integer(8)       :: fwy          !< forward fft plan along y
   integer(8)       :: bwx          !< backward fft plan along x
   integer(8)       :: bwy          !< backward fft plan along y
   complex(8), pointer :: tmp_x(:) !< x fft transform
   complex(8), pointer :: tmp_y(:) !< y fft transform

end type maxwell_pstd


contains

!> Initialize 2d maxwell solver on cartesian mesh with PSTD scheme
subroutine init_pstd(self, dimx, nx, dimy, ny)

   type(maxwell_pstd) :: self         !< maxwell object
   real(8), intent(in)      :: dimx         !< x min
   real(8), intent(in)      :: dimy         !< x max
   integer,  intent(in)     :: nx         !< x cells number
   integer,  intent(in)     :: ny         !< y cells number

   real(8)                  :: dx           !< x space step
   real(8)                  :: dy           !< y space step
   real(8)                  :: kx0
   real(8)                  :: ky0
   integer                  :: i, j
   real(8)                  :: pi

   self%nc_eta1 = nx
   self%nc_eta2 = ny
   self%eta1_min = 0d0
   self%eta2_min = 0d0
   self%eta1_max = dimx
   self%eta2_max = dimy

   pi = 4 * atan(1d0)

   allocate(self%d_dx(nx))
   allocate(self%d_dy(ny))

   allocate(self%tmp_x(nx/2+1))
   allocate(self%tmp_y(ny/2+1))

   call dfftw_plan_dft_r2c_1d(self%fwx, nx, self%d_dx,  self%tmp_x, FFTW_ESTIMATE)
   call dfftw_plan_dft_c2r_1d(self%bwx, nx, self%tmp_x, self%d_dx,  FFTW_ESTIMATE)
   call dfftw_plan_dft_r2c_1d(self%fwy, ny, self%d_dy,  self%tmp_y, FFTW_ESTIMATE)
   call dfftw_plan_dft_c2r_1d(self%bwy, ny, self%tmp_y, self%d_dy,  FFTW_ESTIMATE)

   allocate(self%kx(nx/2+1))
   allocate(self%ky(ny/2+1))

   dx = dimx / nx
   dy = dimy / ny

   kx0 = 2._8*pi/(nx*dx)
   ky0 = 2._8*pi/(ny*dy)

   do i=2,nx/2+1
      self%kx(i) = (i-1)*kx0
   end do
   self%kx(1) = 1.0_8
   do j=2,ny/2+1
      self%ky(j) = (j-1)*ky0
   end do
   self%ky(1) = 1.0_8

end subroutine init_pstd

!> Solve faraday equation (ex,ey,hz)
subroutine faraday_pstd(self, ex, ey, hz, dt)

   type(maxwell_pstd),intent(inout) :: self   !< maxwell object

   real(8), dimension(:,:), intent(inout) :: ex     !< electric field x
   real(8), dimension(:,:), intent(inout) :: ey     !< electric field y
   real(8), dimension(:,:), intent(inout) :: hz     !< magnetic field z
   integer                                :: nx   !< x cells number
   integer                                :: ny   !< y cells number
   real(8), intent(in)                    :: dt     !< time step
   integer                                :: i, j

   nx = self%nc_eta1
   ny = self%nc_eta2

   do i = 1, nx
      self%d_dy = ex(i,1:ny)
      call dfftw_execute_dft_r2c(self%fwy, self%d_dy, self%tmp_y)
      self%tmp_y(2:ny/2+1)=-cmplx(0d0,self%ky(2:ny/2+1),kind=8)*self%tmp_y(2:ny/2+1)
      call dfftw_execute_dft_c2r(self%bwy, self%tmp_y, self%d_dy)
      self%d_dy = self%d_dy / ny
      hz(i,1:ny) = hz(i,1:ny) + dt * self%d_dy
   end do

   do j = 1, ny
      self%d_dx = ey(1:nx,j)
      call dfftw_execute_dft_r2c(self%fwx, self%d_dx, self%tmp_x)
      self%tmp_x(2:nx/2+1)=-cmplx(0d0,self%kx(2:nx/2+1),kind=8)*self%tmp_x(2:nx/2+1)
      call dfftw_execute_dft_c2r(self%bwx, self%tmp_x, self%d_dx)
      self%d_dx = self%d_dx / nx
      hz(1:nx,j) = hz(1:nx,j) - dt * self%d_dx
   end do

   if (size(hz,1) == nx+1) hz(nx+1,:) = hz(1,:)
   if (size(hz,2) == ny+1) hz(:,ny+1) = hz(:,1)

end subroutine faraday_pstd

subroutine ampere_pstd(self, ex, ey, hz, dt, jx, jy)

   type(maxwell_pstd),intent(inout) :: self   !< maxwell equation
   real(8), dimension(:,:)                :: ex     !< electric field x
   real(8), dimension(:,:)                :: ey     !< electric field y
   real(8), dimension(:,:)                :: hz     !< magnetic field z
   real(8)                                :: dt     !< time step
   real(8), dimension(:,:), optional      :: jx     !< current x
   real(8), dimension(:,:), optional      :: jy     !< current y

   integer                                :: nx
   integer                                :: ny
   integer                                :: i, j

   nx = self%nc_eta1
   ny = self%nc_eta2

   do j = 1, ny
      self%d_dx = hz(1:nx,j)
      call dfftw_execute_dft_r2c(self%fwx, self%d_dx, self%tmp_x)
      self%tmp_x(2:nx/2+1)=-cmplx(0d0,self%kx(2:nx/2+1),kind=8)*self%tmp_x(2:nx/2+1)
      call dfftw_execute_dft_c2r(self%bwx, self%tmp_x, self%d_dx)
      self%d_dx = self%d_dx / nx
      ey(1:nx,j) = ey(1:nx,j) - dt * self%d_dx
   end do

   do i = 1, nx
      self%d_dy = hz(i,1:ny)
      call dfftw_execute_dft_r2c(self%fwy, self%d_dy, self%tmp_y)
      self%tmp_y(2:ny/2+1)=-cmplx(0d0,self%ky(2:ny/2+1),kind=8)*self%tmp_y(2:ny/2+1)
      call dfftw_execute_dft_c2r(self%bwy, self%tmp_y, self%d_dy)
      self%d_dy = self%d_dy / ny
      ex(i,1:ny) = ex(i,1:ny) + dt * self%d_dy
   end do

   If (present(jx) .and. present(jy)) then
      ex = ex - dt * jx
      ey = ey - dt * jy
   end if

   if (size(ex,1) == nx+1) ex(nx+1,:) = ex(1,:)
   if (size(ey,2) == ny+1) ey(:,ny+1) = ey(:,1)

end subroutine ampere_pstd

end module pstd
