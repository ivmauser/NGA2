!> Various definitions and tools for initializing NGA2 config
module geometry
   use config_class,   only: config
   use precision,      only: WP
   implicit none
   private
   
   !> Single config
   type(config), public :: cfg

   public :: geometry_init
   
contains
   
   
   !> Initialization of problem geometry
   subroutine geometry_init
      use sgrid_class, only: sgrid
      use param,       only: param_read
      implicit none
      type(sgrid) :: grid
      
      
      ! Create a grid from input params
      create_grid: block
         use sgrid_class, only: cartesian
         integer :: i,j,k,nx,ny,nz,N_r,N_xp,N_xm,N_y,N_z
         real(WP) :: R,elem,Lx,Ly,Lz
         real(WP), dimension(:), allocatable :: x,y,z
         
         ! Read in grid definition
         call param_read('R',R,default=0.5_WP)
         call param_read('N_r',N_r,default=10)
         call param_read('X+ ratio',N_xp,default=6)
         call param_read('X- ratio',N_xm,default=6)
         call param_read('Y ratio',N_y,default=6)
         call param_read('Z ratio',N_z,default=6)

         ! Use the same spacing in every direction.  The cylinder center is
         ! located at x=0, with N_xm radii upstream and N_xp radii downstream.
         elem=R/real(N_r,WP)
         Lx=real(N_xm+N_xp,WP)*R
         Ly=real(N_y,WP)*R
         
         
         nx = N_r*(N_xp + N_xm); allocate(x(nx+1))
         ny = N_r*N_y;           allocate(y(ny+1))
         
         ! Option to make 2D in the z direction
         if (N_z.eq.0) then
            Lz = elem/3.0_WP; nz = 1; allocate(z(nz+1))
         else
            Lz=real(N_z,WP)*R; nz = N_r*N_z; allocate(z(nz+1))
         end if

         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=-real(N_xm,WP)*R+real(i-1,WP)*elem
         end do
         do j=1,ny+1
            y(j)=-0.5_WP*Ly+real(j-1,WP)*elem
         end do
         if (N_z.eq.0) then
            do k=1,nz+1
               z(k)=-0.5_WP*Lz+real(k-1,WP)*elem/3.0_WP
            end do
         else
            do k=1,nz+1
               z(k)=-0.5_WP*Lz+real(k-1,WP)*elem
            end do
         end if
         
         
         ! General serial grid object (no=3 needed to support ghost/image point interpolation/extrapolation)
         grid=sgrid(coord=cartesian,no=2,x=x,y=y,z=z,xper=.false.,yper=.true.,zper=.true.,name='box')
         
      end block create_grid
      
      
      ! Create a config from that grid on our entire group
      create_cfg: block
         use parallel, only: group
         integer, dimension(3) :: partition
         ! Read in partition
         call param_read('Partition',partition,short='p')
         ! Create partitioned grid
         cfg=config(grp=group,decomp=partition,grid=grid)
      end block create_cfg
      
      
      ! Create walls for this config
      create_walls: block
        cfg%VF=1.0_WP
      end block create_walls
      
      
   end subroutine geometry_init
   
   
end module geometry
