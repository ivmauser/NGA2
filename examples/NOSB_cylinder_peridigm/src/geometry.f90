!> Various definitions and tools for initializing NGA2 config
module geometry
   use ibconfig_class,   only: ibconfig
   use precision,      only: WP
   implicit none
   private
   
   !> Single config
   type(ibconfig), public :: cfg

   public :: geometry_init
   
contains
   
   
   !> Initialization of problem geometry
   subroutine geometry_init
      use sgrid_class, only: sgrid
      use param,       only: param_read
      implicit none
      type(sgrid) :: grid
      real(WP) :: Ly,elem
      
      
      ! Create a grid from input params
      create_grid: block
         use sgrid_class, only: cartesian
         integer :: i,j,k,nx,ny,nz,N_r,N_xp,N_xm,N_y,N_z,npad
         real(WP) :: R,Lx,Lz
         real(WP), dimension(:), allocatable :: x,y,z
         
         npad = 1
         ! Read in grid definition
         call param_read('R',R,default=0.5_WP)
         call param_read('N_r',N_r,default=10)
         call param_read('Lx',Lx,default=2.5_WP)      
         call param_read('Ly',Ly,default=0.41_WP)     
         call param_read('Lz',Lz,default=0.41_WP)

         ! Use the same spacing in every direction. 
         elem=R/real(N_r,WP)
         nx = int(Lx/elem); allocate(x(nx+1))
         ny = int(Ly/elem)+ 2*npad; allocate(y(ny+1)) 
         
         ! Option to make 2D in the z direction
         if (Lz.eq.0.0_WP) then
            Lz = elem/3.0_WP; nz = 1; allocate(z(nz+1))
         else
            nz = int(Lz/elem); allocate(z(nz+1))
         end if

         ! Create simple rectilinear grid
         do i=1,nx+1
            x(i)=real(i-1,WP)*elem
         end do
         do j=1,ny+1
            y(j)=real(j-1-npad,WP)*elem ! -2 to account of the extra cell which would make up the wall
         end do
         if (nz.eq.1) then
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
         cfg=ibconfig(grp=group,decomp=partition,grid=grid)
      end block create_cfg
      
      
      ! Create IB walls for this config
      create_walls: block
      use ibconfig_class, only: bigot,sharp
      integer :: i,j,k
      ! Create IB field
      do k=cfg%kmino_,cfg%kmaxo_
         do j=cfg%jmino_,cfg%jmaxo_
            do i=cfg%imino_,cfg%imaxo_
               cfg%Gib(i,j,k)=sqrt((cfg%ym(j) - 0.5_WP*Ly)**2)- 0.5_WP*Ly
            end do
         end do
      end do
      ! Get normal vector
      call cfg%calculate_normal()
      ! Get VF field
      call cfg%calculate_vf(method=sharp,allow_zero_vf=.false.)
    end block create_walls
      
      
   end subroutine geometry_init
   
   
end module geometry
