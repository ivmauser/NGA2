!> Lagrangian solid solver object
!> Attempt at integrating the pdsolver_class without AMR
module lsspd_class
   use precision,      only: WP, I8
   use string,         only: str_medium
   use config_class,   only: config
   use ddadi_class,    only: ddadi
   use mpi_f08,        only: MPI_Datatype,MPI_INTEGER8,MPI_INTEGER,MPI_DOUBLE_PRECISION
   use NOSB_class,     only: pdsolver, PDC_IS_DEAD, PDC_BONDS, PDC_INTEGRATES, PDC_MOVES, pd_partition
   implicit none
   private
   
   
   ! Expose type/constructor/methods
   public :: lss, PDC_MOVES, PDC_IS_DEAD, PDC_BONDS, PDC_INTEGRATES, pd_partition
   
   
   !> Memory adaptation parameter
   real(WP), parameter :: coeff_up=1.3_WP      !< Particle array size increase factor
   real(WP), parameter :: coeff_dn=0.7_WP      !< Particle array size decrease factor
   

   !> I/O chunk size to read at a time
   integer, parameter :: part_chunk_size=1000  !< Read 1000 particles at a time before redistributing

   ! MPI message tags for the NOSB <-> fluid particle bridge.
   integer, parameter :: TAG_PARTICLE_COUNT = 1001
   integer, parameter :: TAG_PARTICLE_GID   = 1002
   integer, parameter :: TAG_PARTICLE_Y     = 1003
   integer, parameter :: TAG_PARTICLE_V     = 1004
   integer, parameter :: TAG_PARTICLE_TEST  = 1005

   integer, parameter :: TAG_FORCE_COUNT    = 1011
   integer, parameter :: TAG_FORCE_GID      = 1012
   integer, parameter :: TAG_FORCE_FF       = 1013

   ! We will use this copy for giving and recieving fluid solver information, since we don't want to pass the entire lss object
   type :: pd_copy
      integer :: nown=0                        !< Number of particles in this copy
      real(WP),    allocatable :: y(:,:)       !< position, (3,nown)
      real(WP),    allocatable :: v(:,:)       !< velocity, (3,nown)
      real(WP),    allocatable :: ff(:,:)      !< filled locally by fluid, (3,nown)
   end type pd_copy

   !> Lagrangian solid solver object definition
   !> Extends the existing pdsolver_class, incorporating the coupling functions
   type, extends(pdsolver) :: lss
   
      ! This config is used for parallelization and for calculating bond/collision forces
      class(config), pointer :: cfg

      type(ddadi) :: implicit                             !< Implicit solver for filtering
      type(pd_copy) :: fluid_copy                         !< Copy of the particle data to send to the fluid solver
      ! Solid volume fraction and momentum
      real(WP), dimension(:,:,:), allocatable :: VF       !< Volume fraction, cell-centered
      real(WP), dimension(:,:,:), allocatable :: VFU      !< Solid velocity, U-face
      real(WP), dimension(:,:,:), allocatable :: VFV      !< Solid velocity, V-face
      real(WP), dimension(:,:,:), allocatable :: VFW      !< Solid velocity, W-face
      
       ! CFL numbers
      real(WP) :: CFLp_x,CFLp_y,CFLp_z,CFLp_a
      
      real(WP) :: VFmax                              !< Volume fraction info
      real(WP), dimension(3) :: ibmForce             !< Total force due to IBM

      ! Filtering operation
      real(WP) :: filter_width                       !< Characteristic filter width
      real(WP), dimension(:,:,:,:), allocatable :: div_x,div_y,div_z    !< Divergence operator
      real(WP), dimension(:,:,:,:), allocatable :: grd_x,grd_y,grd_z    !< Gradient operator

      ! Compatibility with the old non-amr version
      integer, dimension(:,:), allocatable :: icell   !< Index of cell containing the particle 
                                                      !< (this might be unnecessary or already exist somewhere, 
                                                      !<  but not in the pdsolver alone I think)

      ! Moving or not (allow flow to setup)
      real(WP) :: unfreeze_time

      ! Communcation related fields for handling the fluid solver
      integer, allocatable :: fluid_rank(:)           !< Fluid rank associated with the particle
      
      ! Communication testing
      integer, allocatable :: which_rank(:)          !< Which rank owns me (for debugging, comment for runs)
   contains
      procedure :: advance                           !< Step forward the particle ODEs
      procedure :: update_partmesh                   !< Update a partmesh object using current particles
      procedure :: update_VF                         !< Compute volume fraction
      procedure :: filter                            !< Apply volume filtering to field
      ! procedure :: get_cfl
      procedure :: update_fluid_location             !< Update the fluid rank and index for each particle
      procedure :: fluid_sync
      procedure :: compute_fluid_forces
   end type lss

   
   
   
   !> Declare lss constructor
   interface lss
      procedure constructor
   end interface lss
   
contains
   
   
   !> Default constructor for Lagrangian solid solver
   function constructor(cfg,name) result(self)
      implicit none
      type(lss) :: self
      class(config), target, intent(in) :: cfg
      character(len=*), optional :: name
      integer :: i,j,k
      
      ! Set the name for the solver
      if (present(name)) self%name=trim(adjustl(name))
      
      ! Point to pgrid object
      self%cfg=>cfg

      self%Ldom = [self%cfg%xL, self%cfg%yL, self%cfg%zL]

      self%per = [self%cfg%xper, &
                  self%cfg%yper, &
                  self%cfg%zper]

      self%collapsed = [self%cfg%nx.eq.1, &
                        self%cfg%ny.eq.1, &
                        self%cfg%nz.eq.1]

      self%dom_lo = [self%cfg%x(self%cfg%imin), &
                     self%cfg%y(self%cfg%jmin), &
                     self%cfg%z(self%cfg%kmin)]

      self%dom_hi = [self%cfg%x(self%cfg%imax+1), &
                     self%cfg%y(self%cfg%jmax+1), &
                     self%cfg%z(self%cfg%kmax+1)]
      
      ! ! Initialize MPI derived datatype for a particle
      ! call prepare_mpi_part() ! IVM, do we need this still? I think that pdsolver handles communication...

      
      ! Allocate VF array on cfg mesh
      allocate(self%VF(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%VF=0.0_WP
      allocate(self%VFU(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%VFU=0.0_WP
      allocate(self%VFV(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%VFV=0.0_WP
      allocate(self%VFW(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%VFW=0.0_WP

      ! Allocate finite volume divergence operators
      allocate(self%div_x(0:+1,self%cfg%imin_:self%cfg%imax_,self%cfg%jmin_:self%cfg%jmax_,self%cfg%kmin_:self%cfg%kmax_)) !< Cell-centered
      allocate(self%div_y(0:+1,self%cfg%imin_:self%cfg%imax_,self%cfg%jmin_:self%cfg%jmax_,self%cfg%kmin_:self%cfg%kmax_)) !< Cell-centered
      allocate(self%div_z(0:+1,self%cfg%imin_:self%cfg%imax_,self%cfg%jmin_:self%cfg%jmax_,self%cfg%kmin_:self%cfg%kmax_)) !< Cell-centered
      ! Create divergence operator to cell center [xm,ym,zm]
      do k=self%cfg%kmin_,self%cfg%kmax_
         do j=self%cfg%jmin_,self%cfg%jmax_
            do i=self%cfg%imin_,self%cfg%imax_
               self%div_x(:,i,j,k)=self%cfg%dxi(i)*[-1.0_WP,+1.0_WP] !< Divergence from [x ,ym,zm]
               self%div_y(:,i,j,k)=self%cfg%dyi(j)*[-1.0_WP,+1.0_WP] !< Divergence from [xm,y ,zm]
               self%div_z(:,i,j,k)=self%cfg%dzi(k)*[-1.0_WP,+1.0_WP] !< Divergence from [xm,ym,z ]
            end do
         end do
      end do

      ! Allocate finite difference velocity gradient operators
      allocate(self%grd_x(-1:0,self%cfg%imin_:self%cfg%imax_+1,self%cfg%jmin_:self%cfg%jmax_+1,self%cfg%kmin_:self%cfg%kmax_+1)) !< X-face-centered
      allocate(self%grd_y(-1:0,self%cfg%imin_:self%cfg%imax_+1,self%cfg%jmin_:self%cfg%jmax_+1,self%cfg%kmin_:self%cfg%kmax_+1)) !< Y-face-centered
      allocate(self%grd_z(-1:0,self%cfg%imin_:self%cfg%imax_+1,self%cfg%jmin_:self%cfg%jmax_+1,self%cfg%kmin_:self%cfg%kmax_+1)) !< Z-face-centered
      ! Create gradient coefficients to cell faces
      do k=self%cfg%kmin_,self%cfg%kmax_+1
         do j=self%cfg%jmin_,self%cfg%jmax_+1
            do i=self%cfg%imin_,self%cfg%imax_+1
               self%grd_x(:,i,j,k)=self%cfg%dxmi(i)*[-1.0_WP,+1.0_WP] !< Gradient in x from [xm,ym,zm] to [x,ym,zm]
               self%grd_y(:,i,j,k)=self%cfg%dymi(j)*[-1.0_WP,+1.0_WP] !< Gradient in y from [xm,ym,zm] to [xm,y,zm]
               self%grd_z(:,i,j,k)=self%cfg%dzmi(k)*[-1.0_WP,+1.0_WP] !< Gradient in z from [xm,ym,zm] to [xm,ym,z]
            end do
         end do
      end do

      ! Loop over the domain and zero divergence in walls
      do k=self%cfg%kmin_,self%cfg%kmax_
         do j=self%cfg%jmin_,self%cfg%jmax_
            do i=self%cfg%imin_,self%cfg%imax_
               if (self%cfg%VF(i,j,k).eq.0.0_WP) then
                  self%div_x(:,i,j,k)=0.0_WP
                  self%div_y(:,i,j,k)=0.0_WP
                  self%div_z(:,i,j,k)=0.0_WP
               end if
            end do
         end do
      end do
      
      ! Zero out gradient to wall faces
      do k=self%cfg%kmin_,self%cfg%kmax_+1
         do j=self%cfg%jmin_,self%cfg%jmax_+1
            do i=self%cfg%imin_,self%cfg%imax_+1
               if (self%cfg%VF(i,j,k).eq.0.0_WP.or.self%cfg%VF(i-1,j,k).eq.0.0_WP) self%grd_x(:,i,j,k)=0.0_WP
               if (self%cfg%VF(i,j,k).eq.0.0_WP.or.self%cfg%VF(i,j-1,k).eq.0.0_WP) self%grd_y(:,i,j,k)=0.0_WP
               if (self%cfg%VF(i,j,k).eq.0.0_WP.or.self%cfg%VF(i,j,k-1).eq.0.0_WP) self%grd_z(:,i,j,k)=0.0_WP
            end do
         end do
      end do

      ! Adjust metrics to account for lower dimensionality
      if (self%cfg%nx.eq.1) then
         self%div_x=0.0_WP
         self%grd_x=0.0_WP
      end if
      if (self%cfg%ny.eq.1) then
         self%div_y=0.0_WP
         self%grd_y=0.0_WP
      end if
      if (self%cfg%nz.eq.1) then
         self%div_z=0.0_WP
         self%grd_z=0.0_WP
      end if

      ! Create implicit solver object for filtering
      self%implicit=ddadi(cfg=self%cfg,name='Filter',nst=7)
      self%implicit%stc(1,:)=[ 0, 0, 0]
      self%implicit%stc(2,:)=[+1, 0, 0]
      self%implicit%stc(3,:)=[-1, 0, 0]
      self%implicit%stc(4,:)=[ 0,+1, 0]
      self%implicit%stc(5,:)=[ 0,-1, 0]
      self%implicit%stc(6,:)=[ 0, 0,+1]
      self%implicit%stc(7,:)=[ 0, 0,-1]
      call self%implicit%init()

      ! Set default filter width
      self%filter_width=1.0_WP*self%cfg%min_meshsize

      ! Log/screen output
      logging: block
         use, intrinsic :: iso_fortran_env, only: output_unit
         use param,    only: verbose
         use messager, only: log
         use string,   only: str_long
         character(len=str_long) :: message
         if (self%cfg%amRoot) then
            write(message,'("LSS object [",a,"] on partitioned grid [",a,"]")') trim(self%name),trim(self%cfg%name)
            if (verbose.gt.1) write(output_unit,'(a)') trim(message)
            if (verbose.gt.0) call log(message)
         end if
      end block logging
      
   end function constructor

   

   !> Advance the particle equations by a specified time step dt
   subroutine advance(this,dt,unfreeze,div_stress_x,div_stress_y,div_stress_z)
      implicit none
      class(lss), intent(inout) :: this
      real(WP), intent(inout) :: dt  !< Timestep size over which to advance
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: div_stress_x  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: div_stress_y  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: div_stress_z  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: n,i
      logical, intent(in) :: unfreeze
      ! updates particle positions, shares particles to fluid-cell rank-owners, compute fluid force, and return ff (fluid force) 
      call this%fluid_sync(d_stress_x=div_stress_x,d_stress_y=div_stress_y,d_stress_z=div_stress_z) 
      ! now the particles are communicated down to get their fluid forces
      call this%pd_advance(dt) ! use fluid forces and compute bond forces, and update position due to verlet scheme
      call this%update_VF()
      if (unfreeze) then
         do i = 1,this%nown
            this%damping_rate = 0.0005_WP
            if (this%flag(i).eq.(PDC_MOVES+PDC_BONDS)) this%v(:,i)= 0.0_WP
         end do
      end if


      
      ! call this%update_VF() ! now we update the volume fraction

      
      ! Log/screen output (do we need to do this still?)
      ! logging: block
      !    use, intrinsic :: iso_fortran_env, only: output_unit
      !    use param,    only: verbose
      !    use messager, only: log
      !    use string,   only: str_long
      !    character(len=str_long) :: message
      !    if (this%cfg%amRoot) then
      !       write(message,'("Particle solver [",a,"] on partitioned grid [",a,"]: ",i0," particles were advanced")') trim(this%name),trim(this%cfg%name),this%np
      !       if (verbose.gt.1) write(output_unit,'(a)') trim(message)
      !       if (verbose.gt.0) call log(message)
      !    end if
      ! end block logging
      
   end subroutine advance


   !> Update particle volume fraction using our current particles
   subroutine update_VF(this)
      implicit none
      class(lss), intent(inout) :: this
      integer :: i
      integer, dimension(3) :: idx
      ! Reset volume fraction and momentum
      this%VF=0.0_WP; this%VFU=0.0_WP; this%VFV=0.0_WP; this%VFW=0.0_WP
      ! Transfer particle volume
      idx = 0
      do i=1,this%fluid_copy%nown ! we only do the particles that we actually physically have
         ! Transfer volume to mesh
         idx = this%cfg%get_ijk_global(this%fluid_copy%y(:,i),idx)
         call this%cfg%set_scalar(Sp=this%dV,                        pos=this%fluid_copy%y(:,i),i0=idx(1),j0=idx(2),k0=idx(3),S=this%VF ,bc='n')
         call this%cfg%set_scalar(Sp=this%dV*this%fluid_copy%v(1,i), pos=this%fluid_copy%y(:,i),i0=idx(1),j0=idx(2),k0=idx(3),S=this%VFU,bc='n')
         call this%cfg%set_scalar(Sp=this%dV*this%fluid_copy%v(2,i), pos=this%fluid_copy%y(:,i),i0=idx(1),j0=idx(2),k0=idx(3),S=this%VFV,bc='n')
         call this%cfg%set_scalar(Sp=this%dV*this%fluid_copy%v(3,i), pos=this%fluid_copy%y(:,i),i0=idx(1),j0=idx(2),k0=idx(3),S=this%VFW,bc='n')
      end do
      this%VF =this%VF /this%cfg%vol
      this%VFU=this%VFU/this%cfg%vol
      this%VFV=this%VFV/this%cfg%vol
      this%VFW=this%VFW/this%cfg%vol
      ! Sum at boundaries
      call this%cfg%syncsum(this%VF )
      call this%cfg%syncsum(this%VFU)
      call this%cfg%syncsum(this%VFV)
      call this%cfg%syncsum(this%VFW)
      ! Apply volume filter
      call this%filter(this%VF )
      call this%filter(this%VFU)
      call this%filter(this%VFV)
      call this%filter(this%VFW)
      ! Clip
      where (this%VF.lt.0.0_WP) this%VF=0.0_WP
      this%VF=min(this%VF,1.0_WP-epsilon(1.0_WP))

    end subroutine update_VF

   !  subroutine get_cfl(this,dt,cfl)
   !    use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX
   !    use parallel, only: MPI_REAL_WP
   !    implicit none
   !    class(lss), intent(inout) :: this
   !    real(WP), intent(in)  :: dt
   !    real(WP), intent(out) :: cfl
   !    integer :: i,ierr
   !    real(WP) :: my_CFLp_x,my_CFLp_y,my_CFLp_z,kk,mu,a
      
   !    ! Set the CFLs to zero
   !    my_CFLp_x=0.0_WP; my_CFLp_y=0.0_WP; my_CFLp_z=0.0_WP
   !    do i=1,this%nown
   !       my_CFLp_x=max(my_CFLp_x,abs(this%v(1,i))*this%cfg%dxi(this%icell(1,i)))
   !       my_CFLp_y=max(my_CFLp_y,abs(this%v(2,i))*this%cfg%dyi(this%icell(2,i)))
   !       my_CFLp_z=max(my_CFLp_z,abs(this%v(3,i))*this%cfg%dzi(this%icell(3,i)))
   !    end do
   !    my_CFLp_x=my_CFLp_x*dt; my_CFLp_y=my_CFLp_y*dt; my_CFLp_z=my_CFLp_z*dt
      
   !    ! Get the parallel max
   !    call MPI_ALLREDUCE(my_CFLp_x,this%CFLp_x,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
   !    call MPI_ALLREDUCE(my_CFLp_y,this%CFLp_y,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
   !    call MPI_ALLREDUCE(my_CFLp_z,this%CFLp_z,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)

   !    ! CFL based on elastic wave speed in material
   !    kk=this%elastic_modulus/(3.0_WP-6.0_WP*this%poisson_ratio)
   !    mu=this%elastic_modulus/(2.0_WP+2.0_WP*this%poisson_ratio)      
   !    a=sqrt((kk+4.0_WP*mu/3.0_WP)/this%rho)
   !    this%CFLp_a=dt*a/this%delta
      
   !    ! Return the maximum CFL
   !    cfl=max(this%CFLp_x,this%CFLp_y,this%CFLp_z,this%CFLp_a)
      
   ! end subroutine get_cfl
    

    !> Laplacian filtering operation
    subroutine filter(this,A)
      implicit none
      class(lss), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: A     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP) :: filter_coeff
      integer :: i,j,k,n,nstep
      real(WP), dimension(:,:,:), allocatable :: FX,FY,FZ

      ! Recompute filter coefficient
      filter_coeff=max(this%filter_width**2-this%cfg%min_meshsize**2,0.0_WP)/(16.0_WP*log(2.0_WP))
      if (filter_coeff.le.0.0_WP) return

      ! Allocate flux arrays
      allocate(FX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))

      if (.not.this%implicit%setup_done) then
         ! Prepare diffusive operator (only need to do this once)
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  this%implicit%opr(1,i,j,k)=1.0_WP-(this%div_x(+1,i,j,k)*filter_coeff*this%grd_x(-1,i+1,j,k)+&
                  &                                  this%div_x( 0,i,j,k)*filter_coeff*this%grd_x( 0,i  ,j,k)+&
                  &                                  this%div_y(+1,i,j,k)*filter_coeff*this%grd_y(-1,i,j+1,k)+&
                  &                                  this%div_y( 0,i,j,k)*filter_coeff*this%grd_y( 0,i,j  ,k)+&
                  &                                  this%div_z(+1,i,j,k)*filter_coeff*this%grd_z(-1,i,j,k+1)+&
                  &                                  this%div_z( 0,i,j,k)*filter_coeff*this%grd_z( 0,i,j,k  ))
                  this%implicit%opr(2,i,j,k)=      -(this%div_x(+1,i,j,k)*filter_coeff*this%grd_x( 0,i+1,j,k))
                  this%implicit%opr(3,i,j,k)=      -(this%div_x( 0,i,j,k)*filter_coeff*this%grd_x(-1,i  ,j,k))
                  this%implicit%opr(4,i,j,k)=      -(this%div_y(+1,i,j,k)*filter_coeff*this%grd_y( 0,i,j+1,k))
                  this%implicit%opr(5,i,j,k)=      -(this%div_y( 0,i,j,k)*filter_coeff*this%grd_y(-1,i,j  ,k))
                  this%implicit%opr(6,i,j,k)=      -(this%div_z(+1,i,j,k)*filter_coeff*this%grd_z( 0,i,j,k+1))
                  this%implicit%opr(7,i,j,k)=      -(this%div_z( 0,i,j,k)*filter_coeff*this%grd_z(-1,i,j,k  ))
               end do
            end do
         end do
      end if
      ! Explicit step
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               FX(i,j,k)=filter_coeff*sum(this%grd_x(:,i,j,k)*A(i-1:i,j,k))
               FY(i,j,k)=filter_coeff*sum(this%grd_y(:,i,j,k)*A(i,j-1:j,k))
               FZ(i,j,k)=filter_coeff*sum(this%grd_z(:,i,j,k)*A(i,j,k-1:k))
            end do
         end do
      end do
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%implicit%rhs(i,j,k)=sum(this%div_x(:,i,j,k)*FX(i:i+1,j,k))+sum(this%div_y(:,i,j,k)*FY(i,j:j+1,k))+sum(this%div_z(:,i,j,k)*FZ(i,j,k:k+1))
            end do
         end do
      end do
      ! Implicit step
      call this%implicit%setup()
      this%implicit%sol=0.0_WP
      call this%implicit%solve()
      A=A+this%implicit%sol
      call this%cfg%sync(A)

      ! Deallocate flux arrays
      deallocate(FX,FY,FZ)

    end subroutine filter
   
   !> Update particle mesh using our current particles
   subroutine update_partmesh(this,pmesh)
      use partmesh_class, only: partmesh
      implicit none
      class(lss), intent(inout) :: this
      class(partmesh), intent(inout) :: pmesh
      integer :: i
      ! Reset particle mesh storage
      call pmesh%reset()
      ! Nothing else to do if no particle is present
      if (this%nown.eq.0) return
      ! Copy particle info
      call pmesh%set_size(this%nown)
      do i=1,this%nown !< IVM, this might not be good, I think we will get duplicates this way
         pmesh%pos(:,i)=this%y(:,i)
      end do
   end subroutine update_partmesh
   
   
   ! !> Creation of the MPI datatype for particle ! IVM, Maybe we dont need this, since comm is handled by pdsolver?
   ! subroutine prepare_mpi_part()
   !    use mpi_f08
   !    use messager, only: die
   !    implicit none
   !    integer(MPI_ADDRESS_KIND), dimension(part_nblock) :: disp
   !    integer(MPI_ADDRESS_KIND) :: lb,extent
   !    type(MPI_Datatype) :: MPI_PART_TMP
   !    integer :: i,mysize,ierr
   !    ! Prepare the displacement array
   !    disp(1)=0
   !    do i=2,part_nblock
   !       call MPI_Type_size(part_tblock(i-1),mysize,ierr)
   !       disp(i)=disp(i-1)+int(mysize,MPI_ADDRESS_KIND)*int(part_lblock(i-1),MPI_ADDRESS_KIND)
   !    end do
   !    ! Create and commit the new type
   !    call MPI_Type_create_struct(part_nblock,part_lblock,disp,part_tblock,MPI_PART_TMP,ierr)
   !    call MPI_Type_get_extent(MPI_PART_TMP,lb,extent,ierr)
   !    call MPI_Type_create_resized(MPI_PART_TMP,lb,extent,MPI_PART,ierr)
   !    call MPI_Type_commit(MPI_PART,ierr)
   !    ! If a problem was encountered, say it
   !    if (ierr.ne.0) call die('[lss prepare_mpi_part] MPI Particle type creation failed')
   !    ! Get the size of this type
   !    call MPI_type_size(MPI_PART,MPI_PART_SIZE,ierr)
   ! end subroutine prepare_mpi_part

   subroutine update_fluid_location(this)
      implicit none
      class(lss), intent(inout) :: this
      integer :: i
      ! Ensure that we have allocated the fluid_rank array
      if (allocated(this%fluid_rank)) then
         if (size(this%fluid_rank).lt.max(this%nown,1)) then
            deallocate(this%fluid_rank)
         end if
      end if
      if (.not.allocated(this%fluid_rank)) then
         allocate(this%fluid_rank(max(this%nown,1)))
      end if

      if (allocated(this%icell)) then
         if (size(this%icell,dim=2).lt.max(this%nown,1)) then
            deallocate(this%icell)
         end if
      end if
      if (.not.allocated(this%icell)) then
         allocate(this%icell(3,max(this%nown,1)))
         this%icell=0
      end if

      ! Update the fluid_rank and icell for each particle
      do i=1,this%nown
         if (this%flag(i).eq.PDC_IS_DEAD) cycle
         ! Get the cell index for the particle
         this%icell(:,i) = this%cfg%get_ijk_global(this%y(:,i), this%icell(:,i))
         ! Get the fluid rank for the particle
         this%fluid_rank(i) = this%cfg%get_rank(this%icell(:,i))
      end do

   end subroutine update_fluid_location

   ! subroutine share_particles(this)
   !    use parallel, only: MPI_REAL_WP
   !    use mpi_f08
   !    implicit none

   !    class(lss), intent(inout) :: this
   !    integer :: i, ierr, nranks, r
   !    integer, allocatable :: send_count(:)     ! List (nproc) of how many particles we should send
   !    integer, allocatable :: recv_count(:)     ! List (nproc) of how many particles we should recieve
   !    integer, allocatable :: send_disp(:)      ! Displacement for sending particles
   !    integer, allocatable :: recv_disp(:)      ! Displacement for recieving particles
   !    integer :: nsend, nrecv                   ! Total number of particles to send and recieve
   !    integer :: nreq                           ! Number of required messages
   !    integer :: first                          ! first index of the recieve buffer for a given rank, and which message we are on
   !    integer :: q                              ! which message we are on
   !    integer :: slot                           ! index for where the rank information starta
   !    real(WP), allocatable :: sy(:,:), sv(:,:) ! Send buffers for position, velocity
   !    real(WP), allocatable :: ry(:,:), rv(:,:) ! Recieve buffers for position, velocity
   !    integer, allocatable :: stest(:), rtest(:)  ! Send and recieve buffers for testing
   !    integer, allocatable :: next(:)                  ! Next index for sending particles to a given rank
   !    type(MPI_Request), allocatable :: req(:)  ! MPI requests for non-blocking communication
   !    type(MPI_Status),  allocatable :: stat(:) ! Status for the MPI requests
   !    ! First we figure out which ranks we need to communicate with
   !    call this%update_fluid_location() ! We update the fluid location and rank information for each particle nown

   !    nranks=this%cfg%nproc
   !    allocate(send_count(0:nranks-1), recv_count(0:nranks-1))
   !    send_count=0
   !    recv_count=0
   !    do i = 1,this%nown
   !       if (this%flag(i).eq.PDC_IS_DEAD) cycle
   !       send_count(this%fluid_rank(i)) = send_count(this%fluid_rank(i)) + 1 ! count up how many particles need to get passed along 
   !    end do

   !    ! Now we populate the recv_count array by doing an all-to-all communication (only 1)

   !    call MPI_ALLTOALL(send_count,1,MPI_INTEGER, &
   !                    & recv_count,1,MPI_INTEGER, &
   !                    & this%cfg%comm,ierr)

   !    ! now everyone knows who they are recieving from and how many particles to get
   !    ! now we size the recieve buffers to fit the amount of information and number of particles

   !    ! Total number of particle to send and recieve on this rank
   !    nsend = sum(send_count)
   !    nrecv = sum(recv_count)

   !    ! track the displacement needed for each rank for sending and recieving
   !    allocate(send_disp(0:nranks-1), recv_disp(0:nranks-1))
   !    send_disp(0)=0; recv_disp(0)=0
   !    do r=1,nranks-1
   !       ! displacement for the contigous send and recieve buffers
   !       ! we start at 0 on the 0th index, then add the number of particles to send
   !       ! we then start the nexxt rank at the previous displacement plus the number of particles to send from that rank
   !       ! e.g. if rank 0 is sending 3 particles, then rank 1 will start at index 3 in the send buffer, 
   !       ! and if rank 1 is sending 5 particles, then rank 2 will start at index 8 in the send buffer
   !       send_disp(r) = send_disp(r-1) + send_count(r-1)
   !       ! Same process for the recieve buffer
   !       recv_disp(r) = recv_disp(r-1) + recv_count(r-1)
   !    end do

   !    ! since we have different data types we are communicating, we need to allocate buffers for each tyep'
   !    allocate(stest(max(nsend,1)),rtest(max(nrecv,1)))    ! Global id buffers
   !    allocate(sy(3,max(nsend,1)),ry(3,max(nrecv,1)))    ! position buffers
   !    allocate(sv(3,max(nsend,1)),rv(3,max(nrecv,1)))    ! velocity buffers

   !    nreq=3*(count(recv_count.gt.0) + count(send_count.gt.0)) ! number of required messages

   !    if(nreq.gt.0) then
   !       allocate(req(nreq),stat(nreq)) ! allocating status and request arrays for the number of messages
   !       q=0 ! which message we are on

   !       do r=0,nranks-1
   !          if (recv_count(r).eq.0) cycle

   !          first=recv_disp(r)+1 ! first index of the recieve buffer for a given rank

   !          ! Setup recieve buffer for position
   !          q=q+1
   !          call MPI_IRECV(ry(1,first),3*recv_count(r),MPI_REAL_WP, & ! multiply by three since vector
   !          &              r,TAG_PARTICLE_Y,this%cfg%comm,req(q),ierr)

   !          ! Setup recieve buffer for velocity
   !          q=q+1
   !          call MPI_IRECV(rv(1,first),3*recv_count(r),MPI_REAL_WP, & ! mutliply by three since vector
   !          &              r,TAG_PARTICLE_V,this%cfg%comm,req(q),ierr)

   !          ! Setup test buffer for comms testing
   !          q=q+1
   !          call MPI_IRECV(rtest(first),recv_count(r),MPI_INTEGER, & ! mutliply by three since vector
   !          &              r,TAG_PARTICLE_TEST,this%cfg%comm,req(q),ierr)
   !       end do
   !    end if 

   !    ! Now we pack the message
   !   allocate(next(0:nranks-1))
   !   next=send_disp
   !    do i = 1,this%nown
   !       if (this%flag(i).eq.PDC_IS_DEAD) cycle
   !       r=this%fluid_rank(i) !which rank are we sending to

   !       next(r)=next(r)+1 ! we +1 this to put the particles one after 
   !                         ! the next for this rank (if we didn't we would 
   !                         ! overwrite the previous particle for this rank in the send buffer)
   !       slot = next(r) ! we 1 index for the this% but 0 index for the send_disp
   !       sy(:,slot)=this%y(:,i)
   !       sv(:,slot)=this%v(:,i)
   !       stest(slot)=this%which_rank(i)
   !    end do
   !    deallocate(next)

   !    if(nreq.gt.0) then
   !       ! Now we send the messages
   !       do r=0,nranks-1
   !          if (send_count(r).eq.0) cycle

   !          first=send_disp(r)+1 ! first index of the send buffer for a given rank

   !          ! Setup send buffer for position
   !          q=q+1
   !          call MPI_ISEND(sy(1,first),3*send_count(r),MPI_REAL_WP, &
   !          &              r,TAG_PARTICLE_Y,this%cfg%comm,req(q),ierr)

   !          ! Setup send buffer for velocity
   !          q=q+1
   !          call MPI_ISEND(sv(1,first),3*send_count(r),MPI_REAL_WP, &
   !          &              r,TAG_PARTICLE_V,this%cfg%comm,req(q),ierr)

   !          q=q+1
   !          call MPI_ISEND(stest(first),send_count(r),MPI_INTEGER, & 
   !          &              r,TAG_PARTICLE_TEST,this%cfg%comm,req(q),ierr)
   !       end do

         
   !    end if

   !    if (nreq.gt.0) then
   !       if (q.ne.nreq) then
   !          error stop '[share_particles] MPI request count mismatch'
   !       end if

   !       call MPI_WAITALL(nreq,req,stat,ierr) ! make sure everyone has done their sending and recieving

   !       deallocate(req,stat)
   !    end if

   !    do r=0,nranks-1
   !       if (recv_count(r).eq.0) cycle

   !       first=recv_disp(r)+1

   !       ! The block came from rank r, so every value should equal r.
   !       if (any(rtest(first:first+recv_count(r)-1).ne.r)) then
   !          write(*,*) 'Rank ',this%cfg%rank, &
   !          &          ' received incorrect test data from rank ',r
   !          error stop '[share_particles] which_rank test failed'
   !       end if
   !    end do

   !    write(*,*) 'Rank ',this%cfg%rank, &
   !    &          ': particle communication test passed; received ',nrecv

      
    
   ! end subroutine share_particles

   ! Sends out particles accross ranks to those who own them for computing fluid forces on the particles
   ! Also needed to update volume fractions
   subroutine fluid_sync(this,d_stress_x,d_stress_y,d_stress_z)
      use parallel, only: MPI_REAL_WP
      use mpi_f08
      implicit none

      class(lss), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: d_stress_x  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: d_stress_y  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: d_stress_z  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      
      integer :: i, ierr, nranks, r, nsend, nrecv, rank, column
      integer, allocatable :: send_count(:), recv_count(:)
      integer, allocatable :: send_disp(:), recv_disp(:)
      integer, allocatable :: next(:)
      integer, allocatable :: send_lid(:) ! local id of the send (needed for the return trip)
      real(WP), allocatable :: send_yv(:,:), recv_yv(:,:), recv_ff(:,:)
      ! Maybe we do this seperatately instead of tying it in?
      call this%update_fluid_location() ! We update the fluid location and rank information for each particle nown

      nranks=this%cfg%nproc
      allocate(send_count(0:nranks-1), recv_count(0:nranks-1))
      send_count=0
      recv_count=0
      do i = 1,this%nown
         if (this%flag(i).eq.PDC_IS_DEAD) cycle
         send_count(this%fluid_rank(i)) = send_count(this%fluid_rank(i)) + 1 ! count up how many particles need to get passed along 
      end do
      ! Now we populate the recv_count array by doing an all-to-all communication (only 1)
      call MPI_ALLTOALL(send_count,1,MPI_INTEGER, &
                      & recv_count,1,MPI_INTEGER, &
                      & this%cfg%comm,ierr)
      nsend = sum(send_count) ! total number of particles that are being send 
      nrecv = sum(recv_count) ! total number of particles we expect to recieve
      ! Set up the displacement counts based on number of particles
      allocate(send_disp(0:nranks-1), recv_disp(0:nranks-1))
      send_disp(0) = 0
      recv_disp(0) = 0
      do rank=1,nranks-1
         send_disp(rank) = send_disp(rank-1) + send_count(rank-1)
         recv_disp(rank) = recv_disp(rank-1) + recv_count(rank-1)
      end do
      ! Pack message
      allocate(send_yv(6,max(nsend,1)))
      allocate(recv_yv(6,max(nrecv,1)))
      send_yv=0.0_WP
      recv_yv=0.0_WP
      allocate(next(0:nranks-1))
      next=send_disp ! cursor
      ! setup local map
      allocate(send_lid(max(nsend,1)))
      send_lid=0
      ! pack up the send buffer
      do i = 1,this%nown
         if (this%flag(i).eq.PDC_IS_DEAD) cycle
         rank = this%fluid_rank(i)
         next(rank) = next(rank) + 1 ! plus up for each particle that we are including starting at the rank displacement +1, (this will 
                                     ! not overwrite because the next(:) is based on the displacments and ensures the column mappings keeps
                                     ! same-rank particles in blocks of adjacent columns)
         column = next(rank) ! this is the index for the start of the information for that particle
         ! The way that fortran stores information is column major so we associate each particle with a particular column
         ! that way when we send it, it sends all in order the information for a particle, and we can unwrap
         ! the send information in the same way if we construct the recieve buffer the same
         send_yv(1:3,column) = this%y(:,i)
         send_yv(4:6,column) = this%v(:,i)
         ! save the local id for the way back
         send_lid(column)=i
      end do
      !    MPI_Alltoallv(
      !    sendbuf,        Starting address of the send buffer in memory
      !    sendcounts,     Counts for the number of elements to send to each rank
      !    sdispls,        Displacements for the starting address of each rank's data in the send buffer, relative to sendbuf
      !    sendtype,       Data type of the send buffer elements
      !    recvbuf,        Starting address of the receive buffer in memory
      !    recvcounts[],   Counts for the number of elements to receive from each rank
      !    rdispls[],      Displacements for the starting address of each rank's data in the receive buffer, relative to recvbuf
      !    recvtype,       Data type of the receive buffer elements
      !    comm            Communicator handle
      !    )
      ! each entry of send_counts is just 6 times the number of particles we are sending to that rank

      call MPI_ALLTOALLV(send_yv,send_count*6,send_disp*6,MPI_REAL_WP, &
                        recv_yv,recv_count*6,recv_disp*6,MPI_REAL_WP, &
                        this%cfg%comm,ierr)
      
      ! now recv_yv has all the particles that this processor needed, which we need to reconstruct
      ! clear out the existing copies
      if (allocated(this%fluid_copy%y)) then
         deallocate(this%fluid_copy%y)
      end if

      if (allocated(this%fluid_copy%v)) then
         deallocate(this%fluid_copy%v)
      end if
      this%fluid_copy%nown=nrecv
      allocate(this%fluid_copy%y(3,max(nrecv,1)))
      allocate(this%fluid_copy%v(3,max(nrecv,1)))
      this%fluid_copy%y=0.0_WP
      this%fluid_copy%v=0.0_WP

      if (nrecv.gt.0) then
         this%fluid_copy%y(:,1:nrecv)=recv_yv(1:3,1:nrecv)
         this%fluid_copy%v(:,1:nrecv)=recv_yv(4:6,1:nrecv)
      end if

      ! I think if we are careful about the order of nown we send and recieve, and keep it identically the same
      ! we can get away without having to send the particle global id

      if (allocated(this%fluid_copy%ff)) then
         deallocate(this%fluid_copy%ff)
      end if
      allocate(this%fluid_copy%ff(3,max(nrecv,1)))
      this%fluid_copy%ff=0.0_WP
      call this%update_VF()
      call this%compute_fluid_forces(stress_x=d_stress_x,stress_y=d_stress_y,stress_z=d_stress_z)

      ! Back the way we came
      allocate(recv_ff(3,max(nsend,1)))
      recv_ff=0.0_WP
      ! swap send and recieve since we are getting things back (should be the same order?)
      call MPI_ALLTOALLV(this%fluid_copy%ff,recv_count*3,recv_disp*3,MPI_REAL_WP, &
                        & recv_ff,send_count*3,send_disp*3,MPI_REAL_WP, &
                        & this%cfg%comm,ierr)
      do column=1,nsend
         i=send_lid(column)
         if (i.lt.1 .or. i.gt.this%nown) then
            error stop '[sync] Invalid returned-force mapping'
         end if
         this%ff(:,i)=recv_ff(:,column)
      end do
      deallocate(next)
      deallocate(send_yv,recv_yv)
      deallocate(send_count,recv_count)
      deallocate(recv_ff,send_lid)

   end subroutine fluid_sync

   ! Compute fluid forces acting on fluid_copy particles on each rank
   subroutine compute_fluid_forces(this,stress_x,stress_y,stress_z)
      implicit none
      class(lss), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: stress_x  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: stress_y  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(inout) :: stress_z  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: n,ierr
      real(WP), dimension(3) :: stress
      integer, dimension(3) :: idx
      ! Assumes that we have already shared copies and particle locations with ranks
      idx = 0
      ! Compute the fluid forces on each particle using the copied particles we have
      do n=1,this%fluid_copy%nown
         ! Advance with Verlet scheme
         
         idx = this%cfg%get_ijk_global(this%fluid_copy%y(:,n),idx) ! is this slow? should we store it?
         this%fluid_copy%ff(:,n)=this%cfg%get_velocity(pos=this%fluid_copy%y(:,n),i0=idx(1),j0=idx(2),k0=idx(3),U=stress_x,V=stress_y,W=stress_z)
         ! we will divide by rho later

      end do
   
   end subroutine compute_fluid_forces 
   
end module lsspd_class
