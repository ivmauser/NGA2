!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP, I8
   use geometry,          only: cfg
   use fft2d_class,       only: fft2d
   use ddadi_class,       only: ddadi
   use incomp_class,      only: incomp
   use lsspd_class,       only: lss, pd_partition, PDC_MOVES,PDC_INTEGRATES,PDC_BONDS,PDC_SURFACE
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use partmesh_class,    only: partmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   use string,            only: str_medium
   use datafile_class,    only: datafile
   implicit none
   private
   
   !> Get a couple linear solvers, an incompressible flow solver and corresponding time tracker
   type(fft2d),       public :: ps
   type(ddadi),       public :: vs
   type(incomp),      public :: fs
   type(lss),         public :: ls
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   
   type(ensight) :: ens_out
   type(event)   :: ens_evt, save_evt
   type(datafile) :: df
   logical :: restarted = .false.
   type(partmesh),    public :: pmesh
   character(len=str_medium) :: restart_dir
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,sfile,tfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:), allocatable :: div_x,div_y,div_z
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:), allocatable :: Uib,Vib,Wib,srcM
   real(WP), dimension(:,:,:,:,:), allocatable :: gradU

   !> Max timestep size for solid solver
   real(WP) :: ls_dt,ls_dt_max
   integer :: solid_substeps = 0

   !> Steady state evolution boolean
   logical :: steady_state = .false.
   
 contains


   !> Function that localizes the left (x-) of the domain
   function left_of_domain(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
     implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin) isIn=.true.
   end function left_of_domain
   
   
   !> Function that localizes the right (x+) of the domain
   function right_of_domain(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
     implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1) isIn=.true.
   end function right_of_domain


   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read,param_exists
      use parallel, only: amRoot
      implicit none

      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
        time=timetracker(amRoot=cfg%amRoot)
        call param_read('Max timestep size',time%dtmax)
        call param_read('Max cfl number',time%cflmax)
        call param_read('Max time',time%tmax)
        time%dt=time%dtmax
        time%itmax=2
      end block initialize_timetracker

      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(div_x(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(div_y(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(div_z(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Uib (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vib (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wib (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(srcM(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(gradU(1:3,1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays

      ! Handle restart/saves here
      restart_and_save: block
         restart_dir=''
         ! Create event for saving restart files
         save_evt=event(time,'Restart output')
         call param_read('Restart output period',save_evt%tper,default=huge(1.0_WP))
         ! Check if we are restarting
         call param_read(tag='Restart from',val=restart_dir,short='r',default='')
         restarted=.false.; if (len_trim(restart_dir).gt.0) restarted=.true.
         if (restarted) then
            ! If we are, read the name of the directory
            call param_read('Restart from',restart_dir,'r')
            ! Read the datafile
            df = datafile(pg=cfg, fdata=trim(restart_dir)//'/fluid')
         else
            ! Prepare a new directory for storing files for restart
            call execute_command_line('mkdir -p restart')
            ! If we are not restarting, we will still need a datafile for saving restart files
            df=datafile(pg=cfg,filename=trim(cfg%name),nval=4,nvar=4)
            df%valname(1)='t'
            df%valname(2)='dt'
            df%valname(3)='step'
            df%valname(4)='ls_dt'
            df%varname(1)='U'
            df%varname(2)='V'
            df%varname(3)='W'
            df%varname(4)='P'
         end if
      end block restart_and_save

      ! Revisit timetracker to adjust time and time step values if this is a restart
      update_timetracker: block
         if (restarted) then
            call df%pullval(name='t' ,val=time%t )
            call df%pullval(name='dt',val=time%dt)
            time%told=time%t-time%dt
         end if
      end block update_timetracker

      ! Initialize Lagrangian solid solver
      initialize_lss: block
         use mathtools, only: Pi
         
         integer(I8), allocatable :: gids(:),rgid(:)
         real(WP), allocatable :: pos(:,:),vel(:,:),voll(:),rpos(:,:),rvel(:,:),rvol(:)
         integer, allocatable :: flags(:),owner(:),rflag(:)
         integer :: i,j,k,n,nn,nr,nx,ny,nz,N_r,N_z,bnx,bny
         real(WP) :: dx,x0,y0,z0,Beam_L,Beam_H,Lx,Ly,Lz,x_c,y_c,a,b
         real(WP) :: R, R_int
         real(WP) :: rho,E,nu,elem,delta

         ls=lss(cfg=cfg,name='solid')
         
         ! Problem dimensions
         call param_read('R',R,default=0.5_WP)        ! Cylinder radius
         call param_read('N_r',N_r,default=10)        ! Number of particles accross cylinder radius
         call param_read('Lz',Lz,default=6.0_WP)     ! Z dimensionality 
         call param_read('Lx',Lx,default=2.5_WP)      ! Z length of the domain
         call param_read('Ly',Ly,default=0.41_WP)     ! Y length of the domain
         call param_read('x_c',x_c,default=0.2_WP)    ! x-Center of the cylinder in space (measured from lower left corner)
         call param_read('y_c',y_c,default=0.2_WP)    ! y-Center of the cylinder in space (measured from lower left corner)
         call param_read('a',a,default=0.35_WP)       ! Length of beam (measured from end of cylinder)
         call param_read('b',b,default=0.02_WP)       ! Width of beam
         ! need to define an elem value
         elem=R/(3.0_WP*real(N_r,WP))
         ls%delta=elem*3.01_WP ! 3 times the elem value should be the grid spacing
         call param_read('Material density',ls%rho)
         call param_read('Elastic modulus', ls%elastic_modulus,default=2.0e11_WP)
         call param_read('Poisson ratio',   ls%poisson_ratio,default=0.3_WP)
         call param_read('Tau',             ls%tau, default=huge(1.0_WP))
         call param_read('Particle timestep size',ls_dt_max,default=huge(1.0_WP))
         call param_read('Steady state', steady_state)

         ls%damping_rate=0.0_WP
         ls_dt=min(ls_dt_max,time%dtmax)
         ! Configure by field assignment (grid-free: no domain, no periodicity)
         ls%dV=elem**3 ! we treat everybody as a cube
         ! Root builds the whole lattice; pd_partition routes it (gids are
         ! simply 1..n -- any unique positive keys work)
         ! Check dimensionality
         if(cfg%nx.eq.1) ls%collapsed(1)=.true.
         if(cfg%ny.eq.1) ls%collapsed(2)=.true.
         if(cfg%nz.eq.1) ls%collapsed(3)=.true.
         if (restarted) then
            call ls%read_state(dirname=trim(restart_dir))
            ! do i=1,ls%nown
            !    if(ls%flag(i).eq.PDC_BONDS+PDC_MOVES) ls%flag(i)=PDC_BONDS+PDC_MOVES+PDC_INTEGRATES
            ! end do         
         else
            nx = 2*N_r*3
            ny = 2*N_r*3
            ! Z is special for a cylinder, needs to match width of the domain 
            if (Lz.eq.0.0_WP) then
               nz = 1
            else
               nz = N_r*N_z*3
            end if 
            nn=0

            ! Now we set up the beam
            Beam_L = a + R
            Beam_H = b
            bnx = int(Beam_L/elem)
            bny = int(Beam_H/elem)
            R_int = R*(real(N_r-4,WP)/real(N_r,WP))
            if(amRoot) then
               ! Cylinder
               do k=1,nz; do j=1,ny; do i=1,nx
                  
                  x0 = (real(i,WP) - 0.5_WP*real(nx+1,WP))*elem + x_c
                  y0 = (real(j,WP) - 0.5_WP*real(ny+1,WP))*elem + y_c
                  z0 = (real(k,WP) - 0.5_WP*real(nz+1,WP))*elem
                  if (((x0-x_c)*(x0-x_c) + (y0-y_c)*(y0-y_c)).ge.R*R) cycle;
                  if (((x0-x_c)*(x0-x_c) + (y0-y_c)*(y0-y_c)).le.R_int*R_int) cycle;
                  if (((x0-x_c).gt.0.0_WP).and.(abs(y0-y_c).le.Beam_H/2.0_WP)) cycle;
                  nn=nn+1
               end do; end do; end do

               ! Beam
               do k=1,nz; do j=1,bny; do i=1,bnx
                  x0 = (real(i,WP) - 0.5_WP*real(bnx+1,WP))*elem + Beam_L/2.0_WP + x_c
                  y0 = (real(j,WP) - 0.5_WP*real(bny+1,WP))*elem + y_c
                  z0 = (real(k,WP) - 0.5_WP*real(nz+1,WP))*elem
                  ! if (((x0)*(x0) + y0*y0).le.R*R) cycle;
                  nn=nn+1
               end do; end do; end do

               print*, nn
            end if
         
            allocate(gids(max(nn,1)),pos(3,max(nn,1)),vel(3,max(nn,1)),flags(max(nn,1)),voll(max(nn,1)),owner(max(nn,1)))
            ! allocate(ls%icell(3,max(nn,1)))
            n=0
            do k=1,nz; do j=1,ny; do i=1,nx
               if (.not.amRoot) exit
               x0 = (real(i,WP) - 0.5_WP*real(nx+1,WP))*elem + x_c
               y0 = (real(j,WP) - 0.5_WP*real(ny+1,WP))*elem + y_c
               z0 = (real(k,WP) - 0.5_WP*real(nz+1,WP))*elem
               if (((x0-x_c)*(x0-x_c) + (y0-y_c)*(y0-y_c)).ge.R*R) cycle;
               if (((x0-x_c)*(x0-x_c) + (y0-y_c)*(y0-y_c)).le.R_int*R_int) cycle;
               if (((x0-x_c).gt.0.0_WP).and.(abs(y0-y_c).le.Beam_H/2.0_WP)) cycle;
               n=n+1
               pos(:,n)=[x0, y0, z0]
               vel(:,n)=[0.0_WP, 0.0_WP, 0.0_WP]
               flags(n)= PDC_BONDS ! PDC_MOVES+PDC_INTEGRATES+PDC_BONDS !< IVM, bitwise, this should keep it still?
               if (((x0-x_c)*(x0-x_c) + (y0-y_c)*(y0-y_c)).ge.(R-ls%delta)*(R-ls%delta)) flags(n)= PDC_BONDS + PDC_SURFACE
               gids(n)=int(n,I8)
               voll(n)=elem**3
            end do; end do; end do
            print*, n

         
            do k=1,nz; do j=1,bny; do i=1,bnx
               if (.not.amRoot) exit
               x0 = (real(i,WP) - 0.5_WP*real(bnx+1,WP))*elem + Beam_L/2.0_WP + x_c
               y0 = (real(j,WP) - 0.5_WP*real(bny+1,WP))*elem + y_c
               z0 = (real(k,WP) - 0.5_WP*real(nz+1,WP))*elem
               ! if (((x0)*(x0) + y0*y0).le.R*R) cycle;
               n=n+1
               pos(:,n)=[x0, y0, z0]
               vel(:,n)=[0.0_WP, 0.0_WP, 0.0_WP]
               flags(n)= PDC_BONDS !PDC_MOVES+PDC_INTEGRATES+PDC_BONDS !< IVM, bitwise, this should keep it still?
               if ((x0-x_c).gt.R) flags(n)=PDC_MOVES+PDC_INTEGRATES+PDC_BONDS
               if((j.le.3).and.((x0-x_c).gt.(R-ls%delta))) flags(n) = PDC_MOVES+PDC_INTEGRATES+PDC_BONDS+ PDC_SURFACE
               if((j.ge.bny-2).and.((x0-x_c).gt.(R-ls%delta))) flags(n) = PDC_MOVES+PDC_INTEGRATES+PDC_BONDS + PDC_SURFACE
               if((i.ge.bnx-2).and.((x0-x_c).gt.R-ls%delta)) flags(n) = PDC_MOVES+PDC_INTEGRATES+PDC_BONDS + PDC_SURFACE
               gids(n)=int(n,I8)
               voll(n)=elem**3
            end do; end do; end do
            call pd_partition(nn,gids,pos,vel,flags,voll,owner,nr,rgid,rpos,rvel,rflag,rvol)
            call ls%set_nodes(nr,rgid,rpos,rvel,rflag,rvol)
            call ls%detect_families()
         end if
         
         
         

         ! COMMS TEST
         allocate(ls%which_rank(ls%nown))
         ls%which_rank = ls%cfg%rank

         call ls%update_fluid_sync()
         
      end block initialize_lss

     ! Create partmesh object for visualizing Lagrangian particles
      create_pmesh: block
         use lss_class, only: max_bond
         integer :: i,n,nbond
         pmesh=partmesh(nvar=3,nvec=3,name='solid')
         pmesh%varname(1)='damage'
         pmesh%varname(2)='flag'
         pmesh%varname(3)='which_rank'
        
         pmesh%vecname(1)='velocity'
         pmesh%vecname(2)='fluid_force'
         pmesh%vecname(3) = 'displacement'
         call ls%update_partmesh(pmesh)
         
         do i=1,ls%nown 
            pmesh%vec(:,1,i)=ls%v(:,i)
            pmesh%vec(:,2,i)=ls%ff(:,i)
            pmesh%vec(:,3,i)=ls%y(:,i)-ls%x0(:,i)
            pmesh%var(2,i)=ls%flag(i)
            pmesh%var(3,i)=ls%which_rank(i)
         end do
      end block create_pmesh

      ! Create a flow solver with inflow-outflow
      create_flow_solver: block
         use incomp_class, only: dirichlet,clipped_neumann
         real(WP) :: visc
         ! Create flow solver
         fs=incomp(cfg=cfg,name='Incompressible NS')
         ! Set the flow properties
         call param_read('Density',fs%rho)
         call param_read('Dynamic viscosity',visc); fs%visc=visc
         ! Define boundary conditions
         call fs%add_bcond(name='inflow', type=dirichlet      ,locator=left_of_domain ,face='x',dir=-1,canCorrect=.false.)
         call fs%add_bcond(name='outflow',type=clipped_neumann,locator=right_of_domain,face='x',dir=+1,canCorrect=.true. )
         ! Configure pressure solver
         ps=fft2d(cfg=cfg,name='Pressure',nst=7)
         ! Configure implicit velocity solver
         vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
      end block create_flow_solver
      
      
      ! Initialize our velocity field
      initialize_velocity: block
         use random,       only: random_normal
         use incomp_class, only: bcond
         type(bcond), pointer :: mybc
         integer :: n,i,j,k
         real(WP) :: Uin, Ly
         ! Read inflow velocity
         call param_read('Inlet velocity',Uin)
         call param_read('Ly',Ly)
         ! IB arrays
         Uib=0.0_WP; Vib=0.0_WP; Wib=0.0_WP; srcM=0.0_WP
         ! Make initial velocity field random to trigger transition
         if (restarted) then
            restart_fluid: block
            real(WP) :: step_real
            call df%pullval(name='step', val=step_real)
            call df%pullval(name='ls_dt', val=ls_dt)

            time%n = nint(step_real)
            time%told = time%t - time%dt

            call df%pullvar(name='U', var=fs%U)
            call df%pullvar(name='V', var=fs%V)
            call df%pullvar(name='W', var=fs%W)
            call df%pullvar(name='P', var=fs%P)
            end block restart_fluid
         else
            do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        fs%U(i,j,k)=0.0_WP
                        fs%V(i,j,k)=0.0_WP
                        fs%W(i,j,k)=0.0_WP
                  end do
               end do
            end do
         end if
         call fs%cfg%sync(fs%U)
         call fs%cfg%sync(fs%V)
         call fs%cfg%sync(fs%W)
         ! Set inflow velocity
         call fs%get_bcond('inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%U(i,j,k)=6.0_WP * Uin * cfg%ym(j)*(Ly -cfg%ym(j))/Ly**2 
         end do
         ! Compute MFR through all boundary conditions
         call fs%get_mfr()
         ! Adjust MFR for global mass balance
         call fs%correct_mfr(src=srcM)
         ! Compute cell-centered velocity
        call fs%interp_vel(Ui,Vi,Wi)
         ! Compute divergence
         resU=srcM/fs%rho           !< Careful, we need to provide
         call fs%get_div(src=resU)  !< a volume source term to div
         
      end block initialize_velocity

      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='cylinder')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_particle('particles',pmesh)
         call ens_out%add_scalar('divergence',fs%div)
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_vector('velocity_s',Uib,Vib,Wib)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_scalar('VFs',ls%VF)
         call ens_out%add_scalar('SRCM',srcM)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      ! Create monitor files
      create_monitor: block
        real(WP) :: cfl
        ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
        ! Create simulation monitor
        mfile=monitor(fs%cfg%amRoot,'simulation')
        call mfile%add_column(time%n,'Timestep number')
        call mfile%add_column(time%t,'Time')
        call mfile%add_column(time%dt,'Timestep size')
        call mfile%add_column(time%cfl,'Maximum CFL')
        call mfile%add_column(fs%Umax,'Umax')
        call mfile%add_column(fs%Vmax,'Vmax')
        call mfile%add_column(fs%Wmax,'Wmax')
        call mfile%add_column(fs%Pmax,'Pmax')
        call mfile%add_column(fs%divmax,'Maximum divergence')
        call mfile%add_column(fs%psolv%it,'Pressure iteration')
        call mfile%add_column(fs%psolv%rerr,'Pressure error')
        call mfile%add_column(solid_substeps,'Particle Sub-Stpes')
        call mfile%write()

        ! Create CFL monitor
        cflfile=monitor(fs%cfg%amRoot,'cfl')
        call cflfile%add_column(time%n,'Timestep number')
        call cflfile%add_column(time%t,'Time')
        call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
        call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
        call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
        call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
        call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
        call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
        call cflfile%write()

        ! Create solid monitor
        call ls%get_info()
        sfile=monitor(fs%cfg%amRoot,'solid')
        call sfile%add_column(time%n,'Timestep number')
        call sfile%add_column(time%t,'Time')
        call sfile%add_column(ls_dt,'Particle dt')
        call sfile%add_column(ls%CFLe,'CFLe')
        call sfile%add_column(ls%CFLp,'CFLp')
        ! call sfile%add_column(ls%VFmax,'VFmax')
        call sfile%add_column(ls%Umin,'Particle Umin')
        call sfile%add_column(ls%Umax,'Particle Umax')
        call sfile%add_column(ls%Vmin,'Particle Vmin')
        call sfile%add_column(ls%Vmax,'Particle Vmax')
        call sfile%add_column(ls%Wmin,'Particle Wmin')
        call sfile%add_column(ls%Wmax,'Particle Wmax')       
      !   call sfile%add_column(ls%ibmForce(1),'Particle Fx')
      !   call sfile%add_column(ls%ibmForce(2),'Particle Fy')
      !   call sfile%add_column(ls%ibmForce(3),'Particle Fz')
        call sfile%write()

        ! Create solid timing monitor
        tfile=monitor(amRoot=amRoot,name='timing')
        call tfile%add_column(time%n,'Timestep')
        call tfile%add_column(time%t,'Time')
        call tfile%add_column(ls%np,'Nodes')
        call tfile%add_column(ls%nb,'Bonds')
        call tfile%add_column(ls%wtmax_kick,   'kick_max')
        call tfile%add_column(ls%wtmax_halo,   'halo_max')
        call tfile%add_column(ls%wtmax_dil,    'dil_max')
        call tfile%add_column(ls%wtmin_dil,    'dil_min')
        call tfile%add_column(ls%wtmax_force,  'force_max')
        call tfile%add_column(ls%wtmin_force,  'force_min')
        call tfile%add_column(ls%wtmax_reduce, 'reduce_max')
        call tfile%add_column(ls%wtmax_contact,'contact_max')
        call tfile%add_column(ls%wtmax_broad,  'broad_max')
        call tfile%add_column(ls%maxtot_time,  'total_max')
        call tfile%write()
      end block create_monitor

      print *, '================== simulation_init COMPLETE =================='

    end subroutine simulation_init


   !> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
    subroutine simulation_run
      implicit none
      real(WP) :: solid_cfl

      ! ls%damping_rate = 0.00025
      ! print*, "Restart Value check:"
      ! print*,'Material density: ',ls%rho
      ! print*,'Elastic modulus: ', ls%elastic_modulus 
      ! print*,'Poisson ratio: ',   ls%poisson_ratio
      ! print*,'Particle timestep size: ',ls_dt_max
      ! print*,'Particle Volume: ', ls%dV

      ! Perform time integration
      do while (.not.time%done())
         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()
         solid_substeps=0
      
         solid: block
            real(WP) :: dt_done,mydt
            
            if (.not.steady_state) then
               ! Advance solid solver
               
               ! Compute divergence of fluid stress (old way, currently not used)
               ! call fs%get_div_stress(divx=div_x(:,:,:),divy=div_y(:,:,:),divz=div_z(:,:,:))
               ! Sub-iteratore
               call ls%get_cfl(ls_dt,cfl=solid_cfl)
               if (solid_cfl.gt.0.0_WP) ls_dt=min(ls_dt*time%cflmax/solid_cfl,ls_dt_max)
               dt_done=0.0_WP
               
               do while (dt_done.lt.time%dtmid)
                  ! Decide the timestep size
                  mydt=min(ls_dt,time%dtmid-dt_done)
                  ! Advance particles
                  call ls%fluid_sync(dt_f=time%dtmid,fvel_x=fs%U,fvel_y=fs%V,fvel_z=fs%W) 
                  call ls%advance(dt      =mydt,& 
                  &               fluid_dt=time%dtmid,&
                  &               fluid_rho=fs%rho) !,&
                  !   &               div_stress_x=fs%U(:,:,:),&
                  !   &               div_stress_y=fs%V(:,:,:),&
                  !   &               div_stress_z=fs%W(:,:,:))
                  ! Increment
                  dt_done=dt_done+mydt
                  if(cfg%amRoot) solid_substeps=solid_substeps+1
               end do
               call ls%update_fluid_sync()
               
            end if
         end block solid

         ! ! Evaluate IB velocity and mass source
         ! calc_ib_velocity: block
         !    integer :: i,j,k
         !    do k=fs%cfg%kmin_,fs%cfg%kmax_
         !       do j=fs%cfg%jmin_,fs%cfg%jmax_
         !          do i=fs%cfg%imin_,fs%cfg%imax_
         !             ! VF based velocity
         !             Uib(i,j,k)=0.5_WP*(ls%VFU(i-1,j,k)+ls%VFU(i,j,k))/(sum(fs%itpr_x(:,i,j,k)*ls%VF(i-1:i,j,k))+epsilon(1.0_WP))
         !             Vib(i,j,k)=0.5_WP*(ls%VFV(i,j-1,k)+ls%VFV(i,j,k))/(sum(fs%itpr_y(:,i,j,k)*ls%VF(i,j-1:j,k))+epsilon(1.0_WP))
         !             Wib(i,j,k)=0.5_WP*(ls%VFW(i,j,k-1)+ls%VFW(i,j,k))/(sum(fs%itpr_z(:,i,j,k)*ls%VF(i,j,k-1:k))+epsilon(1.0_WP))
         !          end do
         !       end do
         !    end do
         !    call cfg%sync(Uib)
         !    call cfg%sync(Vib)
         !    call cfg%sync(Wib)
         !    ! Compute IB mass source
         !    do k=fs%cfg%kmin_,fs%cfg%kmax_
         !       do j=fs%cfg%jmin_,fs%cfg%jmax_
         !          do i=fs%cfg%imin_,fs%cfg%imax_
         !             srcM(i,j,k)=fs%rho*(ls%VF(i,j,k)*(sum(fs%divp_x(:,i,j,k)*Uib(i:i+1,j,k))+&
         !             &                                sum(fs%divp_y(:,i,j,k)*Vib(i,j:j+1,k))+&
         !             &                                sum(fs%divp_z(:,i,j,k)*Wib(i,j,k:k+1))))

         !          end do
         !       end do
         !    end do
         !    call cfg%sync(srcM)
         ! end block calc_ib_velocity
         
         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W
         
         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            
            ! Build mid-time velocity
            fs%U=0.5_WP*(fs%U+fs%Uold)
            fs%V=0.5_WP*(fs%V+fs%Vold)
            fs%W=0.5_WP*(fs%W+fs%Wold)
            
            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)
            
            ! Assemble explicit residual
            resU=-2.0_WP*(fs%rho*fs%U-fs%rho*fs%Uold)+time%dtmid*resU
            resV=-2.0_WP*(fs%rho*fs%V-fs%rho*fs%Vold)+time%dtmid*resV
            resW=-2.0_WP*(fs%rho*fs%W-fs%rho*fs%Wold)+time%dtmid*resW
            
            ! Form implicit residuals
            call fs%solve_implicit(time%dtmid,resU,resV,resW)

            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW
            
            ! Apply direct IB forcing
            ! ibforcing: block
            !    integer :: i,j,k
            !    do k=fs%cfg%kmin_,fs%cfg%kmax_; do j=fs%cfg%jmin_,fs%cfg%jmax_; do i=fs%cfg%imin_,fs%cfg%imax_
            !       fs%U(i,j,k)=(1.0_WP-sum(fs%itpr_x(:,i,j,k)*ls%VF(i-1:i,j,k)))*fs%U(i,j,k)+0.5_WP*(ls%VFU(i-1,j,k)+ls%VFU(i,j,k))
            !       fs%V(i,j,k)=(1.0_WP-sum(fs%itpr_y(:,i,j,k)*ls%VF(i,j-1:j,k)))*fs%V(i,j,k)+0.5_WP*(ls%VFV(i,j-1,k)+ls%VFV(i,j,k))
            !       fs%W(i,j,k)=(1.0_WP-sum(fs%itpr_z(:,i,j,k)*ls%VF(i,j,k-1:k)))*fs%W(i,j,k)+0.5_WP*(ls%VFW(i,j,k-1)+ls%VFW(i,j,k))
            !       ! Enforcing no slip on the walls
            !       fs%U(i,j,k)=sum(fs%itpr_x(:,i,j,k)*cfg%VF(i-1:i,j,k))*fs%U(i,j,k)
            !       fs%V(i,j,k)=sum(fs%itpr_y(:,i,j,k)*cfg%VF(i,j-1:j,k))*fs%V(i,j,k)
            !       fs%W(i,j,k)=sum(fs%itpr_z(:,i,j,k)*cfg%VF(i,j,k-1:k))*fs%W(i,j,k)
            !    end do; end do; end do
            !    call fs%cfg%sync(fs%U)
            !    call fs%cfg%sync(fs%V)
            !    call fs%cfg%sync(fs%W)
            ! end block ibforcing
            ibm_correction: block
               integer :: i,j,k
               ! Interpolate to the staggered cells and synchronize
               resU=0.0_WP; resV=0.0_WP; resW=0.0_WP
               call ls%fluid_sync(dt_f=time%dtmid,fvel_x=fs%U,fvel_y=fs%V,fvel_z=fs%W) 
               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        resU(i,j,k)=sum(fs%itpr_x(:,i,j,k)*ls%srcU(i-1:i,j,k))
                        resV(i,j,k)=sum(fs%itpr_y(:,i,j,k)*ls%srcV(i,j-1:j,k))
                        resW(i,j,k)=sum(fs%itpr_z(:,i,j,k)*ls%srcW(i,j,k-1:k))
                        fs%U(i,j,k)=sum(fs%itpr_x(:,i,j,k)*cfg%VF(i-1:i,j,k))*fs%U(i,j,k)
                        fs%V(i,j,k)=sum(fs%itpr_y(:,i,j,k)*cfg%VF(i,j-1:j,k))*fs%V(i,j,k)
                        fs%W(i,j,k)=sum(fs%itpr_z(:,i,j,k)*cfg%VF(i,j,k-1:k))*fs%W(i,j,k)
                     end do
                  end do
               end do
               call fs%cfg%sync(resU)
               call fs%cfg%sync(resV)
               call fs%cfg%sync(resW)
               ! Increment velocity field
               fs%U=fs%U+resU
               fs%V=fs%V+resV
               fs%W=fs%W+resW
            end block ibm_correction
            
            ! Apply other boundary conditions
            call fs%apply_bcond(time%t,time%dtmid)

            ! Solve Poisson equation
            ! call fs%correct_mfr(src=srcM)
            call fs%correct_mfr()
            ! resU=srcM/fs%rho           !< Careful, we need to provide
            ! call fs%get_div(src=resU)  !< a volume source term to div
            call fs%get_div()
            fs%psolv%rhs=-fs%cfg%vol*fs%div*fs%rho/time%dtmid
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)
            
            ! Correct velocity
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%U=fs%U-time%dtmid*resU/fs%rho
            fs%V=fs%V-time%dtmid*resV/fs%rho
            fs%W=fs%W-time%dtmid*resW/fs%rho
            
            ! Increment sub-iteration counter
            time%it=time%it+1
            
         end do
         
         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         ! resU=srcM/fs%rho           !< Careful, we need to provide
         ! call fs%get_div(src=resU)  !< a volume source term to div
         call fs%get_div()

         ! Output to ensight
         if (ens_evt%occurs()) then
            update_pmesh: block
              use lss_class, only: max_bond
              integer :: i,n,nbond
              call ls%update_partmesh(pmesh)
               do i=1,ls%nown ! IVM, probably not the right thing 
                  pmesh%vec(:,1,i)=ls%v(:,i)
                  pmesh%vec(:,2,i)=ls%ff(:,i)
                  pmesh%vec(:,3,i)=ls%y(:,i)-ls%x0(:,i)
                  pmesh%var(2,i)=ls%flag(i)
                  pmesh%var(3,i)=ls%which_rank(i)
               end do
            end block update_pmesh
            call ens_out%write_data(time%t)

            
         end if

         if (save_evt%occurs()) then
            save_restart: block
               character(len=str_medium) :: checkpoint_dir
               character(len=10) :: step_string
               ! Prefix for files
               write(step_string,'(i10.10)') time%n
               checkpoint_dir = 'restart/state_'//step_string
               call ls%write_state(dirname=trim(checkpoint_dir))
               ! Populate df and write it
               call df%pushval(name='t',     val=time%t)
               call df%pushval(name='dt',    val=time%dt)
               call df%pushval(name='step',  val=real(time%n,WP))
               call df%pushval(name='ls_dt', val=ls_dt)
               call df%pushvar(name='U' ,var=fs%U       )
               call df%pushvar(name='V' ,var=fs%V       )
               call df%pushvar(name='W' ,var=fs%W       )
               call df%pushvar(name='P' ,var=fs%P       )
               call df%write(fdata=trim(checkpoint_dir)//'/fluid')
               ! Write particle file
               
            end block save_restart
         end if

         ! ! Perform and output monitoring
         call ls%get_info()
         call fs%get_max()
         call mfile%write()
         call cflfile%write()
         call sfile%write()
         call tfile%write()
         
      end do

 end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker
      
      ! Deallocate work arrays
      deallocate(div_x,div_y,div_z,resU,resV,resW,Ui,Vi,Wi,Uib,Vib,Wib,srcM,gradU)
      
   end subroutine simulation_final
   
   
end module simulation
