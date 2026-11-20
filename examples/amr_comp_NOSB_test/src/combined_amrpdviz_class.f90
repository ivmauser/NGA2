!> Particle visualization appended to the matching native AMReX mesh plotfile.
!>
!> AMReX ParticleContainer%WritePlotFile expects its output directory to
!> already be a valid mesh plotfile.  The base amrpdviz class instead creates
!> standalone plt.part.NNNNNN directories, which leaves them without the
!> top-level mesh Header required by ParaView.  This case-local extension keeps
!> the base initialization/component-selection machinery and only redirects
!> particle output into plt.nga2.cell.NNNNNN after amrviz has written it.
module combined_amrpdviz_class
   use precision,       only: WP
   use string,          only: str_long
   use amrpdviz_class,  only: amrpdviz
   use filesys,         only: isfile
   use messager,        only: die
   use mpi_f08,         only: MPI_BCAST,MPI_LOGICAL
   use iso_c_binding,   only: c_ptr,c_char,c_int,c_double,c_null_char
   implicit none
   private

   public :: combined_amrpdviz

   interface
      subroutine amrpd_write_plotfile(pc,basedir,pname,write_real,write_int,time) bind(c)
         import :: c_ptr,c_char,c_int,c_double
         type(c_ptr), value :: pc
         character(kind=c_char) :: basedir(*),pname(*)
         integer(c_int), intent(in) :: write_real(*),write_int(*)
         real(c_double), value :: time
      end subroutine amrpd_write_plotfile
   end interface

   type, extends(amrpdviz) :: combined_amrpdviz
   contains
      procedure :: write => write_combined
   end type combined_amrpdviz

contains

   !> Append particles to the cell-centered plotfile written by amrviz.
   subroutine write_combined(this,time)
      implicit none
      class(combined_amrpdviz), intent(inout) :: this
      real(WP), intent(in) :: time

      character(len=str_long) :: pltdir
      real(WP), allocatable :: tmp(:)
      integer :: i,n,ierr
      logical :: header_ok

      ! Match amrviz rewind/restart indexing exactly.
      if (this%ntime.eq.0) then
         this%ntime=1
         if (allocated(this%time)) deallocate(this%time)
         allocate(this%time(1))
         this%time(1)=time
      else
         n=1
         rewind: do i=this%ntime,1,-1
            if (this%time(i).lt.time-1.0e-6_WP) then
               n=i+1
               exit rewind
            end if
         end do rewind
         this%ntime=n
         allocate(tmp(this%ntime))
         tmp=[this%time(1:this%ntime-1),time]
         call move_alloc(tmp,this%time)
      end if

      ! amrviz writes this native mesh plotfile immediately before this call.
      write(pltdir,'("amrviz/",a,"/plt.nga2.cell.",i6.6)') trim(this%name),this%ntime
      ! Only root probes the root-written Header.  The broadcast also holds
      ! other ranks here until the mesh writer has made the file visible.
      header_ok=.false.
      if (this%pd%amr%amRoot) header_ok=isfile(trim(pltdir)//'/Header')
      call MPI_BCAST(header_ok,1,MPI_LOGICAL,0,this%pd%amr%comm,ierr)
      if (.not.header_ok) &
      &  call die('[combined_amrpdviz] mesh plotfile Header is missing; call amrviz%write before particle write')
      call amrpd_write_plotfile(this%pd%pcp,trim(pltdir)//c_null_char, &
      &                         'particles'//c_null_char, &
      &                         this%write_real,this%write_int, &
      &                         real(time,c_double))

      ! Keep a particle-specific ParaView series, but reference the combined
      ! mesh+particle plotfile directories.
      if (this%pd%amr%amRoot) then
         open(newunit=n,file='amrviz/'//trim(this%name)//'/plt.part.series', &
         &    status='replace',action='write')
         write(n,'(a)') '{ "file-series-version": "1.0",'
         write(n,'(a)') '  "files": ['
         do i=1,this%ntime
            write(pltdir,'("plt.nga2.cell.",i6.6)') i
            if (i.lt.this%ntime) then
               write(n,'(4x,a,a,a,es24.17,a)') '{ "name": "',trim(pltdir),'", "time": ',this%time(i),' },'
            else
               write(n,'(4x,a,a,a,es24.17,a)') '{ "name": "',trim(pltdir),'", "time": ',this%time(i),' }'
            end if
         end do
         write(n,'(a)') '  ]'
         write(n,'(a)') '}'
         close(n)
      end if

   end subroutine write_combined

end module combined_amrpdviz_class
