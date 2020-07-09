!------------------------------------------------------------------------------
! BRIXS: BSE Calculations for RIXS spectra
!------------------------------------------------------------------------------
!
! MODULE: mod_mpi
!
!> @author
!> Christian Vorwerk, Humboldt Universität zu Berlin.
!
! DESCRIPTION: 
!> MPI variables and interface functions
!> In case of compiled without MPI support it defines the
!> mpi specific variables such that the code behaves exactly as
!> the unmodified scalar version
!
! REVISION HISTORY:
! 01 10 2006 - Created (CHM)
! 01 01 2007 - Added wrapper routines. (S. Sagmeister)
! 01 08 2010 - Added allgatherv interface. (S. Sagmeister)
! 01 01 2016 - Added subroutines/functions to documentation scheme. (Aurich)
! 01 01 2016 - Adapted partitioning functions to handle cases with more processes than elements. (Aurich)
! 01 01 2016 - Added proc groups functionality. (Aurich)
! 01 01 2018 - Adapted from exciting code. (Christian Vorwerk)
! 09 07 2020 - Added documentation. (Christian Vorwerk)
!------------------------------------------------------------------------------
module mod_mpi
#ifdef MPI
  use mpi
#endif

  implicit none

  ! MPI info type
  ! Contains basic information regarding 
  ! a mpi communicator 
  type mpiinfo
    integer(4) :: rank
    integer(4) :: procs
    integer :: comm
    integer(4) :: ierr
  end type mpiinfo

  ! Groups of MPI communicators
  ! connected via inter-communicator
  type procgroup
    ! Total number of process groups
    ! this group belongs to
    integer(4) :: ngroups  
    ! Group id
    integer(4) :: id
    ! MPI information for current
    ! process group
    type(mpiinfo) :: mpi
    ! Inter-groups communicator
    type(mpiinfo) :: mpiintercom
  end type procgroup

  ! mpiinfo for global scope
  type(mpiinfo) :: mpiglobal

  ! Nodes as procgroup
  type(procgroup) :: mpinodes

  !! Legacy code:
  ! Variables (contained in mpiglobal type)
  integer(4) :: rank
  integer(4) :: procs
  integer(4) :: ierr
  ! Variables (contained in mpinodes)
  integer(4) :: firstinnode_comm

  ! Some parts use these 
  logical :: splittfile, firstinnode

  contains

    !+++++++++++++++++++++++++++++++++++++++++!
    ! Initialization and finalization  of MPI !
    !+++++++++++++++++++++++++++++++++++++++++!

    !---------------------------------------------------------------------------  
    !> @author 
    !> Christian Vorwerk, Humboldt Universität zu Berlin.
    !
    ! DESCRIPTION: 
    !> Brief description of routine. 
    !> @brief
    !> Initializes MPI and sets procs and rank numbers.
    !> Sets splittfile and initializes the first in node
    !> list. Or if -DMPI is not used sets procs=1 and rank=1
    !> and splittfile=.false.
    !
    ! REVISION HISTORY:
    ! 09 07 2020 - Added documentation 
    !
    !---------------------------------------------------------------------------  
    subroutine initmpi
      integer(4) :: ierr
#ifdef MPI
      ! Initialize MPI and get number of processes 
      ! and current rank
      call mpi_init(ierr)
      call mpi_comm_size(mpi_comm_world, procs, ierr)
      call mpi_comm_rank(mpi_comm_world, rank, ierr)

      ! Set global mpiinfo type
      mpiglobal%procs = procs
      mpiglobal%comm = mpi_comm_world
      mpiglobal%rank = rank
      mpiglobal%ierr = ierr

      ! Each rank writes its own file (use in pars of GS and XS)
      splittfile = .true.
       
#endif
#ifndef MPI
      procs = 1
      rank = 0

      mpiglobal%procs = procs
      mpiglobal%rank = rank
      mpiglobal%comm = 0

      splittfile = .false.
      firstinnode = .true.
#endif
    end subroutine initmpi

    !---------------------------------------------------------------------------  
    !> @author 
    !> Christian Vorwerk, Humboldt Universität zu Berlin.
    !
    ! DESCRIPTION: 
    !> Brief description of routine. 
    !> @brief
    !> If -DMPI calls {\tt mpi\_finalize}, after waiting
    !> for the processes of mpiglobal, mpinodes%mpi and mpinodes%mpiintercom
    !> to finish communicating.
    !> Note: For other communicators you need to make sure that they finished
    !>       communication.
    !
    ! REVISION HISTORY:
    ! 09 07 2020 - Added documentation 
    !
    !---------------------------------------------------------------------------  
    subroutine finitmpi
      integer(4) :: ierr

      ! Wait for everyone to reach this point
      call barrier(mpiglobal)
#ifdef MPI
      call mpi_finalize(ierr)
      if(ierr /= 0) then 
        write(*,*) "Error (finitmpi): ierr =",ierr
      end if
#endif
    end subroutine finitmpi

    !---------------------------------------------------------------------------  
    !> @author 
    !> Christian Vorwerk, Humboldt Universität zu Berlin.
    !
    ! DESCRIPTION: 
    !> Brief description of routine. 
    !> @brief
    !> Kills the program in {\tt MPI} or
    !> single execution.
    !
    ! REVISION HISTORY:
    ! 09 07 2020 - Added documentation 
    !
    !---------------------------------------------------------------------------  
    subroutine terminate
      implicit none
      integer(4) :: ierr
      ! Abort mpi if necessary
#ifdef MPI
      if(mpiglobal%rank == 0) then 
        write(*,*) "Goodbye, cruel world. (terminate)"
      end if
      call mpi_abort(mpi_comm_world, 1, ierr)
      if(ierr .eq. 0) then
         write (*, '(a)') 'MPI abort'
      else
         write (*, '(a)') 'MPI abort with errors - zombie processes might remain!'
      end if
#endif
#ifndef MPI
      write (*, '(a)') 'Abort'
#endif
      ! stop program
      stop
    end subroutine terminate

    !+++++++++++++++++++++++++++++++++++++++++++!
    ! Partitioning of N elements to P processes !
    ! in continuous blocks. Each element is     !
    ! associated to one and only one process.   !
    !+++++++++++++++++++++++++++++++++++++++++++!

    !---------------------------------------------------------------------------  
    !> @author 
    !> Christian Vorwerk, Humboldt Universität zu Berlin.
    !
    ! DESCRIPTION: 
    !> Brief description of routine. 
    !> @brief
    !> This functions helps with distributing a set of \f$ N_\mathrm{el}\f$ elements 
    !> to \f$ N_\mathrm{p} \f$ {\tt MPI} processes in continuous blocks. The function calculates
    !> the number of elements $N_\mathrm{el}(p)$ a given process is responsible for.
    !> Example:\
    !> \f$ N_\mathrm{el}=10, N_\mathrm{p}=3 \rightarrow N_\mathrm{el}(0)=4, N_\mathrm{el}(1)=3, N_\mathrm{el}(2)=3 \f$\
    !> Example:\
    !> \f$N_\mathrm{el}=2, N_\mathrm{p}=3 \rightarrow N_\mathrm{el}(0)=1, N_\mathrm{el}(1)=1, N_\mathrm{el}(2)=0 \f$
    !
    ! REVISION HISTORY:
    ! 09 07 2020 - Added documentation 
    !
    !> @param[in] myrank 
    !> @param[in] set
    !> @param[in] nprocs    
    !---------------------------------------------------------------------------  
    function nofset(myrank, set, nprocs)
      integer(4) :: nofset
      integer(4), intent(in) :: myrank, set
      integer(4), intent(in), optional :: nprocs
      integer(4) :: np

      ! Sanity checks
      if( myrank < 0 ) then 
        write(*,*) "nofset (Error): myrank < 0"
        call terminate
      end if
      if( set < 1 ) then 
        write(*,*) "nofset (Error): set < 1"
        call terminate
      end if
      if(present(nprocs)) then
        np = nprocs
        if(np < 1) then
          write(*,*) "nofset (Error): np < 1"
          call terminate
        end if
        if(np > mpiglobal%procs) then
          write(*,*) "nofset (Error): np > np_max"
          call terminate
        end if
      else
        np = mpiglobal%procs
      end if
      if(np < myrank+1) then
        write(*,*) "nofset (Error): np < myrank+1"
        call terminate
      end if

      ! Compute number of elements on current rank
      nofset = set / np
      if((mod(set, np) > myrank)) nofset = nofset + 1

    end function nofset

    !---------------------------------------------------------------------------  
    !> @author 
    !> Christian Vorwerk, Humboldt Universität zu Berlin.
    !
    ! DESCRIPTION: 
    !> Brief description of routine. 
    !> @brief
    !> This functions helps with distributing a set of \f$ N_\mathrm{el} \f$ elements 
    !> to \f$ N_\mathrm{p} \f$ MPI processes in continuous blocks. The function calculates
    !> the index of the fist element \f$ i_\mathrm{el}(p) \f$ a given process is responsible for.
    !> If there are more processes than elements a process responsible for no element
    !> gets the assignment \f$ i_\mathrm{el}(p > N_\mathrm{el}-1) = 0 \f$.\
    !> Example:\
    !> \f$ N_\mathrm{el}=10, N_\mathrm{p}=3 \rightarrow i_\mathrm{el}(0)=1, i_\mathrm{el}(1)=5, i_\mathrm{el}(2)=8 \f$\
    !> Example:\
    !> \f$ N_\mathrm{el}=2, N_\mathrm{p}=4 \rightarrow i_\mathrm{el}(0)=1, i_\mathrm{el}(1)=2, i_\mathrm{el}(2)=0, i_\mathrm{el}(3)=0 \f$
    !
    ! REVISION HISTORY:
    ! 09 07 2020 - Added documentation 
    !
    !> @param[in] myrank 
    !> @param[in] set
    !> @param[in] nprocs    
    !---------------------------------------------------------------------------  
    function firstofset(myrank, set, nprocs)
      integer(4) :: firstofset
      integer(4), intent(in) :: myrank, set
      integer(4), intent(in), optional :: nprocs
      integer(4) :: i, np

      ! Sanity checks
      if( myrank < 0 ) then 
        write(*,*) "firstofset (Error): myrank < 0"
        call terminate
      end if
      if( set < 1 ) then 
        write(*,*) "firstofset (Error): set < 1"
        call terminate
      end if
      if(present(nprocs)) then
        np = nprocs
        if(np < 1) then
          write(*,*) "firstofset (Error): np < 1"
          call terminate
        end if
        if(np > mpiglobal%procs) then
          write(*,*) "firstofset (Error): np > np_max"
          call terminate
        end if
      else
        np = mpiglobal%procs
      end if
      if(np < myrank+1) then
        write(*,*) "firstofset (Error): np < myrank+1"
        call terminate
      end if

      ! Compute first element index on current rank
      firstofset = 1
      do i = 0, min( myrank-1, set-1)
        firstofset = firstofset + nofset(i, set, nprocs=np)
      end do
      if(set <= myrank) firstofset = 0

    end function firstofset

    !---------------------------------------------------------------------------  
    !> @author 
    !> Christian Vorwerk, Humboldt Universität zu Berlin.
    !
    ! DESCRIPTION: 
    !> Brief description of routine. 
    !> @brief
    !> This functions helps with distributing a set of \f$ N_\mathrm{el} \f$ elements 
    !> to \f$ N_\mathrm{p} \f$ MPI processes in continuous blocks. The function calculates
    !> the index of the last element \f$ j_\mathrm{el}(p) \f$ a given process is responsible for.
    !> If there are more processes than elements, a process responsible for no element
    !> gets the assignment \f$ j_\mathrm{el}(p > N_\mathrm{el}-1) = -1 \f$.\
    !> Example:\
    !> \f$ N_\mathrm{el}=10, N_\mathrm{p}=3 \rightarrow j_\mathrm{el}(0)=4, j_\mathrm{el}(1)=7, j_\mathrm{el}(2)=10 \f$\
    !> Example:\
    !> \f$ N_\mathrm{el}=2, N_\mathrm{p}=4 \rightarrow j_\mathrm{el}(0)=1, j_\mathrm{el}(1)=2, j_\mathrm{el}(2)=-1, j_\mathrm{el}(3)=-1 \f$
    !
    ! REVISION HISTORY:
    ! 09 07 2020 - Added documentation 
    !
    !> @param[in] myrank 
    !> @param[in] set
    !> @param[in] nprocs    
    !---------------------------------------------------------------------------  
    function lastofset(myrank, set, nprocs)
    ! !INPUT/OUTPUT PARAMETERS:
    ! IN:
    ! integer(4) :: myrank  ! MPI rank
    ! integer(4) :: set     ! Total number of elements to distribute
    ! integer(4) :: nprocs  ! Number of processes in commuincator
    ! OUT:
    ! integer(4) :: lastofset  ! Index of the total set for the first index 
    !                          ! of the current subset
    !
    ! !DESCRIPTION:
    !
    ! !REVISION HISTORY:
    !   Added to documentation scheme. (Aurich)
    !   Changed behaviour if there are more processes than elements. (Aurich)
    !   Added sanity checks. (Aurich)
    !EOP
    !BOC
      integer(4) :: lastofset
      integer(4), intent(in) :: myrank, set
      integer(4), intent(in), optional :: nprocs
      integer(4) :: i, np

      ! Sanity checks
      if( myrank < 0 ) then 
        write(*,*) "lastofset (Error): myrank < 0"
        call terminate
      end if
      if( set < 1 ) then 
        write(*,*) "lastofset (Error): set < 1"
        call terminate
      end if
      if(present(nprocs)) then
        np = nprocs
        if(np < 1) then
          write(*,*) "lastofset (Error): np < 1"
          call terminate
        end if
        if(np > mpiglobal%procs) then
          write(*,*) "lastofset (Error): np > np_max"
          call terminate
        end if
      else
        np = mpiglobal%procs
      end if
      if(np < myrank+1) then
        write(*,*) "lastofset (Error): np < myrank+1"
        call terminate
      end if

      ! Compute last element index on this rank
      lastofset = 0
      do i = 0, min(myrank, set-1)
         lastofset = lastofset + nofset(i, set, nprocs=np)
      end do
      if(set <= myrank) lastofset = -1

    end function lastofset

    !---------------------------------------------------------------------------  
    !> @author 
    !> Christian Vorwerk, Humboldt Universität zu Berlin.
    !
    ! DESCRIPTION: 
    !> Brief description of routine. 
    !> @brief
    !> This functions helps with distributing a set of \f$ N_\mathrm{el} \f$ elements 
    !> to \f$ N_\mathrm{p} \f$ MPI processes in continuous blocks. The function calculates
    !> the index of the process \f$ i_\mathrm{p}(k) \f$ that is responsible for the element with index \f$ k \f$.
    !> If $k$ is larger than \f$ N_\mathrm{el} \f$ or smaller than $1$ the routine returns terminates the execution.
    !> Example:\
    !> \f$ N_\mathrm{el}=10, N_\mathrm{p}=3 \rightarrow i_\mathrm{p}(1)=1, i_\mathrm{p}(4)=1, i_\mathrm{p}(5)=2, \dots \f$
    !
    ! REVISION HISTORY:
    ! 09 07 2020 - Added documentation 
    !
    !> @param[in] k 
    !> @param[in] set
    !> @param[in] nprocs    
    !---------------------------------------------------------------------------  
    function procofindex(k, set, nprocs)
      integer(4) :: procofindex
      integer(4), intent(in) :: k, set
      integer(4), intent(in), optional :: nprocs
      integer(4) :: iproc, np

      ! Sanity checks
      if( k < 1 ) then 
        write(*,*) "procofindex (Error): k < 1"
        call terminate
      end if
      if( k > set ) then 
        write(*,*) "procofindex (Error): set < k"
        call terminate
      end if
      if( set < 1 ) then 
        write(*,*) "procofindex (Error): set < 1"
        call terminate
      end if
      if(present(nprocs)) then
        np = nprocs
        if(np < 1) then
          write(*,*) "procofindex (Error): np < 1"
          call terminate
        end if
        if(np > mpiglobal%procs) then
          write(*,*) "procofindex (Error): np > np_max"
          call terminate
        end if
      else
        np = mpiglobal%procs
      end if

      ! Compute rank that holds element k
      procofindex = 0
      do iproc = 0, np - 1
        if(k > lastofset(iproc, set, nprocs=np)&
          & .and. lastofset(iproc, set, nprocs=np) > 0) procofindex = procofindex + 1
      end do

    end function procofindex

    !+++++++++++++++++++++++++++++++++++++++++++!
    ! MPI wrapper for "convenience"             !
    !+++++++++++++++++++++++++++++++++++++++++++!

    !---------------------------------------------------------------------------  
    !> @author 
    !> Christian Vorwerk, Humboldt Universität zu Berlin.
    !
    ! DESCRIPTION: 
    !> Brief description of routine. 
    !> @brief
    !> If -DMPI calls "mpi_barrier", else nothing.
    !
    ! REVISION HISTORY:
    ! 09 07 2020 - Added documentation 
    !
    !> @param[in] mpicom 
    !> @param[in] callername
    !---------------------------------------------------------------------------  
    subroutine barrier(mpicom, callername)
      implicit none
      type(mpiinfo), intent(in), optional :: mpicom
      character(*), intent(in), optional :: callername
      
      character(*), parameter :: thisname = "barrier"

      type(mpiinfo) :: mpinf

      if(present(mpicom)) then 
        mpinf = mpicom
      else
        mpinf = mpiglobal
      end if

      if(present(callername)) then 
        if(.false.) then 
          write(*, '("Info(",a,"): Rank ",i3," of mpicom", i16," called barrier from ", a)')&
            & trim(thisname), mpinf%rank, mpinf%comm, trim(callername)
        end if
      end if

      ! do nothing if only one process
#ifndef MPI
      if(mpinf%procs .eq. 1) return
#endif
      ! call the mpi barrier
#ifdef MPI
      call mpi_barrier(mpinf%comm, mpinf%ierr)
      if(mpinf%ierr /= 0) then 
        write(*,*) "Error (barrier): ierr =", mpinf%ierr
      end if
#endif
    end subroutine barrier

end module mod_mpi
