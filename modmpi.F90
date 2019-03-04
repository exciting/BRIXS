! Copyright (C) 2006-2008 C. Ambrosch-Draxl. C. Meisenbichler S. Sagmeister
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details
!
! !MODULE: modmpi
! !DESCRIPTION:
!   MPI variables and interface functions
!   In case of compiled without MPI support it defines the
!   mpi specific variables such that the code behaves exactly as
!   the unmodified scalar version
!
! !REVISION HISTORY:
!   Created October 2006 (CHM)
!   Added wrapper routines, 2007-2008 (S. Sagmeister)
!   Added allgatherv interface, August 2010 (S. Sagmeister)
!   Added subroutines/functions to documentation scheme. 2016 (Aurich)
!   Adapted partitioning functions to handle cases with more processes than elements. 2016 (Aurich)
!   Added proc groups functionality. 2016 (Aurich)
!
module modmpi
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

    !BOP
    ! !ROUTINE: initmpi
    ! !INTERFACE:
    subroutine initmpi
    ! !DESCRIPTION:
    !   Initializes MPI and sets procs and rank numbers.
    !   Sets splittfile and initializes the first in node
    !   list. Or if -DMPI is not used sets procs=1 and rank=1
    !   and splittfile=.false.
    !
    ! !REVISION HISTORY:
    !   Added to documentation scheme (Aurich)
    !EOP
    !BOC
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
    !EOC

    !BOP
    ! !ROUTINE: finitmpi
    ! !INTERFACE: 
    subroutine finitmpi
    ! !DESCRIPTION:
    !   If -DMPI calls {\tt mpi\_finalize}, after wainting
    !   for the processes of mpiglobal, mpinodes%mpi and mpinodes%mpiintercom
    !   to finish communicating.
    !
    !   Note: For other communicators you need to make shure that they finished
    !         communication.
    !
    ! !REVISION HISTORY:
    !   Added to documentation scheme. (Aurich)
    !   Modified to use new mpi types. (Aurich)
    !EOP
    !BOC
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
    !EOC

    !BOP
    ! !ROUTINE: terminate
    ! !INTERFACE: 
    subroutine terminate
    ! !DESCRIPTION:
    !   Kills the program in {\tt MPI} or
    !   single execution.
    ! 
    ! !REVISION HISTORY:
    !   Added to documentation scheme and moved to 
    !   modmpi. 2016 (Aurich)
    !EOP
    !BOC
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
    !EOC

    !+++++++++++++++++++++++++++++++++++++++++++!
    ! Partitioning of N elements to P processes !
    ! in continuous blocks. Each element is     !
    ! associated to one and only one process.   !
    !+++++++++++++++++++++++++++++++++++++++++++!

    !BOP
    ! !ROUTINE: nofset
    ! !INTERFACE:
    function nofset(myrank, set, nprocs)
    ! !INPUT/OUTPUT PARAMETERS:
    ! IN:
    ! integer(4) :: myrank  ! MPI rank
    ! integer(4) :: set     ! Total number of elements to distribute
    ! integer(4), optional :: nprocs ! Number of processes in communicator
    ! OUT:
    ! integer(4) :: nofset  ! Number of elements for that rank
    !
    ! !DESCRIPTION:
    !   This functions helps with distributing a set of $N_\text{el}$ elements 
    !   to $N_\text{p}$ {\tt MPI} processes in continuous blocks. The function calculates
    !   the number of elements $N_\text{el}(p)$ a given process is responsible for. \\
    !   Example:\\
    !   $N_\text{el}=10, N_\text{p}=3 \rightarrow N_\text{el}(0)=4, N_\text{el}(1)=3, N_\text{el}(2)=3$\\
    !   Example:\\
    !   $N_\text{el}=2, N_\text{p}=3 \rightarrow N_\text{el}(0)=1, N_\text{el}(1)=1, N_\text{el}(2)=0$
    !
    ! !REVISION HISTORY:
    !   Added to documentation scheme. (Aurich)
    !   Added sanity checks. (Aurich)
    !EOP
    !BOC
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
    !EOC

    !BOP
    ! !ROUTINE: firstofset
    ! !INTERFACE:
    function firstofset(myrank, set, nprocs)
    ! !INPUT/OUTPUT PARAMETERS:
    ! IN:
    ! integer(4) :: myrank  ! MPI rank
    ! integer(4) :: set     ! Total number of elements to distribute
    ! integer(4), optional :: nprocs  ! Number of processes in commuincator
    ! OUT:
    ! integer(4) :: firstofset ! Index of the total set for the first index 
    !                          ! of the current subset
    !
    ! !DESCRIPTION:
    !   This functions helps with distributing a set of $N_\text{el}$ elements 
    !   to $N_\text{p}$ {\tt MPI} processes in continuous blocks. The function calculates
    !   the index of the fist element $i_\text{el}(p)$ a given process is responsible for.
    !   If there are more processes than elements a process responsible for no element
    !   gets the assignment $i_\text{el}(p > N_\text{el}-1) = 0$.\\
    !   Example:\\
    !   $N_\text{el}=10, N_\text{p}=3 \rightarrow i_\text{el}(0)=1, i_\text{el}(1)=5, i_\text{el}(2)=8$\\
    !   Example:\\
    !   $N_\text{el}=2, N_\text{p}=4 \rightarrow i_\text{el}(0)=1, i_\text{el}(1)=2, i_\text{el}(2)=0, i_\text{el}(3)=0$
    !
    ! !REVISION HISTORY:
    !   Added to documentation scheme. (Aurich)
    !   Changed behaviour if there are more processes than elements. (Aurich)
    !   Added sanity checks. (Aurich)
    !EOP
    !BOC
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
    !EOC

    !BOP
    ! !ROUTINE: lastofset
    ! !INTERFACE:
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
    !   This functions helps with distributing a set of $N_\text{el}$ elements 
    !   to $N_\text{p}$ {\tt MPI} processes in continuous blocks. The function calculates
    !   the index of the last element $j_\text{el}(p)$ a given process is responsible for.
    !   If there are more processes than elements, a process responsible for no element
    !   gets the assignment $j_\text{el}(p > N_\text{el}-1) = -1$.\\
    !   Example:\\
    !   $N_\text{el}=10, N_\text{p}=3 \rightarrow j_\text{el}(0)=4, j_\text{el}(1)=7, j_\text{el}(2)=10$\\
    !   Example:\\
    !   $N_\text{el}=2, N_\text{p}=4 \rightarrow j_\text{el}(0)=1, j_\text{el}(1)=2, j_\text{el}(2)=-1, j_\text{el}(3)=-1$
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
    !EOC

    !BOP
    ! !ROUTINE: procofindex
    ! !INTERFACE:
    function procofindex(k, set, nprocs)
    ! !INPUT/OUTPUT PARAMETERS:
    ! IN:
    ! integer(4) :: k    ! Element number k
    ! integer(4) :: set  ! Total number of distributed elements 
    ! integer(4), optional :: nprocs ! Number of processes in communicator
    ! OUT:
    ! integer(4) :: procofindex  ! Rank that holds the element
    !
    ! !DESCRIPTION:
    !   This functions helps with distributing a set of $N_\text{el}$ elements 
    !   to $N_\text{p}$ {\tt MPI} processes in continuous blocks. The function calculates
    !   the index of the process $i_\text{p}(k)$ that is responsible for the element with index $k$.
    !   If $k$ is larger than $N_\text{el}$ or smaller than $1$ the routine returns terminates the execution.
    !   Example:\\
    !   $N_\text{el}=10, N_\text{p}=3 \rightarrow i_\text{p}(1)=1, i_\text{p}(4)=1, i_\text{p}(5)=2, \dots $
    !
    ! !REVISION HISTORY:
    !   Added to documentation scheme. (Aurich)
    !   Adapted to changes in lastofset. (Aurich)
    !   Changed behaviour if k is smaller of larger than the set. (Aurich)
    !   Added sanity checks. (Aurich)
    !EOP
    !BOC
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
    !EOC

    !BOP
    ! !ROUTINE: lastproc
    ! !INTERFACE:
    function lastproc(col, set, nprocs)
    ! !INPUT/OUTPUT PARAMETERS:
    ! IN:
    ! integer(4) :: col ! ``Column" of process grid (i.e. an element index of rank 0)
    ! integer(4) :: set ! Total number of distributed elements. 
    ! integer(4), optional :: nprocs  ! Number of processes in communicator
    ! OUT:
    ! integer(4) :: lastproc ! Number of processes active in process column
    !
    ! !DESCRIPTION:
    !   This functions helps with collecting a set of $N_\text{el}$ elements which were
    !   distributed to $N_\text{p}$ {\tt MPI} processes in continuous blocks and is used for example for
    !   the writing of {\tt PMATXS.OUT}, {\tt EMAT.OUT}, {\tt SCCLI.OUT} 
    !   and {\tt EXCLI.OUT}.
    !   For further describing the functionality of this routine, let us consider an example
    !   distribution: \\
    !   Let $N_\text{el} = 13$ and $N_\text{p} = 5$ : \\
    !   \begin{tabular}{c|ccc|c}
    !     rank & firstofset & \dots & lastofset & nofset \\
    !     \hline
    !     0 & 1 & 2 & 3 & 3 \\
    !     1 & 4 & 5 & 6 & 3 \\
    !     2 & 7 & 6 & 9 & 3 \\
    !     3 & 10 & 11 & - & 2 \\
    !     4 & 12 & 13 & - & 2 \\ 
    !   \end{tabular}
    !
    ! For inputs of $\text{col}=\{1,2,3\}, \text{set}=13$ the routine returns $\{4,4,2\}$, i.e. the process index
    ! of the last active process in the respective column. For all other input for col the routine halts execution.
    ! In the pathological case, where we have more processes than elements the
    ! following example depicts the routines behaviour:\\
    !   Let $N_\text{el} = 3$ and $N_\text{p} = 5$ : \\
    !   \begin{tabular}{c|ccc|c}
    !     rank & firstofset & \dots & lastofset & nofset \\
    !     \hline
    !     0 & 1 & 1 & - & 1 \\
    !     1 & 2 & 2 & - & 1 \\
    !     2 & 3 & 3 & - & 1 \\
    !     3 & 0 & -1 & - & 0 \\
    !     4 & 0 & -1 & - & 0 \\ 
    !   \end{tabular}
    ! For inputs of $\text{col}=1, \text{set}=3$ the routine returns $2$.
    ! For all other input for col execution is halted.
    ! In the other case pathological case, where we have only one processes the
    ! following example depicts the routines behaviour:\\
    !   Let $N_\text{el} = 3$ and $N_\text{p} = 1$ : \\
    !   \begin{tabular}{c|ccc|c}
    !     rank & firstofset & \dots & lastofset & nofset \\
    !     \hline
    !     0 & 1 & 2 & 3 & 3 \\
    !   \end{tabular}
    !
    !
    ! !REVISION HISTORY:
    !   Added to documentation scheme. (Aurich)
    !   Added sanity checks. (Aurich)
    !   
    !EOP
    !BOC
      implicit none
      integer(4) :: lastproc
      integer(4), intent(in) :: col, set
      integer(4), intent(in), optional :: nprocs
      integer(4) :: np

      ! Sanity checks
      if( set < 1 ) then 
        write(*,*) "lastproc (Error): set < 1"
        call terminate
      end if
      if(present(nprocs)) then
        np = nprocs
        if(np < 1) then
          write(*,*) "lastproc (Error): np < 1"
          call terminate
        end if
        if(np > mpiglobal%procs) then
          write(*,*) "lastproc (Error): np > np_max"
          call terminate
        end if
      else
        np = mpiglobal%procs
      end if
      if(col > nofset(0, set, nprocs=np) .or. col < 1) then
        write(*,*) "lastproc (Error): col > nofset(0,set,np) or col < 1"
        call terminate
      end if

      ! Only in the last (or only) column less then all 
      ! processes can be active.
      ! nofset(0,set,np) gives the number of columns.
      if(col /= nofset(0, set, nprocs=np)) then
        lastproc = np - 1
      ! Processes fit evenly (includes only one row
      else if(modulo(set,np) == 0) then
        lastproc = np - 1
      ! Dangling processes
      else
        ! Rest elements not filling last column
        lastproc = modulo(set, np) - 1
      end if
      
    end function lastproc
    !EOC

    !+++++++++++++++++++++++++++++++++++++++++++!
    ! MPI wrapper for "convenience"             !
    !+++++++++++++++++++++++++++++++++++++++++++!

    !BOP
    ! !ROUTINE: barrier
    ! !INTERFACE:
    subroutine barrier(mpicom, callername)
    ! !DESCRIPTION:
    !   If -DMPI calls {\tt mpi\_barrier}, else nothing.
    !
    ! !REVISION HISTORY:
    !   Added to documentation scheme. (Aurich)
    !EOP
    !BOC
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
    !EOC

    !BOP
    ! !ROUTINE: mpi_allgatherv_ifc
    ! !INTERFACE:
    subroutine mpi_allgatherv_ifc(set, rlen, rlenv, ibuf, rlpbuf, rbuf, zbuf,&
      & inplace, comm)
    ! !INPUT/OUTPUT PARAMETERS:
    ! In:
    ! integer(4) :: set              ! Number of elements in distributed set 
    ! integer(4), optional :: rlen   ! Number of data elements per element (constant)
    ! integer(4), optional :: rlenv(set) ! Number of data elements per element
    ! logical, optional :: inplace   ! Use mpi_in_place
    ! type(mpiifo), optional :: comm ! MPI communicator type
    ! In/Out:
    ! integer(4), optional :: ibuf(*) ! Buffers to send/recive
    ! real(4), optional :: rlbuf(*)   ! for different data types
    ! real(8), optional :: rbuf(*)    !
    ! complex(8), optional :: zbuf(*) !
    !
    ! !DESCRIPTION:
    !   Wrapper routine for {\tt MPI\_ALLGATHERV} for different 
    !   data types which is adapted for the k-point set 
    !   distribution scheme. That is this works, if {\it set} number
    !   of elements (e.g. k-points) is distributed over all 
    !   processes in the {\tt MPI} communicator {\it comm} using
    !   a continuous distribution as created by the functions
    !   {\tt nofset, firstofset, lastofset}.
    !   The routine can handle a constant number of data elements per 
    !   set element by specifying {\tt rlen} or a set element dependent 
    !   number of data elements by passing {\tt rlenv(set)}. 
    !
    ! !REVISION HISTORY:
    !   Added to documentation scheme. (Aurich)
    !   Added input parameter for communicator and
    !   a switch for inplace allgather. (Aurich)
    !   Added support for set element dependent number 
    !   of data elements. (Aurich)
    !EOP
    !BOC
      implicit none
      integer(4), intent(in) :: set
      integer(4), intent(in), optional :: rlen
      integer(4), intent(in), optional :: rlenv(set)
      logical, intent(in), optional :: inplace
      type(mpiinfo), intent(in), optional :: comm
      integer(4), intent(inout), optional :: ibuf(*)
      real(4), intent(inout), optional :: rlpbuf(*)
      real(8), intent(inout), optional :: rbuf(*)
      complex(8), intent(inout), optional :: zbuf(*)

      ! Arrays for out of place send
      integer(4), allocatable :: bufi(:)
      real(4), allocatable :: bufrlp(:)
      real(8), allocatable :: bufr(:)
      complex(8), allocatable :: bufz(:)

      type(mpiinfo) :: mpicom
      integer(4) :: ierr
      integer(4), allocatable :: buf_n(:), buf_dspls(:)
      integer(4) :: j
      logical :: ti, tr, trlp, tz, tinplace
      integer(4) :: myrank, myprocs, mycomm

      ! Sanity checks
      if( set < 1 ) then 
        write(*,*) "Error (mpi_allgatherv_ifc): set < 1"
        call terminate
      end if
      if(present(inplace)) then
        tinplace = inplace
      else
        tinplace = .false.
      end if
      if(present(comm)) then
        mpicom = comm
      else
        mpicom = mpiglobal
      end if
      ti = present(ibuf)
      tr = present(rbuf)
      trlp = present(rlpbuf)
      tz = present(zbuf)
      if(count((/ti, tr, trlp, tz/)).ne.1) then
        write(*,*)
        write(*,'("Error (mpi_allgatherv_ifc): Exactly one array must be defined.")')
        write(*,*)
        call terminate
      end if
      if(present(rlen) .and. present(rlenv)&
        & .or. .not. present(rlen) .and. .not. present(rlenv)) then
        write(*,*)
        write(*,'("Error (mpi_allgatherv_ifc): Specifiy either rlen or rlenv")')
        write(*,*)
        call terminate
      end if
      if(present(rlen)) then
        if(rlen < 0) then
          write(*,'("Error (mpi_allgatherv_ifc): rlen < 0")')
          call terminate
        end if
      end if
      if(present(rlenv)) then
        if(any(rlenv < 0)) then
          write(*,'("Error (mpi_allgatherv_ifc): rlenv < 0")')
          call terminate
        end if
      end if

      myrank = mpicom%rank
      myprocs = mpicom%procs
      mycomm = mpicom%comm

#ifdef MPI
      allocate(buf_n(myprocs), buf_dspls(myprocs))

      ! Number of elements in send buffer (flattened array)
      if(present(rlen)) then
        buf_n =(/(rlen*nofset(j, set, myprocs), j = 0, myprocs-1)/)
      else
        do j = 0, myprocs-1
          buf_n(j+1) = sum(rlenv(firstofset(j,set,myprocs):lastofset(j,set,myprocs)))
        end do
      end if

      ! Displacements within receive buffer (flattened array)
      if(present(rlen)) then
        buf_dspls =(/(rlen*(firstofset(j, set, myprocs)-1), j = 0, myprocs-1)/)
      else
        do j = 0, myprocs-1
          buf_dspls(j+1) = sum(buf_n(1:j))
        end do
      end if

      ! Integers
      if(ti) then

        if(.not. tinplace) then
          ! Make send buffer
          allocate(bufi(buf_n(myrank+1)))
          if(buf_n(myrank+1) > 0) then
            bufi(:)= ibuf(buf_dspls(myrank+1)+1:buf_dspls(myrank+1)+buf_n(myrank+1))
          end if
          call mpi_allgatherv(bufi, &
            buf_n(myrank+1), &
            mpi_integer, &
            ibuf, &
            buf_n, &
            buf_dspls, &
            mpi_integer, &
            mycomm, &
            ierr)
          deallocate(bufi)
        else
          ! Use receive buffer as sendbuffer
          call mpi_allgatherv(mpi_in_place, &
            buf_n(myrank+1), &
            mpi_integer, &
            ibuf, &
            buf_n, &
            buf_dspls, &
            mpi_integer, &
            mycomm, &
            ierr)
        end if

      end if

      ! Floats
      if(trlp) then

        if(.not. tinplace) then
          ! Make send buffer
          allocate(bufrlp(buf_n(myrank+1)))
          if(buf_n(myrank+1) > 0) then
            bufrlp(:) =&
              & rlpbuf(buf_dspls(myrank+1)+1:buf_dspls(myrank+1)+buf_n(myrank+1))
          end if
          call mpi_allgatherv(bufrlp, &
            buf_n(myrank+1), &
            mpi_real4, &
            rlpbuf, &
            buf_n, &
            buf_dspls, &
            mpi_real4, &
            mycomm, &
            ierr)
          deallocate(bufrlp)
        else
          ! Use receive buffer as sendbuffer
          call mpi_allgatherv(mpi_in_place, &
            buf_n(myrank+1), &
            mpi_real4, &
            rlpbuf, &
            buf_n, &
            buf_dspls, &
            mpi_real4, &
            mycomm, &
            ierr)
        end if

      end if

      ! Doubles
      if(tr) then

        if(.not. tinplace) then
          ! Make send buffer
          allocate(bufr(buf_n(myrank+1)))
          if(buf_n(myrank+1) > 0) then
            bufr(:)= rbuf(buf_dspls(myrank+1)+1:buf_dspls(myrank+1)+buf_n(myrank+1))
          end if
          call mpi_allgatherv(bufr, &
            buf_n(myrank+1), &
            mpi_double_precision, &
            rbuf, &
            buf_n, &
            buf_dspls, &
            mpi_double_precision, &
            mycomm, &
            ierr)
          deallocate(bufr)
        else
          ! Use receive buffer as sendbuffer
          call mpi_allgatherv(mpi_in_place, &
            buf_n(myrank+1), &
            mpi_double_precision, &
            rbuf, &
            buf_n, &
            buf_dspls, &
            mpi_double_precision, &
            mycomm, &
            ierr)
        end if

      end if

      ! Complex doubles
      if(tz) then

        if(.not. tinplace) then
          ! Make send buffer
          allocate(bufz(buf_n(myrank+1)))
          if(buf_n(myrank+1) > 0) then
            bufz(:)= zbuf(buf_dspls(myrank+1)+1:buf_dspls(myrank+1)+buf_n(myrank+1))
          end if
          call mpi_allgatherv(bufz, &
            buf_n(myrank+1), &
            mpi_double_complex, &
            zbuf, &
            buf_n, &
            buf_dspls, &
            mpi_double_complex, &
            mycomm, &
            ierr)
          deallocate(bufz)
        else
          ! Use receive buffer as sendbuffer
          call mpi_allgatherv(mpi_in_place, &
            buf_n(myrank+1), &
            mpi_double_complex, &
            zbuf, &
            buf_n, &
            buf_dspls, &
            mpi_double_complex, &
            mycomm, &
            ierr)
        end if

      end if

      deallocate(buf_n, buf_dspls)

#endif
    end subroutine
    !EOC

end module modmpi
