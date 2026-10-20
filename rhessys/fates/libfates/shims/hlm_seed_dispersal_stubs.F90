! Stand-ins for the CLM/ELM modules used only by FATES's
! DetermineGridCellNeighbors (cross-gridcell seed dispersal, MPI). RHESSys
! runs FATES with seeddisp_cadence = 0 and a single process, so these only
! need to compile; the MPI calls abort if ever reached (rhessys_mpi_stubs.c).
module decompMod
   implicit none
   type :: processor_type
      integer :: begg = 1
      integer :: ncells = 0
   end type processor_type
   type(processor_type), public :: procinfo
end module decompMod

module domainMod
   use shr_kind_mod, only : r8 => shr_kind_r8
   implicit none
   type :: domain_type
      real(r8), pointer :: latc(:) => null()
      real(r8), pointer :: lonc(:) => null()
   end type domain_type
   type(domain_type), public :: ldomain
end module domainMod

module spmdMod
   implicit none
   integer, parameter, public :: MPI_INTEGER = 1
   integer, parameter, public :: MPI_REAL8   = 2
   integer, public :: mpicom     = 0
   integer, public :: npes       = 1
   integer, public :: iam        = 0
   logical, public :: masterproc = .true.
end module spmdMod

module perf_mod
   implicit none
   public :: t_startf, t_stopf
contains
   subroutine t_startf(name)
      character(len=*), intent(in) :: name
   end subroutine t_startf
   subroutine t_stopf(name)
      character(len=*), intent(in) :: name
   end subroutine t_stopf
end module perf_mod
