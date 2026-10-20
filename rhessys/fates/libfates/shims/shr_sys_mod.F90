! Minimal stand-in for CESM share shr_sys_mod: abort prints and stops.
module shr_sys_mod
   use, intrinsic :: iso_fortran_env, only : error_unit
   implicit none
   private
   public :: shr_sys_abort, shr_sys_flush
contains
   subroutine shr_sys_abort(string, rc)
      character(len=*), intent(in), optional :: string
      integer,          intent(in), optional :: rc
      if (present(string)) write(error_unit,'(2a)') 'FATES ABORT: ', trim(string)
      if (present(rc))     write(error_unit,'(a,i0)') 'FATES ABORT rc=', rc
      flush(error_unit)
      error stop 1
   end subroutine shr_sys_abort
   subroutine shr_sys_flush(unit)
      integer, intent(in) :: unit
      flush(unit)
   end subroutine shr_sys_flush
end module shr_sys_mod
