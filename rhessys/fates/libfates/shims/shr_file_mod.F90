! Minimal stand-in for CESM share shr_file_mod (unit management only).
module shr_file_mod
   implicit none
   private
   public :: shr_file_getUnit, shr_file_freeUnit
contains
   integer function shr_file_getUnit(unit)
      integer, intent(in), optional :: unit
      logical :: is_open
      if (present(unit)) then
         shr_file_getUnit = unit
      else
         do shr_file_getUnit = 20, 999
            inquire(unit=shr_file_getUnit, opened=is_open)
            if (.not. is_open) return
         end do
         error stop 'shr_file_getUnit: no free units'
      end if
   end function shr_file_getUnit
   subroutine shr_file_freeUnit(unit)
      integer, intent(in) :: unit
      logical :: is_open
      inquire(unit=unit, opened=is_open)
      if (is_open) close(unit)
   end subroutine shr_file_freeUnit
end module shr_file_mod
