! Minimal stand-in for CESM share shr_log_mod (errMsg helper only).
module shr_log_mod
   use shr_kind_mod, only : shr_kind_cx
   implicit none
   private
   public :: shr_log_errMsg
contains
   pure function shr_log_errMsg(file, line) result(msg)
      character(len=*), intent(in) :: file
      integer,          intent(in) :: line
      character(len=shr_kind_cx)   :: msg
      write(msg,'(a,a,a,i0)') 'ERROR in ', trim(file), ' at line ', line
   end function shr_log_errMsg
end module shr_log_mod
