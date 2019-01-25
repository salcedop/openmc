module crass

  use, intrinsic :: ISO_C_BINDING
  use constants,      only: MAX_WORD_LEN
  use hdf5_interface
  use message_passing,     only: n_procs
  implicit none

  public :: openmc_MG_rates


real(8),allocatable,target :: group_tally_results(:,:,:)
contains

  function openmc_MG_rates(ptr, shape_) result(err) bind(C)
    ! Returns a pointer to a MG results array along with its shape. This
    ! allows a user to obtain in-memory tally results from Python directly.
    type(C_PTR),        intent(out) :: ptr
    integer(C_INT),     intent(out) :: shape_(3)
    integer(C_INT) :: err

    if (allocated(group_tally_results)) then
      ptr = C_LOC(group_tally_results(1,1,1))
      shape_(:) = shape(group_tally_results)
      err = 0
    else
      err = 100!E_ALLOCATE
      !call set_errmsg("Tally results have not been allocated yet.")
    end if
  end function openmc_MG_rates

end module crass
