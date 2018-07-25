module groupr_header

  use constants,      only: MAX_WORD_LEN, DEPLETION_STRINGS
  use hdf5_interface

  implicit none

!===============================================================================
! REACTION contains the cross-section and secondary energy and angle
! distributions for a single reaction in a continuous-energy ACE-format table
!===============================================================================

  type GrouprXS
    integer :: threshold             ! Energy grid index of threshold
    real(8), allocatable :: value(:) ! Cross section values
  end type GrouprXS

  type Groupr
    type(GrouprXS) :: xs
  contains
    procedure :: from_hdf5 => groupr_from_hdf5
  end type Groupr

contains

  subroutine groupr_from_hdf5(this, rx_group)
    class(Groupr), intent(inout) :: this
    integer(HID_T),  intent(in)    :: rx_group

    integer(HID_T) :: xs
    integer(HSIZE_T) :: dims(1)
    integer :: i 
    
    ! Read cross section and threshold_idx data
    xs = open_dataset(rx_group, 'groupr')
    call read_attribute(this % xs % threshold, xs, 'start_point')
    call get_shape(xs, dims)
    allocate(this % xs % value(dims(1)))
    call read_dataset(this % xs % value, xs)
    call close_dataset(xs)
  end subroutine groupr_from_hdf5

end module groupr_header
