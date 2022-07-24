!*****
!***** TIMING MODULE
!*****

MODULE mo_timer
   use mo_kind,       ONLY : ni, nr
   use mo_parameters, ONLY : zero
   use mo_mpi,        ONLY : p_wtime, myid
   IMPLICIT NONE
   PRIVATE

   integer(kind=ni), parameter :: max_length_name = 1024

   type, private :: t_timer
      character(len=max_length_name) :: timer_name
      real(kind=nr)    :: timer1
      real(kind=nr)    :: timer2
      real(kind=nr)    :: timer_value
   end type t_timer

   integer(kind=ni), parameter :: max_timers = 32
   integer(kind=ni) :: current_timers = 0
   type(t_timer), dimension(1:max_timers) :: timers

   integer(kind=ni), public :: timer_filter

   public :: timer_start, timer_stop, timer_init, timer_print

   contains

   function add_timer(timer_name) result(timer_id)
      character(len=*), intent(in) :: timer_name
      integer(kind=ni) :: timer_id

      current_timers = current_timers + 1

      if (current_timers > max_timers) then
         write(*,*) 'Maximum number of timers exceeded'
         STOP 1
      end if

      timers(current_timers)%timer_name = trim(timer_name)
      timers(current_timers)%timer_value = zero
      
      timer_id = current_timers

   end function add_timer

   subroutine timer_start(timer_id)
      integer(kind=ni), intent(in) :: timer_id

      if (timer_id > current_timers) then
         write(*,*) 'Provided timer ID wrong'
         STOP 1
      end if

      timers(timer_id)%timer1 = p_wtime()

   end subroutine timer_start

   subroutine timer_stop(timer_id)
      integer(kind=ni), intent(in) :: timer_id

      if (timer_id > current_timers) then
         write(*,*) 'Provided timer ID wrong'
         STOP 1
      end if

      timers(timer_id)%timer2 = p_wtime()
      timers(timer_id)%timer_value = timers(timer_id)%timer2 - &
                                     timers(timer_id)%timer1

   end subroutine timer_stop

   subroutine timer_init()
      timer_filter = add_timer('Filters')

   end subroutine timer_init

   subroutine timer_print()
      integer(kind=ni) :: i

      if (myid == 0) then
         write(*,*) '----------'
         write(*,*) '--TIMING--'
         write(*,*) '----------'         
         do i = 1,current_timers
            write(*,"(a20, ' = ', f12.5, ' s')") trim(timers(i)%timer_name), timers(i)%timer_value
         end do
      end if

   end subroutine timer_print

END MODULE mo_timer