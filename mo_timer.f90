!*****
!***** TIMING MODULE
!*****

MODULE mo_timer
   use mo_kind,       ONLY : ni, nr
   use mo_parameters, ONLY : zero
   use mo_mpi,        ONLY : p_wtime, p_get_process_ID, p_get_n_processes, &
                           & p_max, p_min, p_sum
   IMPLICIT NONE
   PRIVATE

   integer(kind=ni), parameter :: max_length_name = 1024

   type, private :: t_timer
      character(len=max_length_name) :: timer_name
      real(kind=nr)    :: timer1
      real(kind=nr)    :: timer2
      real(kind=nr)    :: timer_value
      real(kind=nr)    :: timer_max
      real(kind=nr)    :: timer_min
      real(kind=nr)    :: timer_avg
   end type t_timer

   integer(kind=ni), parameter :: max_timers = 32
   integer(kind=ni) :: current_timers = 0
   type(t_timer), dimension(1:max_timers) :: timers

   integer(kind=ni), public :: timer_filter, timer_loop, timer_timestep, &
                               timer_VSstress, timer_fluxes,timer_GCBC,  &
                               timer_averaging, timer_recording,         &
                               timer_tot_output, timer_output,           &
                               timer_total

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
      timer_loop       = add_timer('Time Loop')
      timer_filter     = add_timer('Filters')
      timer_timestep   = add_timer('Calc Time Step')
      timer_VSstress   = add_timer('VS stress')
      timer_fluxes     = add_timer('Fluxes')
      timer_GCBC       = add_timer('GCBC')
      timer_averaging  = add_timer('Averaging')
      timer_recording  = add_timer('Recording')
      timer_tot_output = add_timer('Total Output')
      timer_output     = add_timer('Output')
      timer_total      = add_timer('Total')

   end subroutine timer_init

   subroutine timer_print()
      integer(kind=ni) :: i, myid, npro
      character(len=20) :: t1
      character(len=12) :: t2

      t1 = '--------------------'
      t2 = '------------'

      npro = p_get_n_processes()

      do i = 1,current_timers
         call p_max(timers(i)%timer_value, timers(i)%timer_max)
         call p_min(timers(i)%timer_value, timers(i)%timer_min)
         call p_sum(timers(i)%timer_value, timers(i)%timer_avg)
         timers(i)%timer_avg = timers(i)%timer_avg / npro
      end do

      myid = p_get_process_ID()
      if (myid == 0) then
         write(*,*) '----------'
         write(*,*) '--TIMING--'
         write(*,*) '----------'
         write(*,"(a20, ' | ', a12, ' | ', a12, ' | ', a12, ' | ', a12, ' | ')") &
                 "Name", "Min [s]", "Avg [s]", "Max [s]", "%" 
         write(*,"(a20, '---', a12, '---', a12, '---', a12, '---', a12, '-- ')") &
                 t1, t2, t2, t2, t2
         do i = 1,current_timers
            write(*,"(a20, ' | ', f12.5, ' | ', f12.5, ' | ', f12.5, ' | ', f12.5, ' | ')") &
                      trim(timers(i)%timer_name), timers(i)%timer_min, timers(i)%timer_avg, &
                      timers(i)%timer_max, timers(i)%timer_avg / timers(timer_total)%timer_avg * 100.0_nr
         end do
      end if

   end subroutine timer_print

END MODULE mo_timer