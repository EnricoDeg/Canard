!*****
!***** MPI MODULE
!*****

MODULE mo_mpi

  use mpi
  use mainvar3d
  IMPLICIT NONE
  PRIVATE

  integer(kind=ni),dimension(:,:),allocatable :: ista
  integer(kind=ni),dimension(:),allocatable :: ireq
  integer(kind=ni) :: ir,mpro,npro,myid,info,icom,ierr

  INTEGER :: p_status(MPI_STATUS_SIZE) 

  INTEGER, PARAMETER :: nerr = 0

  PUBLIC :: p_send, p_recv, p_bcast, p_sum, p_isend, p_irecv
  PUBLIC :: p_start, p_stop, p_null_req, p_waitall, p_barrier
  PUBLIC :: mpro, npro, myid

  INTERFACE p_send
    MODULE PROCEDURE p_send_int
    MODULE PROCEDURE p_send_long_int
    MODULE PROCEDURE p_send_real_1d_sp
  END INTERFACE p_send

  INTERFACE p_recv
    MODULE PROCEDURE p_recv_int
    MODULE PROCEDURE p_recv_long_int
    MODULE PROCEDURE p_recv_real_1d_sp
  END INTERFACE p_recv

  INTERFACE p_bcast
    MODULE PROCEDURE p_bcast_int
    MODULE PROCEDURE p_bcast_int_1d
  END INTERFACE p_bcast

  INTERFACE p_sum
    MODULE PROCEDURE p_sum_all_dp_0d
  END INTERFACE

  INTERFACE p_isend
     MODULE PROCEDURE p_isend_real_1d
     MODULE PROCEDURE p_isend_real_1d_sp
     MODULE PROCEDURE p_isend_real_2d
  END INTERFACE

  INTERFACE p_irecv
     MODULE PROCEDURE p_irecv_real_1d
     MODULE PROCEDURE p_irecv_real_1d_sp
     MODULE PROCEDURE p_irecv_real_2d
  END INTERFACE

  CONTAINS

  SUBROUTINE p_start

    ! Basic MPI initialization
    call MPI_INIT(ierr)
    IF (ierr /= MPI_SUCCESS) THEN
      write(nerr,'(a)') ' MPI_INIT failed.'
      write(nerr,'(a,i4)') ' Error =  ', ierr
      STOP 1
    END IF

    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    IF (ierr /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_COMM_RANK failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', ierr
       CALL p_abort
    END IF

    call MPI_COMM_SIZE(MPI_COMM_WORLD,npro,ierr)
    IF (ierr /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_COMM_SIZE failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', ierr
       CALL p_abort
    END IF

    mpro=npro-1
    icom=MPI_COMM_WORLD
    info=MPI_INFO_NULL

    ll=max(npro,12)
    allocate(ista(MPI_STATUS_SIZE,ll),ireq(ll))

  END SUBROUTINE p_start

  SUBROUTINE p_null_req
    ir = 0
  END SUBROUTINE p_null_req

  SUBROUTINE p_barrier(comm)
  
    INTEGER ,INTENT(IN) ,OPTIONAL :: comm
    INTEGER :: com
    
    com = icom
    IF(PRESENT(comm)) com = comm
    CALL MPI_BARRIER (com, ierr)

    IF (ierr /= MPI_SUCCESS) THEN
      WRITE (nerr,'(a,i4,a)') ' MPI_BARRIER on ', myid, ' failed.'
      WRITE (nerr,'(a,i4)') ' Error = ', ierr
      CALL p_abort
    END IF

  END SUBROUTINE p_barrier  

  SUBROUTINE p_send_int(buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: buffer
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = icom
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, MPI_INTEGER4, p_destination, p_tag, &
            p_comm, ierr)
    ELSE
       CALL MPI_SEND (buffer, 1, MPI_INTEGER4, p_destination, p_tag, &
            p_comm, ierr)
    END IF

#ifdef DEBUG
    IF (ierr /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', myid, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', ierr
       CALL p_abort
    END IF
#endif

  END SUBROUTINE p_send_int

  SUBROUTINE p_send_long_int(buffer, p_destination, p_tag, p_count, comm)

    INTEGER(nli), INTENT(in) :: buffer
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = icom
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_SEND (buffer, p_count, MPI_INTEGER8, p_destination, p_tag, &
            p_comm, ierr)
    ELSE
       CALL MPI_SEND (buffer, 1, MPI_INTEGER8, p_destination, p_tag, &
            p_comm, ierr)
    END IF

#ifdef DEBUG
    IF (ierr /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', myid, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', ierr
       CALL p_abort
    END IF
#endif

  END SUBROUTINE p_send_long_int

  SUBROUTINE p_send_real_1d_sp (buffer, p_destination, p_tag, p_count, comm)

    REAL(nsp), INTENT(in) :: buffer(:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = icom
    ENDIF

    IF (PRESENT(p_count)) THEN
      CALL MPI_SEND (buffer, p_count, MPI_REAL4, p_destination, p_tag, &
           p_comm, ierr)
    ELSE
      CALL MPI_SEND (buffer, SIZE(buffer), MPI_REAL4, p_destination, p_tag, &
           p_comm, ierr)
    END IF

#ifdef DEBUG
    IF (ierr /= MPI_SUCCESS) THEN
      WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_SEND from ', myid, &
           ' to ', p_destination, ' for tag ', p_tag, ' failed.'
      WRITE (nerr,'(a,i4)') ' Error = ', ierr
      CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_send_real_1d_sp

  SUBROUTINE p_recv_int(buffer, p_source, p_tag, p_count, comm)
    INTEGER, INTENT(out) :: buffer
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = icom
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, MPI_INTEGER4, p_source, p_tag, &
            p_comm, p_status, ierr)
    ELSE
       CALL MPI_RECV (buffer, 1, MPI_INTEGER4, p_source, p_tag, &
            p_comm, p_status, ierr)
    END IF

#ifdef DEBUG
    IF (ierr /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', myid, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', ierr
       CALL p_abort
    END IF
#endif
  END SUBROUTINE p_recv_int

  SUBROUTINE p_recv_long_int(buffer, p_source, p_tag, p_count, comm)
    INTEGER(nli), INTENT(out) :: buffer
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = icom
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, MPI_INTEGER8, p_source, p_tag, &
            p_comm, p_status, ierr)
    ELSE
       CALL MPI_RECV (buffer, 1, MPI_INTEGER8, p_source, p_tag, &
            p_comm, p_status, ierr)
    END IF

#ifdef DEBUG
    IF (ierr /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', myid, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', ierr
       CALL p_abort
    END IF
#endif
  END SUBROUTINE p_recv_long_int


  SUBROUTINE p_recv_real_1d_sp (buffer, p_source, p_tag, p_count, comm)

    REAL(nsp), INTENT(out) :: buffer(:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = icom
    ENDIF

    IF (PRESENT(p_count)) THEN
       CALL MPI_RECV (buffer, p_count, MPI_REAL4, p_source, p_tag, &
            p_comm, p_status, ierr)
    ELSE
       CALL MPI_RECV (buffer, SIZE(buffer), MPI_REAL4, p_source, p_tag, &
            p_comm, p_status, ierr)
    END IF

#ifdef DEBUG
    IF (ierr /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_RECV on ', myid, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', ierr
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_recv_real_1d_sp

  SUBROUTINE p_bcast_int(buffer, p_source, comm)

    INTEGER, INTENT(inout) :: buffer
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = icom
    ENDIF

    CALL MPI_BCAST (buffer, 1, MPI_INTEGER4, p_source, &
        p_comm, ierr)

#ifdef DEBUG
    IF (ierr /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', ierr
       CALL p_abort
    END IF
#endif

  END SUBROUTINE p_bcast_int

  SUBROUTINE p_bcast_int_1d(buffer, p_source, comm)

    INTEGER, INTENT(inout) :: buffer(:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = icom
    ENDIF

    CALL MPI_BCAST (buffer, SIZE(buffer), MPI_INTEGER4, p_source, &
        p_comm, ierr)

#ifdef DEBUG
     IF (ierr /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a)') ' MPI_BCAST from ', p_source, &
            ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', ierr
       CALL p_abort
    END IF
#endif

  END SUBROUTINE p_bcast_int_1d

  SUBROUTINE p_sum_all_dp_0d(zfield, ssum, comm)

    REAL(nr),          INTENT(out) :: ssum
    REAL(nr),          INTENT(in)  :: zfield
    INTEGER, OPTIONAL, INTENT(in)  :: comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = icom
    ENDIF

    CALL MPI_ALLREDUCE (zfield, ssum, 1, MPI_REAL8, &
         MPI_SUM, p_comm, ierr)
#else
    ssum = zfield
#endif

  END SUBROUTINE p_sum_all_dp_0d

  SUBROUTINE p_isend_real_1d (buffer, p_destination, p_tag, p_count, comm)

    REAL(nr), INTENT(inout) :: buffer(:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = icom
    ENDIF

    IF (PRESENT(p_count)) THEN
       ir = ir + 1
       CALL MPI_ISEND (buffer, p_count, MPI_REAL8, p_destination, p_tag, &
            p_comm, ireq(ir), ierr)
    ELSE
       ir = ir + 1
       CALL MPI_ISEND (buffer, SIZE(buffer), MPI_REAL8, p_destination, p_tag, &
            p_comm, ireq(ir), ierr)
    END IF

#ifdef DEBUG
    IF (ierr /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', myid, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', ierr
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real_1d

  SUBROUTINE p_isend_real_1d_sp (buffer, p_destination, p_tag, p_count, comm)

    REAL(nsp), INTENT(inout) :: buffer(:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = icom
    ENDIF

    IF (PRESENT(p_count)) THEN
       ir = ir + 1
       CALL MPI_ISEND (buffer, p_count, MPI_REAL4, p_destination, p_tag, &
            p_comm, ireq(ir), ierr)
    ELSE
       ir = ir + 1
       CALL MPI_ISEND (buffer, SIZE(buffer), MPI_REAL4, p_destination, p_tag, &
            p_comm, ireq(ir), ierr)
    END IF

#ifdef DEBUG
    IF (ierr /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', myid, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', ierr
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real_1d_sp

  SUBROUTINE p_isend_real_2d (buffer, p_destination, p_tag, p_count, comm)

    REAL(nr), INTENT(inout) :: buffer(:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = icom
    ENDIF

    IF (PRESENT(p_count)) THEN
       ir = ir + 1
       CALL MPI_ISEND (buffer, p_count, MPI_REAL8, p_destination, p_tag, &
            p_comm, ireq(ir), ierr)
    ELSE
       ir = ir + 1
       CALL MPI_ISEND (buffer, SIZE(buffer), MPI_REAL8, p_destination, p_tag, &
            p_comm, ireq(ir), ierr)
    END IF

#ifdef DEBUG
    IF (ierr /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_ISEND from ', myid, &
            ' to ', p_destination, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', ierr
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_isend_real_2d

  SUBROUTINE p_irecv_real_1d (buffer, p_source, p_tag, p_count, comm)

    REAL(nr), INTENT(out) :: buffer(:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = icom
    ENDIF

    IF (PRESENT(p_count)) THEN
       ir = ir + 1
       CALL MPI_IRECV (buffer, p_count, MPI_REAL8, p_source, p_tag, &
            p_comm, ireq(ir), ierr)
    ELSE
       ir = ir + 1
       CALL MPI_IRECV (buffer, SIZE(buffer), MPI_REAL8, p_source, p_tag, &
            p_comm, ireq(ir), ierr)
    END IF

#ifdef DEBUG
    IF (ierr /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', myid, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', ierr
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_real_1d

  SUBROUTINE p_irecv_real_1d_sp (buffer, p_source, p_tag, p_count, comm)

    REAL(nsp), INTENT(out) :: buffer(:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = icom
    ENDIF

    IF (PRESENT(p_count)) THEN
       ir = ir + 1
       CALL MPI_IRECV (buffer, p_count, MPI_REAL4, p_source, p_tag, &
            p_comm, ireq(ir), ierr)
    ELSE
       ir = ir + 1
       CALL MPI_IRECV (buffer, SIZE(buffer), MPI_REAL4, p_source, p_tag, &
            p_comm, ireq(ir), ierr)
    END IF

#ifdef DEBUG
    IF (ierr /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', myid, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', ierr
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_real_1d_sp

  SUBROUTINE p_irecv_real_2d (buffer, p_source, p_tag, p_count, comm)

    REAL(nr), INTENT(out) :: buffer(:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

#ifndef NOMPI
    INTEGER :: p_comm

    IF (PRESENT(comm)) THEN
       p_comm = comm
    ELSE
       p_comm = icom
    ENDIF

    IF (PRESENT(p_count)) THEN
       ir = ir + 1
       CALL MPI_IRECV (buffer, p_count, MPI_REAL8, p_source, p_tag, &
            p_comm, ireq(ir), ierr)
    ELSE
       ir = ir + 1
       CALL MPI_IRECV (buffer, SIZE(buffer), MPI_REAL8, p_source, p_tag, &
            p_comm, ireq(ir), ierr)
    END IF

#ifdef DEBUG
    IF (ierr /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a,i4,a,i4,a,i6,a)') ' MPI_IRECV on ', myid, &
            ' from ', p_source, ' for tag ', p_tag, ' failed.'
       WRITE (nerr,'(a,i4)') ' Error = ', ierr
       CALL p_abort
    END IF
#endif
#endif

  END SUBROUTINE p_irecv_real_2d

  SUBROUTINE p_waitall
#ifndef NOMPI
    if(ir/=0) then
      call MPI_WAITALL(ir,ireq,ista,ierr)
    end if
#endif
  END SUBROUTINE p_waitall

  SUBROUTINE p_stop

    ! finish MPI and clean up all PEs
#ifndef NOMPI
    ! to prevent abort due to unfinished communication
    CALL p_barrier

    CALL MPI_FINALIZE (ierr)

    IF (ierr /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_FINALIZE failed.'
       STOP 1
    END IF
!    DEALLOCATE(p_request)
#endif

  END SUBROUTINE p_stop

  SUBROUTINE p_abort

    CALL MPI_ABORT (MPI_COMM_WORLD, 1, ierr)

    IF (ierr /= MPI_SUCCESS) THEN
       WRITE (nerr,'(a)') ' MPI_ABORT failed.'
       WRITE (nerr,'(a,i4)') ' Error =  ', ierr
       STOP 1
    END IF

    deallocate(ista,ireq)

  END SUBROUTINE p_abort

END MODULE mo_mpi