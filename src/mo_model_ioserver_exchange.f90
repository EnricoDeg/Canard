MODULE mo_model_ioserver_exchange
  USE mo_kind, ONLY : ni, nr, ieee32
  USE mo_mpi,  ONLY : p_barrier, p_get_intercomm, p_model2io, p_recv, &
                    & p_send, p_get_global_comm
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: meeting_point, send_cmd, get_cmd, recv_model2io, send_model2io

  INTERFACE send_model2io
    MODULE PROCEDURE send_model2io_1d_ieee32
  END INTERFACE send_model2io

  INTERFACE recv_model2io
    MODULE PROCEDURE recv_model2io_1d_ieee32
  END INTERFACE recv_model2io

  CONTAINS

  SUBROUTINE meeting_point
    INTEGER(kind=ni) :: iintercomm
    iintercomm = p_get_intercomm()
    CALL p_barrier(iintercomm)
  END SUBROUTINE meeting_point

  SUBROUTINE send_cmd(i_cmd)
    INTEGER(kind=ni), INTENT(IN) :: i_cmd
    INTEGER(kind=ni) :: i_cmd_unused
    call p_model2io(model=i_cmd, server=i_cmd_unused, root=0, lmodel_role=.true.)
  END SUBROUTINE send_cmd

  INTEGER(kind=ni) FUNCTION get_cmd() RESULT(i_cmd)
    call p_model2io(model=i_cmd, server=i_cmd, root=0, lmodel_role=.false.)
  END FUNCTION get_cmd

  SUBROUTINE send_model2io_1d_ieee32(field, mp, ljs, lje)
    REAL(kind=ieee32), INTENT(in)      :: field(:)
    INTEGER(kind=ni), intent(in) :: mp
    INTEGER(kind=ni), intent(in) :: ljs, lje

    INTEGER(kind=ni) :: global_comm, lmpi, itag

    global_comm = p_get_global_comm()
    lmpi = lje - ljs + 1
    itag = 1
    call p_send(field(ljs:lje), mp, itag, p_count=lmpi, comm=global_comm)

  END SUBROUTINE send_model2io_1d_ieee32

  SUBROUTINE recv_model2io_1d_ieee32(field, mps, mpe, lis, lie)
    REAL(kind=ieee32), INTENT(inout)      :: field(:)
    INTEGER(kind=ni), INTENT(IN) :: mps
    INTEGER(kind=ni), INTENT(IN) :: mpe
    INTEGER(kind=ni), dimension(0:mpe-mps), intent(in) :: lis
    INTEGER(kind=ni), dimension(0:mpe-mps), intent(in) :: lie

    INTEGER(kind=ni) :: mp, lmpi, itag, global_comm

    global_comm = p_get_global_comm()

    do mp=mps,mpe
      lmpi = lie(mp-mps) - lis(mp-mps) + 1
      itag = 1
      call p_recv(field(lis(mp-mps):lie(mp-mps)), mp, itag, p_count=lmpi, comm=global_comm)
    end do

  END SUBROUTINE recv_model2io_1d_ieee32

END MODULE mo_model_ioserver_exchange