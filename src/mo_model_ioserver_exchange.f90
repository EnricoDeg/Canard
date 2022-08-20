MODULE mo_model_ioserver_exchange
  USE mo_kind, ONLY : ni
  USE mo_mpi,  ONLY : p_barrier, p_get_intercomm, p_model2io
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: meeting_point, send_cmd, get_cmd

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

END MODULE mo_model_ioserver_exchange