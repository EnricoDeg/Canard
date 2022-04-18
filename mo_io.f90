!*****
!***** INPUT/OUTPUT MODULE
!*****

MODULE mo_io
   use mainvar3d
   use mo_mpi, ONLY : mpro, myid, p_barrier, p_recv, p_send, &
                      p_null_req, p_irecv, p_waitall
   use problemcase
   use mo_numerics
   IMPLICIT NONE
   PUBLIC

   integer(kind=ni),   dimension(:),allocatable :: idsgnl,lsgnl
   integer(kind=int64),dimension(:),allocatable :: lhmb
   
   CONTAINS

   SUBROUTINE read_inputo

      open(9,file='inputo.dat',status='old')
      read(9,*) cinput,mbk
      read(9,*) cinput,nts
      read(9,*) cinput,nscrn,nsgnl
      read(9,*) cinput,ndata,ndatafl,ndataav
      read(9,*) cinput,nkrk
      read(9,*) cinput,nviscous
      read(9,*) cinput,nsmf
      read(9,*) cinput,nrestart
      read(9,*) cinput,nextrabc,nextgcic
      read(9,*) cinput,nnf(:)
      read(9,*) cinput,reoo,tempoo
      read(9,*) cinput,amach1,amach2,amach3
      read(9,*) cinput,wtemp
      read(9,*) cinput,cfl
      read(9,*) cinput,tmax,timf,tsam
      read(9,*) cinput,dto
      close(9)

   END SUBROUTINE read_inputo

   SUBROUTINE read_inputp

      open(9,file='inputp.dat',status='old')
      read(9,*) cinput,lxi0,lxi1,lxi2
      read(9,*) cinput,let0
      read(9,*) cinput,lze0
      read(9,*) cinput,nbody,nthick
      read(9,*) cinput,ngridv
      read(9,*) cinput,nbpc(:,1)
      read(9,*) cinput,nbpc(:,2)
      read(9,*) cinput,nbpc(:,3)
      read(9,*) cinput,smg,smgvr
      read(9,*) cinput,doml0,doml1,domh
      read(9,*) cinput,span
      read(9,*) cinput,wlew,wlea
      read(9,*) cinput,szth0,szth1,szco
      read(9,*) cinput,skew,spx
      close(9)

      lximb(:)=(/lxi0,lxi1,lxi2,lxi0,lxi1,lxi2/)
      letmb(:)=(/let0,let0,let0,let0,let0,let0/)
      lzemb(:)=(/lze0,lze0,lze0,lze0,lze0,lze0/)

      allocate(idsgnl(0:lze0),lsgnl(0:lze0))

   END SUBROUTINE read_inputp

   SUBROUTINE allocate_io_memory
      ll=3+5*(ndata+1)
      allocate(times(0:ndata),cfilet(-1:ndata),     &
               ctecplt(-1:ndata),varm(0:1,0:mpro),  &
               varmin(ll),varmax(ll),cthead(0:mbk), &
               czonet(0:mbk),lhmb(0:mbk))
   END SUBROUTINE allocate_io_memory

   SUBROUTINE output_init

      cfilet(-1)='grid'
      do n=0,ndata
         no(2)=n/100
         no(1)=mod(n,100)/10
         no(0)=mod(n,10)
         cno=achar(no+48)
         cfilet(n)='n'//cno(2)//cno(1)//cno(0)
      end do
      do n=-1,ndata
         ctecplt(n)='data/'//cfilet(n)//'.plt'
      end do
      do mm=0,mbk
         no(2)=mm/100
         no(1)=mod(mm,100)/10
         no(0)=mod(mm,10)
         cno=achar(no+48)
         czonet(mm)='z'//cno(2)//cno(1)//cno(0)
         cthead(mm)='data/'//czonet(mm)//'.plt'
      end do
      cgrid='misc/grid'//czonet(mb)//'.dat'
      crestart='misc/restart'//czonet(mb)//'.dat'

      no(4)=myid/10000
      no(3)=mod(myid,10000)/1000
      no(2)=mod(myid,1000)/100
      no(1)=mod(myid,100)/10
      no(0)=mod(myid,10)
      cno=achar(no+48)
      cnnode=cno(4)//cno(3)//cno(2)//cno(1)//cno(0)
      cdata='misc/data'//cnnode//'.dat'

      do mm=0,mbk
         lpos(mo(mm))=0
         do i=1,nbpc(mm,1)-1
            mp=mo(mm)+i
            lpos(mp)=lpos(mp-1)+lxim(mp-1)+1
         end do
         jp=nbpc(mm,1)
         do j=1,nbpc(mm,2)-1
            do i=0,nbpc(mm,1)-1
               mp=mo(mm)+j*jp+i
               lpos(mp)=lpos(mp-jp)+(lximb(mm)+1)*(letm(mp-jp)+1)
            end do
         end do
         kp=nbpc(mm,1)*nbpc(mm,2)
         do k=1,nbpc(mm,3)-1
            do j=0,nbpc(mm,2)-1
               do i=0,nbpc(mm,1)-1
                  mp=mo(mm)+k*kp+j*jp+i
                  lpos(mp)=lpos(mp-kp)+(lximb(mm)+1)*(letmb(mm)+1)*(lzem(mp-kp)+1)
               end do
            end do
         end do
      end do

   END SUBROUTINE output_init

   SUBROUTINE read_restart_file

      open(9,file=crestart,access='direct',form='unformatted',recl=5*nrecd,status='old')
      read(9,rec=1) cha(:)
      read(9,rec=2) dha(:)
      n=cha(1)
      ndt=cha(2)
      dt=cha(3)
      dts=cha(4)
      dte=cha(5)
      timo=dha(1)
      lp=lpos(myid)+2
      do k=0,lze
         do j=0,let
            lq=lp+lio(j,k)
            do i=0,lxi
               l=indx3(i,j,k,1)
               read(9,rec=lq+i+1) qa(l,:)
            end do
         end do
      end do
      close(9)

   END SUBROUTINE read_restart_file

   SUBROUTINE write_restart_file

      if(myid==mo(mb)) then
         open(9,file=crestart,access='direct',form='unformatted',recl=5*nrecd,status='replace')
         cha(:)=(/real(n,kind=nr),real(ndt,kind=nr),dt,dts,dte/)
         dha(:)=(/timo,zero,zero,zero,zero/)
         write(9,rec=1) cha(:)
         write(9,rec=2) dha(:)
         close(9)
      end if
      call p_barrier
      open(9,file=crestart,access='direct',form='unformatted',recl=5*nrecd,status='old')
      lp=lpos(myid)+2
      do k=0,lze
         do j=0,let
            lq=lp+lio(j,k)
            do i=0,lxi
               l=indx3(i,j,k,1)
               write(9,rec=lq+i+1) qa(l,:)
            end do
         end do
      end do
      close(9)

   END SUBROUTINE write_restart_file

   SUBROUTINE write_output_file

      lje=-1
      do n=-1,ndata
         mq=3+2*min(n+1,1)
         llmb=mq*ltomb-1
         allocate(vara(0:llmb),varb(0:llmb))
         ljs=lje+1
         lje=ljs+mq*(lmx+1)-1
         if(myid==mo(mb)) then !===========================
            mps=mo(mb)
            mpe=mps+nbpc(mb,1)*nbpc(mb,2)*nbpc(mb,3)-1
            lis=0
            lie=mq*(lmx+1)-1
            vara(lis:lie)=vart(ljs:lje)
            do mp=mps+1,mpe
               lis=lie+1
               lie=lis+mq*(lxim(mp)+1)*(letm(mp)+1)*(lzem(mp)+1)-1
               lmpi = lie-lis+1
               itag=1
               call p_recv(vara(lis:lie), mp, itag, lmpi)
            end do
            lis=0
            do mp=mps,mpe
               do m=1,mq
                  do k=0,lzem(mp)
                     do j=0,letm(mp)
                        ljs=lpos(mp)+(m-1)*ltomb+k*(leto+1)*(lxio+1)+j*(lxio+1)
                        varb(ljs:ljs+lxim(mp))=vara(lis:lis+lxim(mp))
                        lis=lis+lxim(mp)+1
                     end do
                  end do
               end do
            end do
            open(9,file=cthead(mb),access='stream',form='unformatted',status='replace')
            call techead(9,n,mb,lh)
            deallocate(vara)
            allocate(vara(0:lh+llmb))
            read(9,pos=1) vara(0:lh-1)
            close(9,status='delete')
            lhmb(mb)=lh+llmb+1
            vara(lh:lh+llmb)=varb(:)
            if(mb==0) then !---------------------
               do mm=1,mbk
                  itag=2
                  call p_recv(lhmb(mm), mo(mm), itag)
               end do
               llmo=sum(lhmb(:))-1
               deallocate(varb)
               allocate(varb(0:llmo))
               lis=0
               lie=lhmb(mb)-1
               varb(lis:lie)=vara(:)
               do mm=1,mbk
                  lis=lie+1
                  lie=lis+lhmb(mm)-1
                  lmpi = lie-lis+1
                  itag=3
                  call p_recv(varb(lis:lie), mo(mm), itag, lmpi)
               end do
               open(0,file=ctecplt(n),status='unknown')
               close(0,status='delete') ! 'replace' not suitable as 'recl' may vary
               open(0,file=ctecplt(n),access='direct',form='unformatted',recl=nrecs*(llmo+1),status='new')
               write(0,rec=1) varb(:)
               close(0)
            else !-------------------
               itag=2
               call p_send(lhmb(mb), mo(0), itag)
               lmpi = lhmb(mb)
               itag=3
               call p_send(vara(:), mo(0), itag, lmpi)
            end if !-----------
         else !===========================
            lmpi = lje-ljs+1
            itag=1
            call p_send(vart(ljs:lje), mo(mb), itag, lmpi)
         end if !=========================
         deallocate(vara,varb)
      end do
 
   END SUBROUTINE write_output_file

!===== SUBROUTINE FOR GENERATING TECPLOT DATA FILE

   subroutine techead(nf,n,mb,lh)

      integer(kind=ni),intent(in) :: nf,n,mb
      integer(kind=ni),intent(inout) :: lh

      lh=0
      if(mb==0) then
         write(nf,pos=4*lh+1) '#!TDV112'
         lh=lh+2
         write(nf,pos=4*lh+1) 1
         lh=lh+1 ! Header Section
         write(nf,pos=4*lh+1) int(min(n+2,2),kind=int32)
         lh=lh+1 ! File Type (0 = Full / 1 = Grid / 2 = Solution)
         cinput=cfilet(n)
         call strio(nf,lh,cinput) ! File Title
         write(nf,pos=4*lh+1) int(mq,kind=int32)
         lh=lh+1 ! Number of Variables
         if(n==-1) then
            cinput='x'
            call strio(nf,lh,cinput)
             cinput='y'
            call strio(nf,lh,cinput)
            cinput='z'
            call strio(nf,lh,cinput)
         else
            cinput='rho'
            call strio(nf,lh,cinput)
            cinput='u'
            call strio(nf,lh,cinput)
            cinput='v'
            call strio(nf,lh,cinput)
            cinput='w'
            call strio(nf,lh,cinput)
            cinput='p'
            call strio(nf,lh,cinput)
         end if
         do mm=0,mbk
            write(nf,pos=4*lh+1) 299.0
            lh=lh+1 ! Zone Marker
            cinput=czonet(mm)
            call strio(nf,lh,cinput) ! Zone Name
            write(nf,pos=4*lh+1) -1
            lh=lh+1 ! Parent Zone
            write(nf,pos=4*lh+1) int(mm+1,kind=int32)
            lh=lh+1 ! Strand ID
            write(nf,pos=4*lh+1) real(times(max(n,0)),kind=ieee64)
            lh=lh+2 ! Solution Time (Double)
            write(nf,pos=4*lh+1) -1
            lh=lh+1 ! (Not used. Set to -1.)
            write(nf,pos=4*lh+1) 0
            lh=lh+1 ! Zone Type
            write(nf,pos=4*lh+1) 0
            lh=lh+1 ! Specify Var Location
            write(nf,pos=4*lh+1) 0
            lh=lh+1 ! Raw Local 1-to-1 Face Neighbours Suppliled
            write(nf,pos=4*lh+1) 0
            lh=lh+1 ! Number of Miscellaneous Face Neighbour Connections
            write(nf,pos=4*lh+1) int(lximb(mm)+1,kind=int32)
            lh=lh+1 ! IMax
            write(nf,pos=4*lh+1) int(letmb(mm)+1,kind=int32)
            lh=lh+1 ! JMax
            write(nf,pos=4*lh+1) int(lzemb(mm)+1,kind=int32)
            lh=lh+1 ! KMax
            write(nf,pos=4*lh+1) 0
            lh=lh+1 ! No Auxillary Data Pairs
         end do
         write(nf,pos=4*lh+1) 357.0
         lh=lh+1 ! End of Header Marker
      end if
      write(nf,pos=4*lh+1) 299.0
      lh=lh+1 ! Zone Marker
      do m=1,mq
         write(nf,pos=4*lh+1) 1
         lh=lh+1 ! 1 = Float / 2 = Double
      end do
      write(nf,pos=4*lh+1) 0
      lh=lh+1 ! No Passive Variables
      write(nf,pos=4*lh+1) 0
      lh=lh+1 ! No Variable Sharing
      write(nf,pos=4*lh+1) -1
      lh=lh+1 ! Zero Based Zone Number to Share
      do m=1,mq
         nn=max(3+5*n,0)+m
         write(nf,pos=4*lh+1) real(varmin(nn),kind=ieee64)
         lh=lh+2 ! Minimum Value (Double) of Variables
         write(nf,pos=4*lh+1) real(varmax(nn),kind=ieee64)
         lh=lh+2 ! Maximum Value (Double) of Variables
      end do

   end subroutine techead

!===== SUBROUTINE FOR CHARACTER STRING CONVERSION

   subroutine strio(nfile,lh,cinput)

      integer(kind=ni),intent(in) :: nfile
      integer(kind=ni),intent(inout) :: lh
      character(16),intent(in) :: cinput
      integer(kind=ni) :: ll

      do ll=1,len_trim(cinput)
         write(nfile,pos=4*lh+1) ichar(cinput(ll:ll))
         lh=lh+1
      end do
      write(nfile,pos=4*lh+1) 0
      lh=lh+1

   end subroutine strio

!===== SUBROUTINE FOR FINDING VARIABLE MIN/MAX VALUES FOR TECPLOT DATA FILE

   subroutine vminmax(nn)

      integer(kind=ni),intent(in) :: nn

      varmin(nn)=minval(varr)
      varmax(nn)=maxval(varr)
      varm(:,myid)=(/varmin(nn),varmax(nn)/)

      call p_null_req
      itag=nn
      if(myid==mo(mb)) then
         mps=mo(mb)
         mpe=mps+nbpc(mb,1)*nbpc(mb,2)*nbpc(mb,3)-1
         do mp=mps+1,mpe
            call p_irecv(varm(:,mp), mp, itag, 2)
         end do
         call p_waitall
         varmin(nn)=minval(varm(0,mps:mpe))
         varmax(nn)=maxval(varm(1,mps:mpe))
      else
         call p_send(varm(:,myid), mo(mb), itag, 2)
      end if

   end subroutine vminmax

END MODULE mo_io