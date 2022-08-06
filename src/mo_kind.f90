MODULE mo_kind

   implicit none
   public
   
   integer,parameter :: int32=selected_int_kind(9),int64=selected_int_kind(18)
   integer,parameter :: ieee32=selected_real_kind(6,37),ieee64=selected_real_kind(15,307)
   integer,parameter :: ni=int32,nr=ieee64,nsp=ieee32, nli=int64

END MODULE mo_kind