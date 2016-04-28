      module constants
      double precision hbc,fmsp,fmsn,piln
      double precision rdmavg,fmbyh2,recmass
      double precision rmup,rmun,rlambdaVD,rlambdan
      double precision rhos,rmus,rlambdaEs
      parameter (hbc=197.327053d0)
      parameter (fmsp=938.27231d0,fmsn=939.56563d0)
      parameter (piln=3.14159265358979312d0)
      parameter (rmup=2.79284736d0,rmun=-1.9130427d0)
*      parameter (rmup=2.79d0,rmun=-1.91d0)
      parameter (rlambdaVD=4.97d0,rlambdan=5.6d0)
      parameter (rhos=-2.d0,rmus=0.23d0,rlambdaEs=rlambdan)
      parameter (sinthetaWsq=0.2315)
*      parameter (sinthetaWsq=0.23)
      parameter (rdmavg=2.d0*(fmsp*fmsn)/(fmsp+fmsn))
      parameter (fmbyh2=rdmavg/hbc**2)
      parameter (recmass=hbc*hbc/(4.d0*rdmavg*rdmavg))
      real(kind=kind(0.0d0)) epi,enu,glp,gln,gsp,gsn
      data epi,enu,glp,gln,gsp,gsn
     + /1.d0,0.d0,1.d0,0.d0,5.586d0,-3.826d0/
      end module constants 

      module paramdef
      integer :: nbit,nbit1
      parameter (nbit=64,nbit1=nbit-1)
      integer :: nblock=1
      integer,allocatable :: global_ist1(:,:),global_ist1_sp(:,:),
     $     global_ist1_sp_map(:)
      integer :: global_ist1_dim,global_ist1_sp_dim
      integer :: N_cluster_max,N_cluster_min,N_sp_min,N_proj_max,
     $     num_neutr_proj,num_prot_proj
      integer :: N_RGM_max
      integer :: jrelm=8
      type jsp
      integer,allocatable :: j(:)
      end type jsp
      type lspjsp
      type(jsp),allocatable :: l(:)
      end type lspjsp
      type(lspjsp),allocatable :: nlj_st(:)
      type mjsp
      integer,allocatable :: mj_sp(:)
      end type mjsp
      type(mjsp),allocatable :: nljmmt_st(:,:)
      end module paramdef

      module gamalog
      integer :: nrmax,lrelm,maxgam,lmax
      real(kind=kind(0.0d0)),allocatable,dimension(:,:) :: dsq
      real(kind=kind(0.0d0)),allocatable,dimension(:) :: gamal
      end module gamalog 

      module initst
      integer(2),pointer :: iloci(:,:),mconfi(:)
      integer(4),pointer :: nconfi,nprtni(:,:),nneuti(:,:),iendconfi(:)
      integer(4),pointer:: jt2i(:),it2i(:)
      integer(4),pointer :: n_spi(:),l_spi(:),
     +     j2_spi(:),m2_spi(:),mt2_spi(:)
      real(kind=kind(0.0d0)),pointer::bmpi(:,:)
      integer(4),pointer:: ibasisi(:,:,:)
      integer(4),pointer :: nsdi,mxnwdi,nhomi,nhom12i,naspsi,mxspsi,
     +     nhom123i,nhom1234i
      integer(4),pointer :: mjtotali,mttotali,iparityi,nshlli,nhwi
      integer(4),pointer :: nucleonsi,nprotonsi,nneutrnsi
      real(kind=kind(0.0d0)),pointer :: Jxi(:),Txi(:),eneri(:)
      real(kind=kind(0.0d0)),pointer :: hboi
      integer(4),pointer :: ibascompi(:)
      integer(4),pointer:: I_state_i(:,:),nsd_p_i,nsd_n_i,
     +        occ_p_i(:,:),occ_n_i(:,:),ibas_p_i(:,:),ibas_n_i(:,:)
      integer :: num_of_in_i
      end module initst

      module finast
      integer(2),pointer :: ilocf(:,:),mconff(:)
      integer(4),pointer :: nconff,nprtnf(:,:),nneutf(:,:),iendconff(:)
      integer(4),pointer:: jt2f(:),it2f(:)
      integer(4),pointer :: n_spf(:),l_spf(:),
     +     j2_spf(:),m2_spf(:),mt2_spf(:)
      real(kind=kind(0.0d0)),pointer::bmpf(:,:)
      integer(4),pointer :: ibasisf(:,:,:)
      integer(4),pointer :: nsdf,mxnwdf,nhomf,nhom12f,naspsf,mxspsf,
     +     nhom123f,nhom1234f
      integer(4),pointer :: mjtotalf,mttotalf,iparityf,nshllf,nhwf
      integer(4),pointer :: nucleonsf,nprotonsf,nneutrnsf
      real(kind=kind(0.0d0)),pointer :: Jxf(:),Txf(:),enerf(:)
      real(kind=kind(0.0d0)),pointer :: hbof
      integer(4),pointer :: ibascompf(:)
      integer(4),pointer:: I_state_f(:,:),nsd_p_f,nsd_n_f,
     +        occ_p_f(:,:),occ_n_f(:,:),ibas_p_f(:,:),ibas_n_f(:,:)
      integer :: num_of_in_f
      end module finast

      module cleb_coef
      private
      public cleb_3j,cleb_init,cleb_destroy,num_of_3j_types_init,
     +     num_of_3j_types_destroy
      integer :: num_of_3j_types
      integer,allocatable :: mshift(:)

      type cleb_value
      real(kind(0.d0)),pointer :: val => null()
      end type cleb_value

      type cleb_j3
      integer :: j3_min,j3_max
      type(cleb_value),allocatable :: j3(:)
      end type cleb_j3

      type cleb_m2
      integer :: m2_min,m2_max
      type(cleb_j3),allocatable :: m2(:)
      end type cleb_m2

      type cleb_j2
      integer :: j2_min,j2_max
      type(cleb_m2),allocatable :: j2(:)
      end type cleb_j2

      type cleb_m1
      integer :: m1_min,m1_max
      type(cleb_j2),allocatable :: m1(:)
      end type cleb_m1

      type cleb_j1
      integer :: j1_min,j1_max
      type(cleb_m1),allocatable :: j1(:)
      end type cleb_j1

      type(cleb_j1),allocatable :: cl_3j(:)
      
      contains

      subroutine num_of_3j_types_init(num)
      implicit none
      integer,intent(IN) :: num
      num_of_3j_types=num
      allocate(cl_3j(num_of_3j_types))
      allocate(mshift(num_of_3j_types))
      end subroutine num_of_3j_types_init

      subroutine cleb_init(num,j1_min,j1_max,m1_min,m1_max,
     +     j2_min,j2_max,m2_min,m2_max,j3_min,j3_max)
      implicit none
      integer,intent(IN) :: num,j1_min,j1_max,m1_min,m1_max,
     +     j2_min,j2_max,m2_min,m2_max,j3_min,j3_max

      integer :: j1,m1,j2,m2,j3,m1shift,m2shift
      if (num>num_of_3j_types) then
         print *,'*** error: numn_of_cleb_types,num=',
     +        num_of_3j_types,num
         stop
      endif
      if (mod(j1_min+j2_min+j3_min,2)/=0.or.
     +     mod(j1_max+j2_max+j3_max,2)/=0.or.
     +     mod(j1_min+m1_min,2)/=0.or.
     +     mod(j1_max+m1_max,2)/=0.or.
     +     mod(j1_min+m1_max,2)/=0.or.
     +     mod(j1_max+m1_min,2)/=0.or.
     +     mod(j2_min+m2_max,2)/=0.or.
     +     mod(j2_max+m2_min,2)/=0.or.
     +     mod(j2_min+m2_min,2)/=0.or.
     +     mod(j2_max+m2_max,2)/=0) then
         print *,'*** error: j_min,j_max incompatible'
         stop
      endif
      mshift(num)=max(abs(m1_min),abs(m2_min))
      cl_3j(num)%j1_min=j1_min
      cl_3j(num)%j1_max=j1_max

c      print *,' j1_min,j1_max',j1_min,j1_max

      allocate(cl_3j(num)%j1(j1_min/2:j1_max/2))

c      print *,' j1 alloc'

      do j1=j1_min,j1_max,2
         cl_3j(num)%j1(j1/2)%m1_min=max(m1_min,-j1)
         cl_3j(num)%j1(j1/2)%m1_max=min(m1_max,j1)

c         print *,' j1,m1_min,m1_max=',j1,max(m1_min,-j1),min(m1_max,j1)

         allocate(cl_3j(num)%j1(j1/2)%m1((mshift(num)+max(m1_min,-j1))/2
     +        :(mshift(num)+min(m1_max,j1))/2))

c         print *,' m1 alloc'

         do m1=max(m1_min,-j1),min(m1_max,j1),2
            m1shift=(mshift(num)+m1)/2
            cl_3j(num)%j1(j1/2)%m1(m1shift)%j2_min=j2_min
            cl_3j(num)%j1(j1/2)%m1(m1shift)%j2_max=j2_max

c            print *,' m1,j2_min,j2_max=',m1,j2_min,j2_max

            allocate(cl_3j(num)%j1(j1/2)%m1(m1shift)
     $           %j2(j2_min/2:j2_max/2))

c            print *,' j2 alloc'

            do j2=j2_min,j2_max,2
               cl_3j(num)%j1(j1/2)%m1(m1shift)%j2(j2/2)%m2_min=
     +              max(m2_min,-j2)
               cl_3j(num)%j1(j1/2)%m1(m1shift)%j2(j2/2)%m2_max=
     +              min(m2_max,j2)

c               print *,' j2,m2_min,m2_max=',j2,max(m2_min,-j2),
c     $              min(m2_max,j2)

               allocate(cl_3j(num)%j1(j1/2)%m1(m1shift)%j2(j2/2)%m2(
     +              (mshift(num)+max(m2_min,-j2))/2:
     $              (mshift(num)+min(m2_max,j2))/2))

c               print *,' m2 alloc'

               do m2=max(m2_min,-j2),min(m2_max,j2),2
                  m2shift=(mshift(num)+m2)/2
                  cl_3j(num)%j1(j1/2)%m1(m1shift)%j2(j2/2)%m2(m2shift)
     $                 %j3_min=max(j3_min,abs(j1-j2),abs(m1+m2))
                  cl_3j(num)%j1(j1/2)%m1(m1shift)%j2(j2/2)%m2(m2shift)
     $                 %j3_max=min(j3_max,j1+j2)

c                  print *,' m2,j3_min,j3_max=',m2,
c     $                 max(j3_min,abs(j1-j2),abs(m1+m2)),
c     $                 min(j3_max,j1+j2)

                  allocate(cl_3j(num)%j1(j1/2)%m1(m1shift)%j2(j2/2)
     +                 %m2(m2shift)%j3(
     $                 max(j3_min,abs(j1-j2),abs(m1+m2))/2:
     +                 min(j3_max,j1+j2)/2))

c                  print *,' j3 alloc'

               end do
            end do
         end do
      end do
      end subroutine cleb_init

      subroutine cleb_destroy(num)
      implicit none
      integer,intent(IN) :: num
      if (num>num_of_3j_types) then
         print *,'*** error in destroy: numn_of_cleb_types,num=',
     +        num_of_3j_types,num
         stop
      endif
      if (allocated(cl_3j(num)%j1)) deallocate(cl_3j(num)%j1)
      end subroutine cleb_destroy

      subroutine num_of_3j_types_destroy
      implicit none
      if (allocated(cl_3j)) deallocate(cl_3j)
      if (allocated(mshift)) deallocate(mshift)
      end subroutine num_of_3j_types_destroy

      real(kind(0.d0)) function cleb_3j(num,j1,m1,j2,m2,j3)
      implicit none
      integer,intent(IN) :: num,j1,m1,j2,m2,j3
      real(kind(0.d0)) :: clebd
      integer :: m1shift,m2shift
      m1shift=(mshift(num)+m1)/2
      m2shift=(mshift(num)+m2)/2
      if (num>num_of_3j_types.or.j1<cl_3j(num)%j1_min
     +     .or.j1>cl_3j(num)%j1_max
     +     .or.m1<cl_3j(num)%j1(j1/2)%m1_min
     +     .or.m1>cl_3j(num)%j1(j1/2)%m1_max
     +     .or.j2<cl_3j(num)%j1(j1/2)%m1(m1shift)%j2_min
     +     .or.j2>cl_3j(num)%j1(j1/2)%m1(m1shift)%j2_max
     +     .or.m2<cl_3j(num)%j1(j1/2)%m1(m1shift)%j2(j2/2)%m2_min
     +     .or.m2>cl_3j(num)%j1(j1/2)%m1(m1shift)%j2(j2/2)%m2_max
     +     .or.j3<cl_3j(num)%j1(j1/2)%m1(m1shift)%j2(j2/2)%m2(m2shift)
     $     %j3_min
     +     .or.j3>cl_3j(num)%j1(j1/2)%m1(m1shift)%j2(j2/2)%m2(m2shift)
     $     %j3_max
     +     ) then
         cleb_3j=clebd(j1,m1,j2,m2,j3,m1+m2)
      else
         if (associated(cl_3j(num)%j1(j1/2)%m1(m1shift)%j2(j2/2)
     +        %m2(m2shift)%j3(j3/2)%val)) then
            cleb_3j=cl_3j(num)%j1(j1/2)%m1(m1shift)%j2(j2/2)
     +           %m2(m2shift)%j3(j3/2)%val
         else
            cleb_3j=clebd(j1,m1,j2,m2,j3,m1+m2)            
            allocate(cl_3j(num)%j1(j1/2)%m1(m1shift)%j2(j2/2)
     +           %m2(m2shift)%j3(j3/2)%val)
            cl_3j(num)%j1(j1/2)%m1(m1shift)%j2(j2/2)
     +           %m2(m2shift)%j3(j3/2)%val=cleb_3j
         endif
      endif
      end function cleb_3j

      end module cleb_coef

      module HO_braket
      private
      public osc_br_unit_mass_ratio,HO_braket_init,HO_braket_destroy
      integer :: N_max_glob,Nr_max_glob

c      type HO_value
c      real(kind(0.d0)),pointer :: val => null()
c      end type HO_value

      type HO_lambda
      integer :: lambda_min,lambda_max
c      type(HO_value),allocatable :: lambda(:)
      real(kind(0.d0)),allocatable :: lambda(:)
      end type HO_lambda

      type HO_n2
      integer :: n2_min,n2_max
      type(HO_lambda),allocatable :: n2(:)
      end type HO_n2

      type HO_l1
      integer :: l1_min,l1_max
      type(HO_n2),allocatable :: l1(:)
      end type HO_l1

      type HO_n1
      integer :: n1_min,n1_max
      type(HO_l1),allocatable :: n1(:)
      end type HO_n1

      type HO_Lc
      integer :: Lc_min,Lc_max
      type(HO_n1),allocatable :: Lc(:)
      end type HO_Lc

      type HO_Nc
      integer :: Nc_min,Nc_max
      type(HO_Lc),allocatable :: Nc(:)
      end type HO_Nc

      type HO_lr
      integer :: lr_min,lr_max
      type(HO_Nc),allocatable :: lr(:)
      end type HO_lr

      type HO_nr
      integer :: nr_min,nr_max
      type(HO_lr),allocatable :: nr(:)
      end type HO_nr

      type(HO_nr),save :: oscbr

      contains

      subroutine HO_braket_init(Nr_max,N_max)
      implicit none
      integer,intent(IN) :: Nr_max,N_max
      integer :: nr,lr,Nc,Lc
      N_max_glob=N_max
      Nr_max_glob=Nr_max
      oscbr%nr_min=0
      oscbr%nr_max=Nr_max/2
      allocate(oscbr%nr(0:Nr_max/2))
      do nr=0,Nr_max/2
         oscbr%nr(nr)%lr_min=0
         oscbr%nr(nr)%lr_max=Nr_max-2*nr
         allocate(oscbr%nr(nr)%lr(0:Nr_max-2*nr))
         do lr=0,Nr_max-2*nr
            oscbr%nr(nr)%lr(lr)%Nc_min=0
            oscbr%nr(nr)%lr(lr)%Nc_max=(N_max-2*nr-lr)/2
            allocate(oscbr%nr(nr)%lr(lr)%Nc(0:(N_max-2*nr-lr)/2))
            do Nc=0,(N_max-2*nr-lr)/2
               oscbr%nr(nr)%lr(lr)%Nc(Nc)%Lc_min=0
               oscbr%nr(nr)%lr(lr)%Nc(Nc)%Lc_max=N_max-2*nr-lr-2*Nc
               allocate(oscbr%nr(nr)%lr(lr)%Nc(Nc)%Lc(0:
     $              N_max-2*nr-lr-2*Nc))
               do Lc=0,N_max-2*nr-lr-2*Nc
                  oscbr%nr(nr)%lr(lr)%Nc(Nc)%Lc(Lc)%n1_min=0
                  oscbr%nr(nr)%lr(lr)%Nc(Nc)%Lc(Lc)%n1_max=
     +                 min(N_max/2,(2*nr+lr+2*Nc+Lc)/2)
                  allocate(oscbr%nr(nr)%lr(lr)%Nc(Nc)%Lc(Lc)%n1(0:
     $                 min(N_max/2,(2*nr+lr+2*Nc+Lc)/2)))
                  call HO_braket_init_int(Nr_max,N_max,nr,lr,Nc,Lc)
               end do
            end do
         end do
      end do
      end subroutine HO_braket_init

      subroutine HO_braket_init_int(Nr_max,N_max,nr,lr,Nc,Lc)
      implicit none
      integer,intent(IN) :: Nr_max,N_max,nr,lr,Nc,Lc
      integer :: n1,l1,n2,l2,lambda
      do n1=0,min(N_max/2,(2*nr+lr+2*Nc+Lc)/2)
         oscbr%nr(nr)%lr(lr)%Nc(Nc)%Lc(Lc)%n1(n1)%l1_min=0
         oscbr%nr(nr)%lr(lr)%Nc(Nc)%Lc(Lc)%n1(n1)%l1_max=
     +        min(N_max-2*n1,2*nr+lr+2*Nc+Lc-2*n1)
         allocate(oscbr%nr(nr)%lr(lr)%Nc(Nc)%Lc(Lc)%n1(n1)%
     +        l1(0:min(N_max-2*n1,2*nr+lr+2*Nc+Lc-2*n1)))
         do l1=0,min(N_max-2*n1,2*nr+lr+2*Nc+Lc-2*n1)
            oscbr%nr(nr)%lr(lr)%Nc(Nc)%Lc(Lc)%n1(n1)%l1(l1)
     +           %n2_min=0
            oscbr%nr(nr)%lr(lr)%Nc(Nc)%Lc(Lc)%n1(n1)%l1(l1)
     +           %n2_max=min(N_max/2,
     +           (2*nr+lr+2*Nc+Lc-2*n1-l1)/2)
            allocate(oscbr%nr(nr)%lr(lr)%Nc(Nc)%Lc(Lc)
     +           %n1(n1)%l1(l1)%n2(0:min(N_max/2,
     +           (2*nr+lr+2*Nc+Lc-2*n1-l1)/2)))
            do n2=0,min(N_max/2,(2*nr+lr+2*Nc+Lc-2*n1-l1)/2)
               l2=2*nr+lr+2*Nc+Lc-2*n1-l1-2*n2
               if (l2<0) exit
               oscbr%nr(nr)%lr(lr)%Nc(Nc)%Lc(Lc)%n1(n1)
     +              %l1(l1)%n2(n2)%lambda_min=
     +              max(abs(l1-l2),abs(lr-Lc))
               oscbr%nr(nr)%lr(lr)%Nc(Nc)%Lc(Lc)%n1(n1)
     +              %l1(l1)%n2(n2)%lambda_max=
     +              min(l1+l2,lr+Lc)
               allocate(oscbr%nr(nr)%lr(lr)%Nc(Nc)%Lc(Lc)
     +              %n1(n1)%l1(l1)%n2(n2)%lambda(
     $              max(abs(l1-l2),abs(lr-Lc)):
     $              min(l1+l2,lr+Lc)))
               oscbr%nr(nr)%lr(lr)%Nc(Nc)%Lc(Lc)
     +              %n1(n1)%l1(l1)%n2(n2)%lambda=-111111.d0
            end do
         end do
      end do
      end subroutine HO_braket_init_int

      subroutine HO_braket_destroy
      implicit none
      if (allocated(oscbr%nr)) deallocate(oscbr%nr)
      end subroutine HO_braket_destroy

      subroutine osc_br_unit_mass_ratio(nr,lr,Nc,Lc,n1,l1,n2,l2,
     +     lambda,osc)
c      use hob
      implicit none
      integer,intent(IN) :: nr,lr,Nc,Lc,n1,l1,n2,l2,
     +     lambda
      real(kind(0.d0)),intent(OUT) :: osc
      integer :: success
      real(kind(0.d0)) :: oscb
      if (2*nr+lr+2*Nc+Lc/=2*n1+l1+2*n2+l2.or.
     +     lambda<max(abs(l1-l2),abs(lr-Lc)).or.
     +     lambda>min(l1+l2,lr+Lc)) then
         osc=0.d0
         return
      endif
      if (2*nr+lr>Nr_max_glob.or.2*Nc+Lc>N_max_glob
     $     .or.2*nr+lr+2*Nc+Lc>N_max_glob
     $     .or.2*n1+l1>N_max_glob
     +     .or.2*n2+l2>N_max_glob) then
         call osclbr(nr,lr,Nc,Lc,n1,l1,n2,l2,lambda,
     +        1.d0,oscb)
            osc=oscb
cc         osc=TMB(2*nr+lr,lr,2*Nc+Lc,Lc,2*n1+l1,l1,2*n2+l2,l2,lambda,
cc     +        real(num,kind(0.d0))/real(den,kind(0.d0)))
      else
c         if (associated(oscbr%nr(nr)%lr(lr)%Nc(Nc)%Lc(Lc)
c     +        %n1(n1)%l1(l1)%n2(n2)
c     +        %lambda(lambda)%val)) then
c            osc=oscbr%nr(nr)%lr(lr)%Nc(Nc)%Lc(Lc)
c     +           %n1(n1)%l1(l1)%n2(n2)
c     +           %lambda(lambda)%val
         if (oscbr%nr(nr)%lr(lr)%Nc(Nc)%Lc(Lc)
     +        %n1(n1)%l1(l1)%n2(n2)
     +        %lambda(lambda)/=-111111.d0) then
            osc=oscbr%nr(nr)%lr(lr)%Nc(Nc)%Lc(Lc)
     +           %n1(n1)%l1(l1)%n2(n2)
     +           %lambda(lambda)
         else
            call osclbr(nr,lr,Nc,Lc,n1,l1,n2,l2,lambda,
     +           1.d0,oscb)
            osc=oscb
cc            osc=TMB(2*nr+lr,lr,2*Nc+Lc,Lc,2*n1+l1,l1,2*n2+l2,l2,lambda,
cc     +           real(num,kind(0.d0))/real(den,kind(0.d0)))
c            allocate(oscbr%nr(nr)%lr(lr)%Nc(Nc)%Lc(Lc)
c     +        %n1(n1)%l1(l1)%n2(n2)
c     +        %lambda(lambda)%val,stat=success)
cc            if (success==0)
c            oscbr%nr(nr)%lr(lr)%Nc(Nc)%Lc(Lc)
c     +           %n1(n1)%l1(l1)%n2(n2)
c     +           %lambda(lambda)%val=osc
C$OMP CRITICAL            
            oscbr%nr(nr)%lr(lr)%Nc(Nc)%Lc(Lc)
     +           %n1(n1)%l1(l1)%n2(n2)
     +           %lambda(lambda)=osc
C$OMP END CRITICAL            
         endif
      endif
      end subroutine osc_br_unit_mass_ratio

      end module HO_braket

      module nodeinfo
      integer :: nproc,iproc,icomm
      end module nodeinfo

      module occ
      type occ_bas
      integer,pointer :: nucl,dim_st,prot,neut,mxsp,nwd
      integer(2),pointer :: occ(:,:)
      end type occ_bas
      type(occ_bas) :: ty
      end module occ

      module hash_tables
      use occ
      private
      public construct_hash_table,get_state_index,get_state_index2

      integer :: dim,maxdisp,mxnw
      integer,allocatable :: linhash(:)
      integer,allocatable :: occ_tmp(:)
      integer(8),allocatable :: intbas_tmp(:)

      contains
      
      subroutine construct_hash_table
      use nodeinfo
      implicit none
      integer,allocatable :: temp_ind(:)
      integer :: i,key,dim_st,nuc,max_el,temp,wd,cr_1,ibit,ip,in,
     $     nbit
      integer(8),allocatable :: intbas(:)

      nbit=64
      dim_st=ty%dim_st
      nuc=ty%nucl
      mxnw=2*ty%nwd

c      if (dim_st<1 431 655 765) then
      if (dim_st<750 000 000) then
c         dim=dim_st*3/2
         dim=2*(2**(nint(log(real(dim_st,kind(0.d0)))/log(2.d0))))
      else
c         dim=dim_st
         dim=2147483646
      endif

      if (iproc==0) print *,' dim,hash table length=',dim_st,dim

      allocate(linhash(0:dim-1))
      linhash=-1
      maxdisp=0

      allocate(intbas(mxnw))
      allocate(occ_tmp(nuc))
      do i=1,dim_st

         if (i==1000000*(i/1000000)) then
            if (iproc==0) print *,' construct_hash_table: i=',i
         endif

         occ_tmp(:)=ty%occ(:,i)
c         nwd=(maxval(occ_tmp)-1)/nbit+1
c         if (nwd>mxnwd) mxnwd=nwd
         intbas=0
         do cr_1=1,ty%nucl
            wd=(occ_tmp(cr_1)-1)/nbit+1
            ibit=mod(occ_tmp(cr_1)-1,nbit)
            intbas(wd)=ibset(intbas(wd),ibit)
         end do
c         print *,' occ_tmp=',occ_tmp

         call get_key(mxnw,intbas,key,dim)

         if (key<0.or.key>dim-1) then
            print *,'*** error in construct_hash_table: key,dim=',
     +           key,dim
            stop
         endif
c         print *,' key=',key

         max_el=0
         do
            if (linhash(key)==-1) then
               linhash(key)=i
               exit
            endif
            key=mod(key+1,dim)
            max_el=max_el+1
         end do
         if (max_el>maxdisp) maxdisp=max_el
      end do
      deallocate(intbas)
      allocate(intbas_tmp(mxnw))
      if (iproc==0) print *,' maximal displacement on proc0=',maxdisp

      end subroutine construct_hash_table

      subroutine get_key(nwd,intbas,key,M)
      implicit none
      integer,intent(IN) :: nwd,M
      integer(8),intent(IN) :: intbas(nwd)
      integer,intent(OUT) :: key
      integer(8) :: intbas_tmp
      integer :: wd,ibit
      integer(8) :: isum,mul
      isum=intbas(1)
      mul=isum
c! ieor,ishft appears to work the best
      isum=ieor(isum,ishft(isum,13))
      isum=ieor(isum,ishft(isum,-17))
      isum=ieor(isum,ishft(isum,5))
      do wd=2,nwd
         isum=isum+intbas(wd)
         isum=ieor(isum,ishft(isum,13))
         isum=ieor(isum,ishft(isum,-17))
         isum=ieor(isum,ishft(isum,5))
         mul=mul*intbas(wd)
      end do
      isum=ieor(isum,mul)
      key=abs(mod(isum,int(M,kind=8)))
      end subroutine get_key

      subroutine get_state_index(nuc,occ,type,index)
c      subroutine get_state_index(nuc,nwd,intbas,index)
      implicit none
      integer,intent(IN) :: nuc,occ(nuc),type
c      integer(8),intent(IN) :: intbas(nwd)
      integer,intent(OUT) :: index
      integer :: key,i,j,wd,ibit,ip,in
      logical :: OK

      intbas_tmp=0
      do i=1,nuc
         wd=(occ(i)-1)/64+1
         ibit=mod(occ(i)-1,64)
         intbas_tmp(wd)=ibset(intbas_tmp(wd),ibit)
      end do

      call get_key(mxnw,intbas_tmp,key,dim)

      index=linhash(key)
      if (index==-1) return
      occ_tmp(:)=ty%occ(:,index)

      OK=.true.
      do i=1,nuc
         if (occ_tmp(i)/=occ(i)) then
            OK=.false.
            exit
         endif
      end do

      if (OK) then
         return
      endif

      do j=1,maxdisp
         key=mod(key+1,dim)
         index=linhash(key)
         if (index==-1) return
         occ_tmp(:)=ty%occ(:,index)

         OK=.true.
         do i=1,nuc
            if (occ_tmp(i)/=occ(i)) then
               OK=.false.
               goto 1000
            endif
         end do

         if (OK) then
            return
         endif

 1000    continue
      end do
      index=-1
      end subroutine get_state_index

      subroutine get_state_index2(nuc,nwd,intbas,index)
      implicit none
      integer,intent(IN) :: nuc,nwd
      integer(8),intent(IN) :: intbas(nwd)
      integer,intent(OUT) :: index
      integer :: key,i,j,wd,ibit !,ip,in
      logical :: OK

      call get_key(nwd,intbas,key,dim)

      index=linhash(key)
      if (index==-1) return
c      ip=ty%occ(1,index)
c      in=ty%occ(2,index)
c      occ_tmp(1:ty%protn)=ty%occp(:,ip)
c      occ_tmp(ty%protn+1:nuc)=ty%occn(:,in)
      occ_tmp(:)=ty%occ(:,index)

      intbas_tmp=0
      do i=1,nuc
         wd=(occ_tmp(i)-1)/64+1
         ibit=mod(occ_tmp(i)-1,64)
c     print *,' i,wd,ibit=',i,wd,ibit
         intbas_tmp(wd)=ibset(intbas_tmp(wd),ibit)
      end do

      OK=.true.
      do wd=1,nwd
         if (intbas(wd)/=intbas_tmp(wd)) then
            OK=.false.
            exit
         endif
      end do

      if (OK) then
         return
      endif

      do j=1,maxdisp
         key=mod(key+1,dim)
         index=linhash(key)
         if (index==-1) return
c         ip=ty%occ(1,index)
c         in=ty%occ(2,index)
c         occ_tmp(1:ty%protn)=ty%occp(:,ip)
c         occ_tmp(ty%protn+1:nuc)=ty%occn(:,in)
         occ_tmp(:)=ty%occ(:,index)

         intbas_tmp=0
         do i=1,nuc
            wd=(occ_tmp(i)-1)/64+1
            ibit=mod(occ_tmp(i)-1,64)
c     print *,' i,wd,ibit=',i,wd,ibit
            intbas_tmp(wd)=ibset(intbas_tmp(wd),ibit)
         end do

         OK=.true.
         do wd=1,nwd
            if (intbas(wd)/=intbas_tmp(wd)) then
               OK=.false.
               goto 1000
            endif
         end do

         if (OK) then
            return
         endif

 1000    continue
      end do
      index=-1
      end subroutine get_state_index2

      end module hash_tables


c*** set the input (iunitvec) output (iunitvout) units ***
c*** mjtotal ... 2*M_J, conserved quantity in the SM calculation
c*** mttotal ... Z-N = nprotons-nneutrns
c*** nucleons = nprotons+nneutrns
c*** iparity ... parity=+ iparity=0, parity=- iparity=1
c*** hbo ... hbar*Omega
c*** anu = 1/b^2, b^2=hbar/m/Omega
c*** Jx, Tx ... angular momentum and isopin of a given state
      module intrface
      integer iunitvout,iunitvec,iunitobout
      data iunitvec/3/
      parameter (iunitvout=2,iunitobout=4)
      character(len=8) :: filename,filename2
      integer ki,nki,kf,nkf,kxxx,nkxxx
      data ki/1/,nki/3/,kf/1/,nkf/8/
      integer irestart,jtotal2max,majortot
      data irestart/0/,jtotal2max/4/,majortot/1/
      logical :: formfcal,twobdcal,select,radial,momdist,cluster,antoine
      logical :: mbpt_ncsm,redstick,mfd_james,threebdcal,OLS
      data formfcal/.true./,twobdcal/.false./,select/.false./
      data radial/.false./,momdist/.false./,cluster/.false./,
     +     antoine/.false./,mbpt_ncsm/.false./,redstick/.false./,
     $     mfd_james/.false./,threebdcal/.false./,OLS/.false./
      integer,allocatable :: iseli(:),iself(:),iselxxx(:)
      integer :: ipn,ilast,ilast_tbd
      data ipn/0/,ilast/8/
      character(len=80) :: intbme,twobdop
      data intbme(1:8)/'TBME.int'/
      character(len=4) :: ext1,ext2
      data ext1/'.tmp'/,ext2/'.bak'/
      integer :: nhw_mod,dim_nhw_mod
      end module intrface 

      module obdens
      integer :: mnop,nasps_full=0
      integer,allocatable :: ndiff(:,:) 
      real(kind=kind(0.d0)),allocatable::tdJp(:,:,:,:,:),
     +   tdJn(:,:,:,:,:)
      integer :: ist1dim
      integer,allocatable :: ist1(:,:),iobsind(:)
      integer :: ist2dim,isum2,ist3dim,ist4dim
      integer,allocatable :: ist2(:,:),itbind(:,:,:)
      integer,allocatable :: ist3(:,:),ad3ind(:,:,:,:)
      integer,allocatable :: ist4(:,:),ad4ind(:,:,:,:)
      real(kind=kind(0.d0)),allocatable::t2bd(:,:,:,:,:,:)
      real(kind=kind(0.d0)),allocatable,dimension(:,:,:)::ad_cl,ad_cl_2
      integer,allocatable :: ist2_Jst(:,:)
      integer :: ist2_Jst_J12min,ist2_Jst_J12max

      type de_in
      real(kind(0.d0)),allocatable :: in(:)
      integer :: in_min,in_max
      end type de_in
      type de_mT
      type(de_in),allocatable :: mT(:)
      integer :: mT_min,mT_max
      end type de_mT
      type de_fi
      type(de_mT),allocatable :: fi(:)
      integer :: dim_fi
      end type de_fi
      type de_Jtr
      type(de_fi),allocatable :: Jtr(:)
      integer :: Jtr_min,Jtr_max
      end type de_Jtr
      type(de_Jtr),allocatable :: tbd(:,:)

      contains

      subroutine write_de_Jtr(unit,tb_tmp)
      implicit none
      integer,intent(IN) :: unit
      type(de_Jtr),intent(IN) :: tb_tmp
      integer :: jtrans,ii,mTt
      do jtrans=tb_tmp%Jtr_min,tb_tmp%Jtr_max
         do ii=1,tb_tmp%Jtr(jtrans)%dim_fi
            if (allocated(tb_tmp%Jtr(jtrans)%fi(ii)%mT)) then
               do mTt=tb_tmp%Jtr(jtrans)%fi(ii)%mT_min,
     $              tb_tmp%Jtr(jtrans)%fi(ii)%mT_max
                  write(unit) tb_tmp%Jtr(jtrans)%fi(ii)%mT(mTt)%in
               end do
            endif
         end do
      end do
      end subroutine write_de_Jtr

      subroutine read_de_Jtr(unit,tb_tmp,re_err)
      implicit none
      integer,intent(IN) :: unit
      logical,intent(OUT) :: re_err
      type(de_Jtr) :: tb_tmp
      integer :: jtrans,ii,mTt
      do jtrans=tb_tmp%Jtr_min,tb_tmp%Jtr_max
         do ii=1,tb_tmp%Jtr(jtrans)%dim_fi
            if (allocated(tb_tmp%Jtr(jtrans)%fi(ii)%mT)) then
               do mTt=tb_tmp%Jtr(jtrans)%fi(ii)%mT_min,
     $              tb_tmp%Jtr(jtrans)%fi(ii)%mT_max
                  read(unit,end=9929,err=9929)
     $                 tb_tmp%Jtr(jtrans)%fi(ii)%mT(mTt)%in
               end do
            endif
         end do
      end do
      re_err=.false.
      return
 9929 continue
      re_err=.true.
      end subroutine read_de_Jtr

      end module obdens

      module harmosc
      real(kind=kind(0.d0)) anu,bsquare,bHO
      double precision,allocatable,dimension(:,:,:,:,:) :: rLme
      double precision,allocatable,dimension(:,:,:):: u 
      double precision, allocatable :: rm000(:,:,:),rm001(:,:,:) 
      double precision,allocatable,dimension(:)::sbf0,sbf1
      double precision rc,rinf,rstep
      integer nstep
      parameter (rc=1.d-8,rinf=20.d0,nstep=2000)
c*      parameter (rstep=(rinf-rc)/real(nstep))     ! defined in subroutine
      end module harmosc

      module formfparam
      double precision rmulp0,rmulp1,rmuln0,rmuln1
      double precision rmulpW0,rmulpW1,rmulnW0,rmulnW1
      double precision rmuls0,rmuls1
      end module formfparam

      module kernels
      logical :: NCSMC_kernels=.false.
      integer :: unit_NCSMC=111
      integer :: mult_nucl2_dim,mult_nucl3_dim,mshift2(0:1),
     $     ist2_SDkern_dim,I_ab_max,mult_nucl4_dim,ist3_SDkern_dim,
     $     I_abc_max
      integer, allocatable :: mult_nucl2(:,:),mult_nucl3(:,:),
     $     mult_nucl4(:,:)
      integer,allocatable :: multiple2_dim(:,:,:),
     $     multiple2_point(:,:,:),multiple2_ist(:,:),map_ji(:,:,:),
     $     ist2_SDkern(:,:),ist3_SDkern(:,:)
      type ncsmc_kernel
      integer :: dim_1,start_1
      real(kind(0.d0)),allocatable :: g(:),h_NN(:),h_3N(:)
      end type ncsmc_kernel
      type state_kernel
      integer :: phdim
      integer,allocatable :: phchan(:,:)
      type(ncsmc_kernel),allocatable :: ovl_chan_i(:,:)
      integer :: dim_i,dim_f
      type(ncsmc_kernel),allocatable :: st_fi(:,:)
      end type state_kernel
      type(state_kernel),allocatable :: JpiT(:,:,:)
      integer :: jmi,jma,tmi,tma,ipamin,ipamax
      end module kernels

      module chan_mat_el_definitions
      type symmetric_matrix
      integer :: dim
      real(kind(0.d0)),allocatable :: mat(:)
      real(kind(0.d0)),allocatable :: mat_gen(:,:)
      end type symmetric_matrix
      type Jpi_mat
      type(symmetric_matrix),allocatable :: T_mat(:)
      integer :: T2_min,T2_max
      end type Jpi_mat
      end module chan_mat_el_definitions

      module interaction
      use chan_mat_el_definitions
      type Tz_2b
      type(symmetric_matrix),allocatable :: Tz(:)
      end type Tz_2b
      type T_2b
      type(Tz_2b) :: T(0:1)
      end type T_2b
      type(T_2b),allocatable :: V_2b(:,:)
      type startend
      integer :: start,end
      end type startend
      type j_r
      type(startend),allocatable :: jr(:)
      end type j_r
      type n_c
      type(j_r),allocatable :: nc(:)
      end type n_c
      type N_tot
      integer :: start,end
      type(n_c),allocatable :: lc(:)
      end type N_tot
      type sp_2
      integer :: sp2min,sp2max
      integer,allocatable :: sp2(:)
      end type sp_2
      type sp_1
      integer :: total,sp1min,sp1max,totalrel
      integer,allocatable :: ist(:,:),istrel(:,:)
      type(sp_2),allocatable :: sp1(:)
      type(N_tot),allocatable :: Ntot(:)
      end type sp_1
      type T_tbd
      type(sp_1) :: T(0:1)
      end type T_tbd
      type(T_tbd),allocatable :: tbdst(:,:),tbdst_SDkern(:,:),
     $     thrbdst_SDkern(:,:)
      type n_n
      real(kind(0.d0)),allocatable :: nn(:,:)
      end type n_n
      type l_l
      integer :: lmin,lmax,lminp,lmaxp
      type(n_n),allocatable :: ll(:,:)
      end type l_l
      type T_z
      type(l_l),allocatable :: tz(:)
      end type T_z
      type T_rel
      type(T_z) :: t(0:1)
      end type T_rel
      type(T_rel),allocatable :: V_rel(:,:)
      integer :: num_of_interaction_files
      character(len=80),allocatable :: interaction_file(:)
      integer :: J_sp_int_max,N_sp_int_1max,N_sp_int_2max,
     $     N_sp_int_12max
      end module interaction

      module v3b
      integer :: N1_max=-1,N12_max=-1,N123_max=-1,
     $     isp1ntot,isp3ntot,dim_abc
      integer(8) :: dim_abcdef
      integer,allocatable :: isp1n_cJ(:,:),isp3n_cJ(:,:),sp3ntot(:,:,:)
      real(kind(0.0)),allocatable :: v3b_cJ(:)   ! must be defined in v3b module
      integer,allocatable :: index_abc(:,:,:) !,nlj_orb(:) ! must be defined in v3b module
      integer(8),allocatable :: index_abcdef(:,:),start_abcdef(:) 
      real(kind=kind(0.0)),allocatable :: cgj12(:,:,:,:,:),
     $     cgj123(:,:,:,:,:)
      real(kind=kind(0.0)) :: cgt12(0:1,0:1,0:1),
     $     cgt123(0:1,0:1,-1:1,0:1)
      character(len=80) :: v3intfile
      logical :: V3Nint=.false.
      end module v3b

      program trdensuni
      use initst
      use finast
      use intrface
      use constants
      use paramdef
      use obdens
      use occ
      use hash_tables
      use nodeinfo
      use kernels
      use interaction
      use v3b
***** 9-Dec-1999 *****
***** Calculates one-body transition densities ****
**** Petr Navratil, University of Arizona *****
**** NCSMC kernels, TRIUMF 1/20/2013 ****
      implicit double precision (a-h,o-z)
      include 'mpif.h'
      integer(4),allocatable :: ibasistemp(:,:,:),ibasistempa(:,:)
      logical :: mfdp,mfdi,mfdf,snsp,sndp,dnuc,obscal,irem_cal
      real(kind=kind(0.d0)) :: denom
      integer :: ncut,mionb,ionb,ierr,p_state,n_state,i1,i2
      integer,allocatable :: occ_temp(:)
      logical :: mbpi,mbpf
*****************
      interface
         subroutine confwavread(nsd,mxnwd,nhom,mjtotal,mttotal,iparity,
     +        nhw,
     +        nucleons,nprotons,nneutrns,jt2,it2,Jx,Tx,hbo,iloc,nshll,
     +        nhom12,nhom123,nhom1234,nasps,mxsps,mconf,nconf,nprtn,
     +        nneut,iendconf,
     +        n_sp,l_sp,j2_sp,m2_sp,mt2_sp,bmp,ibasis,ener)
         integer(2),pointer :: iloc(:,:),mconf(:)
         integer(4),pointer:: jt2(:),it2(:)
         integer(4),pointer :: nconf,nprtn(:,:),nneut(:,:),iendconf(:)
         integer(4),pointer :: n_sp(:),l_sp(:),
     +        j2_sp(:),m2_sp(:),mt2_sp(:)
         real(kind=kind(0.0d0)),pointer::bmp(:,:)
         integer(4),pointer :: ibasis(:,:,:)
         integer(4),pointer :: nsd,mxnwd,nhom,nhom12,nasps,mxsps,
     +        nhom123,nhom1234
         integer(4),pointer :: mjtotal,mttotal,iparity,nshll,nhw
         integer(4),pointer :: nucleons,nprotons,nneutrns
         real(kind=kind(0.0d0)),pointer :: Jx(:),Tx(:),ener(:)
         real(kind=kind(0.0d0)),pointer :: hbo
         end subroutine confwavread
         subroutine confwave_antoine(nsd,mxnwd,nhom,mjtotal,
     +        mttotal,iparity,nhw,
     +        nucleons,nprotons,nneutrns,jt2,Jx,
     +        I_state,nsd_p,nsd_n,occ_p,occ_n,ibas_p,ibas_n,nshll,
     +        nhom12,nhom123,nhom1234,nasps,mxsps,mconf,nconf,nprtn,
     +        nneut,iendconf,
     +        n_sp,l_sp,j2_sp,m2_sp,mt2_sp,bmp,ener,i_case)
         integer(2),pointer :: mconf(:)
         integer(4),pointer:: I_state(:,:),nsd_p,nsd_n,
     +        occ_p(:,:),occ_n(:,:),ibas_p(:,:),ibas_n(:,:)
         integer(4),pointer:: jt2(:)
         integer(4),pointer :: nconf,nprtn(:,:),nneut(:,:),iendconf(:)
         integer(4),pointer :: n_sp(:),l_sp(:),
     +        j2_sp(:),m2_sp(:),mt2_sp(:)
         real(kind=kind(0.0d0)),pointer::bmp(:,:)
         integer(4),pointer :: nsd,mxnwd,nhom,nhom12,nasps,mxsps,
     +        nhom123,nhom1234
         integer(4),pointer :: mjtotal,mttotal,iparity,nshll,nhw
         integer(4),pointer :: nucleons,nprotons,nneutrns
         real(kind=kind(0.0d0)),pointer :: Jx(:),ener(:)
         integer :: i_case
         end subroutine confwave_antoine
         subroutine confwave_mbpt(nsd,mxnwd,nhom,mjtotal,mttotal,
     $        iparity,nhw,
     +        nucleons,nprotons,nneutrns,jt2,it2,Jx,Tx,hbo,iloc,nshll,
     +        nhom12,nhom123,nhom1234,nasps,mxsps,mconf,nconf,nprtn,
     +        nneut,iendconf,
     +        n_sp,l_sp,j2_sp,m2_sp,mt2_sp,bmp,ibasis,ener)
         integer(2),pointer :: iloc(:,:),mconf(:)
         integer(4),pointer:: jt2(:),it2(:)
         integer(4),pointer :: nconf,nprtn(:,:),nneut(:,:),iendconf(:)
         integer(4),pointer :: n_sp(:),l_sp(:),
     +        j2_sp(:),m2_sp(:),mt2_sp(:)
         real(kind=kind(0.0d0)),pointer::bmp(:,:)
         integer(4),pointer :: ibasis(:,:,:)
         integer(4),pointer :: nsd,mxnwd,nhom,nhom12,nasps,mxsps,
     +        nhom123,nhom1234
         integer(4),pointer :: mjtotal,mttotal,iparity,nshll,nhw
         integer(4),pointer :: nucleons,nprotons,nneutrns
         real(kind=kind(0.0d0)),pointer :: Jx(:),Tx(:),ener(:)
         real(kind=kind(0.0d0)),pointer :: hbo
         end subroutine confwave_mbpt
         subroutine confwave_redstick(nsd,mxnwd,nhom,mjtotal,mttotal,
     $        iparity,nhw,
     +        nucleons,nprotons,nneutrns,jt2,it2,Jx,Tx,hbo,iloc,nshll,
     +        nhom12,nhom123,nhom1234,nasps,mxsps,mconf,nconf,nprtn,
     +        nneut,iendconf,
     +        n_sp,l_sp,j2_sp,m2_sp,mt2_sp,bmp,ibasis,ener)
         integer(2),pointer :: iloc(:,:),mconf(:)
         integer(4),pointer:: jt2(:),it2(:)
         integer(4),pointer :: nconf,nprtn(:,:),nneut(:,:),iendconf(:)
         integer(4),pointer :: n_sp(:),l_sp(:),
     +        j2_sp(:),m2_sp(:),mt2_sp(:)
         real(kind=kind(0.0d0)),pointer::bmp(:,:)
         integer(4),pointer :: ibasis(:,:,:)
         integer(4),pointer :: nsd,mxnwd,nhom,nhom12,nasps,mxsps,
     +        nhom123,nhom1234
         integer(4),pointer :: mjtotal,mttotal,iparity,nshll,nhw
         integer(4),pointer :: nucleons,nprotons,nneutrns
         real(kind=kind(0.0d0)),pointer :: Jx(:),Tx(:),ener(:)
         real(kind=kind(0.0d0)),pointer :: hbo
         end subroutine confwave_redstick
         subroutine confwave_mfdj(nsd,mxnwd,nhom,mjtotal,mttotal,
     $        iparity,nhw,
     +        nucleons,nprotons,nneutrns,jt2,it2,Jx,Tx,hbo,iloc,nshll,
     +        nhom12,nhom123,nhom1234,nasps,mxsps,mconf,nconf,nprtn,
     +        nneut,iendconf,
     +        n_sp,l_sp,j2_sp,m2_sp,mt2_sp,bmp,ibasis,ener)
         integer(2),pointer :: iloc(:,:),mconf(:)
         integer(4),pointer:: jt2(:),it2(:)
         integer(4),pointer :: nconf,nprtn(:,:),nneut(:,:),iendconf(:)
         integer(4),pointer :: n_sp(:),l_sp(:),
     +        j2_sp(:),m2_sp(:),mt2_sp(:)
         real(kind=kind(0.0d0)),pointer::bmp(:,:)
         integer(4),pointer :: ibasis(:,:,:)
         integer(4),pointer :: nsd,mxnwd,nhom,nhom12,nasps,mxsps,
     +        nhom123,nhom1234
         integer(4),pointer :: mjtotal,mttotal,iparity,nshll,nhw
         integer(4),pointer :: nucleons,nprotons,nneutrns
         real(kind=kind(0.0d0)),pointer :: Jx(:),Tx(:),ener(:)
         real(kind=kind(0.0d0)),pointer :: hbo
         end subroutine confwave_mfdj
      end interface 

      call MPI_INIT(ierr)
      icomm = MPI_COMM_WORLD
      call MPI_COMM_RANK(icomm,iproc,ierr)
      call MPI_COMM_SIZE(icomm,nproc,ierr)

      inquire(file='trdens.in',exist=obscal)
      if (obscal) then
         open(1,file='trdens.in',status='old')
         read(1,*) select
         read(1,*) ki,nki
         if (select) then
            allocate (iseli(nki))
            read(1,*,end=1998,err=1998) (iseli(ii),ii=1,nki)
            goto 1999
 1998       continue
            iseli(:)=(/(ii,ii=ki,ki+nki-1)/)
 1999       continue
         else
            read(1,*)
         endif
         read(1,*) kf,nkf
         if (select) then
            allocate (iself(nkf))
            read(1,*,end=2998,err=2998) (iself(ii),ii=1,nkf)
            goto 2999
 2998       continue
            iself(:)=(/(ii,ii=kf,kf+nkf-1)/)
 2999       continue
         else
            read(1,*)
         endif
         read(1,*) jtotal2max
         read(1,*) irestart
         if (irestart==1) then
            ext1='.tmp'
            ext2='.bak'
         elseif (irestart==2) then
            ext1='.bak'
            ext2='.tmp'
         endif
         read(1,*) majortot
         read(1,*) irem_cal
         read(1,*) formfcal
         read(1,*) radial
         read(1,*) momdist
         read(1,*) twobdcal
         if (twobdcal) then
            read(1,*) ipn
            read(1,'(a)') intbme
            ilast=index(intbme,' ')-1
            read(1,'(a)') twobdop
            ilast_tbd=index(twobdop,' ')-1
         else
            read(1,*,end=8889)
            read(1,*,end=8889)
            read(1,*,end=8889)
         endif   
         read(1,*,end=8889) cluster
         read(1,*,end=8889) antoine
         if (antoine) then
            allocate(hboi)
            read(1,*) hboi
            read(1,*) num_of_in_i
            allocate(it2i(num_of_in_i))
            allocate(Txi(num_of_in_i))
            read(1,*) (it2i(ii),ii=1,num_of_in_i)
            do ii=1,num_of_in_i
               Txi(ii)=it2i(ii)/2.0d0
            end do
            read(1,*) num_of_in_f
            if (num_of_in_f>0) then
               allocate(it2f(num_of_in_f))
               allocate(Txf(num_of_in_f))
               read(1,*) (it2f(ii),ii=1,num_of_in_f)
               do ii=1,num_of_in_f
                  Txf(ii)=it2f(ii)/2.0d0
               end do
            endif
            irem_cal=.false.
            majortot=3
         else
            read(1,*,end=8889)
            read(1,*,end=8889)
            read(1,*,end=8889)
            read(1,*,end=8889)
            read(1,*,end=8889)
         endif
         read(1,*,end=8889) mbpt_ncsm
         read(1,*,end=8889) redstick
         read(1,*,end=8889) mfd_james
         read(1,*,end=8889) threebdcal
         read(1,*,end=8889) NCSMC_kernels
         if (NCSMC_kernels) then
            read(1,*) num_of_interaction_files
            allocate(interaction_file(num_of_interaction_files))
            do jj=1,num_of_interaction_files
               read(1,'(a)') interaction_file(jj)
            end do
            read(1,*,end=8889) V3Nint
            if (V3Nint) then
               read(1,*,end=8889) N1_max,N12_max,N123_max
               read(1,'(a)') v3intfile
            else
               read(1,*,end=8889)
               read(1,*,end=8889)
            endif
         else
            read(1,*,end=8889)
            read(1,*,end=8889)
            read(1,*,end=8889)
            read(1,*,end=8889)
            read(1,*,end=8889)
            read(1,*,end=8889)
            read(1,*,end=8889)
         endif

         read(1,*,end=8889) OLS
         if (OLS) then
            read(1,*,end=8889) nhw_mod
            read(1,*,end=8889) dim_nhw_mod
         else
            read(1,*,end=8889)
            read(1,*,end=8889)
         endif

 8889    continue
         close(1)
      endif

      inquire(file='trdens.out',exist=obscal)
      call MPI_Barrier(icomm,ierr)
      if (obscal) goto 33333

      if (iproc==0) then
         open(iunitvout,file='trdens.out',status='new')
         WRITE(iunitvout,2000)
 2000    FORMAT(' OBDME calculation')
         if (threebdcal) then
            open(44,file='trdens.out_3bd.bin',status='new',
     $           form='unformatted',action='write')
         endif
         if (NCSMC_kernels) then
            open(unit_NCSMC,file='NCSMC_kernels.dat',status='new',
     $           form='formatted',action='write')
         endif
      endif   

      if (antoine) then
         inquire(file='anto.egv',exist=mfdp)
         inquire(file='anti.egv',exist=mfdi)
         inquire(file='antf.egv',exist=mfdf)
      elseif (mbpt_ncsm) then
         inquire(file='mbpt.egv',exist=mfdp)
         inquire(file='mbpi.egv',exist=mfdi)
         inquire(file='mbpf.egv',exist=mfdf)
         majortot=2
      elseif (redstick) then
         inquire(file='reds.egv',exist=mfdp)
         inquire(file='redi.egv',exist=mfdi)
         inquire(file='redf.egv',exist=mfdf)
      elseif (mfd_james) then
         inquire(file='mfdj.egv',exist=mfdp)
         inquire(file='mfji.egv',exist=mfdi)
         inquire(file='mfjf.egv',exist=mfdf)
      else
         inquire(file='mfdp.egv',exist=mfdp)
         inquire(file='mfdi.egv',exist=mfdi)
         inquire(file='mfdf.egv',exist=mfdf)
      endif

      if (NCSMC_kernels) then
         inquire(file='mbpi.egv',exist=mbpi)
         inquire(file='mfdf.egv',exist=mbpf)
         if (mbpf) mfdf=.true.
      endif

      if (mfdp) then
         snsp=.true.
         sndp=.false.
         dnuc=.false.
      elseif (mfdi.and.mfdf.or.(NCSMC_kernels.and.mbpi.and.mfdf)) then
         snsp=.false.
         sndp=.false.
         dnuc=.false.
      else
          print *,'*** error: input files do not exist'
          stop
      endif

      if (mfdp) then
c*** initialization routine that reads in the mfdp.egv file
c*** with the NC SM wave function
         if (iproc==0) then
            write(iunitvout,*) mfdp
            if (antoine) then
               write(iunitvout,*)
     +              'Wave functions read from anto.egv file'
            elseif (mbpt_ncsm) then
               write(iunitvout,*)
     +              'Wave functions read from mbpt.egv file'
               select=.false.
c               ki=1
c               nki=1
c               kf=1
c               nkf=1
            elseif (redstick) then
               write(iunitvout,*)
     +              'Wave functions read from reds.egv file'
            elseif (mfd_james) then
               write(iunitvout,*)
     +              'Wave functions read from mfdj.egv file'
            else
               write(iunitvout,*)
     +              'Wave functions read from mfdp.egv file'
            endif
         endif
         kxxx=min(ki,kf)
         nkxxx=max(ki+nki,kf+nkf)-kxxx
         if (select) then
            if (ki/=kf) stop 'wrong input'
            allocate(iselxxx(nkxxx))
            if (nki>nkf) then
               iselxxx=iseli
            else
               iselxxx=iself
            endif
         endif
         if (antoine) then
            filename='anto.egv'
            call confwave_antoine(nsdi,mxnwdi,nhomi,mjtotali,
     $           mttotali,iparityi,nhwi,
     +           nucleonsi,nprotonsi,nneutrnsi,jt2i,Jxi,
     +           I_state_i,nsd_p_i,nsd_n_i,occ_p_i,occ_n_i,
     +           ibas_p_i,ibas_n_i,nshlli,
     +           nhom12i,nhom123i,nhom1234i,naspsi,mxspsi,mconfi,nconfi,
     +           nprtni,nneuti,iendconfi,
     +           n_spi,l_spi,j2_spi,m2_spi,mt2_spi,bmpi,eneri,1)
         elseif (mbpt_ncsm) then
            filename='mbpt.egv'
            call confwave_mbpt(nsdi,mxnwdi,nhomi,mjtotali,mttotali,
     $           iparityi,nhwi,
     +           nucleonsi,nprotonsi,nneutrnsi,jt2i,it2i,Jxi,Txi,hboi,
     +           iloci,nshlli,
     +           nhom12i,nhom123i,nhom1234i,naspsi,mxspsi,mconfi,nconfi,
     +           nprtni,nneuti,iendconfi,
     +           n_spi,l_spi,j2_spi,m2_spi,mt2_spi,bmpi,ibasisi,eneri)
         elseif (redstick) then
            filename='reds.egv'
            call confwave_redstick(nsdi,mxnwdi,nhomi,mjtotali,mttotali,
     $           iparityi,nhwi,
     +           nucleonsi,nprotonsi,nneutrnsi,jt2i,it2i,Jxi,Txi,hboi,
     +           iloci,nshlli,
     +           nhom12i,nhom123i,nhom1234i,naspsi,mxspsi,mconfi,nconfi,
     +           nprtni,nneuti,iendconfi,
     +           n_spi,l_spi,j2_spi,m2_spi,mt2_spi,bmpi,ibasisi,eneri)
         elseif (mfd_james) then
            filename='mfdj.egv'
            filename2='mf2j.egv'
            call confwave_mfdj(nsdi,mxnwdi,nhomi,mjtotali,mttotali,
     $           iparityi,nhwi,
     +           nucleonsi,nprotonsi,nneutrnsi,jt2i,it2i,Jxi,Txi,hboi,
     +           iloci,nshlli,
     +           nhom12i,nhom123i,nhom1234i,naspsi,mxspsi,mconfi,nconfi,
     +           nprtni,nneuti,iendconfi,
     +           n_spi,l_spi,j2_spi,m2_spi,mt2_spi,bmpi,ibasisi,eneri)
         else
            filename='mfdp.egv'
            call confwavread(nsdi,mxnwdi,nhomi,mjtotali,mttotali,
     $           iparityi,nhwi,
     +           nucleonsi,nprotonsi,nneutrnsi,jt2i,it2i,Jxi,Txi,hboi,
     +           iloci,nshlli,
     +           nhom12i,nhom123i,nhom1234i,naspsi,mxspsi,mconfi,nconfi,
     +           nprtni,nneuti,iendconfi,
     +           n_spi,l_spi,j2_spi,m2_spi,mt2_spi,bmpi,ibasisi,eneri)
         endif
         nsdf=>nsdi
         mxnwdf=>mxnwdi
         nhomf=>nhomi
         nhom12f=>nhom12i
         nhom123f=>nhom123i
         nhom1234f=>nhom1234i
         naspsf=>naspsi
         mxspsf=>mxspsi
         mjtotalf=>mjtotali
         mttotalf=>mttotali
         iparityf=>iparityi
         nhwf=>nhwi
         nucleonsf=>nucleonsi
         nprotonsf=>nprotonsi
         nneutrnsf=>nneutrnsi
         jt2f=>jt2i
         it2f=>it2i
         Jxf=>Jxi
         Txf=>Txi
         hbof=>hboi
         bmpf=>bmpi
         enerf=>eneri
         n_spf=>n_spi
         l_spf=>l_spi
         j2_spf=>j2_spi
         m2_spf=>m2_spi
         mt2_spf=>mt2_spi
         nshllf=>nshlli
         mconff=>mconfi
         nconff=>nconfi
         nprtnf=>nprtni
         nneutf=>nneuti
         iendconff=>iendconfi
         if (nki>nkxxx) nki=nkxxx
         if (nkf>nkxxx) nkf=nkxxx

         if (iproc==0) then
            write(iunitvout,1347) ki,ki+nki-1
 1347       format(/,' wave functions of the states #',i3,'- #',i3,
     +        ' used')
            write(iunitvout,1347) kf,kf+nkf-1
         endif
         if (majortot==2) then
            ilocf=>iloci
         elseif (majortot==3) then
            I_state_f=>I_state_i
            nsd_p_f=>nsd_p_i
            nsd_n_f=>nsd_n_i
            occ_p_f=>occ_p_i
            occ_n_f=>occ_n_i
            ibas_p_f=>ibas_p_i
            ibas_n_f=>ibas_n_i
         else
            ibasisf=>ibasisi
            if (associated(iloci)) deallocate(iloci)
         endif

      else
         if (antoine) then
            filename='anti.egv'
         elseif (mbpt_ncsm) then
            filename='mbpi.egv'
c            ki=1
c            nki=1
            select=.false.
         elseif (redstick) then
            filename='redi.egv'
         elseif (mfd_james) then
            filename='mfji.egv'
            filename2='mf2i.egv'
         else
            filename='mfdi.egv'
         endif
         if (NCSMC_kernels.and.mbpi.and.mfdf) then
            filename='mbpi.egv'
         endif
         kxxx=ki
         nkxxx=nki
         if (select) then
            if (allocated(iselxxx)) deallocate(iselxxx)
            allocate(iselxxx(nkxxx))
            iselxxx=iseli
            print *, ' iseli=',iselxxx
         endif
         if (iproc==0) then
            write(iunitvout,*) mfdp
            if (antoine) then
               write(iunitvout,*) 
     +           'Wave functions read from anti.egv and antf.egv files'
            elseif (mbpt_ncsm) then
               write(iunitvout,*) 
     +          'Wave functions read from mbpi.egv and mbpf.egv files'
            elseif (redstick) then
               write(iunitvout,*) 
     +        'Wave functions read from redi.egv and redf.egv files'
            elseif (mfd_james) then
               write(iunitvout,*) 
     +          'Wave functions read from mfji.egv and mfjf.egv files'
            else
               write(iunitvout,*) 
     +           'Wave functions read from mfdi.egv and mfdf.egv files'
            endif
            write(iunitvout,*) 
            write(iunitvout,*) 
     +'***************************************************************'
            write(iunitvout,*) '*** Initial nucleus ***'
            write(iunitvout,*) '***********************'
         endif
         if (antoine) then
            call confwave_antoine(nsdi,mxnwdi,nhomi,mjtotali,
     $           mttotali,iparityi,nhwi,
     +           nucleonsi,nprotonsi,nneutrnsi,jt2i,Jxi,
     +           I_state_i,nsd_p_i,nsd_n_i,occ_p_i,occ_n_i,
     +           ibas_p_i,ibas_n_i,nshlli,
     +           nhom12i,nhom123i,nhom1234i,naspsi,mxspsi,mconfi,nconfi,
     +           nprtni,nneuti,iendconfi,
     +           n_spi,l_spi,j2_spi,m2_spi,mt2_spi,bmpi,eneri,1)   
         elseif (mbpt_ncsm.or.(NCSMC_kernels.and.mbpi.and.mfdf)) then
            call confwave_mbpt(nsdi,mxnwdi,nhomi,mjtotali,mttotali,
     $           iparityi,nhwi,
     +           nucleonsi,nprotonsi,nneutrnsi,jt2i,it2i,Jxi,Txi,hboi,
     +           iloci,nshlli,
     +           nhom12i,nhom123i,nhom1234i,naspsi,mxspsi,mconfi,nconfi,
     +           nprtni,nneuti,iendconfi,
     +           n_spi,l_spi,j2_spi,m2_spi,mt2_spi,bmpi,ibasisi,eneri)
         elseif (redstick) then
            call confwave_redstick(nsdi,mxnwdi,nhomi,mjtotali,mttotali,
     $           iparityi,nhwi,
     +           nucleonsi,nprotonsi,nneutrnsi,jt2i,it2i,Jxi,Txi,hboi,
     +           iloci,nshlli,
     +           nhom12i,nhom123i,nhom1234i,naspsi,mxspsi,mconfi,nconfi,
     +           nprtni,nneuti,iendconfi,
     +           n_spi,l_spi,j2_spi,m2_spi,mt2_spi,bmpi,ibasisi,eneri)
         elseif (mfd_james) then
            call confwave_mfdj(nsdi,mxnwdi,nhomi,mjtotali,mttotali,
     $           iparityi,nhwi,
     +           nucleonsi,nprotonsi,nneutrnsi,jt2i,it2i,Jxi,Txi,hboi,
     +           iloci,nshlli,
     +           nhom12i,nhom123i,nhom1234i,naspsi,mxspsi,mconfi,nconfi,
     +           nprtni,nneuti,iendconfi,
     +           n_spi,l_spi,j2_spi,m2_spi,mt2_spi,bmpi,ibasisi,eneri)
         else
            iunitvec=3
            call confwavread(nsdi,mxnwdi,nhomi,mjtotali,
     $           mttotali,iparityi,nhwi,
     +           nucleonsi,nprotonsi,nneutrnsi,jt2i,it2i,Jxi,Txi,hboi,
     +           iloci,nshlli,
     +           nhom12i,nhom123i,nhom1234i,naspsi,mxspsi,mconfi,nconfi,
     +           nprtni,nneuti,iendconfi,
     +           n_spi,l_spi,j2_spi,m2_spi,mt2_spi,bmpi,ibasisi,eneri)   
         endif
         if (nki>nkxxx) nki=nkxxx
         if (iproc==0) then         
            write(iunitvout,1347) ki,ki+nki-1
c*      do i1=1,1000
c*         write(iunitvout,3434) (iloci(ia,i1),ia=1,nucleonsi)
c* 3434    format(10i6)
c*      end do
         endif   
         if (majortot/=2.and.associated(iloci)) deallocate(iloci)

         if (antoine) then
            filename='antf.egv'
         elseif (mbpt_ncsm) then
            filename='mbpf.egv'
c            kf=1
c            nkf=1
            select=.false.
         elseif (redstick) then
            filename='redf.egv'
         elseif (mfd_james) then
            filename='mfjf.egv'
            filename2='mf2f.egv'
         else
            filename='mfdf.egv'
         endif
         kxxx=kf
         nkxxx=nkf
         if (select) then
            if (allocated(iselxxx)) deallocate(iselxxx)
            allocate(iselxxx(nkxxx))
            iselxxx=iself
            print *, ' iself=',iselxxx
         endif
         if (iproc==0) then
            write(iunitvout,*) 
            write(iunitvout,*) 
     +'***************************************************************'
            write(iunitvout,*) '*** Final nucleus ***'
            write(iunitvout,*) '*********************'
         endif   
         if (antoine) then
            call confwave_antoine(nsdf,mxnwdf,nhomf,mjtotalf,
     $           mttotalf,iparityf,nhwf,
     +           nucleonsf,nprotonsf,nneutrnsf,jt2f,Jxf,
     +           I_state_f,nsd_p_f,nsd_n_f,occ_p_f,occ_n_f,
     +           ibas_p_f,ibas_n_f,nshllf,
     +           nhom12f,nhom123f,nhom1234f,naspsf,mxspsf,mconff,nconff,
     +           nprtnf,nneutf,iendconff,
     +           n_spf,l_spf,j2_spf,m2_spf,mt2_spf,bmpf,enerf,2)
            hbof=>hboi
         elseif (mbpt_ncsm) then
            call confwave_mbpt(nsdf,mxnwdf,nhomf,mjtotalf,
     $           mttotalf,iparityf,nhwf,
     +           nucleonsf,nprotonsf,nneutrnsf,jt2f,it2f,Jxf,Txf,
     $           hbof,ilocf,nshllf,
     +           nhom12f,nhom123f,nhom1234f,naspsf,mxspsf,mconff,nconff,
     +           nprtnf,nneutf,iendconff,
     +           n_spf,l_spf,j2_spf,m2_spf,mt2_spf,bmpf,ibasisf,enerf)
         elseif (redstick) then
            call confwave_redstick(nsdf,mxnwdf,nhomf,mjtotalf,
     $           mttotalf,iparityf,nhwf,
     +           nucleonsf,nprotonsf,nneutrnsf,jt2f,it2f,Jxf,Txf,
     $           hbof,ilocf,nshllf,
     +           nhom12f,nhom123f,nhom1234f,naspsf,mxspsf,mconff,nconff,
     +           nprtnf,nneutf,iendconff,
     +           n_spf,l_spf,j2_spf,m2_spf,mt2_spf,bmpf,ibasisf,enerf)
         elseif (mfd_james) then
            call confwave_mfdj(nsdf,mxnwdf,nhomf,mjtotalf,
     $           mttotalf,iparityf,nhwf,
     +           nucleonsf,nprotonsf,nneutrnsf,jt2f,it2f,Jxf,Txf,
     $           hbof,ilocf,nshllf,
     +           nhom12f,nhom123f,nhom1234f,naspsf,mxspsf,mconff,nconff,
     +           nprtnf,nneutf,iendconff,
     +           n_spf,l_spf,j2_spf,m2_spf,mt2_spf,bmpf,ibasisf,enerf)
         else
            iunitvec=378
            call confwavread(nsdf,mxnwdf,nhomf,mjtotalf,
     $           mttotalf,iparityf,nhwf,
     +           nucleonsf,nprotonsf,nneutrnsf,jt2f,it2f,Jxf,Txf,
     $           hbof,ilocf,nshllf,
     +           nhom12f,nhom123f,nhom1234f,naspsf,mxspsf,mconff,nconff,
     +           nprtnf,nneutf,iendconff,
     +           n_spf,l_spf,j2_spf,m2_spf,mt2_spf,bmpf,ibasisf,enerf)
         endif
         if (nkf>nkxxx) nkf=nkxxx
         if (iproc==0) then         
            write(iunitvout,1347) kf,kf+nkf-1
         endif   
         if (majortot/=2.and.associated(ilocf)) deallocate(ilocf)
         nhomi=max(nhomi,nhomf)
         nhomf=>nhomi
         nhom12i=max(nhom12i,nhom12f)
         nhom12f=>nhom12i

         if (hbof/=hboi) then
            if (iproc==0) then
            write(iunitvout,*) ' hboi does not correspond to hbof',
     +             hboi,hbof
            endif
cc            stop
         endif 

         do ii=1,min(naspsi,naspsf)
            if (n_spi(ii)/=n_spf(ii)) then
               print *,'*** error: n_spi,n_spf=',n_spi(ii),n_spf(ii)
               stop
            endif   
            if (l_spi(ii)/=l_spf(ii)) then
               print *,'*** error: l_spi,l_spf=',l_spi(ii),l_spf(ii)
               stop
            endif   
            if (j2_spi(ii)/=j2_spf(ii)) then
               print *,'*** error: j2_spi,j2_spf=',j2_spi(ii),
     +              j2_spf(ii)
               stop
            endif   
            if (m2_spi(ii)/=m2_spf(ii)) then
               print *,'*** error: m2_spi,m2_spf=',m2_spi(ii),
     +              m2_spf(ii)
               stop
            endif   
            if (mt2_spi(ii)/=mt2_spf(ii)) then
               print *,'*** error: mt2_spi,mt2_spf=',mt2_spi(ii),
     +              mt2_spf(ii)
               stop
            endif   
         end do

         if (mxnwdf>mxnwdi) then
            if (majortot<2) then
               allocate(ibasistemp(mxnwdf,2,nsdi))
               ibasistemp=0
               ibasistemp(1:mxnwdi,:,:)=ibasisi(1:mxnwdi,:,:)
               deallocate(ibasisi)
               allocate(ibasisi(mxnwdf,2,nsdi))
               ibasisi=ibasistemp
               deallocate(ibasistemp)
            elseif (majortot==3) then
               allocate(ibasistempa(mxnwdf,nsd_p_i))
               ibasistempa=0
               ibasistempa(1:mxnwdi,:)=ibas_p_i(1:mxnwdi,:)
               deallocate(ibas_p_i)
               allocate(ibas_p_i(mxnwdf,nsd_p_i))
               ibas_p_i=ibasistempa
               deallocate(ibasistempa)
               allocate(ibasistempa(mxnwdf,nsd_n_i))
               ibasistempa=0
               ibasistempa(1:mxnwdi,:)=ibas_n_i(1:mxnwdi,:)
               deallocate(ibas_n_i)
               allocate(ibas_n_i(mxnwdf,nsd_n_i))
               ibas_n_i=ibasistempa
               deallocate(ibasistempa)
            endif
            mxnwdi=>mxnwdf
            n_spi=>n_spf
            l_spi=>l_spf
            j2_spi=>j2_spf
            m2_spi=>m2_spf
            mt2_spi=>mt2_spf
         elseif (mxnwdi>mxnwdf) then
            if (majortot<2) then
               allocate(ibasistemp(mxnwdi,2,nsdf))
               ibasistemp=0
               ibasistemp(1:mxnwdf,:,:)=ibasisf(1:mxnwdf,:,:)
               deallocate(ibasisf)
               allocate(ibasisf(mxnwdi,2,nsdf))
               ibasisf=ibasistemp
               deallocate(ibasistemp)
            elseif (majortot==3) then
               allocate(ibasistempa(mxnwdi,nsd_p_f))
               ibasistempa=0
               ibasistempa(1:mxnwdf,:)=ibas_p_f(1:mxnwdf,:)
               deallocate(ibas_p_f)
               allocate(ibas_p_f(mxnwdi,nsd_p_f))
               ibas_p_f=ibasistempa
               deallocate(ibasistempa)
               allocate(ibasistempa(mxnwdi,nsd_n_f))
               ibasistempa=0
               ibasistempa(1:mxnwdf,:)=ibas_n_f(1:mxnwdf,:)
               deallocate(ibas_n_f)
               allocate(ibas_n_f(mxnwdi,nsd_n_f))
               ibas_n_f=ibasistempa
               deallocate(ibasistempa)
            endif
            mxnwdf=>mxnwdi
            n_spf=>n_spi
            l_spf=>l_spi
            j2_spf=>j2_spi
            m2_spf=>m2_spi
            mt2_spf=>mt2_spi
         endif 

         select case(iabs(nucleonsi-nucleonsf))
      case(0)
         if (nprotonsi==nprotonsf) then
            if (iparityi/=iparityf) sndp=.true.
            if (iparityi==iparityf.and.nsdi==nsdf) snsp=.true.
         elseif (abs(nprotonsi-nprotonsf)==1
     $           .or.abs(nprotonsi-nprotonsf)==2) then   
            dnuc=.true.
         else
            if (iproc==0) then
            write(iunitvout,*) '*** incompatible number of nucleons'
            endif
            stop
         endif   
      case(1:4)
         dnuc=.true.
      case default
         if (iproc==0) then
            write(iunitvout,*) 
     +           ' nucleonsi and nucleonsf out of range',
     +           nucleonsi,nucleonsf
         endif
         stop
         end select
            
         naspsi=max(naspsi,naspsf)
         naspsf=>naspsi

      endif

      call onebodyinit

      if (OLS) then
         call OLS_effective_hamiltonian
      endif

      call MPI_Barrier(icomm,ierr)
      if (NCSMC_kernels) then

         if (nucleonsi==nucleonsf+1) then
            N_cluster_max=nhom12i
            N_cluster_min=nhom12f
            N_proj_max=0
            N_sp_min=nhomi !max(nhomi,N_proj_max)
            num_neutr_proj=nneutrnsi-nneutrnsf
            num_prot_proj=nprotonsi-nprotonsf
            N_RGM_max=min(nhwf+1,nhwi)
     $           -(Nmin_HO(nprotonsf)
     $           +Nmin_HO(nneutrnsf))
     $           -(Nmin_HO(num_prot_proj)
     $           +Nmin_HO(num_neutr_proj))
         elseif (nucleonsi==nucleonsf+2) then
            N_cluster_max=nhom12i
            N_cluster_min=nhom12f
            N_proj_max=nhom12f
            N_sp_min=nhomi !max(nhomi,N_proj_max)
            num_neutr_proj=nneutrnsi-nneutrnsf
            num_prot_proj=nprotonsi-nprotonsf
            N_RGM_max=min(nhwf+1,nhwi)
     $           -(Nmin_HO(nprotonsf)
     $           +Nmin_HO(nneutrnsf))
     $           -(Nmin_HO(num_prot_proj)
     $           +Nmin_HO(num_neutr_proj))
         elseif (nucleonsi==nucleonsf+3) then
            N_cluster_max=nhom12i
            N_cluster_min=nhom12f
            N_proj_max=nhom12f
            N_sp_min=nhomi !max(nhomi,N_proj_max)
            num_neutr_proj=nneutrnsi-nneutrnsf
            num_prot_proj=nprotonsi-nprotonsf
            N_RGM_max=min(nhwf+1,nhwi)
     $           -(Nmin_HO(nprotonsf)
     $           +Nmin_HO(nneutrnsf))
     $           -(Nmin_HO(num_prot_proj)
     $           +Nmin_HO(num_neutr_proj))
         else
            if (iproc==0) print *,
     $           ' not implemented for nucelonsi,nucleonsf=',
     $           nucleonsi,nucleonsf
            stop
         endif
         if (iproc==0) then
            print *,' N_cluster_max=',N_cluster_max
            print *,' N_sp_min=',N_sp_min
            print *,' N_proj_max=',N_proj_max
            print *,' num_prot_proj=',num_prot_proj
            print *,' num_neutr_proj=',num_neutr_proj
            print *,' N_RGM_max=',N_RGM_max
         endif
c**** call gamasub before first use
         lrelm=nhwi+N_proj_max+10 !Targe(1,1)%nhom+Projectile(1)%p_chann(1)%nhom+10
         nrmax=50
         lmax=lrelm
         maxgam=2*nrmax+2*lmax+3
         call gamasub
         if (iproc==0) print *,' gamasub called, iproc=',iproc
c***  MPI
         call MPI_Barrier(icomm,ierr)
c***  MPI
         call global_ist1_setup
         if (iproc==0) print *,' global_ist1_setup called, iproc=',
     $        iproc
         call check_jlevel_consistency
         if (iproc==0) print *,
     $        ' check_jlevel_consistency called, iproc=',
     $        iproc
c*** MPI
         call MPI_Barrier(icomm,ierr)
c*** MPI
         call interaction_read
         if (iproc==0) print *, ' interaction_read called, iproc=',iproc
c*** MPI
         call MPI_Barrier(icomm,ierr)
c*** MPI
         call two_body_state_setup
         if (iproc==0) print *, ' two_body_state_setup called, iproc=',
     $        iproc
c*** MPI
         call MPI_Barrier(icomm,ierr)
c*** MPI
         call V_rel_V_2b_trans
         if (iproc==0) print *, ' V_rel_V_2b_trans called, iproc=',iproc
c*** MPI
         call MPI_Barrier(icomm,ierr)
c*** MPI
         if (V3Nint) then
            call threebodysetup_cJ
         endif
c*** MPI
         call MPI_Barrier(icomm,ierr)
c*** MPI
         if (nucleonsi-nucleonsf>1) then
            call ist2_setup
            if (iproc==0) print *,' ist2_setup called, iproc=',
     $           iproc
            call MPI_Barrier(icomm,ierr)
         endif
         if (nucleonsi==nucleonsf+3) then
            call ist3_setup
            if (iproc==0) print *,' ist3_setup called, iproc=',
     $           iproc
            call MPI_Barrier(icomm,ierr)
         endif
      endif

      if (iabs(nucleonsi-nucleonsf)>=0.and.
     +     iabs(nucleonsi-nucleonsf)<=4) then
         if (twobdcal) then
            mnop=iabs(nucleonsi-nucleonsf)+2
         else
            mnop=iabs(nucleonsi-nucleonsf)
         endif            
         if (abs(nucleonsf-nucleonsi)==1.or.nucleonsf==nucleonsi-2
     +        .or.nucleonsf==nucleonsi-3.or.nucleonsf==nucleonsi-4
     +        .or.nucleonsi==nucleonsf) then

c            deallocate(ibascompf)
c            if (.not.(snsp)) then
c               if (associated(ibascompi)) deallocate(ibascompi)
c            endif
c               deallocate(ibascompi)
c            if (associated(ibascompf)) deallocate(ibascompf)
c            if (associated(ibascompi)) deallocate(ibascompi)

            if (antoine) then
               allocate(ilocf(nucleonsf,nsdf))
               do i1=1,nsdf
                  p_state=I_state_f(1,i1)
                  n_state=I_state_f(2,i1)
                  ilocf(1:nprotonsf,i1)=occ_p_f(:,p_state)
                  do i2=1,nneutrnsf
                     ilocf(nprotonsf+i2,i1)=occ_n_f(i2,n_state)
     +                    +mxnwdi*nbit
                  end do
               end do
cc               open(64,file='occ_test_save.tmp',status='unknown',
cc     +              form='formatted')
cc               write(64,'(i4,2x,i10)') nucleonsf,nsdf
cc               do i1=1,nsdf
cc                  write(64,'(i10)') i1
cc                  write(64,'(10i6)') (ilocf(i2,i1),i2=1,nucleonsf)
cc               end do
cc               close(64)
            else
               if (iproc==0) then
                  print *,' mxspsi,mxspsf,mxnwdi,mxnwdf=',
     +                 mxspsi,mxspsf,mxnwdi,mxnwdf
               endif
               if (mxspsf/=mxnwdi*nbit) then
                  do i1=1,nsdf
                     do i2=1,nneutrnsf
                        ilocf(nprotonsf+i2,i1)=ilocf(nprotonsf+i2,i1)
     +                       -mxspsf+mxnwdi*nbit
                     end do
                  end do
               endif
               if (mxspsi/=mxnwdi*nbit) then
                  do i1=1,nsdi
                     do i2=1,nneutrnsi
                        iloci(nprotonsi+i2,i1)=iloci(nprotonsi+i2,i1)
     +                       -mxspsi+mxnwdi*nbit
                     end do
                  end do
               endif
            endif
            ty%nucl=>nucleonsf
            ty%prot=>nprotonsf
            ty%neut=>nneutrnsf
            ty%dim_st=>nsdf
            ty%occ=>ilocf
            ty%mxsp=>mxspsf
            ty%nwd=>mxnwdf
            call construct_hash_table
            allocate(occ_temp(nucleonsf))
            do i1=1,nsdf
               occ_temp(:)=ilocf(:,i1)
               call get_state_index(nucleonsf,occ_temp,1,ionb)
               if (i1/=ionb) then
                  print *,'*** error:'
cc               if (iproc==0) then
                  print *,' i=',i1,'    index=',ionb
cc                  print *, ilocf(:,i1)
cc               endif
                  stop
               endif
            end do
            deallocate(occ_temp)
cc            stop
cc10101       continue
            if (nucleonsi==nucleonsf) then
               call obdtbdcalc_hash(snsp,sndp,dnuc)
               print *,' obdtbdcalc_hash called: iproc=',iproc
            else
               if (NCSMC_kernels) then
                  call NCSMC_coupling_kernels
                  call MPI_Barrier(icomm,ierr)
                  if (iproc==0) print *,
     $                 ' NCSMC_coupling_kernels called: nproc=',nproc
                  if (nucleonsi==nucleonsf+1) then
                     call NCSMC_coupling_kernels_phys_trans
                     call MPI_Barrier(icomm,ierr)
                     if (iproc==0) print *,
     $           ' NCSMC_coupling_kernels_phys_trans called: nproc=',
     $                 nproc
                  endif
                  if (iproc==0) then
                     call NCSMC_coupling_kernels_write
                     print *,' NCSMC_coupling_kernels_write called'
                  endif
               else
                  call cluster_overlap_calc_hash
                  print *,' cluster_overlap_calc_hash called: iproc=',
     $                 iproc
               endif
            endif
         else
            print *,'***error in input****'
            stop
         endif
      endif   

      call MPI_Barrier(icomm,ierr)
      if (iproc==0) close(iunitvout)
      if (nblock>1.or.NCSMC_kernels) then
         call MPI_Finalize(ierr)
         stop
      endif

      if (associated(nsdf)) deallocate(nsdf)  
      if (associated(mxnwdf)) deallocate(mxnwdf)  
      if (allocated(ndiff)) deallocate(ndiff)
      if (associated(nhomf)) deallocate(nhomf)  
      if (associated(mjtotalf)) deallocate(mjtotalf)  
      if (associated(mttotalf)) deallocate(mttotalf)  
      if (associated(iparityf)) deallocate(iparityf)  
      if (associated(nhwf)) deallocate(nhwf)  
      if (associated(nucleonsf)) deallocate(nucleonsf)  
      if (associated(nprotonsf)) deallocate(nprotonsf)  
      if (associated(nneutrnsf)) deallocate(nneutrnsf)  
      if (associated(hbof)) deallocate(hbof)  
      if (associated(nhom12f)) deallocate(nhom12f)  
      if (associated(nhom123f)) deallocate(nhom123f)  
      if (associated(nhom1234f)) deallocate(nhom1234f)  
      if (associated(naspsf)) deallocate(naspsf)  
      if (associated(mxspsf)) deallocate(mxspsf)  
      if (associated(nshllf)) deallocate(nshllf)  
      if (associated(mconff)) deallocate(mconff)  
      if (associated(nconff)) deallocate(nconff)  
      if (associated(nprtnf)) deallocate(nprtnf)  
      if (associated(nneutf)) deallocate(nneutf)  
      if (associated(iendconff)) deallocate(iendconff)  

      print *,' 1st block deallocated'

      if (.not.(snsp)) then
         if (associated(nsdi)) deallocate(nsdi)  

c         print *,' nsdii deallocated'

         if (associated(mxnwdi)) deallocate(mxnwdi)  

c         print *,' mxnwdi deallocated'

c         if (associated(nhomi)) deallocate(nhomi)  
c         print *,' nhomi deallocated'

         if (associated(mjtotali)) deallocate(mjtotali)  

c         print *,' mjtotali deallocated'

         if (associated(mttotali)) deallocate(mttotali)  

c         print *,' mttotali deallocated'

         if (associated(iparityi)) deallocate(iparityi)  

c         print *,' iparityi deallocated'

         if (associated(nhwi)) deallocate(nhwi)

c         print *,' nhwi deallocated'
  
         if (associated(nucleonsi)) deallocate(nucleonsi)  
         if (associated(nprotonsi)) deallocate(nprotonsi)  
         if (associated(nneutrnsi)) deallocate(nneutrnsi)  
         if (associated(hboi)) deallocate(hboi)

c         print *,' hboi deallocated'
  
c         if (associated(nhom12i)) deallocate(nhom12i)  
         if (associated(nhom123i)) deallocate(nhom123i)  
         if (associated(nhom1234i)) deallocate(nhom1234i)  
         if (associated(naspsi)) deallocate(naspsi)  
         if (associated(mxspsi)) deallocate(mxspsi)  
         if (associated(nshlli)) deallocate(nshlli)  
         if (associated(mconfi)) deallocate(mconfi)  
         if (associated(nconfi)) deallocate(nconfi)  
         if (associated(nprtni)) deallocate(nprtni)  
         if (associated(nneuti)) deallocate(nneuti)  
         if (associated(iendconfi)) deallocate(iendconfi)  
      endif

      print *,' 2nd block deallocated'

      deallocate(enerf)      
      if (.not.(snsp)) then
         if (associated(eneri)) deallocate(eneri)
      endif
      deallocate(jt2f,it2f,Jxf,Txf)
      if (.not.(snsp)) then
         if (associated(jt2i)) deallocate(jt2i)
         if (associated(it2i)) deallocate(it2i)
         if (associated(Jxi)) deallocate(Jxi)
         if (associated(Txi)) deallocate(Txi)
      endif

      print *,' 3rd block deallocated'

      deallocate(n_spf,l_spf,j2_spf,m2_spf,mt2_spf)
      if (.not.(snsp)) then
         if (associated(n_spi)) deallocate(n_spi)
         if (associated(l_spi)) deallocate(l_spi)
         if (associated(j2_spi)) deallocate(j2_spi)
         if (associated(m2_spi)) deallocate(m2_spi)
         if (associated(mt2_spi)) deallocate(mt2_spi)
      endif

      print *,' 4th block deallocated'

      if (associated(ibasisf)) deallocate(ibasisf)
      if (associated(bmpf)) deallocate(bmpf)
      if (associated(ilocf)) deallocate(ilocf)
      if (associated(I_state_f)) deallocate(I_state_f)
      if (associated(occ_p_f)) deallocate(occ_p_f)
      if (associated(occ_n_f)) deallocate(occ_n_f)
      if (associated(nsd_p_f)) deallocate(nsd_p_f)
      if (associated(nsd_n_f)) deallocate(nsd_n_f)
      if (associated(ibas_p_f)) deallocate(ibas_p_f)
      if (associated(ibas_n_f)) deallocate(ibas_n_f)

      print *,' 4.5th block deallocated'

      if (.not.(snsp)) then
         if (associated(ibasisi)) deallocate(ibasisi)
         print *,' ibasisi deallocated'
c         if (associated(bmpi)) deallocate(bmpi)
c         print *,' bmpii deallocated'
c         if (associated(iloci)) deallocate(iloci)
         print *,' iloci deallocated'
         if (associated(I_state_i)) deallocate(I_state_i)
         if (associated(occ_p_i)) deallocate(occ_p_i)
         if (associated(occ_n_i)) deallocate(occ_n_i)
         if (associated(nsd_p_i)) deallocate(nsd_p_i)
         if (associated(nsd_n_i)) deallocate(nsd_n_i)
         if (associated(ibas_p_i)) deallocate(ibas_p_i)
         if (associated(ibas_n_i)) deallocate(ibas_n_i)
      endif

      print *,' 5th block deallocated'

c      deallocate(ibascompf)
c      if (.not.(snsp)) then
c      if (associated(ibascompi)) deallocate(ibascompi)
c      endif

      if (allocated(tdJp)) deallocate(tdJp)
      if (allocated(tdJn)) deallocate(tdJn)
      if (allocated(ad_cl)) deallocate(ad_cl)
      if (allocated(ad_cl_2)) deallocate(ad_cl_2)
      deallocate(ist1,iobsind)
      if (allocated(t2bd)) deallocate(t2bd)
      if (allocated(tbd)) deallocate(tbd)
      if (allocated(ist2)) deallocate(ist2)
      if (allocated(ist2_Jst)) deallocate(ist2_Jst)
      if (allocated(itbind)) deallocate(itbind)
      if (allocated(ist3)) deallocate(ist3)
      if (allocated(ad3ind)) deallocate(ad3ind)
      if (allocated(ist4)) deallocate(ist4)
      if (allocated(ad4ind)) deallocate(ad4ind)

      print *,' 6th block deallocated'

c**********************************************
c*** test read in *** 
33333 continue
      if (iproc==0) then
         if (.not.cluster) then
            call testread(snsp,sndp,dnuc)
            if (.not.snsp) formfcal=.false.
            if (snsp.and.formfcal) call obformf
            if (radial) call radialdist
            if (radial) call radialdist_trinv
         else
            call cluster_overlap_testread
            if (abs(nucleonsi-nucleonsf)==1.and.momdist) call momdistr
         endif
         close(iunitobout)
      endif   

      call MPI_Barrier(icomm,ierr)
      call MPI_Finalize(ierr)
      stop
      end


      subroutine confwavread(nsd,mxnwd,nhom,mjtotal,mttotal,iparity,nhw,
     +        nucleons,nprotons,nneutrns,jt2,it2,Jx,Tx,hbo,iloc,nshll,
     +        nhom12,nhom123,nhom1234,nasps,mxsps,mconf,nconf,
     +        nprtn,nneut,iendconf,
     +        n_sp,l_sp,j2_sp,m2_sp,mt2_sp,bmp,ibasis,ener)
**** 9-August-1999 *****
**** reads in configuration-space many-body wave function ****
**** computed using the many-fermion dynamics code ****
**** Petr Navratil, University of Arizona *****
      use intrface
      use paramdef
      use nodeinfo
      implicit none
      include 'mpif.h'
      integer(2),pointer :: iloc(:,:),mconf(:)
      integer(4),pointer :: nconf,nprtn(:,:),nneut(:,:),iendconf(:)
      integer(4),pointer:: jt2(:),it2(:)
      integer(4),pointer :: n_sp(:),l_sp(:),
     +     j2_sp(:),m2_sp(:),mt2_sp(:)
      real(kind=kind(0.0d0)),pointer::bmp(:,:)
      integer(4),pointer :: ibasis(:,:,:)
      integer(4),pointer :: nsd,mxnwd,nhom,nhom12,nasps,mxsps,nshll,
     +     nhom123,nhom1234
      integer(4),pointer :: mjtotal,mttotal,iparity,nhw
      integer(4),pointer :: nucleons,nprotons,nneutrns
      real(kind=kind(0.0d0)),pointer :: Jx(:),Tx(:),ener(:)
      real(kind=kind(0.0d0)),pointer :: hbo

      integer(4) major,i,ia,ierr,iixxx,ii,ix_ii,ix_iii,ixxx
      integer(4) mxsps2
      integer(4) nhme,k1max,k1,i1,kx,iy,ix,iii,icls,locx,iwd,ibit
      real(4) hboin
      real(4) xj,xt
c*      real(8) xj,xt
c*      integer(2) mcon  
c*      integer(2),allocatable:: mconf(:)
      integer(4),allocatable:: ntemp(:)
      real(kind=kind(0.d0)),allocatable :: Jxin(:),Txin(:)
      real(kind=kind(0.d0)),allocatable :: enerin(:)
c*      real(kind=kind(0.0)),allocatable :: enerin(:)
      real(kind=kind(0.d0)),allocatable :: bmp0(:)
      integer(2),allocatable :: occ_tmp(:)
      integer(4),allocatable :: max_in(:),max_out(:)

*****************

      allocate(nsd,mxnwd,nhom,mjtotal,mttotal,iparity,nhw,nshll,nconf,
     +     nucleons,nprotons,nneutrns,hbo,nhom12,nhom123,nhom1234,
     +     nasps,mxsps)

      if (iproc==0) print *,' nsd...hbo allocated'

      open(iunitvec,file=filename,form='unformatted',
     + status='old',action='read')

      if (iproc==0) print *,' unit iunitvec opened'

      read(iunitvec,end=333,err=333) 
     + nsd,nhme,k1max,nucleons,nprotons,nneutrns,
     + hboin,nhw

      if (iproc==0) print *,' nsd...hboin,nhw read'

      if (kxxx>k1max+2) then
         if (iproc==0) then
            write(iunitvout,*)'**** Error: too large kxxx=',kxxx
         endif   
         stop
      endif   

      hbo=dble(hboin)

      read(iunitvec,end=333,err=333) 
     + mjtotal,mttotal,nshll,mxnwd,mxsps,iparity,major

      if (iproc==0) print *,' mjtotal...iparity,major read'
      if (iproc==0) print *,' mxsps=',mxsps

      if (mxsps/=32*mxnwd.and.mxsps/=64*mxnwd) then
         if (iproc==0) then
            write(iunitvout,*)'**** Error: mxsps,mxnwd=',mxsps,mxnwd
         endif
         stop
      endif
      mxsps2=2*mxsps
*      if (mxsps==32*mxnwd.and.nbit==64) then
*         mxnwd=mxnwd/2
*         if (nbit*mxnwd*2/=mxsps2) then
*            if (iproc==0) then
*               write(iunitvout,*)'**** Problem: mxsps,mxnwd,nbit=',
*     +              mxsps,mxnwd,nbit
*            endif
*         endif
*      endif 
   

      allocate(mconf(nsd))
      read(iunitvec,end=333,err=333) (mconf(i1),i1=1,nsd)
c***      read(iunitvec,end=333,err=333) (mcon,i1=1,nsd)

      if (iproc==0) print *,' mcon read'
      if (iproc==0) print *,' kxxx,nkxxx,kxxx+nkxxx-1:',
     $     kxxx,nkxxx,kxxx+nkxxx-1
      if (iproc==0) print *,' k1max=',k1max

      if (kxxx+nkxxx-1>k1max+2) nkxxx=k1max+3-kxxx
      allocate(bmp(nsd,kxxx:kxxx+nkxxx-1))
      allocate(bmp0(nsd))
      read(iunitvec,end=333,err=333) (bmp0(i1),i1=1,nsd)
      if (select) then
         iixxx=1
         if (iselxxx(1)==1) then
            bmp(:,1)=bmp0(:)
            if (iixxx<nkxxx) iixxx=iixxx+1
         endif
         do k1=2,k1max+2
            read(iunitvec,end=333,err=333) (bmp0(i1),i1=1,nsd)
            if (k1==iselxxx(iixxx)) then
               bmp(:,kxxx+iixxx-1)=bmp0(:)
               if (iixxx<nkxxx) iixxx=iixxx+1
            endif
         end do
      else
         if (kxxx==1) bmp(:,1)=bmp0(:)
         do k1=2,k1max+2
            read(iunitvec,end=333,err=333) (bmp0(i1),i1=1,nsd)
            if (k1>=kxxx.and.k1<=kxxx+nkxxx-1) bmp(:,k1)=bmp0(:)
         end do
      endif

      deallocate(bmp0)
      if (iproc==0) print *,' bmp read'
      if (iproc==0) print *,' mxsps2=',mxsps2

      allocate(n_sp(mxsps2),l_sp(mxsps2),j2_sp(mxsps2),
     +     m2_sp(mxsps2),mt2_sp(mxsps2))

      if (iproc==0) print *,' n_sp..mt2_sp allocated'

      read(iunitvec,end=333,err=333) 
     +(n_sp(i),l_sp(i),j2_sp(i),m2_sp(i),mt2_sp(i),i=1,mxsps2)

      if (iproc==0) print *,' nsp read'
      if (iproc==0) print *,' major=',major
      
      if (major>-1) then
         allocate(iloc(nucleons,nsd)) 
         do ia=1,nucleons
            read(iunitvec,end=333,err=333) 
c     +           (iloc(ia,i1),i1=1,nsd)
     +           (mconf(i1),i1=1,nsd)
            do i1=1,nsd
               iloc(ia,i1)=mconf(i1)
            end do
         end do
         deallocate(mconf)
         allocate(occ_tmp(nucleons))
         nasps=0
         nhom=0
         nhom12=0
         nhom123=0
         nhom1234=0
         do i1=iproc+1,nsd,nproc
            occ_tmp(:)=iloc(:,i1)
            do ia=1,nucleons
               iy=occ_tmp(ia) !iloc(ia,i1)
               if (iy>mxsps) iy=iy-mxsps
               if (iy>nasps) nasps=iy
               if (nhom<2*n_sp(iy)+l_sp(iy))nhom=2*n_sp(iy)+l_sp(iy)
               do i=ia+1,nucleons
                  ix=occ_tmp(i) !iloc(i,i1)
                  if(nhom12<2*n_sp(ix)+l_sp(ix)+2*n_sp(iy)+l_sp(iy))
     +                   nhom12=2*n_sp(ix)+l_sp(ix)+2*n_sp(iy)+l_sp(iy)
                  do ii=i+1,nucleons
                     ix_ii=occ_tmp(ii) !iloc(ii,i1)
                     if(nhom123<2*n_sp(ix_ii)+l_sp(ix_ii)
     +                    +2*n_sp(ix)+l_sp(ix)+2*n_sp(iy)+l_sp(iy))
     +                    nhom123=2*n_sp(ix_ii)+l_sp(ix_ii)
     +                    +2*n_sp(ix)+l_sp(ix)+2*n_sp(iy)+l_sp(iy)
                     do iii=ii+1,nucleons
                        ix_iii=occ_tmp(iii) !iloc(iii,i1)
                        if(nhom1234<2*n_sp(ix_iii)+l_sp(ix_iii)
     +                       +2*n_sp(ix_ii)+l_sp(ix_ii)
     +                       +2*n_sp(ix)+l_sp(ix)+2*n_sp(iy)+l_sp(iy))
     +                       nhom1234=2*n_sp(ix_iii)+l_sp(ix_iii)
     +                       +2*n_sp(ix_ii)+l_sp(ix_ii)
     +                       +2*n_sp(ix)+l_sp(ix)+2*n_sp(iy)+l_sp(iy)
                     end do
                  end do
               end do
            end do   
         end do   
         deallocate(occ_tmp)
         allocate(max_in(5),max_out(5))
         max_in(1)=nasps
         max_in(2)=nhom
         max_in(3)=nhom12
         max_in(4)=nhom123
         max_in(5)=nhom1234
         call MPI_Allreduce(max_in,max_out,5,MPI_INTEGER,
     $        MPI_MAX,icomm,ierr)
         nasps=max_out(1)
         nhom=max_out(2)
         nhom12=max_out(3)
         nhom123=max_out(4)
         nhom1234=max_out(5)
         deallocate(max_in,max_out)
         mxnwd=nasps/nbit
         if (mxnwd*nbit<nasps) mxnwd=mxnwd+1
         if (mxnwd*nbit/=mxsps) then
            allocate(ntemp(2*mxnwd*nbit))
            ntemp=0
            do ia=1,min(mxnwd*nbit,mxsps)
               ntemp(ia)=n_sp(ia)
               ntemp(mxnwd*nbit+ia)=n_sp(mxsps+ia)
            end do
            deallocate(n_sp)
            allocate(n_sp(2*mxnwd*nbit))
            n_sp=ntemp
            ntemp=0
            do ia=1,min(mxnwd*nbit,mxsps)
               ntemp(ia)=l_sp(ia)
               ntemp(mxnwd*nbit+ia)=l_sp(mxsps+ia)
            end do
            deallocate(l_sp)
            allocate(l_sp(2*mxnwd*nbit))
            l_sp=ntemp
            ntemp=0
            do ia=1,min(mxnwd*nbit,mxsps)
               ntemp(ia)=j2_sp(ia)
               ntemp(mxnwd*nbit+ia)=j2_sp(mxsps+ia)
            end do
            deallocate(j2_sp)
            allocate(j2_sp(2*mxnwd*nbit))
            j2_sp=ntemp
            ntemp=0
            do ia=1,min(mxnwd*nbit,mxsps)
               ntemp(ia)=m2_sp(ia)
               ntemp(mxnwd*nbit+ia)=m2_sp(mxsps+ia)
            end do
            deallocate(m2_sp)
            allocate(m2_sp(2*mxnwd*nbit))
            m2_sp=ntemp
            ntemp=0
            do ia=1,min(mxnwd*nbit,mxsps)
               ntemp(ia)=mt2_sp(ia)
               ntemp(mxnwd*nbit+ia)=mt2_sp(mxsps+ia)
            end do
            deallocate(mt2_sp)
            allocate(mt2_sp(2*mxnwd*nbit))
            mt2_sp=ntemp
            deallocate(ntemp)
            do ia=1,nucleons
               do i1=1,nsd
                  iy=iloc(ia,i1)
                  if (iy>mxsps)  then
                     iy=iy-mxsps+mxnwd*nbit
                     iloc(ia,i1)=iy
                  endif   
               end do
            end do   
            mxsps=mxnwd*nbit
         endif           
      else
         if (mxsps/=mxnwd*nbit) then
            if (iproc==0) then
               write(iunitvout,*)'**** Error: mxsps,mxnwd,nbit=',
     +              mxsps,mxnwd,nbit
            endif
            call MPI_Abort(icomm,1011,ierr)
            stop
         endif            
         nhom=0
         do i=1,2*mxnwd*nbit
            if (nhom<2*n_sp(i)+l_sp(i)) nhom=2*n_sp(i)+l_sp(i)
         end do
         nhom12=2*nhom         
         nasps=mxsps
         allocate(ibasis(mxnwd,2,nsd))
         do iwd=1,mxnwd
            do icls=1,2
               read(iunitvec,end=333,err=333) 
     +              (ibasis(iwd,icls,i1),i1=1,nsd)  
            end do
         end do 
      endif 

      if (iproc==0) then 
         print *,' nhom12,nhom123,nhom1234=',nhom12,nhom123,nhom1234
         print *,' iloc read'
         print *,' mxsps=',mxsps
      endif

      allocate(enerin(k1max+2))      
      read(iunitvec,end=333,err=333)(enerin(k1),k1=1,k1max+2)

      if (iproc==0) print *,' ener read'

      allocate(Jxin(k1max+2),Txin(k1max+2))

      do k1 = 1, k1max+2
         if (iproc==0) print *,' k1=',k1
         read(iunitvec,end=333,err=333) kx,xj,xt
c*         read(iunitvec,end=333,err=333) kx,xj,xt,enerin(k1)
         if (iproc==0) print *,' kx=',kx,' xj,xt:',xj,xt

         Jxin(k1)=xj
         Txin(k1)=xt 
      end do 
      if (iproc==0) print *,' Jx, Tx read'

c      read(iunitvec,end=333,err=333) nconf
c      print *,' nconf=',nconf
c      allocate(nprtn(nconf,nshll),nneut(nconf,nshll))
c      allocate(iendconf(nconf))
c      read(iunitvec,end=333,err=333) nprtn,nneut
c      read(iunitvec,end=333,err=333) iendconf

c      print *,' nprtn,nneut,iendconf read'

      close(iunitvec)

      if (select) then
         allocate(jt2(kxxx:kxxx+nkxxx-1))
         allocate(it2(kxxx:kxxx+nkxxx-1))
         allocate(Jx(kxxx:kxxx+nkxxx-1),Tx(kxxx:kxxx+nkxxx-1))
         allocate(ener(kxxx:kxxx+nkxxx-1))
         do iixxx=1,nkxxx
            ener(kxxx+iixxx-1)=enerin(iselxxx(iixxx))
            Jx(kxxx+iixxx-1)=Jxin(iselxxx(iixxx))
            Tx(kxxx+iixxx-1)=Txin(iselxxx(iixxx))
         end do
         k1max=nkxxx-2
      else
         allocate(jt2(k1max+2))
         allocate(it2(k1max+2))
         allocate(Jx(k1max+2),Tx(k1max+2))
         allocate(ener(k1max+2))
         ener=enerin
         Jx=Jxin
         Tx=Txin
      endif
      do ixxx=kxxx,kxxx+nkxxx-1
         if (mod(nucleons,2)==0) then
            jt2(ixxx)=2*nint(Jx(ixxx))
            it2(ixxx)=2*nint(Tx(ixxx))
         else
            jt2(ixxx)=nint(2*Jx(ixxx))
            if (mod(jt2(ixxx),2)/=1) then
               if (2*Jx(ixxx)>jt2(ixxx)) then
                  jt2(ixxx)=nint(2*Jx(ixxx))+1
               else
                  jt2(ixxx)=nint(2*Jx(ixxx))-1
               endif
            endif
            it2(ixxx)=nint(2*Tx(ixxx))
            if (mod(it2(ixxx),2)/=1) then
               if (2*Tx(ixxx)>it2(ixxx)) then
                  it2(ixxx)=nint(2*Tx(ixxx))+1
               else
                  it2(ixxx)=nint(2*Tx(ixxx))-1
               endif
            endif
         endif
         if (dabs(dble(jt2(ixxx))-2.d0*Jx(ixxx))>1.d-4) then
            if (iproc==0) print *, ixxx,'**** not a good J ****'
         endif 
      end do
      deallocate(enerin)
      deallocate(Jxin,Txin)

      if (iproc==0) then      
         write(iunitvout,*)
         write(iunitvout,*)'*** Nuclear states ***'

      if (iparity==0) then
         write(iunitvout,1000) nucleons,nprotons,nneutrns,mjtotal,
     + mttotal,hbo,nhw,nsd,nhme,k1max,
     + mxnwd,mxsps,major,iparity
 1000    format(' Nucleus:',/,' A=',i3,'   Z=',i3,'   N=',i3,
     + /,' 2*MJ=',i3,'   2*MT=',i3,'  parity= +',/,
     + ' hbar Omega=',f8.4,'   Nhw=',i3,'   dimension=',i8,'   nhme=',
     + i10,/,' k1max=',i3,'   mxnwd=',i3,'   mxsps=',i8,
     +                    '   major=',i2,'   iparity=',i2,/)
      elseif(iparity==1) then
         write(iunitvout,1007) nucleons,nprotons,nneutrns,mjtotal,
     + mttotal,hbo,nhw,nsd,nhme,k1max,
     + mxnwd,mxsps,major,iparity
 1007    format(' Nucleus:',/,' A=',i3,'   Z=',i3,'   N=',i3,
     + /,' 2*MJ=',i3,'   2*MT=',i3,'  parity= -',/,
     + ' hbar Omega=',f8.4,'   Nhw=',i3,'   dimension=',i8,'   nhme=',
     + i10,/,' k1max=',i3,'   mxnwd=',i3,'   mxsps=',i8,
     +                    '   major=',i2,'   iparity=',i2,/)
      else
         write(iunitvout,*) ' error: iparity=',iparity
         call MPI_Abort(icomm,1111,ierr)
         stop
      endif

      if (select) then
         write(iunitvout,1100) (Jx(k1),Tx(k1),ener(k1),
     +        ener(k1)-ener(kxxx),k1=kxxx,kxxx+nkxxx-1)
      else
         write(iunitvout,1100) (Jx(k1),Tx(k1),ener(k1),
     +        ener(k1)-ener(1),k1=1,k1max+2)
      endif
 1100 format(' J=',f7.4,'    T=',f7.4,'     Energy=',f12.4,
     + '     Ex=',f12.4) 
      endif

      if (iproc==0) write(iunitvout,1156) nhom,nhom12,nasps
 1156 format(/,' N1_max=',i4,'   N12_max=',i4,'   Nasps=',i4)

c      open(64,file='occ_test_save.tmp',status='unknown',
c     +     form='formatted')
c      write(64,'(i4,2x,i10)') nucleons,nsd
c      do i1=1,nsd
c         write(64,'(i10)') i1
c         write(64,'(10i6)') (iloc(ia,i1),ia=1,nucleons)
c      end do
c      close(64)
c*      do i1=1,1000
c*         write(iunitvout,3434) (iloc(ia,i1),ia=1,nucleons)
c* 3434    format(10i6)
c*      end do

*      write(iunitvout,*)' used state:' 
*      write(iunitvout,3435) (bmp(i1,kxxx),i1=1,nsd)
* 3435 format(10f8.4)
c** state basis setup *******************
      if (majortot<2) then
         allocate(ibasis(mxnwd,2,nsd))
         do iy=1,nsd
            ibasis(:,:,iy)=0
            do ix=1,nucleons
               locx=iloc(ix,iy)-1
               icls=locx/mxsps
               locx=locx-icls*mxsps 
               iii=icls+1            
               iwd=locx/nbit
               ibit=locx-iwd*nbit
               iwd=iwd+1
               ibasis(iwd,iii,iy)=ibset(ibasis(iwd,iii,iy),ibit)
            end do 
         end do   
         print *,' state basis calculated' 
      endif
c**********************************************
      return
 333  continue
      if (iproc==0) write(iunitvout,*) '**** error in reading ****'
      call MPI_Abort(icomm,3331,ierr)
      stop
      end


      subroutine confwave_antoine(nsd,mxnwd,nhom,mjtotal,mttotal,
     $     iparity,nhw,
     +     nucleons,nprotons,nneutrns,jt2,Jx,
     +     I_state,nsd_p,nsd_n,occ_p,occ_n,ibas_p,ibas_n,nshll,
     +     nhom12,nhom123,nhom1234,nasps,mxsps,mconf,nconf,
     +     nprtn,nneut,iendconf,
     +     n_sp,l_sp,j2_sp,m2_sp,mt2_sp,bmp,ener,i_case)
**** 25-May-2004 *****
**** reads in configuration-space many-body wave function ****
**** computed using the Antoine code ****
**** Petr Navratil, LLNL *****
      use intrface
      use paramdef
      use nodeinfo
      use initst, only: hboi,it2i,Txi,num_of_in_i
      use finast, only: it2f,Txf,num_of_in_f
      implicit none
      include 'mpif.h'
      integer(2),pointer :: mconf(:)
      integer(4),pointer :: I_state(:,:),nsd_p,nsd_n,
     +     occ_p(:,:),occ_n(:,:),ibas_p(:,:),ibas_n(:,:)
      integer(4),pointer :: nconf,nprtn(:,:),nneut(:,:),iendconf(:)
      integer(4),pointer:: jt2(:)
      integer(4),pointer :: n_sp(:),l_sp(:),
     +     j2_sp(:),m2_sp(:),mt2_sp(:)
      real(kind=kind(0.0d0)),pointer::bmp(:,:)
      integer(4),pointer :: nsd,mxnwd,nhom,nhom12,nasps,mxsps,nshll,
     +     nhom123,nhom1234
      integer(4),pointer :: mjtotal,mttotal,iparity,nhw
      integer(4),pointer :: nucleons,nprotons,nneutrns
      real(kind=kind(0.0d0)),pointer :: Jx(:),ener(:)
      integer :: i_case

      integer(4) :: major,i,ia,ierr,iixxx,ii,ix_ii,ix_iii,ixxx
      integer(4) :: mxsps2
      integer(4) :: nhme,k1max,k1,i1,kx,iy,ix,iii,locx,iwd,ibit
      integer(4),allocatable:: ntemp(:)
      integer :: num_of_shell,num_of_jm,A(0:2),nslt(2),fon(14),
     $     coul,iprec
      integer,allocatable :: nr(:),ll(:),jj(:),num(:),mpr(:),
     $     occ_1(:,:),occ_2(:,:)
      real(kind=kind(0.d0)) :: en
      integer :: p_state,n_state,nhomp,nhomn,nhom12p,nhom12n,nhom123p,
     +     nhom123n,nhom1234p,nhom1234n,diff
      integer,allocatable :: mxpi(:),mnpi(:),mxnu(:),mnnu(:),
     +     mxpn(:),mnpn(:),ntempi(:),ntemnu(:),npib(:),nnub(:),
     $     conf_p(:,:),conf_n(:,:)
*****************

      allocate(nsd,nsd_p,nsd_n,mxnwd,nhom,mjtotal,mttotal,iparity,nhw,
     +     nshll,nconf,
     +     nucleons,nprotons,nneutrns,nhom12,nhom123,nhom1234,
     +     nasps,mxsps)

      print *,' nsd... allocated in confwave_antoine'

      open(iunitvec,file=filename,form='unformatted',
     +     status='old',action='read')

c      open(33,file=filename//'_form',form='formatted',
c     +     status='unknown',action='write')
c      open(33,file=filename//'_form',form='formatted',
c     +     status='old',action='read')

      print *,' unit iunitvec opened'

      read(iunitvec,end=333,err=333) num_of_shell,num_of_jm,nsd,
     $     A(0),A(1),A(2),nslt(1),nslt(2) !new
c     $     A(1),A(2),nslt(1),nslt(2)  !old

c      write(33,'(2i6,i15,5i10)') num_of_shell,num_of_jm,nsd,
ccc     $     A(0),A(1),A(2),nslt(1),nslt(2)
c     $     A(1),A(2),nslt(1),nslt(2)
c      read(33,'(2i6,i15,5i10)') num_of_shell,num_of_jm,nsd,
ccc     $     A(0),A(1),A(2),nslt(1),nslt(2)
c     $     A(1),A(2),nslt(1),nslt(2)

      print *,' num_of_shell,num_of_jm,nsd,A,nslt read:',
     $     num_of_shell,num_of_jm,nsd,A(0),A(1),A(2),nslt(1),nslt(2)

      allocate(nr(num_of_shell),ll(num_of_shell),jj(num_of_shell))

      read(iunitvec,end=333,err=333)
     $     (nr(i),ll(i),jj(i),i=1,num_of_shell)

c      write(33,'(100(3i4))') (nr(i),ll(i),jj(i),i=1,num_of_shell)
c      read(33,'(100(3i4))') (nr(i),ll(i),jj(i),i=1,num_of_shell)

      print *,' nr,ll,jj read'

*      do i=1,num_of_shell
*         print *, i,nr(i),ll(i),jj(i)
*      end do

      allocate(num(num_of_jm),mpr(num_of_jm))

      read(iunitvec,end=333,err=333) (num(i),mpr(i),i=1,num_of_jm)

c      write(33,'(100(2i4))') (num(i),mpr(i),i=1,num_of_jm)
c      read(33,'(100(2i4))') (num(i),mpr(i),i=1,num_of_jm)

      print *,' num,mpr read'

*      do i=1,num_of_jm
*         print *, i,num(i),mpr(i)
*      end do

      mxsps=num_of_jm
      print *,' num_of_jm=',num_of_jm
      mxsps2=2*mxsps
      print *,' mxsps2=',mxsps2
      allocate(n_sp(mxsps2),l_sp(mxsps2),j2_sp(mxsps2),
     +     m2_sp(mxsps2),mt2_sp(mxsps2))
      print *,' n_sp..mt2_sp allocated'
      do i=1,mxsps
         n_sp(i)=nr(num(i))
         n_sp(i+mxsps)=n_sp(i)
         l_sp(i)=ll(num(i))
         l_sp(i+mxsps)=l_sp(i)
         j2_sp(i)=jj(num(i))
         j2_sp(i+mxsps)=j2_sp(i)
         m2_sp(i)=mpr(i)
         m2_sp(i+mxsps)=m2_sp(i)
         mt2_sp(i)=1
         mt2_sp(i+mxsps)=-1
      end do
      print *,' n_sp... set'
      deallocate(num,mpr)
      deallocate(nr,ll,jj)
      print *,' num..nr.. deallocated'
      
      allocate(I_state(2,nsd))
      print *,' I_state allocated'

      read(iunitvec,end=333,err=333) ((I_state(ii,i),ii=1,2),i=1,nsd)

c      write(33,'(100(2i10))') ((I_state(ii,i),ii=1,2),i=1,nsd)
c      read(33,'(100(2i10))') ((I_state(ii,i),ii=1,2),i=1,nsd)

      print *,' I_state read'

*      do i=1,nsd
*         print *, i,I_state(1,i),I_state(2,i)
*      end do

      allocate(occ_1(A(1),nslt(1)),occ_2(A(2),nslt(2)))

      read(iunitvec,end=333,err=333)
     $     ((occ_1(ii,i),ii=1,A(1)),i=1,nslt(1))

c      write(33,'(100i5)') ((occ_1(ii,i),ii=1,A(1)),i=1,nslt(1))
c      read(33,'(100i5)') ((occ_1(ii,i),ii=1,A(1)),i=1,nslt(1))

**      print *,' occ_1:',occ_1
      print *,' occ_1 read'

*      do i=1,nslt(1)
*         do ii=1,A(1)
*            print *,i,ii,occ_1(ii,i)
*         end do
*      end do

      read(iunitvec,end=333,err=333)
     $     ((occ_2(ii,i),ii=1,A(2)),i=1,nslt(2))

c      write(33,'(100i5)') ((occ_2(ii,i),ii=1,A(2)),i=1,nslt(2))
c      read(33,'(100i5)') ((occ_2(ii,i),ii=1,A(2)),i=1,nslt(2))

**      print *,' occ_2:',occ_2
      print *,' occ_2 read'

*      do i=1,nslt(1)
*         do ii=1,A(1)
*            print *,i,ii,occ_1(ii,i)
*         end do
*      end do

      if (i_case==1) then
         k1max=num_of_in_i-2
      else
         k1max=num_of_in_f-2
      endif
c         k1max=nkxxx-2+kxxx-1
         
      print *,' kxxx,nkxxx,kxxx+nkxxx-1:',kxxx,nkxxx,kxxx+nkxxx-1
      print *,' k1max=',k1max

      if (kxxx+nkxxx-1>k1max+2) nkxxx=k1max+3-kxxx
      allocate(bmp(nsd,kxxx:kxxx+nkxxx-1))
c      allocate(jt2(kxxx:kxxx+nkxxx-1))
c      allocate(Jx(kxxx:kxxx+nkxxx-1))
c      allocate(ener(kxxx:kxxx+nkxxx-1))
      allocate(jt2(k1max+2))
      allocate(Jx(k1max+2))
      allocate(ener(k1max+2))

      read(iunitvec,end=333,err=333) fon,en   !new

c      read(iunitvec,end=333,err=333) fon(1:13),en  !old

c      write(33,'(13i16,e14.7)') fon(1:13),en
c      read(33,'(13i16,e14.7)') fon(1:13),en

c      fon(14)=2 !old
      print *,' en=',en
      print *,' fon=',fon

      if (select) then
         k1=1
         iixxx=1
         if (iselxxx(1)==1) then
            ener(1)=en
            jt2(1)=fon(10)
            Jx(1)=real(jt2(1),kind(0.d0))/2.d0
            coul=fon(12)
            iprec=fon(14)
            mjtotal=fon(7)
            iparity=fon(8)
            print *,' k1,iixxx,iprec=',k1,iixxx,iprec
            call read_vector(.true.,iprec,nsd,1)
            if (iixxx<nkxxx) iixxx=iixxx+1
         else
            ener(1)=en
            jt2(1)=fon(10)
            Jx(1)=real(jt2(1),kind(0.d0))/2.d0
            iprec=fon(14)
            call read_vector(.false.,iprec,nsd,0)
         endif
         do 
            k1=k1+1

            read(iunitvec,end=366,err=366) fon,en

c            write(33,'(13i16,e14.7)') fon(1:13),en
c            read(33,'(13i16,e14.7)') fon(1:13),en

            print *,' en=',en
            print *,' fon=',fon
            if (k1==iselxxx(iixxx)) then
               ener(k1)=en
               jt2(k1)=fon(10)
               Jx(k1)=real(jt2(k1),kind(0.d0))/2.d0
               coul=fon(12)
               iprec=fon(14)
               mjtotal=fon(7)
               iparity=fon(8)
               print *,' k1,iixxx,iprec=',k1,iixxx,iprec
               call read_vector(.true.,iprec,nsd,kxxx+iixxx-1)
               if (iixxx<nkxxx) iixxx=iixxx+1
            else
               ener(k1)=en
               jt2(k1)=fon(10)
               Jx(k1)=real(jt2(k1),kind(0.d0))/2.d0
               iprec=fon(14)
               call read_vector(.false.,iprec,nsd,0)
            endif
         end do
      else
         if (kxxx==1) then
            ener(1)=en
            jt2(1)=fon(10) !new
c            jt2(1)=fon(9)   !old
            Jx(1)=real(jt2(1),kind(0.d0))/2.d0
            coul=fon(12) !new
            iprec=fon(14) 
            print *, ' iprec=',iprec
            mjtotal=fon(7) !new
            iparity=fon(8) !new
c            coul=fon(11)    !old
c            mjtotal=fon(6)  !old
c            iparity=fon(7)  !old
            print *,' en=',ener(1),' J=',jt2(1),' coul=',coul
            print *,' M=',mjtotal,'  pi=',iparity
            call read_vector(.true.,iprec,nsd,1)
*            print *,' bmp(1)=',bmp(:,1)
         else
            ener(1)=en
            jt2(1)=fon(10)
            Jx(1)=real(jt2(1),kind(0.d0))/2.d0
            iprec=fon(14)
            call read_vector(.false.,iprec,nsd,0)
         endif
         do k1=2,k1max+2
            read(iunitvec,end=333,err=333) fon,en  !new

c           read(iunitvec,end=333,err=333) fon(1:13),en !old

c            write(33,'(13i16,e14.7)') fon(1:13),en
c            read(33,'(13i16,e14.7)') fon(1:13),en

            print *,' en=',en
            print *,' fon=',fon
            if (k1>=kxxx.and.k1<=kxxx+nkxxx-1) then
               ener(k1)=en
               jt2(k1)=fon(10) !new
c               jt2(k1)=fon(9)    !old
               Jx(k1)=real(jt2(k1),kind(0.d0))/2.d0
               coul=fon(12) !new
               iprec=fon(14) 
               print *, ' iprec=',iprec
               mjtotal=fon(7) !new
               iparity=fon(8) !new
c               coul=fon(11)     !old
c               mjtotal=fon(6)   !old
c               iparity=fon(7)   !old
               print *,' en=',ener(k1),' J=',jt2(k1),' coul=',coul
               print *,' M=',mjtotal,'  pi=',iparity
               call read_vector(.true.,iprec,nsd,k1)
*               print *,' k1,bmp(k1)=',k1,bmp(:,k1)
            else
               ener(k1)=en
               jt2(k1)=fon(10)
               Jx(k1)=real(jt2(k1),kind(0.d0))/2.d0
               iprec=fon(14)
               call read_vector(.false.,iprec,nsd,0)
            endif
         end do
      endif
 366  continue
      close(iunitvec)
c      close(33)
      print *,' bmp read'
**      print *, bmp(:,1)
c***
c      coul=2   !old
c***
      print *,' coul=',coul
      if (coul==1) then
         nprotons=A(1)
         nneutrns=A(2)
         nsd_p=nslt(1)
         nsd_n=nslt(2)
         allocate(occ_p(nprotons,nsd_p))
         allocate(occ_n(nneutrns,nsd_n))
         occ_p=occ_1
         occ_n=occ_2
      else
         nprotons=A(2)
         nneutrns=A(1)
         nsd_p=nslt(2)
         nsd_n=nslt(1)
         allocate(occ_p(nprotons,nsd_p))
         allocate(occ_n(nneutrns,nsd_n))
         occ_p=occ_2
         occ_n=occ_1
         do i1=1,nsd
            ix=I_state(1,i1)
            I_state(1,i1)=I_state(2,i1)
            I_state(2,i1)=ix
         end do
      endif
      deallocate(occ_1,occ_2)
      mttotal=nprotons-nneutrns
      nucleons=nprotons+nneutrns

      major=3
      if (nprotons<=40.and.nneutrns<=40) then
         nhw=fon(1)+max(0,min(6,nprotons-2))+max(0,min(6,nneutrns-2))
     +        +2*max(0,min(12,nprotons-8))+2*max(0,min(12,nneutrns-8))
     +        +3*max(0,min(20,nprotons-20))+3*max(0,min(20,nneutrns-20))
         print *,' fon(1),nhw=',fon(1),nhw
      else
         print *,' error: protons,neutrons=',nprotons,nneutrns
         stop
      endif
      nhme=0

      nasps=0
      nhom=0
      nhom12=0
      nhom123=0
      nhom1234=0
      do i1=1,nsd
         p_state=I_state(1,i1)
         n_state=I_state(2,i1)
         nhomp=0
         nhomn=0
         nhom12p=0
         nhom12n=0
         nhom123p=0
         nhom123n=0
         nhom1234p=0
         nhom1234n=0
         do ia=1,nprotons
            ix=occ_p(ia,p_state)
            if (ix>nasps) nasps=ix
            if (nhomp<2*n_sp(ix)+l_sp(ix)) nhomp=2*n_sp(ix)+l_sp(ix)
            do i=ia+1,nprotons
               iy=occ_p(i,p_state)
               if (nhom12p<2*n_sp(ix)+l_sp(ix)+2*n_sp(iy)+l_sp(iy))
     +              nhom12p=2*n_sp(ix)+l_sp(ix)+2*n_sp(iy)+l_sp(iy)
               do ii=i+1,nprotons
                  ix_ii=occ_p(ii,p_state)
                  if (nhom123p<2*n_sp(ix_ii)+l_sp(ix_ii)
     +                 +2*n_sp(ix)+l_sp(ix)+2*n_sp(iy)+l_sp(iy))
     +                 nhom123p=2*n_sp(ix_ii)+l_sp(ix_ii)
     +                 +2*n_sp(ix)+l_sp(ix)+2*n_sp(iy)+l_sp(iy)
                  do iii=ii+1,nprotons
                     ix_iii=occ_p(iii,p_state)
                     if(nhom1234p<2*n_sp(ix_iii)+l_sp(ix_iii)
     +                    +2*n_sp(ix_ii)+l_sp(ix_ii)
     +                    +2*n_sp(ix)+l_sp(ix)+2*n_sp(iy)+l_sp(iy))
     +                    nhom1234p=2*n_sp(ix_iii)+l_sp(ix_iii)
     +                    +2*n_sp(ix_ii)+l_sp(ix_ii)
     +                    +2*n_sp(ix)+l_sp(ix)+2*n_sp(iy)+l_sp(iy)
                  end do
               end do
            end do
         end do
         do ia=1,nneutrns
            ix=occ_n(ia,n_state)
            if (ix>nasps) nasps=ix
            if (nhomn<2*n_sp(ix)+l_sp(ix)) nhomn=2*n_sp(ix)+l_sp(ix)
            do i=ia+1,nneutrns
               iy=occ_n(i,n_state)
               if (nhom12n<2*n_sp(ix)+l_sp(ix)+2*n_sp(iy)+l_sp(iy))
     +              nhom12n=2*n_sp(ix)+l_sp(ix)+2*n_sp(iy)+l_sp(iy)
               do ii=i+1,nneutrns
                  ix_ii=occ_n(ii,n_state)
                  if (nhom123n<2*n_sp(ix_ii)+l_sp(ix_ii)
     +                 +2*n_sp(ix)+l_sp(ix)+2*n_sp(iy)+l_sp(iy))
     +                 nhom123n=2*n_sp(ix_ii)+l_sp(ix_ii)
     +                 +2*n_sp(ix)+l_sp(ix)+2*n_sp(iy)+l_sp(iy)
                  do iii=ii+1,nneutrns
                     ix_iii=occ_n(iii,n_state)
                     if(nhom1234n<2*n_sp(ix_iii)+l_sp(ix_iii)
     +                    +2*n_sp(ix_ii)+l_sp(ix_ii)
     +                    +2*n_sp(ix)+l_sp(ix)+2*n_sp(iy)+l_sp(iy))
     +                    nhom1234n=2*n_sp(ix_iii)+l_sp(ix_iii)
     +                    +2*n_sp(ix_ii)+l_sp(ix_ii)
     +                    +2*n_sp(ix)+l_sp(ix)+2*n_sp(iy)+l_sp(iy)
                  end do
               end do
            end do
         end do
         if (nhom<max(nhomp,nhomn)) nhom=max(nhomp,nhomn)
         if (nhom12<max(nhom12p,nhom12n,nhomp+nhomn))
     +        nhom12=max(nhom12p,nhom12n,nhomp+nhomn)
         if (nhom123<max(nhom123p,nhom123n,nhom12p+nhomn,nhomp+nhom12n))
     +        nhom123=max(nhom123p,nhom123n,nhom12p+nhomn,nhomp+nhom12n)
         if (nhom1234<max(nhom1234p,nhom1234n,nhom123p+nhomn,
     +        nhomp+nhom123n,nhom12p+nhom12n))
     +        nhom1234=max(nhom1234p,nhom1234n,nhom123p+nhomn,
     +        nhomp+nhom123n,nhom12p+nhom12n)
      end do
      mxnwd=nasps/nbit
      if (mxnwd*nbit<nasps) mxnwd=mxnwd+1
      if (mxnwd*nbit/=mxsps) then
         allocate(ntemp(2*mxnwd*nbit))
         ntemp=0
         do ia=1,min(mxnwd*nbit,mxsps)
            ntemp(ia)=n_sp(ia)
            ntemp(mxnwd*nbit+ia)=n_sp(mxsps+ia)
         end do
         deallocate(n_sp)
         allocate(n_sp(2*mxnwd*nbit))
         n_sp=ntemp
         ntemp=0
         do ia=1,min(mxnwd*nbit,mxsps)
            ntemp(ia)=l_sp(ia)
            ntemp(mxnwd*nbit+ia)=l_sp(mxsps+ia)
         end do
         deallocate(l_sp)
         allocate(l_sp(2*mxnwd*nbit))
         l_sp=ntemp
         ntemp=0
         do ia=1,min(mxnwd*nbit,mxsps)
            ntemp(ia)=j2_sp(ia)
            ntemp(mxnwd*nbit+ia)=j2_sp(mxsps+ia)
         end do
         deallocate(j2_sp)
         allocate(j2_sp(2*mxnwd*nbit))
         j2_sp=ntemp
         ntemp=0
         do ia=1,min(mxnwd*nbit,mxsps)
            ntemp(ia)=m2_sp(ia)
            ntemp(mxnwd*nbit+ia)=m2_sp(mxsps+ia)
         end do
         deallocate(m2_sp)
         allocate(m2_sp(2*mxnwd*nbit))
         m2_sp=ntemp
         ntemp=0
         do ia=1,min(mxnwd*nbit,mxsps)
            ntemp(ia)=mt2_sp(ia)
            ntemp(mxnwd*nbit+ia)=mt2_sp(mxsps+ia)
         end do
         deallocate(mt2_sp)
         allocate(mt2_sp(2*mxnwd*nbit))
         mt2_sp=ntemp
         deallocate(ntemp)
         mxsps=mxnwd*nbit
      endif           

      nshll=nhom+1
      print *,' nhom_1=',nhom
      print *,' nhom12,nhom123,nhom1234=',nhom12,nhom123,nhom1234
      print *,' nbit,mxnwd=',nbit,mxnwd
      print *,' mxsps=',mxsps

      allocate(ibas_p(mxnwd,nsd_p))
      allocate(ibas_n(mxnwd,nsd_n))
      ibas_p=0
      ibas_n=0
      do ii=1,nsd_p
         do ia=1,nprotons
            locx=occ_p(ia,ii)-1
            iwd=locx/nbit
            ibit=locx-iwd*nbit
            iwd=iwd+1
            ibas_p(iwd,ii)=ibset(ibas_p(iwd,ii),ibit)
         end do
      end do
      do ii=1,nsd_n
         do ia=1,nneutrns
            locx=occ_n(ia,ii)-1
            iwd=locx/nbit
            ibit=locx-iwd*nbit
            iwd=iwd+1
            ibas_n(iwd,ii)=ibset(ibas_n(iwd,ii),ibit)
         end do
      end do

      nconf=0
      goto 12345  ! with hash cluster
12346 continue
      allocate(ntempi(nshll),ntemnu(nshll),npib(nshll),nnub(nshll))
      allocate(mxpi(nshll),mnpi(nshll),mxnu(nshll),mnnu(nshll),
     +     mxpn(nshll),mnpn(nshll))
      mnpi=0
      mnnu=0
      mnpn=0
      mxpi(1)=min(2,nprotons)
      mxnu(1)=min(2,nneutrns)
      mxpn(1)=min(4,nprotons+nneutrns)
      do i=2,nshll
         mxpi(i)=min(nprotons,nhw/(i-1),i*(i+1))
         mxnu(i)=min(nneutrns,nhw/(i-1),i*(i+1))
         mxpn(i)=min(nprotons+nneutrns,nhw/(i-1),2*i*(i+1))
      end do
*      print *,' before loopconf'
      nconf=0
      call loopconf(nshll,0)
      print *,' nconf=',nconf
      allocate(nprtn(nconf,nshll),nneut(nconf,nshll))
      nconf = 0
      call loopconf(nshll,1)
      deallocate(ntempi,ntemnu,npib,nnub)
      deallocate(mxpi,mnpi,mxnu,mnnu,mxpn,mnpn)

      allocate(conf_p(nshll,nsd_p),conf_n(nshll,nsd_n))
      conf_p=0
      conf_n=0
      do ii=1,nsd_p
         do ia=1,nprotons
            iy=occ_p(ia,ii)
            ixxx=0
            do i=1,nshll
               if (iy>ixxx.and.iy<=ixxx+i*(i+1)) then
                  conf_p(i,ii)=conf_p(i,ii)+1
                  exit
               endif
               ixxx=ixxx+i*(i+1)
            end do
         end do
**         if (ii<10.or.ii==1369) then
**         if (occ_p(2,ii)<=2.and.occ_p(3,ii)<=8) then
**            print *,'#',ii,' conf_p=',(conf_p(i,ii),i=1,nshll)
**            print *,'occ_p=',(occ_p(ia,ii),ia=1,nprotons)
**         endif
      end do
      do ii=1,nsd_n
         do ia=1,nneutrns
            iy=occ_n(ia,ii)
            ixxx=0
            do i=1,nshll
               if (iy>ixxx.and.iy<=ixxx+i*(i+1)) then
                  conf_n(i,ii)=conf_n(i,ii)+1
                  exit
               endif
               ixxx=ixxx+i*(i+1)
            end do
         end do
**         if (ii<10.or.ii==1369) then
**         if (occ_n(2,ii)<=2.and.occ_n(3,ii)<=8) then
**            print *,'#',ii,' conf_n=',(conf_n(i,ii),i=1,nshll)
**            print *,'occ_n=',(occ_n(ia,ii),ia=1,nneutrns)
**         endif
      end do

*      print *,' before sorting'
      allocate(mconf(nsd))
      allocate(iendconf(nconf))
      ii=1
      do i=1,nconf
         if (10*(i/10)==i) print *,' sorting conf i=',i
         ixxx=ii
         do i1=ixxx,nsd
            p_state=I_state(1,i1)
            n_state=I_state(2,i1)
            diff=0
            do iy=1,nshll
               diff=diff+abs(conf_p(iy,p_state)-nprtn(i,iy))
     $              +abs(conf_n(iy,n_state)-nneut(i,iy))
            end do
**            if (p_state==3271.and.n_state==3271) then
**               print *,'#',p_state,' conf_p=',
**     +              (conf_p(ix_ii,p_state),ix_ii=1,nshll)
**               print *,'#',n_state,' conf_n=',
**     +              (conf_n(ix_ii,n_state),ix_ii=1,nshll)
**               print *,' i,i1,ii,diff=',i,i1,ii,diff
**            endif
            if (diff==0) then
               I_state(:,i1)=I_state(:,ii)
               I_state(1,ii)=p_state
               I_state(2,ii)=n_state
               mconf(ii)=i
               do ixxx=kxxx,kxxx+nkxxx-1
                  en=bmp(i1,ixxx)
                  bmp(i1,ixxx)=bmp(ii,ixxx)
                  bmp(ii,ixxx)=en
               end do
**               if (i==1) then
**                  print *,'i=',i
**                  print *,' ii=',ii
**                  print *,' I_state=',I_state(1,ii),I_state(2,ii)
**                  print *,' mconf=',mconf(ii)
**               endif
               ii=ii+1
            endif
         end do
         iendconf(i)=ii-1
**         print *,' i,iendconf=',i,iendconf(i)
      end do
      print *,' ii-1,nsd=',ii-1,nsd
      deallocate(conf_p,conf_n)
12345 continue

      do ixxx=kxxx,kxxx+nkxxx-1
         if (mod(nucleons,2)==0) then
            jt2(ixxx)=2*nint(Jx(ixxx))
         else
            jt2(ixxx)=nint(2*Jx(ixxx))
            if (mod(jt2(ixxx),2)/=1) then
               if (2*Jx(ixxx)>jt2(ixxx)) then
                  jt2(ixxx)=nint(2*Jx(ixxx))+1
               else
                  jt2(ixxx)=nint(2*Jx(ixxx))-1
               endif
            endif
         endif
         if (dabs(dble(jt2(ixxx))-2.d0*Jx(ixxx))>1.d-4) then
            print *, ixxx,'**** not a good J ****'
         endif 
      end do

      if (iproc==0) then      
         write(iunitvout,*)
         write(iunitvout,*)'*** Nuclear states ***'

      if (iparity==0) then
         write(iunitvout,1000) nucleons,nprotons,nneutrns,mjtotal,
     + mttotal,hboi,nhw,nsd,nhme,k1max,
     + mxnwd,mxsps,major,iparity
 1000    format(' Nucleus:',/,' A=',i3,'   Z=',i3,'   N=',i3,
     + /,' 2*MJ=',i3,'   2*MT=',i3,'  parity= +',/,
     + ' hbar Omega=',f8.4,'   Nhw=',i3,'   dimension=',i8,'   nhme=',
     + i10,/,' k1max=',i3,'   mxnwd=',i3,'   mxsps=',i8,
     +                    '   major=',i2,'   iparity=',i2,/)
      elseif(iparity==1) then
         write(iunitvout,1007) nucleons,nprotons,nneutrns,mjtotal,
     + mttotal,hboi,nhw,nsd,nhme,k1max,
     + mxnwd,mxsps,major,iparity
 1007    format(' Nucleus:',/,' A=',i3,'   Z=',i3,'   N=',i3,
     + /,' 2*MJ=',i3,'   2*MT=',i3,'  parity= -',/,
     + ' hbar Omega=',f8.4,'   Nhw=',i3,'   dimension=',i8,'   nhme=',
     + i10,/,' k1max=',i3,'   mxnwd=',i3,'   mxsps=',i8,
     +                    '   major=',i2,'   iparity=',i2,/)
      else
         write(iunitvout,*) ' error: iparity=',iparity
         call MPI_Abort(icomm,1111,ierr)
         stop
      endif

      if (i_case==1) then
         write(iunitvout,1100) (Jx(k1),Txi(k1),ener(k1),
     +        ener(k1)-ener(1),
     +        k1=1,k1max+2)
      else
         write(iunitvout,1100) (Jx(k1),Txf(k1),ener(k1),
     +        ener(k1)-ener(1),
     +        k1=1,k1max+2)
      endif
 1100 format(' J=',f7.4,'    T=',f7.4,'     Energy=',f12.4,
     + '     Ex=',f12.4) 
      endif

      if (iproc==0) write(iunitvout,1156) nhom,nhom12,nasps
 1156    format(/,' N1_max=',i4,'   N12_max=',i4,'   Nasps=',i4)

c      do i1=1,2*mxsps
c         write(iunitvout,3434) n_sp(i1),l_sp(i1),j2_sp(i1),
c     +        m2_sp(i1),mt2_sp(i1)
c      end do

c      do i1=1,min(1000,nsd)
c         write(iunitvout,3434) (iloc(ia,i1),ia=1,nucleons)
c 3434    format(10i6)
c      end do

c      write(iunitvout,*)' used state #1:' 
c      write(iunitvout,3435) (bmp(i1,1),i1=1,nsd)
c 3435 format(10f8.4)
      return
 333  continue
      if (iproc==0) write(iunitvout,*) '**** error in reading ****'
      call MPI_Abort(icomm,3331,ierr)
      stop
      contains
      recursive subroutine loopconf(Nin,ipass)
      implicit none
      integer, intent(IN) :: Nin,ipass
      integer np,nn,iloop,nengy,nodd,i,npi,nnu

***      print *,' nshll, Nin, iloop:',nshll,Nin,iloop

      iloop=nshll-Nin+1
      do np = mxpi(iloop),mnpi(iloop),-1
         if (iloop==1) then
            npib(1)=np
         else
            if (npib(iloop-1)+np>nprotons) then
               cycle
            else
               npib(iloop)=npib(iloop-1)+np
            endif
         endif
         ntempi(iloop)=np
         do nn = min(mxnu(iloop),mxpn(iloop)-np),
     +           max(mnnu(iloop),mnpn(iloop)-np),-1
c**         do nn = mxnu(iloop),mnnu(iloop),-1
c**            if ((np+nn)<mnpn(iloop).or.(np+nn)>mxpn(iloop)) cycle 
            if (iloop==1) then
               nnub(1)=nn
            else
               if (nnub(iloop-1)+nn>nneutrns) then
                  cycle
               else
                  nnub(iloop)=nnub(iloop-1)+nn
               endif
            endif
            ntemnu(iloop)=nn
***            print *,' Nin,np,nn:',Nin,np,nn
***            print *,' iloop=',iloop
            if (Nin>1) then
               call loopconf(Nin-1,ipass)
            elseif (Nin==1) then

***            write(8,*)' #',iproc,' ntempi,ntempn:'
***            write(8,*)(ntempi(i),ntemnu(i),i=1,nshll)

               npi = 0
               nnu = 0
               nengy = 0
               do i = 1, nshll
                  npi = npi + ntempi(i)
                  nnu = nnu + ntemnu(i)
                  nengy = nengy + (i-1)*(ntempi(i)+ntemnu(i))
               enddo
               if (npi/=nprotons.or.nnu/=nneutrns.or.nengy>nhw) cycle

               nodd = 0
               do i = 1,nshll
                  if ((-1)**(i-1)==-1) then
                     nodd = nodd + ntempi(i) + ntemnu(i)
                  endif
               enddo
               if ((-1)**(nodd+iparity)==-1) cycle
c     This configuration survives the Nhw test and the parity test.
               nconf = nconf + 1
               if (ipass==0) cycle
cc               negy(nconf) = nengy
               do i = 1, nshll
                  nprtn(nconf,i) = ntempi(i)
                  nneut(nconf,i) = ntemnu(i)
               enddo
               if(iproc==0)then
**                  write(*,3) nconf,(nprtn(nconf,i),i=1,nshll)
**                  write(*,4) (nneut(nconf,i),i=1,nshll)
 3                format(1x,'Conf #',i5,16(i3))
 4                format(12x,16(i3))
               endif
               
            else
               print *,'*** Error in loopconf ***'
               stop
            endif    

         end do 
      end do   
      end subroutine loopconf 

      subroutine read_vector(retrieve,iprec,dim,k)
      implicit none
      logical,intent(IN) :: retrieve
      integer,intent(IN) :: iprec,dim,k
      integer,parameter :: lecvec=1000000
      real(kind=kind(0.d0)),allocatable :: bmp0(:)
      real(kind=kind(0.0)),allocatable :: bmps(:)
      integer :: i1,i,ii
      if (iprec==1) then
         allocate(bmps(min(dim,lecvec)))
         if (dim<=lecvec) then

            read(iunitvec) (bmps(i1),i1=1,dim)

c            write(33,'(10(e14.7))') (bmps(i1),i1=1,dim)
c            read(33,'(10(e14.7))') (bmps(i1),i1=1,dim)

c**            if (retrieve) bmp(:,k)=real(bmps(:),kind(0.d0))
            if (retrieve) then
               do i1=1,dim
                  bmp(i1,k)=real(bmps(i1),kind(0.d0))
               end do
            endif
         else
            ii=1
            do i=1,dim/lecvec

               read(iunitvec) (bmps(i1),i1=1,lecvec)

c               write(33,'(10(e14.7))') (bmps(i1),i1=1,lecvec)
c               read(33,'(10(e14.7))') (bmps(i1),i1=1,lecvec)

               if (retrieve) then
c**                  bmp(ii:min(ii+lecvec-1,dim),k)=
c**     +                 real(bmps(:),kind(0.d0))
                  do i1=1,lecvec
                     bmp(ii-1+i1,k)=real(bmps(i1),kind(0.d0))
                  end do
                  print *,' bmp(i:j)=',ii,min(ii+lecvec-1,dim)
               endif
               ii=ii+lecvec
            end do
            if (mod(dim,lecvec)>0) then

               read(iunitvec)
     +              (bmps(i1),i1=1,mod(dim,lecvec))

c               write(33,'(10(e14.7))') (bmps(i1),i1=1,mod(dim,lecvec))
c               read(33,'(10(e14.7))') (bmps(i1),i1=1,mod(dim,lecvec))

               if (retrieve) then
c**                  bmp(ii:dim,k)=real(bmps(1:mod(dim,lecvec)),kind(0.d0))
                  do i1=1,mod(dim,lecvec)
                     bmp(ii-1+i1,k)=real(bmps(i1),kind(0.d0))
                  end do
                  print *,' bmp(i:j)=',ii,dim,dim-ii+1,mod(dim,lecvec),
     +                 ii-1+mod(dim,lecvec)
               endif
            endif
         endif
         deallocate(bmps)
      else
         allocate(bmp0(min(dim,lecvec)))
         if (dim<=lecvec) then

            read(iunitvec) (bmp0(i1),i1=1,dim)

c            write(33,'(10(e14.7))') (bmp0(i1),i1=1,dim)
c            read(33,'(10(e14.7))') (bmp0(i1),i1=1,dim)

            if (retrieve) bmp(:,k)=bmp0(:)
         else
            ii=1
            do i=1,dim/lecvec

               read(iunitvec) (bmp0(i1),i1=1,lecvec)

c               write(33,'(10(e14.7))') (bmp0(i1),i1=1,lecvec)
c               read(33,'(10(e14.7))') (bmp0(i1),i1=1,lecvec)

               if (retrieve) then
                  bmp(ii:min(ii+lecvec-1,dim),k)=bmp0(:)
                  print *,' bmp(i:j)=',ii,min(ii+lecvec-1,dim)
               endif
               ii=ii+lecvec
            end do
            if (mod(dim,lecvec)>0) then

               read(iunitvec)
     +              (bmp0(i1),i1=1,mod(dim,lecvec))

c               write(33,'(10(e14.7))') (bmp0(i1),i1=1,mod(dim,lecvec))
c               read(33,'(10(e14.7))') (bmp0(i1),i1=1,mod(dim,lecvec))

               if (retrieve) then
                  bmp(ii:dim,k)=bmp0(1:mod(dim,lecvec))
                  print *,' bmp(i:j)=',ii,dim,dim-ii+1,mod(dim,lecvec)
               endif
            endif
         endif
         deallocate(bmp0)
      endif
      end subroutine read_vector
      end


      subroutine confwave_mbpt(nsd,mxnwd,nhom,mjtotal,mttotal,
     $     iparity,nhw,
     +     nucleons,nprotons,nneutrns,jt2,it2,Jx,Tx,hbo,iloc,nshll,
     +     nhom12,nhom123,nhom1234,nasps,mxsps,mconf,nconf,
     +     nprtn,nneut,iendconf,
     +     n_sp,l_sp,j2_sp,m2_sp,mt2_sp,bmp,ibasis,ener)
**** 15-November-2006 *****
**** reads in configuration-space many-body wave function ****
**** computed using the MBPT-NCSM code ****
**** Petr Navratil, LLNL *****
      use intrface
      use paramdef
      use nodeinfo
      use obdens, only: nasps_full 
      implicit none
      include 'mpif.h'
      integer(2),pointer :: iloc(:,:),mconf(:)
      integer(4),pointer :: nconf,nprtn(:,:),nneut(:,:),iendconf(:)
      integer(4),pointer:: jt2(:),it2(:)
      integer(4),pointer :: n_sp(:),l_sp(:),
     +     j2_sp(:),m2_sp(:),mt2_sp(:)
      real(kind=kind(0.0d0)),pointer::bmp(:,:)
      integer(4),pointer :: ibasis(:,:,:)
      integer(4),pointer :: nsd,mxnwd,nhom,nhom12,nasps,mxsps,nshll,
     +     nhom123,nhom1234
      integer(4),pointer :: mjtotal,mttotal,iparity,nhw
      integer(4),pointer :: nucleons,nprotons,nneutrns
      real(kind=kind(0.0d0)),pointer :: Jx(:),Tx(:),ener(:),delta_E(:)
      real(kind=kind(0.0d0)),pointer :: hbo
      character(len=20) :: interaction_id
      integer :: ph_selector,mxsps2,i,s,n,l2,j2,m2,mt2,ii,kk,min_oc,
     +     N_max,i1,iy,ia,nhme,k1max,ix,ierr,k1,ix_ii,ix_iii,iii,major,
     +     phase,N_HO,num_of_st
      real(kind(0.d0)) :: kap,norm
      integer(2),allocatable :: occ_tmp(:),sp_map(:)
      integer(4),allocatable :: ntemp(:)
      integer :: Npi_1min_full,Nnu_1min_full,N_1max,N_12max,
     $     Npim1_1min_full,Nnum1_1min_full
      integer(4),allocatable :: max_in(:),max_out(:)

      allocate(nsd,mxnwd,nhom,mjtotal,mttotal,iparity,nhw,nshll,nconf,
     +     nucleons,nprotons,nneutrns,hbo,nhom12,nhom123,nhom1234,
     +     nasps,mxsps)

      if (iproc==0) print *,' nsd...hbo allocated'

      open(iunitvec,file=filename,form='formatted',
     +     status='old',action='read')

      if (iproc==0) print *,' unit iunitvec opened'

      read(iunitvec,*,end=333,err=333) nprotons
      read(iunitvec,*,end=333,err=333) nneutrns
      read(iunitvec,'(a)',end=333,err=333) interaction_id
      read(iunitvec,*,end=333,err=333) hbo
      read(iunitvec,*,end=333,err=333) nshll
      read(iunitvec,*,end=333,err=333) mxsps2
cc      read(iunitvec,*,end=333,err=333) ph_selector
      read(iunitvec,*,end=333,err=333) N_max
cc      read(iunitvec,*,end=333,err=333) kap
      read(iunitvec,*,end=333,err=333) nsd
      read(iunitvec,*,end=333,err=333) iparity
      iparity=abs(iparity-1)/2
      read(iunitvec,*,end=333,err=333) mjtotal
      read(iunitvec,*,end=333,err=333) num_of_st
cc      num_of_st=1
      k1max=num_of_st-2

      nucleons=nprotons+nneutrns
      mxsps=mxsps2/2
      mttotal=nprotons-nneutrns
cc      mjtotal=0
cc      iparity=0
      nconf=0
      nhme=0
      major=2

      if (nprotons<=40.and.nneutrns<=40) then
         nhw=N_max+max(0,min(6,nprotons-2))+max(0,min(6,nneutrns-2))
     +        +2*max(0,min(12,nprotons-8))+2*max(0,min(12,nneutrns-8))
     +        +3*max(0,min(20,nprotons-20))+3*max(0,min(20,nneutrns-20))
         if (iproc==0) print *,' N_max,nhw=',N_max,nhw

         if (nprotons<=2) then
            Npi_1min_full=0
            Npim1_1min_full=0
         elseif (nprotons==3) then
            Npi_1min_full=1
            Npim1_1min_full=0
         elseif (nprotons<=8) then
            Npi_1min_full=1
            Npim1_1min_full=1
         elseif (nprotons==9) then
            Npi_1min_full=2
            Npim1_1min_full=1
         elseif (nprotons<=20) then
            Npi_1min_full=2
            Npim1_1min_full=2
         elseif (nprotons==21) then
            Npi_1min_full=3
            Npim1_1min_full=2
         elseif (nprotons<=40) then
            Npi_1min_full=3
            Npim1_1min_full=3
         endif
         if (nneutrns<=2) then
            Nnu_1min_full=0
            Nnum1_1min_full=0
         elseif (nneutrns==3) then
            Nnu_1min_full=1
            Nnum1_1min_full=0
         elseif (nneutrns<=8) then
            Nnu_1min_full=1
            Nnum1_1min_full=1
         elseif (nneutrns==9) then
            Nnu_1min_full=2
            Nnum1_1min_full=1
         elseif (nneutrns<=20) then
            Nnu_1min_full=2
            Nnum1_1min_full=2
         elseif (nneutrns==21) then
            Nnu_1min_full=3
            Nnum1_1min_full=2
         elseif (nneutrns<=40) then
            Nnu_1min_full=3
            Nnum1_1min_full=3
         endif
         N_1max=max(Npi_1min_full,Nnu_1min_full)+N_max
         N_12max=max(Npi_1min_full+Npim1_1min_full,
     $        Npi_1min_full+Nnu_1min_full,
     $        Nnu_1min_full+Nnum1_1min_full)+N_max

         if (iproc==0) print *,' N_1max,N_12max=',N_1max,N_12max

      else
         print *,' error: protons,neutrons=',nprotons,nneutrns
         stop
      endif
cc      print *,' ph_selector=',ph_selector
cc      print *,' kap=',kap

cc      k1max=-1
      allocate(jt2(num_of_st))
      allocate(it2(num_of_st))
      allocate(Jx(num_of_st),Tx(num_of_st))
      allocate(ener(num_of_st))
      allocate(delta_E(num_of_st))
cc      Jx(1)=0.d0
cc      Tx(1)=0.d0
cc      jt2(1)=nint(Jx(1))
cc      it2(1)=nint(Tx(1))

      do ii=1,num_of_st
cc         read(iunitvec,*,end=333,err=333) ener(ii)
         read(iunitvec,*,end=333,err=333) ener(ii),Jx(ii),Tx(ii)
     $        ,delta_E(ii)
         ener(ii)=ener(ii)+delta_E(ii)
         if (mod(nucleons,2)==0) then
            jt2(ii)=2*nint(Jx(ii))
            it2(ii)=2*nint(Tx(ii))
         else
            jt2(ii)=nint(2*Jx(ii))
            it2(ii)=nint(2*Tx(ii))
            if (mod(jt2(ii),2)/=1) then
               if (2*Jx(ii)>jt2(ii)) then
                  jt2(ii)=nint(2*Jx(ii))+1
               else
                  jt2(ii)=nint(2*Jx(ii))-1
               endif
            endif
            if (mod(it2(ii),2)/=1) then
               if (2*Tx(ii)>it2(ii)) then
                  it2(ii)=nint(2*Tx(ii))+1
               else
                  it2(ii)=nint(2*Tx(ii))-1
               endif
            endif
         endif
         if (dabs(dble(jt2(ii))-2.d0*Jx(ii))>1.d-4) then
            if (iproc==0) print *, ii,'**** not a good J ****'
         endif 
      end do

      read(iunitvec,*,end=333,err=333)
      read(iunitvec,*,end=333,err=333)
      read(iunitvec,*,end=333,err=333)

      allocate(n_sp(mxsps2),l_sp(mxsps2),j2_sp(mxsps2),
     +     m2_sp(mxsps2),mt2_sp(mxsps2))
      if (iproc==0) print *,' nshll=',nshll
      ii=0
      do mt2=1,-1,-2
         do N_HO=0,nshll-1
            do l2=mod(N_HO,2),N_HO,2
               n=(N_HO-l2)/2
               do j2=abs(2*l2-1),2*l2+1,2
                  do m2=-j2,j2,2
                     ii=ii+1
                     n_sp(ii)=n
                     l_sp(ii)=l2
                     j2_sp(ii)=j2
                     m2_sp(ii)=m2
                     mt2_sp(ii)=mt2      
                  end do
               end do
            end do
            if (mt2==1.and.N_HO==N_1max) nasps_full=ii
         end do
      end do
      if (iproc==0) print *,' mxsps2,ii=',mxsps2,ii
      if (iproc==0) print *,' nasps_full=',nasps_full
      if (mxsps2/=ii) stop
      if (iproc==0) print *,' n_sp..mt2_sp allocated'
      
      allocate(sp_map(0:mxsps2-1))
      do i=1,mxsps2
         read(iunitvec,*,end=333,err=333) s,n,l2,j2,m2,mt2
         do ii=1,mxsps2
            if (n_sp(ii)==n.and.l_sp(ii)==l2/2.and.j2_sp(ii)==j2
     +           .and.m2_sp(ii)==m2.and.mt2_sp(ii)==mt2) then
               sp_map(s)=ii
               exit
            endif
         end do
      end do

      read(iunitvec,*,end=333,err=333)
      read(iunitvec,*,end=333,err=333)
      read(iunitvec,*,end=333,err=333)

      allocate(iloc(nucleons,nsd))
      allocate(bmp(nsd,num_of_st))
      allocate(occ_tmp(nucleons))
      norm=0.d0
      do i=1,nsd
         read(iunitvec,*,end=333,err=333) (occ_tmp(ii),ii=1,nucleons),
     +        (bmp(i,kk),kk=1,num_of_st)
         do ii=1,nucleons
            occ_tmp(ii)=sp_map(occ_tmp(ii))
         end do
         phase=1
         do ii=1,nucleons-1
            do kk=ii+1,nucleons
               if (occ_tmp(kk)<occ_tmp(ii)) then
                  min_oc=occ_tmp(ii)
                  occ_tmp(ii)=occ_tmp(kk)
                  occ_tmp(kk)=min_oc
                  phase=-phase
               end if
            end do
         end do
         iloc(:,i)=occ_tmp(:)
         bmp(i,:)=bmp(i,:)*real(phase,kind(0.d0))
c         print *,iloc(:,i)
c         print *,bmp(i,1)
         norm=norm+bmp(i,1)**2
      end do
      deallocate(sp_map)
      close(iunitvec)
      if (iproc==0) print *,' iunitvec closed'

      if (iproc==0) print *,' bmp norm=',norm

      nasps=0
      nhom=0
      nhom12=0
      nhom123=0
      nhom1234=0
      do i1=iproc+1,nsd,nproc
         occ_tmp(:)=iloc(:,i1)
         do ia=1,nucleons
            iy=occ_tmp(ia)      !iloc(ia,i1)
            if (iy>mxsps) iy=iy-mxsps
            if (iy>nasps) nasps=iy
            if (nhom<2*n_sp(iy)+l_sp(iy))nhom=2*n_sp(iy)+l_sp(iy)
            do i=ia+1,nucleons
               ix=occ_tmp(i)    !iloc(i,i1)
               if(nhom12<2*n_sp(ix)+l_sp(ix)+2*n_sp(iy)+l_sp(iy))
     +              nhom12=2*n_sp(ix)+l_sp(ix)+2*n_sp(iy)+l_sp(iy)
               do ii=i+1,nucleons
                  ix_ii=occ_tmp(ii) !iloc(ii,i1)
                  if(nhom123<2*n_sp(ix_ii)+l_sp(ix_ii)
     +                 +2*n_sp(ix)+l_sp(ix)+2*n_sp(iy)+l_sp(iy))
     +                 nhom123=2*n_sp(ix_ii)+l_sp(ix_ii)
     +                 +2*n_sp(ix)+l_sp(ix)+2*n_sp(iy)+l_sp(iy)
                  do iii=ii+1,nucleons
                     ix_iii=occ_tmp(iii) !iloc(iii,i1)
                     if(nhom1234<2*n_sp(ix_iii)+l_sp(ix_iii)
     +                    +2*n_sp(ix_ii)+l_sp(ix_ii)
     +                    +2*n_sp(ix)+l_sp(ix)+2*n_sp(iy)+l_sp(iy))
     +                    nhom1234=2*n_sp(ix_iii)+l_sp(ix_iii)
     +                    +2*n_sp(ix_ii)+l_sp(ix_ii)
     +                    +2*n_sp(ix)+l_sp(ix)+2*n_sp(iy)+l_sp(iy)
                  end do
               end do
            end do
         end do   
      end do   
      deallocate(occ_tmp)
      allocate(max_in(5),max_out(5))
      max_in(1)=nasps
      max_in(2)=nhom
      max_in(3)=nhom12
      max_in(4)=nhom123
      max_in(5)=nhom1234
      call MPI_Allreduce(max_in,max_out,5,MPI_INTEGER,
     $     MPI_MAX,icomm,ierr)
      nasps=max_out(1)
      nhom=max_out(2)
      nhom12=max_out(3)
      nhom123=max_out(4)
      nhom1234=max_out(5)
      deallocate(max_in,max_out)
      if (iproc==0) then
         print *,' nasps,nhom=',nasps,nhom
         print *,' nhom12,nhom123,nhom1234=',nhom12,nhom123,nhom1234
      endif
c      nasps=nasps_full
      nhom=N_1max
      nhom12=N_12max
      if (iproc==0) then
         print *,' From full space: nasps,nhom=',nasps_full,nhom
         print *,' From full space (12 only): nhom12,nhom123,nhom1234=',
     $        nhom12,nhom123,nhom1234
      endif

      mxnwd=nasps_full/nbit
      if (mxnwd*nbit<nasps_full) mxnwd=mxnwd+1
c      mxnwd=nasps/nbit
c      if (mxnwd*nbit<nasps) mxnwd=mxnwd+1
      if (mxnwd*nbit/=mxsps) then
         allocate(ntemp(2*mxnwd*nbit))
         ntemp=0
         do ia=1,min(mxnwd*nbit,mxsps)
            ntemp(ia)=n_sp(ia)
            ntemp(mxnwd*nbit+ia)=n_sp(mxsps+ia)
         end do
         deallocate(n_sp)
         allocate(n_sp(2*mxnwd*nbit))
         n_sp=ntemp
         ntemp=0
         do ia=1,min(mxnwd*nbit,mxsps)
            ntemp(ia)=l_sp(ia)
            ntemp(mxnwd*nbit+ia)=l_sp(mxsps+ia)
         end do
         deallocate(l_sp)
         allocate(l_sp(2*mxnwd*nbit))
         l_sp=ntemp
         ntemp=0
         do ia=1,min(mxnwd*nbit,mxsps)
            ntemp(ia)=j2_sp(ia)
            ntemp(mxnwd*nbit+ia)=j2_sp(mxsps+ia)
         end do
         deallocate(j2_sp)
         allocate(j2_sp(2*mxnwd*nbit))
         j2_sp=ntemp
         ntemp=0
         do ia=1,min(mxnwd*nbit,mxsps)
            ntemp(ia)=m2_sp(ia)
            ntemp(mxnwd*nbit+ia)=m2_sp(mxsps+ia)
         end do
         deallocate(m2_sp)
         allocate(m2_sp(2*mxnwd*nbit))
         m2_sp=ntemp
         ntemp=0
         do ia=1,min(mxnwd*nbit,mxsps)
            ntemp(ia)=mt2_sp(ia)
            ntemp(mxnwd*nbit+ia)=mt2_sp(mxsps+ia)
         end do
         deallocate(mt2_sp)
         allocate(mt2_sp(2*mxnwd*nbit))
         mt2_sp=ntemp
         deallocate(ntemp)
         do ia=1,nucleons
            do i1=1,nsd
               iy=iloc(ia,i1)
               if (iy>mxsps)  then
                  iy=iy-mxsps+mxnwd*nbit
                  iloc(ia,i1)=iy
               endif   
            end do
         end do   
         mxsps=mxnwd*nbit
      endif           

      if (iproc==0) print *,' mxsps=',mxsps

      if (iproc==0) then      
         write(iunitvout,*)
         write(iunitvout,*)'*** Nuclear states ***'

      if (iparity==0) then
         write(iunitvout,1000) nucleons,nprotons,nneutrns,mjtotal,
     + mttotal,hbo,nhw,nsd,nhme,k1max,
     + mxnwd,mxsps,major,iparity
 1000    format(' Nucleus:',/,' A=',i3,'   Z=',i3,'   N=',i3,
     + /,' 2*MJ=',i3,'   2*MT=',i3,'  parity= +',/,
     + ' hbar Omega=',f8.4,'   Nhw=',i3,'   dimension=',i8,'   nhme=',
     + i10,/,' k1max=',i3,'   mxnwd=',i3,'   mxsps=',i8,
     +                    '   major=',i2,'   iparity=',i2,/)
      elseif(iparity==1) then
         write(iunitvout,1007) nucleons,nprotons,nneutrns,mjtotal,
     + mttotal,hbo,nhw,nsd,nhme,k1max,
     + mxnwd,mxsps,major,iparity
 1007    format(' Nucleus:',/,' A=',i3,'   Z=',i3,'   N=',i3,
     + /,' 2*MJ=',i3,'   2*MT=',i3,'  parity= -',/,
     + ' hbar Omega=',f8.4,'   Nhw=',i3,'   dimension=',i8,'   nhme=',
     + i10,/,' k1max=',i3,'   mxnwd=',i3,'   mxsps=',i8,
     +                    '   major=',i2,'   iparity=',i2,/)
      else
         write(iunitvout,*) ' error: iparity=',iparity
         call MPI_Abort(icomm,1111,ierr)
         stop
      endif

      if (select) then
         write(iunitvout,1100) (Jx(k1),Tx(k1),ener(k1),
     +        ener(k1)-ener(kxxx),k1=kxxx,kxxx+nkxxx-1)
      else
         write(iunitvout,1100) (Jx(k1),Tx(k1),ener(k1),
     +        ener(k1)-ener(1),k1=1,k1max+2)
      endif
 1100 format(' J=',f7.4,'    T=',f7.4,'     Energy=',f12.4,
     + '     Ex=',f12.4) 
      endif

      if (iproc==0) write(iunitvout,1156) nhom,nhom12,nasps
 1156 format(/,' N1_max=',i4,'   N12_max=',i4,'   Nasps=',i4)


c**********************************************
      return
 333  continue
      if (iproc==0) write(iunitvout,*) '**** error in reading ****'
      call MPI_Abort(icomm,3331,ierr)
      stop
      end


      subroutine confwave_redstick(nsd,mxnwd,nhom,mjtotal,mttotal,
     $     iparity,nhw,
     +     nucleons,nprotons,nneutrns,jt2,it2,Jx,Tx,hbo,iloc,nshll,
     +     nhom12,nhom123,nhom1234,nasps,mxsps,mconf,nconf,
     +     nprtn,nneut,iendconf,
     +     n_sp,l_sp,j2_sp,m2_sp,mt2_sp,bmp,ibasis,ener)
**** 15-November-2006 *****
**** reads in configuration-space many-body wave function ****
**** computed using the Redstick code ****
**** Petr Navratil, LLNL *****
      use intrface
      use paramdef
      use nodeinfo
      implicit none
      include 'mpif.h'
      integer(2),pointer :: iloc(:,:),mconf(:)
      integer(4),pointer :: nconf,nprtn(:,:),nneut(:,:),iendconf(:)
      integer(4),pointer:: jt2(:),it2(:)
      integer(4),pointer :: n_sp(:),l_sp(:),
     +     j2_sp(:),m2_sp(:),mt2_sp(:)
      real(kind=kind(0.0d0)),pointer::bmp(:,:)
      integer(4),pointer :: ibasis(:,:,:)
      integer(4),pointer :: nsd,mxnwd,nhom,nhom12,nasps,mxsps,nshll,
     +     nhom123,nhom1234
      integer(4),pointer :: mjtotal,mttotal,iparity,nhw
      integer(4),pointer :: nucleons,nprotons,nneutrns
      real(kind=kind(0.0d0)),pointer :: Jx(:),Tx(:),ener(:)
      real(kind=kind(0.0d0)),pointer :: hbo

      real(kind=kind(0.0d0)),allocatable :: Jxin(:),Txin(:),enerin(:),
     $     bmpin(:)
      character(len=100) :: interaction_id
      integer(2),allocatable :: occ_tmp(:),sp_map(:)
      integer :: mxsps2,N_max,num_of_st,ii,n,l,j2,m2,mt2,s,i,phase,kk
      real(kind(0.d0)) :: norm
      integer :: k1max,nhme,major,N_HO,min_oc,ia,i1,ix,iy,ix_ii,iii,
     $     ix_iii,ierr,k1,
     $     mxsps2_in,nshll_in,nshll_core,nprotons_in,nneutrns_in
      integer(4),allocatable :: ntemp(:)

      allocate(nsd,mxnwd,nhom,mjtotal,mttotal,iparity,nhw,nshll,nconf,
     +     nucleons,nprotons,nneutrns,hbo,nhom12,nhom123,nhom1234,
     +     nasps,mxsps)

      print *,' nsd...hbo allocated'

      open(iunitvec,file=filename,form='formatted',
     +     status='old',action='read')

      print *,' unit iunitvec opened'

      read(iunitvec,*,end=333,err=333) nprotons
      read(iunitvec,*,end=333,err=333) nneutrns
      read(iunitvec,'(a)',end=333,err=333) interaction_id
      read(iunitvec,*,end=333,err=333) hbo
      read(iunitvec,*,end=333,err=333) nshll_in  ! needed
      if (nshll_in==1) then
         read(iunitvec,*,end=333,err=333) mxsps2_in,nshll_core
         nshll=nshll_core+nshll_in
      else
         nshll=nshll_in
         read(iunitvec,*,end=333,err=333) mxsps2_in
      endif
      read(iunitvec,*,end=333,err=333) N_max  ! needed
      read(iunitvec,*,end=333,err=333) nsd
      read(iunitvec,*,end=333,err=333) iparity
      iparity=iparity/2
      read(iunitvec,*,end=333,err=333) mjtotal
      read(iunitvec,*,end=333,err=333) num_of_st
      k1max=num_of_st-2
      nucleons=nprotons+nneutrns
cc      mxsps=mxsps2/2
      mttotal=nprotons-nneutrns

      nconf=0
      nhme=0
      major=2

      if (nshll_in/=nshll) then
         nprotons_in=nprotons
         nneutrns_in=nneutrns
         select case (nshll_core)
         case(1)
            nprotons=nprotons_in+2
            nneutrns=nneutrns_in+2
         case(2)
            nprotons=nprotons_in+8
            nneutrns=nneutrns_in+8
         case(3)
            nprotons=nprotons_in+20
            nneutrns=nneutrns_in+20
         case default
            print *,' error: protons,neutrons=',nprotons,nneutrns
            stop
         end select
      endif

      if (nprotons<=40.and.nneutrns<=40) then
         nhw=N_max+max(0,min(6,nprotons-2))+max(0,min(6,nneutrns-2))
     +        +2*max(0,min(12,nprotons-8))+2*max(0,min(12,nneutrns-8))
     +        +3*max(0,min(20,nprotons-20))+3*max(0,min(20,nneutrns-20))
         print *,' N_max,nhw=',N_max,nhw
      else
         print *,' error: protons,neutrons=',nprotons,nneutrns
         stop
      endif

      if (nshll_in/=nshll) then
         nprotons=nprotons_in
         nneutrns=nneutrns_in
      endif

      if (kxxx+nkxxx-1>num_of_st) nkxxx=num_of_st+1-kxxx
      allocate(jt2(kxxx:kxxx+nkxxx-1))
      allocate(it2(kxxx:kxxx+nkxxx-1))
      allocate(Jx(kxxx:kxxx+nkxxx-1),Tx(kxxx:kxxx+nkxxx-1))
      allocate(ener(kxxx:kxxx+nkxxx-1))
      allocate(Jxin(num_of_st),Txin(num_of_st))
      allocate(enerin(num_of_st))

      do ii=1,num_of_st
         read(iunitvec,*,end=333,err=333) enerin(ii),Jxin(ii),Txin(ii)
      end do
      if (select) then
         do ii=1,nkxxx
            ener(kxxx+ii-1)=enerin(iselxxx(ii))
            Jx(kxxx+ii-1)=Jxin(iselxxx(ii))
            Tx(kxxx+ii-1)=Txin(iselxxx(ii))
         end do
      else
         Jx=Jxin
         Tx=Txin
         ener=enerin
      endif
      do ii=kxxx,kxxx+nkxxx-1
         if (mod(nucleons,2)==0) then
            jt2(ii)=2*nint(Jx(ii))
            it2(ii)=2*nint(Tx(ii))
         else
            jt2(ii)=nint(2*Jx(ii))
            it2(ii)=nint(2*Tx(ii))
            if (mod(jt2(ii),2)/=1) then
               if (2*Jx(ii)>jt2(ii)) then
                  jt2(ii)=nint(2*Jx(ii))+1
               else
                  jt2(ii)=nint(2*Jx(ii))-1
               endif
            endif
            if (mod(it2(ii),2)/=1) then
               if (2*Tx(ii)>it2(ii)) then
                  it2(ii)=nint(2*Tx(ii))+1
               else
                  it2(ii)=nint(2*Tx(ii))-1
               endif
            endif
         endif
         if (dabs(dble(jt2(ii))-2.d0*Jx(ii))>1.d-4) then
            print *, ii,'**** not a good J ****'
         endif 
      end do
      deallocate(Jxin,Txin,enerin)

      ii=0
      do mt2=1,-1,-2
         do N_HO=0,nshll-1
            do l=mod(N_HO,2),N_HO,2
               n=(N_HO-l)/2
               do j2=abs(2*l-1),2*l+1,2
                  do m2=-j2,j2,2
                     ii=ii+1
                  end do
               end do
            end do
         end do
      end do
      mxsps2=ii
      mxsps=mxsps2/2
      allocate(n_sp(mxsps2),l_sp(mxsps2),j2_sp(mxsps2),
     +     m2_sp(mxsps2),mt2_sp(mxsps2))
      print *,' nshll=',nshll
      ii=0
      do mt2=1,-1,-2
         do N_HO=0,nshll-1
            do l=mod(N_HO,2),N_HO,2
               n=(N_HO-l)/2
               do j2=abs(2*l-1),2*l+1,2
                  do m2=-j2,j2,2
                     ii=ii+1
                     n_sp(ii)=n
                     l_sp(ii)=l
                     j2_sp(ii)=j2
                     m2_sp(ii)=m2
                     mt2_sp(ii)=mt2      
                  end do
               end do
            end do
         end do
      end do
      print *,' mxsps2_in,mxsps2,ii=',mxsps2_in,mxsps2,ii
      if (mxsps2/=ii) stop
      print *,' n_sp..mt2_sp allocated'
      
      allocate(sp_map(mxsps2_in))
      do i=1,mxsps2_in
         read(iunitvec,*,end=333,err=333) s,n,l,j2,m2,mt2
         do ii=1,mxsps2
            if (n_sp(ii)==n.and.l_sp(ii)==l.and.j2_sp(ii)==j2
     +           .and.m2_sp(ii)==m2.and.mt2_sp(ii)==mt2) then
               sp_map(s)=ii
               exit
            endif
         end do
      end do

      allocate(iloc(nucleons,nsd))
      allocate(bmp(nsd,kxxx:kxxx+nkxxx-1))
      allocate(bmpin(num_of_st))
      allocate(occ_tmp(nucleons))
      norm=0.d0
      do i=1,nsd
         read(iunitvec,*,end=333,err=333) (occ_tmp(ii),ii=1,nucleons),
     +        (bmpin(kk),kk=1,num_of_st)
         if (select) then
            ii=1
            do kk=1,num_of_st
               if (kk==iselxxx(ii)) then
                  bmp(i,kxxx+ii-1)=bmpin(kk)
                  if (ii<nkxxx) ii=ii+1
               endif
            end do
         else
            bmp(i,:)=bmpin(:)
         endif
         do ii=1,nucleons
            occ_tmp(ii)=sp_map(occ_tmp(ii))
         end do
         phase=1
         do ii=1,nucleons-1
            do kk=ii+1,nucleons
               if (occ_tmp(kk)<occ_tmp(ii)) then
                  min_oc=occ_tmp(ii)
                  occ_tmp(ii)=occ_tmp(kk)
                  occ_tmp(kk)=min_oc
                  phase=-phase
               end if
            end do
         end do
         iloc(:,i)=occ_tmp(:)
         bmp(i,:)=bmp(i,:)*real(phase,kind(0.d0))
         norm=norm+bmp(i,kxxx)**2
      end do
      deallocate(occ_tmp)
      deallocate(sp_map)
      deallocate(bmpin)
      close(iunitvec)
      print *,' iunitvec closed'
      print *,' bmp norm=',norm
      if (select) then
         k1max=nkxxx-2
      endif

      nasps=0
      nhom=0
      nhom12=0
      nhom123=0
      nhom1234=0
      do ia=1,nucleons
         do i1=1,nsd
            iy=iloc(ia,i1)
            if (iy>mxsps) iy=iy-mxsps
            if (iy>nasps) nasps=iy
            if (nhom<2*n_sp(iy)+l_sp(iy))nhom=2*n_sp(iy)+l_sp(iy)
            do i=ia+1,nucleons
               ix=iloc(i,i1)
               if(nhom12<2*n_sp(ix)+l_sp(ix)+2*n_sp(iy)+l_sp(iy))
     +              nhom12=2*n_sp(ix)+l_sp(ix)+2*n_sp(iy)+l_sp(iy)
               do ii=i+1,nucleons
                  ix_ii=iloc(ii,i1)
                  if(nhom123<2*n_sp(ix_ii)+l_sp(ix_ii)
     +                 +2*n_sp(ix)+l_sp(ix)+2*n_sp(iy)+l_sp(iy))
     +                 nhom123=2*n_sp(ix_ii)+l_sp(ix_ii)
     +                 +2*n_sp(ix)+l_sp(ix)+2*n_sp(iy)+l_sp(iy)
                  do iii=ii+1,nucleons
                     ix_iii=iloc(iii,i1)
                     if(nhom1234<2*n_sp(ix_iii)+l_sp(ix_iii)
     +                    +2*n_sp(ix_ii)+l_sp(ix_ii)
     +                    +2*n_sp(ix)+l_sp(ix)+2*n_sp(iy)+l_sp(iy))
     +                    nhom1234=2*n_sp(ix_iii)+l_sp(ix_iii)
     +                    +2*n_sp(ix_ii)+l_sp(ix_ii)
     +                    +2*n_sp(ix)+l_sp(ix)+2*n_sp(iy)+l_sp(iy)
                  end do
               end do
            end do
         end do   
      end do   
      mxnwd=nasps/nbit
      if (mxnwd*nbit<nasps) mxnwd=mxnwd+1
      if (mxnwd*nbit/=mxsps) then
         allocate(ntemp(2*mxnwd*nbit))
         ntemp=0
         do ia=1,min(mxnwd*nbit,mxsps)
            ntemp(ia)=n_sp(ia)
            ntemp(mxnwd*nbit+ia)=n_sp(mxsps+ia)
         end do
         deallocate(n_sp)
         allocate(n_sp(2*mxnwd*nbit))
         n_sp=ntemp
         ntemp=0
         do ia=1,min(mxnwd*nbit,mxsps)
            ntemp(ia)=l_sp(ia)
            ntemp(mxnwd*nbit+ia)=l_sp(mxsps+ia)
         end do
         deallocate(l_sp)
         allocate(l_sp(2*mxnwd*nbit))
         l_sp=ntemp
         ntemp=0
         do ia=1,min(mxnwd*nbit,mxsps)
            ntemp(ia)=j2_sp(ia)
            ntemp(mxnwd*nbit+ia)=j2_sp(mxsps+ia)
         end do
         deallocate(j2_sp)
         allocate(j2_sp(2*mxnwd*nbit))
         j2_sp=ntemp
         ntemp=0
         do ia=1,min(mxnwd*nbit,mxsps)
            ntemp(ia)=m2_sp(ia)
            ntemp(mxnwd*nbit+ia)=m2_sp(mxsps+ia)
         end do
         deallocate(m2_sp)
         allocate(m2_sp(2*mxnwd*nbit))
         m2_sp=ntemp
         ntemp=0
         do ia=1,min(mxnwd*nbit,mxsps)
            ntemp(ia)=mt2_sp(ia)
            ntemp(mxnwd*nbit+ia)=mt2_sp(mxsps+ia)
         end do
         deallocate(mt2_sp)
         allocate(mt2_sp(2*mxnwd*nbit))
         mt2_sp=ntemp
         deallocate(ntemp)
         do ia=1,nucleons
            do i1=1,nsd
               iy=iloc(ia,i1)
               if (iy>mxsps)  then
                  iy=iy-mxsps+mxnwd*nbit
                  iloc(ia,i1)=iy
               endif   
            end do
         end do   
         mxsps=mxnwd*nbit
      endif           

      print *,' nasps,nhom=',nasps,nhom
      print *,' nhom12,nhom123,nhom1234=',nhom12,nhom123,nhom1234
      print *,' mxsps=',mxsps

      if (iproc==0) then      
         write(iunitvout,*)
         write(iunitvout,*)'*** Nuclear states ***'

      if (iparity==0) then
         write(iunitvout,1000) nucleons,nprotons,nneutrns,mjtotal,
     + mttotal,hbo,nhw,nsd,nhme,k1max,
     + mxnwd,mxsps,major,iparity
 1000    format(' Nucleus:',/,' A=',i3,'   Z=',i3,'   N=',i3,
     + /,' 2*MJ=',i3,'   2*MT=',i3,'  parity= +',/,
     + ' hbar Omega=',f8.4,'   Nhw=',i3,'   dimension=',i8,'   nhme=',
     + i10,/,' k1max=',i3,'   mxnwd=',i3,'   mxsps=',i8,
     +                    '   major=',i2,'   iparity=',i2,/)
      elseif(iparity==1) then
         write(iunitvout,1007) nucleons,nprotons,nneutrns,mjtotal,
     + mttotal,hbo,nhw,nsd,nhme,k1max,
     + mxnwd,mxsps,major,iparity
 1007    format(' Nucleus:',/,' A=',i3,'   Z=',i3,'   N=',i3,
     + /,' 2*MJ=',i3,'   2*MT=',i3,'  parity= -',/,
     + ' hbar Omega=',f8.4,'   Nhw=',i3,'   dimension=',i8,'   nhme=',
     + i10,/,' k1max=',i3,'   mxnwd=',i3,'   mxsps=',i8,
     +                    '   major=',i2,'   iparity=',i2,/)
      else
         write(iunitvout,*) ' error: iparity=',iparity
         call MPI_Abort(icomm,1111,ierr)
         stop
      endif

      if (select) then
         write(iunitvout,1100) (Jx(k1),Tx(k1),ener(k1),
     +        ener(k1)-ener(kxxx),k1=kxxx,kxxx+nkxxx-1)
      else
         write(iunitvout,1100) (Jx(k1),Tx(k1),ener(k1),
     +        ener(k1)-ener(1),k1=1,k1max+2)
      endif
 1100 format(' J=',f7.4,'    T=',f7.4,'     Energy=',f12.4,
     + '     Ex=',f12.4) 
      endif

      if (iproc==0) write(iunitvout,1156) nhom,nhom12,nasps
 1156 format(/,' N1_max=',i4,'   N12_max=',i4,'   Nasps=',i4)


c**********************************************
      return
 333  continue
      if (iproc==0) write(iunitvout,*) '**** error in reading ****'
      call MPI_Abort(icomm,3331,ierr)
      stop
      end

c*** MFDn_V13 begin ***
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
ccc      subroutine confwave_MFDn_V13(
      subroutine confwave_MFDj(
     $   nsd, mxnwd, nhom, mjtotal, mttotal,
     $   iparity, nhw, nucleons, nprotons, nneutrns, jt2, it2,
     $   Jx, Tx, hbo, iloc, nshll, nhom12, nhom123, nhom1234,
     $   nasps, mxsps, mconf, nconf, nprtn, nneut, iendconf,
     $   n_sp, l_sp, j2_sp, m2_sp, mt2_sp, bmp, ibasis, ener)
****              March/April 2010                        ****
**** reads in configuration-space many-body wave function ****
**** computed using MFDn Version 13                       ****
**** NOTE: each processor reads/stores entire wave function **
**** Pieter Maris, ISU / Petr Navratil, LLNL              ****
      use intrface              ! contains kxxx, nkxxx, select
      use paramdef              ! contains nbit
      use nodeinfo              ! contains communicator info
      implicit none
      include 'mpif.h'
C
C     Number of particles
      integer(4),pointer :: nucleons, nprotons, nneutrns
C     Single-particle basis parameters
      integer(4),pointer :: nshll, mxnwd, nasps, mxsps
C     Many-Body basis parameters
      integer(4),pointer :: mjtotal, mttotal, iparity, nhw, nsd
      integer(4),pointer :: nhom, nhom12, nhom123, nhom1234
C     Single-Particle basis 
      integer(4),pointer ::
     $   n_sp(:), l_sp(:), j2_sp(:), m2_sp(:), mt2_sp(:)
C     Many-body basis 
      integer(2),pointer :: iloc(:,:)
C     hbar-omega value
      real(kind=kind(0.0d0)),pointer :: hbo
C     Set of wavefunctions
      real(kind=kind(0.0d0)),pointer :: bmp(:,:)
C     Ebinding, J-value, and T-value for set of wavefunctions
      real(kind=kind(0.0d0)),pointer ::  ener(:), Jx(:),Tx(:)
C     Integer values for 2J and 2T for set of wavefunctions
      integer(4),pointer :: jt2(:), it2(:)
C     
C     Not used -- still in the interface because of legacy issues...
      integer(2),pointer :: mconf(:)
      integer(4),pointer :: nconf,nprtn(:,:),nneut(:,:),iendconf(:)
      integer(4),pointer :: ibasis(:,:,:)
C
C     local variables
      character(len=8) :: mbsiname, smwfname
      integer(8) :: totaldim
      integer(kind=MPI_OFFSET_KIND) :: disp_obs, disp_smwf
      integer(4) :: fh_basis, fh_smwf, ierr
      integer(4) :: ndiag, NOSPS, MAXCLS, mxsps2, intfour
      integer(4) :: hbarground, Nmin, deltaN, Nmaxindex
      integer(4) :: i, ia, i1, ix, iy, ii, ix_ii, iii, ix_iii
      integer(4) :: i_TRDENS, i_MFDn, iread, istore
      integer(4) :: numstates, numobs
      integer(4) :: major, nhme, k1max
      integer(4), dimension(:), allocatable ::
     $   numprtcls, numshellscls, numSPstatescls, totalnumstatesN
      integer(4), dimension(:,:), allocatable :: SPstates
      real(4), dimension(:), allocatable :: amp
      real(4) :: Ebinding, totalJ, totalT
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      mbsiname = filename
      smwfname = filename2
C
C     allocate non-array variables?
      allocate(nsd,mxnwd,nhom,mjtotal,mttotal,iparity,nhw,nshll,nconf,
     $   nucleons,nprotons,nneutrns,hbo,nhom12,nhom123,nhom1234,
     $   nasps,mxsps)
      print*,' nsd...hbo allocated'
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     open MBSI file, and read parameters
      call read_MFDn_param(mbsiname, MPI_COMM_WORLD,
     $   ndiag, nucleons, NOSPS, MAXCLS,
     $   hbarground, Nmin, deltaN, Nmaxindex, mjtotal,
     $   totaldim, fh_basis)
C
      if (MAXCLS.ne.2) then
         print*, MAXCLS
         print*, ' This code is specifically for protons AND neutrons'
         call MPI_Abort(MPI_COMM_WORLD, 401)
         stop
      endif
      nsd = intfour(totaldim, 402)
C
      allocate(numprtcls(MAXCLS),
     $   numshellscls(MAXCLS), numSPstatescls(MAXCLS))
      allocate(SPstates(8, NOSPS*MAXCLS))
      allocate(totalnumstatesN(Nmaxindex))
      allocate(iloc(nucleons, nsd))
C
      call read_MFDn_basis(fh_basis,
     $   ndiag, nucleons, NOSPS, MAXCLS, Nmaxindex,
     $   numprtcls, numshellscls, numSPstatescls, SPstates,
     $   nsd, totalnumstatesN, iloc)
C
C     close MBSI file
      call MPI_FILE_CLOSE(fh_basis, ierr)
C
      nprotons = numprtcls(1)
      nneutrns = numprtcls(2)
      mttotal = nprotons - nneutrns
C
      nshll = maxval(numshellscls)
C
      nhw = hbarground + Nmin + (Nmaxindex-1)*deltaN
      iparity = mod(nhw,2)      ! this assumes that deltaN = 2
C
      print*, iproc, ' SP and MB Basis read in from file', mbsiname
C     print to STDOUT from root
      if (iproc.eq.0) then
         print*
         print*, 'Number of particles ', nucleons
         print*, 'MAXCLS ', MAXCLS, 'NOSPS  ', NOSPS
         do i = 1, MAXCLS
            print*, 'Class', i,' fermions'
            print*, '   Number of particles ', numprtcls(i) 
            print*, '   Number of S.P.shells', numshellscls(i) 
            print*, '   Number of S.P.states', numSPstatescls(i) 
         enddo
C     
         print*
         print*, 'Number of H.O. quanta lowest configuration ',
     $      hbarground
         print*, 'Number of H.O. quanta start configuration  ',
     $      hbarground + Nmin
         print*, 'Number of H.O. quanta highest configuration ', nhw
         if ((mod(Nmin,2).eq.0).and.(deltaN.eq.2)) then
            print*, 'Natural parity states only'
         elseif ((mod(Nmin,2).eq.1).and.(deltaN.eq.2)) then
            print*, 'Unnatural parity states only'
         elseif (deltaN.eq.1) then
            print*, 'Both natural and unnatural parity states'
         else
            print*, 'Uncommon deltaN?', deltaN
         endif
         print*, 'Total M_j (spin projection) ', mjtotal
         print*
         do i = 1, Nmaxindex
            print*, 'Number of states up to N=',
     $         hbarground + Nmin + (i-1)*deltaN,' is',totalnumstatesN(i)
         enddo
         print*, 'Total Number of M.B. states ', totaldim
         print*
      endif
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     
C     convert Single-Particle basis from MFDn_V13 to TRDENS
C
C     determine max. number of SPstates in MBbasis
      nasps = 0
      do ia = 1, nucleons
         do i1 = 1, nsd
            iy = iloc(ia,i1)
            if (iy.gt.NOSPS) iy = iy - NOSPS
            if (iy.gt.nasps) nasps = iy
         end do
      end do   
C
      mxnwd = nasps/nbit
      if (mxnwd*nbit.lt.nasps) mxnwd = mxnwd+1
      mxsps = mxnwd*nbit
      mxsps2 = 2*mxsps
C
      allocate(n_sp(mxsps2), l_sp(mxsps2), j2_sp(mxsps2),
     $   m2_sp(mxsps2), mt2_sp(mxsps2))
      n_sp = 0
      l_sp = 0
      j2_sp = 0
      m2_sp = 0
      mt2_sp = 0
C     
C     proton SPstates
      do i = 1, nasps
         n_sp(i) = SPstates(1, i)
         l_sp(i) = SPstates(2, i)
         j2_sp(i) = SPstates(4, i)
         m2_sp(i) = SPstates(5, i)
         mt2_sp(i) = SPstates(7, i)
      enddo
C
C     neutron SPstates
      do i = 1, nasps
         i_TRDENS = i + mxsps
         i_MFDn = i + NOSPS
         n_sp(i_TRDENS) = SPstates(1, i_MFDn)
         l_sp(i_TRDENS) = SPstates(2, i_MFDn)
         j2_sp(i_TRDENS) = SPstates(4, i_MFDn)
         m2_sp(i_TRDENS) = SPstates(5, i_MFDn)
         mt2_sp(i_TRDENS) = SPstates(7, i_MFDn)
      enddo
C
C     adjust neutron SPstates in MB basis if necessary... (?)
      if (NOSPS.ne.mxsps) then
         do i1 = 1, nsd
            do ia = nprotons+1, nucleons
               iloc(ia,i1) = iloc(ia,i1) + mxsps - NOSPS
            enddo
         enddo   
      endif
C
      nhom=0
      nhom12=0
      nhom123=0
      nhom1234=0
      do ia=1,nucleons
         do i1=1,nsd
            iy=iloc(ia,i1)
            if (iy>mxsps) iy=iy-mxsps
            if (nhom<2*n_sp(iy)+l_sp(iy))nhom=2*n_sp(iy)+l_sp(iy)
            do i=ia+1,nucleons
               ix=iloc(i,i1)
               if(nhom12<2*n_sp(ix)+l_sp(ix)+2*n_sp(iy)+l_sp(iy))
     +              nhom12=2*n_sp(ix)+l_sp(ix)+2*n_sp(iy)+l_sp(iy)
               do ii=i+1,nucleons
                  ix_ii=iloc(ii,i1)
                  if(nhom123<2*n_sp(ix_ii)+l_sp(ix_ii)
     +                 +2*n_sp(ix)+l_sp(ix)+2*n_sp(iy)+l_sp(iy))
     +                 nhom123=2*n_sp(ix_ii)+l_sp(ix_ii)
     +                 +2*n_sp(ix)+l_sp(ix)+2*n_sp(iy)+l_sp(iy)
                  do iii=ii+1,nucleons
                     ix_iii=iloc(iii,i1)
                     if(nhom1234<2*n_sp(ix_iii)+l_sp(ix_iii)
     +                    +2*n_sp(ix_ii)+l_sp(ix_ii)
     +                    +2*n_sp(ix)+l_sp(ix)+2*n_sp(iy)+l_sp(iy))
     +                    nhom1234=2*n_sp(ix_iii)+l_sp(ix_iii)
     +                    +2*n_sp(ix_ii)+l_sp(ix_ii)
     +                    +2*n_sp(ix)+l_sp(ix)+2*n_sp(iy)+l_sp(iy)
                  end do
               end do
            end do
         end do   
      end do   
C
      print *,' mxsps,nbit,mxnwd=',mxsps,nbit,mxnwd
      print *,' nasps=',nasps
      print *,' nhom,nhom12,nhom123,nhom1234=',
     $     nhom,nhom12,nhom123,nhom1234
C
      deallocate(numshellscls, numSPstatescls)
      deallocate(SPstates)
      deallocate(totalnumstatesN)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     open smwf file, and read parameters
      call read_MFDn_smwf_info(MPI_COMM_WORLD, smwfname,
     $   MAXCLS, numprtcls, fh_smwf, disp_obs, disp_smwf,
     $   hbo, numobs, numstates)
C
C     reset nkxxx if numstates is insufficient
      if (select) then
         do while (maxval(iselxxx(1:nkxxx)).gt.numstates) 
            nkxxx = nkxxx - 1
            print*, ' Reducing nkxxx by one'
         enddo
      else
         if (kxxx+nkxxx-1.gt.numstates) then
            print*, ' Resetting kxxx+nkxxx-1 to numstates'
            print*, kxxx, kxxx+nkxxx-1, numstates
            nkxxx = numstates+1-kxxx
         endif
         if (numstates.le.0) then
            print*, ' Requesting state that is not in .smwf file'
            print*, kxxx, nkxxx, numstates
            call MPI_Abort(411, ierr)
         endif
      endif
C
      allocate(ener(kxxx:kxxx+nkxxx-1))
      allocate(Jx(kxxx:kxxx+nkxxx-1))
      allocate(Tx(kxxx:kxxx+nkxxx-1))
      allocate(jt2(kxxx:kxxx+nkxxx-1))
      allocate(it2(kxxx:kxxx+nkxxx-1))
C
      allocate(bmp(nsd,kxxx:kxxx+nkxxx-1))
      allocate(amp(nsd))
      do i = 1, nkxxx
         if (select) then
            iread = iselxxx(i)
            istore = i
         else
            iread = kxxx + i - 1
            istore = iread
         endif
C
         if (iread.le.numstates) then
            call read_MFDn_smwf(fh_smwf,
     $         disp_obs, numobs, disp_smwf, nsd,
     $         iread, Ebinding, totalJ, totalT, amp)
         else
            print*, ' Requesting state that is not in .smwf file'
            print*, i, iread, numstates
            call MPI_Abort(412, ierr)
            stop
         endif
C
         bmp(1:nsd, istore) = amp(1:nsd)
         ener(istore) = Ebinding
         Jx(istore) = totalJ
         jt2(istore) = nint(2.0*Jx(istore))
         if (abs(jt2(istore)-Jx(istore)).gt.1.0d-4) then
            print*, istore, '**** not a good J ****'
            print*, 'real J', Jx(istore), ' int 2J', jt2(istore)
         endif
         Tx(istore) = totalT
         it2(istore) = nint(2.0*Tx(istore))
      enddo
C
      call MPI_FILE_CLOSE(fh_smwf, ierr)
C
      deallocate(numprtcls)
      deallocate(amp)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      nconf=0
      nhme=0
      major=2
      if (iproc==0) then      
         write(iunitvout,*)
         write(iunitvout,*)'*** Nuclear states ***'

         if (iparity==0) then
            write(iunitvout,1000) nucleons,nprotons,nneutrns,mjtotal,
     +         mttotal,hbo,nhw,nsd,nhme,k1max,mxnwd,mxsps,major,iparity
 1000       format(' Nucleus:',/,' A=',i3,'   Z=',i3,'   N=',i3,
     + /,' 2*MJ=',i3,'   2*MT=',i3,'  parity= +',/,
     + ' hbar Omega=',f8.4,'   Nhw=',i3,'   dimension=',i8,'   nhme=',
     + i10,/,' k1max=',i3,'   mxnwd=',i3,'   mxsps=',i8,
     +                    '   major=',i2,'   iparity=',i2,/)
         elseif(iparity==1) then
            write(iunitvout,1007) nucleons,nprotons,nneutrns,mjtotal,
     +         mttotal,hbo,nhw,nsd,nhme,k1max,mxnwd,mxsps,major,iparity
 1007       format(' Nucleus:',/,' A=',i3,'   Z=',i3,'   N=',i3,
     + /,' 2*MJ=',i3,'   2*MT=',i3,'  parity= -',/,
     + ' hbar Omega=',f8.4,'   Nhw=',i3,'   dimension=',i8,'   nhme=',
     + i10,/,' k1max=',i3,'   mxnwd=',i3,'   mxsps=',i8,
     +                    '   major=',i2,'   iparity=',i2,/)
         else
            write(iunitvout,*) ' error: iparity=',iparity
            call MPI_Abort(MPI_COMM_WORLD,1111,ierr)
            stop
         endif

         write(iunitvout,1100)
     $      (Jx(i), Tx(i), ener(i), ener(i)-ener(kxxx),
     $      i=kxxx,kxxx+nkxxx-1)
 1100    format(' J=',f7.4,'    T=',f7.4,'     Energy=',f12.4,
     +      '     Ex=',f12.4) 

         write(iunitvout,1156) nhom,nhom12,nasps
 1156    format(/,' N1_max=',i4,'   N12_max=',i4,'   Nasps=',i4)
      endif
C
      return
C
ccc      end subroutine confwave_MFDn_V13
      end subroutine confwave_MFDj
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     
      subroutine read_MFDn_param(MBSIname, communicator, 
     $   ndiag, nparticles, NOSPS, MAXCLS,
     $   hbarground, Nmin, deltaN, Nmaxindex, totMj,
     $   totaldimension, fh)
C
      implicit none
      include "mpif.h"
C
      character*8, intent(in) :: MBSIname     
      integer, intent(in) :: communicator
      integer*4, intent(out) :: NOSPS, MAXCLS, fh
      integer, intent(out) :: ndiag, nparticles, 
     $   hbarground, Nmin, deltaN, Nmaxindex, totMj
      integer*8, intent(out) :: totaldimension
C
C     local variables
      integer*4 :: ierr, filemode
      integer :: status(MPI_STATUS_SIZE)
      integer(kind=MPI_OFFSET_KIND) :: ioffset, disp
      integer*8, dimension(:), allocatable :: totalnumstatesN
      integer*4, dimension(:), allocatable :: itemp
      integer :: icls, numlevels, partitionsize
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     fh = filehandler -- similar to Fortran unit number
      fh = 42
C
C     open file read-only
      filemode = MPI_MODE_RDONLY
      call MPI_FILE_OPEN(communicator, MBSIname,
     $   filemode, MPI_INFO_NULL, fh, ierr)
C
      disp = 0
C     create a file-view 
      call MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER4, MPI_INTEGER4,
     $   'native', MPI_INFO_NULL, ierr)
C
C     read MAXCLS, NOSPS
      ioffset = 0
      call MPI_FILE_READ_AT(fh, ioffset, MAXCLS, 1,
     $   MPI_INTEGER4, status, ierr)
      ioffset = 1
      call MPI_FILE_READ_AT(fh, ioffset, NOSPS, 1,
     $   MPI_INTEGER4, status, ierr)
C 
      ioffset = 2
      allocate(itemp(6*MAXCLS))
      call MPI_FILE_READ_AT(fh, ioffset, itemp, 6*MAXCLS,
     $   MPI_INTEGER4, status, ierr)
      nparticles = 0
      do icls = 1, MAXCLS
         nparticles = nparticles + itemp(1+(icls-1)*6)
      enddo
C
C     increase displacement
      disp = disp + 4 * (2 + 6*MAXCLS + 8*NOSPS*MAXCLS)
C
C     read partitioning size info on all processors
C     (needed to calculate correct disp)
      ioffset = 0
      call MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER4, MPI_INTEGER4,
     $   'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_READ_AT(fh, ioffset, itemp, 2,
     $   MPI_INTEGER4, status, ierr)
      numlevels = itemp(1)
      partitionsize = itemp(2)
C     increase displacement accordingly
      disp = disp + 4 * (2 + numlevels + partitionsize)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     create a file-view for reading of
C     hbarground, Nmin, deltaN, Nmaxindex, and totMj
      call MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER4, MPI_INTEGER4,
     $   'native', MPI_INFO_NULL, ierr)
C     actual read
      ioffset = 0
      call MPI_FILE_READ_AT(fh, ioffset, itemp, 6,
     $   MPI_INTEGER4, status, ierr)
      hbarground = itemp(1) 
      Nmin = itemp(2) 
      deltaN = itemp(3)
      Nmaxindex = itemp(4)
      totMj = itemp(5)
      ndiag = itemp(6)
C
C     allocate and read totalnumstatesN
      allocate(totalnumstatesN(Nmaxindex))
      disp = disp + 4 * 6
      totalnumstatesN = 0
      call MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER8, MPI_INTEGER8,
     $   'native', MPI_INFO_NULL, ierr)
      ioffset = 0
      call MPI_FILE_READ_AT(fh, ioffset, totalnumstatesN,
     $   Nmaxindex, MPI_INTEGER8, status, ierr)
C
C     calculate total dimensions
      totaldimension = totalnumstatesN(Nmaxindex)
C
      deallocate(totalnumstatesN)
      deallocate(itemp)
C
      return
C
      end subroutine read_MFDn_param
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     
      subroutine read_MFDn_basis(fh,
     $   ndiag, nparticles, NOSPS, MAXCLS, Nmaxindex,
     $   numprtcls, numshellscls, numSPstatescls, SPstates,
     $   nsd, totalnumstatesN, MBstates)
C
      implicit none
      include "mpif.h"
C
      integer*4, intent(in) :: NOSPS, MAXCLS, fh
      integer, intent(in) :: ndiag, nparticles,  Nmaxindex, nsd
      integer*4, dimension(MAXCLS), intent(out) ::
     $   numprtcls, numshellscls, numSPstatescls
      integer*4, dimension(8, NOSPS*MAXCLS), intent(out) :: SPstates
      integer*8, dimension(Nmaxindex), intent(out) :: totalnumstatesN
      integer*2, dimension(nparticles, nsd), intent(out) :: MBstates
C
C     local variables
      integer*4 :: ierr
      integer :: status(MPI_STATUS_SIZE)
      integer(kind=MPI_OFFSET_KIND) :: ioffset, disp
      integer*4, dimension(6*MAXCLS) :: itemp
      integer :: icls, numlevels, partitionsize
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      disp = 0
      call MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER4, MPI_INTEGER4,
     $   'native', MPI_INFO_NULL, ierr)
C
      ioffset = 2
      call MPI_FILE_READ_AT(fh, ioffset, itemp, 6*MAXCLS,
     $   MPI_INTEGER4, status, ierr)
      do icls = 1, MAXCLS
         numprtcls = itemp(1+(icls-1)*6)
         numshellscls = itemp(3+(icls-1)*6)
         numSPstatescls = itemp(5+(icls-1)*6)
      enddo
C
      ioffset = 2 + 6*MAXCLS
      call MPI_FILE_READ_AT(fh, ioffset, SPstates, 8*NOSPS*MAXCLS,
     $   MPI_INTEGER4, status, ierr)
C
C     increase displacement
      disp = disp + 4 * (ioffset + 8*NOSPS*MAXCLS)
C
C     read partitioning size info on all processors
C     (needed to calculate correct disp)
      ioffset = 0
      call MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER4, MPI_INTEGER4,
     $   'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_READ_AT(fh, ioffset, itemp, 2,
     $   MPI_INTEGER4, status, ierr)
      numlevels = itemp(1)
      partitionsize = itemp(2)
C     increase displacement accordingly
      disp = disp + 4 * (2 + numlevels + partitionsize)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     read totalnumstatesN
      disp = disp + 4 * 6
      totalnumstatesN = 0
      call MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER8, MPI_INTEGER8,
     $   'native', MPI_INFO_NULL, ierr)
      ioffset = 0
      call MPI_FILE_READ_AT(fh, ioffset, totalnumstatesN,
     $   Nmaxindex, MPI_INTEGER8, status, ierr)
C     increase displacement (need to know ndiag!)
      disp = disp + 8 * (Nmaxindex + (Nmaxindex+1)*ndiag)
C
C     create a file-view for reading of MBstates
      call MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER2, MPI_INTEGER2,
     $   'native', MPI_INFO_NULL, ierr)
C
C     read MB basis states (complete vector on each processor)
      ioffset = 0
      call MPI_FILE_READ_AT(fh, ioffset,
     $   MBstates(1:nparticles, 1:nsd),
     $   nparticles*nsd, MPI_INTEGER2, status, ierr)
C
      return
C
      end subroutine read_MFDn_basis
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C
      subroutine read_MFDn_smwf_info(communicator, smwfname,
     $   MAXCLS, numprt,
     $   fh, disp_obs, disp_smwf, hw, numobs, numstates)
C
      implicit none
      include "mpif.h"
C
      character*8, intent(in) :: smwfname
      integer, intent(in) :: communicator
      integer*4, intent(in) :: MAXCLS
      integer*4, dimension(MAXCLS), intent(in) :: numprt
C
      real*4, intent(out) :: hw
      integer*4, intent(out) :: fh, numobs, numstates
      integer(kind=MPI_OFFSET_KIND) :: disp_obs, disp_smwf
C
C     local variables
      integer :: status(MPI_STATUS_SIZE)
      integer(kind=MPI_OFFSET_KIND) :: ioffset, disp
      integer*4 :: icounter, icls, ierr, filemode
      integer*4, dimension(MAXCLS+3) :: itemp
C
C     fh = filehandler -- similar to Fortran unit number
      fh = 17
C     open .smwf file read-only
      filemode = MPI_MODE_RDONLY
      call MPI_FILE_OPEN(communicator, smwfname,
     $   filemode, MPI_INFO_NULL, fh, ierr)
C     read MAXCLS, numprt(MAXCLS), numstates, numobs (integers)
      disp = 0
      ioffset = 0      
      icounter = MAXCLS + 3
      call MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER4, MPI_INTEGER4,
     $   'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_READ_AT(fh, ioffset, itemp, icounter, 
     $   MPI_INTEGER, status, ierr)
C
      icounter = 1
      if (itemp(icounter).ne.MAXCLS) then
         print*, ' Something wrong with input files'
         print*, ' MAXCLS in .smwf and MBSI files different'
         print*, itemp(icounter), MAXCLS
         call MPI_Abort(MPI_COMM_WORLD, 31, ierr)
         stop
      endif
      do icls = 1, MAXCLS
         icounter = icounter+1
         if (itemp(icounter).ne.numprt(icls)) then
            print*, ' Something wrong with input files'
            print*, ' numprt in .smwf and MBSI files different'
            print*, itemp(icounter), numprt(icls)
            call MPI_Abort(MPI_COMM_WORLD, 32, ierr)
            stop
         endif
      enddo
      icounter = icounter + 1
      numstates = itemp(icounter)
      icounter = icounter + 1
      numobs = itemp(icounter) 
C
C     read hw (real)
      disp = icounter*4
      ioffset = 0
      call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL4, MPI_REAL4,
     $   'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_READ_AT(fh, ioffset, hw, 1,
     $   MPI_REAL4, status, ierr)
C
C     set displacement for observables
      disp_obs = disp
C     set displacement for smwf
      disp_smwf = (icounter + 1)*4 + numstates * numobs * 4
C
      return
C
      end subroutine read_MFDn_smwf_info
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      subroutine read_MFDn_smwf(fh,
     $   disp_obs, numobs, disp_smwf, nsd,
     $   statenumber, Ebinding, totalJ, totalT, amp)
C
      implicit none
      include "mpif.h"
C
      integer*4, intent(in) :: fh, numobs, nsd, statenumber
      integer(kind=MPI_OFFSET_KIND), intent(in) :: disp_obs, disp_smwf
      real*4, intent(out) :: Ebinding, totalJ, totalT
      real*4, dimension(nsd), intent(out) :: amp
C
C     local variables
      integer :: status(MPI_STATUS_SIZE), ierr
      integer(kind=MPI_OFFSET_KIND) :: ioffset
      real*4, dimension(numobs) :: arrayobs
C
C     read observables for statenumber (real*4)
      call MPI_FILE_SET_VIEW(fh, disp_obs, MPI_REAL4, MPI_REAL4,
     $   'native', MPI_INFO_NULL, ierr)
C
      ioffset = 1 + numobs*(statenumber-1)
      call MPI_FILE_READ_AT(fh, ioffset, arrayobs, numobs,
     $   MPI_REAL4, status, ierr)
      Ebinding = arrayobs(1)
      totalJ = arrayobs(2)
      if (numobs.gt.2) then
         totalT = arrayobs(3)
      else
         totalT = -99.0
      endif
C
C     read smwf for statenumber (complete vector on each processor)
      call MPI_FILE_SET_VIEW(fh, disp_smwf, MPI_REAL4, MPI_REAL4,
     $   'native', MPI_INFO_NULL, ierr)
C
      ioffset = (statenumber-1)*nsd
      call MPI_FILE_READ_AT(fh, ioffset, amp, nsd,
     $   MPI_REAL4, status, ierr)
C
      return
C
      end subroutine read_MFDn_smwf
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      integer*4 function intfour(neight, icode)
C
      use nodeinfo
      include "mpif.h"
      integer*8 neight
      integer*4 icode, ierr
C
C     converts integer*8 to integer*4
C     aborts MPI and stops if neight > 2^31
C
      intfour = neight
      if (intfour.ne.neight) then
         if (iproc.eq.0) then         
            print*, ' Fatal error in program MFDn'
            print*, ' integer*8 too big in function intfour', neight
            print*, ' Error code', icode
         endif
         call MPI_Abort(MPI_COMM_WORLD,icode,ierr)
         stop
      endif
C
      return
      end function intfour
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c*** MFDn_V13 end ***

      subroutine onebodyinit
      use initst, only: n_spi,l_spi,j2_spi,m2_spi,mt2_spi,mxnwdi,
     +     naspsi,nhom12i,nhom123i,nhom1234i,nucleonsi,iparityi,jt2i,
     $     it2i
      use finast, only: nhom12f,nhom123f,nhom1234f,nucleonsf,
     +     iparityf,jt2f
      use obdens
      use paramdef
      use nodeinfo
      use intrface, only: twobdcal,ki,nki,kf,nkf,threebdcal
      implicit none
      integer :: ii,mt2spo,ind,ij2x,ind_pi,J12ma
      integer :: isum,ij,J12,iT12,ipi,J123,iT123,ia,ib,na,la,nb,lb,
     +     nc,lc,j2c,J1234,iT1234,nd,ld,j2d,ii2,ii1
      integer :: n_ii,l_ii,j2_ii,n_ij,l_ij,j2_ij
      integer :: nhom12max,nhom123max,nhom1234max,ipimin,ipimax
      integer :: J12min,J12max,jtotal2min,jtotal2,kii,kff,jtot2min,jtot2
      integer :: t12min,t12max,ttotal2min,ttotal2,ttot2min,ttot2
      mt2spo=mt2_spi(1)
      ind=1
      ij2x=1
      ii=1
c*         print *,'#',ii,' n=',n_spi(ii),' l=',l_spi(ii),
c*     +   ' j=',j2_spi(ii),' m=',m2_spi(ii),' mt=',mt2_spi(ii)
c*         print *,' ind=',ind
      do ii=2,mxnwdi*nbit
c         if (ii>naspsi) cycle
         if (ii>max(naspsi,nasps_full)) cycle
         if (mt2_spi(ii)/=1.and.mt2_spi(ii)/=-1) cycle
         if (mt2spo/=mt2_spi(ii)) then
            print *,'*** error in onebodyinit: mt2spo,mt2_spi(ii)=',
     +           ii,mt2spo,mt2_spi(ii)
            print *,'#',ii,' n=',n_spi(ii),' l=',l_spi(ii),
     +           ' j=',j2_spi(ii),' m=',m2_spi(ii),' mt=',mt2_spi(ii)
            print *,' ind=',ind
            stop
         endif
         if (n_spi(ii)/=n_spi(ii-1).or.l_spi(ii)/=l_spi(ii-1)
     +        .or.j2_spi(ii)/=j2_spi(ii-1)) then
            ind=ind+1
            if (ij2x/=j2_spi(ii-1)+1) then
               print *,'*** error: ij2x,2*j+1:',ii,ij2x,
     +              j2_spi(ii-1)+1
               stop
            endif   
            ij2x=1
         else
            ij2x=ij2x+1
         endif           
c*         print *,'#',ii,' n=',n_spi(ii),' l=',l_spi(ii),
c*     +   ' j=',j2_spi(ii),' m=',m2_spi(ii),' mt=',mt2_spi(ii)
c*         print *,' ind=',ind
      end do   
      ist1dim=ind
      allocate(ist1(3,ist1dim))
      allocate(iobsind(2*mxnwdi*nbit))
      ist1(1,1)=n_spi(1)
      ist1(2,1)=l_spi(1)
      ist1(3,1)=j2_spi(1)
      iobsind(1)=1
      ind=1
      do ii=2,mxnwdi*nbit
c         if (ii>naspsi) cycle
         if (ii>max(naspsi,nasps_full)) cycle
         if (mt2_spi(ii)/=1.and.mt2_spi(ii)/=-1) cycle
         if (n_spi(ii)/=n_spi(ii-1).or.l_spi(ii)/=l_spi(ii-1)
     +        .or.j2_spi(ii)/=j2_spi(ii-1)) then
            ind=ind+1
            ist1(1,ind)=n_spi(ii)
            ist1(2,ind)=l_spi(ii)
            ist1(3,ind)=j2_spi(ii)
         endif           
         iobsind(ii)=ind
      end do   
      iobsind(mxnwdi*nbit+1:2*mxnwdi*nbit)=iobsind(1:mxnwdi*nbit)
      if (twobdcal.or.abs(nucleonsi-nucleonsf)>1) then
         nhom12max=max(nhom12i,nhom12f)
         if (abs(nucleonsi-nucleonsf)==2) then
            ipimax=(-1)**(iparityi+iparityf)
            ipimin=ipimax
            jtotal2=(jt2i(ki)+jt2f(kf))/2
            jtotal2min=iabs(jt2i(ki)-jt2f(kf))/2
            do kii=ki,ki+nki-1
               do kff=kf,kf+nkf-1
                  jtot2min=iabs(jt2i(kii)-jt2f(kff))/2
                  jtot2=(jt2i(kii)+jt2f(kff))/2
                  if (jtotal2min>jtot2min) jtotal2min=jtot2min
                  if (jtotal2<jtot2) jtotal2=jtot2
               end do
            end do 
            J12min=max(0,jtotal2min)
            J12max=min(nhom12max+1,jtotal2)
         else
            ipimax=1
            ipimin=-1
            J12min=0
            J12max=nhom12max+1
         endif
         ist2_Jst_J12min=J12min
         ist2_Jst_J12max=J12max
         allocate(itbind(ist1dim*(ist1dim+1)/2,J12min:J12max,0:1))
         allocate(ist2_Jst(0:1,J12min:J12max+1))
         ind=0
         isum2=0
         do ipi=ipimax,ipimin,-2
c            ind_pi=0
            if (ipi/=(-1)**nhom12max) then
               J12ma=J12max-1
            else
               J12ma=J12max
            endif
            do J12=J12min,J12ma
               ist2_Jst((ipi+1)/2,J12)=ind+1
               if (iproc==0) then
                  print *,' ipi,J12,ist2_Jst=',ipi,J12,
     $                 ist2_Jst((ipi+1)/2,J12)
               endif
               do iT12=0,1
                  isum=0
                  do ii=1,ist1dim
                     n_ii=ist1(1,ii)
                     l_ii=ist1(2,ii)
                     j2_ii=ist1(3,ii)
                     do ij=ii,ist1dim
                        n_ij=ist1(1,ij)
                        l_ij=ist1(2,ij)
                        j2_ij=ist1(3,ij)
                        if ((-1)**(l_ii+l_ij)/=ipi) cycle
                        if (2*n_ii+l_ii+2*n_ij+l_ij>nhom12max) cycle
                        if (iabs(j2_ii-j2_ij)>2*J12) cycle
                        if (j2_ii+j2_ij<2*J12) cycle
                        if (ii==ij) then
                           if ((-1)**(J12+iT12)/=-1) cycle
                        endif
                        ind=ind+1
                        isum=isum+1
c                        ind_pi=ind_pi+1
                        itbind(ii+ij*(ij-1)/2,J12,iT12)=ind
                     end do
                  end do   
                  isum2=isum2+isum*(isum+1)/2
               end do   
            end do
            ist2_Jst((ipi+1)/2,J12ma+1)=ind+1
            if (iproc==0) then
               print *,' ipi,J12max+1,ist2_Jst=',
     $              ipi,J12ma+1,ist2_Jst((ipi+1)/2,J12ma+1)
            endif
         end do   
         ist2dim=ind
         if (iproc==0) then
            print *,' ist2dim=',ist2dim
            write(*,3456) isum2
 3456       format(/,' Number of hamiltonian 2-b matrix elements:',i7,/)
         endif
         allocate(ist2(4,ist2dim))
         ind=0
         do ipi=ipimax,ipimin,-2
            if (ipi/=(-1)**nhom12max) then
               J12ma=J12max-1
            else
               J12ma=J12max
            endif
            do J12=J12min,J12ma
               do iT12=0,1
                  do ii=1,ist1dim
                     n_ii=ist1(1,ii)
                     l_ii=ist1(2,ii)
                     j2_ii=ist1(3,ii)
                     do ij=ii,ist1dim
                        n_ij=ist1(1,ij)
                        l_ij=ist1(2,ij)
                        j2_ij=ist1(3,ij)
                        if ((-1)**(l_ii+l_ij)/=ipi) cycle
                        if (2*n_ii+l_ii+2*n_ij+l_ij>nhom12max) cycle
                        if (iabs(j2_ii-j2_ij)>2*J12) cycle
                        if (j2_ii+j2_ij<2*J12) cycle
                        if (ii==ij) then
                           if ((-1)**(J12+iT12)/=-1) cycle
                        endif
                        ind=ind+1
                        ist2(1,ind)=ii
                        ist2(2,ind)=ij
                        ist2(3,ind)=J12
                        ist2(4,ind)=iT12
                     end do
                  end do   
               end do   
            end do   
         end do   
         if (ind==ist2dim) then
            if (iproc==0) write(*,3435) ist2dim
 3435       format(' Number of two-nucleon states:',i5)
         else
            print *,
     +           '*** error in two-body state generation:ind,ist2dim=',
     +           ind,ist2dim
            stop
         endif   
      endif   
      if (threebdcal.or.abs(nucleonsi-nucleonsf)>2) then
         nhom123max=max(nhom123i,nhom123f)
         if (abs(nucleonsi-nucleonsf)==3) then
            ipimax=(-1)**(iparityi+iparityf)
            ipimin=ipimax
            jtotal2=jt2i(ki)+jt2f(kf)
            jtotal2min=iabs(jt2i(ki)-jt2f(kf))
            do kii=ki,ki+nki-1
               do kff=kf,kf+nkf-1
                  jtot2min=iabs(jt2i(kii)-jt2f(kff))
                  jtot2=jt2i(kii)+jt2f(kff)
                  if (jtotal2min>jtot2min) jtotal2min=jtot2min
                  if (jtotal2<jtot2) jtotal2=jtot2
               end do
            end do 
            J12min=max(1,jtotal2min)
            J12max=min(2*nhom123max+3,jtotal2)
            t12min=1
            t12max=3
         elseif (threebdcal.and.nucleonsi==3) then
            ipimax=(-1)**iparityi
            ipimin=ipimax
            jtotal2=jt2i(ki)
            jtotal2min=iabs(jt2i(ki))
            ttotal2=it2i(ki)
            ttotal2min=iabs(it2i(ki))
            do kii=ki,ki+nki-1
               jtot2min=iabs(jt2i(kii))
               jtot2=jt2i(kii)
               if (jtotal2min>jtot2min) jtotal2min=jtot2min
               if (jtotal2<jtot2) jtotal2=jtot2
               ttot2min=iabs(it2i(kii))
               ttot2=it2i(kii)
               if (ttotal2min>ttot2min) ttotal2min=ttot2min
               if (ttotal2<ttot2) ttotal2=ttot2
            end do 
            J12min=max(1,jtotal2min)
            J12max=min(2*nhom123max+3,jtotal2)
            t12min=max(1,ttotal2min)
            t12max=min(3,ttotal2)
         elseif (threebdcal.and.nucleonsi==4) then
            ipimax=1
            ipimin=-1
            j2_ii=ist1(3,1)
            jtotal2=jt2i(ki)+j2_ii
            jtotal2min=iabs(jt2i(ki)-j2_ii)
            ttotal2=it2i(ki)+1
            ttotal2min=iabs(it2i(ki)-1)
            do kii=ki,ki+nki-1
               do ii=1,ist1dim
                  j2_ii=ist1(3,ii)         
                  jtot2min=iabs(jt2i(kii)-j2_ii)
                  jtot2=jt2i(kii)+j2_ii
                  if (jtotal2min>jtot2min) jtotal2min=jtot2min
                  if (jtotal2<jtot2) jtotal2=jtot2
               end do
               ttot2min=iabs(it2i(kii)-1)
               ttot2=it2i(kii)+1
               if (ttotal2min>ttot2min) ttotal2min=ttot2min
               if (ttotal2<ttot2) ttotal2=ttot2
            end do 
            J12min=max(1,jtotal2min)
            J12max=min(2*nhom123max+3,jtotal2)
            t12min=max(1,ttotal2min)
            t12max=min(3,ttotal2)
         else
            ipimax=1
            ipimin=-1
            J12min=1
            J12max=2*nhom123max+3
            t12min=1
            t12max=3
         endif
         allocate(ad3ind(ist2dim,ist1dim,J12min/2:J12max/2,
     $        t12min/2:t12max/2))
         ind=0
         do J123=J12min,J12max,2
            do iT123=t12min,t12max,2
               do ipi=ipimax,ipimin,-2
                  do ii=1,ist2dim
                     ia=ist2(1,ii)
                     ib=ist2(2,ii)
                     J12=ist2(3,ii)
                     iT12=ist2(4,ii)
                     na=ist1(1,ia)
                     la=ist1(2,ia)
                     nb=ist1(1,ib)
                     lb=ist1(2,ib)
                     do ij=1,ist1dim
                        nc=ist1(1,ij)
                        lc=ist1(2,ij)
                        j2c=ist1(3,ij)
                        if ((-1)**(la+lb+lc)/=ipi) cycle
                        if (2*na+la+2*nb+lb+2*nc+lc>nhom123max) cycle
                        if (abs(2*J12-j2c)>J123.or.(2*J12+j2c)<J123) 
     +                       cycle
                        if (abs(2*iT12-1)>iT123.or.(2*iT12+1)<iT123) 
     +                       cycle
                        ind=ind+1
                        ad3ind(ii,ij,J123/2,iT123/2)=ind
                     end do
                  end do
               end do
            end do
         end do
         ist3dim=ind
         allocate(ist3(4,ist3dim))
         ind=0
         do J123=J12min,J12max,2
            do iT123=t12min,t12max,2
               do ipi=ipimax,ipimin,-2
                  do ii=1,ist2dim
                     ia=ist2(1,ii)
                     ib=ist2(2,ii)
                     J12=ist2(3,ii)
                     iT12=ist2(4,ii)
                     na=ist1(1,ia)
                     la=ist1(2,ia)
                     nb=ist1(1,ib)
                     lb=ist1(2,ib)
                     do ij=1,ist1dim
                        nc=ist1(1,ij)
                        lc=ist1(2,ij)
                        j2c=ist1(3,ij)
                        if ((-1)**(la+lb+lc)/=ipi) cycle
                        if (2*na+la+2*nb+lb+2*nc+lc>nhom123max) cycle
                        if (abs(2*J12-j2c)>J123.or.(2*J12+j2c)<J123) 
     +                       cycle
                        if (abs(2*iT12-1)>iT123.or.(2*iT12+1)<iT123) 
     +                       cycle
                        ind=ind+1
                        ist3(1,ind)=ii
                        ist3(2,ind)=ij
                        ist3(3,ind)=J123
                        ist3(4,ind)=iT123
                     end do
                  end do
               end do
            end do
         end do
         if (ind==ist3dim) then
            if (iproc==0)
     $           write(*,"(' Number of three-nucleon states:',i7)")
     $           ist3dim
         else
            print *,
     +           '*** error in three-body state generation:ind,ist3dim='
     +           ,ind,ist3dim
            stop
         endif   
      endif
      if ((threebdcal.and.nucleonsi==4)
     $     .or.abs(nucleonsi-nucleonsf)==4) then
         nhom1234max=max(nhom1234i,nhom1234f)
         if (abs(nucleonsi-nucleonsf)==4) then
            ipimax=(-1)**(iparityi+iparityf)
            ipimin=ipimax
            jtotal2=(jt2i(ki)+jt2f(kf))/2
            jtotal2min=iabs(jt2i(ki)-jt2f(kf))/2
            do kii=ki,ki+nki-1
               do kff=kf,kf+nkf-1
                  jtot2min=iabs(jt2i(kii)-jt2f(kff))/2
                  jtot2=(jt2i(kii)+jt2f(kff))/2
                  if (jtotal2min>jtot2min) jtotal2min=jtot2min
                  if (jtotal2<jtot2) jtotal2=jtot2
               end do
            end do 
            J12min=max(0,jtotal2min)
            J12max=min(nhom1234max+2,jtotal2)
            t12min=0
            t12max=2
         elseif (threebdcal.and.nucleonsi==4) then
            nhom1234max=nhom1234i
            ipimax=(-1)**iparityi
            ipimin=ipimax
            jtotal2=jt2i(ki)
            jtotal2min=iabs(jt2i(ki))
            ttotal2=it2i(ki)
            ttotal2min=iabs(it2i(ki))
            do kii=ki,ki+nki-1
               jtot2min=iabs(jt2i(kii))
               jtot2=jt2i(kii)
               if (jtotal2min>jtot2min) jtotal2min=jtot2min
               if (jtotal2<jtot2) jtotal2=jtot2
               ttot2min=iabs(it2i(kii))
               ttot2=it2i(kii)
               if (ttotal2min>ttot2min) ttotal2min=ttot2min
               if (ttotal2<ttot2) ttotal2=ttot2
            end do 
            J12min=max(0,jtotal2min/2)
            J12max=min(nhom123max+2,jtotal2/2)
            t12min=max(0,ttotal2min/2)
            t12max=min(2,ttotal2/2)
         else
            ipimax=1
            ipimin=-1
            J12min=0
            J12max=nhom1234max+2
            t12min=0
            t12max=2
         endif
         allocate(ad4ind(ist3dim,ist1dim,J12min:J12max,t12min:t12max))
         ind=0
         do J1234=J12min,J12max
            do iT1234=t12min,t12max
               do ipi=ipimax,ipimin,-2
                  do ii=1,ist3dim
                     ii2=ist3(1,ii)
                     ii1=ist3(2,ii)
                     J123=ist3(3,ii)
                     iT123=ist3(4,ii)
                     ia=ist2(1,ii2)
                     ib=ist2(2,ii2)
                     na=ist1(1,ia)
                     la=ist1(2,ia)
                     nb=ist1(1,ib)
                     lb=ist1(2,ib)
                     nc=ist1(1,ii1)
                     lc=ist1(2,ii1)
                     do ij=1,ist1dim
                        nd=ist1(1,ij)
                        ld=ist1(2,ij)
                        j2d=ist1(3,ij)
                        if ((-1)**(la+lb+lc+ld)/=ipi) cycle
                        if (2*na+la+2*nb+lb+2*nc+lc+2*nd+ld
     +                       >nhom1234max) cycle
                        if (abs(J123-j2d)>2*J1234.or.(J123+j2d)<2*J1234) 
     +                       cycle
                        if (abs(iT123-1)>2*iT1234.or.(iT123+1)<2*iT1234) 
     +                       cycle
                        ind=ind+1
                        ad4ind(ii,ij,J1234,iT1234)=ind
                     end do
                  end do
               end do
            end do
         end do
         ist4dim=ind
         allocate(ist4(4,ist4dim))
         ind=0
         do J1234=J12min,J12max
            do iT1234=t12min,t12max
               do ipi=ipimax,ipimin,-2
                  do ii=1,ist3dim
                     ii2=ist3(1,ii)
                     ii1=ist3(2,ii)
                     J123=ist3(3,ii)
                     iT123=ist3(4,ii)
                     ia=ist2(1,ii2)
                     ib=ist2(2,ii2)
                     na=ist1(1,ia)
                     la=ist1(2,ia)
                     nb=ist1(1,ib)
                     lb=ist1(2,ib)
                     nc=ist1(1,ii1)
                     lc=ist1(2,ii1)
                     do ij=1,ist1dim
                        nd=ist1(1,ij)
                        ld=ist1(2,ij)
                        j2d=ist1(3,ij)
                        if ((-1)**(la+lb+lc+ld)/=ipi) cycle
                        if (2*na+la+2*nb+lb+2*nc+lc+2*nd+ld
     +                       >nhom1234max) cycle
                        if (abs(J123-j2d)>2*J1234.or.(J123+j2d)<2*J1234) 
     +                       cycle
                        if (abs(iT123-1)>2*iT1234.or.(iT123+1)<2*iT1234) 
     +                       cycle
                        ind=ind+1
                        ist4(1,ind)=ii
                        ist4(2,ind)=ij
                        ist4(3,ind)=J1234
                        ist4(4,ind)=iT1234
                     end do
                  end do
               end do
            end do
         end do
         if (ind==ist4dim) then
            if (iproc==0)
     $           write(*,"(' Number of four-nucleon states:',i8)")
     $           ist4dim
         else
            print *,
     +           '*** error in four-body state generation:ind,ist4dim=',
     +           ind,ist4dim
            stop
         endif   
      endif
      end

      subroutine obdtbdcalc_hash(snsp,sndp,dnuc)
      use initst
      use finast
      use intrface
      use paramdef
      use obdens
      use nodeinfo
      use occ
      use hash_tables
      use cleb_coef
      implicit none
      include 'mpif.h'
      logical,intent(IN) :: snsp,sndp,dnuc
      integer :: mxsps,jtotal2,jtotal2min,kii,kff,ii,jtot2,jtot2min
      integer :: nuc_min,nuc_max,nuc_min_2,nuc_max_2,cr_st_min,
     +     cr_st_min_2,cr_st_max,cr_st_max_2
      integer :: ia_1,an_st_1,cr_1,sp_st_an,sp_st_cr,ia,i1,i2,
     +     p_state,n_state,cr_st_1
      integer,allocatable :: occi(:),occf(:),occim1(:),occim2(:)
      integer :: phase,phase1,N_HO_i,wd,ibit
      integer(8),allocatable :: intbas(:)
      integer,allocatable :: locz(:)
      integer :: jtotali2,itotali2,jtotalf2,itotalf2,jtrans,jtrans2
      integer :: n1,l1,j1,n2,l2,j2,ierr,i1start
      real(kind(0.d0)) :: cleb,clb1,clebd,tbx(3)
      real(kind=kind(0.d0)),allocatable :: tdJtemp(:,:,:,:,:)
      real(kind=kind(0.d0)),allocatable :: t2bdtemp(:),t2bdtemp_rec(:),
     +     clebjj(:,:,:,:,:),clebj12j34(:,:,:,:,:)
      integer :: jf2,mf2,ji2,mi2,iphase,interm_energy
      integer :: ia_2,an_st_2,cr_st_2,sp_st_an_2,sp_st_cr_2,
     $     cr_st_min_1,cr_st_max_1,cr_st_min2b_2,cr_st_max2b_2
      integer :: i3,j12,it12,j34,it34,nonebody,nonebody2,nonebody3,
     $     nonebody4,l3,l4
      integer(8) :: ih3ind,ih3oibuf,iix,ibuf8,one8,i8
      integer :: ibuf,i
      integer :: Jtr_mi,Jtr_ma,mT_mi,mT_ma,mTt,in_mi,in_ma,nhom12max,pi,
     $     J12ma
      integer :: an_st_3,delta_mj,delta_mt,J123,iT123,
     $     beta,j_beta,m_beta,mt_beta,ia_3,ia_4,an_st_4,J1234,iT1234
      real(kind=kind(0.d0)),allocatable :: adtemp(:,:,:)

      integer :: !nblock=8,  !buffer_size=128,multf=100,
     $     nproc_blockgroup,iproc_blockgroup,wri=1,
     $     ist2dim_block,i_block,iproc_block,i_local,dest,ibl !,final_mes
      integer :: kf_in,ki_in,Jtr,i_loc,mT,in
      real(kind(0.d0)) :: aloc_size,aloc_size_tot
      logical :: re_err,clsun
      character(len=80) :: line
cc      type tbd_update
cc      sequence
cc      integer :: kf,ki,Jtr,i_loc,mT,in
cc      real(kind(0.d0)) :: val
cc      end type tbd_update
cc      type(tbd_update),allocatable :: tbd_upda(:,:),buffer_bsend(:)
cc      integer :: blockcounts(0:1),offsets(0:1),oldtypes(0:1),extent,
cc     $     tbd_update_type
      integer :: stat(MPI_STATUS_SIZE),tag,source,count !,buff_size_bs
cc      logical :: flag
cc      logical,allocatable :: finished_loop(:)
cc      integer,allocatable :: buffer_count(:)

      print *,' obdcalc_hash entered'

      if (irestart==-1) wri=-1
      print *,' irestart,wri=',irestart,wri

      mxsps=nbit*mxnwdi

      print *,' mxsps=',mxsps

      jtotal2=(jt2i(ki)+jt2f(kf))/2
      jtotal2min=iabs(jt2i(ki)-jt2f(kf))/2
      do kii=ki,ki+nki-1
         do kff=kf,kf+nkf-1
            jtot2min=iabs(jt2i(kii)-jt2f(kff))/2
            jtot2=(jt2i(kii)+jt2f(kff))/2
            if (jtotal2min>jtot2min) jtotal2min=jtot2min
            if (jtotal2<jtot2) jtotal2=jtot2
         end do
      end do 

      jtotal2=min(jtotal2,jtotal2max)

      print *,' jtotal2min,jtotal2:',jtotal2min,jtotal2

      print *,' naspsi=',naspsi
      print *,' mxspsi=',mxspsi
c      do i1=1,2*mxspsi
c         print *,' #, n,l,j,m,mt=',i1,n_spi(i1),l_spi(i1),j2_spi(i1),
c     $        m2_spi(i1),mt2_spi(i1)
c      end do

      if (iproc==0) then
         write(iunitvout,*)
         write(iunitvout,7768) ist1dim
 7768    format(' number of single-nucleon states =',i4)
         do ii=1,ist1dim
            write(iunitvout,7770) ii,(ist1(ia,ii),ia=1,3)
 7770       format(' #',i4,'  n=',i3,'  l=',i3,'  j=',i2,'/2')
         end do
      endif   

      call num_of_3j_types_init(1)
      call cleb_init(1,1,2*nhomf+1,-2*nhomf-1,2*nhomf+1,
     $     1,2*nhomi+1,-2*nhomi-1,2*nhomi+1,2*jtotal2min,2*jtotal2)

      print *,' ist1dim',ist1dim

      allocate(tdJp(ist1dim,ist1dim,jtotal2min:jtotal2,
     +     kf:kf+nkf-1,ki:ki+nki-1))          
      allocate(tdJn(ist1dim,ist1dim,jtotal2min:jtotal2,
     +     kf:kf+nkf-1,ki:ki+nki-1))          

      tdJp=0.d0
      tdJn=0.d0 

      if (twobdcal) then
         allocate(clebjj(0:(2*nhomf+1)/2,-(2*nhomf+1):2*nhomf+1,
     +        0:(2*nhomi+1)/2,-(2*nhomi+1):2*nhomi+1,0:nhomi+nhomf+1))
         clebjj=0.d0
         allocate(clebj12j34(0:2*nhomf+1,-(2*nhomf+1):2*nhomf+1,
     +        jtotal2min:jtotal2,0:2*nhomi+1,-(2*nhomi+1):2*nhomi+1))
         clebj12j34=0.d0
         
         do jf2=0,2*nhomf+1
            do mf2=-jf2,jf2
               do ji2=0,2*nhomi+1
                  do mi2=-ji2,ji2
                     do jtrans=iabs(jf2-ji2),(jf2+ji2)
                        if (mod(jf2,2)==1.and.mod(ji2,2)==1
     +                       .and.mod(jtrans,2)==0.and.
     +                       mod(mf2,2)/=0.and.mod(mi2,2)/=0) then
                           clebjj(jf2/2,mf2,ji2/2,mi2,jtrans/2)=
     +                          clebd(jf2,mf2,ji2,mi2,jtrans,mf2+mi2)
                        endif   
                        if (jtrans>=jtotal2min.and.jtrans<=jtotal2)
     +                       then 
                           clebj12j34(jf2,mf2,jtrans,ji2,mi2)=
     +              clebd(2*jf2,2*mf2,2*jtrans,2*(mi2-mf2),2*ji2,2*mi2)
                        endif   
                     end do 
                  end do
               end do             
            end do
         end do 
         print *,' ist2dim',ist2dim
c         allocate(t2bd(ist2dim,ist2dim,
c     +        jtotal2min:jtotal2,-1:1,kf:kf+nkf-1,ki:ki+nki-1))          
c         t2bd=0.d0

c*** MPI tbd distribution
         nblock=min(nblock,nproc)
         print *,' nblock=',nblock
         if (mod(ist2dim,nblock)==0) then
            ist2dim_block=ist2dim/nblock
         else
            ist2dim_block=ist2dim/nblock+1
         endif
         print *,' ist2dim_block=',ist2dim_block
         iproc_block=mod(iproc,nblock)+1
         print *,' iproc,iproc_block=',iproc,iproc_block
c*** MPI derived data type
c         offsets(0)=0
c         oldtypes(0)=MPI_INTEGER
c         blockcounts(0)=6
c         call MPI_TYPE_EXTENT(MPI_INTEGER,extent,ierr)
c         offsets(1)=6*extent
c         oldtypes(1)=MPI_REAL8
c         blockcounts(1)=1
c         call MPI_TYPE_STRUCT(2,blockcounts,offsets,oldtypes,
c     $        tbd_update_type,ierr)
c         call MPI_TYPE_COMMIT(tbd_update_type,ierr)
cc         allocate(tbd_upda(buffer_size,0:nproc-1))
cc         allocate(buffer_count(0:nproc-1))
cc         buffer_count=0
cc         allocate(finished_loop(0:nproc-1))
cc         finished_loop=.false.
c         buff_size_bs=nproc*buffer_size*multf
c         allocate(buffer_bsend(buff_size_bs))
c         call MPI_BUFFER_ATTACH(buffer_bsend,32*buff_size_bs,ierr)

         nhom12max=max(nhom12i,nhom12f)
         allocate(tbd(kf:kf+nkf-1,ki:ki+nki-1))
         aloc_size=0.d0
         do kii=ki,ki+nki-1
            do kff=kf,kf+nkf-1
               Jtr_mi=abs(jt2i(kii)-jt2f(kff))/2
               Jtr_ma=min((jt2i(kii)+jt2f(kff))/2,jtotal2max)
               tbd(kff,kii)%Jtr_min=Jtr_mi
               tbd(kff,kii)%Jtr_max=Jtr_ma
               allocate(tbd(kff,kii)%Jtr(Jtr_mi:Jtr_ma))
               do jtrans=Jtr_mi,Jtr_ma
cc                  tbd(kff,kii)%Jtr(jtrans)%dim_fi=ist2dim
cc                  allocate(tbd(kff,kii)%Jtr(jtrans)%fi(ist2dim))
                  tbd(kff,kii)%Jtr(jtrans)%dim_fi=ist2dim_block
                  allocate(tbd(kff,kii)%Jtr(jtrans)%fi(ist2dim_block))
                  do ii=1,ist2dim
                     i_block=(ii-1)/ist2dim_block+1
                     if (i_block/=iproc_block) cycle
                     if (mod(ii,ist2dim_block)==0) then
                        i_local=ist2dim_block
                     else
                        i_local=mod(ii,ist2dim_block)
                     endif
                     i1=ist2(1,ii)
                     i2=ist2(2,ii)
                     j12=ist2(3,ii)
                     it12=ist2(4,ii)

                     mT_mi=max(-abs(it12),(mttotalf-mttotali)/2-1)
                     mT_ma=min(abs(it12),(mttotalf-mttotali)/2+1)
cc                     tbd(kff,kii)%Jtr(jtrans)%fi(ii)%mT_min=mT_mi
cc                     tbd(kff,kii)%Jtr(jtrans)%fi(ii)%mT_max=mT_ma
                     tbd(kff,kii)%Jtr(jtrans)%fi(i_local)%mT_min=mT_mi
                     tbd(kff,kii)%Jtr(jtrans)%fi(i_local)%mT_max=mT_ma
                     if (mT_mi>mT_ma) cycle
                     l1=ist1(2,i1)
                     l2=ist1(2,i2)
                     pi=(-1)**(l1+l2+iparityi+iparityf)
                     if (pi/=(-1)**nhom12max) then
                        J12ma=nhom12max
                     else
                        J12ma=nhom12max+1
                     endif
                     if (abs(jtrans-j12)>ist2_Jst_J12max
     $                    .or.abs(jtrans-j12)<ist2_Jst_J12min
     $                    .or.min(jtrans+j12,J12ma)+1>ist2_Jst_J12max+1
     $                    .or.min(jtrans+j12,J12ma)+1<ist2_Jst_J12min)
     $                    then
                        in_mi=-1
                        in_ma=-1
                     else
                        in_mi=ist2_Jst((pi+1)/2,abs(jtrans-j12))
                        in_ma=ist2_Jst((pi+1)/2,min(jtrans+j12,J12ma)+1)
     $                       -1
                     endif                     
c                     in_mi=ist2_Jst((pi+1)/2,abs(jtrans-j12))
c                     in_ma=ist2_Jst((pi+1)/2,min(jtrans+j12,J12ma)+1)-1
cc                     allocate(tbd(kff,kii)%Jtr(jtrans)%fi(ii)
cc     $                    %mT(mT_mi:mT_ma))
                     allocate(tbd(kff,kii)%Jtr(jtrans)%fi(i_local)
     $                    %mT(mT_mi:mT_ma))
                     do mTt=mT_mi,mT_ma 
cc                        tbd(kff,kii)%Jtr(jtrans)%fi(ii)%mT(mTt)
cc     $                       %in_min=in_mi
cc                        tbd(kff,kii)%Jtr(jtrans)%fi(ii)%mT(mTt)
cc     $                       %in_max=in_ma
cc                        allocate(tbd(kff,kii)%Jtr(jtrans)%fi(ii)%mT(mTt)
cc     $                       %in(in_mi:in_ma))
cc                        tbd(kff,kii)%Jtr(jtrans)%fi(ii)%mT(mTt)
cc     $                       %in=0.d0
                        tbd(kff,kii)%Jtr(jtrans)%fi(i_local)%mT(mTt)
     $                       %in_min=in_mi
                        tbd(kff,kii)%Jtr(jtrans)%fi(i_local)%mT(mTt)
     $                       %in_max=in_ma
                        allocate(tbd(kff,kii)%Jtr(jtrans)%fi(i_local)
     $                       %mT(mTt)%in(in_mi:in_ma))
                        tbd(kff,kii)%Jtr(jtrans)%fi(i_local)%mT(mTt)
     $                       %in=0.d0
                        aloc_size=aloc_size+real(in_ma-in_mi+1,
     $                       kind(0.d0))

                     end do
                  end do
               end do
            end do
         end do
c         print *,' tbd allocated for iproc=',iproc
         print *,' iproc=',iproc,' allocated tbd of size=',aloc_size
         call MPI_Reduce(aloc_size,aloc_size_tot,1,MPI_REAL8,MPI_SUM,0,
     $        icomm,ierr)
         if (iproc==0) print *,' total size of single allocated tbd=',
     $        aloc_size_tot*real(nblock,kind(0.d0))
     $        /real(nproc,kind(0.d0))
         if (iproc==0) print *,' total size of tbd in MB=',
     $        8.d0*aloc_size_tot*real(nblock,kind(0.d0))
     $        /real(nproc,kind(0.d0))/real(1024**2,kind(0.d0))

         if (nucleonsi>2) allocate(occim2(nucleonsi-2))
         nproc_blockgroup=nproc/nblock
         if (iproc_block<=mod(nproc,nblock))
     $        nproc_blockgroup=nproc_blockgroup+1
         iproc_blockgroup=iproc/nblock
         print *,' N,ipr,ipr_gr,N_gr=',
     $        nproc,iproc,iproc_blockgroup,nproc_blockgroup
      else
         iproc_block=1
         nproc_blockgroup=nproc
         iproc_blockgroup=iproc
      endif

      if (threebdcal.and.nucleonsi==3) then
         allocate(ad_cl(ist3dim,1,ki:ki+nki-1))
         ad_cl=0.d0
         if (mjtotali>0) then
            m_beta=1
         else
            m_beta=-1
         endif
         if (mttotali>0) then
            mt_beta=1
         else
            mt_beta=-1
         endif
         allocate(ad_cl_2(ist2dim,ist1dim,ki:ki+nki-1))
         ad_cl_2=0.d0
      elseif (threebdcal.and.nucleonsi==4) then
         if (mjtotali>=0) then
            m_beta=1
         else
            m_beta=-1
         endif
         if (mttotali>=0) then
            mt_beta=1
         else
            mt_beta=-1
         endif
         allocate(ad_cl(ist3dim,ist1dim,ki:ki+nki-1))
         ad_cl=0.d0
         allocate(ad_cl_2(ist4dim,1,ki:ki+nki-1))
         ad_cl_2=0.d0
      endif

      if (snsp.or.sndp) then
         nuc_min=1
         nuc_max=nprotonsi
         cr_st_min=1
c         cr_st_max=mxsps
         cr_st_max=naspsi
         nuc_min_2=nprotonsi+1
         nuc_max_2=nucleonsi
c         cr_st_max_2=2*mxsps
         cr_st_min_2=mxsps+1
         cr_st_max_2=mxsps+naspsi
      elseif (dnuc) then
         if (nprotonsi==nprotonsf+1) then
            nuc_min=1
            nuc_max=0
            cr_st_min=1
            cr_st_max=0
            nuc_min_2=1
            nuc_max_2=nprotonsi
            cr_st_min_2=mxsps+1
c            cr_st_max_2=2*mxsps
            cr_st_max_2=mxsps+naspsi
         else
            nuc_min=nprotonsi+1
            nuc_max=nucleonsi
            cr_st_min=1
c            cr_st_max=mxsps
            cr_st_max=naspsi
            nuc_min_2=1
            nuc_max_2=0
            cr_st_min_2=1
            cr_st_max_2=0
         endif
      else
         print *,'***error in obdcalc_hash'
         stop
      endif

      allocate(intbas(2*mxnwdi))
      allocate(locz(2*mxsps))
      allocate(occi(nucleonsi))
      allocate(occim1(nucleonsi-1))
      allocate(occf(nucleonsf))

c*** main loop for transition density calculation ***********     
      if (irestart==1.or.irestart==2) then

         write(line,'(i8)') iproc
         line=adjustl(line)
         line='trdens'//trim(line)//ext1
         open(10,file=trim(line),access='sequential',form='unformatted',
     $        status='old')

c         select case(iproc)
c         case(0:9)
c            open(10, file='trdens'//achar(iproc+48)//ext1,
c     +           access='sequential',
c     +           form='unformatted',status='old')
c         case(10:99)
c            open(10, file='trdens'//achar(iproc/10+48)
c     +           //achar(mod(iproc,10)+48)//ext1,            
c     +           access='sequential',
c     +           form='unformatted',status='old')
c         case(100:999)
c            open(10, file='trdens'//achar(iproc/100+48)
c     +           //achar(iproc/10+48)
c     +           //achar(mod(iproc,10)+48)//ext1,            
c     +           access='sequential',
c     +           form='unformatted',status='old')
c         case(1000:9999)
c            open(10, file='trdens'//achar(iproc/1000+48)
c     $           //achar(iproc/100+48)
c     +           //achar(iproc/10+48)
c     +           //achar(mod(iproc,10)+48)//ext1,            
c     +           access='sequential',
c     +           form='unformatted',status='old')
c         case default
c            print *,' invalid number of processors'
c            stop
c         end select   

         read(10,err=3457,end=3457) i1start,tdJp,tdJn
         if (twobdcal)  then   !read(10,err=3457,end=3457) tbd
            do kii=ki,ki+nki-1
               do kff=kf,kf+nkf-1
                  call read_de_Jtr(10,tbd(kff,kii),re_err)
                  if (re_err) goto 3457
               end do
            end do
         endif
         close(10)
         goto 3456
 3457    continue
c         stop '*** error in reading'
         inquire(unit=10,opened=clsun)
         if (clsun) close(10)

         write(line,'(i8)') iproc
         line=adjustl(line)
         line='trdens'//trim(line)//ext2
         open(11,file=trim(line),access='sequential',form='unformatted',
     $        status='old')

c         select case(iproc)
c         case(0:9)
c            open(11, file='trdens'//achar(iproc+48)//ext2,
c     +           access='sequential',
c     +           form='unformatted',status='old')
c         case(10:99)
c            open(11, file='trdens'//achar(iproc/10+48)
c     +           //achar(mod(iproc,10)+48)//ext2,            
c     +           access='sequential',
c     +           form='unformatted',status='old')
c         case(100:999)
c            open(11, file='trdens'//achar(iproc/100+48)
c     +           //achar(iproc/10+48)
c     +           //achar(mod(iproc,10)+48)//ext2,            
c     +           access='sequential',
c     +           form='unformatted',status='old')
c         case(1000:9999)
c            open(11, file='trdens'//achar(iproc/1000+48)
c     $           //achar(iproc/100+48)
c     +           //achar(iproc/10+48)
c     +           //achar(mod(iproc,10)+48)//ext2,            
c     +           access='sequential',
c     +           form='unformatted',status='old')
c         case default
c            print *,' invalid number of processors'
c            stop
c         end select
         read(11) i1start,tdJp,tdJn
         if (twobdcal)  then   !read(11,err=3457,end=3457) tbd
            do kii=ki,ki+nki-1
               do kff=kf,kf+nkf-1
                  call read_de_Jtr(11,tbd(kff,kii),re_err)
               end do
            end do
         endif
         close(11)            
      endif   
 3456 continue
      
c      call MPI_Barrier(icomm,ierr)
c      do i1 = iproc+1, nsdi, nproc
      do i1 = iproc_blockgroup+1, nsdi, nproc_blockgroup

         if ((irestart==1.or.irestart==2).and.i1<i1start) cycle

         if (mod(i1-iproc_blockgroup-1,
     $        (3000000/nproc_blockgroup)*nproc_blockgroup)==0
     $        .and.i1>iproc_blockgroup+1) then
            print *, '#',iproc,' doing i1=',i1

            if (wri==1) then
               
               write(line,'(i8)') iproc
               line=adjustl(line)
               line='trdens'//trim(line)//'.tmp'
               open(10,file=trim(line),access='sequential',
     $              form='unformatted',status='unknown')

c               select case(iproc)
c               case(0:9)
c                  open(10, file='trdens'//achar(iproc+48)//'.tmp',
c     +                 access='sequential',
c     +                 form='unformatted',status='unknown')
c               case(10:99)
c                  open(10, file='trdens'//achar(iproc/10+48)
c     +                 //achar(mod(iproc,10)+48)//'.tmp',            
c     +                 access='sequential',
c     +                 form='unformatted',status='unknown')
c               case(100:999)
c                  open(10, file='trdens'//achar(iproc/100+48)
c     +                 //achar(iproc/10+48)
c     +                 //achar(mod(iproc,10)+48)//'.tmp',            
c     +                 access='sequential',
c     +                 form='unformatted',status='unknown')
c               case(1000:9999)
c                  open(10, file='trdens'//achar(iproc/1000+48)
c     $                 //achar(iproc/100+48)
c     +                 //achar(iproc/10+48)
c     +                 //achar(mod(iproc,10)+48)//'.tmp',            
c     +                 access='sequential',
c     +                 form='unformatted',status='unknown')
c               case default
c                  print *,' invalid number of processors'
c                  stop
c               end select

               wri=0   
            elseif (wri==0) then
               
               write(line,'(i8)') iproc
               line=adjustl(line)
               line='trdens'//trim(line)//'.bak'
               open(10,file=trim(line),access='sequential',
     $              form='unformatted',status='unknown')

c               select case(iproc)
c               case(0:9)
c                  open(10, file='trdens'//achar(iproc+48)//'.bak',
c     +                 access='sequential',
c     +                 form='unformatted',status='unknown')
c               case(10:99)
c                  open(10, file='trdens'//achar(iproc/10+48)
c     +                 //achar(mod(iproc,10)+48)//'.bak',            
c     +                 access='sequential',
c     +                 form='unformatted',status='unknown')
c               case(100:999)
c                  open(10, file='trdens'//achar(iproc/100+48)
c     +                 //achar(iproc/10+48)
c     +                 //achar(mod(iproc,10)+48)//'.bak',            
c     +                 access='sequential',
c     +                 form='unformatted',status='unknown')
c               case(1000:9999)
c                  open(10, file='trdens'//achar(iproc/1000+48)
c     $                 //achar(iproc/100+48)
c     +                 //achar(iproc/10+48)
c     +                 //achar(mod(iproc,10)+48)//'.bak',            
c     +                 access='sequential',
c     +                 form='unformatted',status='unknown')
c               case default
c                  print *,' invalid number of processors'
c                  stop
c               end select
               wri=1   
            endif
            if (wri==0.or.wri==1) then
               write(10) i1,tdJp,tdJn
               if (twobdcal) then !write(10) tbd
                  do kii=ki,ki+nki-1
                     do kff=kf,kf+nkf-1
                        call write_de_Jtr(10,tbd(kff,kii))
                     end do
                  end do
               endif
               close(10)
            endif
         endif  

         if (majortot==2) then
            occi(:)=iloci(:,i1)
         elseif (majortot==3) then
            p_state=I_state_i(1,i1)
            n_state=I_state_i(2,i1)
            occi(1:nprotonsi)=occ_p_i(:,p_state)
            do ii=1,nneutrnsi
               occi(nprotonsi+ii)=occ_n_i(ii,n_state)+mxnwdi*nbit
            end do
         endif
         
         N_HO_i=0
         do ia_1=1,nucleonsi
            an_st_1=occi(ia_1)
            N_HO_i=N_HO_i+2*n_spi(an_st_1)+l_spi(an_st_1)
         end do

         if (.not.twobdcal.or.iproc_block==1) then 
         do ia_1=nuc_min,nuc_max
            an_st_1=occi(ia_1)
            occim1(1:ia_1-1)=occi(1:ia_1-1)
            occim1(ia_1:nucleonsi-1)=occi(ia_1+1:nucleonsi)
            phase1=(-1)**(nucleonsi-ia_1
     +           +(j2_spi(an_st_1)+m2_spi(an_st_1))/2)
            sp_st_an=iobsind(an_st_1)

            locz(1:occim1(1))=0
            do cr_1=2,nucleonsi-1
               locz(occim1(cr_1-1)+1:occim1(cr_1))=cr_1-1
            end do
            locz(occim1(nucleonsi-1)+1:2*mxsps)=nucleonsi-1

            intbas=0
            do cr_1=1,nucleonsi-1
               wd=(occim1(cr_1)-1)/nbit+1
               ibit=mod(occim1(cr_1)-1,nbit)
c               print *,' i,wd,ibit=',i,wd,ibit
               intbas(wd)=ibset(intbas(wd),ibit)
            end do

            do cr_st_1=cr_st_min,cr_st_max
               if (cr_st_1==an_st_1.and.snsp) then
                  phase=phase1*((-1)**(nucleonsi-1-locz(cr_st_1)))
                  i2=i1
                  call one_body(cr_st_1,an_st_1,
     +                 tdJp(sp_st_an,sp_st_an,:,:,:),jtotal2min,
     +                 jtotal2,ki,nki,kf,nkf)
                  cycle
               endif
               wd=(cr_st_1-1)/nbit+1
               ibit=mod(cr_st_1-1,nbit)
               if (btest(intbas(wd),ibit)) cycle
               if (m2_spi(cr_st_1)-m2_spi(an_st_1)+mjtotali
     +              /=mjtotalf) cycle
               if (mod(iparityf+l_spi(cr_st_1)+l_spi(an_st_1)+iparityi,
     +              2)/=0) cycle
               if (2*n_spi(cr_st_1)+l_spi(cr_st_1)
     +              -2*n_spi(an_st_1)-l_spi(an_st_1)+N_HO_i>nhwf) exit
               occf(1:locz(cr_st_1))=occim1(1:locz(cr_st_1))
               occf(locz(cr_st_1)+1)=cr_st_1
               occf(locz(cr_st_1)+2:nucleonsf)=
     +              occim1(locz(cr_st_1)+1:nucleonsi-1)
               phase=phase1*((-1)**(nucleonsi-1-locz(cr_st_1)))
cc               if (twobdcal) call check_incoming
               call get_state_index(nucleonsf,occf,1,i2)
               if (i2==-1) cycle
               sp_st_cr=iobsind(cr_st_1)
               call one_body(cr_st_1,an_st_1,
     +              tdJp(sp_st_cr,sp_st_an,:,:,:),jtotal2min,jtotal2,
     +              ki,nki,kf,nkf)

            end do
         end do

         do ia_1=nuc_min_2,nuc_max_2
            an_st_1=occi(ia_1)
            occim1(1:ia_1-1)=occi(1:ia_1-1)
            occim1(ia_1:nucleonsi-1)=occi(ia_1+1:nucleonsi)
            phase1=(-1)**(nucleonsi-ia_1
     +           +(j2_spi(an_st_1)+m2_spi(an_st_1))/2)
            sp_st_an=iobsind(an_st_1)

            locz(1:occim1(1))=0
            do cr_1=2,nucleonsi-1
               locz(occim1(cr_1-1)+1:occim1(cr_1))=cr_1-1
            end do
            locz(occim1(nucleonsi-1)+1:2*mxsps)=nucleonsi-1

            intbas=0
            do cr_1=1,nucleonsi-1
               wd=(occim1(cr_1)-1)/nbit+1
               ibit=mod(occim1(cr_1)-1,nbit)
               intbas(wd)=ibset(intbas(wd),ibit)
            end do

            do cr_st_1=cr_st_min_2,cr_st_max_2
               if (cr_st_1==an_st_1.and.snsp) then
                  phase=phase1*((-1)**(nucleonsi-1-locz(cr_st_1)))
                  i2=i1
                  call one_body(cr_st_1,an_st_1,
     +                 tdJn(sp_st_an,sp_st_an,:,:,:),jtotal2min,
     +                 jtotal2,ki,nki,kf,nkf)
                  cycle
               endif
               wd=(cr_st_1-1)/nbit+1
               ibit=mod(cr_st_1-1,nbit)
               if (btest(intbas(wd),ibit)) cycle
               if (m2_spi(cr_st_1)-m2_spi(an_st_1)+mjtotali
     +              /=mjtotalf) cycle
               if (mod(iparityf+l_spi(cr_st_1)+l_spi(an_st_1)+iparityi,
     +              2)/=0) cycle
               if (2*n_spi(cr_st_1)+l_spi(cr_st_1)
     +              -2*n_spi(an_st_1)-l_spi(an_st_1)+N_HO_i>nhwf) exit
               occf(1:locz(cr_st_1))=occim1(1:locz(cr_st_1))
               occf(locz(cr_st_1)+1)=cr_st_1
               occf(locz(cr_st_1)+2:nucleonsf)=
     +              occim1(locz(cr_st_1)+1:nucleonsi-1)
               phase=phase1*((-1)**(nucleonsi-1-locz(cr_st_1)))
cc               if (twobdcal) call check_incoming
               call get_state_index(nucleonsf,occf,1,i2)
               if (i2==-1) cycle
               sp_st_cr=iobsind(cr_st_1)
               call one_body(cr_st_1,an_st_1,
     +              tdJn(sp_st_cr,sp_st_an,:,:,:),jtotal2min,jtotal2,
     +              ki,nki,kf,nkf)

            end do
         end do
         endif

         if (twobdcal) then
            do ia_1=1,nucleonsi-1
               an_st_1=occi(ia_1)
               do ia_2=ia_1+1,nucleonsi
                  an_st_2=occi(ia_2)

                  select case(mttotalf-mttotali+mt2_spi(an_st_1)
     $                 +mt2_spi(an_st_2))
                  case(2)
                     cr_st_min_1=1
                     cr_st_min2b_2=2
                     cr_st_max_1=naspsi-1
                     cr_st_max2b_2=naspsi
                  case(-2)
                     cr_st_min_1=mxsps+1
                     cr_st_min2b_2=mxsps+2
                     cr_st_max_1=mxsps+naspsi-1
                     cr_st_max2b_2=mxsps+naspsi
                  case(0)
                     cr_st_min_1=1
                     cr_st_min2b_2=mxsps+1
                     cr_st_max_1=naspsi
                     cr_st_max2b_2=mxsps+naspsi
                  case default
                     cycle
                  end select

                  if (nucleonsi>2) then
                     occim2(1:ia_1-1)=occi(1:ia_1-1)
                     occim2(ia_1:ia_2-2)=occi(ia_1+1:ia_2-1)
                     occim2(ia_2-1:nucleonsi-2)=occi(ia_2+1:nucleonsi)

                     locz(1:occim2(1))=0
                     do cr_1=2,nucleonsi-2
                        locz(occim2(cr_1-1)+1:occim2(cr_1))=cr_1-1
                     end do
                     locz(occim2(nucleonsi-2)+1:2*mxsps)=nucleonsi-2

                     intbas=0
                     do cr_1=1,nucleonsi-2
                        wd=(occim2(cr_1)-1)/nbit+1
                        ibit=mod(occim2(cr_1)-1,nbit)
                        intbas(wd)=ibset(intbas(wd),ibit)
                     end do
                  else
                     intbas=0
                  endif
                  interm_energy=-2*n_spi(an_st_1)-l_spi(an_st_1)
     $                 -2*n_spi(an_st_2)-l_spi(an_st_2)+N_HO_i

cc                  call check_incoming
                  do cr_st_1=cr_st_min_1,cr_st_max_1

                     if (2*n_spi(cr_st_1)+l_spi(cr_st_1)
     +                    +interm_energy>nhwf) exit
                     wd=(cr_st_1-1)/nbit+1
                     ibit=mod(cr_st_1-1,nbit)

                     if (btest(intbas(wd),ibit)) cycle
                     do cr_st_2=max(cr_st_min2b_2,cr_st_1+1),
     $                    cr_st_max2b_2

                        if (2*n_spi(cr_st_1)+l_spi(cr_st_1)
     $                       +2*n_spi(cr_st_2)+l_spi(cr_st_2)
     +                       +interm_energy>nhwf) exit

                        if (m2_spi(cr_st_1)+m2_spi(cr_st_2)
     $                       -m2_spi(an_st_1)-m2_spi(an_st_2)+mjtotali
     +                       /=mjtotalf) cycle

                        if (mod(iparityf+l_spi(cr_st_1)+l_spi(an_st_1)
     $                       +l_spi(cr_st_2)+l_spi(an_st_2)+iparityi,
     +                       2)/=0) cycle
                        wd=(cr_st_2-1)/nbit+1
                        ibit=mod(cr_st_2-1,nbit)
                        
                        if (btest(intbas(wd),ibit)) cycle
                        if (nucleonsi>2) then
                           occf(1:locz(cr_st_1))=occim2(1:locz(cr_st_1))
                           occf(locz(cr_st_1)+1)=cr_st_1
                           occf(locz(cr_st_1)+2:locz(cr_st_2)+1)=
     $                          occim2(locz(cr_st_1)+1:locz(cr_st_2))
                           occf(locz(cr_st_2)+2)=cr_st_2
                           occf(locz(cr_st_2)+3:nucleonsf)=
     +                          occim2(locz(cr_st_2)+1:nucleonsi-2)
                        else
                            occf(1)=cr_st_1
                            occf(2)=cr_st_2
                         endif
                         
cc                        call check_incoming
                        call get_state_index(nucleonsf,occf,1,i2)

                        if (i2==-1) cycle

                        if (nucleonsi>2) then
                           iphase=(-1)**(ia_1+ia_2+locz(cr_st_1)
     $                          +locz(cr_st_2)+1) !exchange of an_1 and an_2 
                        else
                           iphase=(-1)**(ia_1+ia_2+1) !exchange of an_1 and an_2 
                        endif
c                                !   or -1 in overall formula
                        call two_body(cr_st_2,cr_st_1,an_st_2,an_st_1)

                     end do
                  end do
               end do
            end do
c            call check_incoming
c            if (i1+nproc>nsdi) then
c               print *,' final message: iproc,i1=',iproc,i1
c               tag=777
c               final_mes=-777
c               do ii=0,nproc-1
c                  if (ii==iproc) cycle
c                  call MPI_Send(final_mes,1,MPI_INTEGER,ii,tag,icomm,
c     $                 ierr)
c               end do
c            endif

         endif

         if (threebdcal) then
            if (nucleonsi==3) then
               an_st_1=occi(1)
               an_st_2=occi(2)
               an_st_3=occi(3)
               delta_mj=mjtotali
               delta_mt=mttotali
               kff=1      
               phase=-1    ! (-1)**(A-i1+A-i2+A-i3) old: +1
               do kii=ki,ki+nki-1
                  J123=jt2i(kii)
                  iT123=it2i(kii)
                  call three_nucleon(an_st_1,an_st_2,an_st_3,1)
                  call three_nucleon(an_st_1,an_st_3,an_st_2,-1)
                  call three_nucleon(an_st_2,an_st_3,an_st_1,1)
               end do
               do ia_3=1,3
                  cr_st_1=occi(ia_3)
                  if (m2_spi(cr_st_1)/=m_beta) cycle
                  if (mt2_spi(cr_st_1)/=mt_beta) cycle
                  ia_1=mod(ia_3+1,3)
                  if (ia_1==0) ia_1=ia_3+1
                  ia_2=mod(ia_3+2,3)
                  if (ia_2==0) ia_2=ia_3+2
                  an_st_1=occi(ia_1)
                  an_st_2=occi(ia_2)
                  phase=+1     ! (-1)**(A-i1+A-i2+1) for i1<i2  old: (-1)**mod(ia_3+1,2)
                  beta=iobsind(cr_st_1)
                  j_beta=j2_spi(cr_st_1)
                  kff=beta
                  call two_nucleon(an_st_1,an_st_2)
               end do
            elseif (nucleonsi==4) then
               do ia_4=1,4
                  cr_st_1=occi(ia_4)
                  if (m2_spi(cr_st_1)/=m_beta) cycle
                  if (mt2_spi(cr_st_1)/=mt_beta) cycle
                  ia_1=mod(ia_4+1,4)
                  if (ia_1==0) ia_1=ia_4+1
                  ia_2=mod(ia_4+2,4)
                  if (ia_2==0) ia_2=ia_4+2
                  ia_3=mod(ia_4+3,4)
                  if (ia_3==0) ia_3=ia_4+3
                  an_st_1=occi(ia_1)
                  an_st_2=occi(ia_2)
                  an_st_3=occi(ia_3)
                  phase=(-1)**(ia_1+ia_2+ia_3) ! (-1)**(A-i1+A-i2+A-i3) old: (-1)**mod(ia_4+1,2)
                  beta=iobsind(cr_st_1)
                  j_beta=j2_spi(cr_st_1)
                  delta_mj=m2_spi(an_st_1)+m2_spi(an_st_2)
     +                 +m2_spi(an_st_3)
                  delta_mt=mt2_spi(an_st_1)+mt2_spi(an_st_2)
     +                 +mt2_spi(an_st_3)
                  kff=beta
                  do kii=ki,ki+nki-1
                     do iT123=max(abs(delta_mt),abs(it2i(kii)-1)),
     $                    min(it2i(kii)+1,3),2
                        do J123=max(abs(delta_mj),
     $                       abs(jt2i(kii)-j_beta)),jt2i(kii)+j_beta,2
                           
                           call three_nucleon(an_st_1,an_st_2,an_st_3,1)
                           call three_nucleon(an_st_1,an_st_3, an_st_2,
     $                          -1)
                           call three_nucleon(an_st_2,an_st_3,an_st_1,1)

                        end do
                     end do
                  end do
               end do
               an_st_1=occi(1)
               an_st_2=occi(2)
               an_st_3=occi(3)
               an_st_4=occi(4)
               delta_mj=mjtotali
               delta_mt=mttotali
               phase=+1    ! (-1)**(A-i1+A-i2+A-i3+A-i4) 
               do kii=ki,ki+nki-1
                  J1234=jt2i(kii)
                  iT1234=it2i(kii)
                  call four_nucleon(an_st_1,an_st_2,an_st_3,an_st_4,1)
                  call four_nucleon(an_st_1,an_st_3,an_st_2,an_st_4,-1)
                  call four_nucleon(an_st_2,an_st_3,an_st_1,an_st_4,1)
                  call four_nucleon(an_st_1,an_st_2,an_st_4,an_st_3,-1)
                  call four_nucleon(an_st_1,an_st_4,an_st_2,an_st_3,1)
                  call four_nucleon(an_st_2,an_st_4,an_st_1,an_st_3,-1)
                  call four_nucleon(an_st_1,an_st_4,an_st_3,an_st_2,-1)
                  call four_nucleon(an_st_1,an_st_3,an_st_4,an_st_2,1)
                  call four_nucleon(an_st_3,an_st_4,an_st_1,an_st_2,1)
                  call four_nucleon(an_st_2,an_st_4,an_st_3,an_st_1,1)
                  call four_nucleon(an_st_2,an_st_3,an_st_4,an_st_1,-1)
                  call four_nucleon(an_st_3,an_st_4,an_st_2,an_st_1,-1)
               end do
            endif
         endif

      end do                   !i1
      deallocate(occi,occim1,occf)
      deallocate(intbas,locz)
      if (twobdcal) then
         if (nucleonsi>2) deallocate(occim2)
      endif
      call cleb_destroy(1)
      call num_of_3j_types_destroy
      print *,' exit main loop: iproc=',iproc

       call MPI_Barrier(icomm,ierr)
c      print *,' main loop send-receive finished: iproc=',iproc
      if (nproc>1) then
         allocate(tdJtemp(ist1dim,ist1dim,
     +        jtotal2min:jtotal2,kf:kf+nkf-1,ki:ki+nki-1))        
         ii=nkf*nki*(jtotal2-jtotal2min+1)*(ist1dim**2)
         call MPI_Reduce(tdJp(1,1,jtotal2min,kf,ki),tdJtemp,ii,
     +        MPI_REAL8,MPI_SUM,0,icomm,ierr)
         tdJp=tdJtemp
         call MPI_Reduce(tdJn(1,1,jtotal2min,kf,ki),tdJtemp,ii,
     +        MPI_REAL8,MPI_SUM,0,icomm,ierr)
         tdJn=tdJtemp
         deallocate(tdJtemp)
         if (twobdcal) then
c            allocate(t2bdtemp(ist2dim,ist2dim,jtotal2min:jtotal2,
c     +           -1:1,kf:kf+nkf-1,ki:ki+nki-1))          
c            ii=3*nkf*nki*(jtotal2-jtotal2min+1)*(ist2dim**2)
c            call MPI_Reduce(t2bd(1,1,jtotal2min,-1,kf,ki),t2bdtemp,ii,
c     +           MPI_REAL8,MPI_SUM,0,icomm,ierr)
c            t2bd=t2bdtemp
            call MPI_Barrier(icomm,ierr)
c            if (allocated(t2bdtemp)) deallocate(t2bdtemp)
            allocate(t2bdtemp(ist2dim))
            do kii=ki,ki+nki-1
               do kff=kf,kf+nkf-1
                  Jtr_mi=tbd(kff,kii)%Jtr_min
                  Jtr_ma=tbd(kff,kii)%Jtr_max
                  do jtrans=Jtr_mi,Jtr_ma
                     do ii=1,ist2dim
                        i_block=(ii-1)/ist2dim_block+1
                        if (i_block/=iproc_block) goto 997
                        if (mod(ii,ist2dim_block)==0) then
                           i_local=ist2dim_block
                        else
                           i_local=mod(ii,ist2dim_block)
                        endif
                        mT_mi=tbd(kff,kii)%Jtr(jtrans)%fi(i_local)
     $                       %mT_min
                        mT_ma=tbd(kff,kii)%Jtr(jtrans)%fi(i_local)
     $                       %mT_max
                        if (mT_mi>mT_ma) cycle
                        do mTt=mT_mi,mT_ma 
                           in_mi=tbd(kff,kii)%Jtr(jtrans)%fi(i_local)
     $                          %mT(mTt)%in_min
                           in_ma=tbd(kff,kii)%Jtr(jtrans)%fi(i_local)
     $                          %mT(mTt)%in_max
                           if (nblock>1) then
                           if (iproc/=i_block-1) then
                              call MPI_Send(tbd(kff,kii)%Jtr(jtrans)
     $                             %fi(i_local)%mT(mTt)%in(in_mi:in_ma),
     $                             in_ma-in_mi+1,
     +                             MPI_REAL8,i_block-1,mTt+2,icomm,ierr)
                           else
                              do source=iproc+nblock,nproc-1,nblock
                                 call MPI_Recv(t2bdtemp,
     $                                in_ma-in_mi+1,MPI_REAL8,source,
     $                                mTt+2,icomm,stat,ierr)
                                 do i1=in_mi,in_ma
                                    tbd(kff,kii)%Jtr(jtrans)
     $                                   %fi(i_local)%mT(mTt)%in(i1)=
     $                                   tbd(kff,kii)%Jtr(jtrans)
     $                                   %fi(i_local)%mT(mTt)%in(i1)
     $                                   +t2bdtemp(i1-in_mi+1)
                                 end do
                              end do
                           endif
                           else
                           call MPI_Reduce(tbd(kff,kii)%Jtr(jtrans)
     $                          %fi(ii)%mT(mTt)%in(in_mi:in_ma),
     $                          t2bdtemp,in_ma-in_mi+1,
     +                          MPI_REAL8,MPI_SUM,0,icomm,ierr)
                           if (iproc==0) tbd(kff,kii)%Jtr(jtrans)
     $                          %fi(ii)%mT(mTt)%in(in_mi:in_ma)=
     $                          t2bdtemp(1:in_ma-in_mi+1)
                           endif
                        end do
 997                    continue
                        call MPI_Barrier(icomm,ierr)
                     end do
                  end do
               end do
            end do
            deallocate(t2bdtemp)
         endif
         if (threebdcal.and.nucleonsi==3) then
            allocate(adtemp(ist3dim,1,ki:ki+nki-1))
            ii=nki*ist3dim
            adtemp=0.d0
            call MPI_Reduce(ad_cl(1,1,ki),adtemp,ii,
     +           MPI_REAL8,MPI_SUM,0,icomm,ierr)
            ad_cl=adtemp
            deallocate(adtemp)
            allocate(adtemp(ist2dim,1,1))
            ii=ist2dim
            do kii=ki,ki+nki-1
               do i1=1,ist1dim
                  adtemp=0.d0
c                  print *,' before:iproc,i1,ad_cl_2=',iproc,i1,
c     $                 ad_cl_2(1,i1,kii),ad_cl_2(2,i1,kii)
                  call MPI_Reduce(ad_cl_2(1,i1,kii),adtemp,ii,
     +                 MPI_REAL8,MPI_SUM,0,icomm,ierr)
                  ad_cl_2(:,i1,kii)=adtemp(:,1,1)
c                  print *,' after:iproc,i1,ad_cl_2=',iproc,i1,
c     $                 ad_cl_2(1,i1,kii),ad_cl_2(2,i1,kii)
               end do
            end do
            deallocate(adtemp)
         endif
         if (threebdcal.and.nucleonsi==4) then
            allocate(adtemp(ist3dim,1,1))
            ii=ist3dim
            do kii=ki,ki+nki-1
               do i1=1,ist1dim
                  call MPI_Reduce(ad_cl(1,i1,kii),adtemp,ii,
     +                 MPI_REAL8,MPI_SUM,0,icomm,ierr)
                  ad_cl(:,i1,kii)=adtemp(:,1,1)
               end do
            end do
            deallocate(adtemp)
            allocate(adtemp(ist4dim,1,1))
            ii=ist4dim
            do kii=ki,ki+nki-1
               call MPI_Reduce(ad_cl_2(1,1,kii),adtemp,ii,
     +              MPI_REAL8,MPI_SUM,0,icomm,ierr)
               ad_cl_2(:,1,kii)=adtemp(:,1,1)
            end do
            deallocate(adtemp)
         endif
      endif

c      print *,' before return: iproc=',iproc
      if (iproc>nblock-1) return
c      print *,' after return: iproc=',iproc

      if (iproc==0) then
      do kii=ki,ki+nki-1
         jtotali2=jt2i(kii)
         itotali2=it2i(kii)
         if (dabs(dble(jtotali2)-2.d0*Jxi(kii))>1.d-4) then
            if (iproc==0) print *, '**** not a good J ****'
         endif 

         do kff=kf,kf+nkf-1
            jtotalf2=jt2f(kff)
            itotalf2=it2f(kff)
            if (dabs(dble(jtotalf2)-2.d0*Jxf(kff))>1.d-4) then
               if (iproc==0) print *, '**** not a good J ****'
            endif 

            if (iproc==0) write(iunitvout,1200) kff,jtotalf2,itotalf2,
     +           enerf(kff)-eneri(ki),
     +           kii,jtotali2,itotali2,eneri(kii)-eneri(ki)
 1200       format(//,' *** Transition matrix elements for states:',
     +           ' ***',/,' #',i3,
     +           ' [2*(J,T),Ex]_f=',i3,i2,f8.4,3x,
     +           '#',i3,' [2*(J,T),Ex]_i=',i3,i2,f8.4) 
c*1200       format(//,' *** Transition matrix elements for states:',
c*     +           ' ***',/,i3,
c*     +           '   2*J_f=',i3,'   Ex_f=',f9.4,5x,
c*     +           i3,'   2*J_i=',i3,'   Ex_i=',f9.4)

            do jtrans=jtotal2min,jtotal2

               if (iabs(jt2i(kii)-jt2f(kff))>2*jtrans
     +              .or.(jt2i(kii)+jt2f(kff))<2*jtrans) cycle

               cleb=clebd(jtotali2,mjtotali,2*jtrans,mjtotalf-mjtotali,
     +              jtotalf2,mjtotalf)

               if (dabs(cleb)<1.d-8) cycle

               if (iproc==0) write(iunitvout,1953) jtrans
 1953          format(/,' Jtrans=',i3)

               cleb=dsqrt(dble(jtotalf2)+1.d0)/cleb 
c*** Glendenning definition
               clb1=-cleb/dsqrt(dble(2*jtrans)+1.d0) 
c*** Alternative definition
c*               clb1=-1.d0               

               do i1=1,ist1dim
                  n1=ist1(1,i1)
                  l1=ist1(2,i1)
                  j1=ist1(3,i1)
                  do i2=1,ist1dim
                     n2=ist1(1,i2)
                     l2=ist1(2,i2)
                     j2=ist1(3,i2)
                     if((j1+j2)<2*jtrans) cycle
                     if(iabs(j1-j2)>2*jtrans) cycle
                     if (((-1)**l1==(-1)**l2.and.
     +                    iparityi==iparityf).or.
     +           ((-1)**l1/=(-1)**l2.and.iparityi/=iparityf)) then
                        if (iproc==0) write(iunitvout,1962) i1,i2,
     +                       clb1*tdJp(i1,i2,jtrans,kff,kii),
     +                       clb1*tdJn(i1,i2,jtrans,kff,kii)
 1962                   format(i4,1x,i4,1x,e15.8,1x,e15.8)
c 1960                   format(' a+=',i3,'    a-=',i3,
c     +                      '     td(a+,a-): p=',f10.6,'     n=',f10.6)
                     endif
                  end do
               end do    
            end do
         end do
      end do
      endif
      if (twobdcal) then
         if (iproc==0) then
            write(iunitvout,*)
     +'***************************************************************'
            write(iunitvout,*)
            write(iunitvout,*) ' Two-body transition matrix elements'
            write(iunitvout,*) ' ***********************************'
            write(iunitvout,7778) ist2dim
 7778       format(' number of two-nucleon states =',i5)
            if (isum2<10000000) then
               write(iunitvout,7878) isum2
            else
               kii=-1
               write(iunitvout,7878) kii
            endif
 7878       format(' number of two-body Hamiltonian matrix elements =',
     +           i7)
            do ii=1,ist2dim
               write(iunitvout,7780) ii,(ist2(ia,ii),ia=1,4)
 7780          format(' #',i5,'  a b J T=',3i3,i2)
            end do
         endif   

         allocate(t2bdtemp(ist2dim))

         do kii=ki,ki+nki-1

            jtotali2=jt2i(kii)
            itotali2=it2i(kii)
            if (dabs(dble(jtotali2)-2.d0*Jxi(kii))>1.d-4) then
               if (iproc==0) print *,'**** not a good J ****'
            endif 

            do kff=kf,kf+nkf-1

               jtotalf2=jt2f(kff)
               itotalf2=it2f(kff)
               if (dabs(dble(jtotalf2)-2.d0*Jxf(kff))>1.d-4) then
                  if (iproc==0)
     $                 print *, '**** not a good J ****'
               endif 

               if (iproc==0) write(iunitvout,1200) kff,
     $              jtotalf2,
     $              itotalf2,enerf(kff)-eneri(ki),
     +              kii,jtotali2,itotali2,eneri(kii)-eneri(ki)


               Jtr_mi=tbd(kff,kii)%Jtr_min
               Jtr_ma=tbd(kff,kii)%Jtr_max

               do jtrans=Jtr_mi,Jtr_ma

                  cleb=clebd(jtotali2,mjtotali,2*jtrans,
     +                 mjtotalf-mjtotali,jtotalf2,mjtotalf)

                  if (dabs(cleb)<1.d-8) cycle

                  if (iproc==0) write(iunitvout,1953) jtrans

                  clb1=dsqrt(dble(jtotalf2+1))/cleb 
c*** Alternative definition
c*               clb1=-1.d0               

                  do ibl=1,nblock

                     if (iproc/=ibl-1.and.iproc/=0) cycle
cc            call MPI_Barrier(icomm,ierr)

                     if (ibl>1) then
c               write(7+iproc,*) ' ibl,iproc=',ibl,iproc

                        do ii=1,ist2dim
                           i_block=(ii-1)/ist2dim_block+1
                           if (i_block/=ibl) cycle

                           if (mod(ii,ist2dim_block)==0) then
                              i_local=ist2dim_block
                           else
                              i_local=mod(ii,ist2dim_block)
                           endif

                           if (iproc==i_block-1) then
                              mT_mi=tbd(kff,kii)%Jtr(jtrans)%fi(i_local)
     $                             %mT_min
                              mT_ma=tbd(kff,kii)%Jtr(jtrans)%fi(i_local)
     $                             %mT_max
                              if (mT_mi>mT_ma) cycle
                           elseif(iproc==0) then
                              i1=ist2(1,ii)
                              i2=ist2(2,ii)
                              j12=ist2(3,ii)
                              it12=ist2(4,ii)
                              mT_mi=max(-abs(it12),(mttotalf-mttotali)
     $                             /2-1)
                              mT_ma=min(abs(it12),(mttotalf-mttotali)
     $                             /2+1)
                              tbd(kff,kii)%Jtr(jtrans)%fi(i_local)
     $                             %mT_min=mT_mi
                              tbd(kff,kii)%Jtr(jtrans)%fi(i_local)
     $                             %mT_max=mT_ma
                              if (mT_mi>mT_ma) cycle
                              l1=ist1(2,i1)
                              l2=ist1(2,i2)
                              pi=(-1)**(l1+l2+iparityi+iparityf)
                              if (pi/=(-1)**nhom12max) then
                                 J12ma=nhom12max
                              else
                                 J12ma=nhom12max+1
                              endif
                              if (abs(jtrans-j12)>ist2_Jst_J12max
     $                             .or.abs(jtrans-j12)<ist2_Jst_J12min
     $                             .or.min(jtrans+j12,J12ma)+1
     $                             >ist2_Jst_J12max+1
     $                             .or.min(jtrans+j12,J12ma)+1
     $                             <ist2_Jst_J12min)
     $                             then
                                 in_mi=-1
                                 in_ma=-1
                              else
                                 in_mi=ist2_Jst((pi+1)/2,
     $                                abs(jtrans-j12))
                                 in_ma=ist2_Jst((pi+1)/2,
     $                                min(jtrans+j12,J12ma)+1)-1
                              endif
c                              in_mi=ist2_Jst((pi+1)/2,abs(jtrans-j12))
c                              in_ma=ist2_Jst((pi+1)/2,min(jtrans+j12,
c     $                             J12ma)+1)-1
                              deallocate(tbd(kff,kii)%Jtr(jtrans)
     $                             %fi(i_local)
     $                             %mT)
                              allocate(tbd(kff,kii)%Jtr(jtrans)
     $                             %fi(i_local)
     $                             %mT(mT_mi:mT_ma))
                              do mTt=mT_mi,mT_ma 
                                 tbd(kff,kii)%Jtr(jtrans)%fi(i_local)
     $                                %mT(mTt)%in_min=in_mi
                                 tbd(kff,kii)%Jtr(jtrans)%fi(i_local)
     $                                %mT(mTt)%in_max=in_ma
                                 allocate(tbd(kff,kii)%Jtr(jtrans)
     $                                %fi(i_local)
     $                                %mT(mTt)%in(in_mi:in_ma))
                                 tbd(kff,kii)%Jtr(jtrans)%fi(i_local)
     $                                %mT(mTt)%in=0.d0
                              end do
c                           else
c                              i1=ist2(1,ii)
c                              i2=ist2(2,ii)
c                              j12=ist2(3,ii)
c                              it12=ist2(4,ii)
c                              mT_mi=max(-abs(it12),(mttotalf-mttotali)
c     $                             /2-1)
c                              mT_ma=min(abs(it12),(mttotalf-mttotali)
c     $                             /2+1)
c                              if (mT_mi>mT_ma) cycle
                           endif
                           do mTt=mT_mi,mT_ma 

                              if (iproc==i_block-1) then
                                 in_mi=tbd(kff,kii)%Jtr(jtrans)
     $                                %fi(i_local)%mT(mTt)%in_min
                                 in_ma=tbd(kff,kii)%Jtr(jtrans)
     $                                %fi(i_local)%mT(mTt)%in_max
                                 
c                                 write(7+iproc,*) 
c     $                      ' SSend:iproc,i_local,ii,mTt,in_mi,in_ma=',
c     $                                iproc,i_local,ii,mTt,in_mi,in_ma

                                 call MPI_SSend(tbd(kff,kii)%Jtr(jtrans)
     $                                %fi(i_local)%mT(mTt)
     $                                %in(in_mi:in_ma),in_ma-in_mi+1,
     +                                MPI_REAL8,0,i_block,icomm,
     $                                ierr)

c                                 write(7+iproc,*) 
c     $                      ' SSent:iproc,i_local,ii,mTt,in_mi,in_ma=',
c     $                                iproc,i_local,ii,mTt,in_mi,in_ma

                              elseif (iproc==0) then
                                 in_mi=tbd(kff,kii)%Jtr(jtrans)
     $                                %fi(i_local)%mT(mTt)%in_min
                                 in_ma=tbd(kff,kii)%Jtr(jtrans)
     $                                %fi(i_local)%mT(mTt)%in_max

c                                 write(7+iproc,*) 
c     $                      ' Recv:iproc,i_local,ii,mTt,in_mi,in_ma=',
c     $                                iproc,i_local,ii,mTt,in_mi,in_ma

                                 call MPI_Recv(t2bdtemp,
     $                                in_ma-in_mi+1,MPI_REAL8,i_block-1,
     $                                i_block,
     $                                icomm,stat,ierr)

c                                 write(7+iproc,*) 
c     $                    ' Recved:iproc,i_local,ii,mTt,in_mi,in_ma=',
c     $                                iproc,i_local,ii,mTt,in_mi,in_ma

                                 do i1=in_mi,in_ma
                                    tbd(kff,kii)%Jtr(jtrans)
     $                                   %fi(i_local)%mT(mTt)%in(i1)=
     $                                   t2bdtemp(i1-in_mi+1)
                                 end do
                              endif
cc                              call MPI_Barrier(icomm,ierr)
                           end do
                        end do
                     endif

                     if (iproc>0) cycle

cc         do kii=ki,ki+nki-1


cc            do kff=kf,kf+nkf-1


c               do jtrans=jtotal2min,jtotal2

c                  jtrans2=jtrans+jtrans
c                  if (iabs(jt2i(kii)-jt2f(kff))>jtrans2
c     +                 .or.(jt2i(kii)+jt2f(kff))<jtrans2) cycle
               
                     do i1=1,ist2dim
                     
                        i_block=(i1-1)/ist2dim_block+1
                        if (ibl/=i_block) cycle
                        if (mod(i1,ist2dim_block)==0) then
                           i_local=ist2dim_block
                        else
                           i_local=mod(i1,ist2dim_block)
                        endif

                        nonebody=ist2(1,i1)
                        nonebody2=ist2(2,i1)
                        j12=ist2(3,i1)
                        it12=ist2(4,i1)

                        l1=ist1(2,nonebody)
                        l2=ist1(2,nonebody2)

                        do i3=1,ist2dim

                           nonebody3=ist2(1,i3)
                           nonebody4=ist2(2,i3)
                           j34=ist2(3,i3)
                           it34=ist2(4,i3)

                           if((j12+j34)<jtrans) cycle
                           if(iabs(j12-j34)>jtrans) cycle

                           l3=ist1(2,nonebody3)
                           l4=ist1(2,nonebody4)

                           if (((-1)**(l1+l2+l3+l4)==1.and.
     +                          iparityi==iparityf).or.
     +                          ((-1)**(l1+l2+l3+l4)==-1
     +                          .and.iparityi/=iparityf)) then
                              tbx=0.d0
                              if (mttotali==mttotalf) then
                                 tbx(1)=tbd(kff,kii)%Jtr(jtrans)
     $                                %fi(i_local)%mT(0)%in(i3)*clb1
                                 if (it12==1) then
                                    tbx(2)=tbd(kff,kii)%Jtr(jtrans)
     $                                   %fi(i_local)%mT(1)%in(i3)*clb1
                                    tbx(3)=tbd(kff,kii)%Jtr(jtrans)
     $                                   %fi(i_local)%mT(-1)%in(i3)*clb1
                                 endif
                              else
                                 mT_mi=tbd(kff,kii)%Jtr(jtrans)
     $                                %fi(i_local)
     $                                %mT_min
                                 mT_ma=tbd(kff,kii)%Jtr(jtrans)
     $                                %fi(i_local)
     $                                %mT_max
                                 do mTt=mT_mi,mT_ma
                                    tbx(mTt-mT_mi+1)=tbd(kff,kii)
     $                                   %Jtr(jtrans)
     $                                   %fi(i_local)%mT(mTt)%in(i3)
     $                                   *clb1
                                 end do
                              endif

                              if (abs(tbx(1))
     $                             >4.d-9.or.
     $                             abs(tbx(2))
     $                             >4.d-9.or.
     $                             abs(tbx(3))
     $                             >4.d-9) then
                                 write(iunitvout,1968) i1,i3,tbx(1:3)
                              endif
 1968                         format('tb',i6,1x,i6,1x,e15.8,1x,
     $                             e15.8,1x,e15.8)
c 1967                         format(' (a+a+)J=',i5,'  (a-a-)J=',i5,
c     +                             '   td: pn=',f10.6,'   pp=',f10.6,
c     +                             '   nn=',f10.6)
                           endif
                        end do
                     end do    
                  end do
               end do
            end do
         end do
         deallocate(t2bdtemp)
cc         deallocate(tbd_upda)
cc         deallocate(buffer_count)
cc     deallocate(buffer_bsend)
      endif   
      if (threebdcal.and.(nucleonsi==3.or.nucleonsi==4)) then
         if (iproc==0) then
            write(44) nhomi,nhom12i,nhom123i
c            print *,nhomi,nhom12i,nhom123i
            write(44) ist3dim
c            print *,' ist3dim=',ist3dim
            do ii=1,ist3dim
               write(44) ii,(ist3(ia,ii),ia=1,4)
c               print *,ii,(ist3(ia,ii),ia=1,4)
            end do
            do kii=ki,ki+nki-1
               jtotali2=jt2i(kii)
               itotali2=it2i(kii)
               if (nucleonsi==3) then
                  cleb=2.d0*sqrt(real((jtotali2+1)*(itotali2+1),
     $                 kind(0.d0)))
               endif
               write(44) kii,jtotali2,itotali2,eneri(kii)
c               print *,kii,jtotali2,itotali2,eneri(kii)
               ii=0
               if (nucleonsi==3) then
                  do i1=1,ist3dim
                     if (abs(ad_cl(i1,1,kii))>1.d-9) ii=ii+1
                  end do
                  write(44) ii
c                  print *,ii
                  do i1=1,ist3dim
                     if (abs(ad_cl(i1,1,ki))>1.d-9) then
                        write(44) i1,ad_cl(i1,1,kii)*cleb
c                        print *,i1,ad_cl(i1,1,kii)*cleb
                     endif
                  end do
                  ii=0
                  do i2=1,ist1dim
                     j_beta=ist1(3,i2)
                     do i1=1,ist2dim
                        J12=ist2(3,i1)
                        iT12=ist2(4,i1)
                        if (abs(ad_cl_2(i1,i2,kii))>1.d-9) then
                           cleb=clebd(jtotali2,mjtotali,
     $                          2*J12,m_beta-mjtotali,j_beta,m_beta)
     $                          *clebd(itotali2,mttotali,
     $                          2*iT12,mt_beta-mttotali,1,mt_beta)
                           if (abs(cleb)>1.d-9) then
                              ii=ii+1
                           endif
                        endif
                     end do
                  end do
                  write(44) ii
c                  print *,ii
                  do i2=1,ist1dim
                     j_beta=ist1(3,i2)
                     do i1=1,ist2dim
                        J12=ist2(3,i1)
                        iT12=ist2(4,i1)
                        if (abs(ad_cl_2(i1,i2,kii))>1.d-9) then
                           cleb=clebd(jtotali2,mjtotali,
     $                          2*J12,m_beta-mjtotali,j_beta,m_beta)
     $                          *clebd(itotali2,mttotali,
     $                          2*iT12,mt_beta-mttotali,1,mt_beta)
                           if (abs(cleb)>1.d-9) then
                              cleb=2.d0*sqrt(real((j_beta+1)
     $                             *2,kind(0.d0)))/cleb
                              write(44) i1,i2,ad_cl_2(i1,i2,kii)*cleb
c                              print *,i1,i2,ad_cl_2(i1,i2,kii)*cleb
                           endif
                        endif
                     end do
                  end do
               elseif (nucleonsi==4) then
                  do i2=1,ist1dim
                     j_beta=ist1(3,i2)
                     do i1=1,ist3dim
                        J123=ist3(3,i1)
                        iT123=ist3(4,i1)
                        if (abs(ad_cl(i1,i2,kii))>1.d-9) then
                           cleb=clebd(j_beta,m_beta,J123,
     $                          mjtotali-m_beta,jtotali2,mjtotali)
     $                          *clebd(1,mt_beta,iT123,
     $                          mttotali-mt_beta,itotali2,mttotali)
                           if (abs(cleb)>1.d-9) then
                              ii=ii+1
                           endif
                        endif
                     end do
                  end do
                  write(44) ii
                  do i2=1,ist1dim
                     j_beta=ist1(3,i2)
                     do i1=1,ist3dim
                        J123=ist3(3,i1)
                        iT123=ist3(4,i1)
                        if (abs(ad_cl(i1,i2,kii))>1.d-9) then
                           cleb=clebd(j_beta,m_beta,J123,
     $                          mjtotali-m_beta,jtotali2,mjtotali)
     $                          *clebd(1,mt_beta,iT123,
     $                          mttotali-mt_beta,itotali2,mttotali)
                           if (abs(cleb)>1.d-9) then
                              cleb=2.d0*sqrt(real((jtotali2+1)
     $                             *(itotali2+1),kind(0.d0)))/cleb
                              write(44) i1,i2,ad_cl(i1,i2,kii)*cleb
                           endif
                        endif
                     end do
                  end do
               endif
            end do
            if (nucleonsi==4) then
               write(44) nhom1234i
c            print *,nhom1234i
               write(44) ist4dim
c            print *,' ist4dim=',ist4dim
               do ii=1,ist4dim
                  write(44) ii,(ist4(ia,ii),ia=1,4)
c               print *,ii,(ist3(ia,ii),ia=1,4)
               end do
               do kii=ki,ki+nki-1
                  jtotali2=jt2i(kii)
                  itotali2=it2i(kii)
                  cleb=2.d0*sqrt(real((jtotali2+1)*(itotali2+1),
     $                 kind(0.d0)))
                  write(44) kii,jtotali2,itotali2,eneri(kii)
c               print *,kii,jtotali2,itotali2,eneri(kii)
                  ii=0
                  do i1=1,ist4dim
                     if (abs(ad_cl_2(i1,1,kii))>1.d-9) ii=ii+1
                  end do
                  write(44) ii
c                  print *,ii
                  do i1=1,ist4dim
                     if (abs(ad_cl_2(i1,1,ki))>1.d-9) then
                        write(44) i1,ad_cl_2(i1,1,kii)*cleb
c                        print *,i1,ad_cl_2(i1,1,kii)*cleb
                     endif
                  end do
               end do
            endif
         endif
      endif
      contains

      subroutine one_body(c_1,a_1,obd,j_min,j_max,k1,n1,k2,n2)
      implicit none
      integer,intent(IN) :: c_1,a_1,j_min,j_max,k1,n1,k2,n2
      real(kind(0.d0)),intent(INOUT) ::
     +     obd(j_min:j_max,k2:k2+n2-1,k1:k1+n1-1)
      integer :: kii,kff,jtrans
      real(kind(0.d0)) :: cleb
cc      if (twobdcal) call check_incoming
      do jtrans=max(jtotal2min,
     +     iabs(j2_spf(c_1)-j2_spi(a_1))/2),
     +     min(jtotal2,(j2_spf(c_1)+j2_spi(a_1))/2)
         cleb=cleb_3j(1,j2_spf(c_1),m2_spf(c_1),
     +        j2_spi(a_1),-m2_spi(a_1),2*jtrans)
         do kii=ki,ki+nki-1
            do kff=kf,kf+nkf-1
               if (iabs(jt2i(kii)-jt2f(kff))>2*jtrans
     +              .or.(jt2i(kii)+jt2f(kff))<2*jtrans) cycle
               obd(jtrans,kff,kii)=obd(jtrans,kff,kii)+
     +              bmpi(i1,kii)*bmpf(i2,kff)*cleb
     +              *real(phase,kind(0.d0))                           
            end do
         end do 
      end do
      end subroutine one_body

      subroutine two_body(icrstatep,icrstate,ianstatep,ianstate)
      implicit none
      integer,intent(IN) :: icrstatep,icrstate,ianstatep,ianstate
      integer :: nonebody,nonebody2,nonebody3,nonebody4,j12,it12,
     $     j34,it34,ntwobody,ipht12,ntwobody2,ipht34,jtrans,kii,kff
      real(kind(0.d0)) :: fact,fact2,cleb,clebt,clebt2,cleb2,clebjt
      integer :: ntwobody_block,ntwobody_local,count
      real(kind(0.d0)) :: tbd_upd
      integer :: kf_in,ki_in,Jtr,i_loc,mT,in
cc      call check_incoming
      nonebody=iobsind(icrstate)
      nonebody3=iobsind(ianstate)
      nonebody2=iobsind(icrstatep)
      nonebody4=iobsind(ianstatep)

      if (nonebody==nonebody2) then 
         fact=dsqrt(2.d0)
      else
         fact=1.d0
      endif
      if (nonebody3==nonebody4) then 
         fact=fact*dsqrt(2.d0)
      else
         fact=1.d0*fact
      endif

      do kii=ki,ki+nki-1
         do kff=kf,kf+nkf-1

            do jtrans=max(abs(jt2i(kii)-jt2f(kff))/2,jtotal2min),
     +           min((jt2i(kii)+jt2f(kff))/2,jtotal2)

               do j12=max(iabs(m2_spf(icrstate)
     +              +m2_spf(icrstatep))/2,
     +              iabs(j2_spf(icrstate)-j2_spf(icrstatep))/2),
     +              (j2_spf(icrstate)+j2_spf(icrstatep))/2

                  fact2=fact/dsqrt(dble(2*j12+1)) 
                  cleb=clebjj(j2_spf(icrstate)/2,m2_spf(icrstate),
     +                 j2_spf(icrstatep)/2,m2_spf(icrstatep),j12)
c     *                  cleb=clebd(j2_spf(icrstate),m2_spf(icrstate),
c     *     +                 j2_spf(icrstatep),m2_spf(icrstatep),2*j12,
c     *     +                 m2_spf(icrstate)+m2_spf(icrstatep))
                     
                  do it12=iabs(mt2_spf(icrstate)
     +                 +mt2_spf(icrstatep))/2,1
                     
                     if (nonebody==nonebody2) then
                        if ((-1)**(j12+it12)/=-1) cycle
                     end if

                     clebt=clebjj(0,mt2_spf(icrstate),
     +                    0,mt2_spf(icrstatep),it12)
c*                     clebt=clebd(1,mt2_spf(icrstate),
c*     +                    1,mt2_spf(icrstatep),2*it12,
c*     +                    mt2_spf(icrstate)+mt2_spf(icrstatep))
                     if(nonebody<=nonebody2) then
                        ntwobody=itbind(nonebody+
     +                       nonebody2*(nonebody2-1)/2,j12,it12)
                        ipht12=1
                     else
                        ntwobody=itbind(nonebody2+
     +                       nonebody*(nonebody-1)/2,j12,it12)
                        ipht12=(-1)**(j12-(j2_spf(icrstate)
     +                       +j2_spf(icrstatep))/2+it12)
                     endif   

                     ntwobody_block=(ntwobody-1)/ist2dim_block+1
                     if (mod(ntwobody,ist2dim_block)==0) then
                        ntwobody_local=ist2dim_block
                     else
                        ntwobody_local=mod(ntwobody,ist2dim_block)
                     endif

c                     print *,' iproc,ntwobody_local=',iproc,
c     $                    ntwobody_local

                     if (iproc_block==ntwobody_block) then

                        do j34=max(iabs(j2_spi(ianstate)
     +                       -j2_spi(ianstatep))/2,
     +                       iabs(m2_spi(ianstate)
     +                       +m2_spi(ianstatep))/2,abs(j12-jtrans)),
     $                       min((j2_spi(ianstate)
     +                       +j2_spi(ianstatep))/2,jtotal2+j12,
     $                       j12+jtrans) 
                           
                           cleb2=clebjj(j2_spi(ianstate)/2,
     +                          m2_spi(ianstate),j2_spi(ianstatep)/2,
     +                          m2_spi(ianstatep),j34)
c*                        cleb2=clebd(j2_spi(ianstate),
c*     +                       m2_spi(ianstate),j2_spi(ianstatep),
c*     +                       m2_spi(ianstatep),2*j34,
c*     +                       m2_spi(ianstate)+m2_spi(ianstatep))

                           do it34=iabs(mt2_spi(ianstate)
     +                          +mt2_spi(ianstatep))/2,1
                                 
                              if (nonebody3==nonebody4) then
                                 if ((-1)**(j34+it34)/=-1) cycle
                              end if
                           
                              clebt2=clebjj(0,mt2_spi(ianstate),
     +                             0,mt2_spi(ianstatep),it34)
c*                           clebt2=clebd(1,mt2_spi(ianstate),
c*     +                          1,mt2_spi(ianstatep),2*it34,
c*     +                          mt2_spi(ianstate)+mt2_spi(ianstatep))
                              
                              if(nonebody3<=nonebody4) then
                                 ntwobody2=itbind(nonebody3+
     +                                nonebody4*(nonebody4-1)/2,
     $                                j34,it34)
                                 ipht34=1 
                              else
                                 ipht34=(-1)**(j34-(j2_spi(ianstate)
     +                                +j2_spi(ianstatep))/2+it34)
                                 ntwobody2=itbind(nonebody4+
     +                                nonebody3*(nonebody3-1)/2,
     $                                j34,it34)
                              endif   
                              
                                 
                              clebjt=clebj12j34(j34,(m2_spi(ianstate)
     +                             +m2_spi(ianstatep))/2,
     +                             jtrans,
     +                             j12,(m2_spi(icrstate)
     +                             +m2_spi(icrstatep))/2)
c*                              clebjt=clebd(2*j34,m2_spi(ianstate)
c*     +                             +m2_spi(ianstatep),
c*     +                             2*jtrans,mjtotalf-mjtotali,
c*     +                             2*j12,m2_spi(icrstate)
c*     +                             +m2_spi(icrstatep))
                                 
                              tbd_upd=clebjt*cleb*cleb2
     +                             *clebt*clebt2*bmpi(i1,kii)
     $                             *bmpf(i2,kff)*fact2
     +                             *dble(iphase*ipht12*ipht34)

                              mT=(mt2_spi(icrstate)
     +                             +mt2_spi(icrstatep))/2

                              tbd(kff,kii)%Jtr(jtrans)
     $                             %fi(ntwobody_local)
     $                             %mT(mT)%in(ntwobody2)
     $                             =tbd(kff,kii)%Jtr(jtrans)
     $                             %fi(ntwobody_local)
     $                             %mT(mT)%in(ntwobody2)
     $                             +tbd_upd

c                              print *,
c     $             ' direct upd: iproc,ntwobody_local,ntwobody_block=',
c     $                             iproc,ntwobody_local,ntwobody_block

                           end do
                        end do   
            
                     endif
                  end do
               end do
            end do
         end do
      end do
      end subroutine two_body
      subroutine three_nucleon(a_1,a_2,a_3,perm)
      implicit none
      integer,intent(IN) :: a_1,a_2,a_3,perm
      integer :: sp_st_1,sp_st_2,sp_st_12,tb_st
      integer :: iT12,J12,sp_st_3,thr_st
      real(kind=kind(0.d0)) :: pht12
      real(kind=kind(0.d0)) :: cleb1,cleb2,clebt1,clebt2

      sp_st_1=iobsind(a_1)
      sp_st_2=iobsind(a_2)
      sp_st_3=iobsind(a_3)
      if (sp_st_1<=sp_st_2) then
         sp_st_12=sp_st_1+sp_st_2*(sp_st_2-1)/2
      else
         sp_st_12=sp_st_2+sp_st_1*(sp_st_1-1)/2
      endif

      do J12=max(abs(j2_spi(a_1)-j2_spi(a_2)),
     +     abs(m2_spi(a_1)+m2_spi(a_2)),
     +     abs(J123-j2_spi(a_3))),
     +     min((j2_spi(a_1)+j2_spi(a_2)),
     +     (J123+j2_spi(a_3))),2
c         cleb1=cleb_3j(1,j2_spi(a_1),m2_spi(a_1),
c     +        j2_spi(a_2),m2_spi(a_2),J12)
         cleb1=clebd(j2_spi(a_1),m2_spi(a_1),
     +        j2_spi(a_2),m2_spi(a_2),J12,
     +        m2_spi(a_1)+m2_spi(a_2))
         if (cleb1==0.d0) cycle
c         cleb2=cleb_3j(2,J12,m2_spi(a_1)+m2_spi(a_2),
c     +        j2_spi(a_3),m2_spi(a_3),J123)
         cleb2=clebd(J12,m2_spi(a_1)+m2_spi(a_2),
     +        j2_spi(a_3),m2_spi(a_3),J123,delta_mj)
         if (cleb2==0.d0) cycle
         
         do iT12=max(abs(mt2_spi(a_1)+mt2_spi(a_2)),
     +        abs(iT123-1)),2,2
            if (sp_st_1==sp_st_2.and.
     +           (-1)**((J12+iT12)/2)/=-1) cycle
            
c            clebt1=cleb_3j(3,1,mt2_spi(a_1),1,
c     +           mt2_spi(a_2),iT12)
            clebt1=clebd(1,mt2_spi(a_1),1,
     +           mt2_spi(a_2),iT12,mt2_spi(a_1)
     +           +mt2_spi(a_2))
            if (clebt1==0.d0) cycle
c            clebt2=cleb_3j(4,iT12,mt2_spi(a_1)+mt2_spi(a_2),
c     +           1,mt2_spi(a_3),iT123)
            clebt2=clebd(iT12,mt2_spi(a_1)+mt2_spi(a_2),
     +           1,mt2_spi(a_3),iT123,delta_mt)
            if (clebt2==0.d0) cycle
            
            if (sp_st_1<=sp_st_2) then
               pht12=1.d0
            else
               pht12=real((-1)**((J12-j2_spi(a_1)
     +              -j2_spi(a_2)+iT12)/2),kind(0.d0))
            endif
            
            tb_st=itbind(sp_st_12,J12/2,iT12/2)
            thr_st=ad3ind(tb_st,sp_st_3,J123/2,iT123/2)
            
            ad_cl(thr_st,kff,kii)=ad_cl(thr_st,kff,kii)
     +           +cleb1*cleb2*clebt1*clebt2
cc     +           *bmpi(i1,kii)*bmpf(i2,kff)
     +           *bmpi(i1,kii)
     +           *real(phase*perm,kind(0.d0))*pht12
cc     +           *real(perm,kind(0.d0))*pht12

         end do
      end do
      end subroutine three_nucleon
      subroutine two_nucleon(a_1,a_2)
      implicit none
      integer,intent(IN) :: a_1,a_2
      integer :: sp_st_1,sp_st_2,sp_st_12,tb_st,delta_mj,delta_mt
      integer :: iT12,J12
      real(kind=kind(0.d0)) :: pht12

      sp_st_1=iobsind(a_1)
      sp_st_2=iobsind(a_2)
      if (sp_st_1<=sp_st_2) then
         sp_st_12=sp_st_1+sp_st_2*(sp_st_2-1)/2
      else
         sp_st_12=sp_st_2+sp_st_1*(sp_st_1-1)/2
      endif
      
      delta_mj=-m2_spi(a_1)-m2_spi(a_2)
      delta_mt=-mt2_spi(a_1)-mt2_spi(a_2)

      do kii=ki,ki+nki-1
         do iT12=max(abs(delta_mt)/2,abs(it2i(kii)-1)/2),
     $        min(1,(it2i(kii)+1)/2)
            do J12=max(abs(delta_mj)/2,
     +           abs(jt2i(kii)-j_beta)/2,
     +           abs(j2_spi(a_1)-j2_spi(a_2))/2),
     +           min((jt2i(kii)+j_beta)/2,
     +           (j2_spi(a_1)+j2_spi(a_2))/2)
                  
               if (sp_st_1==sp_st_2.and.
     +              (-1)**(J12+iT12)/=-1) cycle
               if (sp_st_1<=sp_st_2) then
                  pht12=1.d0
               else
                  pht12=real((-1)**(J12-(j2_spi(a_1)
     +                 +j2_spi(a_2))/2+iT12),kind(0.d0))
               endif

               tb_st=itbind(sp_st_12,J12,iT12)                  
                  
               ad_cl_2(tb_st,kff,kii)=ad_cl_2(tb_st,kff,kii)
     +              +clebd(j2_spi(a_1),-m2_spi(a_1),j2_spi(a_2),
     +              -m2_spi(a_2),2*J12,delta_mj)
     +              *clebd(1,-mt2_spi(a_1),1,-mt2_spi(a_2),2*iT12,
     +              delta_mt)
     +              *bmpi(i1,kii)
     +              *real(phase*(-1)**((j2_spi(a_1)-m2_spi(a_1)
     $              +j2_spi(a_2)-m2_spi(a_2)+1-mt2_spi(a_1)
     $              +1-mt2_spi(a_2))/2),kind(0.d0))*pht12
                  
            end do
         end do
      end do
      end subroutine two_nucleon
      subroutine four_nucleon(a_1,a_2,a_3,a_4,perm)
      implicit none
      integer,intent(IN) :: a_1,a_2,a_3,a_4,perm
      integer :: sp_st_1,sp_st_2,sp_st_12,tb_st
      integer :: iT12,J12,sp_st_3,thr_st,sp_st_4,
     +     fou_st,iT123,J123
      real(kind=kind(0.d0)) :: pht12
      real(kind=kind(0.d0)) :: cleb1,cleb2,cleb3,clebt1,clebt2,clebt3

      sp_st_1=iobsind(a_1)
      sp_st_2=iobsind(a_2)
      sp_st_3=iobsind(a_3)
      sp_st_4=iobsind(a_4)
      if (sp_st_1<=sp_st_2) then
         sp_st_12=sp_st_1+sp_st_2*(sp_st_2-1)/2
      else
         sp_st_12=sp_st_2+sp_st_1*(sp_st_1-1)/2
      endif

      do iT123=max(abs(mt2_spi(a_1) !     (123)4 
     +     +mt2_spi(a_2)+mt2_spi(a_3)),
     +     abs(iT1234-1)),min(iT1234+1,3),2
c         clebt3=cleb_3j(6,iT123,mt2_spi(a_1)+mt2_spi(a_2)
c     +        +mt2_spi(a_3),1,mt2_spi(a_4),iT1234)
         clebt3=clebd(iT123,mt2_spi(a_1)+mt2_spi(a_2)
     +        +mt2_spi(a_3),1,mt2_spi(a_4),iT1234,delta_mt)
         if (clebt3==0.d0) cycle
         do J123=max(abs(m2_spi(a_1)+m2_spi(a_2)
     +        +m2_spi(a_3)),abs(J1234-j2_spi(a_4))),
     +        J1234+j2_spi(a_4),2
c            cleb3=cleb_3j(3,J123,m2_spi(a_1)+m2_spi(a_2)
c     +           +m2_spi(a_3),j2_spi(a_4),m2_spi(a_4),
c     +           J1234)
            cleb3=clebd(J123,m2_spi(a_1)+m2_spi(a_2)
     +           +m2_spi(a_3),j2_spi(a_4),m2_spi(a_4),
     +           J1234,delta_mj)
            if (cleb3==0.d0) cycle
                        
            do iT12=max(abs(mt2_spi(a_1) !    (12)34
     +           +mt2_spi(a_2)),abs(iT123-1)),2,2
c               clebt1=cleb_3j(4,1,mt2_spi(a_1),1,mt2_spi(a_2),iT12)
               clebt1=clebd(1,mt2_spi(a_1),1,mt2_spi(a_2),
     +              iT12,mt2_spi(a_1)+mt2_spi(a_2))
               if (clebt1==0.d0) cycle
c               clebt2=cleb_3j(5,iT12,mt2_spi(a_1)+mt2_spi(a_2),
c     +              1,mt2_spi(a_3),iT123)
               clebt2=clebd(iT12,mt2_spi(a_1)+mt2_spi(a_2),
     +              1,mt2_spi(a_3),iT123,mt2_spi(a_1)
     +              +mt2_spi(a_2)+mt2_spi(a_3))
               if (clebt2==0.d0) cycle
               do J12=max(abs(j2_spi(a_1)-j2_spi(a_2)),
     +              abs(m2_spi(a_1)+m2_spi(a_2)),
     +              abs(J123-j2_spi(a_3))),min(j2_spi(a_1)
     +              +j2_spi(a_2),J123+j2_spi(a_3)),2
                  
                  if (sp_st_1==sp_st_2.and.
     +                 (-1)**((J12+iT12)/2)/=-1) cycle
                  
c                  cleb1=cleb_3j(1,j2_spi(a_1),m2_spi(a_1),
c     +                 j2_spi(a_2),m2_spi(a_2),J12)
                  cleb1=clebd(j2_spi(a_1),m2_spi(a_1),
     +                 j2_spi(a_2),m2_spi(a_2),J12,
     +                 m2_spi(a_1)+m2_spi(a_2))
                  if (cleb1==0.d0) cycle
c                  cleb2=cleb_3j(2,J12,m2_spi(a_1)+m2_spi(a_2),
c     +                 j2_spi(a_3),m2_spi(a_3),J123)
                  cleb2=clebd(J12,m2_spi(a_1)+m2_spi(a_2),
     +                 j2_spi(a_3),m2_spi(a_3),J123,
     +                 m2_spi(a_1)+m2_spi(a_2)+m2_spi(a_3))
                  if (cleb2==0.d0) cycle
                  
                  if (sp_st_1<=sp_st_2) then
                     pht12=1.d0
                  else
                     pht12=real((-1)**((J12-j2_spi(a_1)
     +                    -j2_spi(a_2)+iT12)/2),kind(0.d0))
                  endif
                  tb_st=itbind(sp_st_12,J12/2,iT12/2)
                  thr_st=ad3ind(tb_st,sp_st_3,J123/2,
     +                 iT123/2)
                  fou_st=ad4ind(thr_st,sp_st_4,
     +                 J1234/2,iT1234/2)

                  ad_cl_2(fou_st,1,kii)=
     +                 ad_cl_2(fou_st,1,kii)
     +                 +cleb1*cleb2*cleb3
     +                 *clebt1*clebt2*clebt3
     +                 *bmpi(i1,kii)  !*bmpf(i2,kff)
     +                 *real(phase*perm,kind(0.d0))*pht12

               end do
            end do
         end do
      end do
      end subroutine four_nucleon
      end

      subroutine copy_ar(arr_S,dim_S,arr_T,dim_T,i1,i2,i3,i4)
      implicit none
      integer(8),intent(IN) :: i1,i2,i3,i4,dim_S,dim_T
      real(kind(0.d0)),intent(IN) :: arr_S(dim_S)
      real(kind(0.d0)),intent(OUT) :: arr_T(dim_T)
      integer(8) :: i
      do i=i1,i2
         arr_T(i3+i-i1)=arr_S(i)
      end do
      end subroutine copy_ar


      subroutine cluster_overlap_calc_hash
      use initst
      use finast
      use intrface
      use paramdef
      use obdens
      use occ
      use hash_tables
      use nodeinfo
      use cleb_coef
      implicit none
      include 'mpif.h'
      integer :: nhom12max,nhom123max,nhom1234max
      integer :: jtotal2,jtotal2min,jtot2,kii,kff,ii,ia
      integer :: jtotali2,jtotalf2,jtot2min,jt2min,jt2max,
     +     itotali2,itotalf2,j1
      logical :: prot,neut
      integer,allocatable :: occi(:),occf(:),occim2(:)
      integer :: an_st_1,an_st_2,an_st_3,an_st_4,nuc_min,nuc_max,
     +     nuc_min_2,nuc_max_2,nuc_min_3,nuc_max_3,nuc_min_4,nuc_max_4
      integer :: i1,i2,n_state,p_state,ia_1,ia_2,ia_3,ia_4
      integer :: phase,ind1,ierr,N_HO_i,j12,it12
      integer :: iT123,J123,delta_mj,delta_mt,iT1234,J1234
      real(kind=kind(0.d0)) :: cleb,clebd,clebt
      real(kind=kind(0.d0)),allocatable :: adtemp(:,:,:)
      integer :: cr_st_1,cr_st_min,cr_st_max,cr_1,wd,ibit
      integer(8),allocatable :: intbas(:)
      integer,allocatable :: locz(:)
      real(kind(0.d0)) :: aloc_size,tbdval
      integer :: in_mi,in_ma,J12ma,n1,l1,pi,jtrans,Jtr_mi,Jtr_ma,
     $     interm_energy,cr_st_min_1,cr_st_max_1,iphase,mxsps,
     $     Ttr_mi,Ttr_ma,ttrans
      real(kind=kind(0.d0)),allocatable :: t2bdtemp(:)

      print *,' cluster_overlap_calc_hash entered'

      prot=.false.
      neut=.false.

      if (iproc==0) then
         write(iunitvout,*)
         write(iunitvout,7768) ist1dim
 7768    format(' number of single-nucleon states =',i4)
         do ii=1,ist1dim
            write(iunitvout,7770) ii,(ist1(ia,ii),ia=1,3)
 7770       format(' #',i4,'  n=',i3,'  l=',i3,'  j=',i2,'/2')
         end do
      endif   

      select case(iabs(nucleonsi-nucleonsf))
      case(1)
         allocate(ad_cl(ist1dim,kf:kf+nkf-1,ki:ki+nki-1))
         ad_cl=0.d0
         if (nprotonsi==nprotonsf) then
            neut=.true.
         elseif (nneutrnsi==nneutrnsf) then
            prot=.true.
         endif
         if (twobdcal) then
            mxsps=nbit*mxnwdi
            allocate(locz(2*mxsps))
            allocate(intbas(2*mxnwdi))
            allocate(occim2(nucleonsi-2))
            nhom12max=max(nhom12i,nhom12f)
            jtotal2=jt2i(ki)+jt2f(kf)
            jtotal2min=iabs(jt2i(ki)-jt2f(kf))
            do kii=ki,ki+nki-1
               do kff=kf,kf+nkf-1
                  jtot2min=iabs(jt2i(kii)-jt2f(kff))
                  jtot2=jt2i(kii)+jt2f(kff)
                  if (jtotal2min>jtot2min) jtotal2min=jtot2min
                  if (jtotal2<jtot2) jtotal2=jtot2
               end do
            end do 
            jtotal2min=max(0,jtotal2min)
            print *,' 2*jtotal2min,2*jtotal2:',jtotal2min,jtotal2

            allocate(tbd(kf:kf+nkf-1,ki:ki+nki-1))
            aloc_size=0.d0
            do kii=ki,ki+nki-1
               do kff=kf,kf+nkf-1
                  Ttr_mi=max(abs(it2i(kii)-it2f(kff)),1)
                  Ttr_ma=min(it2i(kii)+it2f(kff),3)
                  Jtr_mi=abs(jt2i(kii)-jt2f(kff))
                  Jtr_ma=jt2i(kii)+jt2f(kff)
                  tbd(kff,kii)%Jtr_min=Jtr_mi
                  tbd(kff,kii)%Jtr_max=Jtr_ma
                  allocate(tbd(kff,kii)%Jtr(Jtr_mi/2:Jtr_ma/2))
                  do jtrans=Jtr_mi,Jtr_ma,2
                     tbd(kff,kii)%Jtr(jtrans/2)%dim_fi=ist1dim
                     allocate(tbd(kff,kii)%Jtr(jtrans/2)%fi(ist1dim))
                     do ii=1,ist1dim
                        n1=ist1(1,ii)
                        l1=ist1(2,ii)
                        j1=ist1(3,ii)
                        pi=(-1)**(l1+iparityi+iparityf)
                        if (pi/=(-1)**nhom12max) then
                           J12ma=nhom12max
                        else
                           J12ma=nhom12max+1
                        endif
                        if (abs(jtrans-j1)>ist2_Jst_J12max
     $                       .or.abs(jtrans-j1)<ist2_Jst_J12min
     $                       .or.min(jtrans+j1,J12ma)+1
     $                       >ist2_Jst_J12max+1
     $                       .or.min(jtrans+j1,J12ma)+1
     $                       <ist2_Jst_J12min)
     $                       then
                           in_mi=-1
                           in_ma=-1
                        else
                           in_mi=ist2_Jst((pi+1)/2,abs(jtrans-j1)/2)
                           in_ma=ist2_Jst((pi+1)/2,
     $                          min((jtrans+j1)/2,J12ma)+1)-1
                        endif
c                        in_mi=ist2_Jst((pi+1)/2,abs(jtrans-j1)/2)
c                        in_ma=ist2_Jst((pi+1)/2,min((jtrans+j1)/2,J12ma)
c     $                       +1)-1
                        tbd(kff,kii)%Jtr(jtrans/2)%fi(ii)%mT_min=Ttr_mi
                        tbd(kff,kii)%Jtr(jtrans/2)%fi(ii)%mT_max=Ttr_ma
                        allocate(tbd(kff,kii)%Jtr(jtrans/2)%fi(ii)
     $                       %mT(Ttr_mi/2:Ttr_ma/2))
                        do ttrans=Ttr_mi,Ttr_ma,2
                           tbd(kff,kii)%Jtr(jtrans/2)%fi(ii)
     $                          %mT(ttrans/2)%in_min=in_mi
                           tbd(kff,kii)%Jtr(jtrans/2)%fi(ii)
     $                          %mT(ttrans/2)%in_max=in_ma
                           allocate(tbd(kff,kii)%Jtr(jtrans/2)%fi(ii)
     $                          %mT(ttrans/2)%in(in_mi:in_ma))
                           tbd(kff,kii)%Jtr(jtrans/2)%fi(ii)
     $                          %mT(ttrans/2)%in=0.d0
                           aloc_size=aloc_size+real(in_ma-in_mi+1,
     $                          kind(0.d0))
                        end do
                     end do
                  end do
               end do
            end do
            print *,' iproc=',iproc,' allocated tbd of size=',aloc_size
         endif
      case(2)
         allocate(ad_cl(ist2dim,kf:kf+nkf-1,ki:ki+nki-1))
         ad_cl=0.d0
         if (nprotonsi==nprotonsf) then
            neut=.true.
         elseif (nneutrnsi==nneutrnsf) then
            prot=.true.
         endif

         nhom12max=max(nhom12i,nhom12f)
         jtotal2=(jt2i(ki)+jt2f(kf))/2
         jtotal2min=iabs(jt2i(ki)-jt2f(kf))/2
         do kii=ki,ki+nki-1
            do kff=kf,kf+nkf-1
               jtot2min=iabs(jt2i(kii)-jt2f(kff))/2
               jtot2=(jt2i(kii)+jt2f(kff))/2
               if (jtotal2min>jtot2min) jtotal2min=jtot2min
               if (jtotal2<jtot2) jtotal2=jtot2
            end do
         end do 
         jtotal2min=max(0,jtotal2min)
         jtotal2=min(nhom12max+1,jtotal2)
         print *,' jtotal2min,jtotal2:',jtotal2min,jtotal2

         if (iproc==0) then
            write(iunitvout,*)
            write(iunitvout,"(' Ni_max=',i3,'   Ni_12_max=',i3)") 
     +           nhomi,nhom12i
            write(iunitvout,"(' Nf_max=',i3,'   Nf_12_max=',i3)") 
     +           nhomf,nhom12f
            write(iunitvout,*)
            write(iunitvout,"(' number of two-nucleon states =',i5)") 
     +           ist2dim
            write(iunitvout,"(' J12_min,J12_max=',2i4)") 
     +           jtotal2min,jtotal2
            write(iunitvout,"(' #,i1,i2,J12,T12')")
            do ii=1,ist2dim
               write(iunitvout,"(i5,2i4,2i3)") 
     +              ii,(ist2(ia,ii),ia=1,4) 
            end do
         endif
      case(3)
         allocate(ad_cl(ist3dim,kf:kf+nkf-1,ki:ki+nki-1))
         ad_cl=0.d0
         if (nprotonsi==nprotonsf) then
            neut=.true.
         elseif (nneutrnsi==nneutrnsf) then
            prot=.true.
         endif

         nhom123max=max(nhom123i,nhom123f)
         jtotal2=jt2i(ki)+jt2f(kf)
         jtotal2min=iabs(jt2i(ki)-jt2f(kf))
         do kii=ki,ki+nki-1
            do kff=kf,kf+nkf-1
               jtot2min=iabs(jt2i(kii)-jt2f(kff))
               jtot2=jt2i(kii)+jt2f(kff)
               if (jtotal2min>jtot2min) jtotal2min=jtot2min
               if (jtotal2<jtot2) jtotal2=jtot2
            end do
         end do 
         jtotal2min=max(1,jtotal2min)
         jtotal2=min(2*nhom123max+3,jtotal2)
         print *,' jtotal2min,jtotal2:',jtotal2min,jtotal2

         if (iproc==0) then
            write(iunitvout,*)
            write(iunitvout,"(' Ni_max=',i3,'   Ni_12_max=',i3,
     +           '   Ni_123_max=',i3)") 
     +           nhomi,nhom12i,nhom123i
            write(iunitvout,"(' Nf_max=',i3,'   Nf_12_max=',i3,
     +           '   Nf_123_max=',i3)")  
     +           nhomf,nhom12f,nhom123f
            write(iunitvout,*)
            write(iunitvout,"(' number of two-nucleon states =',i5)") 
     +           ist2dim
            write(iunitvout,"(' #,i1,i2,J12,T12')")
            do ii=1,ist2dim
               write(iunitvout,"(i5,2i4,2i3)") 
     +              ii,(ist2(ia,ii),ia=1,4) 
            end do
            write(iunitvout,*)
            write(iunitvout,"(' number of three-nucleon states =',i7)") 
     +           ist3dim
            write(iunitvout,"(' J123_min,J123_max=',2i4)")
     +           jtotal2min,jtotal2
            write(iunitvout,"(' #,i1,i2,J123,T123=')")
            do ii=1,ist3dim
               write(iunitvout,"(i7,i5,i4,2i3)") 
     +              ii,(ist3(ia,ii),ia=1,4) 
            end do
         endif
      case(4)
         nhom1234max=max(nhom1234i,nhom1234f)
         allocate(ad_cl(ist4dim,kf:kf+nkf-1,ki:ki+nki-1))
         ad_cl=0.d0
         if (nprotonsi==nprotonsf) then
            neut=.true.
         elseif (nneutrnsi==nneutrnsf) then
            prot=.true.
         endif

         jtotal2=(jt2i(ki)+jt2f(kf))/2
         jtotal2min=iabs(jt2i(ki)-jt2f(kf))/2
         do kii=ki,ki+nki-1
            do kff=kf,kf+nkf-1
               jtot2min=iabs(jt2i(kii)-jt2f(kff))/2
               jtot2=(jt2i(kii)+jt2f(kff))/2
               if (jtotal2min>jtot2min) jtotal2min=jtot2min
               if (jtotal2<jtot2) jtotal2=jtot2
            end do
         end do 
         jtotal2min=max(0,jtotal2min)
         jtotal2=min(nhom1234max+2,jtotal2)
         print *,' jtotal2min,jtotal2:',jtotal2min,jtotal2         

         if (iproc==0) then
            write(iunitvout,*)
            write(iunitvout,"(' Ni_max=',i3,'   Ni_12_max=',i3,
     +           '   Ni_123_max=',i3,'   Ni_1234_max=',i3)") 
     +           nhomi,nhom12i,nhom123i,nhom1234i
            write(iunitvout,"(' Nf_max=',i3,'   Nf_12_max=',i3,
     +           '   Nf_123_max=',i3,'   Nf_1234_max=',i3)")  
     +        nhomf,nhom12f,nhom123f,nhom1234f
            write(iunitvout,*)
            write(iunitvout,"(' number of two-nucleon states =',i5)") 
     +           ist2dim
            write(iunitvout,"(' #,i1,i2,J12,T12')")
            do ii=1,ist2dim
               write(iunitvout,"(i5,2i4,2i3)") 
     +              ii,(ist2(ia,ii),ia=1,4) 
            end do
            write(iunitvout,*)
            write(iunitvout,"(' number of three-nucleon states =',i7)") 
     +           ist3dim
            write(iunitvout,"(' #,i1,i2,J123,T123=')")
            do ii=1,ist3dim
               write(iunitvout,"(i7,i5,i4,2i3)") 
     +              ii,(ist3(ia,ii),ia=1,4) 
            end do
            write(iunitvout,*)
            write(iunitvout,"(' number of four-nucleon states =',i8)") 
     +           ist4dim
            write(iunitvout,"(' J1234_min,J1234_max=',2i4)") 
     +           jtotal2min,jtotal2
            write(iunitvout,"(' #,i1,i2,J1234,T1234=')")
            do ii=1,ist4dim
               write(iunitvout,"(i8,i7,i5,2i3)") 
     +              ii,(ist4(ia,ii),ia=1,4) 
            end do
         endif
      case default
         if (iproc==0) then
            write(iunitvout,*) 
     +           ' nucleonsi and nucleonsf out of range',
     +           nucleonsi,nucleonsf
         endif
         stop
      end select

      allocate(occi(nucleonsi))
      allocate(occf(nucleonsf))

c*** main loop for anihilation operators matrix element calculation ***********

      select case(nucleonsi-nucleonsf)
      case(1)

         if (prot) then
            nuc_min=1
            nuc_max=nprotonsi
         else
            nuc_min=nprotonsi+1
            nuc_max=nucleonsi
         endif

         do i1=iproc+1,nsdi,nproc

            if (i1+nproc-iproc-1==1000000*((i1+nproc-iproc-1)/1000000))
     +           then
               print *, '#',iproc,' doing i1=',i1
            endif
         
            if (majortot==2) then
               occi(:)=iloci(:,i1)
            elseif (majortot==3) then
               p_state=I_state_i(1,i1)
               n_state=I_state_i(2,i1)
               occi(1:nprotonsi)=occ_p_i(:,p_state)
               do ii=1,nneutrnsi
                  occi(nprotonsi+ii)=occ_n_i(ii,n_state)+mxnwdi*nbit
               end do
            endif

            N_HO_i=0
            do ia_1=1,nucleonsi
               an_st_1=occi(ia_1)
               N_HO_i=N_HO_i+2*n_spi(an_st_1)+l_spi(an_st_1)
            end do

            do ia_1=nuc_min,nuc_max
               an_st_1=occi(ia_1)
               if (m2_spi(an_st_1)+mjtotalf/=mjtotali) cycle
               if (mod(iparityf+l_spi(an_st_1)+iparityi,2)/=0) cycle
               if (N_HO_i-2*n_spi(an_st_1)-l_spi(an_st_1)>nhwf) cycle

               phase=(-1)**(nucleonsi-ia_1)
               occf(1:ia_1-1)=occi(1:ia_1-1)
               occf(ia_1:nucleonsf)=occi(ia_1+1:nucleonsi)
               call get_state_index(nucleonsf,occf,1,i2)
               if (i2==-1) cycle
               call one_nucleon(an_st_1)

            end do

            if (twobdcal) then
               do ia_1=1,nucleonsi-1
                  an_st_1=occi(ia_1)
                  do ia_2=ia_1+1,nucleonsi
                     an_st_2=occi(ia_2)
                     select case(mttotalf-mttotali+mt2_spi(an_st_1)
     $                    +mt2_spi(an_st_2))
                     case(1)
                        cr_st_min_1=1
                        cr_st_max_1=naspsi
                     case(-1)
                        cr_st_min_1=mxsps+1
                        cr_st_max_1=mxsps+naspsi
                     case default
                        cycle
                     end select
                     if (nucleonsi>2) then
                        occim2(1:ia_1-1)=occi(1:ia_1-1)
                        occim2(ia_1:ia_2-2)=occi(ia_1+1:ia_2-1)
                        occim2(ia_2-1:nucleonsi-2)=occi(ia_2+1:
     $                       nucleonsi)

                        locz(1:occim2(1))=0
                        do cr_1=2,nucleonsi-2
                           locz(occim2(cr_1-1)+1:occim2(cr_1))=cr_1-1
                        end do
                        locz(occim2(nucleonsi-2)+1:2*mxsps)=nucleonsi-2
                        
                        intbas=0
                        do cr_1=1,nucleonsi-2
                           wd=(occim2(cr_1)-1)/nbit+1
                           ibit=mod(occim2(cr_1)-1,nbit)
                           intbas(wd)=ibset(intbas(wd),ibit)
                        end do
                     else
                        intbas=0
                     endif
                     interm_energy=-2*n_spi(an_st_1)-l_spi(an_st_1)
     $                    -2*n_spi(an_st_2)-l_spi(an_st_2)+N_HO_i
                     do cr_st_1=cr_st_min_1,cr_st_max_1

                        if (2*n_spi(cr_st_1)+l_spi(cr_st_1)
     +                       +interm_energy>nhwf) exit
                        if (m2_spi(cr_st_1)
     $                       -m2_spi(an_st_1)-m2_spi(an_st_2)+mjtotali
     +                       /=mjtotalf) cycle
                        if (mod(iparityf+l_spi(cr_st_1)+l_spi(an_st_1)
     $                       +l_spi(an_st_2)+iparityi,2)/=0) cycle

                        wd=(cr_st_1-1)/nbit+1
                        ibit=mod(cr_st_1-1,nbit)
                        if (btest(intbas(wd),ibit)) cycle
                        if (nucleonsi>2) then
                           occf(1:locz(cr_st_1))=occim2(1:locz(cr_st_1))
                           occf(locz(cr_st_1)+1)=cr_st_1
                           occf(locz(cr_st_1)+2:nucleonsf)=
     $                          occim2(locz(cr_st_1)+1:nucleonsi-2)
                        else
                            occf(1)=cr_st_1
                         endif
                         
                        call get_state_index(nucleonsf,occf,1,i2)

                        if (i2==-1) cycle
                        if (nucleonsi>2) then
                           iphase=(-1)**(ia_1+ia_2
     $                          +nucleonsi-locz(cr_st_1)) 
                        else
                           iphase=(-1)**(ia_1+ia_2) 
                        endif
                        
c                        print *,' i1,i2=',i1,i2

                        call one_nucleon_two_body(cr_st_1,an_st_2,
     $                       an_st_1)

                     end do
                  end do
               end do
            endif
         end do 

      case(-1)

         if (prot) then
            cr_st_min=1
            cr_st_max=mxnwdi*nbit
         else
            cr_st_min=mxnwdi*nbit+1
            cr_st_max=2*mxnwdi*nbit
         endif

c         print *,' cr_st_min,cr_st_max=',cr_st_min,cr_st_max

         allocate(intbas(2*mxnwdi))
         allocate(locz(2*mxnwdi*nbit))
         do i1=iproc+1,nsdi,nproc

            if (i1+nproc-iproc-1==1000000*((i1+nproc-iproc-1)/1000000))
     +           then
               print *, '#',iproc,' doing i1=',i1
            endif
         
            if (majortot==2) then
               occi(:)=iloci(:,i1)
            elseif (majortot==3) then
               p_state=I_state_i(1,i1)
               n_state=I_state_i(2,i1)
               occi(1:nprotonsi)=occ_p_i(:,p_state)
               do ii=1,nneutrnsi
                  occi(nprotonsi+ii)=occ_n_i(ii,n_state)+mxnwdi*nbit
               end do
            endif
            locz(1:occi(1))=0
            do cr_1=2,nucleonsi
               locz(occi(cr_1-1)+1:occi(cr_1))=cr_1-1
            end do
            locz(occi(nucleonsi)+1:2*mxnwdi*nbit)=nucleonsi

c            print *,' occi=',occi
c            print *,' locz=',locz

            N_HO_i=0
            do cr_1=1,nucleonsi
               cr_st_1=occi(cr_1)
               N_HO_i=N_HO_i+2*n_spi(cr_st_1)+l_spi(cr_st_1)
            end do

c            print *,' N_HO_i=',N_HO_i

            intbas=0
            do cr_1=1,nucleonsi
               wd=(occi(cr_1)-1)/nbit+1
               ibit=mod(occi(cr_1)-1,nbit)
c               print *,' i,wd,ibit=',i,wd,ibit
               intbas(wd)=ibset(intbas(wd),ibit)
            end do
            
c            do cr_1=1,wd
c               print *,' wd=',cr_1
c               print '(b64)',intbas(cr_1)
c            end do

            do cr_st_1=cr_st_min,cr_st_max
               wd=(cr_st_1-1)/nbit+1
               ibit=mod(cr_st_1-1,nbit)
               if (btest(intbas(wd),ibit)) cycle
               if (m2_spi(cr_st_1)+mjtotali/=mjtotalf) cycle
               if (mod(iparityf+l_spi(cr_st_1)+iparityi,2)/=0) cycle
               if (2*n_spi(cr_st_1)+l_spi(cr_st_1)+N_HO_i>nhwf) exit
               occf(1:locz(cr_st_1))=occi(1:locz(cr_st_1))
               occf(locz(cr_st_1)+1)=cr_st_1
               occf(locz(cr_st_1)+2:nucleonsf)=
     +              occi(locz(cr_st_1)+1:nucleonsi)
               phase=(-1)**(nucleonsi-locz(cr_st_1))
               call get_state_index(nucleonsf,occf,1,i2)
               if (i2==-1) cycle
               call one_nucleon(cr_st_1)
               
            end do
         end do
         deallocate(intbas)
         deallocate(locz)

      case(2)

         if (prot) then
            nuc_min=1
            nuc_max=nprotonsi-1
            nuc_min_2=2
            nuc_max_2=nprotonsi
         elseif (neut) then
            nuc_min=nprotonsi+1
            nuc_max=nucleonsi-1
            nuc_min_2=nprotonsi+2
            nuc_max_2=nucleonsi
         else
            nuc_min=1
            nuc_max=nprotonsi
            nuc_min_2=nprotonsi+1
            nuc_max_2=nucleonsi
         endif

         do i1=iproc+1,nsdi,nproc

            if (i1+nproc-iproc-1==1000000*((i1+nproc-iproc-1)/1000000))
     +           then
               print *, '#',iproc,' doing i1=',i1
            endif
         
            if (majortot==2) then
               occi(:)=iloci(:,i1)
            elseif (majortot==3) then
               p_state=I_state_i(1,i1)
               n_state=I_state_i(2,i1)
               occi(1:nprotonsi)=occ_p_i(:,p_state)
               do ii=1,nneutrnsi
                  occi(nprotonsi+ii)=occ_n_i(ii,n_state)+mxnwdi*nbit
               end do
            endif

            N_HO_i=0
            do ia_1=1,nucleonsi
               an_st_1=occi(ia_1)
               N_HO_i=N_HO_i+2*n_spi(an_st_1)+l_spi(an_st_1)
            end do

            do ia_1=nuc_min,nuc_max
               an_st_1=occi(ia_1)
               do ia_2=max(nuc_min_2,ia_1+1),nuc_max_2
                  an_st_2=occi(ia_2)
                  if (m2_spi(an_st_1)+m2_spi(an_st_2)+mjtotalf
     +                 /=mjtotali) cycle
                  if (mod(iparityf+l_spi(an_st_1)+l_spi(an_st_2)
     +                 +iparityi,2)/=0) cycle
                  if (N_HO_i-2*n_spi(an_st_1)-l_spi(an_st_1)
     +                 -2*n_spi(an_st_2)-l_spi(an_st_2)>nhwf) cycle
                  phase=(-1)**(ia_1+ia_2)

                  occf(1:ia_1-1)=occi(1:ia_1-1)
                  occf(ia_1:ia_2-2)=occi(ia_1+1:ia_2-1)
                  occf(ia_2-1:nucleonsf)=occi(ia_2+1:nucleonsi)

                  call get_state_index(nucleonsf,occf,1,i2)
                  if (i2==-1) cycle
                  call two_nucleon(an_st_1,an_st_2)

               end do
            end do

         end do

      case(3)

         call num_of_3j_types_init(4)
c         print *,' num_of_3j_types_init allocated'
         call cleb_init(1,1,2*nhomi+1,-2*nhomi-1,2*nhomi+1,
     $        1,2*nhomi+1,-2*nhomi-1,2*nhomi+1,0,2*nhom12i+2)
c         print *,' cleb_init(1.. allocated'
         call cleb_init(2,0,2*nhom12i+2,-2*nhom12i-2,2*nhom12i+2,
     $        1,2*nhomi+1,-2*nhomi-1,2*nhomi+1,1,2*nhom123i+3)
c         print *,' cleb_init(2.. allocated'
         call cleb_init(3,1,1,-1,1,1,1,-1,1,0,2)
c         print *,' cleb_init(3.. allocated'
         call cleb_init(4,0,2,-2,2,1,1,-1,1,1,3)
c         print *,' cleb_init(4.. allocated'

         if (prot) then
            nuc_min=1
            nuc_max=nprotonsi-2
            nuc_min_2=2
            nuc_max_2=nprotonsi-1
            nuc_min_3=3
            nuc_max_3=nprotonsi
         elseif (neut) then
            nuc_min=nprotonsi+1
            nuc_max=nucleonsi-2
            nuc_min_2=nprotonsi+2
            nuc_max_2=nucleonsi-1
            nuc_min_3=nprotonsi+3
            nuc_max_3=nucleonsi
         elseif (nprotonsf==nprotonsi-2) then
            nuc_min=1
            nuc_max=nprotonsi-1
            nuc_min_2=2
            nuc_max_2=nprotonsi
            nuc_min_3=nprotonsi+1
            nuc_max_3=nucleonsi
         else
            nuc_min=1
            nuc_max=nprotonsi
            nuc_min_2=nprotonsi+1
            nuc_max_2=nucleonsi-1
            nuc_min_3=nprotonsi+2
            nuc_max_3=nucleonsi
         endif

         do i1=iproc+1,nsdi,nproc

            if (i1+nproc-iproc-1==1000000*((i1+nproc-iproc-1)/1000000))
     +           then
               print *, '#',iproc,' doing i1=',i1
            endif
         
            if (majortot==2) then
               occi(:)=iloci(:,i1)
            elseif (majortot==3) then
               p_state=I_state_i(1,i1)
               n_state=I_state_i(2,i1)
               occi(1:nprotonsi)=occ_p_i(:,p_state)
               do ii=1,nneutrnsi
                  occi(nprotonsi+ii)=occ_n_i(ii,n_state)+mxnwdi*nbit
               end do
            endif

            N_HO_i=0
            do ia_1=1,nucleonsi
               an_st_1=occi(ia_1)
               N_HO_i=N_HO_i+2*n_spi(an_st_1)+l_spi(an_st_1)
            end do

            do ia_1=nuc_min,nuc_max
               an_st_1=occi(ia_1)
               do ia_2=max(nuc_min_2,ia_1+1),nuc_max_2
                  an_st_2=occi(ia_2)
                  do ia_3=max(nuc_min_3,ia_2+1),nuc_max_3
                     an_st_3=occi(ia_3)
                     if (m2_spi(an_st_1)+m2_spi(an_st_2)
     +                    +m2_spi(an_st_3)+mjtotalf/=mjtotali) cycle
                     if (mod(iparityf+l_spi(an_st_1)+l_spi(an_st_2)
     +                    +l_spi(an_st_3)+iparityi,2)/=0) cycle
                     if (N_HO_i-2*n_spi(an_st_1)-l_spi(an_st_1)
     +                    -2*n_spi(an_st_2)-l_spi(an_st_2)
     +                    -2*n_spi(an_st_3)-l_spi(an_st_3)>nhwf) cycle
                     phase=(-1)**(ia_1+ia_2+ia_3-nucleonsi)

                     occf(1:ia_1-1)=occi(1:ia_1-1)
                     occf(ia_1:ia_2-2)=occi(ia_1+1:ia_2-1)
                     occf(ia_2-1:ia_3-3)=occi(ia_2+1:ia_3-1)
                     occf(ia_3-2:nucleonsf)=occi(ia_3+1:nucleonsi)

                     call get_state_index(nucleonsf,occf,1,i2)
                     if (i2==-1) cycle

                     delta_mj=m2_spi(an_st_1)+m2_spi(an_st_2)
     +                    +m2_spi(an_st_3)
                     delta_mt=mt2_spi(an_st_1)+mt2_spi(an_st_2)
     +                    +mt2_spi(an_st_3)
                     do kii=ki,ki+nki-1
                        do kff=kf,kf+nkf-1
                           jt2min=abs(jt2i(kii)-jt2f(kff))
                           jt2max=jt2i(kii)+jt2f(kff)
                           do iT123=abs(delta_mt),3,2
                              do J123=max(jtotal2min,abs(delta_mj),
     +                             jt2min),min(jtotal2,jt2max),2

                                 call three_nucleon(an_st_1,an_st_2,
     +                                an_st_3,1)
                                 call three_nucleon(an_st_1,an_st_3,
     +                                an_st_2,-1)
                                 call three_nucleon(an_st_2,an_st_3,
     +                                an_st_1,1)

                              end do
                           end do
                        end do
                     end do

                  end do
               end do
            end do
         end do
         
         call cleb_destroy(1)
         call cleb_destroy(2)
         call cleb_destroy(3)
         call cleb_destroy(4)
         call num_of_3j_types_destroy

      case(4)

         call num_of_3j_types_init(6)
         call cleb_init(1,1,2*nhomi+1,-2*nhomi-1,2*nhomi+1,
     $        1,2*nhomi+1,-2*nhomi-1,2*nhomi+1,0,2*nhom12i+2)
         call cleb_init(2,0,2*nhom12i+2,-2*nhom12i-2,2*nhom12i+2,
     $        1,2*nhomi+1,-2*nhomi-1,2*nhomi+1,1,2*nhom123i+3)
         call cleb_init(3,1,2*nhom123i+3,-2*nhom123i-3,2*nhom123i+3,
     $        1,2*nhomi+1,-2*nhomi-1,2*nhomi+1,0,2*nhom1234i+4)
         call cleb_init(4,1,1,-1,1,1,1,-1,1,0,2)
         call cleb_init(5,0,2,-2,2,1,1,-1,1,1,3)
         call cleb_init(6,1,3,-3,3,1,1,-1,1,0,4)

         if (prot) then
            nuc_min=1
            nuc_max=nprotonsi-3
            nuc_min_2=2
            nuc_max_2=nprotonsi-2
            nuc_min_3=3
            nuc_max_3=nprotonsi-1
            nuc_min_4=4
            nuc_max_3=nprotonsi
         elseif (neut) then
            nuc_min=nprotonsi+1
            nuc_max=nucleonsi-3
            nuc_min_2=nprotonsi+2
            nuc_max_2=nucleonsi-2
            nuc_min_3=nprotonsi+3
            nuc_max_3=nucleonsi-1
            nuc_min_4=nprotonsi+4
            nuc_max_4=nucleonsi
         elseif (nprotonsf==nprotonsi-3) then
            nuc_min=1
            nuc_max=nprotonsi-2
            nuc_min_2=2
            nuc_max_2=nprotonsi-1
            nuc_min_3=3
            nuc_max_3=nprotonsi
            nuc_min_4=nprotonsi+1
            nuc_max_4=nucleonsi
         elseif (nprotonsf==nprotonsi-2) then
            nuc_min=1
            nuc_max=nprotonsi-1
            nuc_min_2=2
            nuc_max_2=nprotonsi
            nuc_min_3=nprotonsi+1
            nuc_max_3=nucleonsi-1
            nuc_min_4=nprotonsi+2
            nuc_max_4=nucleonsi
         elseif (nprotonsf==nprotonsi-1) then
            nuc_min=1
            nuc_max=nprotonsi
            nuc_min_2=nprotonsi+1
            nuc_max_2=nucleonsi-2
            nuc_min_3=nprotonsi+2
            nuc_max_3=nucleonsi-1
            nuc_min_3=nprotonsi+3
            nuc_max_3=nucleonsi
         else
            print *,'*** error: cluster_ case(4)'
            stop
         endif

         do i1=iproc+1,nsdi,nproc

            if (i1+nproc-iproc-1==1000000*((i1+nproc-iproc-1)/1000000))
     +           then
               print *, '#',iproc,' doing i1=',i1
            endif
         
            if (majortot==2) then
               occi(:)=iloci(:,i1)
            elseif (majortot==3) then
               p_state=I_state_i(1,i1)
               n_state=I_state_i(2,i1)
               occi(1:nprotonsi)=occ_p_i(:,p_state)
               do ii=1,nneutrnsi
                  occi(nprotonsi+ii)=occ_n_i(ii,n_state)+mxnwdi*nbit
               end do
            endif

c            if (i1==1) print *,' occi=',occi

            N_HO_i=0
            do ia_1=1,nucleonsi
               an_st_1=occi(ia_1)
               N_HO_i=N_HO_i+2*n_spi(an_st_1)+l_spi(an_st_1)
c               if (i1==1) print *,' ia,st,n,l=',ia_1,an_st_1,
c     $              n_spi(an_st_1),l_spi(an_st_1)
            end do

c            if (i1==1) print *,' N_HO_i=',N_HO_i

            do ia_1=nuc_min,nuc_max
               an_st_1=occi(ia_1)
               do ia_2=max(nuc_min_2,ia_1+1),nuc_max_2
                  an_st_2=occi(ia_2)
                  do ia_3=max(nuc_min_3,ia_2+1),nuc_max_3
                     an_st_3=occi(ia_3)
                     do ia_4=max(nuc_min_4,ia_3+1),nuc_max_4
                        an_st_4=occi(ia_4)

c                        if (i1==1) print *,'ia=',ia_1,ia_2,ia_3,ia_4
c                        if (i1==1) print *,'an_st=',
c     $                       an_st_1,an_st_2,an_st_3,an_st_4

                        if (m2_spi(an_st_1)+m2_spi(an_st_2)
     +                       +m2_spi(an_st_3)+m2_spi(an_st_4)
     +                       +mjtotalf/=mjtotali) cycle
                        if (mod(iparityf+l_spi(an_st_1)+l_spi(an_st_2)
     +                       +l_spi(an_st_3)+l_spi(an_st_4)
     +                       +iparityi,2)/=0) cycle
                        if (N_HO_i-2*n_spi(an_st_1)-l_spi(an_st_1)
     +                       -2*n_spi(an_st_2)-l_spi(an_st_2)
     +                       -2*n_spi(an_st_3)-l_spi(an_st_3)
     +                       -2*n_spi(an_st_4)-l_spi(an_st_4)>nhwf)
     +                       cycle
                        phase=(-1)**(ia_1+ia_2+ia_3+ia_4)

                        occf(1:ia_1-1)=occi(1:ia_1-1)
                        occf(ia_1:ia_2-2)=occi(ia_1+1:ia_2-1)
                        occf(ia_2-1:ia_3-3)=occi(ia_2+1:ia_3-1)
                        occf(ia_3-2:ia_4-4)=occi(ia_3+1:ia_4-1)
                        occf(ia_4-3:nucleonsf)=occi(ia_4+1:nucleonsi)

c                        if (i1==1) print *,' occi=',occi
c                        if (i1==1) print *,' occf=',occf
                        call get_state_index(nucleonsf,occf,1,i2)
c                        if (i1==1) print *,' i2=',i2
                        if (i2==-1) cycle

                        delta_mj=m2_spi(an_st_1)+m2_spi(an_st_2)
     +                       +m2_spi(an_st_3)+m2_spi(an_st_4)
                        delta_mt=mt2_spi(an_st_1)+mt2_spi(an_st_2)
     +                       +mt2_spi(an_st_3) +mt2_spi(an_st_4)

c                        if (i1==1) print *,' delta_mj,delta_mt=',
c     $                       delta_mj,delta_mt

                        do kii=ki,ki+nki-1
                           do kff=kf,kf+nkf-1
                              jt2min=abs(jt2i(kii)-jt2f(kff))
                              jt2max=jt2i(kii)+jt2f(kff)
                              do iT1234=abs(delta_mt),4,2
                                 do J1234=max(2*jtotal2min,
     +                                abs(delta_mj),jt2min),
     +                                min(2*jtotal2,jt2max),2

                                    call four_nucleon(an_st_1,an_st_2,
     +                                   an_st_3,an_st_4,1)
                                    call four_nucleon(an_st_1,an_st_3,
     +                                   an_st_2,an_st_4,-1)
                                    call four_nucleon(an_st_2,an_st_3,
     +                                   an_st_1,an_st_4,1)
                                    call four_nucleon(an_st_1,an_st_2,
     +                                   an_st_4,an_st_3,-1)
                                    call four_nucleon(an_st_1,an_st_4,
     +                                   an_st_2,an_st_3,1)
                                    call four_nucleon(an_st_2,an_st_4,
     +                                   an_st_1,an_st_3,-1)
                                    call four_nucleon(an_st_1,an_st_4,
     +                                   an_st_3,an_st_2,-1)
                                    call four_nucleon(an_st_1,an_st_3,
     +                                   an_st_4,an_st_2,1)
                                    call four_nucleon(an_st_3,an_st_4,
     +                                   an_st_1,an_st_2,1)
                                    call four_nucleon(an_st_2,an_st_4,
     +                                   an_st_3,an_st_1,1)
                                    call four_nucleon(an_st_2,an_st_3,
     +                                   an_st_4,an_st_1,-1)
                                    call four_nucleon(an_st_3,an_st_4,
     +                                   an_st_2,an_st_1,-1)
                                 end do
                              end do
                           end do
                        end do
                        
                     end do
                  end do
               end do
            end do
         end do

         call cleb_destroy(1)
         call cleb_destroy(2)
         call cleb_destroy(3)
         call cleb_destroy(4)
         call cleb_destroy(5)
         call cleb_destroy(6)
         call num_of_3j_types_destroy

      case default
         print *,'*** error: not implemented Ai-Af=',
     +        nucleonsi-nucleonsf
         call MPI_Abort(icomm,112,ierr)
         stop
      end select

      deallocate(occi,occf)

      call MPI_Barrier(icomm,ierr)
      if (nproc>1) then
         select case(iabs(nucleonsi-nucleonsf))
         case(1)
            ind1=ist1dim
         case(2)
            ind1=ist2dim
         case(3)
            ind1=ist3dim
         case(4)
            ind1=ist4dim
         case default
            if (iproc==0) then
               write(iunitvout,*) 
     +              ' nucleonsi and nucleonsf out of range',
     +              nucleonsi,nucleonsf
            endif
            stop
         end select
         allocate(adtemp(ind1,kf:kf+nkf-1,ki:ki+nki-1))
         ii=nkf*nki*ind1
         call MPI_Reduce(ad_cl(1,kf,ki),adtemp,ii,
     +        MPI_REAL8,MPI_SUM,0,icomm,ierr)
         ad_cl=adtemp
         deallocate(adtemp)
         if (nucleonsi==nucleonsf+1.and.twobdcal) then
            call MPI_Barrier(icomm,ierr)
            allocate(t2bdtemp(ist2dim))
            do kii=ki,ki+nki-1
               do kff=kf,kf+nkf-1
                  Ttr_mi=max(abs(it2i(kii)-it2f(kff)),1)
                  Ttr_ma=min(it2i(kii)+it2f(kff),3)
                  Jtr_mi=tbd(kff,kii)%Jtr_min
                  Jtr_ma=tbd(kff,kii)%Jtr_max
                  do jtrans=Jtr_mi,Jtr_ma,2
                     do ii=1,ist1dim
                        do ttrans=Ttr_mi,Ttr_ma,2
                           in_mi=tbd(kff,kii)%Jtr(jtrans/2)%fi(ii)
     $                          %mT(ttrans/2)%in_min
                           in_ma=tbd(kff,kii)%Jtr(jtrans/2)%fi(ii)
     $                          %mT(ttrans/2)%in_max
                           call MPI_Reduce(tbd(kff,kii)%Jtr(jtrans/2)
     $                          %fi(ii)%mT(ttrans/2)%in(in_mi:in_ma),
     $                          t2bdtemp,in_ma-in_mi+1,
     +                          MPI_REAL8,MPI_SUM,0,icomm,ierr)
                           if (iproc==0) tbd(kff,kii)%Jtr(jtrans/2)
     $                          %fi(ii)%mT(ttrans/2)%in(in_mi:in_ma)=
     $                          t2bdtemp(1:in_ma-in_mi+1)
                        end do
                     end do
                  end do
               end do
            end do
         endif
      endif          
      if (iproc>0) return
      
      do kii=ki,ki+nki-1

         jtotali2=jt2i(kii)
         itotali2=it2i(kii)
         if (dabs(dble(jtotali2)-2.d0*Jxi(kii))>1.d-4) then
            print *, '**** not a good J ****'
         endif 

         do kff=kf,kf+nkf-1

            jtotalf2=jt2f(kff)
            itotalf2=it2f(kff)
            if (dabs(dble(jtotalf2)-2.d0*Jxf(kff))>1.d-4) then
               print *, '**** not a good J ****'
            endif 

            write(iunitvout,1200) kff,jtotalf2,itotalf2,
     +           enerf(kff)-eneri(ki),
     +           kii,jtotali2,itotali2,eneri(kii)-eneri(ki)
 1200       format(//,' *** Transition matrix elements for states:',
     +           ' ***',/,' #',i3,
     +           ' [2*(J,T),Ex]_f=',i3,i2,f8.4,3x,
     +           '#',i3,' [2*(J,T),Ex]_i=',i3,i2,f8.4) 
c*1200       format(//,' *** Transition matrix elements for states:',
c*     +           ' ***',/,i3,
c*     +           '   2*J_f=',i3,'   Ex_f=',f9.4,5x,
c*     +           i3,'   2*J_i=',i3,'   Ex_i=',f9.4)

            select case(iabs(nucleonsi-nucleonsf))
            case(1)
               if (nucleonsi==nucleonsf+1) then
                  write(iunitvout,*)
                  do i1=1,ist1dim
                     j1=ist1(3,i1)    
                     cleb=clebd(jtotalf2,mjtotalf,j1,
     +                    mjtotali-mjtotalf,jtotali2,mjtotali)
cc                     if (cleb==0.d0) cycle
                     if (abs(cleb)<1.d-12) cycle
                     clebt=clebd(itotalf2,mttotalf,1,
     +                    mttotali-mttotalf,itotali2,mttotali) 
cc                     if (clebt==0.d0) cycle
                     if (abs(clebt)<1.d-12) cycle
                     ad_cl(i1,kff,kii)=ad_cl(i1,kff,kii)/cleb/clebt
                     write(iunitvout,7960) i1,ad_cl(i1,kff,kii)
 7960                format(' a-=',i3,'        <a->=',f10.6)                  
                  
                  end do   
               elseif (nucleonsi==nucleonsf-1) then
                  write(iunitvout,*)
                  do i1=1,ist1dim
                     j1=ist1(3,i1)    
                     cleb=clebd(jtotali2,mjtotali,j1,
     +                    mjtotalf-mjtotali,jtotalf2,mjtotalf)
cc                     if (cleb==0.d0) cycle
                     if (abs(cleb)<1.d-12) cycle
                     clebt=clebd(itotali2,mttotali,1,
     +                    mttotalf-mttotali,itotalf2,mttotalf) 
cc                     if (clebt==0.d0) cycle
                     if (abs(clebt)<1.d-12) cycle
                     ad_cl(i1,kff,kii)=ad_cl(i1,kff,kii)/cleb/clebt
                     write(iunitvout,7962) i1,ad_cl(i1,kff,kii)
 7962                format(' a+=',i3,'        <a+>=',f10.6)                  
                  end do   
               endif   
            case(2)
               if (nucleonsi==nucleonsf+2) then
                  write(iunitvout,*)
                  do i1=1,ist2dim
                     j12=ist2(3,i1)    
                     cleb=clebd(jtotalf2,mjtotalf,2*j12,
     +                    mjtotali-mjtotalf,jtotali2,mjtotali)
cc                     if (cleb==0.d0) cycle
                     if (abs(cleb)<1.d-12) cycle
                     it12=ist2(4,i1)    
                     clebt=clebd(itotalf2,mttotalf,2*it12,
     +                    mttotali-mttotalf,itotali2,mttotali) 
cc                     if (clebt==0.d0) cycle
                     if (abs(clebt)<1.d-12) cycle
                     ad_cl(i1,kff,kii)=ad_cl(i1,kff,kii)/cleb/clebt*2.d0  !2 for j_a,j_b exchange
                     write(iunitvout,79602) i1,ad_cl(i1,kff,kii)
79602                format(' a-a-=',i5,'     <a-a->=',f10.6)                  
                  end do   
               elseif (nucleonsi==nucleonsf-2) then
                  write(iunitvout,*)
                  do i1=1,ist2dim
                     j12=ist2(3,i1)    
                     cleb=clebd(jtotali2,mjtotali,2*j12,
     +                    mjtotalf-mjtotali,jtotalf2,mjtotalf)
cc                     if (cleb==0.d0) cycle
                     if (abs(cleb)<1.d-12) cycle
                     it12=ist2(4,i1)    
                     clebt=clebd(itotali2,mttotali,2*it12,
     +                    mttotalf-mttotali,itotalf2,mttotalf) 
cc                     if (clebt==0.d0) cycle
                     if (abs(clebt)<1.d-12) cycle
                     ad_cl(i1,kff,kii)=ad_cl(i1,kff,kii)/cleb/clebt*2.d0  !2 for j_a,j_b exchange
                     write(iunitvout,79622) i1,ad_cl(i1,kff,kii)
79622                format(' a+a+=',i5,'     <a+a+>=',f10.6)                  
                  end do
               endif
            case(3)
               if (nucleonsi==nucleonsf+3) then
                  write(iunitvout,*)
                  do i1=1,ist3dim
                     j12=ist3(3,i1)    
                     cleb=clebd(jtotalf2,mjtotalf,j12,
     +                    mjtotali-mjtotalf,jtotali2,mjtotali)
cc                     if (cleb==0.d0) cycle
                     if (abs(cleb)<1.d-12) cycle
                     it12=ist3(4,i1)    
                     clebt=clebd(itotalf2,mttotalf,it12,
     +                    mttotali-mttotalf,itotali2,mttotali) 
cc                     if (clebt==0.d0) cycle
                     if (abs(clebt)<1.d-12) cycle
                     ad_cl(i1,kff,kii)=ad_cl(i1,kff,kii)/cleb/clebt*2.d0  !2 for j_a,j_b exchange
                     write(iunitvout,79603) i1,ad_cl(i1,kff,kii)
79603                format(' a-a-a-=',i7,'     <a-a-a->=',f10.6)                  
                  end do   
               elseif (nucleonsi==nucleonsf-3) then
                  write(iunitvout,*)
                  do i1=1,ist3dim
                     j12=ist3(3,i1)    
                     cleb=clebd(jtotali2,mjtotali,j12,
     +                    mjtotalf-mjtotali,jtotalf2,mjtotalf)
                     if (abs(cleb)<1.d-12) cycle
                     it12=ist3(4,i1)    
                     clebt=clebd(itotali2,mttotali,it12,
     +                    mttotalf-mttotali,itotalf2,mttotalf) 
cc                     if (clebt==0.d0) cycle
                     if (abs(clebt)<1.d-12) cycle
                     ad_cl(i1,kff,kii)=ad_cl(i1,kff,kii)/cleb/clebt*2.d0  !2 for j_a,j_b exchange
                     write(iunitvout,79633) i1,ad_cl(i1,kff,kii)
79633                format(' a+a+a+=',i7,'     <a+a+a+>=',f10.6)                  
                  end do
               endif               
            case(4)
               if (nucleonsi==nucleonsf+4) then
                  write(iunitvout,*)
                  do i1=1,ist4dim
                     j12=ist4(3,i1)    
                     cleb=clebd(jtotalf2,mjtotalf,2*j12,
     +                    mjtotali-mjtotalf,jtotali2,mjtotali)

c                     if (i1==212429) then
c                        print *,' i1=',i1
c                        print *,
c     +                       ' jtotalf2,mjtotalf,j12,jtotali2,mjtotali'
c     +                       ,jtotalf2,mjtotalf,j12,jtotali2,mjtotali
c                        print *,' cleb=',cleb
c                        print *,' kff,kii=',kff,kii
c                        print *,' ad_cl=',
c     +                       ad_cl(i1,kff,kii)
c                     endif

cc                     if (cleb==0.d0) cycle
                     if (abs(cleb)<1.d-12) cycle
                     it12=ist4(4,i1)    
                     clebt=clebd(itotalf2,mttotalf,2*it12,
     +                    mttotali-mttotalf,itotali2,mttotali) 

c                     if (i1==212429) then
c                        print *,' it12=',it12
c                        print *,
c     +                       ' itotalf2,mttotalf,it12,itotali2,mttotali'
c     +                       ,itotalf2,mttotalf,it12,itotali2,mttotali
c                        print *,' clebt=',clebt
c                     endif

cc                     if (clebt==0.d0) cycle
                     if (abs(clebt)<1.d-12) cycle
                     ad_cl(i1,kff,kii)=ad_cl(i1,kff,kii)/cleb/clebt*2.d0  !2 for j_a,j_b exchange
                     write(iunitvout,79604) i1,ad_cl(i1,kff,kii)
79604                format(' a-a-a-a-=',i8,'     <a-a-a-a->=',f10.6)                  
                  end do   
               elseif (nucleonsi==nucleonsf-4) then
                  write(iunitvout,*)
                  do i1=1,ist4dim
                     j12=ist4(3,i1)    
                     cleb=clebd(jtotali2,mjtotali,2*j12,
     +                    mjtotalf-mjtotali,jtotalf2,mjtotalf)
cc                     if (cleb==0.d0) cycle
                     if (abs(cleb)<1.d-12) cycle
                     it12=ist4(4,i1)    
                     clebt=clebd(itotali2,mttotali,2*it12,
     +                    mttotalf-mttotali,itotalf2,mttotalf) 
cc                     if (clebt==0.d0) cycle
                     if (abs(clebt)<1.d-12) cycle
                     ad_cl(i1,kff,kii)=ad_cl(i1,kff,kii)/cleb/clebt*2.d0  !2 for j_a,j_b exchange
                     write(iunitvout,79644) i1,ad_cl(i1,kff,kii)
79644                format(' a+a+a+a+=',i8,'     <a+a+a+a+>=',f10.6)                  
                  end do
               endif
            case default
               if (iproc==0) then
                  write(iunitvout,*) 
     +                 ' nucleonsi and nucleonsf out of range',
     +                 nucleonsi,nucleonsf
               endif
               stop
            end select

         end do
      end do
      
      if (nucleonsi==nucleonsf+1.and.twobdcal) then
         write(iunitvout,*)
         write(iunitvout,"(' number of two-nucleon states =',i6)") 
     +        ist2dim
         write(iunitvout,"(' 2*j_tr_min,2*j_tr_max=',2i4)") 
     +        jtotal2min,jtotal2
         write(iunitvout,"(' #,i1,i2,J12,T12')")
         do ii=1,ist2dim
            write(iunitvout,"(i6,2i4,2i3)") 
     +           ii,(ist2(ia,ii),ia=1,4) 
         end do
         do kii=ki,ki+nki-1
            jtotali2=jt2i(kii)
            itotali2=it2i(kii)
            do kff=kf,kf+nkf-1
               jtotalf2=jt2f(kff)
               itotalf2=it2f(kff)
               write(iunitvout,1200) kff,jtotalf2,itotalf2,
     +              enerf(kff)-eneri(ki),
     +              kii,jtotali2,itotali2,eneri(kii)-eneri(ki)

               Ttr_mi=max(abs(it2i(kii)-it2f(kff)),1)
               Ttr_ma=min(it2i(kii)+it2f(kff),3)
               Jtr_mi=tbd(kff,kii)%Jtr_min
               Jtr_ma=tbd(kff,kii)%Jtr_max
               do jtrans=Jtr_mi,Jtr_ma,2
                  cleb=clebd(jtotalf2,mjtotalf,jtrans,
     +                 mjtotali-mjtotalf,jtotali2,mjtotali)
     $                 /sqrt(real(jtotali2+1,kind(0.d0)))

                  if (abs(cleb)<1.d-9) cycle

                  write(iunitvout,1953) jtrans
 1953             format(/,' 2*Jtrans=',i3)

                  do ttrans=Ttr_mi,Ttr_ma,2
                     clebt=clebd(itotalf2,mttotalf,ttrans,
     +                    mttotali-mttotalf,itotali2,mttotali)
     $                    /sqrt(real(itotali2+1,kind(0.d0)))
                     if (abs(clebt)<1.d-9) cycle

                     write(iunitvout,1954) ttrans
 1954                format(/,' 2*Ttrans=',i3)

                     ii=0
                     do i1=1,ist1dim
                        in_mi=tbd(kff,kii)%Jtr(jtrans/2)%fi(i1)
     $                       %mT(ttrans/2)%in_min
                        in_ma=tbd(kff,kii)%Jtr(jtrans/2)%fi(i1)
     $                       %mT(ttrans/2)%in_max

c                        print *,' i1,in_mi,in_ma=',i1,in_mi,in_ma

                        do i2=in_mi,in_ma
                           tbdval=tbd(kff,kii)%Jtr(jtrans/2)
     $                          %fi(i1)%mT(ttrans/2)%in(i2)
                           if (abs(tbdval)>1.d-9) then

c                              print *,' before ii=',ii
                              ii=ii+1
c                              print *,' after ii=',ii

                           endif
                        end do
                     end do
                     write(iunitvout,'(i10)') ii
                     ia=0
                     do i1=1,ist1dim

c                        print *,' i1,ia=',i1,ia

                        in_mi=tbd(kff,kii)%Jtr(jtrans/2)%fi(i1)
     $                       %mT(ttrans/2)%in_min
                        in_ma=tbd(kff,kii)%Jtr(jtrans/2)%fi(i1)
     $                       %mT(ttrans/2)%in_max
                        do i2=in_mi,in_ma

c                        print *,' i1,i2,jtrans,ttrans,ia=',
c     $                          i1,i2,jtrans,ttrans,ia

                        tbdval=tbd(kff,kii)%Jtr(jtrans/2)
     $                       %fi(i1)%mT(ttrans/2)%in(i2)

                           if (abs(tbdval)>1.d-9) then

c                              print *,' ia,kff,kii,i1,i2=',
c     $                             ia,kff,kii,i1,i2

                              write(iunitvout,'(i4,1x,i6,1x,e16.9)')
     $                             i1,i2,tbdval/cleb/clebt*2.d0 ! for a<->b exchange

c                              print *,' before ia=',ia
                              ia=ia+1
c                              print *,' after ia=',ia
c                              print *,' ia,kff,kii,i1,i2,xxx=',
c     $                             ia,kff,kii,i1,i2,tbdval


                           endif
                        end do
                     end do

c                     print *, ' ia=',ia

                  end do
               end do
            end do
         end do
      endif

      contains

      subroutine one_nucleon(a_1)
      implicit none
      integer,intent(IN) :: a_1
      integer :: sp_st_1

      sp_st_1=iobsind(a_1)
      do kii=ki,ki+nki-1
         do kff=kf,kf+nkf-1
            ad_cl(sp_st_1,kff,kii)=ad_cl(sp_st_1,kff,kii)+
     +           bmpi(i1,kii)*bmpf(i2,kff)*real(phase,kind(0.d0))
         end do
      end do   

      end subroutine one_nucleon

      subroutine one_nucleon_two_body(c_1,a_2,a_1)
      implicit none
      integer,intent(IN) :: c_1,a_2,a_1
      integer :: sp_st_c1,sp_st_a1,sp_st_a2,sp_st_12,kii,kff,jtrans,
     $     j12,it12,ipht12,tb_st,Ttr_mi,Ttr_ma,ttrans
      real(kind(0.d0)) :: clebd,cleb1,cleb2,clebt1,clebt2

      sp_st_c1=iobsind(c_1)
      sp_st_a1=iobsind(a_1)
      sp_st_a2=iobsind(a_2)
c      print *,' sp_st_c1,sp_st_a1,sp_st_a2=',sp_st_c1,sp_st_a1,sp_st_a2
      if (sp_st_a1<=sp_st_a2) then
         sp_st_12=sp_st_a1+sp_st_a2*(sp_st_a2-1)/2
      else
         sp_st_12=sp_st_a2+sp_st_a1*(sp_st_a1-1)/2
      endif

c      print *,' sp_st_12=',sp_st_12

      do kii=ki,ki+nki-1
         do kff=kf,kf+nkf-1

c            print *,' kii,kff=',kii,kff

            Ttr_mi=max(abs(it2i(kii)-it2f(kff)),1,
     $           abs(mt2_spi(a_1)+mt2_spi(a_2)-mt2_spf(c_1)))
            Ttr_ma=min(it2i(kii)+it2f(kff),3)

c            if (abs(it2i(kii)-it2f(kff))>1.or.it2i(kii)+it2f(kff)<1)
c     $           cycle
            do jtrans=abs(jt2i(kii)-jt2f(kff)),jt2i(kii)+jt2f(kff),2

c               print *,' jtrans=',jtrans

               do j12=max(iabs(j2_spi(a_1)-j2_spi(a_2)),
     +              iabs(m2_spi(a_1)+m2_spi(a_2)),
     $              abs(j2_spf(c_1)-jtrans)),
     $              min(j2_spi(a_1)+j2_spi(a_2),j2_spf(c_1)+jtrans),2

c                  print *,' j12=',j12

                  cleb1=clebd(j2_spi(a_1),m2_spi(a_1),j2_spi(a_2),
     +                 m2_spi(a_2),j12,m2_spi(a_1)+m2_spi(a_2))
                  cleb2=clebd(j12,m2_spi(a_1)+m2_spi(a_2),j2_spf(c_1),
     $                 -m2_spf(c_1),jtrans,
     $                 m2_spi(a_1)+m2_spi(a_2)-m2_spf(c_1))

c                  print *,' cleb1,cleb2=',cleb1,cleb2

                  do ttrans=Ttr_mi,Ttr_ma,2

                     do it12=abs(mt2_spi(a_1)+mt2_spi(a_2)),2,2
                        if (sp_st_a1==sp_st_a2) then
                           if ((-1)**((j12+it12)/2)/=-1) cycle
                        endif

c                     print *,' it12=',it12

                        clebt1=clebd(1,mt2_spi(a_1),1,mt2_spi(a_2),it12,
     $                       mt2_spi(a_1)+mt2_spi(a_2))
                        clebt2=clebd(it12,mt2_spi(a_1)+mt2_spi(a_2),
     $                       1,-mt2_spf(c_1),ttrans,
     $                       mt2_spi(a_1)+mt2_spi(a_2)-mt2_spf(c_1))

c                     print *,' clebt1,clebt2=',clebt1,clebt2

                        if (sp_st_a1<=sp_st_a2) then
                           ipht12=1
                        else
                           ipht12=(-1)**((j12-j2_spi(a_1)
     +                          -j2_spi(a_2)+iT12)/2)
                        endif

                        tb_st=itbind(sp_st_12,j12/2,it12/2)                  

c                     print *,' tb_st=',tb_st

                        tbd(kff,kii)%Jtr(jtrans/2)%fi(sp_st_c1)
     $                       %mT(ttrans/2)%in(tb_st)=
     $                       tbd(kff,kii)%Jtr(jtrans/2)%fi(sp_st_c1)
     $                       %mT(ttrans/2)%in(tb_st)
     $                       +cleb1*cleb2*clebt1*clebt2
     $                       *bmpi(i1,kii)*bmpf(i2,kff)
     $                       *real(ipht12*iphase*((-1)**((j2_spf(c_1)
     $                       +m2_spf(c_1)+1+mt2_spf(c_1))/2)),
     $                       kind(0.d0))
                     end do
                  end do
               end do
            end do
         end do
      end do
      end subroutine one_nucleon_two_body

      subroutine two_nucleon(a_1,a_2)
      implicit none
      integer,intent(IN) :: a_1,a_2
      integer :: sp_st_1,sp_st_2,sp_st_12,tb_st,delta_mj,delta_mt
      integer :: iT12,J12
      real(kind=kind(0.d0)) :: pht12

      sp_st_1=iobsind(a_1)
      sp_st_2=iobsind(a_2)
      if (sp_st_1<=sp_st_2) then
         sp_st_12=sp_st_1+sp_st_2*(sp_st_2-1)/2
      else
         sp_st_12=sp_st_2+sp_st_1*(sp_st_1-1)/2
      endif
      
      delta_mj=m2_spi(a_1)+m2_spi(a_2)
      delta_mt=mt2_spi(a_1)+mt2_spi(a_2)

      do kii=ki,ki+nki-1
         do kff=kf,kf+nkf-1
            do iT12=abs(delta_mt)/2,1
               do J12=max(jtotal2min,abs(delta_mj)/2,
     +              abs(jt2i(kii)-jt2f(kff))/2,
     +              abs(j2_spi(a_1)-j2_spi(a_2))/2),
     +              min(jtotal2,(jt2i(kii)+jt2f(kff))/2,
     +              (j2_spi(a_1)+j2_spi(a_2))/2)
                  
                  if (sp_st_1==sp_st_2.and.
     +                 (-1)**(J12+iT12)/=-1) cycle
                  if (sp_st_1<=sp_st_2) then
                     pht12=1.d0
                  else
                     pht12=real((-1)**(J12-(j2_spi(a_1)
     +                    +j2_spi(a_2))/2+iT12),kind(0.d0))
                  endif

                  tb_st=itbind(sp_st_12,J12,iT12)                  
                  
                  ad_cl(tb_st,kff,kii)=ad_cl(tb_st,kff,kii)
     +                 +clebd(j2_spi(a_1),m2_spi(a_1),j2_spi(a_2),
     +                 m2_spi(a_2),2*J12,delta_mj)
     +                 *clebd(1,mt2_spi(a_1),1,mt2_spi(a_2),2*iT12,
     +                 delta_mt)
     +                 *bmpi(i1,kii)*bmpf(i2,kff)
     +                 *real(phase,kind(0.d0))*pht12
                  
               end do
            end do
         end do
      end do

      end subroutine two_nucleon

      subroutine three_nucleon(a_1,a_2,a_3,perm)
      implicit none
      integer,intent(IN) :: a_1,a_2,a_3,perm
      integer :: sp_st_1,sp_st_2,sp_st_12,tb_st
      integer :: iT12,J12,sp_st_3,thr_st
      real(kind=kind(0.d0)) :: pht12
      real(kind=kind(0.d0)) :: cleb1,cleb2,clebt1,clebt2

      sp_st_1=iobsind(a_1)
      sp_st_2=iobsind(a_2)
      sp_st_3=iobsind(a_3)
      if (sp_st_1<=sp_st_2) then
         sp_st_12=sp_st_1+sp_st_2*(sp_st_2-1)/2
      else
         sp_st_12=sp_st_2+sp_st_1*(sp_st_1-1)/2
      endif

      do J12=max(abs(j2_spi(a_1)-j2_spi(a_2)),
     +     abs(m2_spi(a_1)+m2_spi(a_2)),
     +     abs(J123-j2_spi(a_3))),
     +     min((j2_spi(a_1)+j2_spi(a_2)),
     +     (J123+j2_spi(a_3))),2
         cleb1=cleb_3j(1,j2_spi(a_1),m2_spi(a_1),
     +        j2_spi(a_2),m2_spi(a_2),J12)
c         cleb1=clebd(j2_spi(a_1),m2_spi(a_1),
c     +        j2_spi(a_2),m2_spi(a_2),J12,
c     +        m2_spi(a_1)+m2_spi(a_2))
         if (cleb1==0.d0) cycle
         cleb2=cleb_3j(2,J12,m2_spi(a_1)+m2_spi(a_2),
     +        j2_spi(a_3),m2_spi(a_3),J123)
c         cleb2=clebd(J12,m2_spi(a_1)+m2_spi(a_2),
c     +        j2_spi(a_3),m2_spi(a_3),J123,delta_mj)
         if (cleb2==0.d0) cycle
         
         do iT12=max(abs(mt2_spi(a_1)+mt2_spi(a_2)),
     +        abs(iT123-1)),2,2
            if (sp_st_1==sp_st_2.and.
     +           (-1)**((J12+iT12)/2)/=-1) cycle
            
            clebt1=cleb_3j(3,1,mt2_spi(a_1),1,
     +           mt2_spi(a_2),iT12)
c            clebt1=clebd(1,mt2_spi(a_1),1,
c     +           mt2_spi(a_2),iT12,mt2_spi(a_1)
c     +           +mt2_spi(a_2))
            if (clebt1==0.d0) cycle
            clebt2=cleb_3j(4,iT12,mt2_spi(a_1)+mt2_spi(a_2),
     +           1,mt2_spi(a_3),iT123)
c            clebt2=clebd(iT12,mt2_spi(a_1)+mt2_spi(a_2),
c     +           1,mt2_spi(a_3),iT123,delta_mt)
            if (clebt2==0.d0) cycle
            
            if (sp_st_1<=sp_st_2) then
               pht12=1.d0
            else
               pht12=real((-1)**((J12-j2_spi(a_1)
     +              -j2_spi(a_2)+iT12)/2),kind(0.d0))
            endif
            
            tb_st=itbind(sp_st_12,J12/2,iT12/2)
            thr_st=ad3ind(tb_st,sp_st_3,J123/2,iT123/2)
            
            ad_cl(thr_st,kff,kii)=ad_cl(thr_st,kff,kii)
     +           +cleb1*cleb2*clebt1*clebt2
     +           *bmpi(i1,kii)*bmpf(i2,kff)
     +           *real(phase*perm,kind(0.d0))*pht12

         end do
      end do

      end subroutine three_nucleon

      subroutine four_nucleon(a_1,a_2,a_3,a_4,perm)
      implicit none
      integer,intent(IN) :: a_1,a_2,a_3,a_4,perm
      integer :: sp_st_1,sp_st_2,sp_st_12,tb_st
      integer :: iT12,J12,sp_st_3,thr_st,sp_st_4,
     +     fou_st,iT123,J123
      real(kind=kind(0.d0)) :: pht12
      real(kind=kind(0.d0)) :: cleb1,cleb2,cleb3,clebt1,clebt2,clebt3

      sp_st_1=iobsind(a_1)
      sp_st_2=iobsind(a_2)
      sp_st_3=iobsind(a_3)
      sp_st_4=iobsind(a_4)
      if (sp_st_1<=sp_st_2) then
         sp_st_12=sp_st_1+sp_st_2*(sp_st_2-1)/2
      else
         sp_st_12=sp_st_2+sp_st_1*(sp_st_1-1)/2
      endif

      do iT123=max(abs(mt2_spi(a_1) !     (123)4 
     +     +mt2_spi(a_2)+mt2_spi(a_3)),
     +     abs(iT1234-1)),min(iT1234+1,3),2
         clebt3=cleb_3j(6,iT123,mt2_spi(a_1)+mt2_spi(a_2)
     +        +mt2_spi(a_3),1,mt2_spi(a_4),iT1234)
c         clebt3=clebd(iT123,mt2_spi(a_1)+mt2_spi(a_2)
c     +        +mt2_spi(a_3),1,mt2_spi(a_4),iT1234,delta_mt)
         if (clebt3==0.d0) cycle
         do J123=max(abs(m2_spi(a_1)+m2_spi(a_2)
     +        +m2_spi(a_3)),abs(J1234-j2_spi(a_4))),
     +        J1234+j2_spi(a_4),2
            cleb3=cleb_3j(3,J123,m2_spi(a_1)+m2_spi(a_2)
     +           +m2_spi(a_3),j2_spi(a_4),m2_spi(a_4),
     +           J1234)
c            cleb3=clebd(J123,m2_spi(a_1)+m2_spi(a_2)
c     +           +m2_spi(a_3),j2_spi(a_4),m2_spi(a_4),
c     +           J1234,delta_mj)
            if (cleb3==0.d0) cycle
                        
            do iT12=max(abs(mt2_spi(a_1) !    (12)34
     +           +mt2_spi(a_2)),abs(iT123-1)),2,2
               clebt1=cleb_3j(4,1,mt2_spi(a_1),1,mt2_spi(a_2),iT12)
c               clebt1=clebd(1,mt2_spi(a_1),1,mt2_spi(a_2),
c     +              iT12,mt2_spi(a_1)+mt2_spi(a_2))
               if (clebt1==0.d0) cycle
               clebt2=cleb_3j(5,iT12,mt2_spi(a_1)+mt2_spi(a_2),
     +              1,mt2_spi(a_3),iT123)
c               clebt2=clebd(iT12,mt2_spi(a_1)+mt2_spi(a_2),
c     +              1,mt2_spi(a_3),iT123,mt2_spi(a_1)
c     +              +mt2_spi(a_2)+mt2_spi(a_3))
               if (clebt2==0.d0) cycle
               do J12=max(abs(j2_spi(a_1)-j2_spi(a_2)),
     +              abs(m2_spi(a_1)+m2_spi(a_2)),
     +              abs(J123-j2_spi(a_3))),min(j2_spi(a_1)
     +              +j2_spi(a_2),J123+j2_spi(a_3)),2
                  
                  if (sp_st_1==sp_st_2.and.
     +                 (-1)**((J12+iT12)/2)/=-1) cycle
                  
                  cleb1=cleb_3j(1,j2_spi(a_1),m2_spi(a_1),
     +                 j2_spi(a_2),m2_spi(a_2),J12)
c                  cleb1=clebd(j2_spi(a_1),m2_spi(a_1),
c     +                 j2_spi(a_2),m2_spi(a_2),J12,
c     +                 m2_spi(a_1)+m2_spi(a_2))
                  if (cleb1==0.d0) cycle
                  cleb2=cleb_3j(2,J12,m2_spi(a_1)+m2_spi(a_2),
     +                 j2_spi(a_3),m2_spi(a_3),J123)
c                  cleb2=clebd(J12,m2_spi(a_1)+m2_spi(a_2),
c     +                 j2_spi(a_3),m2_spi(a_3),J123,
c     +                 m2_spi(a_1)+m2_spi(a_2)+m2_spi(a_3))
                  if (cleb2==0.d0) cycle
                  
                  if (sp_st_1<=sp_st_2) then
                     pht12=1.d0
                  else
                     pht12=real((-1)**((J12-j2_spi(a_1)
     +                    -j2_spi(a_2)+iT12)/2),kind(0.d0))
                  endif
                  tb_st=itbind(sp_st_12,J12/2,iT12/2)
                  thr_st=ad3ind(tb_st,sp_st_3,J123/2,
     +                 iT123/2)
                  fou_st=ad4ind(thr_st,sp_st_4,
     +                 J1234/2,iT1234/2)
                  
c                  if (fou_st==212429) then
c                     print *,' fou_st,thr_st,tb_st=',fou_st,thr_st,tb_st
c                     print *,' sp_st_1,sp_st_2,sp_st_3,sp_st_4=',
c     +                    sp_st_1,sp_st_2,sp_st_3,sp_st_4
c                     print *,' J12,iT12,J123,iT123,J1234,iT1234=',
c     +                    J12,iT12,J123,iT123,J1234,iT1234
c                     print *,' cleb1,cleb2,cleb3=',cleb1,cleb2,cleb3
c                     print *,' clebt1,clebt2,clebt3=',
c     +                    clebt1,clebt2,clebt3
c                     print *,' i1,i2=',i1,i2
c                     print *,' kii,kff=',kii,kff
c                     print *,' bmpi,bmpf=',bmpi(i1,kii),bmpf(i2,kff)
c                     print *,' phase,perm,pht12=',phase,perm,pht12
c                     print *,' before adding: ad_cl=',
c     +                    ad_cl(fou_st,kff,kii)
c                  endif

                  ad_cl(fou_st,kff,kii)=
     +                 ad_cl(fou_st,kff,kii)
     +                 +cleb1*cleb2*cleb3
     +                 *clebt1*clebt2*clebt3
     +                 *bmpi(i1,kii)*bmpf(i2,kff)
     +                 *real(phase*perm,kind(0.d0))*pht12

c                  if (fou_st==212429) then
c                     print *,' after adding: ad_cl=',
c     +                    ad_cl(fou_st,kff,kii)
c                  endif
                              
               end do
            end do
         end do
      end do

      end subroutine four_nucleon

      end

      integer function Nmin_HO(Z)
      implicit none
      integer,intent(IN) :: Z
      integer :: Nmin,N,Zrem
      Zrem=Z
      Nmin=0
      N=0
      do
         Nmin=Nmin+N*min((N+1)*(N+2),Zrem)
         Zrem=Zrem-(N+1)*(N+2)
         if (Zrem<=0) exit
         N=N+1
      end do
      Nmin_HO=Nmin
      end function Nmin_HO

      subroutine global_ist1_setup
      use paramdef
      use nodeinfo, only:iproc
      implicit none
      integer :: Nmax,n,l,j,ii,N_cluster_ma,kk,mj,mt
      N_cluster_ma=max(N_cluster_max,N_RGM_max+N_proj_max)
      allocate(nlj_st(0:N_cluster_ma/2))
      do n=0,N_cluster_ma/2
         allocate(nlj_st(n)%l(0:N_cluster_ma-2*n))
         do l=0,N_cluster_ma-2*n
            allocate(nlj_st(n)%l(l)%j((2*l-1)/2:(2*l+1)/2))
         end do
      end do
      ii=0
      do Nmax=0,N_cluster_ma
         do l=mod(Nmax,2),Nmax,2
            n=(Nmax-l)/2
            do j=abs(2*l-1),2*l+1,2
               ii=ii+1
            end do
         end do
      end do
      global_ist1_dim=ii
c*** MPI
      if (iproc==0) print *,' global_ist1_dim=',global_ist1_dim
c*** MPI
      allocate(global_ist1(3,global_ist1_dim))
      ii=0
      do Nmax=0,N_cluster_ma
         do l=mod(Nmax,2),Nmax,2
            do j=abs(2*l-1),2*l+1,2
               n=(Nmax-l)/2
               ii=ii+1
               global_ist1(1,ii)=n
               global_ist1(2,ii)=l
               global_ist1(3,ii)=j
               nlj_st(n)%l(l)%j(j/2)=ii
c               print *,' ii,nlj=',ii,n,l,j
            end do
         end do
      end do
      ii=0
      do kk=1,global_ist1_dim
         n=global_ist1(1,kk)
         l=global_ist1(2,kk)
         j=global_ist1(3,kk)
         do mj=-j,j,2
            do mt=1,-1,-2
               ii=ii+1
            end do
         end do
      end do
      global_ist1_sp_dim=ii
      allocate(global_ist1_sp(3,ii))
c*** MPI
      if (iproc==0) print *,' global_ist1_sp_dim=',global_ist1_sp_dim
c*** MPI
      ii=0
      do kk=1,global_ist1_dim
         n=global_ist1(1,kk)
         l=global_ist1(2,kk)
         j=global_ist1(3,kk)
         do mj=-j,j,2
            do mt=1,-1,-2
               ii=ii+1
               global_ist1_sp(1,ii)=kk
               global_ist1_sp(2,ii)=mj
               global_ist1_sp(3,ii)=mt
            end do
         end do
      end do
      allocate(nljmmt_st(global_ist1_dim,0:1))
      do kk=1,global_ist1_dim
         n=global_ist1(1,kk)
         l=global_ist1(2,kk)
         j=global_ist1(3,kk)
         allocate(nljmmt_st(kk,0)%mj_sp(-j:j))
         allocate(nljmmt_st(kk,1)%mj_sp(-j:j))
      end do
      do ii=1,global_ist1_sp_dim
         kk=global_ist1_sp(1,ii)
         mj=global_ist1_sp(2,ii)
         mt=global_ist1_sp(3,ii)
         n=global_ist1(1,kk)
         l=global_ist1(2,kk)
         j=global_ist1(3,kk)
         nljmmt_st(kk,(mt+1)/2)%mj_sp(mj)=ii
      end do
      end

      subroutine check_jlevel_consistency
      use paramdef
      use initst, only: n_spi,l_spi,j2_spi,m2_spi,mt2_spi,mxnwdi,
     +     naspsi
      use obdens
      use nodeinfo, only: iproc
      implicit none
      integer :: ii,n,lj,n1,l1,j1,m1,m,mt,jj,kk,l,j
      if (iproc==0) print *,' ist1dim,global_ist1_dim=',
     $     ist1dim,global_ist1_dim
      do ii=1,min(ist1dim,global_ist1_dim)
         n=global_ist1(1,ii)
         l=global_ist1(2,ii)
         j=global_ist1(3,ii)
         n1=ist1(1,ii)
         l1=ist1(2,ii)
         j1=ist1(3,ii)
         if (n/=n1.or.l/=l1.or.j/=j1) then
            print *,'***error iproc,n,l,j,n1,l1,j1:',
     $           iproc,n,l,j,n1,l1,j1
            stop
         endif
      end do
      allocate(global_ist1_sp_map(2*mxnwdi*nbit))
      do ii=1,max(naspsi,nasps_full)
         n1=n_spi(ii)
         l1=l_spi(ii)
         j1=j2_spi(ii)
         m1=m2_spi(ii)
         do jj=1,global_ist1_sp_dim
            kk=global_ist1_sp(1,jj)
            m=global_ist1_sp(2,jj)
            mt=global_ist1_sp(3,jj)
            n=global_ist1(1,kk)
            l=global_ist1(2,kk)
            j=global_ist1(3,kk)
            if (n==n1.and.l==l1.and.j==j1.and.m==m1) then
               if (mt==1) then
                  global_ist1_sp_map(ii)=jj
               elseif (mt==-1) then
                  global_ist1_sp_map(ii+mxnwdi*nbit)=jj
               else
                  print *,'***error in global_ist1_sp_map: mt=',mt
                  stop
               endif
            endif
         end do
      end do
      end

      subroutine ist2_setup
      use paramdef
      use kernels, only: ist2_SDkern,ist2_SDkern_dim,I_ab_max
      use interaction, only: N_sp_int_1max,N_sp_int_2max,N_sp_int_12max,
     $     tbdst_SDkern
      use v3b, only: N1_max,N12_max     
      use intrface, only: kf,nkf,ki,nki
      use initst, only: nhomi,nhom12i,jt2i,it2i,iparityi
      use finast, only: jt2f,it2f,iparityf
      use nodeinfo, only:iproc
      implicit none
      integer :: ia,ib,na,la,ja,nb,lb,jb,I_ab,t_ab,ii,proj_tz,
     $     ia_max,ib_min,ib_max,total,pi,Iabmax,Iabmin,tabmin,tabmax,
     $     iparmin,iparmax

      if (num_prot_proj+num_neutr_proj==2) then
         proj_tz=num_prot_proj-num_neutr_proj
         Iabmax=(jt2i(ki)+jt2f(kf))/2
         Iabmin=abs(jt2i(ki)-jt2f(kf))/2
         tabmax=(it2i(ki)+it2f(kf))/2
         tabmin=abs(it2i(ki)-it2f(kf))/2
         do ia=ki,ki+nki-1
            do ib=kf,kf+nkf-1
               Iabmax=max(Iabmax,(jt2i(ia)+jt2f(ib))/2)
               Iabmin=min(Iabmin,abs(jt2i(ia)-jt2f(ib))/2)
               tabmax=max(tabmax,(it2i(ia)+it2f(ib))/2)
               tabmin=min(tabmin,abs(it2i(ia)-it2f(ib))/2)
            end do
         end do
         iparmin=mod(iparityi+iparityf,2)
         iparmax=iparmin
         if (iproc==0) then
            print *,' proj_tz=',proj_tz
            print *,' Iabmin,Iabmax=',Iabmin,Iabmax
            print *,' tabmin,tabmax=',tabmin,tabmax
            print *,' iparity=',iparmin
         endif
      else
         proj_tz=0
         Iabmin=0
         tabmin=0
         Iabmax=2**30
         tabmax=2**30
         iparmin=0
         iparmax=1
      endif
      I_ab_max=0
      ii=0
      do ia=1,global_ist1_dim
         na=global_ist1(1,ia)
         la=global_ist1(2,ia)
         ja=global_ist1(3,ia)
         if (2*na+la>max(nhomi,N1_max,N_sp_int_2max)) exit
         do ib=1,global_ist1_dim
            nb=global_ist1(1,ib)
            lb=global_ist1(2,ib)
            jb=global_ist1(3,ib)
            if (2*nb+lb>max(nhomi,N1_max,N_sp_int_2max)) exit
            if (2*na+la+2*nb+lb>
     $           max(nhom12i,N12_max,nhomi+N_sp_int_2max)) exit
            pi=mod(la+lb,2)
            if (pi/=iparmin.and.pi/=iparmax) cycle
            do I_ab=max(abs(ja-jb)/2,Iabmin),min((ja+jb)/2,Iabmax)
               do t_ab=max(abs(proj_tz)/2,tabmin),min(1,tabmax)
                  if (ia==ib.and.mod(I_ab+t_ab,2)==0) cycle
                  ii=ii+1
                  if (I_ab>I_ab_max) I_ab_max=I_ab
               end do
            end do
         end do
      end do
      ist2_SDkern_dim=ii
      if (iproc==0) then
         print *,' ist2_SDkern_dim=',ist2_SDkern_dim
         print *,' I_ab_max=',I_ab_max
      endif

      allocate(tbdst_SDkern(Iabmin:I_ab_max,iparmin:iparmax))
      ia_max=0
      do ia=1,global_ist1_dim
         na=global_ist1(1,ia)
         la=global_ist1(2,ia)
         ja=global_ist1(3,ia)
         if (2*na+la>max(nhomi,N1_max,N_sp_int_2max)) exit
         ia_max=ia_max+1
      end do
c*** MPI
      if (iproc==0) then
         print *,' ia_max=',ia_max
      endif
c*** MPI
      do I_ab=Iabmin,I_ab_max
         do pi=iparmin,iparmax
            do t_ab=max(abs(proj_tz)/2,tabmin),min(1,tabmax)
               total=0
               tbdst_SDkern(I_ab,pi)%T(t_ab)%sp1max=ia_max
               allocate(tbdst_SDkern(I_ab,pi)%T(t_ab)%sp1(ia_max))
               do ia=1,ia_max
                  na=global_ist1(1,ia)
                  la=global_ist1(2,ia)
                  ja=global_ist1(3,ia)
                  if (2*na+la>max(nhomi,N1_max,N_sp_int_2max)) exit
                  ib_min=global_ist1_dim
                  ib_max=1
                  do ib=1,global_ist1_dim
                     nb=global_ist1(1,ib)
                     lb=global_ist1(2,ib)
                     jb=global_ist1(3,ib)
                     if (2*nb+lb>max(nhomi,N1_max,N_sp_int_2max)) exit
                     if (2*na+la+2*nb+lb>
     $                    max(nhom12i,N12_max,nhomi+N_sp_int_2max)) exit
                     if (mod(la+lb,2)/=pi) cycle
                     if (abs(ja-jb)>2*I_ab.or.ja+jb<2*I_ab) cycle
                     if (ia==ib.and.mod(I_ab+t_ab,2)/=1) cycle
                     ib_min=min(ib_min,ib)
                     ib_max=max(ib_max,ib)
                     total=total+1
                  end do
                  tbdst_SDkern(I_ab,pi)%T(t_ab)%sp1(ia)%sp2max=ib_max
                  tbdst_SDkern(I_ab,pi)%T(t_ab)%sp1(ia)%sp2min=ib_min
                  allocate(tbdst_SDkern(I_ab,pi)%T(t_ab)%sp1(ia)%sp2
     $                 (ib_min:ib_max))
                  tbdst_SDkern(I_ab,pi)%T(t_ab)%sp1(ia)%sp2(:)=-1
               end do
c               print *,' J,pi,T,tbdst%total=',J,pi,T,total
               tbdst_SDkern(I_ab,pi)%T(t_ab)%total=total
c               allocate(tbdst_SDkern(J,pi)%T(T)%ist(2,total))
            end do
         end do
      end do

      allocate(ist2_SDkern(4,ist2_SDkern_dim))
      ii=0
      do I_ab=Iabmin,I_ab_max
         do pi=iparmin,iparmax
            do t_ab=max(abs(proj_tz)/2,tabmin),min(1,tabmax)
               do ia=1,global_ist1_dim
                  na=global_ist1(1,ia)
                  la=global_ist1(2,ia)
                  ja=global_ist1(3,ia)
                  if (2*na+la>max(nhomi,N1_max,N_sp_int_2max)) exit
                  do ib=1,global_ist1_dim
                     nb=global_ist1(1,ib)
                     lb=global_ist1(2,ib)
                     jb=global_ist1(3,ib)
                     if (2*nb+lb>max(nhomi,N1_max,N_sp_int_2max)) exit
                     if (2*na+la+2*nb+lb>
     $                    max(nhom12i,N12_max,nhomi+N_sp_int_2max)) exit
                     if (mod(la+lb,2)/=pi) cycle
                     if (abs(ja-jb)>2*I_ab.or.ja+jb<2*I_ab) cycle
                     if (ia==ib.and.mod(I_ab+t_ab,2)==0) cycle
                     ii=ii+1
                     ist2_SDkern(1,ii)=ia
                     ist2_SDkern(2,ii)=ib
                     ist2_SDkern(3,ii)=I_ab
                     ist2_SDkern(4,ii)=t_ab
                     tbdst_SDkern(I_ab,pi)%T(t_ab)%sp1(ia)%sp2(ib)=ii
                  end do
               end do
            end do
         end do
      end do
      if (iproc==0) then
         print *,' ii=',ii
      endif
      end

      subroutine ist3_setup
      use paramdef
      use kernels, only: ist2_SDkern,ist2_SDkern_dim,I_ab_max,
     $     ist3_SDkern,ist3_SDkern_dim,I_abc_max
      use interaction, only: N_sp_int_1max,N_sp_int_2max,N_sp_int_12max,
     $     tbdst_SDkern,thrbdst_SDkern
      use v3b, only: N1_max,N12_max,N123_max      
      use intrface, only: kf,nkf,ki,nki
      use initst, only: nhomi,nhom12i,nhom123i,jt2i,it2i,iparityi
      use finast, only: jt2f,it2f,iparityf
      use nodeinfo, only:iproc
      implicit none
      integer :: ia,ib,na,la,ja,nb,lb,jb,I_ab,t_ab,ii,proj_tz,
     $     ia_max,ib_min,ib_max,total,pi,ic,nc,lc,jc,I_abc,t_abc,ic_max,
     $     ic_min,ij,Iabcmax,Iabcmin,tabcmin,tabcmax,iparmin,iparmax

      if (num_prot_proj+num_neutr_proj==3) then
         proj_tz=num_prot_proj-num_neutr_proj
         Iabcmax=jt2i(ki)+jt2f(kf)
         Iabcmin=abs(jt2i(ki)-jt2f(kf))
         tabcmax=it2i(ki)+it2f(kf)
         tabcmin=abs(it2i(ki)-it2f(kf))
         do ia=ki,ki+nki-1
            do ib=kf,kf+nkf-1
               Iabcmax=max(Iabcmax,jt2i(ia)+jt2f(ib))
               Iabcmin=min(Iabcmin,abs(jt2i(ia)-jt2f(ib)))
               tabcmax=max(tabcmax,it2i(ia)+it2f(ib))
               tabcmin=min(tabcmin,abs(it2i(ia)-it2f(ib)))
            end do
         end do
         iparmin=mod(iparityi+iparityf,2)
         iparmax=iparmin
         if (iproc==0) then
            print *,' proj_tz=',proj_tz
            print *,' Iabcmin,Iabcmax=',Iabcmin,Iabcmax
            print *,' tabcmin,tabcmax=',tabcmin,tabcmax
            print *,' iparity=',iparmin
         endif
      else
         proj_tz=1
         Iabcmin=1
         tabcmin=1
         Iabcmax=2**30
         tabcmax=2**30
         iparmin=0
         iparmax=1
      endif
      I_abc_max=1
      ii=0
      do ij=1,ist2_SDkern_dim
         ia=ist2_SDkern(1,ij)
         ib=ist2_SDkern(2,ij)
         I_ab=ist2_SDkern(3,ij)
         t_ab=ist2_SDkern(4,ij)
         na=global_ist1(1,ia)
         la=global_ist1(2,ia)
         ja=global_ist1(3,ia)
         nb=global_ist1(1,ib)
         lb=global_ist1(2,ib)
         jb=global_ist1(3,ib)
         do ic=1,global_ist1_dim
            nc=global_ist1(1,ic)
            lc=global_ist1(2,ic)
            jc=global_ist1(3,ic)
            if (2*nc+lc>max(nhomi,N1_max,N_sp_int_2max)) exit
            if (2*na+la+2*nb+lb+2*nc+lc>
     $           max(nhom123i,N123_max,nhom12i+N_sp_int_2max)) exit   ! OK for NN-only
c     $           max(nhom123i,N123_max,nhomi+N_sp_int_12max)) exit
            pi=mod(la+lb+lc,2)
            if (pi/=iparmin.and.pi/=iparmax) cycle
            do I_abc=max(abs(jc-2*I_ab),Iabcmin),
     $           min(jc+2*I_ab,Iabcmax),2
               do t_abc=max(abs(2*t_ab-1),abs(proj_tz),tabcmin),
     $              min(2*t_ab+1,tabcmax),2
                  ii=ii+1
                  if (I_abc>I_abc_max) I_abc_max=I_abc
               end do
            end do
         end do
      end do

      ist3_SDkern_dim=ii
      if (iproc==0) then
         print *,' ist3_SDkern_dim=',ist3_SDkern_dim
         print *,' I_abc_max=',I_abc_max
      endif
      allocate(thrbdst_SDkern(Iabcmin/2:I_abc_max/2,iparmin:iparmax)) !parity %T(abs(proj_tz)/2:3/2)
      ii=0
      do I_abc=Iabcmin,I_abc_max,2
         do pi=iparmin,iparmax
            do t_abc=max(tabcmin,abs(proj_tz)),min(tabcmax,3),2
               allocate(thrbdst_SDkern(I_abc/2,pi)%T(t_abc/2)
     $              %sp1(ist2_SDkern_dim))
               do ij=1,ist2_SDkern_dim
                  ia=ist2_SDkern(1,ij)
                  ib=ist2_SDkern(2,ij)
                  I_ab=ist2_SDkern(3,ij)
                  t_ab=ist2_SDkern(4,ij)
                  na=global_ist1(1,ia)
                  la=global_ist1(2,ia)
                  ja=global_ist1(3,ia)
                  nb=global_ist1(1,ib)
                  lb=global_ist1(2,ib)
                  jb=global_ist1(3,ib)
                  ic_min=global_ist1_dim
                  ic_max=0
                  do ic=1,global_ist1_dim
                     nc=global_ist1(1,ic)
                     lc=global_ist1(2,ic)
                     jc=global_ist1(3,ic)
                     if (2*nc+lc>max(nhomi,N1_max,N_sp_int_2max)) exit
                     if (2*na+la+2*nb+lb+2*nc+lc>
     $                max(nhom123i,N123_max,nhom12i+N_sp_int_2max)) exit ! OK for NN-only
c     $               max(nhom123i,N123_max,nhomi+N_sp_int_12max)) exit
                     if (mod(la+lb+lc,2)/=pi) cycle
                     if (abs(jc-2*I_ab)>I_abc.or.jc+2*I_ab<I_abc) cycle
                     if (abs(1-2*t_ab)>t_abc.or.1+2*t_ab<t_abc) cycle
                     if (ic<ic_min) ic_min=ic
                     if (ic>ic_max) ic_max=ic
                  end do
                  thrbdst_SDkern(I_abc/2,pi)%T(t_abc/2)
     $                 %sp1(ij)%sp2max=ic_max
                  thrbdst_SDkern(I_abc/2,pi)%T(t_abc/2)
     $                 %sp1(ij)%sp2min=ic_min
                  allocate(thrbdst_SDkern(I_abc/2,pi)%T(t_abc/2)
     $                 %sp1(ij)%sp2(ic_min:ic_max))
                  thrbdst_SDkern(I_abc/2,pi)%T(t_abc/2)
     $                 %sp1(ij)%sp2(:)=-1
c                  do ic=1,global_ist1_dim
                  do ic=ic_min,ic_max
                     nc=global_ist1(1,ic)
                     lc=global_ist1(2,ic)
                     jc=global_ist1(3,ic)
                     if (2*nc+lc>max(nhomi,N1_max,N_sp_int_2max)) exit
                     if (2*na+la+2*nb+lb+2*nc+lc>
     $                max(nhom123i,N123_max,nhom12i+N_sp_int_2max)) exit ! OK for NN-only
c     $               max(nhom123i,N123_max,nhomi+N_sp_int_12max)) exit
                     if (mod(la+lb+lc,2)/=pi) cycle
                     if (abs(jc-2*I_ab)>I_abc.or.jc+2*I_ab<I_abc) cycle
                     if (abs(1-2*t_ab)>t_abc.or.1+2*t_ab<t_abc) cycle
                     ii=ii+1
                  end do
               end do
            end do
         end do
      end do
      if (iproc==0) print *,' ii=',ii
      allocate(ist3_SDkern(4,ist3_SDkern_dim))
      ii=0
      do I_abc=Iabcmin,I_abc_max,2
         do pi=iparmin,iparmax
            do t_abc=max(tabcmin,abs(proj_tz)),min(tabcmax,3),2
               do ij=1,ist2_SDkern_dim
                  ia=ist2_SDkern(1,ij)
                  ib=ist2_SDkern(2,ij)
                  I_ab=ist2_SDkern(3,ij)
                  t_ab=ist2_SDkern(4,ij)
                  na=global_ist1(1,ia)
                  la=global_ist1(2,ia)
                  ja=global_ist1(3,ia)
                  nb=global_ist1(1,ib)
                  lb=global_ist1(2,ib)
                  jb=global_ist1(3,ib)
                  do ic=1,global_ist1_dim
                     nc=global_ist1(1,ic)
                     lc=global_ist1(2,ic)
                     jc=global_ist1(3,ic)
                     if (2*nc+lc>max(nhomi,N1_max,N_sp_int_2max)) exit
                     if (2*na+la+2*nb+lb+2*nc+lc>
     $                max(nhom123i,N123_max,nhom12i+N_sp_int_2max)) exit ! OK for NN-only
c     $           max(nhom123i,N123_max,nhomi+N_sp_int_12max)) exit
                     if (mod(la+lb+lc,2)/=pi) cycle
                     if (abs(jc-2*I_ab)>I_abc.or.jc+2*I_ab<I_abc) cycle
                     if (abs(1-2*t_ab)>t_abc.or.1+2*t_ab<t_abc) cycle
                     ii=ii+1
                     ist3_SDkern(1,ii)=ij
                     ist3_SDkern(2,ii)=ic
                     ist3_SDkern(3,ii)=I_abc
                     ist3_SDkern(4,ii)=t_abc
                     thrbdst_SDkern(I_abc/2,pi)%T(t_abc/2)
     $                    %sp1(ij)%sp2(ic)=ii
                  end do
               end do
            end do
         end do
      end do
      if (iproc==0) print *,' ii=',ii
      end
      
      subroutine interaction_read
      use paramdef
      use interaction
      use nodeinfo, only: iproc
      implicit none
      integer :: i_file,j,pi,t,tz,s,n,l,np,lp,lmin,lma
      real(kind(0.d0)) :: v

      allocate(V_rel(0:jrelm,0:1))
      do j=0,jrelm
         do pi=0,1
            do t=0,1
               s=1-mod(pi+t,2)
               if ((s==0.or.j==0).and.mod(j+s+pi,2)/=0) cycle
               allocate(V_rel(j,pi)%t(t)%tz(-t:t))
               do tz=-t,t
                  if (mod(abs(j-s),2)==pi) then
                     lmin=abs(j-s)
                     lma=j+s
                  else
                     lmin=abs(j-s)+1
                     lma=lmin
                  endif
                  V_rel(j,pi)%t(t)%tz(tz)%lmin=lmin
                  V_rel(j,pi)%t(t)%tz(tz)%lmax=lma
                  allocate(V_rel(j,pi)%t(t)%tz(tz)%ll(lmin:lma,
     $                 lmin:lma))
                  do l=lmin,lma,2
                     do lp=lmin,lma,2
                        allocate(V_rel(j,pi)%t(t)%tz(tz)%ll(l,lp)
     $                       %nn(0:(N_cluster_max-l)/2,
     $                       0:(N_cluster_max-lp)/2))
                        V_rel(j,pi)%t(t)%tz(tz)%ll(l,lp)%nn=0.d0
                     end do
                  end do
               end do
            end do
         end do
      end do

      do i_file=1,num_of_interaction_files
         select case(i_file)
         case(1)
            tz=0
         case(2)
            tz=1
         case(3)
            tz=-1
         case default
            print *,'*** error: i_file=',i_file
            stop
         end select
         open(37,file=trim(interaction_file(i_file)),form='unformatted',
     $        status='old',action='read') 
!SQ         open(37,file=trim(interaction_file(i_file)),form='formatted',
!SQ     $        status='old',action='read')
c*** MPI
        if (iproc==0) then
c*** MPI
         write(2,"(' Interaction file '(a))")
     $        trim(interaction_file(i_file))
c*** MPI
         endif
c*** MPI
         do
            read(37,end=1000,err=1000) s,j,t,n,l,np,lp,v 
!SQ            read(37,*,end=1000,err=1000) s,j,t,n,l,np,lp,v
c            print *,'b: s,j,t,tz,  n,l,np,lp, v=',s,j,t,tz,n,l,np,lp,v
            if (mod(s+l+t,2)==1.and.mod(l+lp,2)==0.and.(s==0.or.s==1)
     $           .and.(t==0.or.t==1).and.j<=jrelm.and.t>=abs(tz).and.
     $           n<=(N_cluster_max-l)/2.and.np<=(N_cluster_max-lp)/2)
cccc     $           n<=(N_cluster_min-l)/2.and.np<=(N_cluster_min-lp)/2) ! test
     $           then

c            if (t==1) v=0.d0
c               if (n==np.and.l==lp) then      !test
c                  v=1.d0                      !test
c               else                           !test
c                  v=0.d0                      !test
c               endif                          !test
c            if (n/=np.or.l/=lp) v=0.d0
c            if (.not.(n==0.and.l==0.and.np==0.and.lp==2)) v=0.d0
c            if ((n==0.and.l==2.and.np==1.and.lp==0)) v=0.d0
c            if ((n==1.and.l==0.and.np==0.and.lp==2)) v=0.d0

               pi=mod(l,2)
               V_rel(j,pi)%t(t)%tz(tz)%ll(l,lp)%nn(n,np)=v
               V_rel(j,pi)%t(t)%tz(tz)%ll(lp,l)%nn(np,n)=v
               if (num_of_interaction_files==1.and.t==1) then
                  V_rel(j,pi)%t(t)%tz(1)%ll(l,lp)%nn(n,np)=v
                  V_rel(j,pi)%t(t)%tz(1)%ll(lp,l)%nn(np,n)=v
                  V_rel(j,pi)%t(t)%tz(-1)%ll(l,lp)%nn(n,np)=v
                  V_rel(j,pi)%t(t)%tz(-1)%ll(lp,l)%nn(np,n)=v
               elseif (num_of_interaction_files==2.and.i_file==2) then
                  v=V_rel(j,pi)%t(t)%tz(0)%ll(l,lp)%nn(n,np)
                  V_rel(j,pi)%t(t)%tz(-1)%ll(l,lp)%nn(n,np)=v
                  V_rel(j,pi)%t(t)%tz(-1)%ll(lp,l)%nn(np,n)=v
               endif
c               print *,'a: s,j,t,tz,  n,l,np,lp, v=',
c     $              s,j,t,tz,n,l,np,lp,v
            endif
         end do
 1000    continue
         close(37)
      end do
      end

      subroutine two_body_state_setup
      use paramdef
      use interaction
      use nodeinfo, only: iproc
      implicit none
      integer :: J,pi,T,sp1,sp2,Jmax,n1,l1,j1,n2,l2,j2,
     $     sp1_max,sp2_min,sp2_max,total,lc,nc,nr,lr,s,jr,tz,
     $     totalrel,lrm,Ntot
      N_sp_int_1max=max(N_sp_min,(N_RGM_max+N_proj_max)/2)
      N_sp_int_2max=max(N_cluster_max,N_RGM_max+N_proj_max)
      N_sp_int_12max=N_sp_int_2max+N_sp_min
c*** This is a truncation for a>1 that appears safe ****
      N_sp_int_12max=max(N_sp_int_2max,N_cluster_max+N_sp_min)
c***
      J_sp_int_max=N_sp_int_12max+1
c*** MPI
      if (iproc==0) then
c*** MPI
         print *,' N_sp_int_1max,N_sp_int_2max,N_sp_int_12max=',
     $     N_sp_int_1max,N_sp_int_2max,N_sp_int_12max
c*** MPI
      endif
c*** MPI
      Jmax=J_sp_int_max
      allocate(tbdst(0:Jmax,0:1))
      allocate(V_2b(0:Jmax,0:1))
      sp1_max=0
      do sp1=1,global_ist1_dim
         n1=global_ist1(1,sp1)
         l1=global_ist1(2,sp1)
         j1=global_ist1(3,sp1)
         if (2*n1+l1>N_sp_int_1max) exit
         sp1_max=sp1_max+1
      end do
c*** MPI
      if (iproc==0) then
c*** MPI
      print *,' sp1_max=',sp1_max
c*** MPI
      endif
c*** MPI
      do J=0,Jmax
         do pi=0,1
            do T=0,1
               total=0
               tbdst(J,pi)%T(T)%sp1max=sp1_max
               allocate(tbdst(J,pi)%T(T)%sp1(sp1_max))
               do sp1=1,sp1_max
                  n1=global_ist1(1,sp1)
                  l1=global_ist1(2,sp1)
                  j1=global_ist1(3,sp1)
                  if (2*n1+l1>N_sp_int_1max) exit
                  sp2_min=global_ist1_dim
                  sp2_max=1
                  do sp2=sp1,global_ist1_dim
                     n2=global_ist1(1,sp2)
                     l2=global_ist1(2,sp2)
                     j2=global_ist1(3,sp2)
                     if (2*n2+l2>N_sp_int_2max) exit   
                     if (2*n1+l1+2*n2+l2>N_sp_int_12max) exit
                     if (mod(l1+l2,2)/=pi) cycle
                     if (abs(j1-j2)>2*J.or.j1+j2<2*J) cycle
                     if (sp1==sp2.and.mod(J+T,2)/=1) cycle
                     sp2_min=min(sp2_min,sp2)
                     sp2_max=max(sp2_max,sp2)
                     total=total+1
                  end do
                  tbdst(J,pi)%T(T)%sp1(sp1)%sp2max=sp2_max
                  tbdst(J,pi)%T(T)%sp1(sp1)%sp2min=sp2_min
                  allocate(tbdst(J,pi)%T(T)%sp1(sp1)%sp2
     $                 (sp2_min:sp2_max))
                  tbdst(J,pi)%T(T)%sp1(sp1)%sp2(:)=-1
               end do
c               print *,' J,pi,T,tbdst%total=',J,pi,T,total
               tbdst(J,pi)%T(T)%total=total
               allocate(tbdst(J,pi)%T(T)%ist(2,total))
               allocate(V_2b(J,pi)%T(T)%Tz(-T:T))
               do tz=-T,T
                  V_2b(J,pi)%T(T)%Tz(tz)%dim=total
                  allocate(V_2b(J,pi)%T(T)%Tz(tz)
     $                 %mat(total*(total+1)/2))
                  V_2b(J,pi)%T(T)%Tz(tz)%mat=0.d0
               end do
               totalrel=0
               allocate(tbdst(J,pi)%T(T)%Ntot(0:N_sp_int_12max))
               do Ntot=0,N_sp_int_12max
                  tbdst(J,pi)%T(T)%Ntot(Ntot)%start=totalrel+1
                  allocate(tbdst(J,pi)%T(T)%Ntot(Ntot)%lc(0:Ntot))
                  do lc=0,Ntot
                     allocate(tbdst(J,pi)%T(T)%Ntot(Ntot)%lc(lc)
     $                    %nc(0:(Ntot-lc)/2))
                     do nc=0,(Ntot-lc)/2
                        allocate(tbdst(J,pi)%T(T)%Ntot(Ntot)%lc(lc)
     $                       %nc(nc)%jr(abs(J-lc):min(J+lc,jrelm)))
                        do jr=abs(J-lc),min(J+lc,jrelm)
                           tbdst(J,pi)%T(T)%Ntot(Ntot)%lc(lc)%nc(nc)
     $                          %jr(jr)%start=totalrel+1
                           do s=0,1
c                           lrm=1-mod(s+T,2)
c                           if (mod(lrm+lc+pi,2)/=0) cycle
                              do lr=abs(jr-s),min(jr+s,Ntot-2*nc-lc)
                                 if (mod(lr+s+T,2)/=1) cycle
                                 if (mod(lr+lc,2)/=pi) cycle
c                                 do nr=0,min((N_cluster_max+N_sp_min
c     $                                -2*nc-lc-lr)/2,(N_cluster_max-lr)/2)
                                 nr=(Ntot-2*nc-lc-lr)/2
                                 if (nr>(N_cluster_max-lr)/2) cycle
                                 if (2*nr+lr+2*nc+lc/=Ntot) cycle
                                 totalrel=totalrel+1
c                                 end do
                              end do
                           end do
                           tbdst(J,pi)%T(T)%Ntot(Ntot)%lc(lc)%nc(nc)
     $                          %jr(jr)%end=totalrel
c                           print *,
c     $                          ' J,pi,T,Ntot,lc,nc,jr,tbdst%start,end:'
c     $                          ,J,pi,T,Ntot,lc,nc,jr,
c     $                          tbdst(J,pi)%T(T)%Ntot(Ntot)%lc(lc)
c     $                          %nc(nc)%jr(jr)%start,
c     $                          tbdst(J,pi)%T(T)%Ntot(Ntot)%lc(lc)
c     $                          %nc(nc)%jr(jr)%end
                        end do
                     end do
                  end do
                  tbdst(J,pi)%T(T)%Ntot(Ntot)%end=totalrel
c                  print *,' J,pi,T,Ntot,tbdst%start,end:',J,pi,T,Ntot,
c     $                 tbdst(J,pi)%T(T)%Ntot(Ntot)%start,
c     $                 tbdst(J,pi)%T(T)%Ntot(Ntot)%end
               end do
c               print *,' J,pi,T,tbdst%totalrel=',J,pi,T,totalrel
               tbdst(J,pi)%T(T)%totalrel=totalrel
               allocate(tbdst(J,pi)%T(T)%istrel(5,totalrel))
            end do
         end do
      end do

      do J=0,Jmax
         do pi=0,1
            do T=0,1
               total=0
               do sp1=1,sp1_max
                  n1=global_ist1(1,sp1)
                  l1=global_ist1(2,sp1)
                  j1=global_ist1(3,sp1)
                  if (2*n1+l1>N_sp_int_1max) exit
                  do sp2=tbdst(J,pi)%T(T)%sp1(sp1)%sp2min,
     $                 tbdst(J,pi)%T(T)%sp1(sp1)%sp2max
                     n2=global_ist1(1,sp2)
                     l2=global_ist1(2,sp2)
                     j2=global_ist1(3,sp2)
                     if (2*n1+l1+2*n2+l2>N_sp_int_12max) exit
                     if (mod(l1+l2,2)/=pi) cycle
                     if (abs(j1-j2)>2*J.or.j1+j2<2*J) cycle
                     if (sp1==sp2.and.mod(J+T,2)/=1) cycle
                     total=total+1
                     tbdst(J,pi)%T(T)%sp1(sp1)%sp2(sp2)=total
                     tbdst(J,pi)%T(T)%ist(1,total)=sp1
                     tbdst(J,pi)%T(T)%ist(2,total)=sp2
c                     print *,' J,pi,T,sp1,sp2,total=',
c     $                    J,pi,T,sp1,sp2,total
                  end do
               end do
               totalrel=0
               do Ntot=0,N_sp_int_12max
                  do lc=0,Ntot
                     do nc=0,(Ntot-lc)/2
                        do jr=abs(J-lc),min(J+lc,jrelm)
                           do s=0,1
c                           lrm=1-mod(s+T,2)
c                           if (mod(lrm+lc+pi,2)/=0) cycle
                              do lr=abs(jr-s),min(jr+s,Ntot-2*nc-lc)
                                 if (mod(lr+s+T,2)/=1) cycle
                                 if (mod(lr+lc,2)/=pi) cycle
                                 nr=(Ntot-2*nc-lc-lr)/2
                                 if (nr>(N_cluster_max-lr)/2) cycle
                                 if (2*nr+lr+2*nc+lc/=Ntot) cycle
c                              do nr=0,min((N_cluster_max+N_sp_min
c     $                             -2*nc-lc-lr)/2,(N_cluster_max-lr)/2)
                                 totalrel=totalrel+1
                                 tbdst(J,pi)%T(T)%istrel(1,totalrel)=nr
                                 tbdst(J,pi)%T(T)%istrel(2,totalrel)=lr
                                 tbdst(J,pi)%T(T)%istrel(3,totalrel)=jr
                                 tbdst(J,pi)%T(T)%istrel(4,totalrel)=nc
                                 tbdst(J,pi)%T(T)%istrel(5,totalrel)=lc
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
c               print *,' J,pi,T,totalrel=',J,pi,T,totalrel
            end do
         end do
      end do

      end

      subroutine V_rel_V_2b_trans
      use paramdef
      use interaction
      use HO_braket
      use nodeinfo
      implicit none
      include 'mpif.h'
      integer :: J,pi,T,sp1,sp2,Jmax,n1,l1,j1,n2,l2,j2,
     $     lc,nc,nr,lr,s,jr,tz,tbst,relst,relsti,nri,lri,si,jri,nci,lci,
     $     tbsti,lambda,sp3,sp4,Ntot
      real(kind(0.d0)),allocatable :: temp(:,:,:) !,sumtrv(:) !,spreltr(:,:)
      real(kind(0.d0)),allocatable :: tmp_in(:),tmp_allred(:)  ! MPI
      integer :: tbst_inc,tbst_step,ierr ! MPI
      real(kind(0.d0)) :: rnorm,c9j,sixj,bmbr,sum,sumtrv(-1:1)
      type sp_rel_tr
      real(kind(0.d0)),allocatable :: rel(:)
      end type sp_rel_tr
      type(sp_rel_tr),allocatable :: spreltr(:)
      interface
         real(kind(0.d0)) function racad(a,b,c,d,e,f)
         integer :: a,b,c,d,e,f
         end function racad
         real(kind(0.d0)) function coef9d(a,b,c,d,e,f,p,q,r)
         integer :: a,b,c,d,e,f,p,q,r
         end function coef9d
      end interface
      Jmax=J_sp_int_max
      call HO_braket_init(N_sp_int_2max,min(N_sp_int_12max,32))
c*** MPI      
      if (iproc==0) print *,' HO_braket_init called'
c*** MPI
      do J=0,Jmax
         do pi=0,1
            do T=0,1
cc               print *,' J,pi,T=',J,pi,T
c               allocate(spreltr(tbdst(J,pi)%T(T)%totalrel,
c     $              tbdst(J,pi)%T(T)%total))
               allocate(spreltr(tbdst(J,pi)%T(T)%total))
               allocate(temp(-T:T,tbdst(J,pi)%T(T)%totalrel,
     $              tbdst(J,pi)%T(T)%total))
cc               allocate(sumtrv(-T:T))
c               spreltr=0.d0
               do tbst=1,tbdst(J,pi)%T(T)%total
                  sp1=tbdst(J,pi)%T(T)%ist(1,tbst)
                  sp2=tbdst(J,pi)%T(T)%ist(2,tbst)
                  n1=global_ist1(1,sp1)
                  l1=global_ist1(2,sp1)
                  n2=global_ist1(1,sp2)
                  l2=global_ist1(2,sp2)
                  Ntot=2*n1+l1+2*n2+l2
                  allocate(spreltr(tbst)%rel(
     $                 tbdst(J,pi)%T(T)%Ntot(Ntot)%start:
     $                 tbdst(J,pi)%T(T)%Ntot(Ntot)%end))
                  spreltr(tbst)%rel=0.d0
               end do
               temp=0.d0
c*** MPI
               do tbst=iproc+1,tbdst(J,pi)%T(T)%total,nproc
c               do tbst=1,tbdst(J,pi)%T(T)%total
c*** MPI
                  sp1=tbdst(J,pi)%T(T)%ist(1,tbst)
                  sp2=tbdst(J,pi)%T(T)%ist(2,tbst)
                  n1=global_ist1(1,sp1)
                  l1=global_ist1(2,sp1)
                  j1=global_ist1(3,sp1)
                  n2=global_ist1(1,sp2)
                  l2=global_ist1(2,sp2)
                  j2=global_ist1(3,sp2)
                  Ntot=2*n1+l1+2*n2+l2
                  if (sp1==sp2) then
                     rnorm=1.d0
                  else
                     rnorm=sqrt(2.d0)
                  endif
c                  print *,' tbst,sp1,sp2,n1,l1,j1,n2,l2,j2=',
c     $                 tbst,sp1,sp2,n1,l1,j1,n2,l2,j2
c                  do relst=1,tbdst(J,pi)%T(T)%totalrel

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(relst,nr,lr,jr,nc,lc,s,
C$OMP&  sum,lambda,
C$OMP&  c9j,sixj,bmbr,relsti,nri,lri,jri,nci,lci,si,tz)
C$OMP& SCHEDULE(DYNAMIC)
                  do relst=tbdst(J,pi)%T(T)%Ntot(Ntot)%start,
     $                 tbdst(J,pi)%T(T)%Ntot(Ntot)%end
                     nr=tbdst(J,pi)%T(T)%istrel(1,relst)
                     lr=tbdst(J,pi)%T(T)%istrel(2,relst)
                     jr=tbdst(J,pi)%T(T)%istrel(3,relst)
                     nc=tbdst(J,pi)%T(T)%istrel(4,relst)
                     lc=tbdst(J,pi)%T(T)%istrel(5,relst)
                     if (2*nr+lr+2*nc+lc/=2*n1+l1+2*n2+l2) then
                        print *,'***error 1 in in V_rel_V_2b_trans'
                        stop
                     endif
                     s=1-mod(lr+T,2)
                     sum=0.d0
                     do lambda=max(abs(lr-lc),abs(l1-l2),abs(J-s)),
     $                    min(lr+lc,l1+l2,J+s)
                        c9j=coef9d(2*l1,1,2*l2,1,j1,j2,2*lambda,2*s,2*J)
                        sixj=racad(2*jr,2*s,2*lc,2*lambda,2*lr,2*J)
c                        call osclbr(nr,lr,nc,lc,n1,l1,n2,l2,lambda,1.d0,
c     $                       bmbr)
                        call osc_br_unit_mass_ratio(nr,lr,nc,lc,
     $                       n1,l1,n2,l2,lambda,bmbr)
                        sum=sum+
     +                       sqrt(real((j1+1)*(j2+1)*(2*s+1)*(2*jr+1)
     $                       ,kind(0.d0)))
     +                       *real(2*lambda+1,kind(0.d0))*c9j*sixj*rnorm
     $                       *bmbr
                     end do
                     spreltr(tbst)%rel(relst)=sum
                  end do
C$OMP END PARALLEL DO
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(relst,nr,lr,jr,nc,lc,s,
C$OMP&  relsti,nri,lri,jri,nci,lci,si,tz,sumtrv)
C$OMP& SCHEDULE(DYNAMIC)
                  do relst=1,tbdst(J,pi)%T(T)%totalrel
                     nr=tbdst(J,pi)%T(T)%istrel(1,relst)
                     lr=tbdst(J,pi)%T(T)%istrel(2,relst)
                     jr=tbdst(J,pi)%T(T)%istrel(3,relst)
                     nc=tbdst(J,pi)%T(T)%istrel(4,relst)
                     lc=tbdst(J,pi)%T(T)%istrel(5,relst)
                     if (2*nc+lc>Ntot) cycle
                     s=1-mod(lr+T,2)
                     sumtrv=0.d0
c                     do relsti=tbdst(J,pi)%T(T)%lc(lc)%nc(nc)%start,
c     $                    tbdst(J,pi)%T(T)%lc(lc)%nc(nc)%end
                     do relsti=tbdst(J,pi)%T(T)%Ntot(Ntot)%lc(lc)%nc(nc)
     $                    %jr(jr)%start,
     $                    tbdst(J,pi)%T(T)%Ntot(Ntot)%lc(lc)%nc(nc)
     $                    %jr(jr)%end
                        nri=tbdst(J,pi)%T(T)%istrel(1,relsti)
                        lri=tbdst(J,pi)%T(T)%istrel(2,relsti)
                        jri=tbdst(J,pi)%T(T)%istrel(3,relsti)
                        nci=tbdst(J,pi)%T(T)%istrel(4,relsti)
                        lci=tbdst(J,pi)%T(T)%istrel(5,relsti)
                        si=1-mod(lri+T,2)
                        if (nc/=nci.or.lc/=lci.or.jr/=jri.or.s/=si) then
                           print *,'*** error 2 in V_rel_V_2b_trans'
                           stop
                        endif
c                        if (jr/=jri.or.s/=si) cycle
c                        if (s/=si) cycle
                        do tz=-T,T
                           sumtrv(tz)=sumtrv(tz)
     $                          +spreltr(tbst)%rel(relsti) !+spreltr(relsti,tbst)
     $                          *V_rel(jr,mod(lr,2))%t(T)%tz(tz)
     $                          %ll(lri,lr)%nn(nri,nr)
                        end do
                     end do
                     do tz=-T,T
                        temp(tz,relst,tbst)=sumtrv(tz)
                     end do
                  end do
C$OMP END PARALLEL DO
               end do
c*** MPI
               if (nproc>1) then
                  call MPI_Barrier(icomm,ierr)
                  tbst=tbdst(J,pi)%T(T)%total
                  if (tbst>0) then
                     allocate(tmp_in(tbst),tmp_allred(tbst))
                     do relst=1,tbdst(J,pi)%T(T)%totalrel
                        do tz=-T,T
                           tmp_in(:)=temp(tz,relst,:)
c                           print *,' before iproc,tz,relst,temp:',iproc,
c     $                          tz,relst,temp(tz,relst,1:3)
                           call MPI_Allreduce(tmp_in(1),tmp_allred(1),
     $                          tbst,MPI_REAL8,MPI_SUM,icomm,ierr)
                           temp(tz,relst,:)=tmp_allred(:)
c                           print *,' after iproc,tz,relst,temp:',iproc,
c     $                          tz,relst,temp(tz,relst,1:3)
                        end do
                     end do
                     deallocate(tmp_in,tmp_allred)
                  endif
                  call MPI_Barrier(icomm,ierr)
                  do tbst=1,tbdst(J,pi)%T(T)%total
                     sp1=tbdst(J,pi)%T(T)%ist(1,tbst)
                     sp2=tbdst(J,pi)%T(T)%ist(2,tbst)
                     n1=global_ist1(1,sp1)
                     l1=global_ist1(2,sp1)
                     n2=global_ist1(1,sp2)
                     l2=global_ist1(2,sp2)
                     Ntot=2*n1+l1+2*n2+l2
                     sp3=tbdst(J,pi)%T(T)%Ntot(Ntot)%start
                     sp4=tbdst(J,pi)%T(T)%Ntot(Ntot)%end-sp3+1
                     if (sp4<1) cycle
                     tbsti=mod(tbst-1,nproc)
c                     print *,' before iproc,tbsti,tbst,spreltr:',
c     $                    iproc,tbst,spreltr(tbst)%rel(sp3),
c     $                    spreltr(tbst)%rel(sp3+sp4-1)
                     call MPI_Bcast(spreltr(tbst)%rel(sp3),
     $                    sp4,MPI_REAL8,tbsti,icomm,ierr)
c                     print *,' after iproc,tbsti,tbst,spreltr:',
c     $                    iproc,tbst,spreltr(tbst)%rel(sp3),
c     $                    spreltr(tbst)%rel(sp3+sp4-1)
                  end do
               endif
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(tbsti,sp3,sp4,n1,l1,
C$OMP&   n2,l2,Ntot,sumtrv,relst,tz,tbst)
C$OMP& SCHEDULE(DYNAMIC)
               do tbsti=1,tbdst(J,pi)%T(T)%total
c               do tbsti=tbst,tbdst(J,pi)%T(T)%total
                  sp3=tbdst(J,pi)%T(T)%ist(1,tbsti)
                  sp4=tbdst(J,pi)%T(T)%ist(2,tbsti)
                  n1=global_ist1(1,sp3)
                  l1=global_ist1(2,sp3)
c     j1=global_ist1(3,sp3)
                  n2=global_ist1(1,sp4)
                  l2=global_ist1(2,sp4)
c     j2=global_ist1(3,sp4)
                  Ntot=2*n1+l1+2*n2+l2

c*** MPI
c               print *,' iproc=',iproc,' spreltr,temp set'
                  do tbst=iproc+1,tbsti,nproc
c                  do tbst=iproc+1,tbdst(J,pi)%T(T)%total,nproc
c               do tbst=1,tbdst(J,pi)%T(T)%total
c*** MPI

                     sumtrv=0.d0
c                     do relst=1,tbdst(J,pi)%T(T)%totalrel
                     do relst=tbdst(J,pi)%T(T)%Ntot(Ntot)%start,
     $                    tbdst(J,pi)%T(T)%Ntot(Ntot)%end
                        do tz=-T,T
                           sumtrv(tz)=sumtrv(tz)+temp(tz,relst,tbst)
     $                          *spreltr(tbsti)%rel(relst)
c     $                       *spreltr(relst,tbsti)
                        end do
                     end do
                     do tz=-T,T
                        V_2b(J,pi)%T(T)%Tz(tz)%mat(
     $                       tbst+tbsti*(tbsti-1)/2)=sumtrv(tz)
c                        print *,' J,pi,T,Tz,a,b,c,d,V=',
c     $                       J,pi,T,tz,sp1,sp2,sp3,sp4,
c     $                       V_2b(J,pi)%T(T)%Tz(tz)%mat(
c     $                       tbst+tbsti*(tbsti-1)/2)
                     end do
c*** MPI
                  end do
c*** MPI
               end do
C$OMP END PARALLEL DO
               deallocate(temp)
cc               deallocate(sumtrv)
               deallocate(spreltr)
c*** MPI
               if (nproc>1) then
                  call MPI_Barrier(icomm,ierr)
                  tbst=tbdst(J,pi)%T(T)%total
                  if (mod(tbst,2)==0) then
                     tbst_inc=tbst+1
                     tbst_step=tbst/2
                  else
                     tbst_inc=tbst
                     tbst_step=(tbst+1)/2
                  endif
                  if (tbst_inc>0) then
                     allocate(tmp_in(tbst_inc),tmp_allred(tbst_inc))
                     do relst=0,tbst_step-1
                        do tz=-T,T
                           tmp_in(:)=V_2b(J,pi)%T(T)%Tz(tz)%mat(
     $                          tbst_inc*relst+1:
     $                          tbst_inc*relst+tbst_inc)
                           call MPI_Allreduce(tmp_in(1),tmp_allred(1),
     $                          tbst_inc,MPI_REAL8,MPI_SUM,icomm,ierr)
                           V_2b(J,pi)%T(T)%Tz(tz)%mat(tbst_inc*relst+1:
     $                          tbst_inc*relst+tbst_inc)=tmp_allred(:)
                        end do
                     end do
                     deallocate(tmp_in,tmp_allred)
                  endif
               endif
c*** MPI
            end do
         end do
      end do
      deallocate(V_rel)
      call HO_braket_destroy
      end

      subroutine get_V2me_uncoup(ialfa,ibeta,igamma,idelta,V2me_uncoupl) 
!this subroutine calls the subroutine get_V2me
c****** HERE MODULES
      use paramdef
      use interaction
      use nodeinfo
c****** DEFINITIONS
      implicit none
      integer,intent(IN) :: ialfa,ibeta,igamma,idelta
      real(kind(0.d0)),intent(OUT) :: V2me_uncoupl
      real(kind(0.d0)) :: V2me
      integer :: ia,na,la,ja,mja,mta,ib,nb,lb,jb,mjb,mtb,
     $     ic,nc,lc,jc,mjc,mtc,id,nd,ld,jd,mjd,mtd
      integer :: Jtot,M_Jtot,Ttot,M_Ttot
      real(kind(0.d0)) :: clebd,cleb1,cleb2,clebt1,clebt2
      integer :: iaa,ibb,icc,idd,JJ,TT,TTz_in
      interface
         pure subroutine get_V2me(iaa,ibb,icc,idd,JJ,TT,TTz_in,VV2me)
         integer,intent(IN) :: iaa,ibb,icc,idd,JJ,TT,TTz_in
         real(kind(0.d0)),intent(OUT) :: VV2me
         end subroutine get_V2me
      end interface
c******
ccc norm test
c      if (ialfa==igamma.and.ibeta==idelta) then
c         V2me_uncoupl=1.d0
c         return
c      elseif (ialfa==idelta.and.igamma==ibeta) then
c         V2me_uncoupl=-1.d0
c         return
c      else
c         V2me_uncoupl=0.d0
c         return
c      endif
ccc norm test
      ia=global_ist1_sp(1,ialfa)
      mja=global_ist1_sp(2,ialfa)
      mta=global_ist1_sp(3,ialfa)

      na=global_ist1(1,ia)
      la=global_ist1(2,ia)
      ja=global_ist1(3,ia)

      ib=global_ist1_sp(1,ibeta)
      mjb=global_ist1_sp(2,ibeta)
      mtb=global_ist1_sp(3,ibeta)

      nb=global_ist1(1,ib)
      lb=global_ist1(2,ib)
      jb=global_ist1(3,ib)

      ic=global_ist1_sp(1,igamma)
      mjc=global_ist1_sp(2,igamma)
      mtc=global_ist1_sp(3,igamma)

      nc=global_ist1(1,ic)
      lc=global_ist1(2,ic)
      jc=global_ist1(3,ic)

      id=global_ist1_sp(1,idelta)
      mjd=global_ist1_sp(2,idelta)
      mtd=global_ist1_sp(3,idelta)

      nd=global_ist1(1,id)
      ld=global_ist1(2,id)
      jd=global_ist1(3,id)

      if (2*na+la+2*nb+lb>N_sp_int_12max.or.
     $     2*nc+lc+2*nd+ld>N_sp_int_12max
     $     .or.max(2*na+la,2*nb+lb,2*nc+lc,2*nd+ld)>N_sp_int_2max
     $     .or.min(2*na+la,2*nb+lb)>N_sp_int_1max
     $     .or.min(2*nc+lc,2*nd+ld)>N_sp_int_1max
     $     .or.(mja+mjb).ne.(mjc+mjd)  !same M_Jtot on the bra and ket
     $     .or.(mta+mtb).ne.(mtc+mtd))  then !same M_Ttot on the bra and ket
         V2me_uncoupl=0.d0
         return
      endif

      M_Jtot=mja+mjb
      M_Ttot=mta+mtb

!      if (iproc==0) print *,' ia,ib,ic,id,ja,jb,jc,jd,M_Jtot,M_Ttot=',
!     $     ia,ib,ic,id,ja,jb,jc,jd,M_Jtot,M_Ttot

      V2me_uncoupl=0.d0
      do Jtot=max(abs(ja-jb),abs(jc-jd),abs(M_Jtot)),
     $     min(ja+jb,jc+jd),2

         cleb1=clebd(ja,mja,jb,mjb,Jtot,M_Jtot)
         cleb2=clebd(jc,mjc,jd,mjd,Jtot,M_Jtot)

         do Ttot=abs(M_Ttot),2,2

            if ((ia==ib.or.ic==id).and.mod((Jtot+Ttot)/2,2)/=1) cycle

            clebt1=clebd(1,mta,1,mtb,Ttot,M_Ttot)
            clebt2=clebd(1,mtc,1,mtd,Ttot,M_Ttot)

!            if (iproc==0) print *,
!     $           ' Jtot,Ttot,cleb1,cleb2,clebt1,clebt2=',
!     $           Jtot,Ttot,cleb1,cleb2,clebt1,clebt2

            call get_V2me(ia,ib,ic,id,Jtot/2,Ttot/2,M_Ttot/2,V2me)

!            if (iproc==0) print *,' V2me=',V2me

            V2me_uncoupl=V2me_uncoupl
     $           +cleb1*cleb2*clebt1*clebt2*V2me

         end do
      end do
      end subroutine get_V2me_uncoup

      subroutine get_V2me(ia,ib,ic,id,J,T,Tz_in,V2me)
c*** It is assumed that compatibility checks were done before calling this
      use paramdef, only: global_ist1,N_cluster_max,N_sp_min
      use interaction
      implicit none
      integer,intent(IN) :: ia,ib,ic,id,J,T,Tz_in
      real(kind(0.d0)),intent(OUT) :: V2me
      integer :: Tz
      integer :: na,la,ja,nb,lb,jb,nc,lc,jc,nd,ld,jd,tab,tcd,
     $     ii_mat,MTot2
      real(kind(0.d0)) :: phase_ab,phase_cd,sq_ab,sq_cd
      real(kind(0.d0)) :: cpp,cnp,cnn
      Tz=Tz_in
c***test
c      Tz=0
c***test
      if (ia==ib) then
         sq_ab=sqrt(2.d0)
      else
         sq_ab=1.d0
      endif
      if (ic==id) then
         sq_cd=sqrt(2.d0)
      else
         sq_cd=1.d0
      endif
      na=global_ist1(1,ia)
      la=global_ist1(2,ia)
      ja=global_ist1(3,ia)
      nb=global_ist1(1,ib)
      lb=global_ist1(2,ib)
      jb=global_ist1(3,ib)
      nc=global_ist1(1,ic)
      lc=global_ist1(2,ic)
      jc=global_ist1(3,ic)
      nd=global_ist1(1,id)
      ld=global_ist1(2,id)
      jd=global_ist1(3,id)
      if (2*na+la+2*nb+lb>N_sp_int_12max.or.
     $     2*nc+lc+2*nd+ld>N_sp_int_12max
     $     .or.max(2*na+la,2*nb+lb,2*nc+lc,2*nd+ld)>N_sp_int_2max
     $     .or.min(2*na+la,2*nb+lb)>N_sp_int_1max
     $     .or.min(2*nc+lc,2*nd+ld)>N_sp_int_1max) then
         V2me=0.d0
         return
      endif
      if (ia>ib) then
         tab=tbdst(J,mod(la+lb,2))%T(T)%sp1(ib)%sp2(ia)
         phase_ab=real((-1)**(J+T+(ja+jb)/2),kind(0.d0))
      else
         tab=tbdst(J,mod(la+lb,2))%T(T)%sp1(ia)%sp2(ib)
         phase_ab=1.d0
      endif
      if (ic>id) then
         tcd=tbdst(J,mod(lc+ld,2))%T(T)%sp1(id)%sp2(ic)
         phase_cd=real((-1)**(J+T+(jc+jd)/2),kind(0.d0))
      else
         tcd=tbdst(J,mod(lc+ld,2))%T(T)%sp1(ic)%sp2(id)
         phase_cd=1.d0
      endif
      if (tab<tcd) then
         ii_mat=tab+tcd*(tcd-1)/2
      else
         ii_mat=tcd+tab*(tab-1)/2
      endif
      V2me=V_2b(J,mod(la+lb,2))%T(T)%Tz(Tz)%mat(ii_mat)
      V2me=V2me*phase_ab*phase_cd*sq_ab*sq_cd
c***test***
cc      if (.not.(ia==1.and.ib==1
cc     $     .and.((ic==1.and.id==4).or.(ic==4.and.id==1)))) V2me=0.d0      
c      if (tab==tcd) then
c         V2me=phase_ab*phase_cd*sq_ab*sq_cd
c      else
c         V2me=0.d0
c      endif
c***test***
      end subroutine get_V2me

      subroutine threebodysetup_cJ
c** Setup for a three-body interaction ***
c** Petr Navratil, 30 June 2011, TRIUMF *****
      use paramdef
      use nodeinfo
      use intrface, only: iunitvout
      use initst
      use finast
      use v3b
      implicit none
      include 'mpif.h'
      integer :: J_3,T_3,ini,fin
      integer :: nhom,nhom2sp,ierr
      integer :: ntot_a,n_a,l_a,j2_a,isp1n_a
      integer :: ntot_b,n_b,l_b,j2_b,isp1n_b
      integer :: ntot_c,n_c,l_c,j2_c,isp1n_c
      integer :: ntot_d,n_d,l_d,j2_d,isp1n_d
      integer :: ntot_e,n_e,l_e,j2_e,isp1n_e,isp1n_e_max
      integer :: ntot_f,n_f,l_f,j2_f,isp1n_f,isp1n_f_max
      integer :: iii,j_ab,t_ab,j_de,t_de,i_abc,i_def
      integer :: i,ibuf,pi_i,pi_f
      integer :: Nmin_HO
      integer(8) ii,i_abcdef,tot_num_of_3bmatel,iix,ih3oibuf
      logical :: compfile

      nhom2sp=N12_max
      nhom=N123_max
      if (iproc==0) print *,' V3N: N1_max,N12_max,N123_max=',
     $     N1_max,N12_max,N123_max

      call sp1nbas_cJ(N1_max)

      allocate(index_abc(isp1ntot,isp1ntot,isp1ntot))
      iii=0
      do isp1n_a=1,isp1ntot
         n_a=isp1n_cJ(1,isp1n_a)
         l_a=isp1n_cJ(2,isp1n_a)
         j2_a=isp1n_cJ(3,isp1n_a)
         ntot_a=2*n_a+l_a
         if (ntot_a>nhom) cycle
         do isp1n_b=1,isp1n_a
            n_b=isp1n_cJ(1,isp1n_b)
            l_b=isp1n_cJ(2,isp1n_b)
            j2_b=isp1n_cJ(3,isp1n_b)
            ntot_b=2*n_b+l_b
            if (ntot_a+ntot_b>nhom2sp) cycle
            do isp1n_c=1,isp1n_b
               n_c=isp1n_cJ(1,isp1n_c)
               l_c=isp1n_cJ(2,isp1n_c)
               j2_c=isp1n_cJ(3,isp1n_c)
               ntot_c=2*n_c+l_c
               if (ntot_a+ntot_c>nhom2sp) cycle
               if (ntot_b+ntot_c>nhom2sp) cycle
               if (ntot_a+ntot_b+ntot_c>nhom) cycle
               iii=iii+1
               index_abc(isp1n_a,isp1n_b,isp1n_c)=iii
            end do
         end do
      end do
      dim_abc=iii
c      print *,' dim_abc=',dim_abc
      allocate(index_abcdef(dim_abc,dim_abc))
      index_abcdef=0
      ii=0
      do isp1n_a=1,isp1ntot
         n_a=isp1n_cJ(1,isp1n_a)
         l_a=isp1n_cJ(2,isp1n_a)
         j2_a=isp1n_cJ(3,isp1n_a)
         ntot_a=2*n_a+l_a
         if (ntot_a>nhom) cycle
         do isp1n_b=1,isp1n_a
            n_b=isp1n_cJ(1,isp1n_b)
            l_b=isp1n_cJ(2,isp1n_b)
            j2_b=isp1n_cJ(3,isp1n_b)
            ntot_b=2*n_b+l_b
            if (ntot_a+ntot_b>nhom2sp) cycle
            do isp1n_c=1,isp1n_b
               n_c=isp1n_cJ(1,isp1n_c)
               l_c=isp1n_cJ(2,isp1n_c)
               j2_c=isp1n_cJ(3,isp1n_c)
               ntot_c=2*n_c+l_c
               if (ntot_a+ntot_c>nhom2sp) cycle
               if (ntot_b+ntot_c>nhom2sp) cycle
               if (ntot_a+ntot_b+ntot_c>nhom) cycle
               i_abc=index_abc(isp1n_a,isp1n_b,isp1n_c)
c               order_abc=(isp1n_a-1)*((isp1ntot-1)*(isp1ntot+1)+1)
c     $              +(isp1n_b-1)*isp1ntot+isp1n_c-1

c               do isp1n_d=1,isp1ntot
               do isp1n_d=1,isp1n_a
                  n_d=isp1n_cJ(1,isp1n_d)
                  l_d=isp1n_cJ(2,isp1n_d)
                  j2_d=isp1n_cJ(3,isp1n_d)
                  ntot_d=2*n_d+l_d
                  if (ntot_d>nhom) cycle
                  if (isp1n_d==isp1n_a) then
                     isp1n_e_max=isp1n_b
                  else
                     isp1n_e_max=isp1n_d
                  endif
c                  do isp1n_e=1,isp1n_d
                  do isp1n_e=1,isp1n_e_max
                     n_e=isp1n_cJ(1,isp1n_e)
                     l_e=isp1n_cJ(2,isp1n_e)
                     j2_e=isp1n_cJ(3,isp1n_e)
                     ntot_e=2*n_e+l_e
                     if (ntot_d+ntot_e>nhom2sp) cycle
                     if (isp1n_d==isp1n_a.and.isp1n_e==isp1n_b) then
                        isp1n_f_max=isp1n_c
                     else
                        isp1n_f_max=isp1n_e
                     endif
c                     do isp1n_f=1,isp1n_e
                     do isp1n_f=1,isp1n_f_max
                        n_f=isp1n_cJ(1,isp1n_f)
                        l_f=isp1n_cJ(2,isp1n_f)
                        j2_f=isp1n_cJ(3,isp1n_f)
                        ntot_f=2*n_f+l_f
                        if (ntot_d+ntot_f>nhom2sp) cycle
                        if (ntot_e+ntot_f>nhom2sp) cycle
                        if (ntot_d+ntot_e+ntot_f>nhom) cycle
                        if (mod(l_a+l_b+l_c+l_d+l_e+l_f,2)==1) cycle
                        i_def=index_abc(isp1n_d,isp1n_e,isp1n_f)
c                        order_def=(isp1n_d-1)
c     $                       *((isp1ntot-1)*(isp1ntot+1)+1)
c     $                       +(isp1n_e-1)*isp1ntot+isp1n_f-1
c                        if (order_def>order_abc) then
                        if (i_def>i_abc) then
                           print *,'a,b,c=',isp1n_a,isp1n_b,isp1n_c
                           print *,'e,d,f=',isp1n_e,isp1n_d,isp1n_f
c                           print *,' order_abc=',order_abc
c                           print *,' order_def=',order_def
                           print *,' i_abc=',i_abc
                           print *,' i_def=',i_def
                           stop
                        endif
cc                        if (i_def<i_abc) cycle
                        ii=ii+1
                        index_abcdef(i_abc,i_def)=ii
                     end do
                  end do
               end do
            end do
         end do
      end do
      dim_abcdef=ii
c      print *,' dim_abcdef=',dim_abcdef
      allocate(start_abcdef(dim_abcdef))
      start_abcdef=0
      ii=0
      do isp1n_a=1,isp1ntot
         n_a=isp1n_cJ(1,isp1n_a)
         l_a=isp1n_cJ(2,isp1n_a)
         j2_a=isp1n_cJ(3,isp1n_a)
         ntot_a=2*n_a+l_a
         if (ntot_a>nhom) cycle
         do isp1n_b=1,isp1n_a
            n_b=isp1n_cJ(1,isp1n_b)
            l_b=isp1n_cJ(2,isp1n_b)
            j2_b=isp1n_cJ(3,isp1n_b)
            ntot_b=2*n_b+l_b
            if (ntot_a+ntot_b>nhom2sp) cycle
            do isp1n_c=1,isp1n_b
               n_c=isp1n_cJ(1,isp1n_c)
               l_c=isp1n_cJ(2,isp1n_c)
               j2_c=isp1n_cJ(3,isp1n_c)
               ntot_c=2*n_c+l_c
               if (ntot_a+ntot_c>nhom2sp) cycle
               if (ntot_b+ntot_c>nhom2sp) cycle
               if (ntot_a+ntot_b+ntot_c>nhom) cycle
               i_abc=index_abc(isp1n_a,isp1n_b,isp1n_c)
c               do isp1n_d=1,isp1ntot
               do isp1n_d=1,isp1n_a
                  n_d=isp1n_cJ(1,isp1n_d)
                  l_d=isp1n_cJ(2,isp1n_d)
                  j2_d=isp1n_cJ(3,isp1n_d)
                  ntot_d=2*n_d+l_d
                  if (ntot_d>nhom) cycle
                  if (isp1n_d==isp1n_a) then
                     isp1n_e_max=isp1n_b
                  else
                     isp1n_e_max=isp1n_d
                  endif
c                  do isp1n_e=1,isp1n_d
                  do isp1n_e=1,isp1n_e_max
                     n_e=isp1n_cJ(1,isp1n_e)
                     l_e=isp1n_cJ(2,isp1n_e)
                     j2_e=isp1n_cJ(3,isp1n_e)
                     ntot_e=2*n_e+l_e
                     if (ntot_d+ntot_e>nhom2sp) cycle
                     if (isp1n_d==isp1n_a.and.isp1n_e==isp1n_b) then
                        isp1n_f_max=isp1n_c
                     else
                        isp1n_f_max=isp1n_e
                     endif
c                     do isp1n_f=1,isp1n_e
                     do isp1n_f=1,isp1n_f_max
                        n_f=isp1n_cJ(1,isp1n_f)
                        l_f=isp1n_cJ(2,isp1n_f)
                        j2_f=isp1n_cJ(3,isp1n_f)
                        ntot_f=2*n_f+l_f
                        if (ntot_d+ntot_f>nhom2sp) cycle
                        if (ntot_e+ntot_f>nhom2sp) cycle
                        if (ntot_d+ntot_e+ntot_f>nhom) cycle
                        if (mod(l_a+l_b+l_c+l_d+l_e+l_f,2)==1) cycle
                        i_def=index_abc(isp1n_d,isp1n_e,isp1n_f)
cc                        if (i_def<i_abc) cycle
                        i_abcdef=index_abcdef(i_abc,i_def)
                        start_abcdef(i_abcdef)=ii+1
                        do j_ab=abs(j2_a-j2_b)/2,(j2_a+j2_b)/2
                           do j_de=abs(j2_d-j2_e)/2,(j2_d+j2_e)/2
                              do J_3=max(abs(2*j_ab-j2_c),
     $                             abs(2*j_de-j2_f)),
     $                             min(2*j_ab+j2_c,2*j_de+j2_f),2
                                 do t_ab=0,1
c                                    if (isp1n_a==isp1n_b.and.
c     $                                   mod(j_ab+t_ab,2)==0) cycle
                                    do t_de=0,1
c                                       if (isp1n_d==isp1n_e.and.
c     $                                      mod(j_de+t_de,2)==0) cycle

                                       do T_3=max(abs(2*t_ab-1),
     $                                      abs(2*t_de-1)),
     $                                      min(2*t_ab+1,2*t_de+1),2
                                          ii=ii+1
                                       end do
                                    end do
                                 end do
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
      tot_num_of_3bmatel=ii
      if (iproc==0) then
         print *,' tot_num_of_3bmatel=',tot_num_of_3bmatel
c         write(2,1000) tot_num_of_3bmatel
c 1000    format(/,' Number of 3N matrix elements in threebodysetup_cJ:',
c     $        i8)
      endif
      allocate(v3b_cJ(tot_num_of_3bmatel))
      if (iproc==0) then

         inquire(file=trim(v3intfile)//'_comp',exist=compfile)
         if (compfile) then

            open(27,file=trim(v3intfile)//'_comp',status='old',
     $           form='unformatted',action='read')
            write(iunitvout,*)
            write(iunitvout,'(a)')
     $           ' 3N interaction file '//trim(v3intfile)//'_comp'
            write(iunitvout,
     $           "(' N1_max=',i4,'    N12_max=',i4,'    N123_max=',i4)")
     $           N1_max,N12_max,N123_max

            read(27) (v3b_cJ(ii),ii=1,tot_num_of_3bmatel)
            close(27)

         else

            open(27,file=trim(v3intfile),status='old',
     $           form='unformatted',action='read')
            write(iunitvout,*)
            write(iunitvout,'(a)')
     $           ' 3N interaction file '//trim(v3intfile)
            write(iunitvout,
     $           "(' N1_max=',i4,'    N12_max=',i4,'    N123_max=',i4)")
     $           N1_max,N12_max,N123_max
            ii=0
         do isp1n_a=1,isp1ntot
            n_a=isp1n_cJ(1,isp1n_a)
            l_a=isp1n_cJ(2,isp1n_a)
            j2_a=isp1n_cJ(3,isp1n_a)
            ntot_a=2*n_a+l_a
            if (ntot_a>nhom) cycle
            do isp1n_b=1,isp1n_a
               n_b=isp1n_cJ(1,isp1n_b)
               l_b=isp1n_cJ(2,isp1n_b)
               j2_b=isp1n_cJ(3,isp1n_b)
               ntot_b=2*n_b+l_b
               if (ntot_a+ntot_b>nhom2sp) cycle
               do isp1n_c=1,isp1n_b
                  n_c=isp1n_cJ(1,isp1n_c)
                  l_c=isp1n_cJ(2,isp1n_c)
                  j2_c=isp1n_cJ(3,isp1n_c)
                  ntot_c=2*n_c+l_c
                  if (ntot_a+ntot_c>nhom2sp) cycle
                  if (ntot_b+ntot_c>nhom2sp) cycle
                  if (ntot_a+ntot_b+ntot_c>nhom) cycle
                  i_abc=index_abc(isp1n_a,isp1n_b,isp1n_c)
                  pi_i=(-1)**(l_a+l_b+l_c)
c                  do isp1n_d=1,isp1ntot
                  do isp1n_d=1,isp1n_a
                     n_d=isp1n_cJ(1,isp1n_d)
                     l_d=isp1n_cJ(2,isp1n_d)
                     j2_d=isp1n_cJ(3,isp1n_d)
                     ntot_d=2*n_d+l_d
                     if (ntot_d>nhom) cycle
                     if (isp1n_d==isp1n_a) then
                        isp1n_e_max=isp1n_b
                     else
                        isp1n_e_max=isp1n_d
                     endif
c                     do isp1n_e=1,isp1n_d
                     do isp1n_e=1,isp1n_e_max
                        n_e=isp1n_cJ(1,isp1n_e)
                        l_e=isp1n_cJ(2,isp1n_e)
                        j2_e=isp1n_cJ(3,isp1n_e)
                        ntot_e=2*n_e+l_e
                        if (ntot_d+ntot_e>nhom2sp) cycle
                        if (isp1n_d==isp1n_a.and.isp1n_e==isp1n_b) then
                           isp1n_f_max=isp1n_c
                        else
                           isp1n_f_max=isp1n_e
                        endif
c                     do isp1n_f=1,isp1n_e
                        do isp1n_f=1,isp1n_f_max
                           n_f=isp1n_cJ(1,isp1n_f)
                           l_f=isp1n_cJ(2,isp1n_f)
                           j2_f=isp1n_cJ(3,isp1n_f)
                           ntot_f=2*n_f+l_f
                           if (ntot_d+ntot_f>nhom2sp) cycle
                           if (ntot_e+ntot_f>nhom2sp) cycle
                           if (ntot_d+ntot_e+ntot_f>nhom) cycle
                           pi_f=(-1)**(l_d+l_e+l_f)
                           if (pi_i/=pi_f) cycle
                           i_def=index_abc(isp1n_d,isp1n_e,isp1n_f)
cc                           if (i_def<i_abc) cycle
c     i_abcdef=index_abcdef(i_abc,i_def)
c                        start_abcdef(i_abcdef)=ii+1
                           do j_ab=abs(j2_a-j2_b)/2,(j2_a+j2_b)/2
                              do j_de=abs(j2_d-j2_e)/2,(j2_d+j2_e)/2
                                 do J_3=max(abs(2*j_ab-j2_c),
     $                                abs(2*j_de-j2_f)),
     $                                min(2*j_ab+j2_c,2*j_de+j2_f),2
                                    do t_ab=0,1
c                                       if (isp1n_a==isp1n_b.and.
c     $                                      mod(j_ab+t_ab,2)==0) cycle
                                       do t_de=0,1
c                                          if (isp1n_d==isp1n_e.and.
c     $                                      mod(j_de+t_de,2)==0) cycle
                                          
                                          do T_3=max(abs(2*t_ab-1),
     $                                         abs(2*t_de-1)),
     $                                         min(2*t_ab+1,2*t_de+1),2
                                             ii=ii+1
                                             read(27) v3b_cJ(ii)
c     print *, isp1n_a,isp1n_b,
c     $                                         isp1n_c,isp1n_d,isp1n_e,
c     $                                         isp1n_f,j_ab,t_ab,
c     $                                         j_de,t_de,J_3,T_3,
c     $                                         v3b_cJ(ii)
                                          end do
                                       end do
                                    end do
                                 end do
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
            close(27)
         endif
      endif
      if (iproc==0) print *,' v3b_cJ read'
      call MPI_Barrier(icomm,ierr)
      ibuf=3000000
      if (tot_num_of_3bmatel<=ibuf) then
         i=tot_num_of_3bmatel
         call MPI_Bcast(v3b_cJ(1),i,MPI_REAL,0,icomm,ierr) 
      else
         ih3oibuf=tot_num_of_3bmatel/ibuf
         iix=1
         do i=1,ih3oibuf
            call MPI_Bcast(v3b_cJ(iix),ibuf,MPI_REAL,0,icomm,ierr) 
            if (iproc==0) print *,' Broadcasted iix=',iix
            iix=iix+ibuf
         end do
         i=mod(tot_num_of_3bmatel,int(ibuf,kind(8)))
         if (i/=0) then
            call MPI_Bcast(v3b_cJ(iix),i,MPI_REAL,0,icomm,ierr)
            if (iproc==0) print *,' Broadcasted iix=',iix
         endif
      endif
      if (iproc==0) print *,' H3 broadcasted'
      call cg_init(N1_max,nhom2sp,nhom)
      if (iproc==0) print *,' cg_init called'
      end

      subroutine sp1nbas_cJ(nhom)
      use paramdef
      use v3b
      use nodeinfo
      implicit none
      include 'mpif.h'
      integer,intent(IN) :: nhom
      integer :: n,l,j2,j,n1,l1,j1
      integer :: ii,ntot
      ii=0
      do ntot=0,nhom
         do l=mod(ntot,2),ntot,2
            n=(ntot-l)/2
            do j2=iabs(2*l-1),2*l+1,2
               ii=ii+1
            end do
         end do
      end do
      isp1ntot=ii
      if (isp1ntot<global_ist1_dim) then
         if (iproc==0) print *,
     $        '***warning: in input isp1ntot,global_ist1_dim=',
     $        isp1ntot,global_ist1_dim
cc         stop
      endif
      if (allocated(isp1n_cJ)) deallocate(isp1n_cJ)
      allocate(isp1n_cJ(3,isp1ntot))
      if (iproc==0) then
         write(*,1000) isp1ntot
 1000    format(' Number of nlj states:',i5)
      endif
      ii=0
      do ntot=0,nhom
         do l=mod(ntot,2),ntot,2
            n=(ntot-l)/2
            do j2=iabs(2*l-1),2*l+1,2
               ii=ii+1
               isp1n_cJ(1,ii)=n
               isp1n_cJ(2,ii)=l
               isp1n_cJ(3,ii)=j2
c               write(*,1100) ii,n,l,j2
c 1100          format(' #',i5,'  n=',i3,'  l=',i3,
c     +              '  j=',i3,'/2')
            end do
         end do
      end do
      do ii=1,min(global_ist1_dim,isp1ntot)
         n=global_ist1(1,ii)
         l=global_ist1(2,ii)
         j=global_ist1(3,ii)
         n1=isp1n_cJ(1,ii)
         l1=isp1n_cJ(2,ii)
         j1=isp1n_cJ(3,ii)
         if (n/=n1.or.l/=l1.or.j/=j1) then
            print *,'***error iproc,n,l,j,n1,l1,j1:',
     $           iproc,n,l,j,n1,l1,j1
            stop
         endif
      end do
      end

      subroutine cg_init(N1max,N12max,N123max)
      use v3b
      implicit none
      integer,intent(IN) :: N1max,N12max,N123max
      integer :: mt1,mt2,t12,mt12,mt3,T3
      integer :: j1,m1,j2,m2,j12,m12,j3,m3,j123,j1max,j12max,j123max
      real(kind(0.d0)) :: clebd
c      allocate(cgj12())
      cgt12=0.d0
      cgt123=0.d0
      do mt1=-1,1,2
cc         print *, ' mt1,mt1/2=',mt1,mt1/2
         do mt2=-1,1,2
            do t12=abs(mt1+mt2),2,2
               cgt12(t12/2,(mt1+1)/2,(mt2+1)/2)=
     $              clebd(1,mt1,1,mt2,t12,mt1+mt2)
            end do
         end do
      end do
      do t12=0,2,2
         do mt12=-t12,t12,2
            do mt3=-1,1,2
               do T3=max(abs(t12-1),abs(mt12+mt3)),t12+1,2
                  cgt123(t12/2,T3/2,mt12/2,(mt3+1)/2)=
     $                 clebd(t12,mt12,1,mt3,T3,mt12+mt3)
               end do
            end do
         end do
      end do

      j1max=2*N1max+1
      j12max=2*(N12max+1)
      j123max=2*N123max+3

      allocate(cgj12(0:j12max/2,0:j1max/2,0:j1max,0:j1max/2,0:j1max))
      allocate(cgj123(0:j123max/2,0:j12max/2,-j12max/2:j12max/2,
     $     0:j1max/2,0:j1max))
      cgj12=0.d0
      cgj123=0.d0
      do j1=1,j1max,2
         do m1=-j1,j1,2
            do j2=1,j1max,2
               do m2=-j2,j2,2
                  do j12=max(abs(j1-j2),abs(m1+m2)),min(j1+j2,j12max),2
                     cgj12(j12/2,j1/2,(j1+m1)/2,j2/2,(j2+m2)/2)=
     $                    clebd(j1,m1,j2,m2,j12,m1+m2)
                  end do
               end do
            end do
         end do
      end do

      do j12=0,j12max,2
         do m12=-j12,j12,2
            do j3=1,j1max,2
               do m3=-j3,j3,2
                  do j123=max(abs(j12-j3),abs(m12+m3)),
     $                 min(j12+j3,j123max),2
                     cgj123(j123/2,j12/2,m12/2,j3/2,(j3+m3)/2)=
     $                    clebd(j12,m12,j3,m3,j123,m12+m3)
                  end do
               end do
            end do
         end do
      end do
      end

      real(kind(0.d0)) function v3b_cJ_unc(ispcr_1_in,ispcr_2_in,
     $     ispcr_3_in,ispan_1_in,ispan_2_in,ispan_3_in)
      use paramdef
      use nodeinfo
      use v3b
      implicit none
      integer,intent(IN) :: ispcr_1_in,ispcr_2_in,
     +           ispcr_3_in,ispan_1_in,ispan_2_in,ispan_3_in
      integer :: ispcr_1,ispcr_2,
     +           ispcr_3,ispan_1,ispan_2,ispan_3
      integer :: iphase,itemp,mjtot2,mttot2,itemp_i !,ipar
      integer:: i_a,i_b,i_c,i_d,i_e,i_f,i_abc,i_def
c     $     ,order_abc,order_def
      integer(8) :: i_abcdef,iii
      real(kind(0.d0)) :: cg1,cg2,cg3,cg4,cgt1(0:1),cgt2(0:1),
     $     cgt3(0:1,0:1),cgt4(0:1,0:1),clebd,sum3
      integer :: j2_a,j2_b,j2_c,j2_d,j2_e,j2_f,j_ab,j_de,t_ab,t_de,
     $     J_3,T_3,m2_a,m2_b,m2_c,m2_d,m2_e,m2_f,mt2_a,mt2_b,mt2_c,
     $     mt2_d,mt2_e,mt2_f
      integer :: mja,mjb,mjc,mjd,mje,mjf,mta,mtb,mtc,mtd,mte,mtf,
     $     ja,jb,jc,jd,je,jf

      ispcr_1=ispcr_1_in
      ispcr_2=ispcr_2_in
      ispcr_3=ispcr_3_in
      ispan_1=ispan_1_in
      ispan_2=ispan_2_in
      ispan_3=ispan_3_in

c      if (mjtot2/=m2_sp(ispcr_1)+m2_sp(ispcr_2)+m2_sp(ispcr_3)) then
c         ham3b=0.0
c         return
c      endif
c      if (mttot2/=mt2_sp(ispcr_1)+mt2_sp(ispcr_2)+mt2_sp(ispcr_3)) then
c         ham3b=0.0
c         return
c      endif
c      ipar=mod(l_sp(ispan_1)+l_sp(ispan_2)+l_sp(ispan_3),2)

      i_a=global_ist1_sp(1,ispcr_1)
c      mja=global_ist1_sp(2,ispcr_1)
c      mta=global_ist1_sp(3,ispcr_1)

c      ja=global_ist1(3,i_a)

      i_b=global_ist1_sp(1,ispcr_2)
c      mjb=global_ist1_sp(2,ispcr_2)
c      mtb=global_ist1_sp(3,ispcr_2)

c      jb=global_ist1(3,i_b)

      i_c=global_ist1_sp(1,ispcr_3)
c      mjc=global_ist1_sp(2,ispcr_3)
c      mtc=global_ist1_sp(3,ispcr_3)

c      jc=global_ist1(3,i_c)

      i_d=global_ist1_sp(1,ispan_1)
c      mjd=global_ist1_sp(2,ispan_1)
c      mtd=global_ist1_sp(3,ispan_1)

c      jd=global_ist1(3,i_d)

      i_e=global_ist1_sp(1,ispan_2)
c      mje=global_ist1_sp(2,ispan_2)
c      mte=global_ist1_sp(3,ispan_2)

c      je=global_ist1(3,i_e)

      i_f=global_ist1_sp(1,ispan_3)
c      mjf=global_ist1_sp(2,ispan_3)
c      mtf=global_ist1_sp(3,ispan_3)

c      jf=global_ist1(3,i_f)

c      mjtot2=mjd+mje+mjf
c      mttot2=mtd+mte+mtf

      iphase=1
      if (i_b>i_a) then
         itemp=ispcr_2
         itemp_i=i_b
         ispcr_2=ispcr_1
         i_b=i_a
         ispcr_1=itemp
         i_a=itemp_i
         iphase=-iphase
      endif
      if (i_c>i_b) then
         itemp=ispcr_3
         itemp_i=i_c
         ispcr_3=ispcr_2
         i_c=i_b
         ispcr_2=itemp
         i_b=itemp_i
         iphase=-iphase
         if (i_b>i_a) then
            itemp=ispcr_2
            itemp_i=i_b
            ispcr_2=ispcr_1
            i_b=i_a
            ispcr_1=itemp
            i_a=itemp_i
            iphase=-iphase
         endif
      endif
      if (i_e>i_d) then
         itemp=ispan_2
         itemp_i=i_e
         ispan_2=ispan_1
         i_e=i_d
         ispan_1=itemp
         i_d=itemp_i
         iphase=-iphase
      endif
      if (i_f>i_e) then
         itemp=ispan_3
         itemp_i=i_f
         ispan_3=ispan_2
         i_f=i_e
         ispan_2=itemp
         i_e=itemp_i
         iphase=-iphase
         if (i_e>i_d) then
            itemp=ispan_2
            itemp_i=i_e
            ispan_2=ispan_1
            i_e=i_d
            ispan_1=itemp
            i_d=itemp_i
            iphase=-iphase
         endif
      endif

c      i_a=global_ist1_sp(1,ispcr_1)
      mja=global_ist1_sp(2,ispcr_1)
      mta=global_ist1_sp(3,ispcr_1)

      ja=global_ist1(3,i_a)

c      i_b=global_ist1_sp(1,ispcr_2)
      mjb=global_ist1_sp(2,ispcr_2)
      mtb=global_ist1_sp(3,ispcr_2)

      jb=global_ist1(3,i_b)

c      i_c=global_ist1_sp(1,ispcr_3)
      mjc=global_ist1_sp(2,ispcr_3)
      mtc=global_ist1_sp(3,ispcr_3)

      jc=global_ist1(3,i_c)

c      i_d=global_ist1_sp(1,ispan_1)
      mjd=global_ist1_sp(2,ispan_1)
      mtd=global_ist1_sp(3,ispan_1)

      jd=global_ist1(3,i_d)

c      i_e=global_ist1_sp(1,ispan_2)
      mje=global_ist1_sp(2,ispan_2)
      mte=global_ist1_sp(3,ispan_2)

      je=global_ist1(3,i_e)

c      i_f=global_ist1_sp(1,ispan_3)
      mjf=global_ist1_sp(2,ispan_3)
      mtf=global_ist1_sp(3,ispan_3)

      jf=global_ist1(3,i_f)

      mjtot2=mjd+mje+mjf
      mttot2=mtd+mte+mtf

      i_abc=index_abc(i_a,i_b,i_c)
      i_def=index_abc(i_d,i_e,i_f)

      if (i_abc==0.or.i_def==0) then
         print *,' i_abc=',i_abc
         print *,' i_a,i_b,i_c=',i_a,i_b,i_c
         print *,' na,nb,nc=',global_ist1(1,i_a),global_ist1(1,i_b),
     $        global_ist1(1,i_c)
         print *,' la,lb,lc=',global_ist1(2,i_a),global_ist1(2,i_b),
     $        global_ist1(2,i_c)
         print *,' ja,jb,jc=',ja,jb,jc
         print *,' mja,mjb,mjc=',mja,mjb,mjc
         print *,' mta,mtb,mtc=',mta,mtb,mtc

         print *,' i_def=',i_def
         print *,' i_d,i_e,i_f=',i_d,i_e,i_f
         print *,' nd,ne,nf=',global_ist1(1,i_d),global_ist1(1,i_e),
     $        global_ist1(1,i_f)
         print *,' le,ld,lf=',global_ist1(2,i_d),global_ist1(2,i_e),
     $        global_ist1(2,i_f)
         print *,' jd,je,jf=',jd,je,jf
         print *,' mjd,mje,mjf=',mjd,mje,mjf
         print *,' mtd,mte,mtf=',mtd,mte,mtf
      endif

      if (i_def>i_abc) then
c         itemp_i=i_a
c         i_a=i_d
c         i_d=itemp_i
c         itemp_i=i_b
c         i_b=i_e
c         i_e=itemp_i
c         itemp_i=i_c
c         i_c=i_f
c         i_f=itemp_i
         j2_a=jd !j2_sp(ispan_1)
         j2_b=je !j2_sp(ispan_2)
         j2_c=jf !j2_sp(ispan_3)
         m2_a=mjd !m2_sp(ispan_1)
         m2_b=mje !m2_sp(ispan_2)
         m2_c=mjf !m2_sp(ispan_3)
         mt2_a=mtd !mt2_sp(ispan_1)
         mt2_b=mte !mt2_sp(ispan_2)
         mt2_c=mtf !mt2_sp(ispan_3)
         j2_d=ja !j2_sp(ispcr_1)
         j2_e=jb !j2_sp(ispcr_2)
         j2_f=jc !j2_sp(ispcr_3)
         m2_d=mja !m2_sp(ispcr_1)
         m2_e=mjb !m2_sp(ispcr_2)
         m2_f=mjc !m2_sp(ispcr_3)
         mt2_d=mta !mt2_sp(ispcr_1)
         mt2_e=mtb !mt2_sp(ispcr_2)
         mt2_f=mtc !mt2_sp(ispcr_3)
         i_abcdef=index_abcdef(i_def,i_abc)
      else
         j2_d=jd !j2_sp(ispan_1)
         j2_e=je !j2_sp(ispan_2)
         j2_f=jf !j2_sp(ispan_3)
         m2_d=mjd !m2_sp(ispan_1)
         m2_e=mje !m2_sp(ispan_2)
         m2_f=mjf !m2_sp(ispan_3)
         mt2_d=mtd !mt2_sp(ispan_1)
         mt2_e=mte !mt2_sp(ispan_2)
         mt2_f=mtf !mt2_sp(ispan_3)
         j2_a=ja !j2_sp(ispcr_1)
         j2_b=jb !j2_sp(ispcr_2)
         j2_c=jc !j2_sp(ispcr_3)
         m2_a=mja !m2_sp(ispcr_1)
         m2_b=mjb !m2_sp(ispcr_2)
         m2_c=mjc !m2_sp(ispcr_3)
         mt2_a=mta !mt2_sp(ispcr_1)
         mt2_b=mtb !mt2_sp(ispcr_2)
         mt2_c=mtc !mt2_sp(ispcr_3)
         i_abcdef=index_abcdef(i_abc,i_def)
      endif
c      print *,' i_abcdef=',i_abcdef

      cgt1(0:1)=cgt12(0:1,(mt2_a+1)/2,(mt2_b+1)/2)
      cgt2(0:1)=cgt12(0:1,(mt2_d+1)/2,(mt2_e+1)/2)

      cgt3(0:1,0:1)=cgt123(0:1,0:1,(mt2_a+mt2_b)/2,(mt2_c+1)/2)
      cgt4(0:1,0:1)=cgt123(0:1,0:1,(mt2_d+mt2_e)/2,(mt2_f+1)/2)

      sum3=0.d0
      iii=start_abcdef(i_abcdef)
      do j_ab=abs(j2_a-j2_b),(j2_a+j2_b),2
c         if (abs(m2_a+m2_b)>j_ab) then
c            cg1=0.d0
c         else
c            cg1=clebd(j2_a,m2_a,j2_b,m2_b,j_ab,m2_a+m2_b)
            cg1=cgj12(j_ab/2,j2_a/2,(j2_a+m2_a)/2,j2_b/2,(j2_b+m2_b)/2)
c         endif
         do j_de=abs(j2_d-j2_e),(j2_d+j2_e),2
c            if (abs(m2_d+m2_e)>j_de.or.cg1==0.d0) then
c               cg2=0.d0
c            else
c               cg2=clebd(j2_d,m2_d,j2_e,m2_e,j_de,m2_d+m2_e)
               cg2=cgj12(j_de/2,j2_d/2,(j2_d+m2_d)/2,
     $              j2_e/2,(j2_e+m2_e)/2)
c            endif
            do J_3=max(abs(j_ab-j2_c),
     $           abs(j_de-j2_f)),
     $           min(j_ab+j2_c,j_de+j2_f),2
c               if (abs(mjtot2)>J_3.or.cg1==0.d0.or.cg2==0.d0) then
c                  cg3=0.d0
c                  cg4=0.d0
c               else
c                  cg3=clebd(j_ab,m2_a+m2_b,j2_c,m2_c,J_3,mjtot2)
c                  cg4=clebd(j_de,m2_d+m2_e,j2_f,m2_f,J_3,mjtot2)
                  cg3=cgj123(J_3/2,j_ab/2,(m2_a+m2_b)/2,
     $                 j2_c/2,(j2_c+m2_c)/2)
                  cg4=cgj123(J_3/2,j_de/2,(m2_d+m2_e)/2,
     $                 j2_f/2,(j2_f+m2_f)/2)
c               endif

               sum3=sum3+cg1*cg2*cg3*cg4*(
     $              v3b_cJ(iii)*cgt1(0)*cgt2(0)*cgt3(0,0)*cgt4(0,0)
     $              +v3b_cJ(iii+1)*cgt1(0)*cgt2(1)*cgt3(0,0)*cgt4(1,0)
     $              +v3b_cJ(iii+2)*cgt1(1)*cgt2(0)*cgt3(1,0)*cgt4(0,0)
     $              +v3b_cJ(iii+3)*cgt1(1)*cgt2(1)*cgt3(1,0)*cgt4(1,0)
     $              +v3b_cJ(iii+4)*cgt1(1)*cgt2(1)*cgt3(1,1)*cgt4(1,1))

               iii=iii+5

            end do
         end do
      end do
      v3b_cJ_unc=sum3*real(iphase,kind(0.d0))
      end

      subroutine NCSMC_coupling_kernels
      use obdens
      use initst
      use finast
      use intrface
      use paramdef
      use occ
      use hash_tables
      use nodeinfo
      use kernels
      use v3b, only: V3Nint
      implicit none
      include 'mpif.h'
      integer :: j,count,ip,t,k_i,k_f
      integer,allocatable :: occi(:),occf(:),
     $     occim2(:),occim3(:)
      logical :: prot=.false.,neut=.false.
      integer :: nuc_min,nuc_max,i1,i2,p_state,n_state,ii,N_HO_i,
     $     ia_1,ia_2,ia_3,an_st_1,an_st_2,an_st_3,phase,ij,jj,ia_4,
     $     cr_st_min_1,cr_st_max_1,cr_st_min_2,cr_st_max_2,cr_1,an_st_4,
     $     cr_st_1,iphase,ia,mxsps,mt,mj,pi,mt_ind,mj_ind,N_HO_f,cr_st_2
      integer :: ibit,wd,interm_energy,ierr,dimalloc,dimallocmin,
     $     I_ab,t_ab,I_abc,t_abc,numel
      integer(8),allocatable :: intbas(:),intbasim3(:),intbasi(:),
     $     intbasim4(:)
      integer,allocatable :: locz(:)
      real(kind(0.d0)),allocatable :: tmp(:),tmp_allred(:)
      
      if (iproc==0) then
         print *,' NCSMC_kernels_calc entered'
      endif

      if (iproc==0) then
         write(iunitvout,*)
         write(iunitvout,7768) ist1dim
 7768    format(' number of single-nucleon states =',i4)
         do ii=1,ist1dim
            write(iunitvout,7770) ii,(ist1(ia,ii),ia=1,3)
 7770       format(' #',i4,'  n=',i3,'  l=',i3,'  j=',i2,'/2')
         end do
      endif   

      if (nprotonsi==nprotonsf) then
         neut=.true.
      elseif (nneutrnsi==nneutrnsf) then
         prot=.true.
      endif

c      if (iproc==0) print *,' jt2i=',jt2i(ki:ki+nki-1)
c      if (iproc==0) print *,' it2i=',it2i(ki:ki+nki-1)
      jmi=minval(jt2i(ki:ki+nki-1))
      jma=maxval(jt2i(ki:ki+nki-1))
      tmi=minval(it2i(ki:ki+nki-1))
      tma=maxval(it2i(ki:ki+nki-1))
      ipamin=(-1)**iparityi
      ipamax=(-1)**iparityi
c      if (iproc==0) print *,' jmi,jma,ipamin,ipamax,tmi,tma=',
c     $     jmi,jma,ipamin,ipamax,tmi,tma
      allocate(JpiT(jmi/2:jma/2,ipamin:ipamax,tmi/2:tma/2))
c      if (iproc==0) print *,' JpiT allocated'

      allocate(map_ji(jmi/2:jma/2,tmi/2:tma/2,ki:ki+nki-1))
c      if (iproc==0) print *,' map_ji allocated'
      map_ji=0
      do j=jmi,jma,2
         do t=tmi,tma,2
            count=0
            do k_i=ki,ki+nki-1   
               if (j==jt2i(k_i).and.t==it2i(k_i)) then
                  count=count+1
                  map_ji(j/2,t/2,k_i)=count
               endif
            end do
            JpiT(j/2,:,t/2)%dim_i=count
            JpiT(j/2,:,t/2)%dim_f=nkf
            do ip=ipamin,ipamax
c               if (iproc==0) print *,' j,ip,t,count=',j,ip,t,count
               allocate(JpiT(j/2,ip,t/2)%st_fi(nkf,count))
c               if (iproc==0) print *,' JpiT%st_fi allocated'
            end do
         end do
      end do

      select case(nucleonsi-nucleonsf)
      case(1)
         dimalloc=global_ist1_dim
         dimallocmin=1
      case(2)
         dimalloc=ist2_SDkern_dim
         dimallocmin=1
      case(3)
         dimalloc=ist3_SDkern_dim
         dimallocmin=1
      case default
         if (iproc==0) then
            write(iunitvout,*) 
     +           ' Not implemented for nucleonsi and nucleonsf',
     +           nucleonsi,nucleonsf
         endif
         stop
      end select

      do j=jmi,jma,2
         do ip=ipamin,ipamax
            do t=tmi,tma,2
               do k_i=1,JpiT(j/2,ip,t/2)%dim_i
                  do k_f=1,JpiT(j/2,ip,t/2)%dim_f
                     select case(nucleonsi-nucleonsf)
                     case(2)
                        dimallocmin=ist2_SDkern_dim+1
                        dimalloc=0
                        do ii=1,ist2_SDkern_dim
                           I_ab=ist2_SDkern(3,ii)
                           t_ab=ist2_SDkern(4,ii)
                           if (I_ab>=abs(j-jt2f(k_f))/2
     $                          .and.I_ab<=(j+jt2f(k_f))/2
     $                          .and.t_ab>=abs(t-it2f(k_f))/2
     $                          .and.t_ab<=(t+it2f(k_f))/2) then
                              dimallocmin=min(dimallocmin,ii)
                              dimalloc=max(dimalloc,ii)
                           endif
                        end do
                        if (iproc==0) then
                           print *,' j,ip,t,k_i,k_f,min,max=',
     $                          j,ip,t,k_i,k_f,dimallocmin,dimalloc
                        endif
                     case(3)
                        dimallocmin=ist3_SDkern_dim+1
                        dimalloc=0
                        do ii=1,ist3_SDkern_dim
                           I_abc=ist3_SDkern(3,ii)
                           t_abc=ist3_SDkern(4,ii)
                           if (I_abc>=abs(j-jt2f(k_f))
     $                          .and.I_abc<=(j+jt2f(k_f))
     $                          .and.t_abc>=abs(t-it2f(k_f))
     $                          .and.t_abc<=(t+it2f(k_f))) then
                              dimallocmin=min(dimallocmin,ii)
                              dimalloc=max(dimalloc,ii)
                           endif
                        end do
                        if (iproc==0) then
                           print *,' j,ip,t,k_i,k_f,min,max=',
     $                          j,ip,t,k_i,k_f,dimallocmin,dimalloc
                        endif
                     end select
                     JpiT(j/2,ip,t/2)%st_fi(k_f,k_i)%dim_1=
     $                    dimalloc
                     JpiT(j/2,ip,t/2)%st_fi(k_f,k_i)%start_1=
     $                    dimallocmin
                     allocate(JpiT(j/2,ip,t/2)%st_fi(k_f,k_i)
     $                    %g(dimallocmin:dimalloc))     
                     allocate(JpiT(j/2,ip,t/2)%st_fi(k_f,k_i)
     $                    %h_NN(dimallocmin:dimalloc))
                     JpiT(j/2,ip,t/2)%st_fi(k_f,k_i)%g=0.d0
                     JpiT(j/2,ip,t/2)%st_fi(k_f,k_i)%h_NN=0.d0
                     if (nucleonsi-nucleonsf==1.or.V3Nint) then
                        allocate(JpiT(j/2,ip,t/2)%st_fi(k_f,k_i)
     $                       %h_3N(dimallocmin:dimalloc))
                        JpiT(j/2,ip,t/2)%st_fi(k_f,k_i)%h_3N=0.d0
                     endif
                  end do
               end do
            end do
         end do
      end do
      if (iproc==0) then
         print *,' ncsmc coupling kernels allocated'
      endif

      allocate(occi(nucleonsi))
      allocate(occf(nucleonsf))
      if (nucleonsi>1) then
         allocate(occim2(nucleonsi-2))
      endif
      if (nucleonsi>2) then
         allocate(occim3(nucleonsi-3))
      endif

      call mult_init
c      if (iproc==0) print *,' mult_init called'

      mxsps=nbit*mxnwdi
      allocate(locz(2*mxsps))
      allocate(intbas(2*mxnwdi))
      allocate(intbasi(2*mxnwdi))
      allocate(intbasim3(2*mxnwdi))
      allocate(intbasim4(2*mxnwdi))

      select case(nucleonsi-nucleonsf)
      case(1)

         if (prot) then
            nuc_min=1
            nuc_max=nprotonsi
         else
            nuc_min=nprotonsi+1
            nuc_max=nucleonsi
         endif

         main_i1_loop_a_1: do i1=iproc+1,nsdi,nproc

            if (i1+nproc-iproc-1==10000*((i1+nproc-iproc-1)/10000))
     +           then
               print *, '#',iproc,' doing i1=',i1
            endif
         
            if (majortot==2) then
               occi(:)=iloci(:,i1)
            elseif (majortot==3) then
               p_state=I_state_i(1,i1)
               n_state=I_state_i(2,i1)
               occi(1:nprotonsi)=occ_p_i(:,p_state)
               do ii=1,nneutrnsi
                  occi(nprotonsi+ii)=occ_n_i(ii,n_state)+mxnwdi*nbit
               end do
            endif

c            if (iproc==0) print *,' i1,occi set:',i1

            intbasi=0
            do cr_1=1,nucleonsi
               wd=(occi(cr_1)-1)/nbit+1
               ibit=mod(occi(cr_1)-1,nbit)
               intbasi(wd)=ibset(intbasi(wd),ibit)
            end do

            N_HO_i=0
            do ia_1=1,nucleonsi
               an_st_1=occi(ia_1)
               N_HO_i=N_HO_i+2*n_spi(an_st_1)+l_spi(an_st_1)
            end do

c            if (iproc==0) print *,' i1,intbasi,N_HO_i set:',i1

            do ia_1=nuc_min,nuc_max
               an_st_1=occi(ia_1)
               if (m2_spi(an_st_1)+mjtotalf/=mjtotali) cycle
               if (mod(iparityf+l_spi(an_st_1)+iparityi,2)/=0) cycle
               if (N_HO_i-2*n_spi(an_st_1)-l_spi(an_st_1)>nhwf) cycle

               phase=(-1)**(nucleonsi-ia_1)
               occf(1:ia_1-1)=occi(1:ia_1-1)
               occf(ia_1:nucleonsf)=occi(ia_1+1:nucleonsi)
               call get_state_index(nucleonsf,occf,1,i2)
               if (i2==-1) cycle
c               if (iproc==0) print *,' i1,ia_1,i2=',i1,ia_1,i2
               call g_a_1(an_st_1)
c               if (iproc==0) print *,' i1,g_a_1 called:',i1

            end do

            do ij=1,mult_nucl2_dim
               ia_1=mult_nucl2(1,ij)
               ia_2=mult_nucl2(2,ij)
               an_st_1=occi(ia_1)
               an_st_2=occi(ia_2)

               select case(mttotalf-mttotali+mt2_spi(an_st_1)
     $              +mt2_spi(an_st_2))
               case(1)
                  cr_st_min_1=1
                  cr_st_max_1=naspsi
               case(-1)
                  cr_st_min_1=mxsps+1
                  cr_st_max_1=mxsps+naspsi
               case default
                  cycle
               end select
               if (nucleonsi>2) then
                  occim2(1:ia_1-1)=occi(1:ia_1-1)
                  occim2(ia_1:ia_2-2)=occi(ia_1+1:ia_2-1)
                  occim2(ia_2-1:nucleonsi-2)=occi(ia_2+1:
     $                 nucleonsi)
                  
                  locz(1:occim2(1))=0
                  do cr_1=2,nucleonsi-2
                     locz(occim2(cr_1-1)+1:occim2(cr_1))=cr_1-1
                  end do
                  locz(occim2(nucleonsi-2)+1:2*mxsps)=nucleonsi-2
                  
                  intbas=0
                  do cr_1=1,nucleonsi-2
                     wd=(occim2(cr_1)-1)/nbit+1
                     ibit=mod(occim2(cr_1)-1,nbit)
                     intbas(wd)=ibset(intbas(wd),ibit)
                  end do
               else
                  intbas=0
               endif
               interm_energy=-2*n_spi(an_st_1)-l_spi(an_st_1)
     $              -2*n_spi(an_st_2)-l_spi(an_st_2)+N_HO_i
               do cr_st_1=cr_st_min_1,cr_st_max_1
                  
                  if (2*n_spi(cr_st_1)+l_spi(cr_st_1)
     +                 +interm_energy>nhwf) exit
                  if (m2_spi(cr_st_1)
     $                 -m2_spi(an_st_1)-m2_spi(an_st_2)+mjtotali
     +                 /=mjtotalf) cycle
                  if (mod(iparityf+l_spi(cr_st_1)+l_spi(an_st_1)
     $                 +l_spi(an_st_2)+iparityi,2)/=0) cycle
                  
                  wd=(cr_st_1-1)/nbit+1
                  ibit=mod(cr_st_1-1,nbit)
                  if (btest(intbas(wd),ibit)) cycle
                  if (nucleonsi>2) then
                     occf(1:locz(cr_st_1))=occim2(1:locz(cr_st_1))
                     occf(locz(cr_st_1)+1)=cr_st_1
                     occf(locz(cr_st_1)+2:nucleonsf)=
     $                    occim2(locz(cr_st_1)+1:nucleonsi-2)
                  else
                     occf(1)=cr_st_1
                  endif
                  
                  call get_state_index(nucleonsf,occf,1,i2)
                  
c                  if (iproc==0) print *,' i1,ia_1,ia_2,cr_st_1,i2=',
c     $                 i1,ia_1,ia_2,cr_st_1,i2

                  if (i2==-1) cycle
                  if (nucleonsi>2) then
                     phase=(-1)**(ia_1+ia_2
     $                    +nucleonsi-locz(cr_st_1)) 
                  else
                     phase=(-1)**(ia_1+ia_2) 
                  endif
                  
                  call h_NN_a_1(cr_st_1,an_st_2,
     $                 an_st_1)
c                  if (iproc==0) print *,' i1,h_NN_a_1 called:',i1
                  
               end do
            end do

            if (.not.V3Nint) cycle

            locz(1:occi(1))=0
            do cr_1=2,nucleonsi
               locz(occi(cr_1-1)+1:occi(cr_1))=cr_1-1
            end do
            locz(occi(nucleonsi)+1:2*mxsps)=nucleonsi

            do ij=1,mult_nucl3_dim

               ia_1=mult_nucl3(1,ij)
               ia_2=mult_nucl3(2,ij)
               ia_3=mult_nucl3(3,ij)

               an_st_1=occi(ia_1)
               an_st_2=occi(ia_2)
               an_st_3=occi(ia_3)

               intbasim3=intbasi
               wd=(an_st_1-1)/nbit+1
               ibit=mod(an_st_1-1,nbit)
               intbasim3(wd)=ibclr(intbasim3(wd),ibit)
               wd=(an_st_2-1)/nbit+1
               ibit=mod(an_st_2-1,nbit)
               intbasim3(wd)=ibclr(intbasim3(wd),ibit)
               wd=(an_st_3-1)/nbit+1
               ibit=mod(an_st_3-1,nbit)
               intbasim3(wd)=ibclr(intbasim3(wd),ibit)

               mt=mttotalf-mttotali+mt2_spi(an_st_1)
     $              +mt2_spi(an_st_2)+mt2_spi(an_st_3)
               mj=mjtotalf-mjtotali+m2_spi(an_st_1)
     $              +m2_spi(an_st_2)+m2_spi(an_st_3)
               pi=mod(iparityf+iparityi+l_spi(an_st_1)+l_spi(an_st_2)
     $              +l_spi(an_st_3),2)

               if (abs(mt)>2) cycle

c               if (abs(mt)>2) then
c                  print *,' iproc,i1,mt,mttotalf,mttotali=',
c     $                 iproc,i1,mt,mttotalf,mttotali
c                  print *,' an_st_1,an_st_2,an_st_3=',
c     $                  an_st_1,an_st_2,an_st_3
c                  print *,' mt2:an_st1,an_st_2,an_st_3=',
c     $                 mt2_spi(an_st_1),
c     $                 mt2_spi(an_st_2),mt2_spi(an_st_3) 
c                  print *,' mj,pi=',mj,pi
c                  stop
c               endif

               mt_ind=(mt+2)/2
               mj_ind=(mshift2(abs(mt)/2)+mj)/2

               interm_energy=-2*n_spi(an_st_1)-l_spi(an_st_1)
     $              -2*n_spi(an_st_2)-l_spi(an_st_2)
     $              -2*n_spi(an_st_3)-l_spi(an_st_3)+N_HO_i

               if (mj_ind>maxval(mshift2).or.mj_ind<0) then
                  print *,' iproc: pi,mt,mj=',iproc,pi,mt,mj
                  print *,' n_spi,l_spi:',n_spi(an_st_1),l_spi(an_st_1),
     $                 n_spi(an_st_2),l_spi(an_st_2),
     $                 n_spi(an_st_3),l_spi(an_st_3)
                  print *,' interm_energy=',interm_energy
                  cycle
               endif
               
               jj=multiple2_point(pi,mt_ind,mj_ind)

               do ii=1,multiple2_dim(pi,mt_ind,mj_ind)

                  jj=jj+1
                  cr_st_1=multiple2_ist(1,jj)
                  cr_st_2=multiple2_ist(2,jj)

                  N_HO_f=2*n_spi(cr_st_1)+l_spi(cr_st_1)
     $                 +2*n_spi(cr_st_2)+l_spi(cr_st_2)
     +                 +interm_energy

                  if (N_HO_f>nhwf) exit

                  intbas=intbasim3
                  wd=(cr_st_1-1)/nbit+1
                  ibit=mod(cr_st_1-1,nbit)
                  if (btest(intbas(wd),ibit)) cycle
                  intbas(wd)=ibset(intbas(wd),ibit)
                  wd=(cr_st_2-1)/nbit+1
                  ibit=mod(cr_st_2-1,nbit)
                  if (btest(intbas(wd),ibit)) cycle
                  intbas(wd)=ibset(intbas(wd),ibit)

                  call get_state_index2(nucleonsf,2*mxnwdi,intbas,i2)
                  if (i2==-1) cycle

                  iphase=ia_1+ia_2+ia_3+nucleonsi+1
     $                 +locz(cr_st_1)+locz(cr_st_2)
                  if (cr_st_1>an_st_1) iphase=iphase+1
                  if (cr_st_1>an_st_2) iphase=iphase+1
                  if (cr_st_1>an_st_3) iphase=iphase+1
                  if (cr_st_2>an_st_1) iphase=iphase+1
                  if (cr_st_2>an_st_2) iphase=iphase+1
                  if (cr_st_2>an_st_3) iphase=iphase+1
                  phase=(-1)**iphase

                  call h_3N_a_1(cr_st_1,cr_st_2,an_st_3,an_st_2,
     $                 an_st_1)
c                  if (iproc==0) print *,' i1,h_3N_a_1 called:',i1

               end do
            end do

         end do main_i1_loop_a_1

      case(2)

         main_i1_loop_a_2: do i1=iproc+1,nsdi,nproc

            if (i1+nproc-iproc-1==1000000*((i1+nproc-iproc-1)/1000000))
     +           then
               print *, '#',iproc,' doing i1=',i1
            endif
         
            if (majortot==2) then
               occi(:)=iloci(:,i1)
            elseif (majortot==3) then
               p_state=I_state_i(1,i1)
               n_state=I_state_i(2,i1)
               occi(1:nprotonsi)=occ_p_i(:,p_state)
               do ii=1,nneutrnsi
                  occi(nprotonsi+ii)=occ_n_i(ii,n_state)+mxnwdi*nbit
               end do
            endif

c            if (iproc==0) print *,' i1,occi set:',i1

            intbasi=0
            do cr_1=1,nucleonsi
               wd=(occi(cr_1)-1)/nbit+1
               ibit=mod(occi(cr_1)-1,nbit)
               intbasi(wd)=ibset(intbasi(wd),ibit)
            end do

            N_HO_i=0
            do ia_1=1,nucleonsi
               an_st_1=occi(ia_1)
               N_HO_i=N_HO_i+2*n_spi(an_st_1)+l_spi(an_st_1)
            end do

c            if (iproc==0) print *,' i1,intbasi,N_HO_i set:',i1

            do ij=1,mult_nucl2_dim
               ia_1=mult_nucl2(1,ij)
               ia_2=mult_nucl2(2,ij)
               an_st_1=occi(ia_1)
               an_st_2=occi(ia_2)

               if (mjtotalf-mjtotali+m2_spi(an_st_1)
     $              +m2_spi(an_st_2)/=0) cycle
               if (mttotalf-mttotali+mt2_spi(an_st_1)
     $              +mt2_spi(an_st_2)/=0) cycle
               if (mod(iparityf+iparityi+l_spi(an_st_1)+l_spi(an_st_2)
     $              ,2)==1)
     $              cycle

               occim2(1:ia_1-1)=occi(1:ia_1-1)
               occim2(ia_1:ia_2-2)=occi(ia_1+1:ia_2-1)
               occim2(ia_2-1:nucleonsi-2)=occi(ia_2+1:
     $              nucleonsi)

               call get_state_index(nucleonsf,occim2,1,i2)

               if (i2==-1) cycle
               phase=(-1)**(ia_1+ia_2)

               call g_a_2(an_st_2,an_st_1) 

            end do

            locz(1:occi(1))=0
            do cr_1=2,nucleonsi
               locz(occi(cr_1-1)+1:occi(cr_1))=cr_1-1
            end do
            locz(occi(nucleonsi)+1:2*mxsps)=nucleonsi

            if (nucleonsi<3) cycle

            do ij=1,mult_nucl3_dim

               ia_1=mult_nucl3(1,ij)
               ia_2=mult_nucl3(2,ij)
               ia_3=mult_nucl3(3,ij)

               an_st_1=occi(ia_1)
               an_st_2=occi(ia_2)
               an_st_3=occi(ia_3)

               intbasim3=intbasi
               wd=(an_st_1-1)/nbit+1
               ibit=mod(an_st_1-1,nbit)
               intbasim3(wd)=ibclr(intbasim3(wd),ibit)
               wd=(an_st_2-1)/nbit+1
               ibit=mod(an_st_2-1,nbit)
               intbasim3(wd)=ibclr(intbasim3(wd),ibit)
               wd=(an_st_3-1)/nbit+1
               ibit=mod(an_st_3-1,nbit)
               intbasim3(wd)=ibclr(intbasim3(wd),ibit)

               mt=mttotalf-mttotali+mt2_spi(an_st_1)
     $              +mt2_spi(an_st_2)+mt2_spi(an_st_3)
               mj=mjtotalf-mjtotali+m2_spi(an_st_1)
     $              +m2_spi(an_st_2)+m2_spi(an_st_3)
               pi=mod(iparityf+iparityi+l_spi(an_st_1)+l_spi(an_st_2)
     $              +l_spi(an_st_3),2)

               interm_energy=-2*n_spi(an_st_1)-l_spi(an_st_1)
     $              -2*n_spi(an_st_2)-l_spi(an_st_2)
     $              -2*n_spi(an_st_3)-l_spi(an_st_3)+N_HO_i

               if (mt==1) then
                  cr_st_min_1=1
                  cr_st_max_1=naspsi
               elseif (mt==-1) then
                  cr_st_min_1=mxsps+1
                  cr_st_max_1=mxsps+naspsi
               else
                  cycle
c                  print *,'***error: mt=',mt
c                  stop
               endif

               do cr_st_1=cr_st_min_1,cr_st_max_1
                  if (2*n_spi(cr_st_1)+l_spi(cr_st_1)
     +                 +interm_energy>nhwf) exit
                  if (m2_spi(cr_st_1)/=mj) cycle
                  if (mod(l_spi(cr_st_1),2)/=pi) cycle

                  intbas=intbasim3
                  wd=(cr_st_1-1)/nbit+1
                  ibit=mod(cr_st_1-1,nbit)
                  if (btest(intbas(wd),ibit)) cycle
                  intbas(wd)=ibset(intbas(wd),ibit)

                  call get_state_index2(nucleonsf,2*mxnwdi,intbas,i2)
                  if (i2==-1) cycle

                  iphase=ia_1+ia_2+ia_3+1+locz(cr_st_1)
                  if (cr_st_1>an_st_1) iphase=iphase+1
                  if (cr_st_1>an_st_2) iphase=iphase+1
                  if (cr_st_1>an_st_3) iphase=iphase+1                  
                  phase=(-1)**iphase

                  call h_NN_a_2(cr_st_1,an_st_3,an_st_2,an_st_1,1)
                  call h_NN_a_2(cr_st_1,an_st_3,an_st_1,an_st_2,-1)
                  call h_NN_a_2(cr_st_1,an_st_2,an_st_1,an_st_3,1)
                  
                  if (V3Nint) then
                     call h_3N_a_2b(cr_st_1,an_st_3,an_st_2,an_st_1)
                  endif

               end do
            end do

            if (.not.V3Nint.or.nucleonsi<4) cycle

            do ij=1,mult_nucl4_dim

               ia_1=mult_nucl4(1,ij)
               ia_2=mult_nucl4(2,ij)
               ia_3=mult_nucl4(3,ij)
               ia_4=mult_nucl4(4,ij)

               an_st_1=occi(ia_1)
               an_st_2=occi(ia_2)
               an_st_3=occi(ia_3)
               an_st_4=occi(ia_4)

               intbasim4=intbasi
               wd=(an_st_1-1)/nbit+1
               ibit=mod(an_st_1-1,nbit)
               intbasim4(wd)=ibclr(intbasim4(wd),ibit)
               wd=(an_st_2-1)/nbit+1
               ibit=mod(an_st_2-1,nbit)
               intbasim4(wd)=ibclr(intbasim4(wd),ibit)
               wd=(an_st_3-1)/nbit+1
               ibit=mod(an_st_3-1,nbit)
               intbasim4(wd)=ibclr(intbasim4(wd),ibit)
               wd=(an_st_4-1)/nbit+1
               ibit=mod(an_st_4-1,nbit)
               intbasim4(wd)=ibclr(intbasim4(wd),ibit)

               mt=mttotalf-mttotali+mt2_spi(an_st_1)
     $              +mt2_spi(an_st_2)+mt2_spi(an_st_3)+mt2_spi(an_st_4)
               mj=mjtotalf-mjtotali+m2_spi(an_st_1)
     $              +m2_spi(an_st_2)+m2_spi(an_st_3)+m2_spi(an_st_4)
               pi=mod(iparityf+iparityi+l_spi(an_st_1)+l_spi(an_st_2)
     $              +l_spi(an_st_3)+l_spi(an_st_4),2)

               if (abs(mt)>2) cycle

               mt_ind=(mt+2)/2
               mj_ind=(mshift2(abs(mt)/2)+mj)/2

               interm_energy=-2*n_spi(an_st_1)-l_spi(an_st_1)
     $              -2*n_spi(an_st_2)-l_spi(an_st_2)
     $              -2*n_spi(an_st_3)-l_spi(an_st_3)
     $              -2*n_spi(an_st_4)-l_spi(an_st_4)+N_HO_i

               if (mj_ind>maxval(mshift2).or.mj_ind<0) then
                  print *,' iproc: pi,mt,mj=',iproc,pi,mt,mj
                  print *,' n_spi,l_spi:',n_spi(an_st_1),l_spi(an_st_1),
     $                 n_spi(an_st_2),l_spi(an_st_2),
     $                 n_spi(an_st_3),l_spi(an_st_3),
     $                 n_spi(an_st_4),l_spi(an_st_4)
                  print *,' interm_energy=',interm_energy
                  cycle
               endif

               jj=multiple2_point(pi,mt_ind,mj_ind)

               do ii=1,multiple2_dim(pi,mt_ind,mj_ind)

                  jj=jj+1
                  cr_st_1=multiple2_ist(1,jj)
                  cr_st_2=multiple2_ist(2,jj)

                  N_HO_f=2*n_spi(cr_st_1)+l_spi(cr_st_1)
     $                 +2*n_spi(cr_st_2)+l_spi(cr_st_2)
     +                 +interm_energy

                  if (N_HO_f>nhwf) exit

                  intbas=intbasim4
                  wd=(cr_st_1-1)/nbit+1
                  ibit=mod(cr_st_1-1,nbit)
                  if (btest(intbas(wd),ibit)) cycle
                  intbas(wd)=ibset(intbas(wd),ibit)
                  wd=(cr_st_2-1)/nbit+1
                  ibit=mod(cr_st_2-1,nbit)
                  if (btest(intbas(wd),ibit)) cycle
                  intbas(wd)=ibset(intbas(wd),ibit)

                  call get_state_index2(nucleonsf,2*mxnwdi,intbas,i2)
                  if (i2==-1) cycle
                  
                  iphase=ia_1+ia_2+ia_3+ia_4+1
     $                 +locz(cr_st_1)+locz(cr_st_2)
                  if (cr_st_1>an_st_1) iphase=iphase+1
                  if (cr_st_1>an_st_2) iphase=iphase+1
                  if (cr_st_1>an_st_3) iphase=iphase+1
                  if (cr_st_1>an_st_4) iphase=iphase+1
                  if (cr_st_2>an_st_1) iphase=iphase+1
                  if (cr_st_2>an_st_2) iphase=iphase+1
                  if (cr_st_2>an_st_3) iphase=iphase+1
                  if (cr_st_2>an_st_4) iphase=iphase+1
                  phase=(-1)**iphase

                  call h_3N_a_2a(cr_st_1,cr_st_2,
     $                 an_st_4,an_st_3,an_st_2,an_st_1,1)
                  call h_3N_a_2a(cr_st_1,cr_st_2,
     $                 an_st_4,an_st_3,an_st_1,an_st_2,-1)
                  call h_3N_a_2a(cr_st_1,cr_st_2,
     $                 an_st_4,an_st_2,an_st_1,an_st_3,1)
                  call h_3N_a_2a(cr_st_1,cr_st_2,
     $                 an_st_3,an_st_2,an_st_1,an_st_4,-1)

               end do
            end do

         end do main_i1_loop_a_2

      case(3)

         main_i1_loop_a_3: do i1=iproc+1,nsdi,nproc

            if (i1+nproc-iproc-1==100000*((i1+nproc-iproc-1)/100000))
c            if (i1+nproc-iproc-1==1000000*((i1+nproc-iproc-1)/1000000))
     +           then
               print *, '#',iproc,' doing i1=',i1
            endif
         
            if (majortot==2) then
               occi(:)=iloci(:,i1)
            elseif (majortot==3) then
               p_state=I_state_i(1,i1)
               n_state=I_state_i(2,i1)
               occi(1:nprotonsi)=occ_p_i(:,p_state)
               do ii=1,nneutrnsi
                  occi(nprotonsi+ii)=occ_n_i(ii,n_state)+mxnwdi*nbit
               end do
            endif

c            if (iproc==0) print *,' i1,occi set:',i1

            intbasi=0
            do cr_1=1,nucleonsi
               wd=(occi(cr_1)-1)/nbit+1
               ibit=mod(occi(cr_1)-1,nbit)
               intbasi(wd)=ibset(intbasi(wd),ibit)
            end do

            N_HO_i=0
            do ia_1=1,nucleonsi
               an_st_1=occi(ia_1)
               N_HO_i=N_HO_i+2*n_spi(an_st_1)+l_spi(an_st_1)
            end do

            do ij=1,mult_nucl3_dim

               ia_1=mult_nucl3(1,ij)
               ia_2=mult_nucl3(2,ij)
               ia_3=mult_nucl3(3,ij)

               an_st_1=occi(ia_1)
               an_st_2=occi(ia_2)
               an_st_3=occi(ia_3)

               if (mjtotalf-mjtotali+m2_spi(an_st_1)
     $              +m2_spi(an_st_2)+m2_spi(an_st_3)/=0) cycle
               if (mttotalf-mttotali+mt2_spi(an_st_1)
     $              +mt2_spi(an_st_2)+mt2_spi(an_st_3)/=0) cycle
               if (mod(iparityf+iparityi+l_spi(an_st_1)+l_spi(an_st_2)
     $              +l_spi(an_st_3),2)==1) cycle

               intbasim3=intbasi
               wd=(an_st_1-1)/nbit+1
               ibit=mod(an_st_1-1,nbit)
               intbasim3(wd)=ibclr(intbasim3(wd),ibit)
               wd=(an_st_2-1)/nbit+1
               ibit=mod(an_st_2-1,nbit)
               intbasim3(wd)=ibclr(intbasim3(wd),ibit)
               wd=(an_st_3-1)/nbit+1
               ibit=mod(an_st_3-1,nbit)
               intbasim3(wd)=ibclr(intbasim3(wd),ibit)

               call get_state_index2(nucleonsf,2*mxnwdi,intbasim3,i2)
               if (i2==-1) cycle
               phase=(-1)**(ia_1+ia_2+ia_3+nucleonsi)
               call g_a_3(an_st_3,an_st_2,an_st_1,1) 
               call g_a_3(an_st_2,an_st_3,an_st_1,-1)
               call g_a_3(an_st_1,an_st_3,an_st_2,1)  

            end do

            locz(1:occi(1))=0
            do cr_1=2,nucleonsi
               locz(occi(cr_1-1)+1:occi(cr_1))=cr_1-1
            end do
            locz(occi(nucleonsi)+1:2*mxsps)=nucleonsi

            do ij=1,mult_nucl4_dim

               ia_1=mult_nucl4(1,ij)
               ia_2=mult_nucl4(2,ij)
               ia_3=mult_nucl4(3,ij)
               ia_4=mult_nucl4(4,ij)

               an_st_1=occi(ia_1)
               an_st_2=occi(ia_2)
               an_st_3=occi(ia_3)
               an_st_4=occi(ia_4)

               intbasim4=intbasi
               wd=(an_st_1-1)/nbit+1
               ibit=mod(an_st_1-1,nbit)
               intbasim4(wd)=ibclr(intbasim4(wd),ibit)
               wd=(an_st_2-1)/nbit+1
               ibit=mod(an_st_2-1,nbit)
               intbasim4(wd)=ibclr(intbasim4(wd),ibit)
               wd=(an_st_3-1)/nbit+1
               ibit=mod(an_st_3-1,nbit)
               intbasim4(wd)=ibclr(intbasim4(wd),ibit)
               wd=(an_st_4-1)/nbit+1
               ibit=mod(an_st_4-1,nbit)
               intbasim4(wd)=ibclr(intbasim4(wd),ibit)

               mt=mttotalf-mttotali+mt2_spi(an_st_1)
     $              +mt2_spi(an_st_2)+mt2_spi(an_st_3)+mt2_spi(an_st_4)
               mj=mjtotalf-mjtotali+m2_spi(an_st_1)
     $              +m2_spi(an_st_2)+m2_spi(an_st_3)+m2_spi(an_st_4)
               pi=mod(iparityf+iparityi+l_spi(an_st_1)+l_spi(an_st_2)
     $              +l_spi(an_st_3)+l_spi(an_st_4),2)

               if (mt==1) then
                  cr_st_min_1=1
                  cr_st_max_1=naspsi
               elseif (mt==-1) then
                  cr_st_min_1=mxsps+1
                  cr_st_max_1=mxsps+naspsi
               else
                  cycle
c                  print *,'***error: mt=',mt
c                  stop
               endif

               interm_energy=-2*n_spi(an_st_1)-l_spi(an_st_1)
     $              -2*n_spi(an_st_2)-l_spi(an_st_2)
     $              -2*n_spi(an_st_3)-l_spi(an_st_3)
     $              -2*n_spi(an_st_4)-l_spi(an_st_4)+N_HO_i

               do cr_st_1=cr_st_min_1,cr_st_max_1
                  if (2*n_spi(cr_st_1)+l_spi(cr_st_1)
     +                 +interm_energy>nhwf) exit
                  if (m2_spi(cr_st_1)/=mj) cycle
                  if (mod(l_spi(cr_st_1),2)/=pi) cycle

                  intbas=intbasim4
                  wd=(cr_st_1-1)/nbit+1
                  ibit=mod(cr_st_1-1,nbit)
                  if (btest(intbas(wd),ibit)) cycle
                  intbas(wd)=ibset(intbas(wd),ibit)

                  call get_state_index2(nucleonsf,2*mxnwdi,intbas,i2)
                  if (i2==-1) cycle

                  iphase=ia_1+ia_2+ia_3+ia_4+nucleonsi-locz(cr_st_1)
                  if (cr_st_1>an_st_1) iphase=iphase+1
                  if (cr_st_1>an_st_2) iphase=iphase+1
                  if (cr_st_1>an_st_3) iphase=iphase+1  
                  if (cr_st_1>an_st_4) iphase=iphase+1                  
                  phase=(-1)**iphase

                  call h_NN_a_3(cr_st_1,an_st_4,an_st_3,an_st_2,an_st_1
     $                 ,1)
                  call h_NN_a_3(cr_st_1,an_st_4,an_st_2,an_st_1,an_st_3
     $                 ,1)
                  call h_NN_a_3(cr_st_1,an_st_2,an_st_1,an_st_4,an_st_3
     $                 ,1)
                  call h_NN_a_3(cr_st_1,an_st_4,an_st_1,an_st_2,an_st_3
     $                 ,-1)
                  call h_NN_a_3(cr_st_1,an_st_3,an_st_2,an_st_4,an_st_1
     $                 ,1)
                  call h_NN_a_3(cr_st_1,an_st_3,an_st_1,an_st_4,an_st_2
     $                 ,-1)

               end do
            end do

         end do main_i1_loop_a_3

      case default
         if (iproc==0) then
            write(iunitvout,*) 
     +           ' Not implemented for nucleonsi and nucleonsf',
     +           nucleonsi,nucleonsf
         endif
         stop
      end select

      call MPI_Barrier(icomm,ierr)
      if (nproc>1) then

         do j=jmi,jma,2
            do ip=ipamin,ipamax
               do t=tmi,tma,2
                  do k_i=1,JpiT(j/2,ip,t/2)%dim_i
                     do k_f=1,JpiT(j/2,ip,t/2)%dim_f
                        dimalloc=JpiT(j/2,ip,t/2)%st_fi(k_f,k_i)%dim_1
                        dimallocmin=JpiT(j/2,ip,t/2)
     $                       %st_fi(k_f,k_i)%start_1
                        numel=dimalloc-dimallocmin+1
                        allocate(tmp(numel))
                        allocate(tmp_allred(numel))
                        tmp(1:numel)=JpiT(j/2,ip,t/2)%st_fi(k_f,k_i)
     $                       %g(dimallocmin:dimalloc)
                        call MPI_Allreduce(tmp(1),tmp_allred(1),
     $                       numel,MPI_REAL8,MPI_SUM,
     $                       icomm,ierr)
                        JpiT(j/2,ip,t/2)%st_fi(k_f,k_i)
     $                       %g(dimallocmin:dimalloc)
     $                       =tmp_allred(1:numel)
                        tmp(1:numel)=JpiT(j/2,ip,t/2)%st_fi(k_f,k_i)
     $                       %h_NN(dimallocmin:dimalloc)
                        call MPI_Allreduce(tmp(1),tmp_allred(1),
     $                       numel,MPI_REAL8,MPI_SUM,
     $                       icomm,ierr)
                        JpiT(j/2,ip,t/2)%st_fi(k_f,k_i)
     $                       %h_NN(dimallocmin:dimalloc)
     $                       =tmp_allred(1:numel)
                        if (V3Nint) then
                           tmp(1:numel)=JpiT(j/2,ip,t/2)%st_fi(k_f,k_i)
     $                          %h_3N(dimallocmin:dimalloc)
                           call MPI_Allreduce(tmp(1),tmp_allred(1),
     $                          numel,MPI_REAL8,MPI_SUM,
     $                          icomm,ierr)
                           JpiT(j/2,ip,t/2)%st_fi(k_f,k_i)
     $                          %h_3N(dimallocmin:dimalloc)
     $                          =tmp_allred(1:numel)
                        endif
                        deallocate(tmp,tmp_allred)
                     end do
                  end do
               end do
            end do
         end do
      endif

      contains
      
      subroutine mult_init
      implicit  none
      integer :: ia_1,ia_2,ia_3,ij,nucleons,mult_nucl_dim
      integer :: mshift(0:1),mj,pi,mt,N_1,N_2,N_3,mj_ind,mt_ind,
     $     cr_st_1,cr_st_2,nasps,N12_max,N123,mxsps,
     $     cr_st_min_1,cr_st_min_2,cr_st_max_1,cr_st_max_2,
     $     num_of_mult,ia_4
      
      nucleons=nucleonsi

      mult_nucl_dim=0
      do ia_1=1,nucleons-1
         do ia_2=ia_1+1,nucleons
            mult_nucl_dim=mult_nucl_dim+1
         end do
      end do
      if (iproc==0) print *,' mult_nucl2_dim=',mult_nucl_dim
      mult_nucl2_dim=mult_nucl_dim
      allocate(mult_nucl2(2,mult_nucl_dim))
      ij=0
      do ia_1=1,nucleons-1
         do ia_2=ia_1+1,nucleons
            ij=ij+1
            mult_nucl2(1,ij)=ia_1
            mult_nucl2(2,ij)=ia_2
         end do
      end do
      if (ij/=mult_nucl_dim) then
         print *,'*** error: ii,mult_nucl_dim=',
     $        ij,mult_nucl_dim
         stop
      endif

      if (nucleonsi-nucleonsf==1.and..not.V3Nint) return

      mult_nucl_dim=0
      do ia_1=1,nucleons-2
         do ia_2=ia_1+1,nucleons-1
            do ia_3=ia_2+1,nucleons
               mult_nucl_dim=mult_nucl_dim+1
            end do
         end do
      end do
      if (iproc==0) print *,' mult_nucl3_dim=',mult_nucl_dim
      mult_nucl3_dim=mult_nucl_dim
      allocate(mult_nucl3(3,mult_nucl_dim))
      ij=0
      do ia_1=1,nucleons-2
         do ia_2=ia_1+1,nucleons-1
            do ia_3=ia_2+1,nucleons
               ij=ij+1
               mult_nucl3(1,ij)=ia_1
               mult_nucl3(2,ij)=ia_2
               mult_nucl3(3,ij)=ia_3
            end do
         end do
      end do
      if (ij/=mult_nucl_dim) then
         print *,'*** error: ii,mult_nucl_dim=',
     $        ij,mult_nucl_dim
         stop
      endif

      nasps=naspsf
      N12_max=nhom12f
      mxsps=mxnwdi*nbit

      mshift=0
      do cr_st_1=1,nasps-1
         N_1=2*n_spi(cr_st_1)+l_spi(cr_st_1)
         do cr_st_2=cr_st_1+1,nasps
            N_2=2*n_spi(cr_st_2)+l_spi(cr_st_2)
            if (N_1+N_2>N12_max) cycle
            mj=m2_spi(cr_st_1)+m2_spi(cr_st_2)
            if (mj>mshift(1)) mshift(1)=mj
         end do
      end do
      do cr_st_1=1,nasps
         N_1=2*n_spi(cr_st_1)+l_spi(cr_st_1)
         do cr_st_2=mxsps+1,mxsps+nasps
            N_2=2*n_spi(cr_st_2)+l_spi(cr_st_2)
            if (N_1+N_2>N12_max) cycle
            mj=m2_spi(cr_st_1)+m2_spi(cr_st_2)
            if (mj>mshift(0)) mshift(0)=mj
         end do
      end do
      
      if (iproc==0) print *,' mshift=',mshift
      
      if (iproc==0) print *,' maxval(mshift)=',maxval(mshift)
      allocate(multiple2_dim(0:1,0:2,0:maxval(mshift)))
      multiple2_dim=0
      allocate(multiple2_point(0:1,0:2,0:maxval(mshift)))
      multiple2_point=0

      mshift2=mshift

      ij=0
      do pi=0,1
         do mt=-2,2,2

            select case(mt)
            case(2)
               cr_st_min_1=1
               cr_st_min_2=2
               cr_st_max_1=nasps-1
               cr_st_max_2=nasps
            case(0)
               cr_st_min_1=1
               cr_st_max_1=nasps
               cr_st_min_2=mxsps+1
               cr_st_max_2=mxsps+nasps
            case(-2)
               cr_st_min_1=mxsps+1
               cr_st_min_2=mxsps+2
               cr_st_max_1=mxsps+nasps-1
               cr_st_max_2=mxsps+nasps
            case default
               cycle
            end select
            mt_ind=(mt+2)/2

            do mj=-mshift(abs(mt)/2),mshift(abs(mt)/2),2
               mj_ind=(mshift(abs(mt)/2)+mj)/2
               multiple2_point(pi,mt_ind,mj_ind)=ij
               num_of_mult=0
               do N123=0,N12_max
                  do cr_st_1=cr_st_min_1,cr_st_max_1
                     N_1=2*n_spi(cr_st_1)+l_spi(cr_st_1)
                     do cr_st_2=max(cr_st_1+1,cr_st_min_2),cr_st_max_2
                        N_2=2*n_spi(cr_st_2)+l_spi(cr_st_2)
                        if (N_1+N_2>N12_max) cycle
                        if (mt2_spi(cr_st_1)+mt2_spi(cr_st_2)
     $                       /=mt) cycle
                        if (m2_spi(cr_st_1)+m2_spi(cr_st_2)
     $                       /=mj) cycle
                        if (mod(l_spi(cr_st_1)+l_spi(cr_st_2)
     $                       ,2)/=pi) cycle
                        if (N_1+N_2/=N123) cycle
                        num_of_mult=num_of_mult+1
                     end do
                  end do
               end do

               ij=ij+num_of_mult
               multiple2_dim(pi,mt_ind,mj_ind)=num_of_mult

               if (num_of_mult==0) cycle
               num_of_mult=0
               do N123=0,N12_max
                  do cr_st_1=cr_st_min_1,cr_st_max_1
                     N_1=2*n_spi(cr_st_1)+l_spi(cr_st_1)
                     do cr_st_2=max(cr_st_1+1,cr_st_min_2),cr_st_max_2
                        N_2=2*n_spi(cr_st_2)+l_spi(cr_st_2)
                        if (N_1+N_2>N12_max) cycle
                        if (mt2_spi(cr_st_1)+mt2_spi(cr_st_2)
     $                       /=mt) cycle
                        if (m2_spi(cr_st_1)+m2_spi(cr_st_2)
     $                       /=mj) cycle
                        if (mod(l_spi(cr_st_1)+l_spi(cr_st_2)
     $                       ,2)/=pi) cycle
                        if (N_1+N_2/=N123) cycle
                        num_of_mult=num_of_mult+1
                     end do
                  end do
               end do
            end do
         end do
      end do
         
      if (iproc==0) print *,' total num_of_mult=',ij
      allocate(multiple2_ist(2,ij))
      multiple2_ist=0
      ij=0

      do pi=0,1
         do mt=-2,2,2
            select case(mt)
            case(2)
               cr_st_min_1=1
               cr_st_min_2=2
               cr_st_max_1=nasps-1
               cr_st_max_2=nasps
            case(0)
               cr_st_min_1=1
               cr_st_max_1=nasps
               cr_st_min_2=mxsps+1
               cr_st_max_2=mxsps+nasps
            case(-2)
               cr_st_min_1=mxsps+1
               cr_st_min_2=mxsps+2
               cr_st_max_1=mxsps+nasps-1
               cr_st_max_2=mxsps+nasps
            case default
               cycle
            end select
            mt_ind=(mt+2)/2
            do mj=-mshift(abs(mt)/2),mshift(abs(mt)/2),2
               mj_ind=(mshift(abs(mt)/2)+mj)/2
               do N123=0,N12_max
                  do cr_st_1=cr_st_min_1,cr_st_max_1
                     N_1=2*n_spi(cr_st_1)+l_spi(cr_st_1)
                     do cr_st_2=max(cr_st_1+1,cr_st_min_2),cr_st_max_2
                        N_2=2*n_spi(cr_st_2)+l_spi(cr_st_2)
                        if (N_1+N_2>N12_max) cycle
                        if (mt2_spi(cr_st_1)+mt2_spi(cr_st_2)
     $                       /=mt) cycle
                        if (m2_spi(cr_st_1)+m2_spi(cr_st_2)
     $                       /=mj) cycle
                        if (mod(l_spi(cr_st_1)+l_spi(cr_st_2)
     $                       ,2)/=pi) cycle
                        if (N_1+N_2/=N123) cycle
                        ij=ij+1
                        multiple2_ist(1,ij)
     $                       =cr_st_1
                        multiple2_ist(2,ij)
     $                       =cr_st_2
                     end do
                  end do
               end do
            end do
         end do
      end do

      if (nucleonsi-nucleonsf==1) return
      if (nucleonsi-nucleonsf==2.and..not.V3Nint) return
      
      mult_nucl_dim=0
      do ia_1=1,nucleons-3
         do ia_2=ia_1+1,nucleons-2
            do ia_3=ia_2+1,nucleons-1
               do ia_4=ia_3+1,nucleons
                  mult_nucl_dim=mult_nucl_dim+1
               end do
            end do
         end do
      end do
      if (iproc==0) print *,' mult_nucl4_dim=',mult_nucl_dim
      mult_nucl4_dim=mult_nucl_dim
      allocate(mult_nucl4(4,mult_nucl_dim))
      ij=0
      do ia_1=1,nucleons-3
         do ia_2=ia_1+1,nucleons-2
            do ia_3=ia_2+1,nucleons-1
               do ia_4=ia_3+1,nucleons
                  ij=ij+1
                  mult_nucl4(1,ij)=ia_1
                  mult_nucl4(2,ij)=ia_2
                  mult_nucl4(3,ij)=ia_3
                  mult_nucl4(4,ij)=ia_4
               end do
            end do
         end do
      end do
      if (ij/=mult_nucl_dim) then
         print *,'*** error: ii,mult_nucl_dim=',
     $        ij,mult_nucl_dim
         stop
      endif

      end subroutine mult_init

      subroutine g_a_1(an_1)
      implicit none
      integer,intent(IN) :: an_1
      real(kind(0.d0)) :: cleb,clebt,clebd
      integer :: k_i,k_f,sp_st_a1,ji,ipi,ti,ki_kernel

      sp_st_a1=iobsind(an_1)

c      if (iproc==0) print *,' an_1,sp_st_a1=',an_1,sp_st_a1

      do k_i=ki,ki+nki-1
         ji=jt2i(k_i)
         ti=it2i(k_i)
         ipi=(-1)**iparityi
         ki_kernel=map_ji(ji/2,ti/2,k_i)
         do k_f=kf,kf+nkf-1
            cleb=clebd(jt2f(k_f),mjtotalf,j2_spi(an_1),
     $           mjtotali-mjtotalf,ji,mjtotali)
            clebt=clebd(it2f(k_f),mttotalf,1,
     $           mttotali-mttotalf,ti,mttotali)
            JpiT(ji/2,ipi,ti/2)%st_fi(k_f,ki_kernel)%g(sp_st_a1)=
     $           JpiT(ji/2,ipi,ti/2)%st_fi(k_f,ki_kernel)%g(sp_st_a1)
     $           +cleb*clebt*bmpi(i1,k_i)*bmpf(i2,k_f)
     $           *real(phase,kind(0.d0))
         end do
      end do
      end subroutine g_a_1

      subroutine h_NN_a_1(cr_1,an_2,an_1)
      use interaction, only: N_sp_int_12max,N_sp_int_2max
      implicit none
      integer,intent(IN) :: an_1,an_2,cr_1
      real(kind(0.d0)) :: cleb,clebt,clebd,V2me_unc
      integer :: k_i,k_f,ji,ipi,ti,ki_kernel,
     $     cr_1_ke,an_2_ke,an_1_ke,N_cr_1,i_nljm,
     $     i_nlj,mj,mt,n,l,m2an,mt2an

      an_1_ke=global_ist1_sp_map(an_1)
      an_2_ke=global_ist1_sp_map(an_2)
      cr_1_ke=global_ist1_sp_map(cr_1)

      m2an=m2_spi(an_1)+m2_spi(an_2)
      mt2an=mt2_spi(an_1)+mt2_spi(an_2)
      N_cr_1=2*n_spi(cr_1)+l_spi(cr_1)

c      if (iproc==0) print *,' an_1,an_2,cr_1,an_1_ke,an_2_ke,cr_1_ke=',
c     $     an_1,an_2,cr_1,an_1_ke,an_2_ke,cr_1_ke

      do i_nlj=1,global_ist1_dim
         mj=mjtotali-mjtotalf
         mt=mttotali-mttotalf
         n=global_ist1(1,i_nlj)
         l=global_ist1(2,i_nlj)
         j=global_ist1(3,i_nlj)
         if (2*n+l+N_cr_1>N_sp_int_12max.or.2*n+l>N_sp_int_2max) exit
         if (mod(l+l_spi(an_1)+l_spi(an_2)+l_spi(cr_1),2)==1) cycle
         if (mj+m2_spi(cr_1)/=m2an.or.mt+mt2_spi(cr_1)/=mt2an) cycle
         if (j<abs(mj).or.abs(mt)>1) cycle
c         if (iproc==0) print *,' i_nlj,n,l,j,mj,mt=',i_nlj,n,l,j,mj,mt
         i_nljm=nljmmt_st(i_nlj,(mt+1)/2)%mj_sp(mj)
         if (i_nljm==cr_1_ke) cycle
c         if (iproc==0) print *,' i_nljm=',i_nljm
         call get_V2me_uncoup(an_1_ke,an_2_ke,cr_1_ke,i_nljm,V2me_unc)
c         if (iproc==0) print *,' V2me_unc=',V2me_unc

         do k_i=ki,ki+nki-1
            ji=jt2i(k_i)
            ti=it2i(k_i)
            ipi=(-1)**iparityi
            ki_kernel=map_ji(ji/2,ti/2,k_i)
            do k_f=kf,kf+nkf-1
               cleb=clebd(jt2f(k_f),mjtotalf,j,
     $              mj,ji,mjtotali)
               clebt=clebd(it2f(k_f),mttotalf,1,
     $              mt,ti,mttotali)
               JpiT(ji/2,ipi,ti/2)%st_fi(k_f,ki_kernel)%h_NN(i_nlj)=
     $              JpiT(ji/2,ipi,ti/2)%st_fi(k_f,ki_kernel)%h_NN(i_nlj)
     $              -cleb*clebt*V2me_unc*bmpi(i1,k_i)*bmpf(i2,k_f)
     $              *real(phase,kind(0.d0))
c               if (iproc==0) print *,' k_i,ki_kernel,k_f,JpiT=',
c     $              k_i,ki_kernel,k_f,
c     $              JpiT(ji/2,ipi,ti/2)%st_fi(k_f,ki_kernel)%h_NN(i_nlj)
            end do
         end do
      end do
      end subroutine h_NN_a_1

      subroutine h_3N_a_1(cr_1,cr_2,an_3,an_2,an_1)
      use v3b, only: N1_max,N12_max,N123_max
      implicit none
      integer,intent(IN) :: an_1,an_2,an_3,cr_1,cr_2
      real(kind(0.d0)) :: cleb,clebt,clebd,V3N,v3b_cJ_unc
      integer :: an_1_ke,an_2_ke,an_3_ke,cr_1_ke,cr_2_ke,m2an,mt2an,N_cr
      integer :: i_nlj,mj,mt,n,l,j,i_nljm,ji,ti,ipi,ki_kernel,k_i,k_f
      integer :: N_cr_1,N_cr_2,N_an_1,N_an_2,N_an_3

      N_cr_1=2*n_spi(cr_1)+l_spi(cr_1)
      N_cr_2=2*n_spi(cr_2)+l_spi(cr_2)
      N_cr=N_cr_1+N_cr_2
      if (N_cr>N12_max) return

      N_an_1=2*n_spi(an_1)+l_spi(an_1)
      N_an_2=2*n_spi(an_2)+l_spi(an_2)
      N_an_3=2*n_spi(an_3)+l_spi(an_3)

      if (N_an_1+N_an_2+N_an_3>N123_max) return
      if (N_an_1+N_an_2>N12_max) return
      if (N_an_1+N_an_3>N12_max) return
      if (N_an_2+N_an_3>N12_max) return

      if (N_an_1>N1_max) return
      if (N_an_2>N1_max) return
      if (N_an_3>N1_max) return
      if (N_cr_1>N1_max) return
      if (N_cr_2>N1_max) return

      m2an=m2_spi(an_1)+m2_spi(an_2)+m2_spi(an_3)
      mt2an=mt2_spi(an_1)+mt2_spi(an_2)+mt2_spi(an_3)

      an_1_ke=global_ist1_sp_map(an_1)
      an_2_ke=global_ist1_sp_map(an_2)
      an_3_ke=global_ist1_sp_map(an_3)
      cr_1_ke=global_ist1_sp_map(cr_1)
      cr_2_ke=global_ist1_sp_map(cr_2)

      do i_nlj=1,global_ist1_dim
         mj=mjtotali-mjtotalf
         mt=mttotali-mttotalf
         n=global_ist1(1,i_nlj)
         l=global_ist1(2,i_nlj)
         j=global_ist1(3,i_nlj)
         if (2*n+l+N_cr>N123_max.or.2*n+l>N1_max) exit
         if (2*n+l+N_cr_1>N12_max.or.2*n+l+N_cr_2>N12_max) exit
         if (mod(l+N_an_1+N_an_2+N_an_3+N_cr_1+N_cr_2,2)==1) cycle
         if (mj+m2_spi(cr_1)+m2_spi(cr_2)/=m2an.or.mt+mt2_spi(cr_1)
     $        +mt2_spi(cr_2)/=mt2an) cycle
         if (j<abs(mj).or.abs(mt)>1) cycle
         i_nljm=nljmmt_st(i_nlj,(mt+1)/2)%mj_sp(mj)
         if (i_nljm==cr_1_ke.or.i_nljm==cr_2_ke) cycle
         V3N=v3b_cJ_unc(an_3_ke,an_2_ke,an_1_ke,cr_2_ke,cr_1_ke,i_nljm)
         do k_i=ki,ki+nki-1
            ji=jt2i(k_i)
            ti=it2i(k_i)
            ipi=(-1)**iparityi
            ki_kernel=map_ji(ji/2,ti/2,k_i)
            do k_f=kf,kf+nkf-1
               cleb=clebd(jt2f(k_f),mjtotalf,j,
     $              mj,ji,mjtotali)
               clebt=clebd(it2f(k_f),mttotalf,1,
     $              mt,ti,mttotali)
               JpiT(ji/2,ipi,ti/2)%st_fi(k_f,ki_kernel)%h_3N(i_nlj)=
     $              JpiT(ji/2,ipi,ti/2)%st_fi(k_f,ki_kernel)%h_3N(i_nlj)
     $              +cleb*clebt*V3N*bmpi(i1,k_i)*bmpf(i2,k_f)
     $              *real(phase,kind(0.d0))
            end do
         end do
      end do
      end subroutine h_3N_a_1

      subroutine g_a_2(an_2,an_1)
      use interaction, only: tbdst_SDkern
      implicit none
      integer,intent(IN) :: an_2,an_1
      real(kind(0.d0)) :: cleb,clebt,cleb2,clebt2,clebd,contr
      integer :: k_i,k_f,ia,ib,ji,ipi,ti,ki_kernel,
     $     I_ab,t_ab,tbin,pi

      ia=iobsind(an_1)
      ib=iobsind(an_2)
      pi=mod(l_spi(an_1)+l_spi(an_2),2)

      do k_i=ki,ki+nki-1
         ji=jt2i(k_i)
         ti=it2i(k_i)
         ipi=(-1)**iparityi
         ki_kernel=map_ji(ji/2,ti/2,k_i)
         do k_f=kf,kf+nkf-1
            do I_ab=max(abs(ji-jt2f(k_f))/2,abs(mjtotali-mjtotalf)/2,
     $           abs(j2_spi(an_1)-j2_spi(an_2))/2),min((ji+jt2f(k_f))/2,
     $           (j2_spi(an_1)+j2_spi(an_2))/2)
               cleb=clebd(jt2f(k_f),mjtotalf,2*I_ab,mjtotali-mjtotalf,
     $              ji,mjtotali)
               cleb2=clebd(j2_spi(an_1),m2_spi(an_1),
     $              j2_spi(an_2),m2_spi(an_2),2*I_ab,mjtotali-mjtotalf)
!SQ               do t_ab=abs(mttotali-mttotalf)/2,1
                  do t_ab=max(abs(ti-it2f(k_f))/2, 
     $                 abs(mttotali-mttotalf)/2),
     $                 (ti+it2f(k_f))/2
                  if (ia==ib.and.mod(I_ab+t_ab,2)==0) cycle
                  clebt=clebd(it2f(k_f),mttotalf,2*t_ab,
     $                 mttotali-mttotalf,ti,mttotali)
                  clebt2=clebd(1,mt2_spi(an_1),1,mt2_spi(an_2),2*t_ab,
     $                 mttotali-mttotalf)
                  contr=cleb*cleb2*clebt*clebt2
     $                 *bmpi(i1,k_i)*bmpf(i2,k_f)*real(phase,kind(0.d0))
     $                 /sqrt(2.d0)
                  tbin=tbdst_SDkern(I_ab,pi)%T(t_ab)%sp1(ia)%sp2(ib)
                  JpiT(ji/2,ipi,ti/2)%st_fi(k_f,ki_kernel)%g(tbin)=
     $                 JpiT(ji/2,ipi,ti/2)%st_fi(k_f,ki_kernel)%g(tbin)
     $                 +contr
                  tbin=tbdst_SDkern(I_ab,pi)%T(t_ab)%sp1(ib)%sp2(ia)
                  JpiT(ji/2,ipi,ti/2)%st_fi(k_f,ki_kernel)%g(tbin)=
     $                 JpiT(ji/2,ipi,ti/2)%st_fi(k_f,ki_kernel)%g(tbin)
     $                 +contr
     $                 *(-1)**(I_ab+t_ab-(j2_spi(an_1)+j2_spi(an_2))/2)
               end do
            end do
         end do
      end do
      end subroutine g_a_2

      subroutine h_NN_a_2(cr_1,an_3,an_2,an_1,iph)
      use interaction, only: tbdst_SDkern,N_sp_int_12max,N_sp_int_2max
      implicit none
      integer,intent(IN) :: cr_1,an_3,an_2,an_1,iph
      real(kind(0.d0)) :: cleb,clebt,cleb2,clebt2,clebd,contr,V2me_unc
      integer :: k_i,k_f,ia,ib,ji,ipi,ti,ki_kernel,
     $     I_ab,t_ab,tbin,pi,cr_1_ke,an_2_ke,an_3_ke,
     $     mj,mt,nb,lb,jb,N_cr_1,ib_nljm

      ia=iobsind(an_1)
      an_2_ke=global_ist1_sp_map(an_2)
      an_3_ke=global_ist1_sp_map(an_3)
      cr_1_ke=global_ist1_sp_map(cr_1)

      N_cr_1=2*n_spi(cr_1)+l_spi(cr_1)

      do ib=1,global_ist1_dim
         mj=m2_spi(an_3)+m2_spi(an_2)-m2_spi(cr_1)
         mt=mt2_spi(an_3)+mt2_spi(an_2)-mt2_spi(cr_1)
         nb=global_ist1(1,ib)
         lb=global_ist1(2,ib)
         jb=global_ist1(3,ib)
         if (2*nb+lb+N_cr_1>N_sp_int_12max.or.2*nb+lb>N_sp_int_2max)
     $        exit
         if (2*nb+lb+2*n_spi(an_1)+l_spi(an_1)>N_sp_int_2max+nhomi) exit
         if (mod(lb+l_spi(an_3)+l_spi(an_2)+l_spi(cr_1),2)==1) cycle
         if (jb<abs(mj).or.abs(mt)>1) cycle
         ib_nljm=nljmmt_st(ib,(mt+1)/2)%mj_sp(mj)
         if (ib_nljm==cr_1_ke) cycle
         call get_V2me_uncoup(an_3_ke,an_2_ke,cr_1_ke,ib_nljm,V2me_unc)
         pi=mod(l_spi(an_1)+lb,2)

         do k_i=ki,ki+nki-1
            ji=jt2i(k_i)
            ti=it2i(k_i)
            ipi=(-1)**iparityi
            ki_kernel=map_ji(ji/2,ti/2,k_i)
            do k_f=kf,kf+nkf-1
               do I_ab=max(abs(ji-jt2f(k_f))/2,abs(mjtotali-mjtotalf)/2,
     $              abs(j2_spi(an_1)-jb)/2),min((ji+jt2f(k_f))/2,
     $              (j2_spi(an_1)+jb)/2)
                  cleb=clebd(jt2f(k_f),mjtotalf,2*I_ab,
     $                 mjtotali-mjtotalf,ji,mjtotali)
                  cleb2=clebd(j2_spi(an_1),m2_spi(an_1),
     $                 jb,mj,2*I_ab,mjtotali-mjtotalf)
!SQ                  do t_ab=abs(mttotali-mttotalf)/2,1
                  do t_ab=max(abs(ti-it2f(k_f))/2, 
     $                 abs(mttotali-mttotalf)/2),
     $                 (ti+it2f(k_f))/2
                     if (ia==ib.and.mod(I_ab+t_ab,2)==0) cycle
                     clebt=clebd(it2f(k_f),mttotalf,2*t_ab,
     $                    mttotali-mttotalf,ti,mttotali)
                     clebt2=clebd(1,mt2_spi(an_1),1,mt,2*t_ab,
     $                    mttotali-mttotalf)
                     contr=cleb*cleb2*clebt*clebt2*V2me_unc
     $                    *bmpi(i1,k_i)*bmpf(i2,k_f)
     $                    *real(phase*iph,kind(0.d0))
     $                    /sqrt(2.d0)
                     tbin=tbdst_SDkern(I_ab,pi)%T(t_ab)%sp1(ia)%sp2(ib)
                     JpiT(ji/2,ipi,ti/2)%st_fi(k_f,ki_kernel)%h_NN(tbin)
     $                    =JpiT(ji/2,ipi,ti/2)%st_fi(k_f,ki_kernel)
     $                    %h_NN(tbin)+contr
                     tbin=tbdst_SDkern(I_ab,pi)%T(t_ab)%sp1(ib)%sp2(ia)
                     JpiT(ji/2,ipi,ti/2)%st_fi(k_f,ki_kernel)%h_NN(tbin)
     $                    =JpiT(ji/2,ipi,ti/2)%st_fi(k_f,ki_kernel)
     $                    %h_NN(tbin)+contr
     $                    *real((-1)**(I_ab+t_ab-(j2_spi(an_1)+jb)/2),
     $                    kind(0.d0))
                  end do
               end do
            end do
         end do
      end do
      end subroutine h_NN_a_2

      subroutine h_3N_a_2a(cr_1,cr_2,an_4,an_3,an_2,an_1,iph)
      use interaction, only: tbdst_SDkern
      use v3b, only: N1_max,N12_max,N123_max
      implicit none
      integer,intent(IN) :: cr_1,cr_2,an_4,an_3,an_2,an_1,iph
      real(kind(0.d0)) :: cleb,clebt,cleb2,clebt2,clebd,contr,
     $     V3N,v3b_cJ_unc
      integer :: an_1_ke,an_2_ke,an_3_ke,cr_1_ke,m2an,mt2an,
     $     cr_2_ke,an_4_ke,m2cr,mt2cr
      integer :: N_cr_1,N_an_1,N_an_2,N_an_3,I_ab,t_ab,ia,ib,iab,
     $     na,la,ja,nb,lb,jb,ma,mb,mta,mtb,pi_ancr,ia_nljm,ib_nljm,pi,
     $     tbin,ji,ti,ipi,ki_kernel,N_cr_2,N_an_4

      N_cr_1=2*n_spi(cr_1)+l_spi(cr_1)
      N_cr_2=2*n_spi(cr_2)+l_spi(cr_2)
      if (N_cr_1>N1_max) return
      if (N_cr_2>N1_max) return
      if (N_cr_1+N_cr_2>N12_max) return

      N_an_1=2*n_spi(an_1)+l_spi(an_1)
      N_an_2=2*n_spi(an_2)+l_spi(an_2)
      N_an_3=2*n_spi(an_3)+l_spi(an_3)
      N_an_4=2*n_spi(an_4)+l_spi(an_4)

      if (N_an_2+N_an_3+N_an_4>N123_max) return
      if (N_an_4+N_an_2>N12_max) return
      if (N_an_4+N_an_3>N12_max) return
      if (N_an_2+N_an_3>N12_max) return

      if (N_an_4>N1_max) return
      if (N_an_2>N1_max) return
      if (N_an_3>N1_max) return

      an_2_ke=global_ist1_sp_map(an_2)
      an_3_ke=global_ist1_sp_map(an_3)
      an_4_ke=global_ist1_sp_map(an_4)
      cr_1_ke=global_ist1_sp_map(cr_1)
      cr_2_ke=global_ist1_sp_map(cr_2)

      m2an=m2_spi(an_4)+m2_spi(an_2)+m2_spi(an_3)
      m2cr=m2_spi(cr_1)+m2_spi(cr_2)
      mt2an=mt2_spi(an_4)+mt2_spi(an_2)+mt2_spi(an_3)
      mt2cr=mt2_spi(cr_1)+mt2_spi(cr_2)
      pi_ancr=l_spi(an_4)+l_spi(an_2)+l_spi(an_3)
     $     +l_spi(cr_1)+l_spi(cr_2)
      mb=m2an-m2cr
      mtb=mt2an-mt2cr

      ia=iobsind(an_1)
      la=l_spi(an_1)

      do ib=1,global_ist1_dim

         nb=global_ist1(1,ib)
         lb=global_ist1(2,ib)
         jb=global_ist1(3,ib)
         if (N_cr_1+N_cr_2+2*nb+lb>N123_max) exit
         if (N_cr_1+2*nb+lb>N12_max) exit
         if (N_cr_2+2*nb+lb>N12_max) exit
         if (2*nb+lb>N1_max) exit
         if (jb<abs(mb).or.abs(mtb)>1) cycle
         if (mod(pi_ancr+lb,2)==1) cycle
         ib_nljm=nljmmt_st(ib,(mtb+1)/2)%mj_sp(mb)
         if (ib_nljm==cr_1_ke.or.ib_nljm==cr_2_ke) cycle
         V3N=v3b_cJ_unc(an_4_ke,an_3_ke,an_2_ke,cr_2_ke,
     $        cr_1_ke,ib_nljm)
         pi=mod(la+lb,2)

         do k_i=ki,ki+nki-1
            ji=jt2i(k_i)
            ti=it2i(k_i)
            ipi=(-1)**iparityi
            ki_kernel=map_ji(ji/2,ti/2,k_i)
            do k_f=kf,kf+nkf-1
               do I_ab=max(abs(ji-jt2f(k_f))/2,abs(mjtotali-mjtotalf)/2,
     $              abs(j2_spi(an_1)-jb)/2),min((ji+jt2f(k_f))/2,
     $              (j2_spi(an_1)+jb)/2)
                  cleb=clebd(jt2f(k_f),mjtotalf,2*I_ab,
     $                 mjtotali-mjtotalf,ji,mjtotali)
                  cleb2=clebd(j2_spi(an_1),m2_spi(an_1),jb,mb,
     $                 2*I_ab,m2_spi(an_1)+mb)
!SQ                  do t_ab=abs(mttotali-mttotalf)/2,1
                  do t_ab=max(abs(ti-it2f(k_f))/2, 
     $                 abs(mttotali-mttotalf)/2),
     $                 (ti+it2f(k_f))/2
                     if (ia==ib.and.mod(I_ab+t_ab,2)==0) cycle
                     clebt=clebd(it2f(k_f),mttotalf,2*t_ab,
     $                    mttotali-mttotalf,ti,mttotali)
                     clebt2=clebd(1,mt2_spi(an_1),1,mtb,
     $                    2*t_ab,mt2_spi(an_1)+mtb)
                     contr=cleb*cleb2*clebt*clebt2*V3N
     $                    *bmpi(i1,k_i)*bmpf(i2,k_f)
     $                    *real(phase*iph,kind(0.d0))
     $                    /sqrt(2.d0)                  
                     tbin=tbdst_SDkern(I_ab,pi)%T(t_ab)%sp1(ia)%sp2(ib)
                     JpiT(ji/2,ipi,ti/2)%st_fi(k_f,ki_kernel)%h_3N(tbin)
     $                    =JpiT(ji/2,ipi,ti/2)%st_fi(k_f,ki_kernel)
     $                    %h_3N(tbin)+contr !*0.5d0
                     tbin=tbdst_SDkern(I_ab,pi)%T(t_ab)%sp1(ib)%sp2(ia)
                     JpiT(ji/2,ipi,ti/2)%st_fi(k_f,ki_kernel)%h_3N(tbin)
     $                    =JpiT(ji/2,ipi,ti/2)%st_fi(k_f,ki_kernel)
     $                    %h_3N(tbin)+contr !*0.5d0
     $                    *real((-1)**(
     $                    (2*I_ab-j2_spi(an_1)-jb+2*t_ab)/2),kind(0.d0))
                  end do
               end do
            end do
         end do
      end do

      end subroutine h_3N_a_2a

      subroutine h_3N_a_2b(cr_1,an_3,an_2,an_1)
      use interaction, only: tbdst_SDkern
      use v3b, only: N1_max,N12_max,N123_max
      implicit none
      integer,intent(IN) :: cr_1,an_3,an_2,an_1
      real(kind(0.d0)) :: cleb,clebt,cleb2,clebt2,clebd,contr,
     $     V3N,v3b_cJ_unc
      integer :: an_1_ke,an_2_ke,an_3_ke,cr_1_ke,m2an,mt2an
      integer :: N_cr_1,N_an_1,N_an_2,N_an_3,I_ab,t_ab,ia,ib,iab,
     $     na,la,ja,nb,lb,jb,ma,mb,mta,mtb,pi_ancr,ia_nljm,ib_nljm,pi,
     $     tbina,tbinb,ji,ti,ipi,ki_kernel

      N_cr_1=2*n_spi(cr_1)+l_spi(cr_1)
      if (N_cr_1>N1_max) return

      N_an_1=2*n_spi(an_1)+l_spi(an_1)
      N_an_2=2*n_spi(an_2)+l_spi(an_2)
      N_an_3=2*n_spi(an_3)+l_spi(an_3)

      if (N_an_1+N_an_2+N_an_3>N123_max) return
      if (N_an_1+N_an_2>N12_max) return
      if (N_an_1+N_an_3>N12_max) return
      if (N_an_2+N_an_3>N12_max) return

      if (N_an_1>N1_max) return
      if (N_an_2>N1_max) return
      if (N_an_3>N1_max) return

      m2an=m2_spi(an_1)+m2_spi(an_2)+m2_spi(an_3)
      mt2an=mt2_spi(an_1)+mt2_spi(an_2)+mt2_spi(an_3)
      pi_ancr=l_spi(an_1)+l_spi(an_2)+l_spi(an_3)+l_spi(cr_1)

      an_1_ke=global_ist1_sp_map(an_1)
      an_2_ke=global_ist1_sp_map(an_2)
      an_3_ke=global_ist1_sp_map(an_3)
      cr_1_ke=global_ist1_sp_map(cr_1)

      do iab=1,ist2_SDkern_dim
         ia=ist2_SDkern(1,iab)
         ib=ist2_SDkern(2,iab)
         I_ab=ist2_SDkern(3,iab)
         t_ab=ist2_SDkern(4,iab)
         na=global_ist1(1,ia)
         la=global_ist1(2,ia)
         ja=global_ist1(3,ia)
         if (N_cr_1+2*na+la>N12_max) cycle
         if (2*na+la>N1_max) cycle
         nb=global_ist1(1,ib)
         lb=global_ist1(2,ib)
         jb=global_ist1(3,ib)
         if (N_cr_1+2*na+la+2*nb+lb>N123_max) cycle
         if (2*na+la+2*nb+lb>N12_max) cycle
         if (N_cr_1+2*nb+lb>N12_max) cycle
         if (2*nb+lb>N1_max) cycle
         if (mod(la+lb+pi_ancr,2)/=0) cycle
         pi=mod(la+lb,2)

         tbina=tbdst_SDkern(I_ab,pi)%T(t_ab)%sp1(ia)%sp2(ib)
c         tbinb=tbdst_SDkern(I_ab,pi)%T(t_ab)%sp1(ib)%sp2(ia)

         do ma=-ja,ja,2
c            do mb=-jb,jb,2
c               if (ma+mb/=mjtotali-mjtotalf) cycle
            mb=mjtotali-mjtotalf-ma
            if (abs(mb)>jb) cycle
            cleb2=clebd(ja,ma,jb,mb,2*I_ab,ma+mb)
            do mta=-1,1,2
               mtb=mttotali-mttotalf-mta
               if (abs(mtb)>1) cycle
               if (m2an/=m2_spi(cr_1)+ma+mb) cycle
               if (mt2an/=mt2_spi(cr_1)+mta+mtb) cycle
               clebt2=clebd(1,mta,1,mtb,2*t_ab,mta+mtb)
               ia_nljm=nljmmt_st(ia,(mta+1)/2)%mj_sp(ma)
               if (ia_nljm==cr_1_ke) cycle
               ib_nljm=nljmmt_st(ib,(mtb+1)/2)%mj_sp(mb)
               if (ib_nljm==cr_1_ke) cycle
               if (ia_nljm==ib_nljm) cycle
               V3N=v3b_cJ_unc(an_3_ke,an_2_ke,an_1_ke,cr_1_ke,
     $              ib_nljm,ia_nljm)

               do k_i=ki,ki+nki-1
                  ji=jt2i(k_i)
                  ti=it2i(k_i)
                  ipi=(-1)**iparityi
                  ki_kernel=map_ji(ji/2,ti/2,k_i)
                  do k_f=kf,kf+nkf-1
                     cleb=clebd(jt2f(k_f),mjtotalf,2*I_ab,
     $                    mjtotali-mjtotalf,ji,mjtotali)
                     clebt=clebd(it2f(k_f),mttotalf,2*t_ab,
     $                    mttotali-mttotalf,ti,mttotali)
                     contr=cleb*cleb2*clebt*clebt2*V3N
     $                    *bmpi(i1,k_i)*bmpf(i2,k_f)
     $                    *real(phase,kind(0.d0))
     $                    /sqrt(2.d0)
                     JpiT(ji/2,ipi,ti/2)%st_fi(k_f,ki_kernel)
     $                    %h_3N(tbina)
     $                    =JpiT(ji/2,ipi,ti/2)%st_fi(k_f,ki_kernel)
     $                    %h_3N(tbina)+contr
c                     JpiT(ji/2,ipi,ti/2)%st_fi(k_f,ki_kernel)
c     $                    %h_3N(tbinb)
c     $                    =JpiT(ji/2,ipi,ti/2)%st_fi(k_f,ki_kernel)
c     $                    %h_3N(tbinb)+contr
c     $                    *real((-1)**(I_ab+t_ab-(ja+jb)/2),kind(0.d0))
                  end do
               end do

            end do
c            end do
         end do
      end do
      end subroutine h_3N_a_2b

      subroutine g_a_3(an_3,an_2,an_1,iph)
      use interaction, only: tbdst_SDkern,thrbdst_SDkern
      implicit none
      integer,intent(IN) :: an_3,an_2,an_1,iph
      real(kind(0.d0)) :: cleb,clebt,cleb2,clebt2,clebd,contr,
     $     cleb3,clebt3
      integer :: k_i,k_f,ia,ib,ji,ipi,ti,ki_kernel,
     $     I_ab,t_ab,tbin,pi,I_abc,t_abc,ic,pi2,thrbin

      ia=iobsind(an_1)
      ib=iobsind(an_2)
      ic=iobsind(an_3)
      pi=mod(l_spi(an_1)+l_spi(an_2)+l_spi(an_3),2)
      pi2=mod(l_spi(an_1)+l_spi(an_2),2)

      do k_i=ki,ki+nki-1
         ji=jt2i(k_i)
         ti=it2i(k_i)
         ipi=(-1)**iparityi
         ki_kernel=map_ji(ji/2,ti/2,k_i)
         do k_f=kf,kf+nkf-1
            do I_abc=max(abs(ji-jt2f(k_f)),abs(mjtotali-mjtotalf)),
     $           ji+jt2f(k_f),2
               cleb=clebd(jt2f(k_f),mjtotalf,I_abc,mjtotali-mjtotalf,
     $              ji,mjtotali)
               do I_ab=max(abs(I_abc-j2_spi(an_3))/2,
     $              abs(j2_spi(an_1)-j2_spi(an_2))/2,
     $              abs(m2_spi(an_1)+m2_spi(an_2))/2),
     $              min((I_abc+j2_spi(an_3))/2,
     $              (j2_spi(an_1)+j2_spi(an_2))/2)
                  cleb2=clebd(j2_spi(an_1),m2_spi(an_1),
     $                 j2_spi(an_2),m2_spi(an_2),2*I_ab,
     $                 m2_spi(an_1)+m2_spi(an_2))
                  cleb3=clebd(2*I_ab,m2_spi(an_1)+m2_spi(an_2),
     $                 j2_spi(an_3),m2_spi(an_3),
     $                 I_abc,mjtotali-mjtotalf)
                  do t_abc=max(abs(mttotali-mttotalf),abs(ti-it2f(k_f)))
     $                 ,min(ti+it2f(k_f),3),2
                     clebt=clebd(it2f(k_f),mttotalf,t_abc,
     $                    mttotali-mttotalf,ti,mttotali)
                     do t_ab=max(abs(mt2_spi(an_1)+mt2_spi(an_2))/2,
     $                    abs(t_abc-1)/2),min(t_abc+1,2)/2
                        if (ia==ib.and.mod(I_ab+t_ab,2)==0) cycle

                        clebt2=clebd(1,mt2_spi(an_1),1,mt2_spi(an_2),
     $                       2*t_ab,mt2_spi(an_1)+mt2_spi(an_2))
                        clebt3=clebd(2*t_ab,mt2_spi(an_1)+mt2_spi(an_2),
     $                       1,mt2_spi(an_3),t_abc,mttotali-mttotalf)
                        contr=cleb*cleb2*cleb3*clebt*clebt2*clebt3  ! no minus sign here!!! 
     $                       *bmpi(i1,k_i)*bmpf(i2,k_f)
     $                       *real(phase*iph,kind(0.d0))
     $                       /sqrt(6.d0)
                        tbin=tbdst_SDkern(I_ab,pi2)%T(t_ab)%sp1(ia)
     $                       %sp2(ib)
                        thrbin=thrbdst_SDkern(I_abc/2,pi)%T(t_abc/2)
     $                       %sp1(tbin)%sp2(ic)
                        JpiT(ji/2,ipi,ti/2)%st_fi(k_f,ki_kernel)
     $                       %g(thrbin)=
     $                       JpiT(ji/2,ipi,ti/2)%st_fi(k_f,ki_kernel)
     $                       %g(thrbin)+contr
                        tbin=tbdst_SDkern(I_ab,pi2)%T(t_ab)%sp1(ib)
     $                       %sp2(ia)
                        thrbin=thrbdst_SDkern(I_abc/2,pi)%T(t_abc/2)
     $                       %sp1(tbin)%sp2(ic)
                        JpiT(ji/2,ipi,ti/2)%st_fi(k_f,ki_kernel)
     $                       %g(thrbin)=
     $                       JpiT(ji/2,ipi,ti/2)%st_fi(k_f,ki_kernel)
     $                       %g(thrbin)+contr*(-1)**(I_ab+t_ab
     $                       -(j2_spi(an_1)+j2_spi(an_2))/2)
                     end do
                  end do
               end do
            end do
         end do
      end do
      end subroutine g_a_3

      subroutine h_NN_a_3(cr_1,an_4,an_3,an_2,an_1,iph)
      use interaction, only: tbdst_SDkern,thrbdst_SDkern,N_sp_int_12max,
     $     N_sp_int_2max
      implicit none
      integer,intent(IN) :: cr_1,an_4,an_3,an_2,an_1,iph
      real(kind(0.d0)) :: cleb,clebt,clebd,contr,V2me_unc,
     $     clebt_ab,clebt_abc,clebt_ca,clebt_cab,clebt_bc,clebt_bca,
     $     cleb_ab,cleb_abc,cleb_ca,cleb_cab,cleb_bc,cleb_bca
      integer :: k_i,k_f,ia,ib,ji,ipi,ti,ki_kernel,mtc,N_cr_1,ic_nljm,
     $     I_ab,t_ab,tbin,pi,I_abc,t_abc,ic,nc,lc,jc,mc,id2,id3,id3p,
     $     an_2_ke,an_3_ke,cr_1_ke,I_ab_min,I_ab_max,I_ca_min,I_ca_max,
     $     I_bc_min,I_bc_max,pi2,thrbin

      ia=iobsind(an_1)
      ib=iobsind(an_2)
      id2=iobsind(an_3)
      id3=iobsind(an_4)
      id3p=iobsind(cr_1)
      an_2_ke=global_ist1_sp_map(an_3)
      an_3_ke=global_ist1_sp_map(an_4)
      cr_1_ke=global_ist1_sp_map(cr_1)
      N_cr_1=2*n_spi(cr_1)+l_spi(cr_1)

      do ic=1,global_ist1_dim
         mc=m2_spi(an_3)+m2_spi(an_4)-m2_spi(cr_1)
         mtc=mt2_spi(an_3)+mt2_spi(an_4)-mt2_spi(cr_1)
         nc=global_ist1(1,ic)
         lc=global_ist1(2,ic)
         jc=global_ist1(3,ic)
         if (2*nc+lc+N_cr_1>N_sp_int_12max.or.2*nc+lc>N_sp_int_2max)
     $        exit
cc         if (2*nc+lc+2*n_spi(an_1)+l_spi(an_1)>N_sp_int_2max+nhomi) exit
         if (mod(lc+l_spi(an_3)+l_spi(an_4)+l_spi(cr_1),2)==1) cycle
         if (jc<abs(mc).or.abs(mtc)>1) cycle
         ic_nljm=nljmmt_st(ic,(mtc+1)/2)%mj_sp(mc)
         if (ic_nljm==cr_1_ke) cycle
         call get_V2me_uncoup(an_3_ke,an_2_ke,cr_1_ke,ic_nljm,V2me_unc)
         if (abs(V2me_unc)<1.d-12) cycle 

         pi=mod(l_spi(an_1)+l_spi(an_2)+lc,2)

         do k_i=ki,ki+nki-1
            ji=jt2i(k_i)
            ti=it2i(k_i)
            ipi=(-1)**iparityi
            ki_kernel=map_ji(ji/2,ti/2,k_i)
            do k_f=kf,kf+nkf-1

               do I_abc=max(abs(ji-jt2f(k_f)),abs(mjtotali-mjtotalf)),
     $              ji+jt2f(k_f),2
                  cleb=clebd(jt2f(k_f),mjtotalf,I_abc,mjtotali-mjtotalf,
     $                 ji,mjtotali)

                  I_ab_min=max(abs(I_abc-jc)/2,
     $                 abs(m2_spi(an_1)+m2_spi(an_2))/2,
     $                 abs(j2_spi(an_1)-j2_spi(an_2))/2)
                  I_ab_max=min((I_abc+jc)/2,
     $                 (j2_spi(an_1)+j2_spi(an_2))/2)
                  I_ca_min=max(abs(I_abc-j2_spi(an_2))/2,
     $                 abs(m2_spi(an_1)+mc)/2,
     $                 abs(j2_spi(an_1)-jc)/2)
                  I_ca_max=min((I_abc+j2_spi(an_2))/2,
     $                 (j2_spi(an_1)+jc)/2)
                  I_bc_min=max(abs(I_abc-j2_spi(an_1))/2,
     $                 abs(m2_spi(an_2)+mc)/2,
     $                 abs(j2_spi(an_2)-jc)/2)
                  I_bc_max=min((I_abc+j2_spi(an_1))/2,
     $                 (j2_spi(an_2)+jc)/2)

                  do t_abc=max(abs(mttotali-mttotalf),abs(ti-it2f(k_f)))
     $                 ,min(ti+it2f(k_f),3),2
                     clebt=clebd(it2f(k_f),mttotalf,t_abc,
     $                    mttotali-mttotalf,ti,mttotali)

                     do t_ab=abs(t_abc-1)/2,min(t_abc+1,2)/2

                        if (abs(mt2_spi(an_1)+mt2_spi(an_2))/2<=t_ab)
     $                       then
                           clebt_ab=
     $                          clebd(1,mt2_spi(an_1),1,mt2_spi(an_2),
     $                          2*t_ab,mt2_spi(an_1)+mt2_spi(an_2))
                           clebt_abc=clebd(2*t_ab,
     $                          mt2_spi(an_1)+mt2_spi(an_2),
     $                          1,mtc,t_abc,mttotali-mttotalf)
                        else
                           clebt_ab=0.d0
                           clebt_abc=0.d0
                        endif
                        if (abs(mt2_spi(an_1)+mtc)/2<=t_ab)
     $                       then
                           clebt_ca=
     $                          clebd(1,mtc,1,mt2_spi(an_1),
     $                          2*t_ab,mt2_spi(an_1)+mtc)
                           clebt_cab=clebd(2*t_ab,
     $                          mt2_spi(an_1)+mtc,
     $                          1,mt2_spi(an_2),t_abc,mttotali-mttotalf)
                        else
                           clebt_ca=0.d0
                           clebt_cab=0.d0
                        endif
                        if (abs(mt2_spi(an_2)+mtc)/2<=t_ab)
     $                       then
                           clebt_bc=
     $                          clebd(1,mt2_spi(an_2),1,mtc,
     $                          2*t_ab,mt2_spi(an_2)+mtc)
                           clebt_bca=clebd(2*t_ab,
     $                          mt2_spi(an_2)+mtc,
     $                          1,mt2_spi(an_1),t_abc,mttotali-mttotalf)
                        else
                           clebt_bc=0.d0
                           clebt_bca=0.d0
                        endif

                        do I_ab=min(I_ab_min,I_ca_min,I_bc_min),
     $                       max(I_ab_max,I_ca_max,I_bc_max)

                           if (I_ab>=I_ab_min.and.I_ab<=I_ab_max.and.
     $                          .not.(ia==ib.and.mod(I_ab+t_ab,2)==0)
     $                          .and.clebt_ab/=0.and.clebt_abc/=0.d0)
     $                          then

                              pi2=mod(l_spi(an_1)+l_spi(an_2),2)

                              cleb_ab=clebd(j2_spi(an_1),m2_spi(an_1),
     $                             j2_spi(an_2),m2_spi(an_2),2*I_ab,
     $                             m2_spi(an_1)+m2_spi(an_2))

                              cleb_abc=clebd(2*I_ab,m2_spi(an_1)
     $                             +m2_spi(an_2),jc,mc,
     $                             I_abc,mjtotali-mjtotalf)

                              contr=cleb*cleb_ab*cleb_abc*clebt*clebt_ab
     $                             *clebt_abc*V2me_unc
     $                             *bmpi(i1,k_i)*bmpf(i2,k_f)
     $                             *real(phase*iph,kind(0.d0))
     $                             /sqrt(6.d0)

                              tbin=tbdst_SDkern(I_ab,pi2)%T(t_ab)
     $                             %sp1(ia)%sp2(ib)
                              thrbin=thrbdst_SDkern(I_abc/2,pi)
     $                             %T(t_abc/2)%sp1(tbin)%sp2(ic)

                              JpiT(ji/2,ipi,ti/2)%st_fi(k_f,ki_kernel)
     $                             %h_NN(thrbin)=JpiT(ji/2,ipi,ti/2)
     $                             %st_fi(k_f,ki_kernel)
     $                             %h_NN(thrbin)+contr

                              tbin=tbdst_SDkern(I_ab,pi2)%T(t_ab)
     $                             %sp1(ib)%sp2(ia)
                              thrbin=thrbdst_SDkern(I_abc/2,pi)
     $                             %T(t_abc/2)%sp1(tbin)%sp2(ic)

                              JpiT(ji/2,ipi,ti/2)%st_fi(k_f,ki_kernel)
     $                             %h_NN(thrbin)=JpiT(ji/2,ipi,ti/2)
     $                             %st_fi(k_f,ki_kernel)
     $                             %h_NN(thrbin)+contr*(-1)**(I_ab+t_ab
     $                             -(j2_spi(an_1)+j2_spi(an_2))/2)
                           endif

                           if (I_ab>=I_ca_min.and.I_ab<=I_ca_max.and.
     $                          .not.(ia==ic.and.mod(I_ab+t_ab,2)==0)
     $                          .and.clebt_ca/=0.and.clebt_cab/=0.d0)
     $                          then

                              pi2=mod(l_spi(an_1)+lc,2)

                              cleb_ca=clebd(jc,mc,
     $                             j2_spi(an_1),m2_spi(an_1),
     $                             2*I_ab,m2_spi(an_1)+mc)

                              cleb_cab=clebd(2*I_ab,m2_spi(an_1)+mc,
     $                             j2_spi(an_2),m2_spi(an_2),
     $                             I_abc,mjtotali-mjtotalf)

                              contr=cleb*cleb_ca*cleb_cab*clebt*clebt_ca
     $                             *clebt_cab*V2me_unc
     $                             *bmpi(i1,k_i)*bmpf(i2,k_f)
     $                             *real(phase*iph,kind(0.d0))
     $                             /sqrt(6.d0)

                              tbin=tbdst_SDkern(I_ab,pi2)%T(t_ab)
     $                             %sp1(ic)%sp2(ia)
                              thrbin=thrbdst_SDkern(I_abc/2,pi)
     $                             %T(t_abc/2)%sp1(tbin)%sp2(ib)

                              JpiT(ji/2,ipi,ti/2)%st_fi(k_f,ki_kernel)
     $                             %h_NN(thrbin)=JpiT(ji/2,ipi,ti/2)
     $                             %st_fi(k_f,ki_kernel)
     $                             %h_NN(thrbin)+contr

                              tbin=tbdst_SDkern(I_ab,pi2)%T(t_ab)
     $                             %sp1(ia)%sp2(ic)
                              thrbin=thrbdst_SDkern(I_abc/2,pi)
     $                             %T(t_abc/2)%sp1(tbin)%sp2(ib)

                              JpiT(ji/2,ipi,ti/2)%st_fi(k_f,ki_kernel)
     $                             %h_NN(thrbin)=JpiT(ji/2,ipi,ti/2)
     $                             %st_fi(k_f,ki_kernel)
     $                             %h_NN(thrbin)+contr*(-1)**(I_ab+t_ab
     $                             -(j2_spi(an_1)+jc)/2)
                           endif

                           if (I_ab>=I_bc_min.and.I_ab<=I_bc_max.and.
     $                          .not.(ib==ic.and.mod(I_ab+t_ab,2)==0)
     $                          .and.clebt_bc/=0.and.clebt_bca/=0.d0)
     $                          then

                              pi2=mod(l_spi(an_2)+lc,2)

                              cleb_bc=clebd(j2_spi(an_2),m2_spi(an_2),
     $                             jc,mc,2*I_ab,m2_spi(an_2)+mc)

                              cleb_bca=clebd(2*I_ab,m2_spi(an_2)+mc,
     $                             j2_spi(an_1),m2_spi(an_1),
     $                             I_abc,mjtotali-mjtotalf)

                              contr=cleb*cleb_bc*cleb_bca*clebt*clebt_bc
     $                             *clebt_bca*V2me_unc
     $                             *bmpi(i1,k_i)*bmpf(i2,k_f)
     $                             *real(phase*iph,kind(0.d0))
     $                             /sqrt(6.d0)

                              tbin=tbdst_SDkern(I_ab,pi2)%T(t_ab)
     $                             %sp1(ib)%sp2(ic)
                              thrbin=thrbdst_SDkern(I_abc/2,pi)
     $                             %T(t_abc/2)%sp1(tbin)%sp2(ia)

                              JpiT(ji/2,ipi,ti/2)%st_fi(k_f,ki_kernel)
     $                             %h_NN(thrbin)=JpiT(ji/2,ipi,ti/2)
     $                             %st_fi(k_f,ki_kernel)
     $                             %h_NN(thrbin)+contr

                              tbin=tbdst_SDkern(I_ab,pi2)%T(t_ab)
     $                             %sp1(ic)%sp2(ib)
                              thrbin=thrbdst_SDkern(I_abc/2,pi)
     $                             %T(t_abc/2)%sp1(tbin)%sp2(ia)

                              JpiT(ji/2,ipi,ti/2)%st_fi(k_f,ki_kernel)
     $                             %h_NN(thrbin)=JpiT(ji/2,ipi,ti/2)
     $                             %st_fi(k_f,ki_kernel)
     $                             %h_NN(thrbin)+contr*(-1)**(I_ab+t_ab
     $                             -(j2_spi(an_2)+jc)/2)
                           endif

                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
      end subroutine h_NN_a_3

      end

      subroutine NCSMC_coupling_kernels_phys_trans
      use paramdef
      use initst
      use finast
      use intrface, only: ki,nki,kf,nkf
      use kernels
      implicit none
      integer :: ii,k_f,s,l,n,I1,T1,j,Ttot,Jtot,ip,k_i,ji,ti,ipi,
     $     i_nlj,n1,l1,j1,ki_kernel
      real(kind(0.d0)) :: sum1,sum2,sum3,sixj,racad,cms_cor

      do Jtot=jmi,jma,2
         do ip=ipamin,ipamax
            do Ttot=tmi,tma,2
               ii=0
               do k_f=kf,kf+nkf-1
                  I1=jt2f(k_f)
                  T1=it2f(k_f)
                  if (abs(T1-1)>Ttot.or.T1+1<Ttot) cycle
                  do s=abs(I1-1),I1+1,2
                     do l=abs(Jtot-s)/2,min((Jtot+s)/2,N_RGM_max)
                        if ((-1)**(l+iparityf)/=ip) cycle
                        do n=0,(N_RGM_max-l)/2
                           ii=ii+1
                        end do
                     end do
                  end do
               end do
               JpiT(Jtot/2,ip,Ttot/2)%phdim=ii
               allocate(JpiT(Jtot/2,ip,Ttot/2)%phchan(4,ii))
               ii=0
               do k_f=kf,kf+nkf-1
                  I1=jt2f(k_f)
                  T1=it2f(k_f)
                  if (abs(T1-1)>Ttot.or.T1+1<Ttot) cycle
                  do s=abs(I1-1),I1+1,2
                     do l=abs(Jtot-s)/2,min((Jtot+s)/2,N_RGM_max)
                        if ((-1)**(l+iparityf)/=ip) cycle
                        do n=0,(N_RGM_max-l)/2
                           ii=ii+1
                           JpiT(Jtot/2,ip,Ttot/2)%phchan(1,ii)=k_f
                           JpiT(Jtot/2,ip,Ttot/2)%phchan(2,ii)=s
                           JpiT(Jtot/2,ip,Ttot/2)%phchan(3,ii)=l
                           JpiT(Jtot/2,ip,Ttot/2)%phchan(4,ii)=n
                        end do
                     end do
                  end do
               end do
               allocate(JpiT(Jtot/2,ip,Ttot/2)%ovl_chan_i(
     $              JpiT(Jtot/2,ip,Ttot/2)%phdim,
     $              JpiT(Jtot/2,ip,Ttot/2)%dim_i))
               do k_i=ki,ki+nki-1
                  ji=jt2i(k_i)
                  ti=it2i(k_i)
                  ipi=(-1)**iparityi
                  if (ji/=Jtot.or.ipi/=ip.or.ti/=Ttot) cycle
                  ki_kernel=map_ji(ji/2,ti/2,k_i)
                  do ii=1,JpiT(Jtot/2,ip,Ttot/2)%phdim
                     allocate(JpiT(Jtot/2,ip,Ttot/2)
     $                    %ovl_chan_i(ii,ki_kernel)%g(1))
                     allocate(JpiT(Jtot/2,ip,Ttot/2)
     $                    %ovl_chan_i(ii,ki_kernel)%h_NN(1))
                     allocate(JpiT(Jtot/2,ip,Ttot/2)
     $                    %ovl_chan_i(ii,ki_kernel)%h_3N(1))
                     JpiT(Jtot/2,ip,Ttot/2)
     $                    %ovl_chan_i(ii,ki_kernel)%g=0.d0
                     JpiT(Jtot/2,ip,Ttot/2)
     $                    %ovl_chan_i(ii,ki_kernel)%h_NN=0.d0
                     JpiT(Jtot/2,ip,Ttot/2)
     $                    %ovl_chan_i(ii,ki_kernel)%h_3N=0.d0
                     k_f=JpiT(Jtot/2,ip,Ttot/2)%phchan(1,ii)
                     s=JpiT(Jtot/2,ip,Ttot/2)%phchan(2,ii)
                     l=JpiT(Jtot/2,ip,Ttot/2)%phchan(3,ii)
                     n=JpiT(Jtot/2,ip,Ttot/2)%phchan(4,ii)
                     cms_cor=sqrt((real(nucleonsi,kind(0.d0))
     $                    /real(nucleonsf,kind(0.d0)))**(2*n+l))
     $                    *real((-1)**l,kind(0.d0))
                     sum1=0.d0
                     sum2=0.d0
                     sum3=0.d0
                     do i_nlj=1,global_ist1_dim
                        n1=global_ist1(1,i_nlj)
                        l1=global_ist1(2,i_nlj)
                        j1=global_ist1(3,i_nlj)
                        if (n1/=n.or.l1/=l) cycle
                        sixj=racad(jt2f(k_f),1,ji,2*l,s,j1)
     $                       *real((-1)**((j1-1)/2-l),kind(0.d0))
     $                       *sqrt(real((s+1)*(j1+1),kind(0.d0)))
                        sum1=sum1+cms_cor*sixj*JpiT(Jtot/2,ip,Ttot/2)
     $                       %st_fi(k_f,ki_kernel)%g(i_nlj)
                        sum2=sum2+cms_cor*sixj*JpiT(Jtot/2,ip,Ttot/2)
     $                       %st_fi(k_f,ki_kernel)%h_NN(i_nlj)
                        sum3=sum3+cms_cor*sixj*JpiT(Jtot/2,ip,Ttot/2)
     $                       %st_fi(k_f,ki_kernel)%h_3N(i_nlj)
                     end do
                     JpiT(Jtot/2,ip,Ttot/2)
     $                    %ovl_chan_i(ii,ki_kernel)%g(1)=sum1
                     JpiT(Jtot/2,ip,Ttot/2)
     $                    %ovl_chan_i(ii,ki_kernel)%h_NN(1)=sum2
                     JpiT(Jtot/2,ip,Ttot/2)
     $                    %ovl_chan_i(ii,ki_kernel)%h_3N(1)=sum3
                  end do
               end do
            end do
         end do
      end do
      end

      subroutine NCSMC_coupling_kernels_write
      use constants, only: rdmavg,hbc
      use paramdef
      use initst
      use finast
      use intrface, only: ki,nki,kf,nkf
      use kernels
      use interaction, only: N_sp_int_1max
      use v3b, only: N1_max,N12_max,N123_max,V3Nint
      implicit none
      integer :: ii,kk,j,ip,t,Nkernelmax,kernelist1dim,i_nlj,n,l,
     $     ki_kernel,k_i,k_f,ji,ti,ipi,s,dimalloc,dimnonzero
      real(kind(0.d0)) :: mass

      write(unit_NCSMC,*) hboi
      ii=1
      write(unit_NCSMC,*) ii
      write(unit_NCSMC,*) N_cluster_max,N_cluster_min,
     $     N_sp_min,N_proj_max
      write(unit_NCSMC,*) N1_max,N12_max,N123_max
      write(unit_NCSMC,*) N_RGM_max

      write(unit_NCSMC,*) 
      write(unit_NCSMC,*) nprotonsi,nneutrnsi
      write(unit_NCSMC,*) nki
      write(unit_NCSMC,*) jmi,jma,ipamin,ipamax,tmi,tma
      write(unit_NCSMC,*) 
      do j=jmi,jma,2
         do ip=ipamin,ipamax
            do t=tmi,tma,2
               write(unit_NCSMC,*) j,ip,t,JpiT(j/2,ip,t/2)%dim_i
               do k_i=1,JpiT(j/2,ip,t/2)%dim_i
                  do kk=ki,ki+nki-1
                     ji=jt2i(kk)
                     ti=it2i(kk)
                     ipi=(-1)**iparityi
                     if (ji/=j.or.ipi/=ip.or.ti/=t) cycle
                     ki_kernel=map_ji(ji/2,ti/2,kk)
                     if (k_i/=ki_kernel) cycle
                     write(unit_NCSMC,*) jt2i(kk),it2i(kk),eneri(kk)
                  end do
               end do
            end do
         end do
      end do
      write(unit_NCSMC,*) nhomi,nhom12i,nhwi
      write(unit_NCSMC,*) mjtotali,mttotali

      mass=rdmavg*hboi/(hbc**2)
     $     *(real(nucleonsf,kind(0.d0))
     $     *real(num_prot_proj+num_neutr_proj,kind(0.d0)))
     $     /real(nucleonsi,kind(0.d0))
      write(unit_NCSMC,*) 
      write(unit_NCSMC,*) ii,mass
      write(unit_NCSMC,*) num_prot_proj,num_neutr_proj

      write(unit_NCSMC,*) 
      write(unit_NCSMC,*) nprotonsf,nneutrnsf
      write(unit_NCSMC,*) nkf
      do kk=kf,kf+nkf-1
         write(unit_NCSMC,*) jt2f(kk),it2f(kk),enerf(kk)
      end do
      write(unit_NCSMC,*) nhomf,nhom12f,nhwf
      write(unit_NCSMC,*) mjtotalf,mttotalf

      write(unit_NCSMC,*) 
      write(unit_NCSMC,*) global_ist1_dim
      do ii=1,global_ist1_dim
         write(unit_NCSMC,*) ii,(global_ist1(kk,ii),kk=1,3)
      end do
      write(unit_NCSMC,*) 

      select case(nucleonsi-nucleonsf)
      case(1)
         dimalloc=global_ist1_dim
      case(2)
         dimalloc=ist2_SDkern_dim
         write(unit_NCSMC,*) ist2_SDkern_dim
         do ii=1,dimalloc
            write(unit_NCSMC,*) ii,(ist2_SDkern(kk,ii),kk=1,4)
         end do
         write(unit_NCSMC,*) 
      case(3)
         dimalloc=ist3_SDkern_dim
         write(unit_NCSMC,*) ist2_SDkern_dim
         do ii=1,ist2_SDkern_dim
            write(unit_NCSMC,*) ii,(ist2_SDkern(kk,ii),kk=1,4)
         end do
         write(unit_NCSMC,*) 
         write(unit_NCSMC,*) ist3_SDkern_dim
         do ii=1,ist3_SDkern_dim
            write(unit_NCSMC,*) ii,(ist3_SDkern(kk,ii),kk=1,4)
         end do
         write(unit_NCSMC,*) 
      case default
         write(unit_NCSMC,*) 
     +        ' Not implemented for nucleonsi and nucleonsf',
     +        nucleonsi,nucleonsf
         stop
      end select

c      Nkernelmax=max(2*n_spi(naspsi)+l_spi(naspsi),
c     $     N_sp_int_1max,N1_max)
c      ii=0
c      do i_nlj=1,global_ist1_dim
c         n=global_ist1(1,i_nlj)
c         l=global_ist1(2,i_nlj)
c         if (2*n+l>Nkernelmax) exit
c         ii=ii+1
c      end do
c      kernelist1dim=ii
c      write(unit_NCSMC,*) 
c      write(unit_NCSMC,*) Nkernelmax,kernelist1dim

      do j=jmi,jma,2
         do ip=ipamin,ipamax
            do t=tmi,tma,2

               write(unit_NCSMC,*) j,ip,t
               if (nucleonsi==nucleonsf+1) then
                  write(unit_NCSMC,*)
     $                 JpiT(j/2,ip,t/2)%dim_i,JpiT(j/2,ip,t/2)%phdim
               else
                  write(unit_NCSMC,*)
     $                 JpiT(j/2,ip,t/2)%dim_i,JpiT(j/2,ip,t/2)%dim_f
               endif

               do k_i=ki,ki+nki-1
                  ji=jt2i(k_i)
                  ti=it2i(k_i)
                  ipi=(-1)**iparityi
                  if (ji/=j.or.ipi/=ip.or.ti/=t) cycle
                  ki_kernel=map_ji(ji/2,ti/2,k_i)
                  if (nucleonsi==nucleonsf+1) then
                     write(unit_NCSMC,*)
                     write(unit_NCSMC,*) JpiT(j/2,ip,t/2)%phdim
                     do ii=1,JpiT(j/2,ip,t/2)%phdim
                        k_f=JpiT(j/2,ip,t/2)%phchan(1,ii)
                        s=JpiT(j/2,ip,t/2)%phchan(2,ii)
                        l=JpiT(j/2,ip,t/2)%phchan(3,ii)
                        n=JpiT(j/2,ip,t/2)%phchan(4,ii)
                        write(unit_NCSMC,*) ki_kernel,k_f,s,l,n
                        write(unit_NCSMC,*) JpiT(j/2,ip,t/2)
     $                       %ovl_chan_i(ii,ki_kernel)%g(1),
     $                       JpiT(j/2,ip,t/2)
     $                       %ovl_chan_i(ii,ki_kernel)%h_NN(1),
     $                       JpiT(j/2,ip,t/2)
     $                       %ovl_chan_i(ii,ki_kernel)%h_3N(1)
                     end do
                  else
                     do k_f=1,JpiT(j/2,ip,t/2)%dim_f
                        
                        dimnonzero=0
                        do i_nlj=JpiT(ji/2,ipi,ti/2)
     $                       %st_fi(k_f,ki_kernel)%start_1,
     $                       JpiT(ji/2,ipi,ti/2)
     $                       %st_fi(k_f,ki_kernel)%dim_1
                           if (V3Nint) then
                              if (abs(JpiT(ji/2,ipi,ti/2)
     $                             %st_fi(k_f,ki_kernel)
     $                             %g(i_nlj))>1.d-9.or.
     $                             abs(JpiT(ji/2,ipi,ti/2)
     $                             %st_fi(k_f,ki_kernel)
     $                             %h_NN(i_nlj))>1.d-9.or.
     $                             abs(JpiT(ji/2,ipi,ti/2)
     $                             %st_fi(k_f,ki_kernel)
     $                             %h_3N(i_nlj))>1.d-9) then
                                 dimnonzero=dimnonzero+1
                              endif
                           else
                              if (abs(JpiT(ji/2,ipi,ti/2)
     $                             %st_fi(k_f,ki_kernel)
     $                             %g(i_nlj))>1.d-9.or.
     $                             abs(JpiT(ji/2,ipi,ti/2)
     $                             %st_fi(k_f,ki_kernel)
     $                             %h_NN(i_nlj))>1.d-9) then
                                 dimnonzero=dimnonzero+1
                              endif
                           endif
                        end do

                        write(unit_NCSMC,*)
                        write(unit_NCSMC,*) dimnonzero,
     $                       JpiT(ji/2,ipi,ti/2)
     $                       %st_fi(k_f,ki_kernel)%start_1,
     $                       JpiT(ji/2,ipi,ti/2)
     $                       %st_fi(k_f,ki_kernel)%dim_1
                        do i_nlj=JpiT(ji/2,ipi,ti/2)
     $                       %st_fi(k_f,ki_kernel)%start_1,
     $                       JpiT(ji/2,ipi,ti/2)
     $                       %st_fi(k_f,ki_kernel)%dim_1
                           if (V3Nint) then
                              if (abs(JpiT(ji/2,ipi,ti/2)
     $                             %st_fi(k_f,ki_kernel)
     $                             %g(i_nlj))>1.d-9.or.
     $                             abs(JpiT(ji/2,ipi,ti/2)
     $                             %st_fi(k_f,ki_kernel)
     $                             %h_NN(i_nlj))>1.d-9.or.
     $                             abs(JpiT(ji/2,ipi,ti/2)
     $                             %st_fi(k_f,ki_kernel)
     $                             %h_3N(i_nlj))>1.d-9) then
                                 write(unit_NCSMC,*) k_f,ki_kernel,
     $                                i_nlj,JpiT(ji/2,ipi,ti/2)
     $                                %st_fi(k_f,ki_kernel)%g(i_nlj),
     $                                JpiT(ji/2,ipi,ti/2)
     $                                %st_fi(k_f,ki_kernel)%h_NN(i_nlj),
     $                                JpiT(ji/2,ipi,ti/2)
     $                                %st_fi(k_f,ki_kernel)%h_3N(i_nlj)
                              endif
                           else
                              if (abs(JpiT(ji/2,ipi,ti/2)
     $                             %st_fi(k_f,ki_kernel)
     $                             %g(i_nlj))>1.d-9.or.
     $                             abs(JpiT(ji/2,ipi,ti/2)
     $                             %st_fi(k_f,ki_kernel)
     $                             %h_NN(i_nlj))>1.d-9) then
                                 write(unit_NCSMC,*) k_f,ki_kernel,
     $                                i_nlj,JpiT(ji/2,ipi,ti/2)
     $                                %st_fi(k_f,ki_kernel)%g(i_nlj),
     $                                JpiT(ji/2,ipi,ti/2)
     $                                %st_fi(k_f,ki_kernel)%h_NN(i_nlj),
     $                                0.d0
                              endif
                           endif
                        end do
                     end do
                  endif

               end do
               write(unit_NCSMC,*)
            end do
         end do
      end do
      close(unit_NCSMC)
      end

      subroutine OLS_effective_hamiltonian
      use nodeinfo
      use initst
      use intrface
      use obdens, only: ist1,ist1dim
      implicit none
      logical :: leesuz=.true.
      integer,parameter :: unitolseff=342,unitolsout=341
      real(kind(0.d0)),allocatable :: heff(:,:),root(:),eigv(:,:),
     $     coutr(:,:),coutrtr(:,:),fv1(:),fv2(:)
      integer :: ii,jj,j2min,j2max,kk,num_of_orb,orbmin
      integer,allocatable :: ist_val(:,:)
      integer :: n_i1,l_i1,j_i1,n_i2,l_i2,j_i2,ipi,T12,J12,ind,nu1,nu2
      integer :: i1,i2,jerr,matz,nul1,nul2
      real(kind(0.d0)) :: clebd

      if (iproc/=0) return

      open(unit=unitolsout,file='info_OLS.out',status='unknown',
     $     form='formatted',action='write')
      write(unitolsout,"(' *** OLS effective interaction ***')")
      write(unitolsout,"(/,' A,Z=',2i4,'   full space dim=',i8,
     $'   model space dim=',i6)")
     $     nucleonsi,nprotonsi,nsdi,dim_nhw_mod
      write(unitolsout,
     $     "(' Nmax_tot=',i4,'   model space Nmax_tot=',i4)")
     $     nhwi,nhw_mod
      write(unitolsout,"(' total number of eigenvectors=',i6)") nki
      if (nki<dim_nhw_mod) then
         write(unitolsout,
     $        "('*** error: insufficient number of eigenvectors')")
         print *,
     $        '*** error in OLS: insufficient number of eigenvectors'
         return
      endif

      allocate(heff(dim_nhw_mod,dim_nhw_mod))
      heff=0.d0
      allocate(root(dim_nhw_mod))
      root=0.d0
      allocate(eigv(dim_nhw_mod,dim_nhw_mod))
      eigv=0.d0

      num_of_orb=0
      orbmin=2**30
      do i1=1,ist1dim 
         n_i1=ist1(1,i1)
         l_i1=ist1(2,i1)
         j_i1=ist1(3,i1)
         if (2*n_i1+l_i1/=nhw_mod/2) cycle
         if (2*n_i1+l_i1>nhw_mod/2) exit
         num_of_orb=num_of_orb+1
         if (i1<orbmin) orbmin=i1
      end do

c      print *,' num_of_orb,orbmin=',num_of_orb,orbmin

      j2min=jt2i(1)
      j2max=j2min
      do ii=2,dim_nhw_mod
         if (j2min>jt2i(ii)) j2min=jt2i(ii)
         if (j2max<jt2i(ii)) j2max=jt2i(ii)
      end do

c      print *,' j2min,j2max=',j2min,j2max

      ipi=(-1)**iparityi
      ind=0
      do J12=j2min/2,j2max/2
         do T12=abs(mttotali)/2,1
           do i1=1,ist1dim 
              n_i1=ist1(1,i1)
              l_i1=ist1(2,i1)
              j_i1=ist1(3,i1)
              if (2*n_i1+l_i1/=nhw_mod/2) cycle
              if (2*n_i1+l_i1>nhw_mod/2) exit
              do i2=i1,ist1dim
                 n_i2=ist1(1,i2)
                 l_i2=ist1(2,i2)
                 j_i2=ist1(3,i2)
                 if (2*n_i2+l_i2/=nhw_mod/2) cycle
                 if (2*n_i2+l_i2>nhw_mod/2) exit
                 if ((-1)**(l_i1+l_i2)/=ipi) cycle
                 if (2*n_i1+l_i1+2*n_i2+l_i2/=nhw_mod) cycle
                 if (iabs(j_i1-j_i2)>2*J12) cycle
                 if (j_i1+j_i2<2*J12) cycle
                 if (i1==i2) then
                    if ((-1)**(J12+T12)/=-1) cycle
                 endif
                 ind=ind+1
              end do
           end do
         end do
      end do
      if (ind/=dim_nhw_mod) then
         write(unitolsout,"('*** error: dim_nhw_mod,ind=',2i6)")
     $        dim_nhw_mod,ind
         return
      endif
      allocate(ist_val(4,dim_nhw_mod))
      ind=0
      do J12=j2min/2,j2max/2
         do T12=abs(mttotali)/2,1
           do i1=1,ist1dim 
              n_i1=ist1(1,i1)
              l_i1=ist1(2,i1)
              j_i1=ist1(3,i1)
              if (2*n_i1+l_i1/=nhw_mod/2) cycle
              if (2*n_i1+l_i1>nhw_mod) exit
              do i2=i1,ist1dim
                 n_i2=ist1(1,i2)
                 l_i2=ist1(2,i2)
                 j_i2=ist1(3,i2)
                 if (2*n_i2+l_i2/=nhw_mod/2) cycle
                 if (2*n_i1+l_i1+2*n_i2+l_i2>nhw_mod) exit
                 if ((-1)**(l_i1+l_i2)/=ipi) cycle
                 if (2*n_i1+l_i1+2*n_i2+l_i2/=nhw_mod) cycle
                 if (iabs(j_i1-j_i2)>2*J12) cycle
                 if (j_i1+j_i2<2*J12) cycle
                 if (i1==i2) then
                    if ((-1)**(J12+T12)/=-1) cycle
                 endif
                 ind=ind+1
                 ist_val(1,ind)=i1
                 ist_val(2,ind)=i2
                 ist_val(3,ind)=J12
                 ist_val(4,ind)=T12
              end do
           end do
         end do
      end do
      
c      print *,' ind=',ind

cc      kk=0
cc      do jj=j2min,j2max,2
      do ii=1,dim_nhw_mod
cc            if (jj==jt2i(ii)) then
cc               kk=kk+1
cc              root(kk)=eneri(ii)
         root(ii)=eneri(ii)
         eigv(1:dim_nhw_mod,ii)=bmpi(1:dim_nhw_mod,ii)
cc              eigv(1:dim_nhw_mod,kk)=bmpi(1:dim_nhw_mod,ii)
cc            endif
      end do
cc      end do

      call OLSeff(dim_nhw_mod,dim_nhw_mod,root,eigv,heff,unitolsout,
     $     leesuz)

      write(unitolsout,"(/,' OLSeff called')")

c      print *,' dim_nhw_mod=',dim_nhw_mod
      allocate(coutr(dim_nhw_mod,dim_nhw_mod))
      coutr=0.d0
c      print *,' coutr allocated'
      do ii=1,dim_nhw_mod
         i1=ist_val(1,ii)
         i2=ist_val(2,ii)
         J12=ist_val(3,ii)
         T12=ist_val(4,ii)
         n_i1=ist1(1,i1)
         l_i1=ist1(2,i1)
         j_i1=ist1(3,i1)
         n_i2=ist1(1,i2)
         l_i2=ist1(2,i2)
         j_i2=ist1(3,i2)
         do jj=1,dim_nhw_mod
            select case(mttotali/2)
            case(0)
               nu1=nprotonsi
               nu2=nprotonsi+nneutrnsi
            case(1)
               nu1=nprotonsi-1
               nu2=nprotonsi
            case(-1)
               nu1=nprotonsi+nneutrnsi-1
               nu2=nprotonsi+nneutrnsi
            case default
               write(unitolsout,"('*** error: mttotali=',i4)") mttotali
               return
            end select
            nul1=iloci(nu1,jj)
            nul2=iloci(nu2,jj)
            if ((n_i1==n_spi(nul1).and.l_i1==l_spi(nul1)
     $           .and.j_i1==j2_spi(nul1).and.
     $           n_i2==n_spi(nul2).and.l_i2==l_spi(nul2)
     $           .and.j_i2==j2_spi(nul2)).or.
     $           (n_i1==n_spi(nul2).and.l_i1==l_spi(nul2)
     $           .and.j_i1==j2_spi(nul2).and.
     $           n_i2==n_spi(nul1).and.l_i2==l_spi(nul1)
     $           .and.j_i2==j2_spi(nul1))) then
               coutr(jj,ii)=clebd(j2_spi(nul1),m2_spi(nul1),
     $              j2_spi(nul2),m2_spi(nul2),
     $              2*J12,m2_spi(nul1)+m2_spi(nul2))
     $              *clebd(1,mt2_spi(nul1),1,mt2_spi(nul2),
     $              2*T12,mt2_spi(nul1)+mt2_spi(nul2))
     $              *(krd(j_i1,j2_spi(nul1))*krd(j_i2,j2_spi(nul2))
     $              *krd(l_i1,l_spi(nul1))*krd(l_i2,l_spi(nul2))
     $              *krd(n_i1,n_spi(nul1))*krd(n_i2,n_spi(nul2))
     $              +(-1)**(J12+T12-(j_i1+j_i2)/2)
     $              *krd(j_i1,j2_spi(nul2))*krd(j_i2,j2_spi(nul1))
     $              *krd(l_i1,l_spi(nul2))*krd(l_i2,l_spi(nul1))
     $              *krd(n_i1,n_spi(nul2))*krd(n_i2,n_spi(nul1)))
     $              /sqrt(1.d0
     $              +krd(j_i1,j_i2)*krd(n_i1,n_i2)*krd(l_i1,l_i2))
            endif
c            print *,' ii,jj,coutr=',jj,ii,coutr(jj,ii)
         end do
      end do
c      print *,' coutr calculated'
      eigv=matmul(heff,coutr)
c      print *,' eigv calculated'

      allocate(coutrtr(dim_nhw_mod,dim_nhw_mod))
      coutrtr=0.d0
c      print *,' coutrtr allocated'
      coutrtr=transpose(coutr)
c      print *,' coutrtr calculated'
      heff=matmul(coutrtr,eigv)
c      print *,' heff calculated'

      open(unit=unitolseff,file='Heff_OLS.dat',status='unknown',
     $     form='formatted',action='write')

      write(unitolseff,'(i6,5i4,i10,i4,f8.4)') dim_nhw_mod,nhw_mod,
     $     num_of_orb,nucleonsi,nprotonsi,nneutrnsi,nsdi,nhwi,hboi
      do ii=1,dim_nhw_mod
         i1=ist_val(1,ii)
         i2=ist_val(2,ii)
         J12=ist_val(3,ii)
         T12=ist_val(4,ii)
         n_i1=ist1(1,i1)
         l_i1=ist1(2,i1)
         j_i1=ist1(3,i1)
         n_i2=ist1(1,i2)
         l_i2=ist1(2,i2)
         j_i2=ist1(3,i2)
         write(unitolseff,'(11i4)') ii,i1-orbmin+1,i2-orbmin+1,
     $        n_i1,l_i1,j_i1,n_i2,l_i2,j_i2,J12,T12
      end do
      do ii=1,dim_nhw_mod
         write(unitolseff,'(1x,100(e14.7,1x))')
     $        (heff(jj,ii),jj=1,dim_nhw_mod)
      end do
      close(unit=unitolseff)

      write(unitolsout,"(/,' Heff saved in Heff_OLS.dat')")


      write(unitolsout,"(/,' Heff diagonalization')")
      allocate(fv1(dim_nhw_mod),fv2(dim_nhw_mod))
      matz=1
      CALL rs(dim_nhw_mod,dim_nhw_mod,heff,root,matz,eigv,fv1,fv2,jerr)
      write(unitolsout,"(' Eigenvalues')")
      write(unitolsout,'(8f10.4)') (root(i1),i1=1,dim_nhw_mod)
      write(unitolsout,"(' Ground-state eigenvector')")
      write(unitolsout,'(8f10.4)') (eigv(i1,1),i1=1,dim_nhw_mod) 

      close(unit=unitolsout)

      deallocate(fv1,fv2)
      deallocate(heff,root,eigv,coutr,coutrtr)
      deallocate(ist_val)

      contains

      real(kind(0.d0)) function krd(i,j)
      implicit none
      integer,intent(IN) :: i,j
      if (i==j) then
         krd=1.d0
      else
         krd=0.d0
      endif
      end function krd

      end subroutine OLS_effective_hamiltonian

      subroutine OLSeff(imodel,ifull,root,eigv,heff,iunit,leesuz)
      implicit none
      integer,intent(IN) :: imodel,ifull,iunit
      logical,intent(IN) :: leesuz
      real(kind=kind(0.d0)),intent(IN) :: root(imodel),
     +                                   eigv(ifull,imodel)
      real(kind=kind(0.d0)),intent(OUT) :: heff(imodel,imodel)
      real(kind=kind(0.d0)) :: heffnh(imodel,imodel),
     +                         heffold(imodel,imodel)
      real(kind=kind(0.d0)) :: sum,ovlpmax,tmproot,rootexac(imodel)
      real(kind=kind(0.d0)) :: rootp(imodel),eigvp(imodel,imodel),
     +                      fv1(imodel),fv2(imodel)
      integer :: i1,i2,matz,ierr,l,m,imax,ixchange,liga(imodel),
     +           ligb(imodel)
      logical :: save_wave=.false.

      rootexac=root
      ixchange=0
      do l=1,imodel
         ovlpmax=rootexac(l)
         imax=l
         do m=l+1,imodel
            if (ovlpmax>rootexac(m)) then
               ovlpmax=rootexac(m)
               imax=m
            endif
         end do
         if (imax/=l) then
            ixchange=ixchange+1
            tmproot=rootexac(l)
            rootexac(l)=rootexac(imax)
            rootexac(imax)=tmproot
         endif
      end do
      if (ixchange>0) then
         write(iunit,
     +        '(" *** number of exchanges for re-ordering root=",i3)')
     +        ixchange
      endif

      if (leesuz) then
         heffold(1:imodel,1:imodel)=eigv(1:imodel,1:imodel)
         heffnh=heffold
         call inversg(heffnh,imodel)
cc         call dminv(heffnh,imodel,sum,liga,ligb)
         
         do i1=1,imodel
            do i2=1,imodel
               sum=0.d0
               do l=1,imodel
                  sum=sum+heffnh(i1,l)*heffold(l,i2)
               end do
               if ((i1==i2.and.dabs(sum-1.d0)>1.d-6)
     +              .or.(i1/=i2.and.dabs(sum)>1.d-6)) then
                  write(iunit,'(2i3,f10.6)') i1,i2,sum
               endif
            end do
         end do
         
         do i1=1,imodel
            do i2=1,imodel
               sum=0.d0
               do l=1,imodel
                  sum=sum+heffnh(l,i1)*heffnh(l,i2)
               end do
               eigvp(i1,i2)=sum
            end do
         end do

         write(iunit,2279)
 2279    format(/,' metric diagonalization')
         heffnh=eigvp

         heff=heffnh

         matz=1
         CALL rs(imodel,imodel,heffnh,rootp,matz,eigvp,fv1,fv2,ierr)
         write(iunit,1336) (rootp(i1),i1=1,imodel)
 1336    format(6d12.5)

         do i1=1,imodel
            do i2=1,imodel
               sum=0.d0
               do l=1,imodel
                  sum=sum+heff(i1,l)*eigvp(l,i2)
               end do
               if (dabs(sum-eigvp(i1,i2)*rootp(i2))>1.d-6) then
                  write(iunit,"(' root:',2i4,1x,f12.6)") i1,i2,sum
               endif
            end do
         end do

         do i1=1,imodel
            do i2=1,imodel
               sum=0.d0
               do l=1,imodel
                  sum=sum+eigvp(i1,l)*dsqrt(rootp(l))*eigvp(i2,l)
               end do
               heffnh(i1,i2)=sum
            end do
         end do         
         eigvp=heffnh
         do i1=1,imodel
            do i2=1,imodel
               sum=0.d0
               do l=1,imodel
                  sum=sum+eigvp(i1,l)*heffold(l,i2)
               end do
               heffnh(i1,i2)=sum
            end do
         end do         

      else
         do l=1,imodel
            fv1(1:imodel)=eigv(1:imodel,l)
            do m=1,l-1
               fv2(1:imodel)=heffnh(1:imodel,m)
               sum=0.d0
               do i1=1,imodel
                  sum=sum+fv1(i1)*fv2(i1)
               end do
               do i1=1,imodel
                  fv1(i1)=fv1(i1)-sum*fv2(i1)
               end do
            end do
            sum=0.d0
            do i1=1,imodel
               sum=sum+fv1(i1)*fv1(i1)
            end do
            sum=1.d0/dsqrt(sum)
            fv1=fv1*sum
            heffnh(:,l)=fv1(:)
         end do
      endif
      
      do i1=1,imodel
         do i2=1,imodel
            sum=0.d0
            do l=1,imodel
               sum=sum+heffnh(l,i1)*heffnh(l,i2)
            end do
            if ((i1==i2.and.dabs(sum-1.d0)>1.d-4)
     +           .or.(i1/=i2.and.dabs(sum)>1.d-4)) then
               write(iunit,"(' ortho:',2i4,1x,f12.6)") i1,i2,sum
            endif
         end do
      end do

      do i1=1,imodel
         do i2=1,imodel      
            sum=0.d0
            do l=1,imodel
               sum=sum+heffnh(i1,l)*root(l)*heffnh(i2,l)
            end do
            heff(i1,i2)=sum
         end do
      end do

      sum=0.d0
      do i1=1,imodel
         do i2=i1,imodel
            sum=sum+dabs(heff(i1,i2)-heff(i2,i1))
         end do
      end do
      write (iunit,4561) sum
 4561 format(/,' non-hermiticity:',d14.7)      

      heffold=heff
      if (save_wave) then
         matz=imodel
      else
         matz=1
      endif
      CALL rs(imodel,imodel,heffold,rootp,matz,eigvp,fv1,fv2,ierr)

      if (leesuz) then
         write(iunit,2678)
 2678    format(/,
     +' model-space diagonalization:
     $Okubo-Lee-Suzuki effective interaction')
      else
         write(iunit,2679)
 2679    format(/,
     + ' model-space diagonalization: CORE effective interaction')
      endif
      write(iunit,1235) (rootp(i1),i1=1,imodel)
      write(iunit,1235) (eigvp(i1,1),i1=1,imodel) 
 1235 format(8f10.4)
      if (save_wave) then
         write(49,'(6f12.6)') (rootp(i1),i1=1,imodel)
         write(49,*)
         do i2=1,imodel
            write(49,'(6f13.7)') (eigvp(i1,i2),i1=1,imodel)
         end do
         write(49,*)
      endif

c**      do i2=1,imodel
c**         sum=0.d0
c**         do i1=1,imodel
c**            sum=sum+dabs(dabs(eigvp(i1,i2))-dabs(heffnh(i1,i2)))
c**         end do
c**         if (sum>1.d-6) then
c**            write(iunit,1236) rootp(i2)
c**            write(iunit,1236) (eigvp(i1,i2),heffnh(i1,i2),i1=1,imodel)
c** 1236       format(2f10.4)
c**            write(iunit,*)
c**         endif
c**      end do

      sum=0.d0
      do i1=1,imodel
         sum=sum+dabs(rootexac(i1)-rootp(i1))
      end do
      write(iunit,3451) sum
 3451 format(/,' diff=',f10.5,/)

      end subroutine OLSeff

      subroutine testread(snsp,sndp,dnuc)
      use obdens
      use initst
      use finast
      use intrface
      use constants
      use paramdef
      use gamalog
      use harmosc, only: bsquare
      implicit none
      logical mfdp,snsp,sndp,dnuc,obscal,rad2b,rad2beff,rad2beffx,e2_2b
      integer(4) nhme,k1max,mxsps,major
      real(kind=kind(0.d0)) xxx,cleb,clb,clb1,sumpl,sumnl,sumps,sumns
      real(kind=kind(0.d0)) dsqj,sumto,clebd,oeff1,sumf,sumgt,racad
      real(kind=kind(0.d0)),allocatable:: sumN(:,:),occup(:,:),
     +     sumrule(:,:)
      real(kind=kind(0.d0)) ELp,ELn,obme,clebt,E1sp,E1sn,E1s,sumco,
     $     Ecore,Rcore,Eobd,Robd
      integer jtrans,jtrans2,jtotali2,jtotalf2,kii,kff
      integer jtot2,jtotal2,jtotal2min,itotali2,itotalf2,jtx2
      integer jtransin,N1,N2,N1in,N2in,j1,j2,k1,ii,jtrmin,jtrmax
      integer jtot2min,kiiin,kffin,l1,l2,i1,i3,ipi,i2
      integer nonebody,nonebody2,nonebody3,nonebody4,l3,l4,isum
      integer n3,n4,k2,j3,j4,ob1,ob2,ob3,ob4,k,ids,ip1,ip2,iq
      real(kind=kind(0.d0)),allocatable :: ham(:,:,:),tkin(:,:),
     +                                     vcoul(:,:)
      real(kind=kind(0.d0)),allocatable :: rreleff(:,:,:),rrel(:,:),
     +     e2_op_2b(:,:,:),spin_spin(:,:,:,:),Esp(:,:),Rsp(:,:)
      integer :: j12,j34,it12,it34,i1x,i2x,i3x,i4x,ixxx,numsp,numspR,
     +     n_ii,l_ii,j2_ii,n_ij,l_ij,j2_ij,ij,ind,J_e2_2b,T_e2_2b
      real(kind=kind(0.d0)) vpn,vpp,vnn,trel,hrel,coul,rel,releff
      real(kind=kind(0.d0)) rad0,radpp,radeff0,radeffpp,radmass,
     +     radneut,radprot,releffpp,releffnn,twobsum,e2_op_sum,
     +     pht12,pht34,t20,t21,t2m1,SpSp,SpSn,SnSn,releffnp
      character(len=2) :: text
      integer :: nhom12max,Jtr_mi,Jtr_ma,mT_mi,mT_ma,in_mi,in_ma,mTt,
     $     J12min,J12max,ind_pi,J12ma
      integer :: i1x_r,i2x_r,i3x_r,i4x_r,j34_r,it34_r

      open(iunitvout,file='trdens.out',status='old',action='read')
      open(iunitobout,file='observ.out',status='unknown')
      read(iunitvout,*)
      read(iunitvout,*) mfdp
      read(iunitvout,*)

      allocate(nsdi,mxnwdi,nhomi,mjtotali,mttotali,iparityi,nhwi,
     +        nucleonsi,nprotonsi,nneutrnsi,hboi,nhom12i,naspsi)

      sndp=.false.
      dnuc=.false.

      if (mfdp) then

         snsp=.true.

         read(iunitvout,*)
         read(iunitvout,*)
         read(iunitvout,1001) nucleonsi,nprotonsi,nneutrnsi,mjtotali,
     + mttotali,hboi,nhwi,nsdi,nhme,k1max,
     + mxnwdi,mxsps,major,iparityi 
 1000    format(' Nucleus:',/,' A=',i3,'   Z=',i3,'   N=',i3,
     + /,' 2*MJ=',i3,'   2*MT=',i3,'  parity= +',/,
     + ' hbar Omega=',f8.4,'   Nhw=',i3,'   dimension=',i8,'   nhme=',
     + i10,/,' k1max=',i3,'   mxnwd=',i3,'   mxsps=',i8,
     +                    '   major=',i2,'   iparity=',i2,/)
 1007    format(' Nucleus:',/,' A=',i3,'   Z=',i3,'   N=',i3,
     + /,' 2*MJ=',i3,'   2*MT=',i3,'  parity= -',/,
     + ' hbar Omega=',f8.4,'   Nhw=',i3,'   dimension=',i8,'   nhme=',
     + i10,/,' k1max=',i3,'   mxnwd=',i3,'   mxsps=',i8,
     +                    '   major=',i2,'   iparity=',i2,/)
 1001    format(9x,/,3x,i3,5x,i3,5x,i3,
     + /,6x,i3,8x,i3,11x,/,
     + 12x, f8.4, 7x, i3, 13x, i8, 8x,
     + i10,/, 7x, i3, 9x, i3, 9x, i8,
     + 9x, i2, 11x, i2,/)

         if (iparityi==0) then
         write(iunitobout,1000) nucleonsi,nprotonsi,nneutrnsi,mjtotali,
     + mttotali,hboi,nhwi,nsdi,nhme,k1max,
     + mxnwdi,mxsps,major,iparityi 
         else
         write(iunitobout,1007) nucleonsi,nprotonsi,nneutrnsi,mjtotali,
     + mttotali,hboi,nhwi,nsdi,nhme,k1max,
     + mxnwdi,mxsps,major,iparityi 
         endif
         allocate(Jxi(k1max+2),Txi(k1max+2),jt2i(k1max+2),
     +        it2i(k1max+2),
     +        eneri(k1max+2))
         read(iunitvout,1101) (Jxi(k1),Txi(k1),eneri(k1),xxx,
     + k1=1,k1max+2)         
 1100    format(' J=',f7.4,'    T=',f7.4,'     Energy=',f12.4,
     + '     Ex=',f12.4) 
 1101    format(3x, f7.4, 6x, f7.4, 12x, f12.4, 8x, f12.4) 

         write(iunitobout,1100) (Jxi(k1),Txi(k1),eneri(k1),
     +        eneri(k1)-eneri(1),k1=1,k1max+2)         

         do ixxx=1,k1max+2
            if (mod(nucleonsi,2)==0) then
               jt2i(ixxx)=2*nint(Jxi(ixxx))
               it2i(ixxx)=2*nint(Txi(ixxx))
            else
               jt2i(ixxx)=nint(2*Jxi(ixxx))
               if (mod(jt2i(ixxx),2)/=1) then
                  if (2*Jxi(ixxx)>jt2i(ixxx)) then
                     jt2i(ixxx)=nint(2*Jxi(ixxx))+1
                  else
                     jt2i(ixxx)=nint(2*Jxi(ixxx))-1
                  endif
               endif
               it2i(ixxx)=nint(2*Txi(ixxx))
               if (mod(it2i(ixxx),2)/=1) then
                  if (2*Txi(ixxx)>it2i(ixxx)) then
                     it2i(ixxx)=nint(2*Txi(ixxx))+1
                  else
                     it2i(ixxx)=nint(2*Txi(ixxx))-1
                  endif
               endif
            endif
            if (dabs(dble(jt2i(ixxx))-2.d0*Jxi(ixxx))>1.d-4) then
               print *, ixxx,'**** not a good J ****'
            endif 
         end do

         read(iunitvout,1157) nhomi,nhom12i,naspsi
 1156    format(/,' N1_max=',i4,'   N12_max=',i4,'   Nasps=',i4)
 1157    format(/, 8x, i4,11x,i4,9x,i4)
         write(iunitobout,1156) nhomi,nhom12i,naspsi

         read(iunitvout,1348) ki,nki
 1347    format(/,' wave functions of the states #',i3,'- #',i3,
     +        ' used')
 1348    format(/, 31x, i3, 3x, i3, 5x)
         nki=nki+1-ki
         write(iunitobout,1347) ki,nki

         nsdf=>nsdi
         mxnwdf=>mxnwdi
         nhomf=>nhomi
         nhom12f=>nhom12i
         naspsf=>naspsi
         mjtotalf=>mjtotali
         mttotalf=>mttotali
         iparityf=>iparityi
         nucleonsf=>nucleonsi
         nprotonsf=>nprotonsi
         nneutrnsf=>nneutrnsi
         jt2f=>jt2i
         it2f=>it2i
         Jxf=>Jxi
         Txf=>Txi
         hbof=>hboi
         enerf=>eneri

      else

         snsp=.false.

         allocate(nsdf,mxnwdf,nhomf,mjtotalf,mttotalf,iparityf,nhwf,
     +        nucleonsf,nprotonsf,nneutrnsf,hbof,nhom12f,naspsf)

         read(iunitvout,*) 
         read(iunitvout,*) 
         read(iunitvout,*) 
         read(iunitvout,*) 

         read(iunitvout,*)
         read(iunitvout,*)
         read(iunitvout,1001) nucleonsi,nprotonsi,nneutrnsi,mjtotali,
     +        mttotali,hboi,nhwi,nsdi,nhme,k1max,
     +        mxnwdi,mxsps,major,iparityi 

         if (iparityi==0) then
         write(iunitobout,1000) nucleonsi,nprotonsi,nneutrnsi,mjtotali,
     + mttotali,hboi,nhwi,nsdi,nhme,k1max,
     + mxnwdi,mxsps,major,iparityi 
         else
         write(iunitobout,1007) nucleonsi,nprotonsi,nneutrnsi,mjtotali,
     + mttotali,hboi,nhwi,nsdi,nhme,k1max,
     + mxnwdi,mxsps,major,iparityi 
         endif
         allocate(Jxi(k1max+2),Txi(k1max+2),jt2i(k1max+2),
     +        it2i(k1max+2),
     +        eneri(k1max+2))
         read(iunitvout,1101) (Jxi(k1),Txi(k1),eneri(k1),xxx,
     +        k1=1,k1max+2)         

         write(iunitobout,1100) (Jxi(k1),Txi(k1),eneri(k1),
     +        eneri(k1)-eneri(1),k1=1,k1max+2)         

         do ixxx=1,k1max+2
            if (mod(nucleonsi,2)==0) then
               jt2i(ixxx)=2*nint(Jxi(ixxx))
               it2i(ixxx)=2*nint(Txi(ixxx))
            else
               jt2i(ixxx)=nint(2*Jxi(ixxx))
               if (mod(jt2i(ixxx),2)/=1) then
                  if (2*Jxi(ixxx)>jt2i(ixxx)) then
                     jt2i(ixxx)=nint(2*Jxi(ixxx))+1
                  else
                     jt2i(ixxx)=nint(2*Jxi(ixxx))-1
                  endif
               endif
               it2i(ixxx)=nint(2*Txi(ixxx))
               if (mod(it2i(ixxx),2)/=1) then
                  if (2*Txi(ixxx)>it2i(ixxx)) then
                     it2i(ixxx)=nint(2*Txi(ixxx))+1
                  else
                     it2i(ixxx)=nint(2*Txi(ixxx))-1
                  endif
               endif
            endif
            if (dabs(dble(jt2i(ixxx))-2.d0*Jxi(ixxx))>1.d-4) then
               print *, ixxx,'**** not a good J ****'
            endif 
         end do

         read(iunitvout,1157) nhomi,nhom12i,naspsi
         write(iunitobout,1156) nhomi,nhom12i,naspsi

         read(iunitvout,1348) ki,nki
         nki=nki+1-ki
         write(iunitobout,1347) ki,nki

         read(iunitvout,*) 
         read(iunitvout,*) 
         read(iunitvout,*) 
         read(iunitvout,*) 

         read(iunitvout,*)
         read(iunitvout,*)
         read(iunitvout,1001) nucleonsf,nprotonsf,nneutrnsf,mjtotalf,
     +        mttotalf,hbof,nhwf,nsdf,nhme,k1max,
     +        mxnwdf,mxsps,major,iparityf 

         if (iparityf==0) then
         write(iunitobout,1000) nucleonsf,nprotonsf,nneutrnsf,mjtotalf,
     +        mttotalf,hbof,nhwf,nsdf,nhme,k1max,
     +        mxnwdf,mxsps,major,iparityf 
         else
         write(iunitobout,1007) nucleonsf,nprotonsf,nneutrnsf,mjtotalf,
     +        mttotalf,hbof,nhwf,nsdf,nhme,k1max,
     +        mxnwdf,mxsps,major,iparityf 
         endif

         allocate(Jxf(k1max+2),Txf(k1max+2),jt2f(k1max+2),
     +        it2f(k1max+2),
     +        enerf(k1max+2))
         read(iunitvout,1101) (Jxf(k1),Txf(k1),enerf(k1),xxx,
     +        k1=1,k1max+2)         

         write(iunitobout,1100) (Jxf(k1),Txf(k1),enerf(k1),
     +        enerf(k1)-enerf(1),k1=1,k1max+2)         

         do ixxx=1,k1max+2
            if (mod(nucleonsf,2)==0) then
               jt2f(ixxx)=2*nint(Jxf(ixxx))
               it2f(ixxx)=2*nint(Txf(ixxx))
            else
               jt2f(ixxx)=nint(2*Jxf(ixxx))
               if (mod(jt2f(ixxx),2)/=1) then
                  if (2*Jxf(ixxx)>jt2f(ixxx)) then
                     jt2f(ixxx)=nint(2*Jxf(ixxx))+1
                  else
                     jt2f(ixxx)=nint(2*Jxf(ixxx))-1
                  endif
               endif
               it2f(ixxx)=nint(2*Txf(ixxx))
               if (mod(it2f(ixxx),2)/=1) then
                  if (2*Txf(ixxx)>it2f(ixxx)) then
                     it2f(ixxx)=nint(2*Txf(ixxx))+1
                  else
                     it2f(ixxx)=nint(2*Txf(ixxx))-1
                  endif
               endif
            endif
            if (dabs(dble(jt2f(ixxx))-2.d0*Jxf(ixxx))>1.d-4) then
               print *, ixxx,'**** not a good J ****'
            endif 
         end do

         read(iunitvout,1157) nhomf,nhom12f,naspsf
         write(iunitobout,1156) nhomf,nhom12f,naspsf

         nhomi=max(nhomi,nhomf)
         nhomf=>nhomi
         nhom12i=max(nhom12i,nhom12f)
         nhom12f=>nhom12i

         select case(iabs(nucleonsi-nucleonsf))
      case(0)
         if (nprotonsi==nprotonsf) then
            if (iparityi/=iparityf) sndp=.true.
            if (iparityi==iparityf.and.nsdi==nsdf) snsp=.true.
         elseif (iabs(nprotonsi-nprotonsf)==1) then   
            dnuc=.true.
         else
            write(iunitobout,*) '*** incompatible number of nucleons'
            stop
         endif   
      case default
         write(iunitobout,*) 
     +        ' nucleonsi and nucleonsf out of range',
     +        nucleonsi,nucleonsf
         stop
         end select

      endif
   
      read(iunitvout,1348) kf,nkf
      nkf=nkf+1-kf
      write(iunitobout,1347) kf,nkf

      jtotal2=(jt2i(ki)+jt2f(kf))/2
      jtotal2min=iabs(jt2i(ki)-jt2f(kf))/2
      do kii=ki,ki+nki-1
         do kff=kf,kf+nkf-1
            jtot2min=iabs(jt2i(kii)-jt2f(kff))/2
            jtot2=(jt2i(kii)+jt2f(kff))/2
            if (jtotal2min>jtot2min) jtotal2min=jtot2min
            if (jtotal2<jtot2) jtotal2=jtot2
         end do
      end do 

      jtotal2=min(jtotal2,jtotal2max)

      print *,' jtotal2min,jtotal2:',jtotal2min,jtotal2

      read(iunitvout,*)
      write(iunitobout,*)
      read(iunitvout,7769) ist1dim
      write(iunitobout,7768) ist1dim
 7768 format(' number of single-nucleon states =',i4)
 7769 format(34x,i4)
      allocate(ist1(3,ist1dim))
      do i1=1,ist1dim
         read(iunitvout,7771) ii,(ist1(i3,i1),i3=1,3)
 7771    format(2x,i4,4x,i3,4x,i3,4x,i2,2x)
         if (ii/=i1) then
            print *,'*** error: i1,ii=',i1,ii
            stop
         endif   
      end do

      allocate(tdJp(ist1dim,ist1dim,jtotal2min:jtotal2,
     +     kf:kf+nkf-1,ki:ki+nki-1))          
      allocate(tdJn(ist1dim,ist1dim,jtotal2min:jtotal2,
     +     kf:kf+nkf-1,ki:ki+nki-1))          
      tdJp=0.d0
      tdJn=0.d0

      do kii=ki,ki+nki-1

         jtotali2=jt2i(kii)
         itotali2=it2i(kii)
         if (dabs(dble(jtotali2)-2.d0*Jxi(kii))>1.d-4) then
            print *, '**** not a good J ****'
         endif 

         do kff=kf,kf+nkf-1

            jtotalf2=jt2f(kff)
            itotalf2=it2f(kff)
            if (dabs(dble(jtotalf2)-2.d0*Jxf(kff))>1.d-4) then
               print *, '**** not a good J ****'
            endif 

            read(iunitvout,1202) kffin,jtotalf2,itotalf2,xxx,
     +           kiiin,jtotali2,itotali2,sumto
 1200       format(//,' *** Transition matrix elements for states:',
     +           ' ***',/,' #',i3,
     +           ' [2*(J,T),Ex]_f=',i3,i2,f8.4,3x,
     +           '#',i3,' [2*(J,T),Ex]_i=',i3,i2,f8.4) 
 1202       format(//, 47x,/, 2x, i3,
     +           16x, i3, i2, f8.4, 3x,
     +           1x, i3, 16x, i3, i2, f8.4) 
c* 1200       format(//,' *** Transition matrix elements for states:',
c*     +           ' ***',/,' #',i2,
c*     +           '  [2*(J,T),Ex]_f=',2i3,f8.4,4x,
c*     +           '#',i2,'  [2*(J,T),Ex]_i=',2i3,f8.4) 
c* 1200       format(//,' *** Transition matrix elements for states:',
c*     +           ' ***',/,' #',i2,
c*     +           '  2*(J,T)_f=',2i3,'  Ex_f=',f8.4,4x,
c*     +           '#',i2,'   2*(J,T)_i=',2i3,'  Ex_i=',f8.4) 
c* 1200       format(//,' *** Transition matrix elements for states:',
c*     +           ' ***',/,' #',i3,
c*     +           '   2*J_f,2*T_f=',2i3,'   Ex_f=',f9.4,5x,
c*     +           ' #',i3,'   2*J_i,2*T_i=',2i3,'   Ex_i=',f9.4)

            if (kffin/=kff.or.kiiin/=kii) then
               print *,'*** error: kiiin,kffin=',kiiin,kffin,
     +              '  kii,kff=',kii,kff
               stop
            endif   

            do jtrans=jtotal2min,jtotal2

               jtrans2=jtrans+jtrans
               if (iabs(jt2i(kii)-jt2f(kff))>jtrans2
     +              .or.(jt2i(kii)+jt2f(kff))<jtrans2) cycle
               
               cleb=clebd(jtotali2,mjtotali,jtrans2,mjtotalf-mjtotali,
     +              jtotalf2,mjtotalf)

               if (dabs(cleb)<1.d-8) cycle

               read(iunitvout,1954) jtransin
 1954          format(/, 8x, i3)

               if (jtransin/=jtrans) then
                  print *,'*** error: jtransin=',jtransin,
     +                 '  jtrans=',jtrans
                  stop
               endif   

               do i1=1,ist1dim
                  n1=ist1(1,i1)
                  l1=ist1(2,i1)
                  j1=ist1(3,i1)
                  do i3=1,ist1dim
                     n2=ist1(1,i3)
                     l2=ist1(2,i3)
                     j2=ist1(3,i3)

                     if((j1+j2)<2*jtrans) cycle
                     if(iabs(j1-j2)>2*jtrans) cycle

                     if (((-1)**l1==(-1)**l2.and.
     +                    iparityi==iparityf).or.
     +           ((-1)**l1/=(-1)**l2.and.iparityi/=iparityf)) then
                        read(iunitvout,1962) N1in,N2in,
     +                       tdJp(i1,i3,jtrans,kff,kii),
     +                       tdJn(i1,i3,jtrans,kff,kii)
 1962                   format(i4,1x,i4,1x,e15.8,1x,e15.8)
c 1960                   format(' a+=',i3,'    a-=',i3,
c     +                      '     td(a+,a-): p=',f10.6,'     n=',f10.6)
c 1961                   format(4x, i3, 7x, i3,
c     +                      18x, f10.6, 7x, f10.6)
                        if (N1in/=i1.or.N2in/=i3) then
                           print *,'*** error: N1in,N2in=',
     +                                N1in,N2in,'  i1,i3=',i1,i3
                           stop
                        endif   
                     endif
                  end do
               end do    
            end do   
         end do
      end do

      if (twobdcal) then
         jtrmin=0
         jtrmax=max(4,jtotal2)
      else
         jtrmin=jtotal2min
         jtrmax=jtotal2
      endif

      if (twobdcal) then
         read(iunitvout,*)
         read(iunitvout,*)
         read(iunitvout,*) 
         read(iunitvout,*)
         read(iunitvout,7777) ist2dim
         write(iunitobout,7778) ist2dim
 7777    format(31x,i5)
 7778    format(' number of two-nucleon states =',i5)
         read(iunitvout,7877) isum2
         write(iunitobout,7878) isum2
 7877    format(49x,i7)
 7878    format(' number of two-body Hamiltonian matrix elements =',
     +           i7)
         allocate(ist2(4,ist2dim))
         do i1=1,ist2dim
            read(iunitvout,7779) ii,(ist2(i3,i1),i3=1,4)
 7779       format(2x,i5,10x,3i3,i2)
            if (ii/=i1) then
               print *,'*** error: i1,ii=',i1,ii
               stop
            endif   
         end do

         if (allocated(sumrule)) deallocate(sumrule)
         allocate(sumrule(jtrmin:jtrmax,0:3))

         print *,'ist2dim,jtotal2min,jtotal2,kf,kf+nkf-1,ki,ki+nki-1'
         print *,ist2dim,jtotal2min,jtotal2,kf,kf+nkf-1,ki,ki+nki-1

c         allocate(t2bd(ist2dim,ist2dim,jtotal2min:jtotal2,
c     +          -1:1,kf:kf+nkf-1,ki:ki+nki-1))          
c         t2bd=0.d0

         nhom12max=max(nhom12i,nhom12f)
         print *,' nhom12max=',nhom12max
         J12min=0
         J12max=nhom12max+1
         ist2_Jst_J12min=J12min
         ist2_Jst_J12max=J12max
         allocate(ist2_Jst(0:1,J12min:J12max+1))
         ind=0
         do ipi=1,-1,-2
c            ind_pi=0
            if (ipi/=(-1)**nhom12max) then
               J12ma=J12max-1
            else
               J12ma=J12max
            endif
            do J12=J12min,J12ma
               ist2_Jst((ipi+1)/2,J12)=ind+1
c            print *,' J12,ist2_Jst=',J12,ist2_Jst(J12)
               do iT12=0,1
                  do ii=1,ist1dim
                     n_ii=ist1(1,ii)
                     l_ii=ist1(2,ii)
                     j2_ii=ist1(3,ii)
                     do ij=ii,ist1dim
                        n_ij=ist1(1,ij)
                        l_ij=ist1(2,ij)
                        j2_ij=ist1(3,ij)
                        if ((-1)**(l_ii+l_ij)/=ipi) cycle
                        if (2*n_ii+l_ii+2*n_ij+l_ij>nhom12max) cycle
                        if (iabs(j2_ii-j2_ij)>2*J12) cycle
                        if (j2_ii+j2_ij<2*J12) cycle
                        if (ii==ij) then
                           if ((-1)**(J12+iT12)/=-1) cycle
                        endif
                        ind=ind+1
c                        ind_pi=ind_pi+1
                     end do
                  end do   
               end do   
            end do
            ist2_Jst((ipi+1)/2,J12ma+1)=ind+1   
         end do   

         allocate(tbd(kf:kf+nkf-1,ki:ki+nki-1))
         do kii=ki,ki+nki-1
            do kff=kf,kf+nkf-1
               Jtr_mi=abs(jt2i(kii)-jt2f(kff))/2
               Jtr_ma=min((jt2i(kii)+jt2f(kff))/2,jtotal2max)
               tbd(kff,kii)%Jtr_min=Jtr_mi
               tbd(kff,kii)%Jtr_max=Jtr_ma
               allocate(tbd(kff,kii)%Jtr(Jtr_mi:Jtr_ma))
               do jtrans=Jtr_mi,Jtr_ma
                  tbd(kff,kii)%Jtr(jtrans)%dim_fi=ist2dim
                  allocate(tbd(kff,kii)%Jtr(jtrans)%fi(ist2dim))
                  do ii=1,ist2dim
                     i1=ist2(1,ii)
                     i2=ist2(2,ii)
                     j12=ist2(3,ii)
                     it12=ist2(4,ii)
                     mT_mi=max(-abs(it12),(mttotalf-mttotali)/2-1)
                     mT_ma=min(abs(it12),(mttotalf-mttotali)/2+1)
                     tbd(kff,kii)%Jtr(jtrans)%fi(ii)%mT_min=mT_mi
                     tbd(kff,kii)%Jtr(jtrans)%fi(ii)%mT_max=mT_ma
                     if (mT_mi>mT_ma) cycle
                     l1=ist1(2,i1)
                     l2=ist1(2,i2)
                     ipi=(-1)**(l1+l2+iparityi+iparityf)
                     if (ipi/=(-1)**nhom12max) then
                        J12ma=nhom12max
                     else
                        J12ma=nhom12max+1
                     endif
                     if (abs(jtrans-j12)>ist2_Jst_J12max
     $                    .or.abs(jtrans-j12)<ist2_Jst_J12min
     $                    .or.min(jtrans+j12,J12ma)+1>ist2_Jst_J12max+1
     $                    .or.min(jtrans+j12,J12ma)+1<ist2_Jst_J12min)
     $                    then
                        in_mi=-1
                        in_ma=-1
                     else
                        in_mi=ist2_Jst((ipi+1)/2,abs(jtrans-j12))
                        in_ma=ist2_Jst((ipi+1)/2,
     $                       min(jtrans+j12,J12ma)+1)-1
                     endif                     
c                     in_mi=ist2_Jst((ipi+1)/2,abs(jtrans-j12))
c                     in_ma=ist2_Jst((ipi+1)/2,min(jtrans+j12,J12ma)+1)-1
                     allocate(tbd(kff,kii)%Jtr(jtrans)%fi(ii)
     $                    %mT(mT_mi:mT_ma))
                     do mTt=mT_mi,mT_ma 
                        tbd(kff,kii)%Jtr(jtrans)%fi(ii)%mT(mTt)
     $                       %in_min=in_mi
                        tbd(kff,kii)%Jtr(jtrans)%fi(ii)%mT(mTt)
     $                       %in_max=in_ma
                        allocate(tbd(kff,kii)%Jtr(jtrans)%fi(ii)%mT(mTt)
     $                       %in(in_mi:in_ma))
                        tbd(kff,kii)%Jtr(jtrans)%fi(ii)%mT(mTt)
     $                       %in=0.d0
                     end do
                  end do
               end do
            end do
         end do
         print *,' tbd allocated'

         do kii=ki,ki+nki-1

            jtotali2=jt2i(kii)
            itotali2=it2i(kii)
            if (dabs(dble(jtotali2)-2.d0*Jxi(kii))>1.d-4) then
               print *, '**** not a good J ****'
            endif 

            do kff=kf,kf+nkf-1

               jtotalf2=jt2f(kff)
               itotalf2=it2f(kff)
               if (dabs(dble(jtotalf2)-2.d0*Jxf(kff))>1.d-4) then
                  print *, '**** not a good J ****'
               endif 

               read(iunitvout,1202) kffin,jtotalf2,itotalf2,xxx,
     +              kiiin,jtotali2,itotali2,sumto

               if (kffin/=kff.or.kiiin/=kii) then
                  print *,'*** error: kiiin,kffin=',kiiin,kffin,
     +                 '  kii,kff=',kii,kff
                  stop
               endif   

               do jtrans=jtotal2min,jtotal2

                  if (iabs(jt2i(kii)-jt2f(kff))>2*jtrans
     +                 .or.(jt2i(kii)+jt2f(kff))<2*jtrans) cycle
               
                  cleb=clebd(jtotali2,mjtotali,2*jtrans,
     +                 mjtotalf-mjtotali,jtotalf2,mjtotalf)

*               if (jtrans>2) print *,' jtotali2,jtotalf2,jtrans:',
*     +              jtotali2,jtotalf2,jtrans,' cleb=',cleb

                  if (dabs(cleb)<1.d-8) cycle

                  read(iunitvout,1954) jtransin
                  if (jtransin/=jtrans) then
                     print *,'*** error: jtransin=',jtransin,
     +                    '  jtrans=',jtrans
                     stop
                  endif   

                  cleb=dsqrt(dble(jtotalf2)+1.d0)/cleb 
c*** Glendenning definition
                  clb1=-cleb/dsqrt(dble(2*jtrans+1)) 
c*** Alternative definition
c*               clb1=-1.d0               

*                if (jtrans>2) write(2,*) ' clb1=',clb1

                  print *,' kii,kff,jtrans=',kii,kff,jtrans                  

                  do
                     read(iunitvout,'(a)',end=12366) text
                     backspace iunitvout
cc                     if (text==' (a+a+)') then
                     if (text=='tb') then
                        read(iunitvout,1965) i1,i3,t20,t21,t2m1
 1965                   format(2x,i6,1x,i6,1x,e15.8,1x,e15.8,1x,e15.8)
c 1966                   format(9x,i5,10x,i5,10x,f10.6,6x,f10.6,
c     +                       6x,f10.6)
                        tbd(kff,kii)%Jtr(jtrans)%fi(i1)
     $                             %mT(0)%in(i3)=t20

c                        print *,' i1,i3,tbd=',i1,i3,
c     $                       tbd(kff,kii)%Jtr(jtrans)%fi(i1)
c     $                          %mT(0)%in(i3)

                        it12=ist2(4,i1)
                        if (it12==1) then
                           tbd(kff,kii)%Jtr(jtrans)%fi(i1)
     $                          %mT(1)%in(i3)=t21
                           tbd(kff,kii)%Jtr(jtrans)%fi(i1)
     $                          %mT(-1)%in(i3)=t2m1

c                        print *,' i1,i3,tbd=',i1,i3,
c     $                       tbd(kff,kii)%Jtr(jtrans)%fi(i1)
c     $                          %mT(1)%in(i3),
c     $                       tbd(kff,kii)%Jtr(jtrans)%fi(i1)
c     $                          %mT(-1)%in(i3)

                        endif

c                        t2bd(i1,i3,jtrans,0,kff,kii)=t20
c                        t2bd(i1,i3,jtrans,1,kff,kii)=t21
c                        t2bd(i1,i3,jtrans,-1,kff,kii)=t2m1
                     else
                        exit
                     endif
                  end do

c                  do i1=1,ist2dim
c                     nonebody=ist2(1,i1)
c                     nonebody2=ist2(2,i1)
c                     j12=ist2(3,i1)
c                     it12=ist2(4,i1)
c
c                     l1=ist1(2,nonebody)
c                     l2=ist1(2,nonebody2)
c
c                     do i3=1,ist2dim
c                        nonebody3=ist2(1,i3)
c                        nonebody4=ist2(2,i3)
c                        j34=ist2(3,i3)
c                        it34=ist2(4,i3)
c
c                        if((j12+j34)<jtrans) cycle
c                        if(iabs(j12-j34)>jtrans) cycle
c
c                        l3=ist1(2,nonebody3)
c                        l4=ist1(2,nonebody4)
c
c                        if (((-1)**(l1+l2+l3+l4)==1.and.
c     +                       iparityi==iparityf).or.
c     +                       ((-1)**(l1+l2+l3+l4)==-1
c     +                       .and.iparityi/=iparityf)) then
c                           read(iunitvout,1966) N1in,N2in,
c     +                          t2bd(i1,i3,jtrans,0,kff,kii),
c     +                          t2bd(i1,i3,jtrans,1,kff,kii),
c     +                          t2bd(i1,i3,jtrans,-1,kff,kii)
c 1966                      format(9x,i4,10x,i4,10x,f10.6,6x,f10.6,
c     +                          6x,f10.6)
c                           if (N1in/=i1.or.N2in/=i3) then
c                              print *,'*** error: N1in,N2in=',
c     +                             N1in,N2in,'  i1,i3=',i1,i3
c                              stop
c                           endif   
c                        endif
c                     end do
c                  end do

               end do
            end do
         end do

      endif   
12366 continue
      close(iunitvout)

      ipi=((-1)**iparityi)*((-1)**iparityf)
      call rLmatel(hboi,nhomi,jtrmin,jtrmax,ipi)

c      if (twobdcal.and.snsp) then
      if (twobdcal) then

         allocate(itbind(ist1dim*(ist1dim+1)/2,0:J12max,0:1))
         itbind=0
         ind=0
         do ipi=1,-1,-2
            if (ipi/=(-1)**nhom12max) then
               J12ma=J12max-1
            else
               J12ma=J12max
            endif
            do J12=0,J12ma
               do iT12=0,1
                  do ii=1,ist1dim
                     n_ii=ist1(1,ii)
                     l_ii=ist1(2,ii)
                     j2_ii=ist1(3,ii)
                     do ij=ii,ist1dim
                        n_ij=ist1(1,ij)
                        l_ij=ist1(2,ij)
                        j2_ij=ist1(3,ij)
                        if ((-1)**(l_ii+l_ij)/=ipi) cycle
                        if (2*n_ii+l_ii+2*n_ij+l_ij>nhom12i) cycle
                        if (iabs(j2_ii-j2_ij)>2*J12) cycle
                        if (j2_ii+j2_ij<2*J12) cycle
                        if (ii==ij) then
                           if ((-1)**(J12+iT12)/=-1) cycle
                        endif
                        ind=ind+1
                        itbind(ii+ij*(ij-1)/2,J12,iT12)=ind
                     end do
                  end do   
               end do   
            end do   
         end do   
         if (ist2dim/=ind) then
            print *,' error: ind,ist2dim=',ind,ist2dim
            stop
         endif

         J_e2_2b=0   
         inquire(file=twobdop(1:ilast_tbd),exist=e2_2b)
         if (e2_2b) then
            allocate(e2_op_2b(ist2dim,ist2dim,-1:1))    
            e2_op_2b=0.d0
            open(93,file=twobdop(1:ilast_tbd),status='old',
     +           form='formatted',action='read')
            print *,twobdop(1:ilast_tbd),' file opened'

            read(93,*)
            read(93,*) J_e2_2b,T_e2_2b
            read(93,*)
            do
               read(93,*,end=76789) nonebody,nonebody2,j12,it12,
     +              nonebody3,nonebody4,j34,it34,rel

               if (nonebody>ist1dim) then
                  print *,' error: nonebody,ist1dim=',nonebody,ist1dim
                  stop
               endif
               if (nonebody2>ist1dim) then
                  print *,' error: nonebody,ist1dim=',nonebody2,ist1dim
                  stop
               endif
               if (nonebody3>ist1dim) then
                  print *,' error: nonebody,ist1dim=',nonebody3,ist1dim
                  stop
               endif
               if (nonebody4>ist1dim) then
                  print *,' error: nonebody,ist1dim=',nonebody4,ist1dim
                  stop
               endif
               if (j12/2>nhom12i+1.or.it12/2>1.or.j34/2>nhom12i+1
     +              .or.it34/2>1) then
                  print *,' error: j12,it12,j34,it34=',
     +                 j12/2,it12/2,j34/2,it34/2
                  stop
               endif

               if (nonebody<=nonebody2) then
                  i1=itbind(nonebody+nonebody2*(nonebody2-1)/2,
     +                 j12/2,it12/2)
                  pht12=1.d0
               else
                  i1=itbind(nonebody2+nonebody*(nonebody-1)/2,
     +                 j12/2,it12/2)
                  j1=ist1(3,nonebody)
                  j2=ist1(3,nonebody2)
                  pht12=real((-1)**((j12+it12-j1-j2)/2),kind(0.d0))
               endif
               if (nonebody3<=nonebody4) then
                  i3=itbind(nonebody3+nonebody4*(nonebody4-1)/2,
     +                 j34/2,it34/2)
                  pht34=1.d0
               else
                  i3=itbind(nonebody4+nonebody3*(nonebody3-1)/2,
     +                 j34/2,it34/2)
                  j1=ist1(3,nonebody3)
                  j2=ist1(3,nonebody4)
                  pht34=real((-1)**((j34+it34-j1-j2)/2),kind(0.d0))
               endif
               if (i1>ist2dim.or.i3>ist2dim) then
                  print *,' error: i1,i3,ist2dim=',i1,i3,ist2dim
                  stop
               endif

c               e2_op_2b(i1,i3,0)=rel*pht12*pht34
c               e2_op_2b(i1,i3,1)=rel*pht12*pht34
c               e2_op_2b(i1,i3,-1)=rel*pht12*pht34

c               e2_op_2b(i1,i3,0)=rel*pht12*pht34
c     +              /dsqrt(it12+1.d0)
c     +              *clebd(it34,0,2*t_e2_2b,0,it12,0)
c               e2_op_2b(i1,i3,1)=rel*pht12*pht34
c     +              /dsqrt(it12+1.d0)
c     +              *clebd(it34,2,2*t_e2_2b,0,it12,2)
c               e2_op_2b(i1,i3,-1)=rel*pht12*pht34
c     +              /dsqrt(it12+1.d0)
c     +              *clebd(it34,-2,2*t_e2_2b,0,it12,-2)
c               e2_op_2b(i3,i1,:)=e2_op_2b(i1,i3,:)
c     +              *(-1)**((j12-j34+it12-it34)/2)

               ip1=(-1)**(t_e2_2b-(it34-it12)/2)
               ip2=(-1)**(t_e2_2b+(j12-j34)/2)

               do iq=-1,1
               e2_op_2b(i1,i3,iq)=ip1*rel*pht12*pht34
     +              *clebd(2*t_e2_2b,0,it34,2*iq,it12,2*iq)
     &               /dsqrt(it12+1.d0)
               e2_op_2b(i3,i1,iq)=ip2*rel*pht12*pht34
     +              *clebd(2*t_e2_2b,0,it12,2*iq,it34,2*iq)
     &               /dsqrt(it34+1.d0)
               enddo

            end do
76789       continue
            close(93)
            print *,' e2_op_2b read from file'
         endif

         allocate(spin_spin(ist2dim,ist2dim,-1:1,-1:1))
         call spindotspin

         inquire(file=intbme(1:ilast),exist=obscal)
         if (obscal) then
            allocate(ham(ist2dim,ist2dim,-1:1))    
            allocate(tkin(ist2dim,ist2dim))    
            allocate(vcoul(ist2dim,ist2dim))
            ham=0.d0
            tkin=0.d0
            vcoul=0.d0
            open(66,file=intbme(1:ilast),status='old')
            read(66,*) isum
            if (isum/=isum2) write(iunitobout,*)
     +           ' *** Warning: isum2, isum:',isum2,isum     
            inquire(file='rrel2b.int',exist=rad2b)
            if (ipn/=1) rad2b=.false.
            if (rad2b) then
               allocate(rreleff(ist2dim,ist2dim,-1:1))    
c               allocate(rrel(ist2dim,ist2dim))    
               rreleff=0.d0
c               rrel=0.d0
               open(76,file='rrel2b.int',status='old')
               read(76,*) isum
               if (isum/=isum2) write(iunitobout,*)
     +              ' *** Warning: isum2, isum:',isum2,isum     
               read(76,*) numspR
               allocate(Rsp(numspR,-1:1))
               Rsp=0.d0
               do ii=1,numspR
                  read(76,*) i1,Rsp(ii,1),Rsp(ii,-1)
                  if (i1/=ii) print *,'***error: ii,i1=',ii,i1
               end do
               read(76,*) Rcore

               read(66,*) numsp
               if (numsp/=numspR) print *,' numsp,numspR=',numsp,numspR
               allocate(Esp(numsp,-1:1))
               Esp=0.d0
               do ii=1,numsp
                  read(66,*) i1,Esp(ii,1),Esp(ii,-1)
                  if (i1/=ii) print *,'***error: ii,i1=',ii,i1
               end do
               read(66,*) Ecore

            endif   
            if (ipn==1) then
               inquire(file='rreleffnp2b.int',exist=rad2beff)
               if (rad2beff) then
                  open(86,file='rreleffnp2b.int',status='old')
                  read(86,*) isum
                  if (isum/=isum2) write(iunitobout,*)
     +                 ' *** Warning: isum2, isum:',isum2,isum   
               endif   
               inquire(file='rreleffpp2b.int',exist=rad2beffx)
               rad2beff=rad2beff.and.rad2beffx
               if (rad2beff) then
                  open(87,file='rreleffpp2b.int',status='old')
                  read(87,*) isum
                  if (isum/=isum2) write(iunitobout,*)
     +                 ' *** Warning: isum2, isum:',isum2,isum
               endif
               inquire(file='rreleffnn2b.int',exist=rad2beffx)
               rad2beff=rad2beff.and.rad2beffx     
               if (rad2beff) then
                  open(88,file='rreleffnn2b.int',status='old')
                  read(88,*) isum
                  if (isum/=isum2) write(iunitobout,*)
     +                 ' *** Warning: isum2, isum:',isum2,isum     
               endif
            else   
               inquire(file='rreleff2b.int',exist=rad2beff)
               if (rad2beff) then
                  open(86,file='rreleff2b.int',status='old')
                  read(86,*) isum
                  if (isum/=isum2) write(iunitobout,*)
     +                 ' *** Warning: isum2, isum:',isum2,isum     
               endif   
            endif   
            do
c            do i1=1,ist2dim
c               nonebody=ist2(1,i1)
c               nonebody2=ist2(2,i1)
c               j12=ist2(3,i1)
c               it12=ist2(4,i1)
               
c               l1=ist1(2,nonebody)
c               l2=ist1(2,nonebody2)

c               do i3=i1,ist2dim
c                  nonebody3=ist2(1,i3)
c                  nonebody4=ist2(2,i3)
c                  j34=ist2(3,i3)
c                  it34=ist2(4,i3)
                  
c                  if (j12/=j34) cycle
c                  if (it12/=it34) cycle
                  
c                  l3=ist1(2,nonebody3)
c                  l4=ist1(2,nonebody4)

c                  if ((-1)**(l1+l2+l3+l4)/=1) cycle

c                  if (rad2b) then
c15878                continue
c                     read(76,*,end=5667)i1x,i2x,i3x,i4x,j34,it34,
c     +                    rel
c                     if (i1x/=nonebody.or.i2x/=nonebody2.or.
c     +                    i3x/=nonebody3.or.i4x/=nonebody4.or.
c     +                    j34/=j12.or.it34/=it12) goto 15878
c                  endif
                     
                  vpn=0.d0
                  vpp=0.d0
                  vnn=0.d0
                  if (ipn==1) then
15876                continue
                     read(66,*,end=5667)i1x,i2x,i3x,i4x,j34,it34,
     +                    trel,hrel,coul,vpn,vpp,vnn
                     coul=vpp-vnn

                     if (rad2b) then
                        read(76,*,end=5667) i1x_r,i2x_r,i3x_r,i4x_r,
     $                       j34_r,it34_r,releffnp,releffpp,releffnn
                        if (i1x/=i1x_r.or.i2x/=i2x_r.or.
     +                       i3x/=i3x_r.or.i4x/=i4x_r.or.
     +                       j34/=j34_r.or.it34/=it34_r) then
                           print *,'***incompatible radius file***'
                           releffnp=0.d0;releffpp=0.d0;releffnn=0.d0
                        endif
                     endif

c     if (i1x/=nonebody.or.i2x/=nonebody2.or.
c     +                    i3x/=nonebody3.or.i4x/=nonebody4.or.
c     +                    j34/=j12.or.it34/=it12) goto 15876
c                     if (rad2beff) then
c75879                   continue
c                        read(86,*,end=5667)i1x,i2x,i3x,i4x,j34,it34,
c     +                       releff
c                        if (i1x/=nonebody.or.i2x/=nonebody2.or.
c     +                       i3x/=nonebody3.or.i4x/=nonebody4.or.
c     +                       j34/=j12.or.it34/=it12) goto 75879
c85879                   continue
c                        read(87,*,end=5667)i1x,i2x,i3x,i4x,j34,it34,
c     +                       releffpp
c                        if (i1x/=nonebody.or.i2x/=nonebody2.or.
c     +                       i3x/=nonebody3.or.i4x/=nonebody4.or.
c     +                       j34/=j12.or.it34/=it12) goto 85879
c95879                   continue
c                        read(88,*,end=5667)i1x,i2x,i3x,i4x,j34,it34,
c     +                       releffnn
c                        if (i1x/=nonebody.or.i2x/=nonebody2.or.
c     +                       i3x/=nonebody3.or.i4x/=nonebody4.or.
c     +                       j34/=j12.or.it34/=it12) goto 95879
c                        rreleff(i1,i3,0)=releff*2.d0*bsquare
c                        rreleff(i1,i3,1)=releffpp*2.d0*bsquare
c                        rreleff(i1,i3,-1)=releffnn*2.d0*bsquare
c                        rreleff(i3,i1,:)=rreleff(i1,i3,:)
c                     endif   
                  else
15877                continue
                     read(66,*,end=5667)i1x,i2x,i3x,i4x,j34,it34,
     +                    trel,hrel,coul,vpn
c                     if (i1x/=nonebody.or.i2x/=nonebody2.or.
c     +                    i3x/=nonebody3.or.i4x/=nonebody4.or.
c     +                    j34/=j12.or.it34/=it12) goto 15877

c                     if (rad2beff) then
c15879                   continue
c                        read(86,*,end=5667)i1x,i2x,i3x,i4x,j34,it34,
c     +                       releff
c                        if (i1x/=nonebody.or.i2x/=nonebody2.or.
c     +                       i3x/=nonebody3.or.i4x/=nonebody4.or.
c     +                       j34/=j12.or.it34/=it12) goto 15879
c                        rreleff(i1,i3,:)=releff*2.d0*bsquare
c                        rreleff(i3,i1,:)=rreleff(i1,i3,:)
c                     endif   
                     if (it34==1) then
                        coul=coul*sqrt(938.9185*hboi)/197.327
                        vpp=vpn+coul
c*                        vpp=vpn
                        vnn=vpn
                     endif   
                  endif
                  if (i1x>ist1dim.or.i2x>ist1dim
     $                 .or.i3x>ist1dim.or.i4x>ist1dim) cycle
                  i1=itbind(i1x+i2x*(i2x-1)/2,j34,it34)
                  i3=itbind(i3x+i4x*(i4x-1)/2,j34,it34)
                  ham(i1,i3,0)=2.d0*trel*hboi/dble(nucleonsi)
     +                 + vpn 
                  ham(i3,i1,0)=ham(i1,i3,0)  
                  ham(i1,i3,1)=2.d0*trel*hboi/dble(nucleonsi)
     +                 + vpp 
                  ham(i3,i1,1)=ham(i1,i3,1)
                  ham(i1,i3,-1)=2.d0*trel*hboi/dble(nucleonsi)
     +                 + vnn 
                  ham(i3,i1,-1)=ham(i1,i3,-1)
                  tkin(i1,i3)=2.d0*trel*hboi/dble(nucleonsi)
                  tkin(i3,i1)=tkin(i1,i3)                     
                  vcoul(i1,i3)=coul
                  vcoul(i3,i1)=vcoul(i1,i3)                     
                  if (rad2b) then
                     ham(i1,i3,0)=vpn 
                     ham(i3,i1,0)=ham(i1,i3,0)  
                     ham(i1,i3,1)=vpp 
                     ham(i3,i1,1)=ham(i1,i3,1)
                     ham(i1,i3,-1)=vnn 
                     ham(i3,i1,-1)=ham(i1,i3,-1)
                     rreleff(i1,i3,0)=releffnp !*2.d0*bsquare
                     rreleff(i1,i3,1)=releffpp !*2.d0*bsquare
                     rreleff(i1,i3,-1)=releffnn !*2.d0*bsquare
                     rreleff(i3,i1,:)=rreleff(i1,i3,:)
c                     rrel(i1,i3)=rel*2.d0*bsquare
c                     rrel(i3,i1)=rrel(i1,i3)
                  endif
c               end do
            end do
 5667       continue
            close(66)
            if (rad2b) close(76)
       
         else
            write(iunitobout,*) ' *** Interaction file ',
     +           intbme(1:ilast),' does not exist'
         endif   
      endif   

      allocate(sumN(0:nhomi,2))
      allocate(occup(ist1dim,2))

      do kii=ki,ki+nki-1

         jtotali2=jt2i(kii)
         itotali2=it2i(kii)
         if (dabs(dble(jtotali2)-2.d0*Jxi(kii))>1.d-4) then
            print *, '**** not a good J ****'
c            cycle
         endif 

         do kff=kf,kf+nkf-1
            
            jtotalf2=jt2f(kff)
            itotalf2=it2f(kff)
            if (dabs(dble(jtotalf2)-2.d0*Jxf(kff))>4.d-4) then
               write(iunitobout,*) '**** not a good J ****'
c               cycle
            endif 

            write(iunitobout,1201) kff,jtotalf2,itotalf2,
     +           enerf(kff)-eneri(ki),
     +           kii,jtotali2,itotali2,eneri(kii)-eneri(ki)
 1201       format(/,' *** Transition matrix elements for states:',
     +           ' ***',/,' #',i3,
     +           ' [2*(J,T),Ex]_f=',i3,i2,f8.4,3x,
     +           '#',i3,' [2*(J,T),Ex]_i=',i3,i2,f8.4) 

            do jtrans=jtotal2min,jtotal2

               jtrans2=jtrans+jtrans
               if (iabs(jt2i(kii)-jt2f(kff))>jtrans2
     +              .or.(jt2i(kii)+jt2f(kff))<jtrans2) cycle
               
               cleb=clebd(jtotali2,mjtotali,jtrans2,
     +              mjtotalf-mjtotali,jtotalf2,mjtotalf)
               if (dabs(cleb)<1.d-8) then
                  cycle
               else 
                  cleb=dsqrt(dble(jtotalf2)+1.d0)/cleb 
                  clb1=-cleb/dsqrt(dble(jtrans2)+1.d0) 
               endif

c               if (twobdcal.and.obscal.and.snsp.and.
               if (twobdcal.and.
     +              (jtrans==0.or.jtrans==J_e2_2b)) then
c*     +              .and.kii==kff) then
c***** Non-energy weighted sum rules ***
                  sumrule=0.d0
c***************************************
                  e2_op_sum=0.d0
                  xxx=0.d0
                  sumto=0.d0
                  sumco=0.d0
                  rad0=0.d0
                  radeff0=0.d0
                  radpp=0.d0
                  radeffpp=0.d0
                  SpSp=0.d0
                  SpSn=0.d0
                  SnSn=0.d0
                  do i1=1,ist2dim

                     ob1=ist2(1,i1)
                     ob2=ist2(2,i1)
                     j12=ist2(3,i1)
                     it12=ist2(4,i1)
                     
                     l1=ist1(2,ob1)
                     l2=ist1(2,ob2)

                     do i3=1,ist2dim

                        ob3=ist2(1,i3)
                        ob4=ist2(2,i3)
                        j34=ist2(3,i3)
                        it34=ist2(4,i3)

                        l3=ist1(2,ob3)
                        l4=ist1(2,ob4)
                  
                        if (mod(l1+l2+l3+l4+iparityi+iparityf,2)==1)
     $                       cycle

                        if (jtrans==0.and.j12==j34) then
                           if (it12==0.or.it34==0) then
                              SpSp=SpSp
     $                             +sqrt(2.d0*j12+1.d0)
     $                             *spin_spin(i1,i3,0,1)
     $                             *tbd(kff,kii)%Jtr(jtrans)%fi(i1)
     $                             %mT(0)%in(i3)
                              SpSn=SpSn+sqrt(2.d0*j12+1.d0)
     $                             *spin_spin(i1,i3,0,0)
     $                             *tbd(kff,kii)%Jtr(jtrans)%fi(i1)
     $                             %mT(0)%in(i3)
                              SnSn=SnSn
     $                             +sqrt(2.d0*j12+1.d0)
     $                             *spin_spin(i1,i3,0,-1)
     $                             *tbd(kff,kii)%Jtr(jtrans)%fi(i1)
     $                             %mT(0)%in(i3)
                           else
                              SpSp=SpSp+sqrt(2.d0*j12+1.d0)
     $                             *spin_spin(i1,i3,1,1)
     $                             *tbd(kff,kii)%Jtr(jtrans)%fi(i1)
     $                             %mT(1)%in(i3)
     $                             +sqrt(2.d0*j12+1.d0)
     $                             *spin_spin(i1,i3,0,1)
     $                             *tbd(kff,kii)%Jtr(jtrans)%fi(i1)
     $                             %mT(0)%in(i3)
                              SpSn=SpSn+sqrt(2.d0*j12+1.d0)
     $                             *spin_spin(i1,i3,0,0)
     $                             *tbd(kff,kii)%Jtr(jtrans)%fi(i1)
     $                             %mT(0)%in(i3)
                              SnSn=SnSn+sqrt(2.d0*j12+1.d0)
     $                             *spin_spin(i1,i3,-1,-1)
     $                             *tbd(kff,kii)%Jtr(jtrans)%fi(i1)
     $                             %mT(-1)%in(i3)
     $                             +sqrt(2.d0*j12+1.d0)
     $                             *spin_spin(i1,i3,0,-1)
     $                             *tbd(kff,kii)%Jtr(jtrans)%fi(i1)
     $                             %mT(0)%in(i3)
                           endif
                        endif

c                        print *,' i1,i3,xxx=',i1,i3,xxx

                        if (e2_2b.and.jtrans==J_e2_2b) then
                           if (it12==1) then
                              e2_op_sum=e2_op_sum+e2_op_2b(i1,i3,1)
c     +                          *t2bd(i1,i3,jtrans,1,kff,kii)
     +                             *tbd(kff,kii)%Jtr(jtrans)%fi(i1)
     $                             %mT(1)%in(i3)
     +                          +e2_op_2b(i1,i3,0)
     +                             *tbd(kff,kii)%Jtr(jtrans)%fi(i1)
     $                             %mT(0)%in(i3)
c     +                          *t2bd(i1,i3,jtrans,0,kff,kii)
     +                          +e2_op_2b(i1,i3,-1)
     +                             *tbd(kff,kii)%Jtr(jtrans)%fi(i1)
     $                             %mT(-1)%in(i3)
c     +                          *t2bd(i1,i3,jtrans,-1,kff,kii)
                           else
                              e2_op_sum=e2_op_sum+e2_op_2b(i1,i3,0)
     +                             *tbd(kff,kii)%Jtr(jtrans)%fi(i1)
     $                             %mT(0)%in(i3)
                           endif
                        endif

                        if (.not.snsp) cycle
                        if (jtrans/=0) cycle
                        if (j12/=j34) cycle
                        if (it12/=it34) cycle

                        if (it12==0) then
                           xxx=xxx
     +                          +ham(i1,i3,0)*dsqrt(dble(2*j12+1))
     +                          *tbd(kff,kii)%Jtr(0)%fi(i1)
     $                          %mT(0)%in(i3)
c     +                          *t2bd(i1,i3,0,0,kff,kii)
                           sumto=sumto+tkin(i1,i3)*dsqrt(dble(2*j12+1))
     +                          *tbd(kff,kii)%Jtr(0)%fi(i1)
     $                          %mT(0)%in(i3)
c     +                          *t2bd(i1,i3,0,0,kff,kii)
                           if (rad2beff) then
                              radeffpp=radeffpp
     +                          +rreleff(i1,i3,0)*dsqrt(dble(2*j12+1))
     +                          *tbd(kff,kii)%Jtr(0)%fi(i1)
     $                          %mT(0)%in(i3)
c     +                          *t2bd(i1,i3,0,0,kff,kii)
                              radeff0=radeff0
     +                          +rreleff(i1,i3,0)*dsqrt(dble(2*j12+1))
     +                          *tbd(kff,kii)%Jtr(0)%fi(i1)
     $                          %mT(0)%in(i3)
c     +                          *t2bd(i1,i3,0,0,kff,kii)
                           endif
                           if (rad2b) then
                              rad0=rad0+rreleff(i1,i3,0) !rrel(i1,i3)
     +                             *dsqrt(dble(2*j12+1))
     +                             *tbd(kff,kii)%Jtr(0)%fi(i1)
     $                             %mT(0)%in(i3)
c                              radpp=radpp+rreleff(i1,i3,0) !rrel(i1,i3)
c     +                             *dsqrt(dble(2*j12+1))
c     +                             *tbd(kff,kii)%Jtr(0)%fi(i1)
c     $                             %mT(0)%in(i3)
                           endif
c***** Non-energy weighted sum rules ***
                           do k=jtrmin,jtrmax
                              do ids=0,3
                                 sumrule(k,ids)=sumrule(k,ids)+
     +                                twobsum(ob1,ob2,
     +                    j12,it12,0,k,ids,ob3,ob4,j34,it34,0)
     +                          *tbd(kff,kii)%Jtr(0)%fi(i1)
     $                          %mT(0)%in(i3)
c     +                             *t2bd(i1,i3,0,0,kff,kii)     
                              end do
                           end do
c***************************************
                        elseif(it12==1) then
                           xxx=xxx
     +                          +dsqrt(dble(2*j12+1))*(ham(i1,i3,0)*
     +                          tbd(kff,kii)%Jtr(0)%fi(i1)
     $                          %mT(0)%in(i3)
c     +                          t2bd(i1,i3,0,0,kff,kii)
     +                          +ham(i1,i3,1)*  !t2bd(i1,i3,0,1,kff,kii)
     +                          tbd(kff,kii)%Jtr(0)%fi(i1)
     $                          %mT(1)%in(i3)
     +                         +ham(i1,i3,-1)*  !t2bd(i1,i3,0,-1,kff,kii)
     $                          tbd(kff,kii)%Jtr(0)%fi(i1)
     $                          %mT(-1)%in(i3))
                           sumto=sumto
     +                          +dsqrt(dble(2*j12+1))*tkin(i1,i3)*
     $                          (tbd(kff,kii)%Jtr(0)%fi(i1)%mT(0)%in(i3)
     $                          +tbd(kff,kii)%Jtr(0)%fi(i1)%mT(1)%in(i3)
     $                          +tbd(kff,kii)%Jtr(0)%fi(i1)%mT(-1)
     $                          %in(i3))
c     +                          (t2bd(i1,i3,0,0,kff,kii)
c     +                          +t2bd(i1,i3,0,1,kff,kii)
c     +                          +t2bd(i1,i3,0,-1,kff,kii))
                           sumco=sumco
     +                          +dsqrt(dble(2*j12+1))*vcoul(i1,i3)*
     $                          tbd(kff,kii)%Jtr(0)%fi(i1)%mT(1)%in(i3)
c     +                          t2bd(i1,i3,0,1,kff,kii)
                           if (rad2beff) then
                              radeffpp=radeffpp
     +                         +dsqrt(dble(2*j12+1))*(rreleff(i1,i3,0)*
     $                          tbd(kff,kii)%Jtr(0)%fi(i1)%mT(0)%in(i3)
c     +                          t2bd(i1,i3,0,0,kff,kii)
     +                   +2.d0*rreleff(i1,i3,1)*  !t2bd(i1,i3,0,1,kff,kii)
     $                         tbd(kff,kii)%Jtr(0)%fi(i1)%mT(1)%in(i3))
                              radeff0=radeff0
     +                          +dsqrt(dble(2*j12+1))*(rreleff(i1,i3,0)*
     $                           tbd(kff,kii)%Jtr(0)%fi(i1)%mT(0)%in(i3)
c     +                          t2bd(i1,i3,0,0,kff,kii)
     +                             +rreleff(i1,i3,1)* !t2bd(i1,i3,0,1,kff,kii)
     $                           tbd(kff,kii)%Jtr(0)%fi(i1)%mT(1)%in(i3)
     +                             +rreleff(i1,i3,-1)* !t2bd(i1,i3,0,-1,kff,kii)
     $                         tbd(kff,kii)%Jtr(0)%fi(i1)%mT(-1)%in(i3))
                           endif
                           if (rad2b) then
                              rad0=rad0
     +                          +dsqrt(dble(2*j12+1))*(rreleff(i1,i3,0)*
     +                          tbd(kff,kii)%Jtr(0)%fi(i1)
     $                          %mT(0)%in(i3)
     +                          +rreleff(i1,i3,1)*  
     +                          tbd(kff,kii)%Jtr(0)%fi(i1)
     $                          %mT(1)%in(i3)
     +                         +rreleff(i1,i3,-1)*  
     $                          tbd(kff,kii)%Jtr(0)%fi(i1)
     $                          %mT(-1)%in(i3))
                              
c                              rad0=rad0
c     +                          +dsqrt(dble(2*j12+1))*rrel(i1,i3)*
c     $                          (tbd(kff,kii)%Jtr(0)%fi(i1)%mT(0)%in(i3)
c     $                          +tbd(kff,kii)%Jtr(0)%fi(i1)%mT(1)%in(i3)
c     $                          +tbd(kff,kii)%Jtr(0)%fi(i1)%mT(-1)
c     $                          %in(i3))
c                              radpp=radpp
c     +                          +dsqrt(dble(2*j12+1))*rrel(i1,i3)*(
c     $                          tbd(kff,kii)%Jtr(0)%fi(i1)%mT(0)%in(i3)
c     $                         +2.d0*tbd(kff,kii)%Jtr(0)%fi(i1)
c     $                             %mT(1)%in(i3))
                           endif
c***** Non-energy weighted sum rules ***
                           do k=jtrmin,jtrmax
                              do ids=0,3
                                 sumrule(k,ids)=sumrule(k,ids)+
     +                                twobsum(ob1,ob2,
     +                    j12,it12,0,k,ids,ob3,ob4,j34,it34,0)
     $                    *tbd(kff,kii)%Jtr(0)%fi(i1)%mT(0)%in(i3)
c     +                                *t2bd(i1,i3,0,0,kff,kii)     
     +                                +twobsum(ob1,ob2,
     +                    j12,it12,1,k,ids,ob3,ob4,j34,it34,1)
     $                    *tbd(kff,kii)%Jtr(0)%fi(i1)%mT(1)%in(i3)
c     +                             *t2bd(i1,i3,0,1,kff,kii)     
     +                                +twobsum(ob1,ob2,
     +                  j12,it12,-1,k,ids,ob3,ob4,j34,it34,-1)
     $                    *tbd(kff,kii)%Jtr(0)%fi(i1)%mT(-1)%in(i3)
c     +                             *t2bd(i1,i3,0,-1,kff,kii)     
                              end do
                           end do
c***************************************
                        else
                           write(2,*)'**** error: it12=',it12
                        endif
*                        print *,'i1,i3,sumto:',i1,i3,sumto
*                        print *,' tkin=',tkin(i1,i3)
                     end do
                  end do

                  if (rad2b.and.jtrans==0) then
                     Eobd=0.d0
                     Robd=0.d0
                     do i1=1,ist1dim
                        n1=ist1(1,i1)
                        l1=ist1(2,i1)
                        j1=ist1(3,i1)
                        dsqj=dsqrt(j1+1.d0)
                        Eobd=Eobd+Esp(i1,1)*dsqj*tdJp(i1,i1,0,kff,kii)
     $                       +Esp(i1,-1)*dsqj*tdJn(i1,i1,0,kff,kii)
                        Robd=Robd+Rsp(i1,1)*dsqj*tdJp(i1,i1,0,kff,kii)
     $                       +Rsp(i1,-1)*dsqj*tdJn(i1,i1,0,kff,kii)
                     end do
                     print *,' Ecore,Eobd,V2b=',
     $                    Ecore*krd(kff,kii),
     $                    Eobd/dsqrt(dble(jtotalf2+1)),
     $                    xxx/dsqrt(dble(jtotalf2+1))
                     print *,' Rcore,Robd,R2b=',
     $                    Rcore*krd(kff,kii),
     $                    Robd/dsqrt(dble(jtotalf2+1)),
     $                    rad0/dsqrt(dble(jtotalf2+1))
                  else
                     if (.not.rad2b) then
                        Eobd=0.d0
                        Robd=0.d0
                        Ecore=0.d0
                        Rcore=0.d0
                     endif
                  endif
                  
                  if (jtrans==0) write(iunitobout,6656)
     $                 xxx/dsqrt(dble(jtotalf2+1))
     $                 +Eobd/dsqrt(dble(jtotalf2+1))+Ecore*krd(kff,kii),
     +                 sumto/dsqrt(dble(jtotalf2+1)),
     +                 sumco/dsqrt(dble(jtotalf2+1))
 6656             format(/,
     + ' Hamiltonian, kinetic and Coulomb operator mean values',/,
     +                 ' < H > = ',f10.5,
     +                 '     < T > =',f9.4,'     < V_Coul > =',f8.4) 

                  if (jtrans==0) then
                     SpSp=SpSp/dsqrt(dble(jtotalf2+1))
                     SpSn=SpSn/dsqrt(dble(jtotalf2+1))
                     SnSn=SnSn/dsqrt(dble(jtotalf2+1))
                     write(iunitobout,
     $       "(/,' <SpSp>=',f8.4,'   <SnSn>=',f8.4,'   <SpSn>=',f8.4)")
     $                    SpSp,SnSn,SpSn
                     write(iunitobout,"(' <S^2>=',f8.4)")
     $                    SpSp+SnSn+2.d0*SpSn
                  endif

                  if (e2_2b.and.jtrans==0.and.J_e2_2b==0) then
                     write(iunitobout,"(/,' <O_2b(T)>=',f8.4,
     +                    ' <O_2b(Trel)>=',f8.4)")
c     +                    e2_op_sum/dble(nucleonsi-1)
     +                    e2_op_sum/dble(nucleonsi)
     +                    /dsqrt(dble(jtotalf2+1)),
c     +                    e2_op_sum/dble(nucleonsi-1)
c     +                    /dsqrt(dble(jtotalf2+1))-0.75d0*hboi
     +                    e2_op_sum/dble(nucleonsi)
     +                    /dsqrt(dble(jtotalf2+1))
                  endif
                  if (e2_2b.and.jtrans==2.and.J_e2_2b==2) then
                     write(iunitobout,"(/,' <O_2b(E2)>=',f8.4,
     +                    ' B(E2)=',f8.4)")
     +                    e2_op_sum
     +                    /dble(nucleonsi-1),
     +                    e2_op_sum**2/dble(jtotali2+1)
     +                    /dble(nucleonsi-1)**2
                  endif
                  if (rad2b.and.jtrans==0) then
c                     if (nprotonsi>0) then
c                        radprot=(radpp/(dble(nucleonsi*nprotonsi))
c     +               -rad0/dble(nucleonsi**2))/dsqrt(dble(jtotalf2+1))
c                     else
c                        radprot=0.d0
c                     endif   
c                     radmass=rad0/dble(nucleonsi**2)
c     +                    /dsqrt(dble(jtotalf2+1))
c                     if (nneutrnsi>0) then
c                        radneut=(dble(nucleonsi)*radmass
c     +                 -dble(nprotonsi)*radprot)/dble(nneutrnsi)
c                     else
c                        radneut=0.d0
c                     endif   
c                     write(iunitobout,6956) dsqrt(dabs(radprot)),
c     +                 dsqrt(dabs(radneut)),
c     +                 dsqrt(dabs(radmass))
c 6956                format(/,
c     +               ' Point rms radius operator mean values',/,
c     +                 ' < r_p > = ',f10.5,
c     +             '      < r_n > =',f10.5,'      < r_m > =',f10.5) 


                     write(iunitobout,"(/,' Rcore,Robd,R2b=',3f10.5)")
     $                    Rcore*krd(kff,kii),
     $                    Robd/dsqrt(dble(jtotalf2+1)),
     $                    rad0/dsqrt(dble(jtotalf2+1))
                     
                     rad0=rad0/dsqrt(dble(jtotalf2+1))
     $                    +Robd/dsqrt(dble(jtotalf2+1))
     $                    +Rcore*krd(kff,kii)
                     
                     write(iunitobout,6956) rad0,dsqrt(dabs(rad0))
 6956                format(/,
     +               ' Point rms radius operator mean values',/,
     +                 ' < r^2 > = ',f10.5,
     +             '      < r > =',f10.5)
                  endif   
                  if (rad2beff.and.jtrans==0) then
                     if (nprotonsi>0) then
                        radprot=(radeffpp/(dble(nucleonsi*nprotonsi))
     +             -radeff0/dble(nucleonsi**2))/dsqrt(dble(jtotalf2+1))
                     else
                        radprot=0.d0
                     endif   
                     radmass=radeff0/dble(nucleonsi**2)
     +                 /dsqrt(dble(jtotalf2+1))
                     if (nneutrnsi>0) then
                        radneut=(dble(nucleonsi)*radmass
     +                 -dble(nprotonsi)*radprot)/dble(nneutrnsi)
                     else
                        radneut=0.d0
                     endif   
                     write(iunitobout,6976) dsqrt(dabs(radprot)),
     +                    dsqrt(dabs(radneut)),
     +                    dsqrt(dabs(radmass))
 6976                format(/,
     +            ' Point rms radius effective operator mean values',/,
     +                 ' < r_p > = ',f10.5,
     +             '      < r_n > =',f10.5,'      < r_m > =',f10.5) 
                  endif
c***** Non-energy weighted sum rules ***
                  if (jtrans==0) then
                     write(iunitobout,
     +                    "(/,' Non-energy weighted sum rules')")
                     sumrule=sumrule/dsqrt(dble(jtotalf2+1))
                     do k=jtrmin,jtrmax
                        write(iunitobout,"(' E',i2,'=',f8.4)") 
     +                       k,sumrule(k,0)
                        select case(k)
                        case(0)
                           write(iunitobout,"(' (Y1 sigma)^',i1,
     +                          '=',f8.4)") 
     +                          k,sumrule(k,2)
                        case(1)
                           write(iunitobout,"(' M1=',f8.4)")
     +                          sumrule(k,1)
                           write(iunitobout,"(' (Y1 sigma)^',i1,
     +                          '=',f8.4)") 
     +                          k,sumrule(k,2)
                           write(iunitobout,"(' GT=',f8.4)")
     +                          sumrule(k,3)
                        case(2)
                           write(iunitobout,"(' (Y1 sigma)^',i1,
     +                          '=',f8.4)") 
     +                          k,sumrule(k,2)
                        case default
                           continue
                        end select
                     end do
                  endif
c***************************************
               endif   
c*** end of twobdcal

               sumN=0.d0
               occup=0.d0
               sumpl=0.d0
               sumnl=0.d0
               sumps=0.d0
               sumns=0.d0
               sumgt=0.d0
               ELp=0.d0
               ELn=0.d0
               E1sp=0.d0
               E1sn=0.d0
               E1s=0.d0

               do i1=1,ist1dim
                  n1=ist1(1,i1)
                  l1=ist1(2,i1)
                  j1=ist1(3,i1)
                  do i3=1,ist1dim
                     n2=ist1(1,i3)
                     l2=ist1(2,i3)
                     j2=ist1(3,i3)
                     if((j1+j2)<2*jtrans) cycle
                     if(iabs(j1-j2)>2*jtrans) cycle
                     if (2*n1+l1==2*n2+l2.and.j1==j2
     +                    .and.jtrans==0) then
                        dsqj=dsqrt(j1+1.d0)
           sumN(2*n1+l1,1)=sumN(2*n1+l1,1)+dsqj*tdJp(i1,i3,0,kff,kii)
           sumN(2*n1+l1,2)=sumN(2*n1+l1,2)+dsqj*tdJn(i1,i3,0,kff,kii)
           occup(i1,1)=occup(i1,1)+dsqj*tdJp(i1,i3,0,kff,kii)
           occup(i1,2)=occup(i1,2)+dsqj*tdJn(i1,i3,0,kff,kii)
                     endif
 
                     if (snsp.and.jtrans==1.and.2*n1+l1==2*n2+l2) then
                        sumpl=sumpl+obme(n1,l1,j1,n2,l2,j2,1,2)
     +                         *tdJp(i1,i3,jtrans,kff,kii)
                        sumnl=sumnl+obme(n1,l1,j1,n2,l2,j2,1,2)
     +                         *tdJn(i1,i3,jtrans,kff,kii)
                        sumps=sumps+obme(n1,l1,j1,n2,l2,j2,1,3)
     +                         *tdJp(i1,i3,jtrans,kff,kii)
                        sumns=sumns+obme(n1,l1,j1,n2,l2,j2,1,3)
     +                         *tdJn(i1,i3,jtrans,kff,kii)
                     endif   

                     if (dnuc.and.iparityi==iparityf.and.jtrans==1) then
                        if (n1==n2.and.l1==l2) then
                           sumgt=sumgt+dsqrt((j1+1.d0)
     +                         *(j2+1.d0)*6.d0)*racad(j2,2,2*l2,1,j1,1)
     +                           *(tdJp(i3,i1,jtrans,kff,kii)
     +                           +tdJn(i3,i1,jtrans,kff,kii))
                        endif   
                     endif               
            
                     ELp=ELp+
     +                    obme(n2,l2,j2,n1,l1,j1,jtrans,0)*
     +                    tdJp(i3,i1,jtrans,kff,kii)

                     ELn=ELn+
     +                    obme(n2,l2,j2,n1,l1,j1,jtrans,0)*
     +                    tdJn(i3,i1,jtrans,kff,kii)

                     if (sndp.and.jtrans<3) then
                        E1sp=E1sp+
     +                       obme(n2,l2,j2,n1,l1,j1,jtrans,4)*
     +                       tdJp(i3,i1,jtrans,kff,kii)
                        
                        E1sn=E1sn+
     +                       obme(n2,l2,j2,n1,l1,j1,jtrans,4)*
     +                       tdJn(i3,i1,jtrans,kff,kii)
                     endif

                     if (dnuc.and.jtrans<3.and.iparityi/=iparityf) then
                        E1s=E1s+
     +                       obme(n2,l2,j2,n1,l1,j1,jtrans,4)*
     +                       (tdJp(i3,i1,jtrans,kff,kii)
     +                       +tdJn(i3,i1,jtrans,kff,kii))
                     endif

                  end do
               end do

               if (snsp.and.kii==kff.and.jtrans==0) then
                  sumN=-sumN/clb1
                  occup=-occup/clb1
                  write(iunitobout,1970)
 1970             format(/,' Occupation Numbers for N=2n+l:')
                  do N1=0,nhomi
                     if (sumN(N1,1)/=0.d0.or.sumN(N1,2)/=0.d0) then
                        write(iunitobout,1980) N1,sumN(N1,1),sumN(N1,2)
                     endif
 1980                format(' N=',i3,'    p=',f8.4,'    n=',f8.4)
                  end do 
                  write(iunitobout,9970)
 9970             format(/,' Occupation Numbers for nlj-levels:')
                  do i1=1,ist1dim
                     n1=ist1(1,i1)
                     l1=ist1(2,i1)
                     j1=ist1(3,i1)
                     write(iunitobout,7980) n1,l1,j1,
     +                    occup(i1,1),occup(i1,2)
 7980                format(' n,l,2j=',3i3,'    p=',f8.4,'    n=',f8.4)
                  end do 
                  if (dabs(sum(sumN(:,1))-nprotonsi)>1.d-5) then
                     write(iunitobout,*)'*** error: sumN,protons=',
     +                    sum(sumN(:,1)),nprotonsi
                  endif   
                  if (dabs(sum(sumN(:,2))-nneutrnsi)>1.d-5) then
                     write(iunitobout,*)'*** error: sumN,neutrons=',
     +                    sum(sumN(:,2)),nneutrnsi
                  endif   
                  if (dabs(sum(occup(:,1))-nprotonsi)>1.d-5) then
                     write(iunitobout,*)'*** error: occup,protons=',
     +                    sum(occup(:,1)),nprotonsi
                  endif   
                  if (dabs(sum(occup(:,2))-nneutrnsi)>1.d-5) then
                     write(iunitobout,*)'*** error: occup,neutrons=',
     +                    sum(occup(:,2)),nneutrnsi
                  endif   
               endif

               if (dnuc.and.jtrans==0.and.iparityi==iparityf) then
                  sumf=sum(sumN)
                  sumf=sumf/dsqrt(dble(jtotalf2)+1.d0)
                  write(iunitobout,6754) sumf,sumf**2
 6754           format(' Fermi matrix element=',f15.8,'   B(F)=',f10.5)
               endif  

               if (dnuc.and.iparityi==iparityf.and.jtrans==1) then
                  sumgt=sumgt/dsqrt(dble(jtotali2)+1.d0)
                  write(iunitobout,6755) sumgt,sumgt**2
 6755             format(' Gamow-Teller matrix element=',f15.8,
     +                 '    B(GT)=',f10.5)
               endif    

               if (snsp.and.jtrans==1) then
                  sumto=sumpl+sumps+sumnl+sumns
                  jtx2=jt2i(kii) 
c*                  if (kii==kff) then
c*                     dsqj=dsqrt(jtx2/2.d0*(jtx2/2.d0+1.d0)
c*     +                    *(jtx2+1.d0))
c*                  else
c*                     dsqj=0.d0
c*                  endif
                  
c*                  write(iunitobout,1971)
c* 1971             format(/,' M1 matrix elements:')
c*                  write(iunitobout,1979) sumpl,sumnl,sumps,sumns,
c*     +                 sumto,dsqj
c* 1979             format(' pl=',f8.4,'   nl=',f8.4,'   ps=',f8.4,
c*     +                 '   ns=',f8.4,/,' sum=',f8.4,'    <J>=',f8.4)   

                  sumto=glp*sumpl+gln*sumnl+gsp*sumps+gsn*sumns
                  
                  if (kii==kff) then
                     clb=dsqrt(0.5d0*jtx2/(jtx2+1.d0)
     +                    /(0.5d0*jtx2+1.d0))
                     write(iunitobout,1975)
 1975                format(/,' M1 moments:')
                     write(iunitobout,1976) sumpl*clb,sumnl*clb,
     +                    sumps*clb,sumns*clb,sumto*clb
 1976                format(' pl=',f8.4,'   nl=',f8.4,'   ps=',f8.4,
     +                    '   ns=',f8.4,'     M1=',f8.4)   
                  endif   

                  clb=dsqrt(3.d0/4.d0/piln/
     +                 (dble(jtx2)+1.d0))
                  sumpl=sumpl*clb
                  sumnl=sumnl*clb
                  sumps=sumps*clb
                  sumns=sumns*clb  
                  sumto=sumto*clb
                  sumto=sumto*sumto
                  write(iunitobout,1973)
 1973             format(/,' BM1 matrix elements:')
                  write(iunitobout,1974) sumpl,sumnl,sumps,sumns,
     +                 sumto
 1974             format(' pl=',f8.4,'   nl=',f8.4,'   ps=',f8.4,
     +                    '   ns=',f8.4,' B(M1)=',f8.4)   
                     
               endif

               if (.not.dnuc.and.
     +              ((-1)**(iparityi+iparityf)==(-1)**jtrans)) then
                  jtx2=jt2i(kii) 
                  sumto=epi*ELp+enu*ELn

                  if (kii==kff.and.jtrans==2) then
                     clb=dsqrt(4.d0*jtx2*(jtx2-1.d0)/(jtx2+2.d0)
     +                    /(jtx2+1.d0)/(jtx2+3.d0))
                     write(iunitobout,1977)
 1977                format(/,' E2 moments:')
                     write(iunitobout,1978) ELp*clb,ELn*clb,sumto*clb
 1978                format(' Q2p=',f8.4,'   Q2n=',f8.4,
     +                                  '     Q2=',f8.4)   
                  endif   

                  clb=dsqrt((2*jtrans+1.d0)/4.d0/piln/
     +                 (dble(jtx2)+1.d0))
                  ELp=ELp*clb
                  ELn=ELn*clb
                  sumto=sumto*clb
                  sumto=sumto*sumto
                  write(iunitobout,1998) jtrans
 1998             format(/,' BE',i1,' matrix elements:')
                  write(iunitobout,1999) jtrans,jtrans,ELp,jtrans,ELn,
     +                 jtrans,sumto
 1999             format(' L=',i2,' E',i1,'p=',f8.4,
     +                          '   E',i1,'n=',f8.4,
     +                          '     B(E',i1,')=',f8.4)
               endif

               if (dnuc.and.
     +              ((-1)**(iparityi+iparityf)==(-1)**jtrans)) then
                  jtx2=jt2i(kii) 
                  sumto=ELp+ELn

                  clb=dsqrt((2*jtrans+1.d0)/4.d0/piln/
     +                 (dble(jtx2)+1.d0))
                  ELp=ELp*clb
                  ELn=ELn*clb
                  sumto=sumto*clb
                  sumto=sumto*sumto
                  write(iunitobout,1998) jtrans

                  write(iunitobout,1999) jtrans,jtrans,ELp,jtrans,ELn,
     +                 jtrans,sumto

               endif

               if (sndp.and.jtrans<3) then
                  E1sp=E1sp/dsqrt(4.d0*piln*
     +                 dble(jt2i(kii)+1))
                  E1sn=E1sn/dsqrt(4.d0*piln*
     +                 dble(jt2i(kii)+1))
                  sumto=epi*E1sp+enu*E1sn
                  sumto=sumto*sumto
                  write(iunitobout,2998) jtrans
 2998             format(/,' BE1s',i1,' matrix elements:')
                  write(iunitobout,2999) jtrans,jtrans,E1sp,
     +                 jtrans,E1sn,
     +                 jtrans,sumto
 2999             format(' k=',i2,' E1sp',i1,'=',f8.4,
     +                 '   E1sn',i1,'=',f8.4,
     +                 '     B(E1s',i1,')=',f8.4)
               endif
                  
               if (dnuc.and.jtrans<3.and.iparityi/=iparityf) then
                  E1s=E1s/dsqrt(4.d0*piln*dble(jt2i(kii)+1))
                  sumto=E1s*E1s
                  write(iunitobout,3998) jtrans
 3998             format(/,' BE1s',i1,' matrix elements:')
                  write(iunitobout,3999) jtrans,jtrans,E1s,
     +                 jtrans,sumto
 3999             format(' k=',i2,' E1s',i1,'=',f8.4,
     +                 '     B(E1s',i1,')=',f8.4)
               endif

            end do

            write(iunitobout,*)
     +'***************************************************************'
         end do
      end do

      contains
      subroutine spindotspin
      implicit none
      integer :: i1,i2,a,b,c,d,J,Jin,Tab,Tcd,MT
      integer :: na,la,ja,nb,lb,jb,nc,lc,jc,nd,ld,jd
      real(kind(0.d0)) :: term1,term2,clebd,racad
      spin_spin=0.d0
      do i1=1,ist2dim
         a=ist2(1,i1)
         b=ist2(2,i1)
         J=ist2(3,i1)
         Tab=ist2(4,i1)
         na=ist1(1,a)
         la=ist1(2,a)
         ja=ist1(3,a)
         nb=ist1(1,b)
         lb=ist1(2,b)
         jb=ist1(3,b)
         do i2=1,ist2dim
            c=ist2(1,i2)
            d=ist2(2,i2)
            Jin=ist2(3,i2)
            Tcd=ist2(4,i2)
            if (J/=Jin) cycle
            nc=ist1(1,c)
            lc=ist1(2,c)
            jc=ist1(3,c)
            nd=ist1(1,d)
            ld=ist1(2,d)
            jd=ist1(3,d)
            
            if (na==nc.and.la==lc.and.nb==nd.and.lb==ld) then
               term1=(1.d0+(-1)**(Tab+Tcd))
     $              *sqrt((ja+1.d0)*(jb+1.d0)*(jc+1.d0)*(jd+1.d0))
     $              *racad(2*la,1,ja,2,jc,1)*racad(2*lb,1,jb,2,jd,1)
     $              *racad(jc,2,2*J,jb,ja,jd)*(-1)
            else
               term1=0.d0
            endif
            if (na==nd.and.la==ld.and.nb==nc.and.lb==lc) then
               term2=(1.d0+(-1)**(Tab+Tcd))*(-1)**(Tcd-1)
     $              *sqrt((ja+1.d0)*(jb+1.d0)*(jc+1.d0)*(jd+1.d0))
     $              *racad(2*la,1,ja,2,jd,1)*racad(2*lb,1,jb,2,jc,1)
     $              *racad(jd,2,2*J,jb,ja,jc)*(-1)**((jc-jd)/2-J)
            else
               term2=0.d0
            endif

            if (Tab==1.and.Tcd==1) then
               spin_spin(i1,i2,1,1)= ! Sp Sp
     $              1.5d0/sqrt((1.d0+krd(a,b))*(1.d0+krd(c,d)))
     $              *clebd(1,1,1,1,2*Tab,2)*clebd(1,1,1,1,2*Tcd,2)
     $              *((krd(b,d)*krd(a,c)
     $              -krd(b,c)*krd(a,d)*(-1)**(J-(jc+jd)/2+Tcd-1))
     $              +term1-term2)
            endif

            spin_spin(i1,i2,0,1)=            ! Sp Sp
     $           1.5d0/sqrt((1.d0+krd(a,b))*(1.d0+krd(c,d)))
     $           *clebd(1,1,1,-1,2*Tab,0)*clebd(1,1,1,-1,2*Tcd,0)
     $           *(krd(b,d)*krd(a,c)
     $           -krd(b,c)*krd(a,d)*(-1)**(J-(jc+jd)/2+Tcd-1))

            spin_spin(i1,i2,0,0)=            ! Sp Sn
     $           1.5d0/sqrt((1.d0+krd(a,b))*(1.d0+krd(c,d)))
     $           *clebd(1,1,1,-1,2*Tab,0)*clebd(1,1,1,-1,2*Tcd,0)
     $           *(term1-term2)

            if (Tab==1.and.Tcd==1) then
               spin_spin(i1,i2,-1,-1)= ! Sn Sn
     $              1.5d0/sqrt((1.d0+krd(a,b))*(1.d0+krd(c,d)))
     $              *clebd(1,-1,1,-1,2*Tab,-2)*clebd(1,-1,1,-1,2*Tcd,-2)
     $              *((krd(b,d)*krd(a,c)
     $              -krd(b,c)*krd(a,d)*(-1)**(J-(jc+jd)/2+Tcd-1))
     $              +term1-term2)
            endif

            spin_spin(i1,i2,0,-1)=            ! Sn Sn
     $           1.5d0/sqrt((1.d0+krd(a,b))*(1.d0+krd(c,d)))
     $           *clebd(1,-1,1,1,2*Tab,0)*clebd(1,-1,1,1,2*Tcd,0)
     $           *(krd(b,d)*krd(a,c)
     $           -krd(b,c)*krd(a,d)*(-1)**(J-(jc+jd)/2+Tcd-1))

         end do
      end do
      end subroutine spindotspin

      real(kind(0.d0)) function krd(i,j)
      implicit none
      integer,intent(IN) :: i,j
      if (i==j) then
         krd=1.d0
      else
         krd=0.d0
      endif
      end function krd

      end

      subroutine cluster_overlap_testread
      use obdens
      use initst
      use finast
      use intrface
      use constants
      use paramdef
      use gamalog
      implicit none
      logical :: mfdp
      integer(4) nhme,k1max,mxsps,major
      real(kind=kind(0.d0)) :: xxx,sumto
      real(kind=kind(0.d0)) :: clebd,cleb,clebt,racad
      real(kind=kind(0.d0)),allocatable:: sumN(:,:)
      integer :: jtotali2,jtotalf2,kii,kff
      integer :: jtotal2,jtotal2min,itotali2,itotalf2
      integer :: N1in,j1,j2,k1,ii,ij
      integer :: kiiin,kffin,i1,i3
      integer :: j12,it12,ia,ixxx

      open(iunitvout,file='trdens.out',status='old',action='read')
      open(iunitobout,file='observ.out',status='unknown')
      read(iunitvout,*)
      read(iunitvout,*) mfdp
      read(iunitvout,*)

      allocate(nsdi,mxnwdi,nhomi,mjtotali,mttotali,iparityi,nhwi,
     +     nucleonsi,nprotonsi,nneutrnsi,hboi,nhom12i,
     +     nhom123i,nhom1234i,naspsi)

      if (mfdp) then
         print *,' *** input error in cluster_overlap_testread'
         stop
      else

         allocate(nsdf,mxnwdf,nhomf,mjtotalf,mttotalf,iparityf,nhwf,
     +        nucleonsf,nprotonsf,nneutrnsf,hbof,nhom12f,
     +        nhom123f,nhom1234f,naspsf)

         read(iunitvout,*) 
         read(iunitvout,*) 
         read(iunitvout,*) 
         read(iunitvout,*) 

         read(iunitvout,*)
         read(iunitvout,*)
         read(iunitvout,1001) nucleonsi,nprotonsi,nneutrnsi,mjtotali,
     +        mttotali,hboi,nhwi,nsdi,nhme,k1max,
     +        mxnwdi,mxsps,major,iparityi 

 1001    format(9x,/,3x,i3,5x,i3,5x,i3,
     + /,6x,i3,8x,i3,11x,/,
     + 12x, f8.4, 7x, i3, 13x, i8, 8x,
     + i10,/, 7x, i3, 9x, i3, 9x, i8,
     + 9x, i2, 11x, i2,/)

         if (iparityi==0) then
         write(iunitobout,1000) nucleonsi,nprotonsi,nneutrnsi,mjtotali,
     + mttotali,hboi,nhwi,nsdi,nhme,k1max,
     + mxnwdi,mxsps,major,iparityi 
 1000    format(' Nucleus:',/,' A=',i3,'   Z=',i3,'   N=',i3,
     + /,' 2*MJ=',i3,'   2*MT=',i3,'  parity= +',/,
     + ' hbar Omega=',f8.4,'   Nhw=',i3,'   dimension=',i8,'   nhme=',
     + i10,/,' k1max=',i3,'   mxnwd=',i3,'   mxsps=',i8,
     +                    '   major=',i2,'   iparity=',i2,/)
         else
         write(iunitobout,1007) nucleonsi,nprotonsi,nneutrnsi,mjtotali,
     + mttotali,hboi,nhwi,nsdi,nhme,k1max,
     + mxnwdi,mxsps,major,iparityi 
 1007    format(' Nucleus:',/,' A=',i3,'   Z=',i3,'   N=',i3,
     + /,' 2*MJ=',i3,'   2*MT=',i3,'  parity= -',/,
     + ' hbar Omega=',f8.4,'   Nhw=',i3,'   dimension=',i8,'   nhme=',
     + i10,/,' k1max=',i3,'   mxnwd=',i3,'   mxsps=',i8,
     +                    '   major=',i2,'   iparity=',i2,/)
         endif
         allocate(Jxi(k1max+2),Txi(k1max+2),jt2i(k1max+2),
     +        it2i(k1max+2),eneri(k1max+2))
         read(iunitvout,1101) (Jxi(k1),Txi(k1),eneri(k1),xxx,
     +        k1=1,k1max+2)         
 1101    format(3x, f7.4, 6x, f7.4, 12x, f12.4, 8x, f12.4) 
         write(iunitobout,1100) (Jxi(k1),Txi(k1),eneri(k1),
     +        eneri(k1)-eneri(1),k1=1,k1max+2)         
 1100    format(' J=',f7.4,'    T=',f7.4,'     Energy=',f12.4,
     + '     Ex=',f12.4) 

         do ixxx=1,k1max+2
            if (mod(nucleonsi,2)==0) then
               jt2i(ixxx)=2*nint(Jxi(ixxx))
               it2i(ixxx)=2*nint(Txi(ixxx))
            else
               jt2i(ixxx)=nint(2*Jxi(ixxx))
               if (mod(jt2i(ixxx),2)/=1) then
                  if (2*Jxi(ixxx)>jt2i(ixxx)) then
                     jt2i(ixxx)=nint(2*Jxi(ixxx))+1
                  else
                     jt2i(ixxx)=nint(2*Jxi(ixxx))-1
                  endif
               endif
               it2i(ixxx)=nint(2*Txi(ixxx))
               if (mod(it2i(ixxx),2)/=1) then
                  if (2*Txi(ixxx)>it2i(ixxx)) then
                     it2i(ixxx)=nint(2*Txi(ixxx))+1
                  else
                     it2i(ixxx)=nint(2*Txi(ixxx))-1
                  endif
               endif
            endif
            if (dabs(dble(jt2i(ixxx))-2.d0*Jxi(ixxx))>1.d-4) then
               print *, ixxx,'**** not a good J ****'
            endif 
         end do

         read(iunitvout,1157) nhomi,nhom12i,naspsi
 1157    format(/, 8x, i4,11x,i4,9x,i4)
         write(iunitobout,1156) nhomi,nhom12i,naspsi
 1156    format(/,' N1_max=',i4,'   N12_max=',i4,'   Nasps=',i4)

         read(iunitvout,1348) ki,nki
 1348    format(/, 31x, i3, 3x, i3, 5x)
         nki=nki+1-ki
         write(iunitobout,1347) ki,ki+nki-1
 1347    format(/,' wave functions of the states #',i3,'- #',i3,
     +        ' used')

         read(iunitvout,*) 
         read(iunitvout,*) 
         read(iunitvout,*) 
         read(iunitvout,*) 

         read(iunitvout,*)
         read(iunitvout,*)
         read(iunitvout,1001) nucleonsf,nprotonsf,nneutrnsf,mjtotalf,
     +        mttotalf,hbof,nhwf,nsdf,nhme,k1max,
     +        mxnwdf,mxsps,major,iparityf 

         if (iparityf==0) then
         write(iunitobout,1000) nucleonsf,nprotonsf,nneutrnsf,mjtotalf,
     +        mttotalf,hbof,nhwf,nsdf,nhme,k1max,
     +        mxnwdf,mxsps,major,iparityf 
         else
         write(iunitobout,1007) nucleonsf,nprotonsf,nneutrnsf,mjtotalf,
     +        mttotalf,hbof,nhwf,nsdf,nhme,k1max,
     +        mxnwdf,mxsps,major,iparityf 
         endif

         allocate(Jxf(k1max+2),Txf(k1max+2),jt2f(k1max+2),
     +        it2f(k1max+2),enerf(k1max+2))
         read(iunitvout,1101) (Jxf(k1),Txf(k1),enerf(k1),xxx,
     +        k1=1,k1max+2)         

         write(iunitobout,1100) (Jxf(k1),Txf(k1),enerf(k1),
     +        enerf(k1)-enerf(1),k1=1,k1max+2)         

         do ixxx=1,k1max+2
            if (mod(nucleonsf,2)==0) then
               jt2f(ixxx)=2*nint(Jxf(ixxx))
               it2f(ixxx)=2*nint(Txf(ixxx))
            else
               jt2f(ixxx)=nint(2*Jxf(ixxx))
               if (mod(jt2f(ixxx),2)/=1) then
                  if (2*Jxf(ixxx)>jt2f(ixxx)) then
                     jt2f(ixxx)=nint(2*Jxf(ixxx))+1
                  else
                     jt2f(ixxx)=nint(2*Jxf(ixxx))-1
                  endif
               endif
               it2f(ixxx)=nint(2*Txf(ixxx))
               if (mod(it2f(ixxx),2)/=1) then
                  if (2*Txf(ixxx)>it2f(ixxx)) then
                     it2f(ixxx)=nint(2*Txf(ixxx))+1
                  else
                     it2f(ixxx)=nint(2*Txf(ixxx))-1
                  endif
               endif
            endif
            if (dabs(dble(jt2f(ixxx))-2.d0*Jxf(ixxx))>1.d-4) then
               print *, ixxx,'**** not a good J ****'
            endif 
         end do

         read(iunitvout,1157) nhomf,nhom12f,naspsf
         write(iunitobout,1156) nhomf,nhom12f,naspsf

         nhomi=max(nhomi,nhomf)
         nhomf=>nhomi
         nhom12i=max(nhom12i,nhom12f)
         nhom12f=>nhom12i

         select case(iabs(nucleonsi-nucleonsf))
         case(0)
            print *,' *** input error in cluster_overlap_testread'
            write(iunitobout,*) '*** incompatible number of nucleons'
            stop
         case(1:4)
            continue
         case default
            write(iunitobout,*) 
     +           ' nucleonsi and nucleonsf out of range',
     +           nucleonsi,nucleonsf
            stop
         end select

      endif
   
      read(iunitvout,1348) kf,nkf
      nkf=nkf+1-kf
      write(iunitobout,1347) kf,kf+nkf-1

      read(iunitvout,*)
      write(iunitobout,*)
      read(iunitvout,7769) ist1dim
      write(iunitobout,7768) ist1dim
 7768 format(' number of single-nucleon states =',i4)
 7769 format(34x,i4)
      allocate(ist1(3,ist1dim))
      do i1=1,ist1dim
         read(iunitvout,7771) ii,(ist1(i3,i1),i3=1,3)
 7771    format(2x,i4,4x,i3,4x,i3,4x,i2,2x)
         if (ii/=i1) then
            print *,'*** error: i1,ii=',i1,ii
            stop
         endif   
      end do
      do ii=1,ist1dim
         write(iunitobout,7770) ii,(ist1(i3,ii),i3=1,3)
 7770    format(' #',i4,'  n=',i3,'  l=',i3,'  j=',i2,'/2')
      end do

      select case(iabs(nucleonsi-nucleonsf))
      case(1)
         allocate(ad_cl(ist1dim,kf:kf+nkf-1,ki:ki+nki-1))
         ad_cl=0.d0
      case(2)
         read(iunitvout,*)
         write(iunitobout,*)
         read(iunitvout,"(8x,i3,13x,i3)") 
     +        nhomi,nhom12i
         write(iunitobout,"(' Ni_max=',i3,'   Ni_12_max=',i3)") 
     +           nhomi,nhom12i
         read(iunitvout,"(8x,i3,13x,i3)") 
     +        nhomf,nhom12f
         write(iunitobout,"(' Nf_max=',i3,'   Nf_12_max=',i3)") 
     +        nhomf,nhom12f
         read(iunitvout,*)
         write(iunitobout,*)
         read(iunitvout,"(31x,i5)") 
     +        ist2dim
         write(iunitobout,"(' number of two-nucleon states =',i5)") 
     +        ist2dim
         allocate(ist2(4,ist2dim))
         allocate(ad_cl(ist2dim,kf:kf+nkf-1,ki:ki+nki-1))
         ad_cl=0.d0
         read(iunitvout,"(17x,2i4)") 
     +        jtotal2min,jtotal2
         write(iunitobout,"(' J12_min,J12_max=',2i4)") 
     +        jtotal2min,jtotal2
         read(iunitvout,*)
         do ij=1,ist2dim
            read(iunitvout,"(i5,2i4,2i3)") 
     +           ii,(ist2(ia,ij),ia=1,4) 
         end do
      case(3)
         read(iunitvout,*)
         write(iunitobout,*)
         read(iunitvout,"(8x,i3,13x,i3,14x,i3)") 
     +        nhomi,nhom12i,nhom123i
         write(iunitobout,"(' Ni_max=',i3,'   Ni_12_max=',i3,
     +        '   Ni_123_max=',i3)") 
     +        nhomi,nhom12i,nhom123i
         read(iunitvout,"(8x,i3,13x,i3,14x,i3)") 
     +        nhomf,nhom12f,nhom123f
         write(iunitobout,"(' Nf_max=',i3,'   Nf_12_max=',i3,
     +        '   Nf_123_max=',i3)")  
     +        nhomf,nhom12f,nhom123f
         read(iunitvout,*)
         write(iunitobout,*)
         read(iunitvout,"(31x,i5)") ist2dim
         write(iunitobout,"(' number of two-nucleon states =',i5)") 
     +        ist2dim
         allocate(ist2(4,ist2dim))
         read(iunitvout,*)
         do ij=1,ist2dim
            read(iunitvout,"(i5,2i4,2i3)") 
     +           ii,(ist2(ia,ij),ia=1,4) 
         end do
         read(iunitvout,*)
         write(iunitobout,*)
         read(iunitvout,"(33x,i7)") ist3dim
         write(iunitobout,"(' number of three-nucleon states =',i7)") 
     +        ist3dim
         allocate(ist3(4,ist3dim))
         allocate(ad_cl(ist3dim,kf:kf+nkf-1,ki:ki+nki-1))
         ad_cl=0.d0
         read(iunitvout,"(19x,2i4)") 
     +        jtotal2min,jtotal2
         write(iunitobout,"(' J123_min,J123_max=',2i4)")
     +        jtotal2min,jtotal2
         read(iunitvout,*)
         do ij=1,ist3dim
            read(iunitvout,"(i7,i5,i4,2i3)") 
     +           ii,(ist3(ia,ij),ia=1,4) 
         end do
      case(4)
         read(iunitvout,*)
         write(iunitobout,*)
         read(iunitvout,"(8x,i3,13x,i3,14x,i3,15x,i3)") 
     +        nhomi,nhom12i,nhom123i,nhom1234i
         write(iunitobout,"(' Ni_max=',i3,'   Ni_12_max=',i3,
     +        '   Ni_123_max=',i3,'   Ni_1234_max=',i3)") 
     +        nhomi,nhom12i,nhom123i,nhom1234i
         read(iunitvout,"(8x,i3,13x,i3,14x,i3,15x,i3)") 
     +        nhomf,nhom12f,nhom123f,nhom1234f
         write(iunitobout,"(' Nf_max=',i3,'   Nf_12_max=',i3,
     +        '   Nf_123_max=',i3,'   Nf_1234_max=',i3)")  
     +        nhomf,nhom12f,nhom123f,nhom1234f
         read(iunitvout,*)
         write(iunitobout,*)
         read(iunitvout,"(31x,i5)") ist2dim
         write(iunitobout,"(' number of two-nucleon states =',i5)") 
     +        ist2dim
         allocate(ist2(4,ist2dim))
         read(iunitvout,*)
         do ij=1,ist2dim
            read(iunitvout,"(i5,2i4,2i3)") 
     +           ii,(ist2(ia,ij),ia=1,4) 
         end do
         read(iunitvout,*)
         write(iunitobout,*)
         read(iunitvout,"(33x,i7)") ist3dim
         write(iunitobout,"(' number of three-nucleon states =',i7)") 
     +        ist3dim
         allocate(ist3(4,ist3dim))
         read(iunitvout,*)
         do ij=1,ist3dim
            read(iunitvout,"(i7,i5,i4,2i3)") 
     +           ii,(ist3(ia,ij),ia=1,4) 
         end do
         read(iunitvout,*)
         write(iunitobout,*)
         read(iunitvout,"(32x,i8)") ist4dim
         write(iunitobout,"(' number of four-nucleon states =',i8)") 
     +        ist4dim
         allocate(ist4(4,ist4dim))
         allocate(ad_cl(ist4dim,kf:kf+nkf-1,ki:ki+nki-1))
         ad_cl=0.d0
         read(iunitvout,"(21x,2i4)") 
     +        jtotal2min,jtotal2
         write(iunitobout,"(' J1234_min,J1234_max=',2i4)") 
     +        jtotal2min,jtotal2
         read(iunitvout,*)
         do ij=1,ist4dim
            read(iunitvout,"(i8,i7,i5,2i3)") 
     +           ii,(ist4(ia,ij),ia=1,4) 
         end do
      case default
         write(iunitobout,*) 
     +        ' nucleonsi and nucleonsf out of range',
     +        nucleonsi,nucleonsf
         stop
      end select

      do kii=ki,ki+nki-1

         jtotali2=jt2i(kii)
         itotali2=it2i(kii)
         if (dabs(dble(jtotali2)-2.d0*Jxi(kii))>1.d-4) then
            print *, '**** not a good J ****'
         endif 

         do kff=kf,kf+nkf-1

            jtotalf2=jt2f(kff)
            itotalf2=it2f(kff)
            if (dabs(dble(jtotalf2)-2.d0*Jxf(kff))>1.d-4) then
               print *, '**** not a good J ****'
            endif 

            read(iunitvout,1202) kffin,jtotalf2,itotalf2,xxx,
     +           kiiin,jtotali2,itotali2,sumto
 1200       format(//,' *** Transition matrix elements for states:',
     +           ' ***',/,' #',i3,
     +           ' [2*(J,T),Ex]_f=',i3,i2,f8.4,3x,
     +           '#',i3,' [2*(J,T),Ex]_i=',i3,i2,f8.4) 
 1202       format(//, 47x,/, 2x, i3,
     +           16x, i3, i2, f8.4, 3x,
     +           1x, i3, 16x, i3, i2, f8.4) 
c* 1200       format(//,' *** Transition matrix elements for states:',
c*     +           ' ***',/,' #',i2,
c*     +           '  [2*(J,T),Ex]_f=',2i3,f8.4,4x,
c*     +           '#',i2,'  [2*(J,T),Ex]_i=',2i3,f8.4) 
c* 1200       format(//,' *** Transition matrix elements for states:',
c*     +           ' ***',/,' #',i2,
c*     +           '  2*(J,T)_f=',2i3,'  Ex_f=',f8.4,4x,
c*     +           '#',i2,'   2*(J,T)_i=',2i3,'  Ex_i=',f8.4) 
c* 1200       format(//,' *** Transition matrix elements for states:',
c*     +           ' ***',/,' #',i3,
c*     +           '   2*J_f,2*T_f=',2i3,'   Ex_f=',f9.4,5x,
c*     +           ' #',i3,'   2*J_i,2*T_i=',2i3,'   Ex_i=',f9.4)

            if (kffin/=kff.or.kiiin/=kii) then
               print *,'*** error: kiiin,kffin=',kiiin,kffin,
     +              '  kii,kff=',kii,kff
               stop
            endif   

            select case(iabs(nucleonsi-nucleonsf))
            case(1)
               if (nucleonsi==nucleonsf+1) then
                  read(iunitvout,*)
                  do i1=1,ist1dim
                     j1=ist1(3,i1)    
                     cleb=clebd(jtotalf2,mjtotalf,j1,
     +                    mjtotali-mjtotalf,jtotali2,mjtotali)
cc                     if (cleb==0.d0) cycle
                     if (abs(cleb)<1.d-12) cycle
                     clebt=clebd(itotalf2,mttotalf,1,
     +                    mttotali-mttotalf,itotali2,mttotali) 
cc                     if (clebt==0.d0) cycle
                     if (abs(clebt)<1.d-12) cycle
                     read(iunitvout,'(4x,i3,13x,f10.6)') 
     +                    N1in,ad_cl(i1,kff,kii)
                     if (N1in/=i1) then
                        print *,'*** error: N1in=',
     +                       N1in,'  i1=',i1
                        stop
                     endif  
                  end do   
               elseif (nucleonsi==nucleonsf-1) then
                  read(iunitvout,*)
                  do i1=1,ist1dim
                     j1=ist1(3,i1)    
                     cleb=clebd(jtotali2,mjtotali,j1,
     +                    mjtotalf-mjtotali,jtotalf2,mjtotalf)
cc                     if (cleb==0.d0) cycle
                     if (abs(cleb)<1.d-12) cycle
                     clebt=clebd(itotali2,mttotali,1,
     +                    mttotalf-mttotali,itotalf2,mttotalf) 
cc                     if (clebt==0.d0) cycle
                     if (abs(clebt)<1.d-12) cycle
                     read(iunitvout,'(4x,i3,13x,f10.6)') 
     +                    N1in,ad_cl(i1,kff,kii)
                     if (N1in/=i1) then
                        print *,'*** error: N1in=',
     +                       N1in,'  i1=',i1
                        stop
                     endif  
                  end do   
               endif   
            case(2)
               if (nucleonsi==nucleonsf+2) then
                  read(iunitvout,*)
                  do i1=1,ist2dim
                     j12=ist2(3,i1)    
                     cleb=clebd(jtotalf2,mjtotalf,2*j12,
     +                    mjtotali-mjtotalf,jtotali2,mjtotali)
cc                     if (cleb==0.d0) cycle
                     if (abs(cleb)<1.d-12) cycle
                     it12=ist2(4,i1)    
                     clebt=clebd(itotalf2,mttotalf,2*it12,
     +                    mttotali-mttotalf,itotali2,mttotali) 
cc                     if (clebt==0.d0) cycle
                     if (abs(clebt)<1.d-12) cycle
                     read(iunitvout,'(6x,i5,12x,f10.6)') 
     +                    N1in,ad_cl(i1,kff,kii)
                     if (N1in/=i1) then
                        print *,'*** error: N1in=',
     +                       N1in,'  i1=',i1
                        stop
                     endif  
                  end do   
               elseif (nucleonsi==nucleonsf-2) then
                  read(iunitvout,*)
                  do i1=1,ist2dim
                     j12=ist2(3,i1)    
                     cleb=clebd(jtotali2,mjtotali,2*j12,
     +                    mjtotalf-mjtotali,jtotalf2,mjtotalf)
cc                     if (cleb==0.d0) cycle
                     if (abs(cleb)<1.d-12) cycle
                     it12=ist2(4,i1)    
                     clebt=clebd(itotali2,mttotali,2*it12,
     +                    mttotalf-mttotali,itotalf2,mttotalf) 
cc                     if (clebt==0.d0) cycle
                     if (abs(clebt)<1.d-12) cycle
                     read(iunitvout,'(6x,i5,12x,f10.6)') 
     +                    N1in,ad_cl(i1,kff,kii)
                     if (N1in/=i1) then
                        print *,'*** error: N1in=',
     +                       N1in,'  i1=',i1
                        stop
                     endif  
                  end do
               endif
            case(3)
               if (nucleonsi==nucleonsf+3) then
                  read(iunitvout,*)
                  do i1=1,ist3dim
                     j12=ist3(3,i1)    
                     cleb=clebd(jtotalf2,mjtotalf,j12,
     +                    mjtotali-mjtotalf,jtotali2,mjtotali)
cc                     if (cleb==0.d0) cycle
                     if (abs(cleb)<1.d-12) cycle
                     it12=ist3(4,i1)    
                     clebt=clebd(itotalf2,mttotalf,it12,
     +                    mttotali-mttotalf,itotali2,mttotali) 
cc                     if (clebt==0.d0) cycle
                     if (abs(clebt)<1.d-12) cycle
                     read(iunitvout,'(8x,i7,14x,f10.6)') 
     +                    N1in,ad_cl(i1,kff,kii)
                     if (N1in/=i1) then
                        print *,'*** error: N1in=',
     +                       N1in,'  i1=',i1
                        stop
                     endif
                  end do   
               elseif (nucleonsi==nucleonsf-3) then
                  read(iunitvout,*)
                  do i1=1,ist3dim
                     j12=ist3(3,i1)    
                     cleb=clebd(jtotali2,mjtotali,j12,
     +                    mjtotalf-mjtotali,jtotalf2,mjtotalf)
cc                     if (cleb==0.d0) cycle
                     if (abs(cleb)<1.d-12) cycle
                     it12=ist3(4,i1)    
                     clebt=clebd(itotali2,mttotali,it12,
     +                    mttotalf-mttotali,itotalf2,mttotalf) 
cc                     if (clebt==0.d0) cycle
                     if (abs(clebt)<1.d-12) cycle
                     read(iunitvout,'(8x,i7,14x,f10.6)') 
     +                    N1in,ad_cl(i1,kff,kii)
                     if (N1in/=i1) then
                        print *,'*** error: N1in=',
     +                       N1in,'  i1=',i1
                        stop
                     endif
                  end do
               endif               
            case(4)
               if (nucleonsi==nucleonsf+4) then
                  read(iunitvout,*)
                  do i1=1,ist4dim
                     j12=ist4(3,i1)    
                     cleb=clebd(jtotalf2,mjtotalf,2*j12,
     +                    mjtotali-mjtotalf,jtotali2,mjtotali)
cc                     if (cleb==0.d0) cycle
                     if (abs(cleb)<1.d-12) cycle
                     it12=ist4(4,i1)    
                     clebt=clebd(itotalf2,mttotalf,2*it12,
     +                    mttotali-mttotalf,itotali2,mttotali) 
cc                     if (clebt==0.d0) cycle
                     if (abs(clebt)<1.d-12) cycle
                     read(iunitvout,'(10x,i8,16x,f10.6)') 
     +                    N1in,ad_cl(i1,kff,kii)
                     if (N1in/=i1) then
                        print *,'*** error: N1in=',
     +                       N1in,'  i1=',i1
                        stop
                     endif
                  end do   
               elseif (nucleonsi==nucleonsf-4) then
                  read(iunitvout,*)
                  do i1=1,ist4dim
                     j12=ist4(3,i1)    
                     cleb=clebd(jtotali2,mjtotali,2*j12,
     +                    mjtotalf-mjtotali,jtotalf2,mjtotalf)
cc                     if (cleb==0.d0) cycle
                     if (abs(cleb)<1.d-12) cycle
                     it12=ist4(4,i1)    
                     clebt=clebd(itotali2,mttotali,2*it12,
     +                    mttotalf-mttotali,itotalf2,mttotalf) 
cc                     if (clebt==0.d0) cycle
                     if (abs(clebt)<1.d-12) cycle
                     read(iunitvout,'(10x,i8,16x,f10.6)') 
     +                    N1in,ad_cl(i1,kff,kii)
                     if (N1in/=i1) then
                        print *,'*** error: N1in=',
     +                       N1in,'  i1=',i1
                        stop
                     endif
                  end do
               endif
            case default
               write(iunitobout,*) 
     +              ' nucleonsi and nucleonsf out of range',
     +              nucleonsi,nucleonsf
               stop
            end select
         end do
      end do

      close(iunitvout)

      if (iabs(nucleonsi-nucleonsf)==1) then
         allocate(sumN(ist1dim,1))
      else
         write(iunitobout,
     +        "(' Reading of <a+-...a+-> matrix elements OK')")
         write(iunitobout,"(
     +'***************************************************************'
     +        )")
         return
      endif

      do kii=ki,ki+nki-1

         jtotali2=jt2i(kii)
         itotali2=it2i(kii)
         if (dabs(dble(jtotali2)-2.d0*Jxi(kii))>1.d-4) then
            print *, '**** not a good J ****'
c            cycle
         endif 

         sumN=0.d0

         do kff=kf,kf+nkf-1
            
            jtotalf2=jt2f(kff)
            itotalf2=it2f(kff)
            if (dabs(dble(jtotalf2)-2.d0*Jxf(kff))>4.d-4) then
               write(iunitobout,*) '**** not a good J ****'
c               cycle
            endif 

            write(iunitobout,1201) kff,jtotalf2,itotalf2,
     +           enerf(kff)-eneri(ki),
     +           kii,jtotali2,itotali2,eneri(kii)-eneri(ki)
 1201       format(/,' *** Transition matrix elements for states:',
     +           ' ***',/,' #',i3,
     +           ' [2*(J,T),Ex]_f=',i3,i2,f8.4,3x,
     +           '#',i3,' [2*(J,T),Ex]_i=',i3,i2,f8.4) 

            if (nucleonsi==nucleonsf+1) then
               write(iunitobout,*)
               do i1=1,ist1dim
                  j1=ist1(3,i1)    
                  cleb=clebd(jtotalf2,mjtotalf,j1,
     +              mjtotali-mjtotalf,jtotali2,mjtotali)
                  if (cleb==0.d0) cycle
                  clebt=clebd(itotalf2,mttotalf,1,
     +                 mttotali-mttotalf,itotali2,mttotali) 
                  if (clebt==0.d0) cycle
                  write(iunitobout,7960) i1,ad_cl(i1,kff,kii)**2
 7960             format(' a-=',i3,'     spf(lj)^2=',f10.6)   
                  sumN(i1,1)=sumN(i1,1)+ad_cl(i1,kff,kii)**2
               end do   
               cycle
            elseif (nucleonsi==nucleonsf-1) then
               write(iunitobout,*)
               do i1=1,ist1dim
                  j1=ist1(3,i1)    
                  cleb=clebd(jtotali2,mjtotali,j1,
     +                 mjtotalf-mjtotali,jtotalf2,mjtotalf)
                  if (cleb==0.d0) cycle
                  clebt=clebd(itotali2,mttotali,1,
     +                 mttotalf-mttotali,itotalf2,mttotalf) 
                  if (clebt==0.d0) cycle
                  write(iunitobout,7962) i1,ad_cl(i1,kff,kii)**2
 7962             format(' a+=',i3,'     spf(lj)^2=',f10.6)     
                  sumN(i1,1)=sumN(i1,1)+dble(jtotalf2+1)
     +                 /dble(jtotali2+1)*ad_cl(i1,kff,kii)**2
               end do   
               cycle
            endif   

         end do

         if (nucleonsi==nucleonsf+1) then
            write(iunitobout,8759) 
 8759       format(/,' Number operator mean values from sum rule')
            do i1=1,ist1dim
               write(iunitobout,2396)i1, sumN(i1,1)
 2396          format(' <n(',i2,')>=',f8.4)
            end do   
            write(iunitobout,*)
     +'***************************************************************'   
         elseif (nucleonsi==nucleonsf-1) then
            write(iunitobout,8759) 
            do i1=1,ist1dim
               j1=ist1(3,i1)    
               write(iunitobout,2396) i1,dble(j1+1)-sumN(i1,1)
            end do
            write(iunitobout,*)
     +'***************************************************************'   
         endif   
      end do
      
      end

      subroutine radialdist
      use obdens
      use initst
      use finast
      use intrface
      use constants
      use harmosc
      implicit none
      real(kind=kind(0.d0)) :: rr,wave,sump,sumn,pintegral,nintegral,
     +     rpintegral,rnintegral
      integer :: kii,kff,jtot2min,jtotal2,jtotal2min,jtotali2,jtotalf2
      integer :: itotali2,itotalf2,jtrans,ii,i1,i2,n1,l1,j1,n2,l2,j2
      integer :: jtot2,numpoint
      real(kind=kind(0.d0)) :: obme,fact,spart,clb,clebd

      rstep=(rinf-rc)/dble(nstep)

      write(iunitobout,"(/,' Radial density')")
      open(42,file='radialdens.dat',status='unknown')
      write(iunitobout,2010) rc,rinf,nstep,rstep
 2010 format(/,' Integration limits:',/,
     + ' rc=',d16.8,'   rinf=',f8.4,'   nstep=',i5,'   rstep=',f8.4)
c***      write(42,"(d16.8,f8.4,i5,f8.4)") rc,rinf,nstep,rstep
      
      jtotal2=(jt2i(ki)+jt2f(kf))/2
      jtotal2min=iabs(jt2i(ki)-jt2f(kf))/2
      do kii=ki,ki+nki-1
         do kff=kf,kf+nkf-1
            jtot2min=iabs(jt2i(kii)-jt2f(kff))/2
            jtot2=(jt2i(kii)+jt2f(kff))/2
            if (jtotal2min>jtot2min) jtotal2min=jtot2min
            if (jtotal2<jtot2) jtotal2=jtot2
         end do
      end do 

      jtotal2=min(jtotal2,jtotal2max)

      print *,' jtotal2min,jtotal2:',jtotal2min,jtotal2

      numpoint=0
      do ii=0,nstep
         rr=rc+ii*rstep
         if (rr>=0.01d0.and.rr<=13.01d0) then
c**            if (ii==10*(ii/10)) numpoint=numpoint+1
            numpoint=numpoint+1
         endif
      end do
      write(42,"(' Number of r-points=',i8)") numpoint

      do kii=ki,ki+nki-1

         jtotali2=jt2i(kii)
         itotali2=it2i(kii)
         if (dabs(dble(jtotali2)-2.d0*Jxi(kii))>1.d-4) then
            print *, '**** not a good J ****'
            cycle
         endif 

         do kff=kf,kf+nkf-1
            
            jtotalf2=jt2f(kff)
            itotalf2=it2f(kff)
            if (dabs(dble(jtotalf2)-2.d0*Jxf(kff))>4.d-4) then
               write(iunitobout,*) '**** not a good J ****'
               cycle
            endif 

            if (iabs(jtotalf2-jtotali2)>2*jtotal2) cycle

            write(iunitobout,1201) kff,jtotalf2,itotalf2,
     +           enerf(kff)-eneri(ki),
     +           kii,jtotali2,itotali2,eneri(kii)-eneri(ki)
 1201       format(/,' *** Transition matrix elements for states:',
     +           ' ***',/,' #',i3,
     +           ' [2*(J,T),Ex]_f=',i3,i2,f8.4,3x,
     +           '#',i3,' [2*(J,T),Ex]_i=',i3,i2,f8.4) 

            do jtrans=max(jtotal2min,iabs(jtotalf2-jtotali2)/2),
     +                min(jtotal2,(jtotalf2+jtotali2)/2)

               if ((-1)**(iparityi+iparityf)/=(-1)**jtrans) cycle

               clb=clebd(jtotali2,mjtotali,2*jtrans,
     +              mjtotalf-mjtotali,jtotalf2,mjtotalf)
               if (dabs(clb)<1.d-8) cycle

               write(iunitobout,1953) jtrans
 1953          format(/,' Jtrans=',i3)
               write(42,4979) kff,jtotalf2,itotalf2,
     +                        kii,jtotali2,itotali2,jtrans 
 4979          format(' Fin st#',i3,'   2*(J,T)_f=',i3,i2,
     +             '     Init st#',i3,'   2*(J,T)_i=',i3,i2,/,
     +             ' Jtrans=',i2)
               write(42,
     +             "('     r      radial proton    radial neutron')")
               if (jtrans==0) then
                  clb=dsqrt((2*jtrans+1.d0)/
     +                 (dble(jtotali2)+1.d0))
                  rpintegral=0.d0
                  rnintegral=0.d0
               else
                  clb=dsqrt((2*jtrans+1.d0)/4.d0/piln/
     +                 (dble(jtotali2)+1.d0))
               endif
               pintegral=0.d0
               nintegral=0.d0
               fact=4.d0
               do ii=0,nstep
                  rr=rc+ii*rstep
                  sump=0.d0
                  sumn=0.d0
                  do i1=1,ist1dim
                     n1=ist1(1,i1)
                     l1=ist1(2,i1)
                     j1=ist1(3,i1)
                     do i2=1,ist1dim
                        n2=ist1(1,i2)
                        l2=ist1(2,i2)
                        j2=ist1(3,i2)

                        if((j1+j2)<2*jtrans) cycle
                        if(iabs(j1-j2)>2*jtrans) cycle
                        if (mod(l1+l2+jtrans,2)/=0) cycle
                        
                        spart=u(ii,n1,l1)*u(ii,n2,l2)
     +                       *obme(n1,l1,j1,n2,l2,j2,jtrans,5)
                        sump=sump+spart
     +                       *tdJp(i1,i2,jtrans,kff,kii)
                        sumn=sumn+spart
     +                       *tdJn(i1,i2,jtrans,kff,kii)

                     end do
                  end do
                  sump=sump*clb
                  sumn=sumn*clb

                  pintegral=pintegral+sump*(rr**jtrans)*fact
                  nintegral=nintegral+sumn*(rr**jtrans)*fact
                  if (jtrans==0) then
                     rpintegral=rpintegral+sump*(rr**2)*fact
                     rnintegral=rnintegral+sumn*(rr**2)*fact
                  endif
                  fact=6.d0-fact
                  if (rr>=0.01d0.and.rr<=13.01d0) then
                     if (ii==50*(ii/50)) write(iunitobout,
     +                    "(' r=',f10.5,'   p=',e14.7,'   n=',e14.7)") 
     +                    rr,sump,sumn
c**                     if (ii==10*(ii/10)) 
c**     +                    write(42,"(f10.5,1x,e14.7,1x,e14.7)") 
                     write(42,"(f10.5,1x,e14.7,1x,e14.7)") 
     +                    rr,sump,sumn
                  endif

               end do
               pintegral=pintegral*rstep/3.d0
               nintegral=nintegral*rstep/3.d0
               write(iunitobout,"(' <E',i1,'>=',f10.5,
     +                          ' <N',i1,'>=',f10.5)") jtrans,
     +              pintegral,jtrans,nintegral
               if (jtrans==0) then
                  rpintegral=rpintegral*rstep/3.d0
                  rnintegral=rnintegral*rstep/3.d0
                  if (nprotonsi/=0.and.nneutrnsi/=0) then
                     if (kii==kff) then
                        write(iunitobout,"('  r_p=',f9.4,
     +                       '   r_n=',f9.4)") 
     +                       dsqrt(rpintegral/dble(nprotonsi)),
     +                       dsqrt(rnintegral/dble(nneutrnsi))
                     else
                        write(iunitobout,"(' <r_p^2>=',f9.4,
     +                       ' <r_n^2>=',f9.4)") 
     +                       rpintegral/dble(nprotonsi),
     +                       rnintegral/dble(nneutrnsi)
                     endif
                  endif
               endif
               write(42,*)
            end do
         end do
      end do

      end


      subroutine radialdist_trinv
      use obdens
      use initst
      use finast
      use intrface
      use constants
      use harmosc
      use gamalog
      implicit none
      real(kind=kind(0.d0)) :: rr,wave,sump,sumn,pintegral,nintegral,
     +     rpintegral,rnintegral
      integer :: kii,kff,jtot2min,jtotal2,jtotal2min,jtotali2,jtotalf2
      integer :: itotali2,itotalf2,jtrans,ii,i1,i2,n1,l1,j1,n2,l2,j2
      integer :: jtot2,numpoint,i,n,l,np,lp,i1_inv_M_K
      real(kind=kind(0.d0)) :: obme,fact,spart,clb,clebd,sqtA
      integer,allocatable :: ist_M_K(:,:),point_ist_M_K(:,:,:,:)
      integer :: dim_M_K
      real(kind=kind(0.d0)),allocatable :: inv_M_K(:,:)

      write(iunitobout,
     +     "(/,' Translationally invariant radial density')")
      open(42,file='radialdens_trinv.dat',status='unknown')
      write(iunitobout,2010) rc,rinf,nstep,rstep
 2010 format(/,' Integration limits:',/,
     + ' rc=',d16.8,'   rinf=',f8.4,'   nstep=',i5,'   rstep=',f8.4)
c***      write(42,"(d16.8,f8.4,i5,f8.4)") rc,rinf,nstep,rstep

      sqtA=sqrt(real(nucleonsi,kind(0.d0))
     +     /real(nucleonsi-1,kind(0.d0)))
      allocate(point_ist_M_K(0:nhom12i/2,0:nhom12i,0:nhom12i/2,
     +     0:nhom12i))
      lrelm=nhom12i
      nrmax=nhom12i/2
c**** call gamasub when lrelm,nrmax changes
      call gamasub
      if (allocated(u)) deallocate(u)
      allocate(u(0:nstep,0:nrmax,0:lrelm))
      rstep=(rinf-rc)/dble(nstep)
      do l=0,lrelm
         do n=0,nrmax
*C$DOACROSS LOCAL(ii,rr,wave),SHARE(n,l,anu,rc,rstep,nstep)
            do ii=0,nstep
               rr=rc+ii*rstep
               call waver(n,l,anu,sqtA*rr,wave)
               u(ii,n,l)=wave
            end do
         end do
      end do   
      print *,' u calculated'
      
      jtotal2=(jt2i(ki)+jt2f(kf))/2
      jtotal2min=iabs(jt2i(ki)-jt2f(kf))/2
      do kii=ki,ki+nki-1
         do kff=kf,kf+nkf-1
            jtot2min=iabs(jt2i(kii)-jt2f(kff))/2
            jtot2=(jt2i(kii)+jt2f(kff))/2
            if (jtotal2min>jtot2min) jtotal2min=jtot2min
            if (jtotal2<jtot2) jtotal2=jtot2
         end do
      end do 

      jtotal2=min(jtotal2,jtotal2max)

      print *,' jtotal2min,jtotal2:',jtotal2min,jtotal2

      numpoint=0
      do ii=0,nstep
         rr=rc+ii*rstep
         if (rr>=0.01d0.and.rr<=13.01d0) then
c**            if (ii==10*(ii/10)) numpoint=numpoint+1
            numpoint=numpoint+1
         endif
      end do
      write(42,"(' Number of r-points=',i8)") numpoint

      do kii=ki,ki+nki-1

         jtotali2=jt2i(kii)
         itotali2=it2i(kii)
         if (dabs(dble(jtotali2)-2.d0*Jxi(kii))>1.d-4) then
            print *, '**** not a good J ****'
c            cycle
         endif 

         do kff=kf,kf+nkf-1
            
            jtotalf2=jt2f(kff)
            itotalf2=it2f(kff)
            if (dabs(dble(jtotalf2)-2.d0*Jxf(kff))>4.d-4) then
               write(iunitobout,*) '**** not a good J ****'
c               cycle
            endif 

            if (iabs(jtotalf2-jtotali2)>2*jtotal2) cycle

            write(iunitobout,1201) kff,jtotalf2,itotalf2,
     +           enerf(kff)-eneri(ki),
     +           kii,jtotali2,itotali2,eneri(kii)-eneri(ki)
 1201       format(/,' *** Transition matrix elements for states:',
     +           ' ***',/,' #',i3,
     +           ' [2*(J,T),Ex]_f=',i3,i2,f8.4,3x,
     +           '#',i3,' [2*(J,T),Ex]_i=',i3,i2,f8.4) 

            do jtrans=max(jtotal2min,iabs(jtotalf2-jtotali2)/2),
     +                min(jtotal2,(jtotalf2+jtotali2)/2)

               if ((-1)**(iparityi+iparityf)/=(-1)**jtrans) cycle

               clb=clebd(jtotali2,mjtotali,2*jtrans,
     +              mjtotalf-mjtotali,jtotalf2,mjtotalf)
               if (dabs(clb)<1.d-8) cycle

               write(iunitobout,1953) jtrans
 1953          format(/,' Jtrans=',i3)
               write(42,4979) kff,jtotalf2,itotalf2,
     +                        kii,jtotali2,itotali2,jtrans 
 4979          format(' Fin st#',i3,'   2*(J,T)_f=',i3,i2,
     +             '     Init st#',i3,'   2*(J,T)_i=',i3,i2,/,
     +             ' Jtrans=',i2)
               write(42,
     +             "('     r      radial proton    radial neutron')")

               call ist_ M_K_init(jtrans)
               if (allocated(inv_M_K)) deallocate(inv_M_K)
               allocate(inv_M_K(dim_M_K,dim_M_K))
               call M_K_matrix_inv(jtrans)

               if (jtrans==0) then
                  clb=dsqrt((2*jtrans+1.d0)/
     +                 (dble(jtotali2)+1.d0))
                  rpintegral=0.d0
                  rnintegral=0.d0
               else
                  clb=((-1)**jtrans)*dsqrt((2*jtrans+1.d0)/4.d0/piln/
     +                 (dble(jtotali2)+1.d0))
               endif
               pintegral=0.d0
               nintegral=0.d0
               fact=4.d0
               do ii=0,nstep
                  rr=rc+ii*rstep
                  sump=0.d0
                  sumn=0.d0
                  do i1=1,ist1dim
                     n1=ist1(1,i1)
                     l1=ist1(2,i1)
                     j1=ist1(3,i1)
                     do i2=1,ist1dim
                        n2=ist1(1,i2)
                        l2=ist1(2,i2)
                        j2=ist1(3,i2)

                        if((j1+j2)<2*jtrans) cycle
                        if(iabs(j1-j2)>2*jtrans) cycle
                        if (mod(l1+l2+jtrans,2)/=0) cycle
                        
                        i1_inv_M_K=point_ist_M_K(n1,l1,n2,l2)

                        do i=1,dim_M_K
                           n=ist_M_K(1,i)
                           l=ist_M_K(2,i)
                           np=ist_M_K(3,i)
                           lp=ist_M_K(4,i)

                           spart=u(ii,n,l)*u(ii,np,lp)
     +                          *inv_M_K(i,i1_inv_M_K)
     +                          *obme(n1,l1,j1,n2,l2,j2,jtrans,5)
                           sump=sump+spart
     +                          *tdJp(i1,i2,jtrans,kff,kii)
                           sumn=sumn+spart
     +                          *tdJn(i1,i2,jtrans,kff,kii)
                           
                        end do

                     end do
                  end do
                  sump=sump*clb*sqtA
                  sumn=sumn*clb*sqtA

                  pintegral=pintegral+sump*(rr**jtrans)*fact
                  nintegral=nintegral+sumn*(rr**jtrans)*fact
                  if (jtrans==0) then
                     rpintegral=rpintegral+sump*(rr**2)*fact
                     rnintegral=rnintegral+sumn*(rr**2)*fact
                  endif
                  fact=6.d0-fact
                  if (rr>=0.01d0.and.rr<=13.01d0) then
                     if (ii==50*(ii/50)) write(iunitobout,
     +                    "(' r=',f10.5,'   p=',e14.7,'   n=',e14.7)") 
     +                    rr,sump,sumn
c**                     if (ii==10*(ii/10)) 
c**     +                    write(42,"(f10.5,1x,e14.7,1x,e14.7)") 
                     write(42,"(f10.5,1x,e14.7,1x,e14.7)") 
     +                    rr,sump,sumn
                  endif

               end do
               pintegral=pintegral*rstep/3.d0
               nintegral=nintegral*rstep/3.d0
               write(iunitobout,"(' <E',i1,'>=',f10.5,
     +                          ' <N',i1,'>=',f10.5)") jtrans,
     +              pintegral,jtrans,nintegral
               if (jtrans==0) then
                  rpintegral=rpintegral*rstep/3.d0
                  rnintegral=rnintegral*rstep/3.d0
                  if (nprotonsi/=0.and.nneutrnsi/=0) then
                     if (kii==kff) then
                        write(iunitobout,"('  r_p=',f9.4,
     +                       '   r_n=',f9.4)") 
     +                       dsqrt(rpintegral/dble(nprotonsi)),
     +                       dsqrt(rnintegral/dble(nneutrnsi))
                     else
                        write(iunitobout,"(' <r_p^2>=',f9.4,
     +                       ' <r_n^2>=',f9.4)") 
     +                       rpintegral/dble(nprotonsi),
     +                       rnintegral/dble(nneutrnsi)
                     endif
                  endif
               endif
               write(42,*)
               if (allocated(inv_M_K)) deallocate(inv_M_K)
            end do
         end do
      end do
      contains
      subroutine ist_M_K_init(K)
      implicit none
      integer,intent(IN):: K
      integer :: n1,l1,n2,l2,ii
      dim_M_K=0
      do n1=0,nhom12i/2
         do l1=0,nhom12i-2*n1
            do n2=0,nhom12i/2
               do l2=0,nhom12i-2*n2
                  if (mod(l1+l2+K,2)==1) cycle
                  if (abs(l1-l2)>K.or.l1+l2<K) cycle
                  dim_M_K=dim_M_K+1
               end do
            end do
         end do
      end do

c**      print *,'K=',K
c**      print *,' dim_M_K=',dim_M_K

      if (allocated(ist_M_K)) deallocate(ist_M_K)
      allocate(ist_M_K(4,dim_M_K))
      point_ist_M_K=0
      ii=0
      do n1=0,nhom12i/2
         do l1=0,nhom12i-2*n1
            do n2=0,nhom12i/2
               do l2=0,nhom12i-2*n2
                  if (mod(l1+l2+K,2)==1) cycle
                  if (abs(l1-l2)>K.or.l1+l2<K) cycle
                  ii=ii+1
                  ist_M_K(1,ii)=n1
                  ist_M_K(2,ii)=l1
                  ist_M_K(3,ii)=n2
                  ist_M_K(4,ii)=l2

c**                  print *,'ii=',ii
c**                  print *,' n1,l1,n2,l2=',n1,l1,n2,l2

                  point_ist_M_K(n1,l1,n2,l2)=ii
               end do
            end do
         end do
      end do
      if (ii/=dim_M_K) then
         print *,'*** error in ist_M_K: dim_M_K,ii=',dim_M_K,ii
         stop
      endif
      end subroutine ist_M_K_init
      subroutine M_K_matrix_inv(K)
      implicit none
      integer,intent(IN):: K
c      integer,allocatable :: liga(:),ligb(:)
      real(kind=kind(0.d0)) :: racad,oscb1,oscb2,dratio,sum,clebd
cc     +     ,sum_test            !,deter
      integer :: n1,l1,n2,l2,n,l,np,lp,cN1,cL1,i1,i,i2
      real(kind=kind(0.d0)),allocatable :: temp(:,:)
      dratio=1.d0/dble(nucleonsi-1)
      inv_M_K=0.d0
      do i1=1,dim_M_K
         n1=ist_M_K(1,i1)
         l1=ist_M_K(2,i1)
         n2=ist_M_K(3,i1)
         l2=ist_M_K(4,i1)
         do i=1,dim_M_K
            n=ist_M_K(1,i)
            l=ist_M_K(2,i)
            np=ist_M_K(3,i)
            lp=ist_M_K(4,i)
            if (2*n+l-2*n1-l1/=2*np+lp-2*n2-l2) cycle

cc            if (K==0) sum_test=0.d0

            sum=0.d0
            do cL1=max(abs(l1-l),abs(l2-lp)),min(l1+l,l2+lp,nhom12i)
               if (mod(l1+l+cL1,2)==1) cycle
               cN1=(2*n+l-2*n1-l1-cL1)/2
               if (cN1<0) exit
               call osclbr(n,l,0,0,cN1,cL1,n1,l1,l,dratio,oscb1)
               call osclbr(np,lp,0,0,cN1,cL1,n2,l2,lp,dratio,oscb2)
               sum=sum+racad(2*l1,2*cL1,2*K,2*lp,2*l,2*l2)
     +              *((-1)**(l2-lp))
     +              *sqrt(real((2*l+1)*(2*lp+1),kind(0.d0)))
     +              *oscb1*oscb2

cc               if (K==0) sum_test=sum_test+oscb1*oscb2

            end do
            inv_M_K(i1,i)=sum
c            print *,' i1,i,M_K=',i1,i,inv_M_K(i1,i)            
cc            if (K==0) then
cc               print *,' i1: n1,l1,n2,l2=',i1,n1,l1,n2,l2
cc               print *,'  i:  n, l,np,lp=',i,n,l,np,lp
cc               print *,' sum_test=',sum_test
cc            endif

         end do
      end do
c      allocate(liga(dim_M_K),ligb(dim_M_K))
      allocate(temp(dim_M_K,dim_M_K))
      temp=inv_M_K
      call inversg(inv_M_K,dim_M_K)
c      call dminv(inv_M_K,dim_M_K,deter,liga,ligb)
c      print *,' deter=',deter
      do i1=1,dim_M_K
         do i2=1,dim_M_K
            sum=0.d0
            do i=1,dim_M_K
               sum=sum+temp(i1,i)*inv_M_K(i,i2)
            end do
            if ((i1==i2.and.dabs(sum-1.d0)>1.d-6)
     +           .or.(i1/=i2.and.dabs(sum)>1.d-6)) then
               write(iunitobout,'(2i3,f10.6)') i1,i2,sum
            endif
         end do
      end do
c      deallocate(liga,ligb)
      deallocate(temp)
      do i1=1,dim_M_K
         n1=ist_M_K(1,i1)
         l1=ist_M_K(2,i1)
         n2=ist_M_K(3,i1)
         l2=ist_M_K(4,i1)
         do i=1,dim_M_K
            n=ist_M_K(1,i)
            l=ist_M_K(2,i)
            np=ist_M_K(3,i)
            lp=ist_M_K(4,i)
            if (2*n+l-2*n1-l1/=2*np+lp-2*n2-l2) cycle
            inv_M_K(i,i1)=inv_M_K(i,i1)
     +           *sqrt(real((2*l+1)*(2*lp+1),kind(0.d0))
     +           /real((2*l1+1)*(2*l2+1),kind(0.d0)))
     +           *clebd(2*l,0,2*lp,0,2*K,0)/clebd(2*l1,0,2*l2,0,2*K,0)
            
cc            if (K==0) then
cc               print *,' i1: n1,l1,n2,l2=',i1,n1,l1,n2,l2
cc               print *,'  i:  n, l,np,lp=',i,n,l,np,lp
cc               print *, ' i,i1,inv_MK=',i,i1,inv_M_K(i,i1)
cc            endif

         end do
      end do
      end subroutine M_K_matrix_inv
      end


      subroutine obformf
      use obdens
      use initst
      use finast
      use intrface
      use constants
      use paramdef
      use gamalog
      use formfparam
      use harmosc
      implicit none
      real(kind=kind(0.d0)) taumax,taustep,qtr,corr,ratio
      real(kind=kind(0.d0)) formC0,formC0W,formC0s,formC,formCW,formCs
      integer iarrq,iq,j1,j2,n1,n2,kii,kff,l1,l2
      integer jtrans,jtotali2,jtotalf2,itotali2,itotalf2
      integer jtotal2,jtotal2min,i1,i2

      taumax=6.0d0
      taustep=0.02d0 
      iarrq=int(taumax/taustep)

      open(22,file='formf.dat',status='unknown')
      WRITE(iunitobout,2000)
 2000 FORMAT(/,' Charge form factor calculation')
      write(iunitobout,2010) rc,rinf,nstep,rstep
 2010 format(/,' Integration limits:',/,
     + ' rc=',d16.8,'   rinf=',f8.4,'   nstep=',i5,'   rstep=',f8.4)

      jtotal2min=0
      jtotal2=0

      allocate(sbf0(0:nstep))
      allocate(sbf1(0:nstep))
      allocate(rm000(0:nhomf,0:nhomi,0:nhomi))
      allocate(rm001(0:nhomf,0:nhomi,0:nhomi))

      do kii=ki,ki+nki-1

         jtotali2=jt2i(kii)
         itotali2=it2i(kii)
         if (dabs(dble(jtotali2)-2.d0*Jxi(kii))>1.d-4) then
            print *, '**** not a good J ****'
            cycle
         endif 

         do kff=kf,kf+nkf-1
            
            jtotalf2=jt2f(kff)
            itotalf2=it2f(kff)
            if (dabs(dble(jtotalf2)-2.d0*Jxf(kff))>4.d-4) then
               write(iunitobout,*) '**** not a good J ****'
               cycle
            endif 

            if (iabs(jtotalf2-jtotali2)>2*jtotal2) cycle

            write(iunitobout,1201) kff,jtotalf2,itotalf2,
     +           enerf(kff)-eneri(ki),
     +           kii,jtotali2,itotali2,eneri(kii)-eneri(ki)
 1201       format(/,' *** Transition matrix elements for states:',
     +           ' ***',/,' #',i3,
     +           ' [2*(J,T),Ex]_f=',i3,i2,f8.4,3x,
     +           '#',i3,' [2*(J,T),Ex]_i=',i3,i2,f8.4) 

            do jtrans=jtotal2min,jtotal2


               write(22,4979) kii,kff 
 4979          format('  Initial state #',i3,'     Final state #',i3)
               write(22,2236)
 2236          format(
     +'  q        (F_C/Z)^2      |F_C|/Z        |F_Cs|',
     +                                  '        F_Cs/F_C',
     +    '        |F_CW|        Gamma')
               do iq=1,iarrq
                  qtr=iq*taustep
                  call m00q(qtr,nhomi)
                  call formfpar(qtr,nucleonsi)
                  formC0=0.d0
                  formC0W=0.d0
                  formC0s=0.d0
                  do i1=1,ist1dim
                     n1=ist1(1,i1)
                     l1=ist1(2,i1)
                     j1=ist1(3,i1)
                     do i2=1,ist1dim
                        n2=ist1(1,i2)
                        l2=ist1(2,i2)
                        j2=ist1(3,i2)

                        if((j1+j2)<2*jtrans) cycle
                        if(iabs(j1-j2)>2*jtrans) cycle
                        if ((-1)**l1/=(-1)**l2) cycle
             
*                  rmulp1=0.d0   
*                  rmuln1=0.d0
*                  rmuls1=0.d0      
 
                        formC0=formC0
     +                       +(rmulp0*rm000(2*n1+l1,2*n2+l2,j1/2)
     +                       +rmulp1*rm001(2*n1+l1,2*n2+l2,j1/2))
     +                       *tdJp(i1,i2,jtrans,kff,kii)
     +                       +(rmuln0*rm000(2*n1+l1,2*n2+l2,j1/2)
     +                       +rmuln1*rm001(2*n1+l1,2*n2+l2,j1/2))
     +                       *tdJn(i1,i2,jtrans,kff,kii)

                        formC0W=formC0W
     +                       +(rmulpW0*rm000(2*n1+l1,2*n2+l2,j1/2)
     +                       +rmulpW1*rm001(2*n1+l1,2*n2+l2,j1/2))
     +                       *tdJp(i1,i2,jtrans,kff,kii)
     +                       +(rmulnW0*rm000(2*n1+l1,2*n2+l2,j1/2)
     +                       +rmulnW1*rm001(2*n1+l1,2*n2+l2,j1/2))
     +                       *tdJn(i1,i2,jtrans,kff,kii)

                        formC0s=formC0s
     +                       +(rmuls0*rm000(2*n1+l1,2*n2+l2,j1/2)
     +                       +rmuls1*rm001(2*n1+l1,2*n2+l2,j1/2))
     +                       *(tdJp(i1,i2,jtrans,kff,kii)
     +                       +tdJn(i1,i2,jtrans,kff,kii))
                     end do
                  end do
                  formC0  = formC0/dsqrt(dble(jtotali2+1))
                  formC0s = formC0s/dsqrt(dble(jtotali2+1))
                  formC0W = formC0W/dsqrt(dble(jtotali2+1))
                  formC  = 2.d0*dsqrt(piln)*formC0
                  formCs = 2.d0*dsqrt(piln)*formC0s
                  formCW = 2.d0*dsqrt(piln)*formC0W
                  corr=-(1.d0+formC0W/(4.d0*sinthetaWsq*formC0))
                  ratio=formC0s/formC0 
                  write(iunitobout,3000) qtr,formC0,formC0s,formC0W
 3000             format(' q=',f8.4,'     F_C0 =',e16.7,'     F_C0s=',
     +                 e16.7,'     F_C0W=',e16.7  )
                  write(22,3002) qtr,(formC/dble(nprotonsi))**2,
     +                 dabs(formC)/dble(nprotonsi),dabs(formCs),ratio,
     +                 dabs(formCW),corr 
 3002             format(f8.4,7e15.6)
                         
               end do       
            end do
         end do
      end do   
      deallocate(rm000)
      deallocate(rm001)
      deallocate(sbf0,sbf1)

      write(iunitobout,*)
     +'***************************************************************'
      close(22)
      end

      subroutine m00q(qtr,nhom)
      use constants
      use harmosc
      implicit double precision (a-h,o-z)

      do i=0,nstep
         r=rc+i*rstep
         rho=qtr*r
         sbf0(i)=sbfj(0,rho)
         sbf1(i)=sbfj(1,rho)/rho
      end do 
      

      nrm=nhom/2
      lrm=nhom
    
      do nrb=0,nrm
         do lrb=0,lrm
            Nb=nrb+nrb+lrb 
            if (Nb>nhom) cycle
            do nra=nrb,nrm
c*               do lra=lrb,lrb
                  lra=lrb
                  Na=nra+nra+lra
                  if (Na>nhom) cycle
                  overx0=0.d0
                  overx1=0.d0
                  fact=4.d0
                  do i=1,nstep-1
                     r=rc+i*rstep
c*                     rho=qtr*r
c*                     wavb=wave(nrb,lrb,anu,r)
c*                     wava=wave(nra,lra,anu,r)
                     wavb=u(i,nrb,lrb)
                     wava=u(i,nra,lra)
c*                     bes0=sbfj(0,rho)
c*                     bes1=sbfj(1,rho)
                     bes0=sbf0(i)
                     bes1=sbf1(i)
                     overx0=overx0+wava*bes0*wavb*fact
                     overx1=overx1+wava*bes1*wavb*fact
                     fact=6.d0-fact
                  enddo
                  over0=(overx0
     +      +u(0,nra,lra)*u(0,nrb,lrb)*sbf0(0)
     +      +u(nstep,nra,lra)*u(nstep,nrb,lrb)*sbf0(nstep))
     +      *rstep/3.d0 
                  over1=(overx1
     +      +u(0,nra,lra)*u(0,nrb,lrb)*sbf1(0)
     +      +u(nstep,nra,lra)*u(nstep,nrb,lrb)*sbf1(nstep))
     +      *rstep/3.d0 
                  do ja=iabs(2*lra-1),2*lra+1,2
                     rja=dble(ja)/2.d0
                     dsqja=dsqrt(dble(ja+1))
                     indj=ja/2
                     rm000(Na,Nb,indj)=dsqja*over0
*                     rm001(Na,Nb,indj)=2.d0*0.75d0*over1   !test for possible error
                     rm001(Na,Nb,indj)=dsqja*over1
     +              *(rja*(rja+1.d0)-lra*(lra+1.d0)-0.75d0)
                      rm000(Nb,Na,indj)=rm000(Na,Nb,indj) 
                      rm001(Nb,Na,indj)=rm001(Na,Nb,indj) 
*                     write(iunitobout,1000) nra,lra,nrb,lrb,ja,Na,Nb,indj,
*     +  over0,over1,rm000(Na,Nb,indj),rm001(Na,Nb,indj)
                  end do
* 1000             format(8i3,4f16.8)
c*               end do
            end do
         end do
      end do 

      end 


      subroutine formfpar(qtr,nucleons)
      use constants
      use harmosc
      use formfparam
      implicit none 
      double precision qtr,tau,GVD,rksin,rksis,GEp,GMp,GEn,GMn
c*      double precision GE0,GM0,GE1,GM1 
      double precision dex,GEs,GMs,GEWp,GMWp,GEWn,GMWn
      integer nucleons

***********
      tau=qtr*qtr*recmass
      GVD=1.d0/((1.d0+rlambdaVD*tau)**2)
      rksin=1.d0/(1.d0+rlambdan*tau)
      rksis=1.d0/(1.d0+rlambdaEs*tau)
      GEp=GVD
      GMp=rmup*GVD
      GEn=-rmun*tau*GVD*rksin
      GMn=rmun*GVD
      GEWp=(1.d0-4.d0*sinthetaWsq)*GEp-GEn
      GEWn=-GEp+(1.d0-4.d0*sinthetaWsq)*GEn
      GMWp=(1.d0-4.d0*sinthetaWsq)*GMp-GMn
      GMWn=-GMp+(1.d0-4.d0*sinthetaWsq)*GMn
      GEs=rhos*tau*GVD*rksis
      GMs=rmus*GVD 

c*      GE0=0.5d0*(GEp+GEn)  
c*      GM0=0.5d0*(GMp+GMn)  
c*      GE1=0.5d0*(GEp-GEn)  
c*      GM1=0.5d0*(GMp-GMn)  

      rmulp0=0.5d0*GEp/dsqrt(piln*(1.d0+tau))
      rmuln0=0.5d0*GEn/dsqrt(piln*(1.d0+tau))
      rmulp1=tau*(GEp-2.d0*GMp)/dsqrt(piln)
      rmuln1=tau*(GEn-2.d0*GMn)/dsqrt(piln)

      rmuls0=0.5d0*GEs/dsqrt(piln*(1.d0+tau))
      rmuls1=tau*(GEs-2.d0*GMs)/dsqrt(piln)

      rmulpW0=0.5d0*GEWp/dsqrt(piln*(1.d0+tau))
      rmulnW0=0.5d0*GEWn/dsqrt(piln*(1.d0+tau))
      rmulpW1=tau*(GEWp-2.d0*GMWp)/dsqrt(piln)
      rmulnW1=tau*(GEWn-2.d0*GMWn)/dsqrt(piln)

c** COM correction: **********
      dex=dexp(bsquare*qtr*qtr/(4.d0*dble(nucleons)))

      rmulp0=rmulp0*dex
      rmuln0=rmuln0*dex
      rmulp1=rmulp1*dex
      rmuln1=rmuln1*dex

      rmuls0=rmuls0*dex
      rmuls1=rmuls1*dex

      rmulpW0=rmulpW0*dex
      rmulnW0=rmulnW0*dex
      rmulpW1=rmulpW1*dex
      rmulnW1=rmulnW1*dex

***********
*      write(2,*) ' COM=',dex
*      write(2,*) ' rmultp0,rmultn0,rmultp1,rmultn1:',
*     + rmultp0,rmultn0,rmultp1,rmultn1
*******
*      rmultp0=1.d0
*      rmultn0=0.d0
*      rmultp1=0.d0
*      rmultn1=0.d0
********

      end  


      double precision function sbfj(k,rho)
      implicit double precision(a-h,o-z)
      select case (k)
      case (0)
         sbfj=dsin(rho)/rho
      case (1)
         xxx=dsin(rho)/rho-dcos(rho)
         sbfj=xxx/rho
      case (2)
         xx1=dsin(rho)
         xxx=xx1/rho-dcos(rho)
         xxx=xxx/rho
         sbfj=(3.d0*xxx-xx1)/rho
*         sbfj=(3.d0/(rho*rho*rho)-1.d0/rho)*dsin(rho)
*     +                      -3.d0*dcos(rho)/(rho*rho)
      case default
         xx1=dsin(rho)
         sbfjlm1=xx1/rho
         xxx=xx1/rho-dcos(rho)
         sbfjl=xxx/rho
*         sbfjl=xx1/(rho*rho)-dcos(rho)/rho
         do l=1,k-1
            sbfjlp1=(2.d0*dble(l)+1.d0)*sbfjl/rho-sbfjlm1
            sbfjlm1=sbfjl
            sbfjl=sbfjlp1
         end do
         sbfj=sbfjl
      end select    
      end 

      double precision function obme(nA13,l3,j23,nA11,l1,j21,kk,ik)
      use harmosc
      implicit double precision (a-h,o-z)

      N23=2*nA13+l3
      N21=2*nA11+l1

      if (ik==0) then
         rL=rLme(nA13,l3,nA11,l1,kk)
         rn=0.5d0*dble((1+(-1)**(l1+l3+kk))*(-1)**kk)
         clb=dsqrt(dble(j23+1))*clebd(j23,1,2*kk,0,j21,1)
         obme=rn*rL*clb
         return
      elseif(ik==2) then   
         if (kk==1) then
            if (N23/=N21.or.l3/=l1) then
               obme=0.d0
               return
            endif
            rn=dsqrt((j23+1.d0)*(j21+1.d0))
            l1l1=l1+l1
            rn=rn*dsqrt(dble(l1*(l1+1)*(l1l1+1)))
            rac=racad(j21,1,2,l1l1,l1l1,j23)
            obme=rn*rac
            return 
         endif
      elseif (ik==3) then
         if (kk==1) then
            if (N23/=N21.or.l3/=l1) then
               obme=0.d0
               return
            endif
            rn=dsqrt((j23+1.d0)*(j21+1.d0))
            l1l1=l1+l1
            rn=rn*dsqrt(1.5d0)
            rac=racad(j23,2,l1l1,1,j21,1)
            obme=rn*rac
            return 
         endif          
      elseif (ik==4) then
         if (kk<3) then
            rL=rLme(nA13,l3,nA11,l1,1)
            rn=dsqrt(dble((j23+1)*(j21+1)*(2*kk+1)*(2*l1+1)*3*6))
            clb=clebd(2*l1,0,2,0,2*l3,0)
     +           *coef9d(2*l1,2,1,2,2*l3,1,j21,2*kk,j23)
            obme=rL*rn*clb
            return
         endif
      elseif (ik==5) then
         rn=0.5d0*dble((1+(-1)**(l1+l3+kk))*(-1)**kk)
         clb=dsqrt(dble(j23+1))*clebd(j23,1,2*kk,0,j21,1)
         obme=rn*clb
         return
      else
         obme=0.d0
         return
      endif
      end

      real(kind=kind(0.d0)) function twobsum(ob3,ob4,
     +              j34,it34,mt34,k,ids,ob1,ob2,j12,it12,mt12)
      use obdens, only: ist1,ist1dim
      implicit none
      integer, intent(IN):: ob3,ob4,j34,it34,mt34,k,ids,
     +                      ob1,ob2,j12,it12,mt12
      integer :: ob6,mt1,mt2,mt3,mt4,mt6
      real(kind=kind(0.d0)) :: sum,clebd,clebt12,clebt34,
     +           obmat1,obmat2,obmat3,obmat4,fac12,fac34,obmat,
     +           racad,rac1,rac2,obme
      integer :: n1,l1,j1,n2,l2,j2,n3,l3,j3,n4,l4,j4,n6,l6,j6

      if (j12/=j34.or.it12/=it34.or.mt12/=mt34) then
         twobsum=0.d0
         return
      endif
      
      n1=ist1(1,ob1)
      l1=ist1(2,ob1)
      j1=ist1(3,ob1)
      n2=ist1(1,ob2)
      l2=ist1(2,ob2)
      j2=ist1(3,ob2)

      n3=ist1(1,ob3)
      l3=ist1(2,ob3)
      j3=ist1(3,ob3)
      n4=ist1(1,ob4)
      l4=ist1(2,ob4)
      j4=ist1(3,ob4)

      if (ob1==ob2) then
         fac12=dsqrt(2.d0)
      else
         fac12=1
      endif
      if (ob3==ob4) then
         fac34=dsqrt(2.d0)
      else
         fac34=1
      endif

      sum=0.d0
      do mt1=-1,1,2
         mt2=2*mt12-mt1
         clebt12=clebd(1,mt1,1,mt2,it12,mt12)
         if (abs(mt2)/=1) cycle
         if (ob1==ob3.and.ob2==ob4) then
            do ob6=1,ist1dim
               n6=ist1(1,ob6)
               l6=ist1(2,ob6)
               j6=ist1(3,ob6)
               do mt6=-1,1,2
                  obmat=0.d0
                  select case(ids)
                  case(0)
                     if (mt1==1.and.mt6==1) then
                        obmat1=obme(n1,l1,j1,n6,l6,j6,k,0)
                     endif
                     if (mt2==1.and.mt6==1) then
                        obmat2=obme(n2,l2,j2,n6,l6,j6,k,0)
                     endif
                  case default
                     continue
                  end select
                  sum=sum+(obmat1/dble(j1+1)+obmat2/dble(j2+1))
     +                     *clebt12**2
               end do
            end do
         endif
         if (ob1==ob4.and.ob2==ob3) then
            do ob6=1,ist1dim
               n6=ist1(1,ob6)
               l6=ist1(2,ob6)
               j6=ist1(3,ob6)
               do mt6=-1,1,2
                  obmat=0.d0
                  select case(ids)
                  case(0)
                     if (mt1==1.and.mt6==1) then
                        obmat1=obme(n1,l1,j1,n6,l6,j6,k,0)
                     endif
                     if (mt2==1.and.mt6==1) then
                        obmat2=obme(n2,l2,j2,n6,l6,j6,k,0)
                     endif
                  case default
                     continue
                  end select
                  sum=sum+(obmat1/dble(j1+1)+obmat2/dble(j2+1))
     +                     *dble((-1)**(j12+it12-(j1+j2)/2))*clebt12**2
               end do
            end do
         endif
         do mt3=-1,1,2
            mt4=2*mt34-mt3
            clebt34=clebd(1,mt3,1,mt4,it34,mt34)
            if (abs(mt2)/=1) cycle
            rac1=racad(j1,j4,j2,j3,2*k,2*j12)*dble((-1)**((j2-j3)/2))
            rac2=racad(j1,j3,j2,j4,2*k,2*j12)
     +           *dble((-1)**(j12-(j2+j3)/2))

         end do
      end do
      twobsum=sum
      end

      subroutine rLmatel(hbo,nhom,jtrmin,jtrmax,ipi)
      use constants
      use harmosc
      use gamalog
      implicit double precision (a-h,o-z)
c*      parameter(nstep=6000,rc=1.d-6,rinf=22.d0)  ! defined in module
      integer,intent(IN) :: jtrmin,jtrmax,ipi
      integer nhom
      real(kind=kind(0.d0)) hbo
c*      double precision,allocatable,dimension(:,:,:):: u   ! declared in module
      lrelm=nhom
      nrmax=nhom/2
c*      rdmavg=2.d0*(fmsp*fmsn)/(fmsp+fmsn) ! defined in module
c*      fmbyh2=rdmavg/hbc**2
      anu=fmbyh2*hbo
      bsquare=1.d0/anu
      bHO=dsqrt(bsquare)
c**** call gamasub before first use
      call gamasub
      allocate(u(0:nstep,0:nrmax,0:lrelm))
      rstep=(rinf-rc)/dble(nstep)
      do lr=0,lrelm
         do nr=0,nrmax
*C$DOACROSS LOCAL(i,rr,wave),SHARE(nr,lr,anu,rc,rstep,nstep)
            do i=0,nstep
               rr=rc+i*rstep
               call waver(nr,lr,anu,rr,wave)
               u(i,nr,lr)=wave
            end do
         end do
      end do   
      print *,' u calculated'

      allocate(rLme(0:nrmax,0:lrelm,0:nrmax,0:lrelm,jtrmin:jtrmax))
      rLme=0.d0
      do lra=0,lrelm
         do nra=0,nrmax
            if (2*nra+lra>nhom) cycle
            do lrb=0,lrelm
c***               if ((-1)**(lra+lrb)/=ipi) cycle
               do nrb=nra,nrmax
                  if (2*nrb+lrb>nhom) cycle
                  do ltr=max(jtrmin,iabs(lra-lrb)),min(jtrmax,lra+lrb)
                     over=0.d0
                     fact=4.d0
                     do i=1,nstep-1
                        rr=rc+i*rstep
                        over=over+(rr**ltr)
     +                       *u(i,nra,lra)*u(i,nrb,lrb)*fact
                        fact=6.d0-fact
                     end do
                     rLme(nra,lra,nrb,lrb,ltr)=(over
     +  +(rc**ltr)*u(0,nra,lra)*u(0,nrb,lrb)
     +  +(rinf**ltr)*u(nstep,nra,lra)*u(nstep,nrb,lrb))*rstep/3.d0
                    rLme(nrb,lrb,nra,lra,ltr)=rLme(nra,lra,nrb,lrb,ltr)
*                    if (ltr==2) then
*                       print *,'na,la,nb,lb,L,<r^L>=',
*     +                      nrb,lrb,nra,lra,ltr,
*     +                      rLme(nra,lra,nrb,lrb,ltr)/bsquare
*                    endif
                  end do
               end do
            end do
         end do
      end do   
         
      end

      subroutine momdistr
      use obdens
      use initst
      use finast
      use intrface
      use constants
      use paramdef
      use gamalog
      use harmosc
      implicit none
      real(kind=kind(0.d0)) taumax,taustep,qtr,scalfac
      real(kind=kind(0.d0)) formp,wave1,wave2,cleb,clebt,clebd
      integer iarrq,iq,j1,j2,n1,n2,kii,kff,l1,l2
      integer jtotali2,jtotalf2,itotali2,itotalf2
      integer i1,i2

      taumax=6.0d0
      taustep=0.02d0 
      iarrq=int(taumax/taustep)

      open(22,file='momdist.dat',status='unknown')
      WRITE(iunitobout,2000)
 2000 FORMAT(/,' Momentum distribution calculation')


      do kii=ki,ki+nki-1

         jtotali2=jt2i(kii)
         itotali2=it2i(kii)
         if (dabs(dble(jtotali2)-2.d0*Jxi(kii))>1.d-4) then
            print *, '**** not a good J ****'
            cycle
         endif 

         do kff=kf,kf+nkf-1
            
            jtotalf2=jt2f(kff)
            itotalf2=it2f(kff)
            if (dabs(dble(jtotalf2)-2.d0*Jxf(kff))>4.d-4) then
               write(iunitobout,*) '**** not a good J ****'
               cycle
            endif 

            write(iunitobout,1201) kff,jtotalf2,itotalf2,
     +           enerf(kff)-eneri(ki),
     +           kii,jtotali2,itotali2,eneri(kii)-eneri(ki)
 1201       format(/,' *** Transition matrix elements for states:',
     +           ' ***',/,' #',i3,
     +           ' [2*(J,T),Ex]_f=',i3,i2,f8.4,3x,
     +           '#',i3,' [2*(J,T),Ex]_i=',i3,i2,f8.4) 

               write(22,4979) kii,kff 
 4979          format('  Initial state #',i3,'     Final state #',i3)
               write(22,2236)
 2236          format('   q         Formp')

               do iq=1,iarrq
                  qtr=iq*taustep

                  formp=0.d0

                  do i1=1,ist1dim
                     n1=ist1(1,i1)
                     l1=ist1(2,i1)
                     j1=ist1(3,i1)
            if (nucleonsi==nucleonsf+1) then
                  cleb=clebd(jtotalf2,mjtotalf,j1,
     +              mjtotali-mjtotalf,jtotali2,mjtotali)
                  if (cleb==0.d0) cycle
                  clebt=clebd(itotalf2,mttotalf,1,
     +                 mttotali-mttotalf,itotali2,mttotali) 
                  if (clebt==0.d0) cycle
                  scalfac=dble(jtotali2+1)
            elseif (nucleonsi==nucleonsf-1) then
                  cleb=clebd(jtotali2,mjtotali,j1,
     +                 mjtotalf-mjtotali,jtotalf2,mjtotalf)
                  if (cleb==0.d0) cycle
                  clebt=clebd(itotali2,mttotali,1,
     +                 mttotalf-mttotali,itotalf2,mttotalf) 
                  if (clebt==0.d0) cycle
                  scalfac=dble(jtotalf2+1)
            endif   
                     call waver(n1,l1,bsquare,qtr,wave1) 
                     wave1=wave1*dble((-1)**n1)/qtr
                     do i2=1,ist1dim
                        n2=ist1(1,i2)
                        l2=ist1(2,i2)
                        j2=ist1(3,i2)

                        if(j1/=j2.or.l1/=l2) cycle
                        call waver(n2,l2,bsquare,qtr,wave2) 
                        wave2=wave2*dble((-1)**n2)/qtr
                        formp=formp+ad_cl(i1,kff,kii)*ad_cl(i2,kff,kii)
     +                       *wave1*wave2*clebt**2

                     end do
                  end do
                  formp  = formp*scalfac/dble(jtotali2+1)/4.d0/piln
*                  formC  = 2.d0*dsqrt(piln)*formC0
                  write(iunitobout,3000) qtr,formp
 3000             format(' q=',f8.4,'     Formp =',e16.7)
                  write(22,3002) qtr*hbc,formp/(hbc**3)
 3002             format(f10.4,e15.6)
                         
            end do
         end do
      end do   

      write(iunitobout,*)
     +'***************************************************************'
      close(22)
      end

      subroutine gamasub
**      use paramdef
      use gamalog
      implicit none
      integer :: i1,il,in
      lmax=max(30,lrelm)
      maxgam=2*nrmax+2*lmax+3
      if (allocated(gamal)) deallocate(gamal)
      if (allocated(dsq)) deallocate(dsq)
      allocate(dsq(nrmax,0:max(lrelm,lmax)))
      allocate(gamal(maxgam))
      gamal(2)=0.d0
      gamal(1)=0.5d0*dlog(3.14159265358979312d0)
      do i1=3,maxgam
         gamal(i1)=dlog(dble(float(i1))/2.d0-1.d0)+gamal(i1-2)
      end do
      do il=0,max(lrelm,lmax)
         do in=1,nrmax
            dsq(in,il)=dsqrt(dble(in)*(dble(il+in)+0.5d0))
         end do
      end do
      end


      subroutine waver(n,l,anu,r,wave)
      use gamalog
      implicit double precision (a-h,o-z) 
c     the radial part wave function (wave=rr(r)) int wave(r)**2 dr =1)
c     anu is the size parameter of nuclear well
c     anu=mass times omega over h bar

      dlanu=dlog(anu)
      zz=anu*r*r
      wavel=0.25d0*dlanu-zz/2.d0
     &     +dble(l+1)*(0.5d0*dlanu+dlog(r))
      if (n==0) then
         guerp=dexp(0.5d0*(dlog(2.d0)-gamal(2*l+3)))
      elseif (n==1) then
         guerp=dexp(0.5d0*(dlog(2.d0)-gamal(2*l+5)))
     +                   *(dble(l)+1.5d0-zz)         
      else   
         guerp=dexp(0.5d0*(dlog(2.d0)-gamal(2*l+5)))
     +                   *(dble(l)+1.5d0-zz)         
         a=dexp(0.5d0*(dlog(2.d0)-gamal(2*l+3)))
         do nnn=2,n
            b=((dble(l+2*nnn)-0.5d0-zz)*guerp
     +           -a*dsq(nnn-1,l))
     +           /dsq(nnn,l)
*     +           -a*dsqrt(dble(nnn-1)*(dble(l+nnn)-0.5d0)))
*     +           /dsqrt(dble(nnn)*(dble(l+nnn)+0.5d0))
            a=guerp
            guerp=b
         end do   
      endif 
      wave=dexp(wavel)*guerp
      end 


      subroutine osclbr(nr,lr,n,l,n1,l1,n2,l2,lambda,dratio,oscb)
**      use paramdef
      use gamalog
c**** calculates generalized ocillator bracket according to
c**** L. Trlifaj, Phys. Rev. C5 (72) 1534.
c**** nr, lr relative quantum numbers
c**** n, l cms quantum numbers
c**** n1, l1, n2, l2 sp quantum numbers
c**** lambda total orbital momentum
c**** dratio=m2/m1 mass ratio of the particles 2 and 1
c**** Petr Navratil, University of Arizona *******
c**** September 1996 *****************************
      implicit double precision (a-h,o-z)
c**** call gamasub before first use
      n12=n1+n1
      n22=n2+n2
      nrnr=nr+nr
      nn=n+n
      ntotr=nrnr+lr+nn+l
      ntotsp=n12+l1+n22+l2
      if (ntotr.ne.ntotsp) then
         oscb=0.d0
         return
      elseif (lambda.lt.iabs(lr-l).or.lambda.gt.(lr+l)) then
         oscb=0.d0
         return
      elseif (lambda.lt.iabs(l1-l2).or.lambda.gt.(l1+l2)) then
         oscb=0.d0
         return
      endif
      l12=l1+l1 
      l22=l2+l2
      lambda2=lambda+lambda
      ll=l+l
      lrlr=lr+lr 
      sq1dra=dsqrt(1.d0+dratio)
      sum=0.d0
      do 100 la1=0,lr
         la2=lr-la1   
         la12=la1+la1
         la22=la2+la2
         sqbin=0.5d0*(gamal(4*lr+2)-gamal(4*la1+2)-gamal(4*la2+2))
         ip1min=iabs(l1-la1)
         ip1max=l1+la1
         ip2min=iabs(l2-la2)
         ip2max=l2+la2
         do 200 ip1=ip1min,ip1max,2
            ip12=ip1+ip1
*            if ((ip1+la1+l1).ne.((ip12+la12+l12)/2)) goto 200
            cp1la1l1=clebd(ip12,0,la12,0,l12,0)
            it1max=n1-(la1+ip1-l1)/2
            do 300 ip2=ip2min,ip2max,2
            ip22=ip2+ip2
*            if ((ip2+la2+l2).ne.((ip22+la22+l22)/2)) goto 300   
            if ((ip1+ip2+l).ne.((ip12+ip22+ll)/2)) goto 300   
            cp2la2l2=clebd(ip22,0,la22,0,l22,0)
            cp1p2l=clebd(ip12,0,ip22,0,ll,0)
            if (cp1p2l.eq.0.d0) goto 300
            coe9=coef9d(ip12,ip22,la12,la22,ll,lrlr,l12,l22,lambda2)
            if (coe9.eq.0.d0) goto 300
            it2max=n2-(la2+ip2-l2)/2
            do it1=0,it1max
               it1it1=it1+it1
               gamt1=gamal(ip12+it1it1+3)+gamal(it1max+it1max-it1it1+2)
     +   +gamal(it1it1+2)
               do 400 it2=0,it2max
                  it2it2=it2+it2
                  iptl=ip1+ip2-l+it1it1+it2it2
                  if (iptl.lt.0) goto 400
                  iptln=iptl-nn
                  if (iptln.lt.0) goto 400
                gamt2=gamal(ip22+it2it2+3)+gamal(it2max+it2max-it2it2+2)     
     +   +gamal(it2it2+2)
                  iupt1t2=ip1+ip2+l+it1it1+it2it2+3
                  if (iupt1t2.gt.maxgam) then 
                   print *,' error: iupt1t2=',iupt1t2,' > maxgam',maxgam
                     stop
                  endif
                  gamtot=gamal(iupt1t2)-gamt1-gamt2+sqbin
     +   +gamal(iptl+2)-gamal(iptln+2)
                  gamtot=dexp(gamtot)
                  iph=la1+it1+it2+(ip1+ip2+l)/2
                  phase=(-1)**iph
                  dfactor=(dsqrt(dratio))**(la1+ip2+it2it2)
                  dfactor=dfactor/(sq1dra**(ip1+ip2+it1it1+it2it2))
               sum=sum+(ip12+1.d0)*(ip22+1.d0)*cp1la1l1*cp2la2l2*cp1p2l
     +   *coe9*gamtot*phase*dfactor
  400          continue
            end do 
  300       continue
  200    continue
  100 continue
      gamtot=0.5d0*(gamal(n12+l12+3)+gamal(n22+l22+3)-gamal(nn+ll+3)
     + -gamal(nrnr+lrlr+3)+gamal(nrnr+2)+gamal(n12+2)+gamal(n22+2)
     + -gamal(nn+2))+gamal(1)
      gamtot=dexp(gamtot)
      phase=(-1)**(n1+n2+nr+lambda)
      dfactor=sq1dra**lr
      oscb=phase/dfactor*0.5d0*(lrlr+1.d0)*gamtot*sum
      end 


      double precision FUNCTION CLEBRD(A,B,C,D,E,F)
C
C      ARGUMENTS ARE REAL AND OF TRUE VALUE; J1,M1,J2,M2,J3,M3
C
      IMPLICIT double precision (A-H,O-Z)
      COMMON / LOGFAD / FIRST,LFACT,FACLOG(200)
      LOGICAL FIRST
*      REAL A,B,C,D,E,F
C      REAL*8 FACLOG
CCCCCC      IA=2*J1,ID=2*M1 ETC.(J1 IS OF TRUE VALUE)
      IA=idINT(2.d0*A)
      IB=idINT(2.d0*C)
      IC=idINT(2.d0*E)
      ID=idINT(2.d0*B)
      IE=idINT(2.d0*D)
      IF=idINT(2.d0*F)
      GOTO 7000
C ...............  CLEBI  ...........................
C
C      ARGUMENTS ARE INTEGER AND OF TRUE VALUE
C
      ENTRY CLEBID(LL1,LM1,LL2,LM2,LL3,LM3)
      IA=2*LL1
      IB=2*LL2
      IC=2*LL3
      ID=2*LM1
      IE=2*LM2
      IF=2*LM3
      GOTO 7000
C ..............  CLEB  ..............................
C
C      ARGUMENTS ARE INTEGER AND REPRESENT TWICE THE REAL VALUE
C
      ENTRY CLEBD(I2J1,I2M1,I2J2,I2M2,I2J3,I2M3)
      IA=I2J1
      IB=I2J2
      IC=I2J3
      ID=I2M1
      IE=I2M2
      IF=I2M3
 7000 IF (FIRST) CALL THFACD
      RAC=0.0D0
      IF(ID+IE-IF) 1000,105,1000
  105 K1=IA+IB+IC
      IF((-1)**K1) 1000,110,110
  110 K1=IA+IB-IC
      K2=IC+IA-IB
      K3=IB+IC-IA
      K4=IA-IABS (IB-IC)
      K5=IB-IABS (IC-IA)
      K6=IC-IABS (IA-IB)
      K7= MIN0 (K1,K2,K3,K4,K5,K6)
      IF(K7) 1000,120,120
  120 IF((-1)**(IA+ID)) 1000,1000,130
  130 IF((-1)**(IB+IE)) 1000,1000,140
  140 IF((-1)**(IC+IF)) 1000,1000,150
  150 IF(IA-IABS (ID)) 1000,152,152
  152 IF(IB-IABS (IE)) 1000,154,154
  154 IF(IC-IABS (IF)) 1000,160,160
  160 SIGNFC=1.0D0
      IAM=IA
      IBM=IB
      ICM=IC
      IDM=ID
      IEM=IE
      IFM=IF
      IF(IA-IB) 210,220,220
  210 IF(IA-IC) 215,225,225
  215 IT=IA
      IA=IB
      IB=IT
      IT=ID
      ID=IE
      IE=IT
      SIGNFC=(-1.0D0)**((IA+IB-IC)/2)
      GO TO 235
  220 IF(IC-IB) 225,235,235
  225 IT=IC
      IC=IB
      IB=IT
      IT=IF
      IF=-IE
      IE=-IT
      FIBM=IBM+1
      FICM=ICM+1
      SIGNFC=(-1.D0)**((IAM-IDM)/2)*DSQRT (FICM/FIBM)
  235 IF(IB) 237,236,237
  236 RAC=SIGNFC
      GO TO 900
  237 IF(IE) 250,250,240
  240 SIGNFC=SIGNFC*((-1.D0)**((IA+IB-IC)/2))
      ID=-ID
      IE=-IE
      IF=-IF
  250 FC2=IC+1
      IABCP=(IA+IB+IC)/2+1
      IABC=IABCP-IC
      ICAB=IABCP-IB
      IBCA=IABCP-IA
      IAPD=(IA+ID)/2+1
      IAMD=IAPD-ID
      IBPE=(IB+IE)/2+1
      IBME=IBPE-IE
      ICPF=(IC+IF)/2+1
      ICMF=ICPF-IF
      SQFCLG=0.5D0*(DLOG(FC2)-FACLOG(IABCP+1)
     1      +FACLOG(IABC)+FACLOG(ICAB)+FACLOG(IBCA)
     2      +FACLOG(IAPD)+FACLOG(IAMD)+FACLOG(IBPE)
     3      +FACLOG(IBME)+FACLOG(ICPF)+FACLOG(ICMF))
      NZMIC2=(IB-IC-ID)/2
      NZMIC3=(IA-IC+IE)/2
      NZMI= MAX0 (0,NZMIC2,NZMIC3)+1
      NZMX= MIN0 (IABC,IAMD,IBPE)
      IF(NZMI-NZMX) 310,310,900
  310 SS=0.D0
      S1=(-1.D0)**(NZMI-1)
      DO 400 NZ=NZMI,NZMX
      NZM1=NZ-1
      NZT1=IABC-NZM1
      NZT2=IAMD-NZM1
      NZT3=IBPE-NZM1
      NZT4=NZ-NZMIC2
      NZT5=NZ-NZMIC3
      TERMLG=SQFCLG-FACLOG(NZ)-FACLOG(NZT1)-FACLOG(NZT2)
     1           -FACLOG(NZT3)-FACLOG(NZT4)-FACLOG(NZT5)
      SSTERM=S1*DEXP (TERMLG)
      SS=SS+SSTERM
  400 S1=-S1
      RAC=SIGNFC*SS
  900 IA=IAM
      IB=IBM
      IC=ICM
      ID=IDM
      IE=IEM
      IF=IFM
 1000 CLEBRD=RAC
      RETURN
      END


      SUBROUTINE THFACD
C ...  SET UP LOG OF FACTORIALS
      PARAMETER (LFACTC=200)
      LOGICAL FIRST
      double precision FACLOG,FN
      COMMON / LOGFAD / FIRST,LFACT,FACLOG(200)
      FIRST=.FALSE.
      LFACT=200
      FACLOG(1)=0.D0
      FACLOG(2)=0.D0
      FN=1.D0
      DO 10 I=3,LFACT
      FN=FN+1.D0
      FACLOG(I)=FACLOG(I-1)+DLOG(FN)
   10 CONTINUE
      RETURN
      END

      block data
      PARAMETER (LFACTC=200)
      LOGICAL FIRST
      double precision FACLOG
      COMMON / LOGFAD / FIRST,LFACT,FACLOG(200)
      DATA FIRST/.TRUE./,LFACT/LFACTC/
      end

      double precision FUNCTION RACAD(JAD,JBD,JCD,JDD,JED,JFD)
C
C        CALCULATES RACAH COEFFICIENTS
C
C        ORIGINAL SOURCE : UNKNOWN
C        SOURCE : IBA_PROGRAM LIBRARY
C        MODIFIED : MARCH 1982 , OLAF SCHOLTEN
C              RUN TIME OPTIMIZED FOR VAX780 MACHINE
C
C	 ENTRIES : RACAD , RACADI , RACADR
C            RACAH  : INTEGER ARGUMENTS = 2*J
C            RACAHI : ARGUMENTS = TRUE INTEGER VALUE
C            RACAHR : ARGUMENTS = TRUE REAL VALUE
C	 EXTERNAL : THFACD , GENERATES FACTORIAL TABLE
C
      IMPLICIT double precision(A-H,O-Z)
      DIMENSION I(16)
      LOGICAL FIRST
*      double precision G,S,PH,H
*      double precision A,B,C,D,E,F
      COMMON / LOGFAD / FIRST,LFACT,G(200)
      EQUIVALENCE(I(1),I1),(I(2),I2),(I(3),I3),(I(4),I4),(I(5),I5),
     1 (I(6),I6),(I(7),I7),(I(8),I8),(I(9),I9),(I(10),I10),(I(11),I11),
     2 (I(12),I12),(I(13),I13),(I(14),I14),(I(15),I15),(I(16),I16)
C        MAKE USEFULL COMBINATIONS
      K=JAD+JBD-JED+2
      I1=K/2
      IF((2*I1).NE.K) GOTO 300
      K=JCD+JDD-JED+2
      I4=K/2
      IF((2*I4).NE.K) GOTO 300
      K=JAD+JCD-JFD+2
      I7=K/2
      IF((2*I7).NE.K) GOTO 300
      K=JBD+JDD-JFD+2
      I10=K/2
      IF((2*I10).NE.K) GOTO 300
      I13=I1+JED
      I14=I4+JED
      I15=I7+JFD
      I16=I10+JFD
      I2=I13-JAD
      I3=I13-JBD
      I5=I14-JCD
      I6=I14-JDD
      I8=I15-JAD
      I9=I15-JCD
      I11=I16-JBD
      I12=I16-JDD
C       CHECK TRIANGULAR INEQUALITIES,FIND NO. OF TERMS IN SUM
      N=MIN(I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12)-1
      IF(N) 300,2,2
C       FIND MINIMUM VALUE OF SUMMATION INDEX
    2 IL=MAX(I13,I14,I15,I16)
      IF(MIN(JAD,JBD,JCD,JDD,JED,JFD)) 300,20,1
C    ..............
      ENTRY RACADI(JA1,JB1,JC1,JD1,JE1,JF1)
C        MAKE USEFULL COMBINATIONS
      I13=JA1+JB1+JE1+1
      I14=JC1+JD1+JE1+1
      I15=JA1+JC1+JF1+1
      I16=JB1+JD1+JF1+1
      I1=I13-JE1*2
      I2=I13-JA1*2
      I3=I13-JB1*2
      I4=I14-JE1*2
      I5=I14-JC1*2
      I6=I14-JD1*2
      I7=I15-JF1*2
      I8=I15-JA1*2
      I9=I15-JC1*2
      I10=I16-JF1*2
      I11=I16-JB1*2
      I12=I16-JD1*2
C       CHECK TRIANGULAR INEQUALITIES,FIND NO. OF TERMS IN SUM
      N=MIN(I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12)-1
      IF(N) 300,4,4
C       FIND MINIMUM VALUE OF SUMMATION INDEX
    4 IL=MAX(I13,I14,I15,I16)
      LMIN=MIN(JA1,JB1,JC1,JD1,JE1,JF1)
      IF(LMIN)300,20,1
C     ............
      ENTRY RACADR(A,B,C,D,E,F)
C     CONVERT ARGUMENTS TO INTEGER
      JA=idINT(2.d0*A)
      JB=idINT(2.d0*B)
      JC=idINT(2.d0*C)
      JD=idINT(2.d0*D)
      JE=idINT(2.d0*E)
      JF=idINT(2.d0*F)
C        MAKE USEFULL COMBINATIONS
      K=JA+JB-JE+2
      I1=K/2
      IF((2*I1-K).NE.0) GOTO 300
      K=JC+JD-JE+2
      I4=K/2
      IF((2*I4-K).NE.0) GOTO 300
      K=JA+JC-JF+2
      I7=K/2
      IF((2*I7-K).NE.0) GOTO 300
      K=JB+JD-JF+2
      I10=K/2
      IF((2*I10-K).NE.0) GOTO 300
      I13=I1+JE
      I14=I4+JE
      I15=I7+JF
      I16=I10+JF
      I2=I13-JA
      I3=I13-JB
      I5=I14-JC
      I6=I14-JD
      I8=I15-JA
      I9=I15-JC
      I11=I16-JB
      I12=I16-JD
C       CHECK TRIANGULAR INEQUALITIES,FIND NO. OF TERMS IN SUM
      N=MIN(I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12)-1
      IF(N) 300,3,3
C       FIND MINIMUM VALUE OF SUMMATION INDEX
    3 IL=MAX(I13,I14,I15,I16)
      LMIN=MIN(JA,JB,JC,JD,JE,JF)
      IF(LMIN)300,20,1
C      ------------
    1 IF(FIRST) CALL THFACD
      IF(IL.GE.LFACT) STOP 'RACAH: LENGTH FACTORIAL TABLE INSUFFICIENT'
      J1=IL-I13+1
      J2=IL-I14+1
      J3=IL-I15+1
      J4=IL-I16+1
      J5=I13+I4-IL
      J6=I15+I5-IL
      J7=I16+I6-IL
      PH=1.D0
      IF(2*(J5/2).EQ.J5) PH=-1.D0
      H=PH*DEXP ((G(I1)+G(I2)+G(I3)-G(I13+1)+G(I4)+G(I5)+G(I6)-
     1G(I14+1)+G(I7)+G(I8)+G(I9)-G(I15+1)+G(I10)+G(I11)+G(I12)-G(I16+1))
     2*.5D0+G(IL+1)-G(J1)-G(J2)-G(J3)-G(J4)-G(J5)-G(J6)-G(J7))
      IF(N)300,110,120
C
  110 RACAD=H
      RETURN
C
  120 S=1.D0
      K=N-1
      KL=IL+1
      J5=J5-1
      J6=J6-1
      J7=J7-1
      DO 130 J=1,N 
C  ! K=N-J
      S=1.D0-((KL+K)*(J5-K)*(J6-K)*(J7-K))*S/((J1+K)*(J2+K)*(J3+K)
     + *(J4+K))
      K=K-1
  130 CONTINUE
      RACAD=H*S
      RETURN
C
C      ONE OF THE ARGUMENTS =0
   20 IAD=IL
      IBD=IL
      DO 21 J=13,16
      IF(IAD.LT.I(J)) GOTO 22
      IF(IAD.LT.IBD) IBD=IAD
      IAD=I(J)
      GOTO 21
   22 IF(IBD.GT.I(J)) IBD=I(J)
   21 CONTINUE
      J5=I13+I4-IL
      PH=1.D0
      IF(2*(J5/2).EQ.J5) PH=-1.D0
      RACAD=PH/DSQRT(DFLOAT(IAD*IBD))
      RETURN
C
C      IMPOSSIBLE COMBINATION OF ARGUMENTS
  300 RACAD=0.D0
      RETURN
      END


      double precision FUNCTION COEF9D(J1,J2,J3,J4,J5,J6,J7,J8,J9)
      IMPLICIT double precision (A-H,O-Z) 
C
CCCC  TAMURAS NUMBERING CONVENTION IS USED HERE.
CCCCC THE ARGUMENTS OF COEF9J, WHEN NUMBERED SEQUENTIALLY 1 THROUGH 9,
CCCC    CORRESPOND TO THE ARRAY
CCCC                               1  2  5
CCCC                               3  4  6
CCCC                               7  8  9
C
      DIMENSION LT(9)
CCCCCC
CCCCCC      ALL L9 MUST BE TWICE AS LARGE AS TRUE ARGUMENTS
CCCCCC
C
C               CHANGED FOR THEORY LIBRARY 4/8/82   HK
C
      U9=0.D0
      LT(1)=J1
      LT(2)=J2
      LT(3)=J3
      LT(4)=J4
      LT(5)=J5
      LT(6)=J6
      LT(7)=J7
      LT(8)=J8
      LT(9)=J9
      LMIN=LT(1)
      IMIN=1
      DO 20 I=2,9
      IF(LT(I)-LMIN) 15,20,20
   15 LMIN=LT(I)
      IMIN=I
   20 CONTINUE
      KEX=0
      GO TO (110,110,110,110,150,150,170,170,190),IMIN
  110 MM=(IMIN-1)/2+1
      M1=MM+MM-1
      M2=M1+1
      M3=MM+4
      L1=LT(7)
      LT(7)=LT(M1)
      LT(M1)=L1
      L1=LT(8)
      LT(8)=LT(M2)
      LT(M2)=L1
      L1=LT(9)
      LT(9)=LT(M3)
      LT(M3)=L1
      IMIN=IMIN+(7-M1)
      GO TO 175
  150 KEX=1
      M1=7
      M2=8
      M3=IMIN+IMIN-9
      M4=M3+1
      GO TO 180
  170 KEX=1
  175 M1=5
      M2=6
      M3=IMIN-6
      M4=M3+2
  180 L1=LT(M1)
      L1=LT(M1)
      LT(M1)=LT(M3)
      LT(M3)=L1
      L1=LT(M2)
      LT(M2)=LT(M4)
      LT(M4)=L1
      L1=LT(9)
      LT(9)=LT(IMIN)
      LT(IMIN)=L1
  190 IF(LT(9)) 200,200,300
  200 IF(LT(5)-LT(6)) 1000,210,1000
  210 IF(LT(7)-LT(8)) 1000,220,1000
  220 RT=(LT(5)+1)*(LT(7)+1)
      K=(LT(5)+LT(7)-LT(1)-LT(4))/2
      RAC= RACAD(LT(1),LT(2),LT(3),LT(4),LT(5),LT(7))
      PH=1.D0
      IF (2*(K/2) .NE. K) PH=-1.D0
      U9=(RAC/DSQRT(RT))*PH
      GO TO 370
  300 K1=IABS(LT(2)-LT(7))
      K2=IABS(LT(3)-LT(5))
      K3=IABS(LT(4)-LT(9))
      NMIN=MAX0(K1,K2,K3)
      K1=LT(2)+LT(7)
      K2=LT(3)+LT(5)
      K3=LT(4)+LT(9)
      NMAX=MIN0(K1,K2,K3)
      IF (NMIN-NMAX) 320, 320, 1000
  320 DO 350 N=NMIN,NMAX,2
      W1=N+1
      RAC= RACAD(LT(2),LT(5),LT(7),LT(3),LT(1),N)
      IF (RAC) 321, 350, 321
  321 W1=W1*RAC
      RAC= RACAD(LT(2),LT(4),LT(7),LT(9),LT(8),N)
      IF (RAC) 322, 350, 322
  322 W1=W1*RAC
      RAC= RACAD(LT(3),LT(4),LT(5),LT(9),LT(6),N)
      IF (RAC) 323, 350, 323
  323 U9=U9+W1*RAC
  350 CONTINUE
  370 IF(KEX) 400,1000,400
  400 KP=0
      DO 410 I=1,9
  410 KP=KP+LT(I)
      K=KP/2
      PH=1.d0
      IF (2*(K/2) .NE. K) PH=-1.d0
      U9=U9*PH
 1000 COEF9D=U9
      RETURN
      END
c********************************************************************
      subroutine inversg(a,n)
      implicit none
      integer, intent(IN) :: n
      real(kind=kind(0.d0)) a(n,n)

      integer ipvt(n),info,job
      real(kind=kind(0.d0)) work(n),det(2)

      call dgefa(a,n,n,ipvt,info)
      job=01
      call dgedi(a,n,n,ipvt,det,work,job)

      end


      subroutine dgefa(a,lda,n,ipvt,info)
      integer lda,n,ipvt(1),info
      double precision a(lda,1)
c
c     dgefa factors a double precision matrix by gaussian elimination.
c
c     dgefa is usually called by dgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgesl or dgedi will divide by zero
c                     if called.  use  rcond  in dgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c
c     internal variables
c
      double precision t
      integer idamax,j,k,kp1,l,nm1
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0d0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end


      subroutine dgedi(a,lda,n,ipvt,det,work,job)
      integer lda,n,ipvt(1),job
      double precision a(lda,1),det(2),work(1)
c
c     dgedi computes the determinant and inverse of a matrix
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        work    double precision(n)
c                work vector.  contents destroyed.
c
c        job     integer
c                = 11   both determinant and inverse.
c                = 01   inverse only.
c                = 10   determinant only.
c
c     on return
c
c        a       inverse of original matrix if requested.
c                otherwise unchanged.
c
c        det     double precision(2)
c                determinant of original matrix if requested.
c                otherwise not referenced.
c                determinant = det(1) * 10.0**det(2)
c                with  1.0 .le. dabs(det(1)) .lt. 10.0
c                or  det(1) .eq. 0.0 .
c
c     error condition
c
c        a division by zero will occur if the input factor contains
c        a zero on the diagonal and the inverse is requested.
c        it will not occur if the subroutines are called correctly
c        and if dgeco has set rcond .gt. 0.0 or dgefa has set
c        info .eq. 0 .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,dswap
c     fortran dabs,mod
c
c     internal variables
c
      double precision t
      double precision ten
      integer i,j,k,kb,kp1,l,nm1
c
c
c     compute determinant
c
      if (job/10 .eq. 0) go to 70
         det(1) = 1.0d0
         det(2) = 0.0d0
         ten = 10.0d0
         do 50 i = 1, n
            if (ipvt(i) .ne. i) det(1) = -det(1)
            det(1) = a(i,i)*det(1)
c        ...exit
            if (det(1) .eq. 0.0d0) go to 60
   10       if (dabs(det(1)) .ge. 1.0d0) go to 20
               det(1) = ten*det(1)
               det(2) = det(2) - 1.0d0
            go to 10
   20       continue
   30       if (dabs(det(1)) .lt. ten) go to 40
               det(1) = det(1)/ten
               det(2) = det(2) + 1.0d0
            go to 30
   40       continue
   50    continue
   60    continue
   70 continue
c
c     compute inverse(u)
c
      if (mod(job,10) .eq. 0) go to 150
         do 100 k = 1, n
            a(k,k) = 1.0d0/a(k,k)
            t = -a(k,k)
            call dscal(k-1,t,a(1,k),1)
            kp1 = k + 1
            if (n .lt. kp1) go to 90
            do 80 j = kp1, n
               t = a(k,j)
               a(k,j) = 0.0d0
               call daxpy(k,t,a(1,k),1,a(1,j),1)
   80       continue
   90       continue
  100    continue
c
c        form inverse(u)*inverse(l)
c
         nm1 = n - 1
         if (nm1 .lt. 1) go to 140
         do 130 kb = 1, nm1
            k = n - kb
            kp1 = k + 1
            do 110 i = kp1, n
               work(i) = a(i,k)
               a(i,k) = 0.0d0
  110       continue
            do 120 j = kp1, n
               t = work(j)
               call daxpy(n,t,a(1,j),1,a(1,k),1)
  120       continue
            l = ipvt(k)
            if (l .ne. k) call dswap(n,a(1,k),1,a(1,l),1)
  130    continue
  140    continue
  150 continue
      return
      end

      subroutine daxpy(n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dy(*),da
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end

      subroutine  dscal(n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end

      integer function idamax(n,dx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dmax
      integer i,incx,ix,n
c
      idamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end


      subroutine  dswap (n,dx,incx,dy,incy)
c
c     interchanges two vectors.
c     uses unrolled loops for increments equal one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
c
c       clean-up loop
c
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
        dtemp = dx(i + 1)
        dx(i + 1) = dy(i + 1)
        dy(i + 1) = dtemp
        dtemp = dx(i + 2)
        dx(i + 2) = dy(i + 2)
        dy(i + 2) = dtemp
   50 continue
      return
      end
c***end of inversg *********************************************

c*** For a machine that does not support popcnt ***      
      integer function popcnt(integ)
      use paramdef
      implicit none
      integer(4), intent(IN) :: integ
      integer :: i
      popcnt=0
      do i=0,nbit1
         if (btest(integ,i)) popcnt=popcnt+1
      end do      
      end


      subroutine rs(nm,n,a,w,matz,z,fv1,fv2,ierr)
c
      integer n,nm,ierr,matz
      double precision a(nm,n),w(n),z(nm,*),fv1(n),fv2(n)
c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a real symmetric matrix.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrix  a.
c
c        a  contains the real symmetric matrix.
c
c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        w  contains the eigenvalues in ascending order.
c
c        z  contains the eigenvectors if matz is not zero.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for tqlrat
c           and tql2.  the normal completion code is zero.
c
c        fv1  and  fv2  are temporary storage arrays.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  tred1(nm,n,a,w,fv1,fv2)
*  tqlrat encounters catastrophic underflow on the Vax
      call  tqlrat(n,w,fv2,ierr)
*      call  tql1(n,w,fv1,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 call  tred2(nm,n,a,w,fv1,z)
      call  tql2(nm,n,w,fv1,z,ierr)
   50 return
      end


      subroutine rsm(nm,n,a,w,m,z,fwork,iwork,ierr)
c 
      integer n,nm,m,iwork(n),ierr
      integer k1,k2,k3,k4,k5,k6,k7
      double precision a(nm,n),w(n),z(nm,m),fwork(1)
c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find all of the eigenvalues and some of the eigenvectors
c     of a real symmetric matrix.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrix  a.
c
c        a  contains the real symmetric matrix.
c
c        m  the eigenvectors corresponding to the first m eigenvalues
c           are to be computed.
c           if m = 0 then no eigenvectors are computed.
c           if m = n then all of the eigenvectors are computed.
c
c     on output
c
c        w  contains all n eigenvalues in ascending order.
c
c        z  contains the orthonormal eigenvectors associated with
c           the first m eigenvalues.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for tqlrat,
c           imtqlv and tinvit.  the normal completion code is zero.
c
c        fwork  is a temporary storage array of dimension 8*n.
c
c        iwork  is an integer temporary storage array of dimension n.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 10 * n
      if (n .gt. nm .or. m .gt. nm) go to 50
      k1 = 1
      k2 = k1 + n
      k3 = k2 + n
      k4 = k3 + n
      k5 = k4 + n
      k6 = k5 + n
      k7 = k6 + n
      k8 = k7 + n
      if (m .gt. 0) go to 10
c     .......... find eigenvalues only ..........
      call  tred1(nm,n,a,w,fwork(k1),fwork(k2))
      call  tqlrat(n,w,fwork(k2),ierr)
      go to 50
c     .......... find all eigenvalues and m eigenvectors ..........
   10 call  tred1(nm,n,a,fwork(k1),fwork(k2),fwork(k3))
      call  imtqlv(n,fwork(k1),fwork(k2),fwork(k3),w,iwork,
     x             ierr,fwork(k4))
      call  tinvit(nm,n,fwork(k1),fwork(k2),fwork(k3),m,w,iwork,z,ierr,
     x             fwork(k4),fwork(k5),fwork(k6),fwork(k7),fwork(k8))
      call  trbak1(nm,n,a,fwork(k2),m,z)
   50 return
      end

      subroutine tred1(nm,n,a,d,e,e2)
c
      integer i,j,k,l,n,ii,nm,jp1
      double precision a(nm,n),d(n),e(n),e2(n)
      double precision f,g,h,scale
c
c     this subroutine is a translation of the algol procedure tred1,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a real symmetric matrix
c     to a symmetric tridiagonal matrix using
c     orthogonal similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        a contains the real symmetric input matrix.  only the
c          lower triangle of the matrix need be supplied.
c
c     on output
c
c        a contains information about the orthogonal trans-
c          formations used in the reduction in its strict lower
c          triangle.  the full upper triangle of a is unaltered.
c
c        d contains the diagonal elements of the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        e2 contains the squares of the corresponding elements of e.
c          e2 may coincide with e if the squares are not needed.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      do 100 i = 1, n
         d(i) = a(n,i)
         a(n,i) = a(i,i)
  100 continue
c     .......... for i=n step -1 until 1 do -- ..........
      do 300 ii = 1, n
         i = n + 1 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 1) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + dabs(d(k))
c
         if (scale .ne. 0.0d0) go to 140
c
         do 125 j = 1, l
            d(j) = a(l,j)
            a(l,j) = a(i,j)
            a(i,j) = 0.0d0
  125    continue
c
  130    e(i) = 0.0d0
         e2(i) = 0.0d0
         go to 300
c
  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue
c
         e2(i) = scale * scale * h
         f = d(l)
         g = -dsign(dsqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
         if (l .eq. 1) go to 285
c     .......... form a*u ..........
         do 170 j = 1, l
  170    e(j) = 0.0d0
c
         do 240 j = 1, l
            f = d(j)
            g = e(j) + a(j,j) * f
            jp1 = j + 1
            if (l .lt. jp1) go to 220
c
            do 200 k = jp1, l
               g = g + a(k,j) * d(k)
               e(k) = e(k) + a(k,j) * f
  200       continue
c
  220       e(j) = g
  240    continue
c     .......... form p ..........
         f = 0.0d0
c
         do 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    continue
c
         h = f / (h + h)
c     .......... form q ..........
         do 250 j = 1, l
  250    e(j) = e(j) - h * d(j)
c     .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)
c
            do 260 k = j, l
  260       a(k,j) = a(k,j) - f * e(k) - g * d(k)
c
  280    continue
c
  285    do 290 j = 1, l
            f = d(j)
            d(j) = a(l,j)
            a(l,j) = a(i,j)
            a(i,j) = f * scale
  290    continue
c
  300 continue
c
      return
      end



      subroutine tql1(n,d,e,ierr)
c
      integer i,j,l,m,n,ii,l1,l2,mml,ierr
      double precision d(n),e(n)
      double precision c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
c
c     this subroutine is a translation of the algol procedure tql1,
c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
c     wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
c
c     this subroutine finds the eigenvalues of a symmetric
c     tridiagonal matrix by the ql method.
c
c     on input
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c      on output
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct and
c          ordered for indices 1,2,...ierr-1, but may not be
c          the smallest eigenvalues.
c
c        e has been destroyed.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e(i-1) = e(i)
c
      f = 0.0d0
      tst1 = 0.0d0
      e(n) = 0.0d0
c
      do 290 l = 1, n
         j = 0
         h = dabs(d(l)) + dabs(e(l))
         if (tst1 .lt. h) tst1 = h
c     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) go to 120
c     .......... e(n) is always zero, so there is no exit
c                through the bottom of the loop ..........
  110    continue
c
  120    if (m .eq. l) go to 210
  130    if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r = pythag(p,1.0d0)
         d(l) = e(l) / (p + dsign(r,p))
         d(l1) = e(l) * (p + dsign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
c
         do 140 i = l2, n
  140    d(i) = d(i) - h
c
  145    f = f + h
c     .......... ql transformation ..........
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
  200    continue
c
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + dabs(e(l))
         if (tst2 .gt. tst1) go to 130
  210    p = d(l) + f
c     .......... order eigenvalues ..........
         if (l .eq. 1) go to 250
c     .......... for i=l step -1 until 2 do -- ..........
         do 230 ii = 2, l
            i = l + 2 - ii
            if (p .ge. d(i-1)) go to 270
            d(i) = d(i-1)
  230    continue
c
  250    i = 1
  270    d(i) = p
  290 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end



      subroutine tred2(nm,n,a,d,e,z)
c
      integer i,j,k,l,n,ii,nm,jp1
      double precision a(nm,n),d(n),e(n),z(nm,*)
      double precision f,g,h,hh,scale
c
c     this subroutine is a translation of the algol procedure tred2,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a real symmetric matrix to a
c     symmetric tridiagonal matrix using and accumulating
c     orthogonal similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        a contains the real symmetric input matrix.  only the
c          lower triangle of the matrix need be supplied.
c
c     on output
c
c        d contains the diagonal elements of the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        z contains the orthogonal transformation matrix
c          produced in the reduction.
c
c        a and z may coincide.  if distinct, a is unaltered.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      do 100 i = 1, n
c
         do 80 j = i, n
   80    z(j,i) = a(j,i)
c
         d(i) = a(n,i)
  100 continue
c
      if (n .eq. 1) go to 510
c     .......... for i=n step -1 until 2 do -- ..........
      do 300 ii = 2, n
         i = n + 2 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 2) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + dabs(d(k))
c
         if (scale .ne. 0.0d0) go to 140
  130    e(i) = d(l)
c
         do 135 j = 1, l
            d(j) = z(l,j)
            z(i,j) = 0.0d0
            z(j,i) = 0.0d0
  135    continue
c
         go to 290
c
  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue
c
         f = d(l)
         g = -dsign(dsqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
c     .......... form a*u ..........
         do 170 j = 1, l
  170    e(j) = 0.0d0
c
         do 240 j = 1, l
            f = d(j)
            z(j,i) = f
            g = e(j) + z(j,j) * f
            jp1 = j + 1
            if (l .lt. jp1) go to 220
c
            do 200 k = jp1, l
               g = g + z(k,j) * d(k)
               e(k) = e(k) + z(k,j) * f
  200       continue
c
  220       e(j) = g
  240    continue
c     .......... form p ..........
         f = 0.0d0
c
         do 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    continue
c
         hh = f / (h + h)
c     .......... form q ..........
         do 250 j = 1, l
  250    e(j) = e(j) - hh * d(j)
c     .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)
c
            do 260 k = j, l
  260       z(k,j) = z(k,j) - f * e(k) - g * d(k)
c
            d(j) = z(l,j)
            z(i,j) = 0.0d0
  280    continue
c
  290    d(i) = h
  300 continue
c     .......... accumulation of transformation matrices ..........
      do 500 i = 2, n
         l = i - 1
         z(n,l) = z(l,l)
         z(l,l) = 1.0d0
         h = d(i)
         if (h .eq. 0.0d0) go to 380
c
         do 330 k = 1, l
  330    d(k) = z(k,i) / h
c
         do 360 j = 1, l
            g = 0.0d0
c
            do 340 k = 1, l
  340       g = g + z(k,i) * z(k,j)
c
            do 360 k = 1, l
               z(k,j) = z(k,j) - g * d(k)
  360    continue
c
  380    do 400 k = 1, l
  400    z(k,i) = 0.0d0
c
  500 continue
c
  510 do 520 i = 1, n
         d(i) = z(n,i)
         z(n,i) = 0.0d0
  520 continue
c
      z(n,n) = 1.0d0
      e(1) = 0.0d0
      return
      end


      subroutine tql2(nm,n,d,e,z,ierr)
c
      integer i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr
      double precision d(n),e(n),z(nm,n)
      double precision c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
c
c     this subroutine is a translation of the algol procedure tql2,
c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
c     wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a symmetric tridiagonal matrix by the ql method.
c     the eigenvectors of a full symmetric matrix can also
c     be found if  tred2  has been used to reduce this
c     full matrix to tridiagonal form.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c        z contains the transformation matrix produced in the
c          reduction by  tred2, if performed.  if the eigenvectors
c          of the tridiagonal matrix are desired, z must contain
c          the identity matrix.
c
c      on output
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct but
c          unordered for indices 1,2,...,ierr-1.
c
c        e has been destroyed.
c
c        z contains orthonormal eigenvectors of the symmetric
c          tridiagonal (or full) matrix.  if an error exit is made,
c          z contains the eigenvectors associated with the stored
c          eigenvalues.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e(i-1) = e(i)
c
      f = 0.0d0
      tst1 = 0.0d0
      e(n) = 0.0d0
c
      do 240 l = 1, n
         j = 0
         h = dabs(d(l)) + dabs(e(l))
         if (tst1 .lt. h) tst1 = h
c     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) go to 120
c     .......... e(n) is always zero, so there is no exit
c                through the bottom of the loop ..........
  110    continue
c
  120    if (m .eq. l) go to 220
  130    if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r = pythag(p,1.0d0)
         d(l) = e(l) / (p + dsign(r,p))
         d(l1) = e(l) * (p + dsign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
c
         do 140 i = l2, n
  140    d(i) = d(i) - h
c
  145    f = f + h
c     .......... ql transformation ..........
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
c     .......... form vector ..........
            do 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
  180       continue
c
  200    continue
c
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + dabs(e(l))
         if (tst2 .gt. tst1) go to 130
  220    d(l) = d(l) + f
  240 continue
c     .......... order eigenvalues and eigenvectors ..........
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
c
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
c
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
c
         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue
c
  300 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end


      double precision function pythag(a,b)
      double precision a,b
c
c     finds dsqrt(a**2+b**2) without overflow or destructive underflow
c
      double precision p,r,s,t,u
      p = dmax1(dabs(a),dabs(b))
      if (p .eq. 0.0d0) go to 20
      r = (dmin1(dabs(a),dabs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end


      double precision function epslon (x)
      double precision x
c
c     estimate unit roundoff in quantities of size x.
c
      double precision a,b,c,eps
c
c     this program should function properly on all systems
c     satisfying the following two assumptions,
c        1.  the base used in representing floating point
c            numbers is not a power of three.
c        2.  the quantity  a  in statement 10 is represented to 
c            the accuracy used in floating point variables
c            that are stored in memory.
c     the statement number 10 and the go to 10 are intended to
c     force optimizing compilers to generate code satisfying 
c     assumption 2.
c     under these assumptions, it should be true that,
c            a  is not exactly equal to four-thirds,
c            b  has a zero for its last bit or digit,
c            c  is not exactly equal to one,
c            eps  measures the separation of 1.0 from
c                 the next larger floating point number.
c     the developers of eispack would appreciate being informed
c     about any systems where these assumptions do not hold.
c
c     this version dated 4/6/83.
c
      a = 4.0d0/3.0d0
   10 b = a - 1.0d0
      c = b + b + b
      eps = dabs(c-1.0d0)
      if (eps .eq. 0.0d0) go to 10
      epslon = eps*dabs(x)
      return
      end



**** for old version, "send otqlrat from eispack"
** From dana!moler Tue, 1 Sep 87 10:15:40 PDT
** New TQLRAT
      SUBROUTINE TQLRAT(N,D,E2,IERR)
C
      INTEGER I,J,L,M,N,II,L1,MML,IERR
      DOUBLE PRECISION D(N),E2(N)
      DOUBLE PRECISION B,C,F,G,H,P,R,S,T,EPSLON,PYTHAG
C
C     This subroutine is a translation of the Algol procedure tqlrat,
C     Algorithm 464, Comm. ACM 16, 689(1973) by Reinsch.
C
C     This subroutine finds the eigenvalues of a symmetric
C     tridiagonal matrix by the rational QL method.
C
C     On input
C
C        N is the order of the matrix.
C
C        D contains the diagonal elements of the input matrix.
C
C        E2 contains the squares of the subdiagonal elements of the
C          input matrix in its last N-1 positions.  E2(1) is arbitrary.
C
C      On output
C
C        D contains the eigenvalues in ascending order.  If an
C          error exit is made, the eigenvalues are correct and
C          ordered for indices 1,2,...IERR-1, but may not be
C          the smallest eigenvalues.
C
C        E2 has been destroyed.
C
C        IERR is set to
C          zero       for normal return,
C          J          if the J-th eigenvalue has not been
C                     determined after 30 iterations.
C
C     Calls PYTHAG for  DSQRT(A*A + B*B) .
C
C     Questions and comments should be directed to Burton S. Garbow,
C     Mathematics and Computer Science Div, Argonne National Laboratory
C
C     This version dated August 1987.
C     Modified by C. Moler to fix underflow/overflow difficulties,
C     especially on the VAX and other machines where epslon(1.0d0)**2
C     nearly underflows.  See the loop involving statement 102 and
C     the two statements just before statement 200.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E2(I-1) = E2(I)
C
      F = 0.0D0
      T = 0.0D0
      E2(N) = 0.0D0
C
      DO 290 L = 1, N
         J = 0
         H = DABS(D(L)) + DSQRT(E2(L))
         IF (T .GT. H) GO TO 105
         T = H
         B = EPSLON(T)
         C = B * B
         if (c .ne. 0.0d0) go to 105
C        Spliting tolerance underflowed.  Look for larger value.
         do 102 i = l, n
            h = dabs(d(i)) + dsqrt(e2(i))
            if (h .gt. t) t = h
  102    continue
         b = epslon(t)
         c = b * b
C     .......... LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT ..........
  105    DO 110 M = L, N
            IF (E2(M) .LE. C) GO TO 120
C     .......... E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
C
  120    IF (M .EQ. L) GO TO 210
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     .......... FORM SHIFT ..........
         L1 = L + 1
         S = DSQRT(E2(L))
         G = D(L)
         P = (D(L1) - G) / (2.0D0 * S)
         R = PYTHAG(P,1.0D0)
         D(L) = S / (P + DSIGN(R,P))
         H = G - D(L)
C
         DO 140 I = L1, N
  140    D(I) = D(I) - H
C
         F = F + H
C     .......... RATIONAL QL TRANSFORMATION ..........
         G = D(M)
         IF (G .EQ. 0.0D0) G = B
         H = G
         S = 0.0D0
         MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            I = M - II
            P = G * H
            R = P + E2(I)
            E2(I+1) = S * R
            S = E2(I) / R
            D(I+1) = H + S * (H + D(I))
            G = D(I) - E2(I) / G
C           Avoid division by zero on next pass
            if (g .eq. 0.0d0) g = epslon(d(i))
            h = g * (p / r)
  200    CONTINUE
C
         E2(L) = S * G
         D(L) = H
C     .......... GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST ..........
         IF (H .EQ. 0.0D0) GO TO 210
         IF (DABS(E2(L)) .LE. DABS(C/H)) GO TO 210
         E2(L) = H * E2(L)
         IF (E2(L) .NE. 0.0D0) GO TO 130
  210    P = D(L) + F
C     .......... ORDER EIGENVALUES ..........
         IF (L .EQ. 1) GO TO 250
C     .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........
         DO 230 II = 2, L
            I = L + 2 - II
            IF (P .GE. D(I-1)) GO TO 270
            D(I) = D(I-1)
  230    CONTINUE
C
  250    I = 1
  270    D(I) = P
  290 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END



      subroutine imtqlv(n,d,e,e2,w,ind,ierr,rv1)
c
      integer i,j,k,l,m,n,ii,mml,tag,ierr
      double precision d(n),e(n),e2(n),w(n),rv1(n)
      double precision b,c,f,g,p,r,s,tst1,tst2,pythag
      integer ind(n)
c
c     this subroutine is a variant of  imtql1  which is a translation of
c     algol procedure imtql1, num. math. 12, 377-383(1968) by martin and
c     wilkinson, as modified in num. math. 15, 450(1970) by dubrulle.
c     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
c
c     this subroutine finds the eigenvalues of a symmetric tridiagonal
c     matrix by the implicit ql method and associates with them
c     their corresponding submatrix indices.
c
c     on input
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c        e2 contains the squares of the corresponding elements of e.
c          e2(1) is arbitrary.
c
c     on output
c
c        d and e are unaltered.
c
c        elements of e2, corresponding to elements of e regarded
c          as negligible, have been replaced by zero causing the
c          matrix to split into a direct sum of submatrices.
c          e2(1) is also set to zero.
c
c        w contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct and
c          ordered for indices 1,2,...ierr-1, but may not be
c          the smallest eigenvalues.
c
c        ind contains the submatrix indices associated with the
c          corresponding eigenvalues in w -- 1 for eigenvalues
c          belonging to the first submatrix from the top,
c          2 for those belonging to the second submatrix, etc..
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c        rv1 is a temporary storage array.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      k = 0
      tag = 0
c
      do 100 i = 1, n
         w(i) = d(i)
         if (i .ne. 1) rv1(i-1) = e(i)
  100 continue
c
      e2(1) = 0.0d0
      rv1(n) = 0.0d0
c
      do 290 l = 1, n
         j = 0
c     .......... look for small sub-diagonal element ..........
  105    do 110 m = l, n
            if (m .eq. n) go to 120
            tst1 = dabs(w(m)) + dabs(w(m+1))
            tst2 = tst1 + dabs(rv1(m))
            if (tst2 .eq. tst1) go to 120
c     .......... guard against underflowed element of e2 ..........
            if (e2(m+1) .eq. 0.0d0) go to 125
  110    continue
c
  120    if (m .le. k) go to 130
         if (m .ne. n) e2(m+1) = 0.0d0
  125    k = m
         tag = tag + 1
  130    p = w(l)
         if (m .eq. l) go to 215
         if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         g = (w(l+1) - p) / (2.0d0 * rv1(l))
         r = pythag(g,1.0d0)
         g = w(m) - p + rv1(l) / (g + dsign(r,g))
         s = 1.0d0
         c = 1.0d0
         p = 0.0d0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            i = m - ii
            f = s * rv1(i)
            b = c * rv1(i)
            r = pythag(f,g)
            rv1(i+1) = r
            if (r .eq. 0.0d0) go to 210
            s = f / r
            c = g / r
            g = w(i+1) - p
            r = (w(i) - g) * s + 2.0d0 * c * b
            p = s * r
            w(i+1) = g + p
            g = c * r - b
  200    continue
c
         w(l) = w(l) - p
         rv1(l) = g
         rv1(m) = 0.0d0
         go to 105
c     .......... recover from underflow ..........
  210    w(i+1) = w(i+1) - p
         rv1(m) = 0.0d0
         go to 105
c     .......... order eigenvalues ..........
  215    if (l .eq. 1) go to 250
c     .......... for i=l step -1 until 2 do -- ..........
         do 230 ii = 2, l
            i = l + 2 - ii
            if (p .ge. w(i-1)) go to 270
            w(i) = w(i-1)
            ind(i) = ind(i-1)
  230    continue
c
  250    i = 1
  270    w(i) = p
         ind(i) = tag
  290 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end


      subroutine tinvit(nm,n,d,e,e2,m,w,ind,z,
     x                  ierr,rv1,rv2,rv3,rv4,rv6)
c
      integer i,j,m,n,p,q,r,s,ii,ip,jj,nm,its,tag,ierr,group
      double precision d(n),e(n),e2(n),w(m),z(nm,m),
     x       rv1(n),rv2(n),rv3(n),rv4(n),rv6(n)
      double precision u,v,uk,xu,x0,x1,eps2,eps3,eps4,norm,order,epslon,
     x       pythag
      integer ind(m)
c
c     this subroutine is a translation of the inverse iteration tech-
c     nique in the algol procedure tristurm by peters and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 418-439(1971).
c
c     this subroutine finds those eigenvectors of a tridiagonal
c     symmetric matrix corresponding to specified eigenvalues,
c     using inverse iteration.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c        e2 contains the squares of the corresponding elements of e,
c          with zeros corresponding to negligible elements of e.
c          e(i) is considered negligible if it is not larger than
c          the product of the relative machine precision and the sum
c          of the magnitudes of d(i) and d(i-1).  e2(1) must contain
c          0.0d0 if the eigenvalues are in ascending order, or 2.0d0
c          if the eigenvalues are in descending order.  if  bisect,
c          tridib, or  imtqlv  has been used to find the eigenvalues,
c          their output e2 array is exactly what is expected here.
c
c        m is the number of specified eigenvalues.
c
c        w contains the m eigenvalues in ascending or descending order.
c
c        ind contains in its first m positions the submatrix indices
c          associated with the corresponding eigenvalues in w --
c          1 for eigenvalues belonging to the first submatrix from
c          the top, 2 for those belonging to the second submatrix, etc.
c
c     on output
c
c        all input arrays are unaltered.
c
c        z contains the associated set of orthonormal eigenvectors.
c          any vector which fails to converge is set to zero.
c
c        ierr is set to
c          zero       for normal return,
c          -r         if the eigenvector corresponding to the r-th
c                     eigenvalue fails to converge in 5 iterations.
c
c        rv1, rv2, rv3, rv4, and rv6 are temporary storage arrays.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (m .eq. 0) go to 1001
      tag = 0
      order = 1.0d0 - e2(1)
      q = 0
c     .......... establish and process next submatrix ..........
  100 p = q + 1
c
      do 120 q = p, n
         if (q .eq. n) go to 140
         if (e2(q+1) .eq. 0.0d0) go to 140
  120 continue
c     .......... find vectors by inverse iteration ..........
  140 tag = tag + 1
      s = 0
c
      do 920 r = 1, m
         if (ind(r) .ne. tag) go to 920
         its = 1
         x1 = w(r)
         if (s .ne. 0) go to 510
c     .......... check for isolated root ..........
         xu = 1.0d0
         if (p .ne. q) go to 490
         rv6(p) = 1.0d0
         go to 870
  490    norm = dabs(d(p))
         ip = p + 1
c
         do 500 i = ip, q
  500    norm = dmax1(norm, dabs(d(i))+dabs(e(i)))
c     .......... eps2 is the criterion for grouping,
c                eps3 replaces zero pivots and equal
c                roots are modified by eps3,
c                eps4 is taken very small to avoid overflow ..........
         eps2 = 1.0d-3 * norm
         eps3 = epslon(norm)
         uk = q - p + 1
         eps4 = uk * eps3
         uk = eps4 / dsqrt(uk)
         s = p
  505    group = 0
         go to 520
c     .......... look for close or coincident roots ..........
  510    if (dabs(x1-x0) .ge. eps2) go to 505
         group = group + 1
         if (order * (x1 - x0) .le. 0.0d0) x1 = x0 + order * eps3
c     .......... elimination with interchanges and
c                initialization of vector ..........
  520    v = 0.0d0
c
         do 580 i = p, q
            rv6(i) = uk
            if (i .eq. p) go to 560
            if (dabs(e(i)) .lt. dabs(u)) go to 540
c     .......... warning -- a divide check may occur here if
c                e2 array has not been specified correctly ..........
            xu = u / e(i)
            rv4(i) = xu
            rv1(i-1) = e(i)
            rv2(i-1) = d(i) - x1
            rv3(i-1) = 0.0d0
            if (i .ne. q) rv3(i-1) = e(i+1)
            u = v - xu * rv2(i-1)
            v = -xu * rv3(i-1)
            go to 580
  540       xu = e(i) / u
            rv4(i) = xu
            rv1(i-1) = u
            rv2(i-1) = v
            rv3(i-1) = 0.0d0
  560       u = d(i) - x1 - xu * v
            if (i .ne. q) v = e(i+1)
  580    continue
c
         if (u .eq. 0.0d0) u = eps3
         rv1(q) = u
         rv2(q) = 0.0d0
         rv3(q) = 0.0d0
c     .......... back substitution
c                for i=q step -1 until p do -- ..........
  600    do 620 ii = p, q
            i = p + q - ii
            rv6(i) = (rv6(i) - u * rv2(i) - v * rv3(i)) / rv1(i)
            v = u
            u = rv6(i)
  620    continue
c     .......... orthogonalize with respect to previous
c                members of group ..........
         if (group .eq. 0) go to 700
         j = r
c
         do 680 jj = 1, group
  630       j = j - 1
            if (ind(j) .ne. tag) go to 630
            xu = 0.0d0
c
            do 640 i = p, q
  640       xu = xu + rv6(i) * z(i,j)
c
            do 660 i = p, q
  660       rv6(i) = rv6(i) - xu * z(i,j)
c
  680    continue
c
  700    norm = 0.0d0
c
         do 720 i = p, q
  720    norm = norm + dabs(rv6(i))
c
         if (norm .ge. 1.0d0) go to 840
c     .......... forward substitution ..........
         if (its .eq. 5) go to 830
         if (norm .ne. 0.0d0) go to 740
         rv6(s) = eps4
         s = s + 1
         if (s .gt. q) s = p
         go to 780
  740    xu = eps4 / norm
c
         do 760 i = p, q
  760    rv6(i) = rv6(i) * xu
c     .......... elimination operations on next vector
c                iterate ..........
  780    do 820 i = ip, q
            u = rv6(i)
c     .......... if rv1(i-1) .eq. e(i), a row interchange
c                was performed earlier in the
c                triangularization process ..........
            if (rv1(i-1) .ne. e(i)) go to 800
            u = rv6(i-1)
            rv6(i-1) = rv6(i)
  800       rv6(i) = u - rv4(i) * rv6(i-1)
  820    continue
c
         its = its + 1
         go to 600
c     .......... set error -- non-converged eigenvector ..........
  830    ierr = -r
         xu = 0.0d0
         go to 870
c     .......... normalize so that sum of squares is
c                1 and expand to full order ..........
  840    u = 0.0d0
c
         do 860 i = p, q
  860    u = pythag(u,rv6(i))
c
         xu = 1.0d0 / u
c
  870    do 880 i = 1, n
  880    z(i,r) = 0.0d0
c
         do 900 i = p, q
  900    z(i,r) = rv6(i) * xu
c
         x0 = x1
  920 continue
c
      if (q .lt. n) go to 100
 1001 return
      end


      subroutine trbak1(nm,n,a,e,m,z)
c
      integer i,j,k,l,m,n,nm
      double precision a(nm,n),e(n),z(nm,m)
      double precision s
c
c     this subroutine is a translation of the algol procedure trbak1,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine forms the eigenvectors of a real symmetric
c     matrix by back transforming those of the corresponding
c     symmetric tridiagonal matrix determined by  tred1.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        a contains information about the orthogonal trans-
c          formations used in the reduction by  tred1
c          in its strict lower triangle.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is arbitrary.
c
c        m is the number of eigenvectors to be back transformed.
c
c        z contains the eigenvectors to be back transformed
c          in its first m columns.
c
c     on output
c
c        z contains the transformed eigenvectors
c          in its first m columns.
c
c     note that trbak1 preserves vector euclidean norms.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (m .eq. 0) go to 200
      if (n .eq. 1) go to 200
c
      do 140 i = 2, n
         l = i - 1
         if (e(i) .eq. 0.0d0) go to 140
c
         do 130 j = 1, m
            s = 0.0d0
c
            do 110 k = 1, l
  110       s = s + a(i,k) * z(k,j)
c     .......... divisor below is negative of h formed in tred1.
c                double division avoids possible underflow ..........
            s = (s / a(i,l)) / e(i)
c
            do 120 k = 1, l
  120       z(k,j) = z(k,j) + s * a(i,k)
c
  130    continue
c
  140 continue
c
  200 return
      end


c---------------------------------------------------------------------
c     *******  The following subroutines and functions are     *******
c     *******  to be used in the sequential environment.       *******
c---------------------------------------------------------------------


c      subroutine MPI_INIT(ierr)
c      return
c      end

c      subroutine MPI_COMM_RANK(icomm,iproc,ierr)
c      iproc=0
c      return
c      end

c      subroutine MPI_COMM_SIZE(icomm,nproc,ierr)
c      nproc=1
c      return
c      end

c      subroutine MPI_Barrier(icomm,ierr)
c      return
c      end


c      subroutine MPI_Finalize(ierr)
c      return
c      end

c      subroutine MPI_Bcast()
c      return
c      end

c      subroutine MPI_Reduce()
c      return
c      end

c      subroutine MPI_Allreduce()
c      return
c      end

c      subroutine MPI_Abort()
c      return
c      end



