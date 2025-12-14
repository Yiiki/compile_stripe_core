program stripe_core
implicit none

type rho_data
  integer*4 :: n1,n2,n3,nnodes_tmp
  real*8 :: AL(3,3)
  real*8,allocatable,dimension(:,:,:) :: rho
end type rho_data
type(rho_data) :: a,b,c
character*256 frho, fout
real*8,parameter :: A_AU=0.529177d0, E_EU=27.211386d0
real*8 :: AL(3,3)

integer*4 :: m1,m2,m3,i,j,k
real*8 :: va,vb,aa,bb
integer*4 :: ai,aj,ak,bi,bj,bk

call getarg(1,frho)
call read_rho(frho,a)

write(fout,191) trim(frho), ".dscore"
191 format(A,A)

call getarg(2,frho)
call read_rho(frho,b)

call check_match(a,b)
! we actually assume a and b are concentrated at (1/2 1/2 1/2)
! and share mesh resolution
! grid number be even
m1=min(a%n1,b%n1)
m2=min(a%n2,b%n2)
m3=min(a%n3,b%n3)

va=1.d0*a%n1*a%n2*a%n3
vb=1.d0*b%n1*b%n2*b%n3

! output c%rho on the common grid
!
!     va                 vb
!  ------- * a%rho  + ------- * b%rho
!  va - vb            vb - va

aa=va/(va-vb)
bb=vb/(vb-va)

c%n1=m1
c%n2=m2
c%n3=m3
c%nnodes_tmp=1
AL=0.d0
AL(1,1)=min(a%AL(1,1),b%AL(1,1))
AL(2,2)=min(a%AL(2,2),b%AL(2,2))
AL(3,3)=min(a%AL(3,3),b%AL(3,3))
c%AL=AL
allocate(c%rho(m1,m2,m3))
c%rho=0.d0
do k=1,m3
   ak=a%n3/2+k-m3/2
   bk=b%n3/2+k-m3/2
if(ak<1.or.ak>a%n3.or.&
   bk<1.or.bk>b%n3) cycle
do j=1,m2
   aj=a%n2/2+j-m2/2
   bj=b%n2/2+j-m2/2
if(aj<1.or.aj>a%n2.or.&
   bj<1.or.bj>b%n2) cycle
do i=1,m1
   ai=a%n1/2+i-m1/2
   bi=b%n1/2+i-m1/2
if(ai<1.or.ai>a%n1.or.&
   bi<1.or.bi>b%n1) cycle
c%rho(i,j,k)=aa*a%rho(ai,aj,ak)+bb*b%rho(bi,bj,bk)
end do
end do
end do

call write_rho(fout,c)
write(6,*) trim(fout), " has been written"

contains
subroutine check_match(a,b)
implicit none
type(rho_data) :: a,b
logical :: bad
bad=&
mod(a%n1,2).ne.0.or.&
mod(a%n2,2).ne.0.or.&
mod(a%n3,2).ne.0.or.&
mod(b%n1,2).ne.0.or.&
mod(b%n2,2).ne.0.or.&
mod(b%n3,2).ne.0
if(bad) then
  write(6,*) "even n123 both required for a and b"
  write(6,*) "a%n123"
  write(6,*) a%n1, a%n2, a%n3
  write(6,*) "b%n123"
  write(6,*) b%n1, b%n2, b%n3
  stop
end if
bad=&
norm2(a%AL(:,1)/a%n1-b%AL(:,1)/b%n1).gt.1.d-10.or.&
norm2(a%AL(:,2)/a%n2-b%AL(:,2)/b%n2).gt.1.d-10.or.&
norm2(a%AL(:,3)/a%n3-b%AL(:,3)/b%n3).gt.1.d-10
if(bad) then
  write(6,*) "mismatch AL resolution, stop"
  write(6,*) a%AL
  write(6,*) a%n1, a%n2, a%n3
  write(6,*) b%AL
  write(6,*) b%n1, b%n2, b%n3
  stop
end if
end subroutine check_match
subroutine check_file(frho)
implicit none
character(len=*) :: frho
logical :: iqst
inquire(file=trim(frho), exist=iqst)
if(.not.iqst) then
  write(6,*) "missing file: ", trim(frho)
  stop
end if
end subroutine check_file
subroutine read_rho(file,p)
implicit none
character(len=*) :: file
type(rho_data) :: p
integer*4 :: n1,n2,n3,nnodes_tmp,nr_n,iread,i,j,k,ii,jj
real*8 :: AL(3,3)
real*8,allocatable :: vr0(:,:,:),vr_tmp(:)
call check_file(file)
open(11,file=trim(adjustl(file)),form="unformatted")
rewind(11)
read(11) n1,n2,n3,nnodes_tmp
read(11) AL
nr_n=(n1*n2*n3)/nnodes_tmp
allocate(vr0(n1,n2,n3))
allocate(vr_tmp(nr_n))
do iread=1,nnodes_tmp
   read(11) vr_tmp
   if(iread.eq.1) then
   write(6,*) "read file1 ..."
   write(6,*) "vr_tmp first 9 "
   write(6,'(3(E14.7,1x))') vr_tmp(1:3)
   write(6,'(3(E14.7,1x))') vr_tmp(4:6)
   write(6,'(3(E14.7,1x))') vr_tmp(7:9)
  endif
  do ii=1,nr_n
      jj=ii+(iread-1)*nr_n
      i=(jj-1)/(n2*n3)+1
      j=(jj-1-(i-1)*n2*n3)/n3+1
      k=jj-(i-1)*n2*n3-(j-1)*n3
      vr0(i,j,k)=vr_tmp(ii)
   enddo
enddo
close(11)
deallocate(vr_tmp)
p%n1=n1
p%n2=n2
p%n3=n3
p%nnodes_tmp=nnodes_tmp
p%AL=AL
allocate(p%rho(n1,n2,n3))
p%rho=vr0
deallocate(vr0)
end subroutine read_rho
subroutine write_rho(file,p)
implicit none
character(len=*) :: file
type(rho_data) :: p
integer*4 :: n1,n2,n3,nnodes_tmp,nr_n,iread,i,j,k,ii,jj
real*8 :: AL(3,3)
real*8,allocatable :: vr0(:,:,:),vr_tmp(:)
nnodes_tmp=p%nnodes_tmp
n1=p%n1
n2=p%n2
n3=p%n3
AL=p%AL     
allocate(vr0(n1,n2,n3))
vr0=p%rho
open(13,file=trim(adjustl(file)),form="unformatted")
rewind(13)
write(13) n1,n2,n3,nnodes_tmp
write(13) AL
nr_n=(n1*n2*n3)/nnodes_tmp
allocate(vr_tmp(nr_n))
do iread=1,nnodes_tmp
   do ii=1,nr_n
      jj=ii+(iread-1)*nr_n
      i=(jj-1)/(n2*n3)+1
      j=(jj-1-(i-1)*n2*n3)/n3+1
      k=jj-(i-1)*n2*n3-(j-1)*n3
      vr_tmp(ii)=vr0(i,j,k)
   enddo
   if(iread.eq.1) then
   write(6,*) "vr_tmp first 9 "
   write(6,'(3(E14.7,1x))') vr_tmp(1:3)
   write(6,'(3(E14.7,1x))') vr_tmp(4:6)
   write(6,'(3(E14.7,1x))') vr_tmp(7:9)
   endif
   write(13) vr_tmp
enddo
close(13)
deallocate(vr0)
deallocate(vr_tmp)
end subroutine write_rho

end
