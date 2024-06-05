subroutine jspd(m, n, k, Adata, jb, Bdata, d, tv, na, nb, nd)
  implicit none
  integer, intent(in) :: m, n, k, na, nb, nd
  integer, intent(out) :: jb, d(nd)
  integer :: i, j, kend, q, mn, kk1, js, je, qx
  double precision, intent(inout) :: Adata(na), Bdata(nb), tv(k)
  logical :: not0

  jb = 0 
  kend=0
  q = (m-1)+n 
  mn = m*n
  Bdata(1:mn) = Adata(1:mn)
  Adata(1:(k-1)*k) = 0.0
  Adata((k-1)*k+1:(k-1)*k+mn) = Bdata(1:mn)
  Adata((k-1)*k+mn+1:(k-1)*k+mn+(k-1)*k) = 0.0

  do i=1,q
    qx = k*(i-1)+1
    tv(1) = Adata(qx)
    tv(2:k) = 0.0
    qx = min(q-i, k-1) 
    if (qx > 0) then
      do j=1,qx
        tv(j+1) = Adata(k*(i+j-1)+j+1)
      end do
    endif
    not0 = any(tv(1:k) /= 0.0)
    if (not0) then 
      jb = jb+1
      d(jb) = (i-k)
      Bdata(kend+1:kend+k) = tv(1:k)
      kend = kend + k
    endif
  end do
  return
end subroutine jspd