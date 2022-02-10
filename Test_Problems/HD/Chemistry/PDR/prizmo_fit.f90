module prizmo_fit
contains

   ! ********************
   subroutine load_data_1d(fname, fdata, xdata, nx, xmin, xfact, invdx, do_log)
      implicit none
      character(len=*),intent(in)::fname
      integer,intent(in)::nx
      integer::unit, i
      real*8::dx
      real*8,intent(out)::fdata(nx), xmin, xfact, invdx
      real*8,intent(out)::xdata(nx)
      logical,intent(in),optional::do_log
      logical::dolog

      dolog = .false.
      if(present(do_log)) then
         dolog = do_log
      end if

      open(newunit=unit, file=trim(fname), status="old")
      do i=1,nx
         read(unit, *) xdata(i), fdata(i)
      end do
      close(unit)

      if(dolog) then
         xdata(:) = log10(xdata(:))
         fdata(:) = log10(fdata(:))
      end if

      call check_regular_grid(xdata(:), fname)

      xmin = xdata(1)
      dx = xdata(2) - xdata(1)
      xfact = (nx - 1) / (xdata(nx) - xdata(1))
      invdx = 1d0 / dx

   end subroutine load_data_1d

   ! *******************
   ! given xdata return variables to be used as arguments for the 1d fit routine
   subroutine fit_prepare(xdata, nx, xmin, xfact, invdx)
      implicit none
      integer,intent(in)::nx
      real*8,intent(in)::xdata(nx)
      real*8::xmin, dx, xfact, invdx

      ! check if grid is uniform
      call check_regular_grid(xdata(:), "unknown")

      xmin = xdata(1)
      dx = xdata(2) - xdata(1)
      xfact = (nx - 1) / (xdata(nx) - xdata(1))
      invdx = 1d0 / dx

   end subroutine fit_prepare

   ! ********************
   subroutine load_data_2d(fname, fdata, &
      xdata, nx, xmin, xfact, invdx, &
      ydata, ny, ymin, yfact, invdy, mode, do_log)
      implicit none
      character(len=*),intent(in)::fname
      character(len=2),intent(in),optional::mode
      character(len=2)::mode_read
      integer,intent(in)::nx, ny
      integer::unit, i, j
      real*8::dx, dy
      real*8,intent(out)::fdata(nx, ny), xmin, xfact, invdx
      real*8,intent(out)::ymin, yfact, invdy
      real*8,intent(out)::xdata(nx), ydata(ny)
      logical,intent(in),optional::do_log
      logical::dolog

      mode_read = "xy"
      if(present(mode)) then
         mode_read = mode
      end if

      dolog = .false.
      if(present(do_log)) then
         dolog = do_log
      end if

      open(newunit=unit, file=trim(fname), status="old")
      if(mode_read=="xy") then
         do i=1,nx
            do j=1,ny
               read(unit, *) xdata(i), ydata(j), fdata(i, j)
            end do
         end do
      else if(mode_read=="yx") then
         do j=1,ny
            do i=1,nx
               read(unit, *) xdata(i), ydata(j), fdata(i, j)
            end do
         end do
      else
         print *, "ERROR: read mode in load interp_2d can be only xy or yx, found "//mode_read
         stop
      end if
      close(unit)

      if(dolog) then
         xdata(:) = log10(xdata(:))
         ydata(:) = log10(ydata(:))
         fdata(:, :) = log10(fdata(:, :))
      end if

      call check_regular_grid(xdata(:), fname)
      call check_regular_grid(ydata(:), fname)

      xmin = xdata(1)
      ymin = ydata(1)
      dx = xdata(2) - xdata(1)
      dy = ydata(2) - ydata(1)
      xfact = (nx - 1) / (xdata(nx) - xdata(1))
      yfact = (ny - 1) / (ydata(ny) - ydata(1))
      invdx = 1d0 / dx
      invdy = 1d0 / dy

   end subroutine load_data_2d

   ! ********************
   subroutine load_data_3d(fname, fdata, &
      xdata, nx, xmin, xfact, invdx, &
      ydata, ny, ymin, yfact, invdy, &
      zdata, nz, zmin, zfact, invdz, mode, do_log)
      implicit none
      character(len=*),intent(in)::fname
      integer,intent(in)::nx, ny, nz
      integer::unit, i, j, k
      real*8::dx, dy, dz
      real*8,intent(out)::fdata(nx, ny, nz), xdata(nx), xmin, xfact, invdx
      real*8,intent(out)::ydata(ny), ymin, yfact, invdy
      real*8,intent(out)::zdata(nz), zmin, zfact, invdz
      character(len=3),intent(in),optional::mode
      character(len=3)::mode_read
      logical,intent(in),optional::do_log
      logical::dolog

      mode_read = "xyz"
      if(present(mode)) then
         mode_read = mode
      end if

      dolog = .false.
      if(present(do_log)) then
         dolog = do_log
      end if

      open(newunit=unit, file=trim(fname), status="old")
      if(mode_read=="xyz") then
         do i=1,nx
            do j=1,ny
               do k=1,nz
                  read(unit, *) xdata(i), ydata(j), zdata(k), fdata(i, j, k)
               end do
            end do
         end do
      elseif(mode_read=="zyx") then
         do k=1,nz
            do j=1,ny
               do i=1,nx
                  read(unit, *) xdata(i), ydata(j), zdata(k), fdata(i, j, k)
               end do
            end do
         end do
      else
         print *, "ERROR: read mode in load interp_3d can be only xyz or zyx, found "//mode_read
         stop
      end if
      close(unit)

      if(dolog) then
         xdata(:) = log10(xdata(:))
         ydata(:) = log10(ydata(:))
         zdata(:) = log10(zdata(:))
         fdata(:, :, :) = log10(fdata(:, :, :))
      end if

      call check_regular_grid(xdata(:), fname)
      call check_regular_grid(ydata(:), fname)
      call check_regular_grid(zdata(:), fname)

      xmin = xdata(1)
      ymin = ydata(1)
      zmin = zdata(1)

      dx = xdata(2) - xdata(1)
      dy = ydata(2) - ydata(1)
      dz = zdata(2) - zdata(1)

      xfact = (nx - 1) / (xdata(nx) - xdata(1))
      yfact = (ny - 1) / (ydata(ny) - ydata(1))
      zfact = (nz - 1) / (zdata(nz) - zdata(1))

      invdx = 1d0 / dx
      invdy = 1d0 / dy
      invdz = 1d0 / dz

   end subroutine load_data_3d

   ! ********************
   subroutine check_regular_grid(xdata, fname)
      implicit none
      real*8,intent(in)::xdata(:)
      real*8::dx, err
      integer::n, i
      character(len=*),optional::fname
      character(len=200)::f_name

      f_name = "unknown file"
      if(present(fname)) then
         f_name = fname
      end if

      dx = xdata(2) - xdata(1)
      n = size(xdata)

      do i=3,n
         err = abs(dx - (xdata(i) - xdata(i-1))) / dx
         if(err > 1d-4) then
            print *, "ERROR: grid is not regular in "//trim(f_name), err, i
            print *, xdata(i), xdata(i-1)
            print *, xdata(i) - xdata(i-1), dx
            stop
         end if
      end do

   end subroutine check_regular_grid

   ! ********************
   function interp_1d(x, fdata, xmin, xfact, invdx) result(f)
      implicit none
      real*8,intent(in)::x, xmin, xfact, invdx, fdata(:)
      real*8::f, xlow
      integer::i

      i = floor((x - xmin) * xfact + 1)

      xlow = (i - 1) / xfact + xmin

      f = (x - xlow) * (fdata(i+1) - fdata(i)) * invdx + fdata(i)

   end function interp_1d

   ! ********************
   function interp_2d(x, y, fdata, xmin, xfact, invdx, &
      ymin, yfact, invdy) result(f)
      implicit none
      real*8,intent(in)::x, xmin, xfact, invdx, fdata(:, :)
      real*8,intent(in)::y, ymin, yfact, invdy
      real*8::f, xlow, ylow, f0, f1, dxi
      integer::i, j

      i = floor((x - xmin) * xfact + 1)
      j = floor((y - ymin) * yfact + 1)

      xlow = (i - 1) / xfact + xmin
      ylow = (j - 1) / yfact + ymin

      dxi = (x - xlow) * invdx

      f0 = dxi * (fdata(i+1, j) - fdata(i, j)) + fdata(i, j)
      f1 = dxi * (fdata(i+1, j+1) - fdata(i, j+1)) + fdata(i, j+1)

      f = (y - ylow) * (f1 - f0) * invdy + f0

   end function interp_2d

   ! ********************
   function interp_3d(x, y, z, fdata, &
      xmin, xfact, invdx, &
      ymin, yfact, invdy, &
      zmin, zfact, invdz) result(f)
      implicit none
      real*8,intent(in)::x, xmin, xfact, invdx, fdata(:, :, :)
      real*8,intent(in)::y, ymin, yfact, invdy
      real*8,intent(in)::z, zmin, zfact, invdz
      real*8::f, xlow, ylow, zlow, f00, f10, f01, f11, f0, f1, dxi, dyi
      integer::i, j, k

      i = floor((x - xmin) * xfact + 1)
      j = floor((y - ymin) * yfact + 1)
      k = floor((z - zmin) * zfact + 1)

      xlow = (i - 1) / xfact + xmin
      ylow = (j - 1) / yfact + ymin
      zlow = (k - 1) / zfact + zmin

      dxi = (x - xlow) * invdx
      dyi = (y - ylow) * invdy

      f00 = dxi * (fdata(i+1, j, k) - fdata(i, j, k)) + fdata(i, j, k)
      f10 = dxi * (fdata(i+1, j+1, k) - fdata(i, j+1, k)) + fdata(i, j+1, k)
      f01 = dxi * (fdata(i+1, j, k+1) - fdata(i, j, k+1)) + fdata(i, j, k+1)
      f11 = dxi * (fdata(i+1, j+1, k+1) - fdata(i, j+1, k+1)) + fdata(i, j+1, k+1)

      f0 = dyi * (f10 - f00) + f00
      f1 = dyi * (f11 - f01) + f01

      f = (z - zlow) * (f1 - f0) * invdz + f0

   end function interp_3d


end module prizmo_fit
