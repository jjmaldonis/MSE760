
! Compile with:
!	ifort -c hw1.f90
!	ifort hw1.f90 -o hw1
! Or just use the makefile I probably also sent in, which does the same thing.

! There are multiple lines of output code commented out.
! Feel free to uncomment any that you wish to use.
! These usually occur at the end of my subroutines.
! I only left in the cohesive energy as ouputs so I could easily
! collect them for the graphs.
! The main program loops over multiple different number of unit cells
! to generate all the cohesive energies needed. If you are using PBCs,
! it takes ~ a minute to run.

! First in this file you will see a module containing all my subroutines.
! Last is the main program.

module lattice
contains

subroutine lat_save(x,y,z,natoms,filename,sigma)
    ! x, y, and z are assumed to be in reduced Angstroms
    ! outputs are saved in Angstroms
    implicit none
    ! Atomic coordinates
    real, dimension(:), allocatable, intent(in) :: x
    real, dimension(:), allocatable, intent(in) :: y
    real, dimension(:), allocatable, intent(in) :: z
    integer, intent(in) :: natoms
    character (len=*), intent(in) :: filename
    real, intent(in) :: sigma
    integer :: istat ! For error tracking
    integer :: i ! iterables

    open(20, file=trim(filename), iostat=istat, status='unknown')
    write(20,*) natoms
    write(20,*) "Ar lattice xyz file for vesta"
    do i=1,natoms
        write(20,*) "Ar", x(i)*sigma, y(i)*sigma, z(i)*sigma
        !write(20,*) "Ar", x(i), y(i), z(i)
    enddo
    close(20)
end subroutine lat_save


subroutine generate_lat(nuc, latcon, x, y, z, natoms, sigma)
    ! NOTE: We assuming a single element model and ignore znum
    implicit none
    integer, intent(in) :: nuc ! Number of unit cells
    real, intent(in) :: latcon ! Lattice constant
    ! Atomic coordinates
    real, dimension(:), allocatable, intent(out) :: x
    real, dimension(:), allocatable, intent(out) :: y
    real, dimension(:), allocatable, intent(out) :: z
    real, intent(in) :: sigma
    integer, intent(inout) :: natoms ! Number of atoms
    integer :: i,j,k,n ! iterables

    if(allocated(x)) deallocate(x)
    if(allocated(y)) deallocate(y)
    if(allocated(z)) deallocate(z)

    ! Use an FCC lattice
    ! That means natoms per unit cell = 4
    natoms = 4*(nuc**3) ! natoms = 4 atoms / unit cell * nuc^3 unit cells.
    allocate(x(natoms))
    allocate(y(natoms))
    allocate(z(natoms))
    n = 1 
    do i=0, nuc-1
        do j=0, nuc-1
            do k=0, nuc-1
                ! We need to add 4 atoms each time
                ! I leave the 0.0 multiplication in for clarity, not for
                ! speed
                x(n) = (0.0+i)*latcon
                y(n) = (0.0+j)*latcon
                z(n) = (0.0+k)*latcon
                n = n+1
                x(n) = (0.5+i)*latcon
                y(n) = (0.5+j)*latcon
                z(n) = (0.0+k)*latcon
                n = n+1
                x(n) = (0.0+i)*latcon
                y(n) = (0.5+j)*latcon
                z(n) = (0.5+k)*latcon
                n = n+1
                x(n) = (0.5+i)*latcon
                y(n) = (0.0+j)*latcon
                z(n) = (0.5+k)*latcon
                n = n+1
            enddo
        enddo
    enddo
    ! Recenter from - to + coords
    x = x - (latcon*nuc)/2.0
    y = y - (latcon*nuc)/2.0
    z = z - (latcon*nuc)/2.0
    ! Change to reduced units
    x = x/sigma
    y = y/sigma
    z = z/sigma
    !write(*,*) "Box size = ",nuc, "unit cells =", latcon*nuc
end subroutine generate_lat


subroutine etot(x,y,z,natoms,sigma,nuc,eps,latcon)
    implicit none
    ! Atomic coordinates
    real, dimension(:), allocatable, intent(in) :: x
    real, dimension(:), allocatable, intent(in) :: y
    real, dimension(:), allocatable, intent(in) :: z
    real, intent(in) :: sigma, latcon, eps
    integer, intent(in) :: natoms, nuc
    integer :: i,j,k,n ! iterables
    real :: cutoff
    real :: rx, ry, rz, r2, boxsize, etotal
    real, dimension(:), allocatable :: e1

    cutoff = 2.5*sigma
    boxsize = nuc*latcon/sigma
    allocate(e1(natoms))
    e1 = 0; etotal = 0

    do i=1, natoms
        do j=1, natoms
            if(j > i) then
                rx = x(i) - x(j)
                ry = y(i) - y(j)
                rz = z(i) - z(j)
                rx = rx - boxsize*anint(rx/boxsize)
                ry = ry - boxsize*anint(ry/boxsize)
                rz = rz - boxsize*anint(rz/boxsize)
                r2 = rx**2+ry**2+rz**2
                if(r2 .le. cutoff**2) then
                    e1(i) = e1(i) + LJ(sqrt(r2))
                endif
            endif
        enddo
    enddo
    e1 = e1*eps
    etotal = sum(e1)
    !write(*,*) "Total energy of system:", etotal
    !write(*,*) "Cohesive energy:", etotal/natoms
    write(*,*) etotal/natoms
end subroutine etot


subroutine etot_nopbc(x,y,z,natoms,sigma,nuc,eps,latcon)
    implicit none
    ! Atomic coordinates
    real, dimension(:), allocatable, intent(in) :: x
    real, dimension(:), allocatable, intent(in) :: y
    real, dimension(:), allocatable, intent(in) :: z
    real, intent(in) :: sigma, latcon, eps
    integer, intent(in) :: natoms, nuc
    integer :: i,j,k,n ! iterables
    real :: cutoff
    real :: rx, ry, rz, r2, boxsize, etotal
    real, dimension(:), allocatable :: e1

    cutoff = 2.5*sigma
    boxsize = nuc*latcon/sigma
    allocate(e1(natoms))
    e1 = 0; etotal = 0

    do i=1, natoms
        do j=1, natoms
            if(j > i) then
                rx = x(i) - x(j)
                ry = y(i) - y(j)
                rz = z(i) - z(j)
                !if(boxsize*anint(rx/boxsize) .ne. 0.0) then
                !    write(*,*) "Warning!  Needed pbcs on x:", rx, x(i), x(j)
                !endif
                !if(boxsize*anint(ry/boxsize) .ne. 0.0) then
                !    write(*,*) "Warning!  Needed pbcs on y:", ry, y(i), y(j)
                !endif
                !if(boxsize*anint(rz/boxsize) .ne. 0.0) then
                !    write(*,*) "Warning!  Needed pbcs on z:", rz, z(i), z(j)
                !endif
                !rx = rx - boxsize*anint(rx/boxsize)
                !ry = ry - boxsize*anint(ry/boxsize)
                !rz = rz - boxsize*anint(rz/boxsize)
                r2 = rx**2+ry**2+rz**2
                if(r2 .le. cutoff**2) then
                    e1(i) = e1(i) + LJ(sqrt(r2))
                endif
            endif
        enddo
    enddo
    e1 = e1*eps
    etotal = sum(e1)
    !write(*,*) "Total energy of system:", etotal
    !write(*,*) "Cohesive energy:", etotal/natoms
    write(*,*) etotal/natoms
end subroutine etot_nopbc


function LJ(rij)
    ! Done in reduced coords
    implicit none
    real, intent(in) :: rij
    real :: LJ
    LJ = 4*( (1/rij)**12 - (1/rij)**6 )
end function LJ
end module

program hw1
    use lattice
    implicit none
    !integer, parameter :: nuc = 5 ! Number of unit cells
    integer :: nuc ! Number of unit cells
    real, parameter :: latcon = 5.260 ! Lattice constant in Angstroms
    !real, parameter :: latcon = 1.0 ! Debugging
    real, parameter :: sigma = 3.4 ! For LJ, in Angstroms
    real, parameter :: eps = 0.0104 ! For LJ, in eV
    ! Atomic coordinates
    real, dimension(:), allocatable :: x
    real, dimension(:), allocatable :: y
    real, dimension(:), allocatable :: z
    integer :: natoms ! Number of atoms
    integer :: i,j,k,n ! iterables

    do nuc=1, 20
    call generate_lat(nuc, latcon, x, y, z, natoms, sigma)
    call lat_save(x,y,z,natoms,"lattice.xyz", sigma)
    call etot(x,y,z,natoms,sigma,nuc,eps,latcon)
    enddo
end program


