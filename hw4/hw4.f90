
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
implicit none

type hutch
  integer, dimension(:), allocatable :: at
  integer :: nat = 0
end type hutch
type hutch_array
    type(hutch), dimension(:,:,:), allocatable :: h
    integer :: nhutch1d
    double precision:: hutchsize
end type hutch_array

contains

subroutine lat_save(x,y,z,natoms,filename,sigma)
    ! Saves the given lattice positions in outfile.
    ! x, y, and z are assumed to be in reduced Angstroms
    ! outputs are saved in Angstroms
    implicit none
    ! Atomic coordinates
    double precision, dimension(:), allocatable, intent(in) :: x
    double precision, dimension(:), allocatable, intent(in) :: y
    double precision, dimension(:), allocatable, intent(in) :: z
    integer, intent(in) :: natoms
    character (len=*), intent(in) :: filename
    double precision, intent(in) :: sigma
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


subroutine generate_lat(nuc, latcon, x, y, z, natoms)
    ! Generates a lattice based on the inputs, stores in outputs.
    ! NOTE: We assuming a single element model and ignore znum
    implicit none
    integer, intent(in) :: nuc ! Number of unit cells
    double precision, intent(in) :: latcon ! Lattice constant
    ! Atomic coordinates
    double precision, dimension(:), allocatable, intent(out) :: x
    double precision, dimension(:), allocatable, intent(out) :: y
    double precision, dimension(:), allocatable, intent(out) :: z
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
                ! I leave the 0.0 addition in for clarity, not for
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
    !write(*,*) "Box size = ",nuc, "unit cells =", latcon*nuc
end subroutine generate_lat


!subroutine etot_atom(x,y,z,fx,fy,fz,ha,atom,x0,y0,z0,natoms,boxsize,cutoff,e1,etotal,vir)
subroutine etot_atom(x,y,z,atom,x0,y0,z0,natoms,boxsize,cutoff,e1,etotal,vir,vir_atom)
    ! Calculates the energy of the lattice
    ! Calculates the forces within the lattice
    implicit none
    ! Atomic coordinates
    double precision, dimension(:), allocatable, intent(in) :: x,y,z
    double precision, intent(in) :: x0,y0,z0
    !double precision, dimension(:), allocatable, intent(inout) :: fx,fy,fz
    !double precision, dimension(:), intent(inout) :: fx,fy,fz
    integer, intent(in) :: natoms, atom
    double precision, intent(in) :: boxsize
    !double precision, dimension(:), allocatable, intent(out) :: e1
    double precision, dimension(:), intent(out) :: e1, vir_atom
    double precision, intent(out) :: etotal
    double precision, intent(inout) :: vir
    integer :: i,j,k,n ! iterables
    integer :: istat
    double precision, intent(in) :: cutoff
    double precision :: rx, ry, rz, r2, r, lj_r, lj_rc, dudr, dudrc
    double precision :: rx0, ry0, rz0, r02
    double precision :: xi, yi, zi, cutoff2, rx12r2, ry12r2, rz12r2, temp
    integer,  dimension(:), allocatable :: atoms
    !type(hutch_array), intent(inout) :: ha ! hutch array
    integer :: nlist

    cutoff2 = cutoff**2
    !if(.not. allocated(e1)) allocate(e1(natoms))
    !write(*,*) e1(1)
    !if(size(e1) .ne. natoms) then
    !    deallocate(e1)
    !    allocate(e1(natoms))
    !endif
    etotal = 0.0
    !fx=0.0; fy=0.0; fz=0.0;
    !vir = 0.0

    lj_rc = LJ(cutoff)
    dudrc = derivLJ(cutoff)

    e1(atom) = 0.0
    vir_atom(atom) = 0.0
    xi = x(atom)
    yi = y(atom)
    zi = z(atom)
    !call hutch_list_3D(ha, xi, yi, zi, cutoff, boxsize, natoms, atoms, nlist)
    do j=1, natoms
        rx = xi - x(j)
        ry = yi - y(j)
        rz = zi - z(j)
        rx = rx - boxsize*anint(rx/boxsize)
        ry = ry - boxsize*anint(ry/boxsize)
        rz = rz - boxsize*anint(rz/boxsize)
        r2 = rx**2+ry**2+rz**2
        if(j .gt. atom) then ! count part only toward info about 'atom'
            ! This part recalulates E_atom
            if(r2 .le. cutoff2) then ! new position within cutoff of j
                r = sqrt(r2)
                lj_r = LJ(r)
                dudr = derivLJ(r)
                temp = (dudrc-dudr)/r
                e1(atom) = e1(atom) + (lj_r - lj_rc - dudrc*(r-cutoff))
                !fx(atom) = fx(atom) + rx*temp
                !fy(atom) = fy(atom) + ry*temp
                !fz(atom) = fz(atom) + rz*temp
                !fx(j) = fx(j) - rx*temp
                !fy(j) = fy(j) - ry*temp
                !fz(j) = fz(j) - rz*temp
                vir_atom(atom) = vir_atom(atom) + r2*temp
                !vir = vir + (rx**2 + ry**2 + rz**2)*temp/3.0
                !write(*,*) "DEBUG1", r,dudrc, dudr
            endif
        else if(j .lt. atom) then ! count part only toward 'j' < 'atom'
            ! This part correctly modifies the changes to E_j due to 'atom'
            ! being moved
            rx0 = x0 - x(j)
            ry0 = y0 - y(j)
            rz0 = z0 - z(j)
            rx0 = rx0 - boxsize*anint(rx0/boxsize)
            ry0 = ry0 - boxsize*anint(ry0/boxsize)
            rz0 = rz0 - boxsize*anint(rz0/boxsize)
            r02 = rx0**2+ry0**2+rz0**2
            if(r02 .le. cutoff2) then ! old position within cutoff of j
                r = sqrt(r02)
                lj_r = LJ(r)
                dudr = derivLJ(r)
                temp = (dudrc-dudr)/r
                e1(j) = e1(j) - (lj_r - lj_rc - dudrc*(r-cutoff))
                !fx(j) = fx(j) - rx0*temp 
                !fy(j) = fy(j) - ry0*temp
                !fz(j) = fz(j) - rz0*temp
                !fx(atom) = fx(atom) + rx0*temp
                !fy(atom) = fy(atom) + ry0*temp
                !fz(atom) = fz(atom) + rz0*temp
                vir_atom(j) = vir_atom(j) - r2*temp
                !vir = vir - (rx0**2 + ry0**2 + rz0**2)*temp/3.0
            endif
            if(r2 .le. cutoff2) then ! new position within cutoff of j
                r = sqrt(r2)
                lj_r = LJ(r)
                dudr = derivLJ(r)
                temp = (dudrc-dudr)/r
                e1(j) = e1(j) + (lj_r - lj_rc - dudrc*(r-cutoff))
                !fx(j) = fx(j) + rx*temp
                !fy(j) = fy(j) + ry*temp
                !fz(j) = fz(j) + rz*temp
                !fx(atom) = fx(atom) - rx*temp
                !fy(atom) = fy(atom) - ry*temp
                !fz(atom) = fz(atom) - rz*temp
                vir_atom(j) = vir_atom(j) + r2*temp
                !vir = vir + (rx**2 + ry**2 + rz**2)*temp/3.0
            endif
        endif
        !if( .not. e1(j) < 0) then
        !    write(*,*) "HEREATOM2", e1(j),atom,j, r
        !    write(*,*) xi,yi,zi
        !    write(*,*) x(j),y(j),z(j)
        !    stop
        !endif
    enddo
    if( allocated(atoms)) deallocate(atoms)
    etotal = sum(e1)
    vir = sum(vir_atom)/3.0
end subroutine etot_atom


!subroutine etot(x,y,z,fx,fy,fz,ha,natoms,boxsize,cutoff,e1,etotal,vir)
subroutine etot(x,y,z,natoms,boxsize,cutoff,e1,etotal,vir,vir_atom)
    ! Calculates the energy of the lattice
    ! Calculates the forces within the lattice
    implicit none
    ! Atomic coordinates
    double precision, dimension(:), allocatable, intent(in) :: x,y,z
    !double precision, dimension(:), allocatable, intent(inout) :: fx,fy,fz
    integer, intent(in) :: natoms
    double precision, intent(in) :: boxsize
    double precision, dimension(:), allocatable, intent(out) :: e1, vir_atom
    double precision, intent(out) :: etotal, vir
    integer :: i,j,k,n ! iterables
    integer :: istat
    double precision, intent(in) :: cutoff
    double precision :: rx, ry, rz, r2, r, lj_r, lj_rc, devlj_r, devlj_rc
    double precision :: xi, yi, zi, cutoff2, rx12r2, ry12r2, rz12r2, temp
    integer,  dimension(:), allocatable :: atoms
    !type(hutch_array), intent(inout) :: ha ! hutch array
    integer :: nlist

    cutoff2 = cutoff**2
    if(.not. allocated(e1)) allocate(e1(natoms))
    if(size(e1) .ne. natoms) then
        deallocate(e1)
        allocate(e1(natoms))
    endif
    if(.not. allocated(vir_atom)) allocate(vir_atom(natoms))
    if(size(vir_atom) .ne. natoms) then
        deallocate(vir_atom)
        allocate(vir_atom(natoms))
    endif
    e1 = 0.0; etotal = 0.0
    !fx=0.0; fy=0.0; fz=0.0;
    vir = 0.0; vir_atom = 0.0

    lj_rc = LJ(cutoff)
    devlj_rc = derivLJ(cutoff)

    do i=1, natoms
        xi = x(i)
        yi = y(i)
        zi = z(i)
        !call hutch_list_3D(ha, xi, yi, zi, cutoff, boxsize, natoms, atoms, nlist)
        do j=1, natoms
        !do j=1, nlist
            if(j > i) then
            !if(atoms(j) .ne. i) then
                rx = xi - x(j)
                ry = yi - y(j)
                rz = zi - z(j)
                !rx = xi - x(atoms(j))
                !ry = yi - y(atoms(j))
                !rz = zi - z(atoms(j))
                rx = rx - boxsize*anint(rx/boxsize)
                ry = ry - boxsize*anint(ry/boxsize)
                rz = rz - boxsize*anint(rz/boxsize)
                r2 = rx**2+ry**2+rz**2
                if(r2 .le. cutoff2) then
        !write(*,*) i,atoms(j)
        !write(*,*) i,j
                    r = sqrt(r2)
                    lj_r = LJ(r)
                    devlj_r = derivLJ(r)
                    e1(i) = e1(i) + lj_r - lj_rc - devlj_rc*(r-cutoff)
                    temp = (devlj_rc-devlj_r)/r
                    !fx(i) = fx(i) + rx*temp
                    !fy(i) = fy(i) + ry*temp
                    !fz(i) = fz(i) + rz*temp
                    !fx(j) = fx(j) - rx*temp
                    !fy(j) = fy(j) - ry*temp
                    !fz(j) = fz(j) - rz*temp
                    vir_atom(i) = vir_atom(i) + (rx**2 + ry**2 + rz**2)*temp
                    !vir = vir + (rx**2 + ry**2 + rz**2)*temp
                    !fx(atoms(j)) = fx(atoms(j)) - rx*(devlj_rc-devlj_r/r)
                    !fy(atoms(j)) = fy(atoms(j)) - ry*(devlj_rc-devlj_r/r)
                    !fz(atoms(j)) = fz(atoms(j)) - rz*(devlj_rc-devlj_r/r)
                    !if( .not. e1(i) < 0) then
                    !    write(*,*) "HERE1", e1(i),i,j, r2
                    !    write(*,*) xi,yi,zi
                    !    write(*,*) x(j),y(j),z(j)
                    !endif
                endif
            endif
        enddo
        if( allocated(atoms)) deallocate(atoms)
    enddo
    etotal = sum(e1)
    vir = sum(vir_atom)/3.0
end subroutine etot


function LJ(rij)
    ! Done in reduced coords
    implicit none
    double precision, intent(in) :: rij
    double precision :: LJ
    double precision :: ri, ri2, ri6
    ri = 1.0/rij
    ri2 = ri**2
    ri6 = ri2**3
    LJ = 4.0*( ri6**2 - ri6 )
    ! I am going to leave the *4.0 here because the 4 is inherently part of the
    ! LJ potential, and that means if I ever want to replace the potential then
    ! I would have to be really careful to realize that I have an extra *4.0
    ! somewhere else that I need to remove. A single extra calculation in this
    ! isn't the worst thing in the world, and the speed decrease shouldn't be
    ! drastic.
    ! If I was only ever going to use the LJ potential then I wouldn't worry
    ! about it, but there is a good chance I will want to look back at this
    ! code and that could screw me up for a while.
end function LJ


function derivLJ(rij)
    ! Done in reduced coords
    implicit none
    double precision, intent(in) :: rij
    double precision :: derivLJ
    double precision :: ri, ri2, ri6
    ri = 1.0/rij
    ri2 = ri**2
    ri6 = ri2**3
    derivLJ = -48.0*( ri6**2*ri - ri6*ri*0.5 )
    ! I am going to leave the *4.0 here because the 4 is inherently part of the
    ! LJ potential, and that means if I ever want to replace the potential then
    ! I would have to be really careful to realize that I have an extra *4.0
    ! somewhere else that I need to remove. A single extra calculation in this
    ! isn't the worst thing in the world, and the speed decrease shouldn't be
    ! drastic.
    ! If I was only ever going to use the LJ potential then I wouldn't worry
    ! about it, but there is a good chance I will want to look back at this
    ! code and that could screw me up for a while.
end function derivLJ


subroutine gr(xx,yy,zz,natoms,boxsize,g,nhis,sigma)
    implicit none
    double precision, dimension(:), allocatable, intent(in) :: xx
    double precision, dimension(:), allocatable, intent(in) :: yy
    double precision, dimension(:), allocatable, intent(in) :: zz
    integer, intent(in) :: natoms
    double precision, intent(in) :: boxsize
    double precision, intent(in) :: sigma !TODO
    double precision, dimension(:), allocatable, intent(out) :: g ! output histogram
    integer, intent(in) :: nhis ! total number of bins
    double precision :: rho ! density
    integer :: i,j,ig
    double precision :: rx, ry, rz, r, delg, vb, nid
    double precision, parameter :: pi = 3.14159265358979323846264
    if(allocated(g)) deallocate(g)
    allocate(g(nhis)); g = 0
    delg = boxsize/(2*nhis)

    rho = natoms/(boxsize**3)
    
    do i=1,natoms-1
        do j=i+1,natoms
            rx = xx(i)-xx(j)
            ry = yy(i)-yy(j)
            rz = zz(i)-zz(j)
            rx = rx - boxsize*anint(rx/boxsize)
            ry = ry - boxsize*anint(ry/boxsize)
            rz = rz - boxsize*anint(rz/boxsize)
            r = sqrt(rx**2 + ry**2 + rz**2)
            if(r .lt. boxsize/2.0) then
                ig = int(r/delg)
                g(ig) = g(ig) + 2 ! Contribution for particle i and j
            endif
        enddo
    enddo

    ! Normalize gr, ie g
    do i=1,nhis
        r = delg*(i+0.5)
        vb = ((i+1)**3-i**3)*(delg**3)
        nid = (4.0/3.0)*pi*vb*rho
        g(i) = g(i) / (nid*natoms)
    enddo

end subroutine gr


subroutine rescale_vol(xx,yy,zz,natoms,boxsize,rho)
    implicit none
    double precision, dimension(:), allocatable, intent(inout) :: xx
    double precision, dimension(:), allocatable, intent(inout) :: yy
    double precision, dimension(:), allocatable, intent(inout) :: zz
    integer, intent(in) :: natoms
    double precision, intent(inout) :: boxsize
    double precision, intent(in) :: rho ! density
    integer :: i,j
    double precision :: rhoi ! intial density of current model
    double precision :: ratio, scaler

    rhoi = float(natoms)/boxsize**3
    ratio = rhoi/rho
    scaler = ratio**(1.0/3.0)
    xx = xx*scaler
    yy = yy*scaler
    zz = zz*scaler
    boxsize = boxsize*scaler
end subroutine rescale_vol


subroutine randommove(x,y,z,ha,natoms,boxsize,cutoff,maxmove,atom,xo,yo,zo)
    implicit none
    ! Atomic coordinates
    double precision, dimension(:), allocatable, intent(inout) :: x,y,z
    integer, intent(in) :: natoms ! Number of atoms
    double precision, intent(in) :: boxsize
    double precision, intent(in) :: cutoff
    double precision, intent(in) :: maxmove
    integer, intent(out) :: atom ! Atom to be moved
    double precision, intent(inout) :: xo,yo,zo ! 'Old' positions
    double precision :: rand1, rand2, rand3, rand4
    integer :: iseed = 0, i
    double precision :: rx, ry, rz, r2, r, delx, dely, delz, cutoff2
    logical :: good
    type(hutch_array), intent(inout) :: ha ! hutch array

    cutoff2 = cutoff**2

    good = .false. ! To enter loop
    do while( .not. good)
        rand1 = ran(iseed)
        rand2 = ran(iseed)
        rand3 = ran(iseed)
        rand4 = ran(iseed)
        
        atom = int(natoms*rand1)+1
        delx = maxmove*(rand2 - 0.5)
        dely = maxmove*(rand3 - 0.5)
        delz = maxmove*(rand4 - 0.5)

        xo = x(atom)
        yo = y(atom)
        zo = z(atom)

        x(atom) = x(atom) + delx
        y(atom) = y(atom) + dely
        z(atom) = z(atom) + delz
        x(atom) = x(atom) - boxsize*anint(x(atom)/boxsize)
        y(atom) = y(atom) - boxsize*anint(y(atom)/boxsize)
        z(atom) = z(atom) - boxsize*anint(z(atom)/boxsize)

        ! Check if move was within cutoff of all atoms
        good = .true.
        do i=1,natoms
            if(atom .ne. i) then
                rx = x(i)-x(atom)
                ry = y(i)-y(atom)
                rz = z(i)-z(atom)
                rx = rx - boxsize*anint(rx/boxsize)
                ry = ry - boxsize*anint(ry/boxsize)
                rz = rz - boxsize*anint(rz/boxsize)
                r2 = rx**2 + ry**2 + rz**2
                if(r2 .lt. cutoff2) then
                    good = .false.
                    !write(*,*) "Bad one: tried to move atom", atom, i, sqrt(r)
                    exit
                endif
            endif
        enddo

        if(.not. good) then ! Move back
            !write(*,*) "Bad move"
            x(atom) = xo
            y(atom) = yo
            z(atom) = zo
        else
            !write(*,*) "Good move!"
            !call hutch_move_atom(ha, atom, xo,yo,zo, x(atom), y(atom), z(atom), boxsize)
        endif
    enddo

end subroutine randommove

subroutine moveatom(x,y,z,ha,atom,boxsize,xn,yn,zn)
    implicit none
    ! Atomic coordinates
    double precision, dimension(:), allocatable, intent(inout) :: x,y,z
    integer, intent(in) :: atom ! Atom to be moved
    double precision, intent(in) :: xn, yn, zn, boxsize
    type(hutch_array), intent(inout) :: ha ! hutch array
    !call hutch_move_atom(ha, atom, x(atom), y(atom), z(atom), xn, yn, zn, boxsize)
    x(atom) = xn
    y(atom) = yn
    z(atom) = zn
    x(atom) = x(atom) - boxsize*anint(x(atom)/boxsize)
    y(atom) = y(atom) - boxsize*anint(y(atom)/boxsize)
    z(atom) = z(atom) - boxsize*anint(z(atom)/boxsize)
end subroutine moveatom


subroutine init_hutch(ha,x,y,z,natoms,boxsize)
    implicit none
    type(hutch_array), intent(out) :: ha ! hutch array
    double precision, dimension(:), allocatable, intent(in) :: x,y,z
    integer, intent(in) :: natoms ! Atom to be moved
    double precision, intent(in) :: boxsize
    integer :: i,j,k, hx,hy,hz
    ha%nhutch1d = anint( natoms**(1./3.) )
    ha%hutchsize = boxsize/float(ha%nhutch1d)
    allocate(ha%h(ha%nhutch1d,ha%nhutch1d,ha%nhutch1d))
    do i=1, ha%nhutch1d
        do j=1, ha%nhutch1d
            do k=1, ha%nhutch1d
                if(allocated(ha%h(i,j,k)%at)) deallocate(ha%h(i,j,k)%at)
                allocate(ha%h(i,j,k)%at(10))
                ha%h(i,j,k)%nat = 0
            enddo
        enddo
    enddo

    do i=1,natoms
        call hutch_position(ha, x(i), y(i), z(i), boxsize, hx, hy, hz)
        call hutch_add_atom(ha, i, hx, hy, hz)
        write(*,*) "hutch,nat", hx,hy,hz,ha%h(hx,hy,hz)%nat
    enddo
end subroutine init_hutch

subroutine hutch_position(ha, xx, yy, zz,boxsize,hx, hy, hz)
    implicit none
    double precision, intent(in) :: xx, yy, zz
    type(hutch_array), intent(in) :: ha ! hutch array
    double precision, intent(in) :: boxsize
    integer, intent(out) :: hx, hy, hz
    hx = mod(ceiling( (xx + 0.5*boxsize) / ha%hutchsize ), ha%nhutch1d+1)
    hy = mod(ceiling( (yy + 0.5*boxsize) / ha%hutchsize ), ha%nhutch1d+1)
    hz = mod(ceiling( (zz + 0.5*boxsize) / ha%hutchsize ), ha%nhutch1d+1)
    if (hx == 0) hx = 1
    if (hy == 0) hy = 1
    if (hz == 0) hz = 1
end subroutine hutch_position

subroutine hutch_add_atom(ha,atom,hx,hy,hz)
    implicit none
    type(hutch_array), intent(inout) :: ha ! hutch array
    integer, intent(in) :: atom, hx, hy, hz
    ha%h(hx,hy,hz)%nat = ha%h(hx,hy,hz)%nat + 1
    ha%h(hx,hy,hz)%at(ha%h(hx,hy,hz)%nat) = atom
end subroutine hutch_add_atom

subroutine hutch_remove_atom(ha, atom, hx, hy, hz)
    implicit none
    type(hutch_array), intent(inout) :: ha ! hutch array
    integer, intent(in) :: atom, hx, hy, hz
    integer :: i
    do i=1, ha%h(hx,hy,hz)%nat
        if (ha%h(hx,hy,hz)%at(i) .eq. atom) then
            ha%h(hx,hy,hz)%at( i:ha%h(hx,hy,hz)%nat-1 ) = ha%h(hx,hy,hz)%at( i+1:ha%h(hx,hy,hz)%nat )
        endif
    enddo
    ha%h(hx,hy,hz)%nat = ha%h(hx,hy,hz)%nat-1
end subroutine hutch_remove_atom

subroutine hutch_move_atom(ha, atom, xo, yo, zo, xx, yy, zz, boxsize)
    implicit none
    type(hutch_array), intent(inout) :: ha ! hutch array
    integer, intent(in) :: atom
    double precision, intent(in) :: xo, yo, zo
    double precision, intent(in) :: xx, yy, zz
    double precision, intent(in) :: boxsize
    integer :: hx, hy, hz
    call hutch_position(ha, xo, yo, zo, boxsize, hx,hy,hz)
    call hutch_remove_atom(ha, atom, hx, hy, hz)

    call hutch_position(ha, xx, yy, zz, boxsize, hx,hy,hz)
    call hutch_add_atom(ha, atom, hx, hy, hz)
end subroutine hutch_move_atom

subroutine hutch_list_3D(ha, px, py, pz, radius, boxsize, natoms, atoms, nlist)
    implicit none
    type(hutch_array), intent(inout) :: ha ! hutch array
    double precision, intent(in) :: px, py, pz, radius, boxsize
    integer, allocatable, dimension(:) :: atoms
    integer, intent(out) :: nlist
    integer, intent(in) :: natoms
    integer :: i_start, i_end, j_start, j_end, k_start, k_end
    double precision:: x_start, x_end, y_start, y_end, z_start, z_end
    double precision, dimension(3) :: hcenter
    integer, dimension(:), allocatable, target :: temp_atoms
    integer :: i, j, k, hx,hy,hz
    real :: dist2

        allocate(temp_atoms(natoms))
        x_start = px-radius
        x_end = px+radius
        y_start = py-radius
        y_end = py+radius
        z_start = pz-radius
        z_end = pz+radius
        !write(*,*) "radius=", radius
        !write(*,*) "x_start, x_end=", x_start, x_end
        !write(*,*) "y_start, y_end=", y_start, y_end
        !write(*,*) "z_start, z_end=", z_start, z_end
        if(x_start < -boxsize/2.0) x_start = x_start + boxsize !PBC
        if(x_end > boxsize/2.0) x_end = x_end - boxsize !PBC
        if(y_start < -boxsize/2.0) y_start = y_start + boxsize !PBC
        if(y_end > boxsize/2.0) y_end = y_end - boxsize !PBC
        if(z_start < -boxsize/2.0) z_start = z_start + boxsize !PBC
        if(z_end > boxsize/2.0) z_end = z_end - boxsize !PBC
        !write(*,*) "x_start, x_end=", x_start, x_end
        !write(*,*) "y_start, y_end=", y_start, y_end
        !write(*,*) "z_start, z_end=", z_start, z_end
        call hutch_position(ha, x_start, y_start, z_start, boxsize, i_start, j_start, k_start)
        call hutch_position(ha, x_end, y_end, z_end, boxsize, i_end, j_end, k_end)
        !write(*,*) "i_start, i_end=", i_start, i_end
        !write(*,*) "j_start, j_end=", j_start, j_end
        !write(*,*) "k_start, k_end=", k_start, k_end

        !nh = 0
        nlist = 1
        ! There will be a major problem here if i_start > i_end due to the pixel
        ! being out of bounds of the model. Same with j and k.
        do i=1, ha%nhutch1d
            if(i_start <= i_end) then ! This takes care of pbc. It's complicated but it works.
                if(i < i_start .or. i > i_end) cycle
            else
                if(i < i_start .and. i > i_end) cycle
            endif
            do j=1, ha%nhutch1d
                if(j_start <= j_end) then ! This takes care of pbc. It's complicated but it works.
                    if(j < j_start .or. j > j_end) cycle
                else
                    if(j < j_start .and. j > j_end) cycle
                endif
                do k=1, ha%nhutch1d
                    if(k_start <= k_end) then ! This takes care of pbc. It's complicated but it works.
                        if(k < k_start .or. k > k_end) cycle
                    else
                        if(k < k_start .and. k > k_end) cycle
                    endif
                    ! Calculate hutch centers.
                    hcenter(1) = -boxsize/2.0 + ha%hutchsize/2.0 + (i-1)*ha%hutchsize
                    hcenter(2) = -boxsize/2.0 + ha%hutchsize/2.0 + (j-1)*ha%hutchsize
                    hcenter(3) = -boxsize/2.0 + ha%hutchsize/2.0 + (k-1)*ha%hutchsize
                    ! Calculate distance.
                    dist2 = (px-hcenter(1))**2 + (py-hcenter(2))**2 + (pz-hcenter(3))**2
                    if( dist2 < (radius + ha%hutchsize/sqrt(2.0))**2 .or.  dist2 > ((boxsize-radius)*sqrt(3.0))**2 ) then ! The 2nd part is for PBC. It only works if the world is a cube.
                        call hutch_position(ha, hcenter(1), hcenter(2), hcenter(3), boxsize, hx, hy, hz)
                        !used_hutches(hx,hy,hz) = 1
                        if(ha%h(hx, hy, hz)%nat /= 0) then
                            temp_atoms(nlist:nlist+ha%h(hx, hy, hz)%nat-1) = ha%h(hx, hy, hz)%at(1:ha%h(hx, hy, hz)%nat)
                            nlist = nlist + ha%h(hx, hy, hz)%nat
                        endif
                        !nh = nh + 1
                        if(i .ne. hx .or. j .ne. hy .or. k .ne. hz) then
                            write(*,*) "ERROR Hutches:"
                            write(*,*) i,j,k
                            write(*,*) hx, hy, hz
                        endif
                    endif
                enddo
            enddo
        enddo

        ! Copy all the atoms we found in the previous loop into atoms.
        if (nlist > 0) then
            !if( allocated(atoms)) deallocate(atoms)
            allocate(atoms(nlist))
            atoms = temp_atoms(1:nlist)
        else
            error stop
        endif
        deallocate(temp_atoms)
end subroutine hutch_list_3D

subroutine jitter(x,y,z,ha,natoms,boxsize,cutoff,maxmove)
    implicit none
    double precision, dimension(:), allocatable, intent(inout) :: x,y,z
    integer, intent(in) :: natoms ! Number of atoms
    double precision, intent(in) :: boxsize
    type(hutch_array), intent(inout) :: ha ! hutch array
    double precision, intent(in) :: cutoff
    real, intent(in) :: maxmove
    double precision :: xo,yo,zo
    integer :: atom
    double precision :: rand1, rand2, rand3, rand4, delx, dely, delz
    integer :: iseed

    iseed = 2
    do atom=1,natoms
        !rand1 = ran(iseed)
        rand2 = ran(iseed)
        rand3 = ran(iseed)
        rand4 = ran(iseed)
        
        !atom = int(natoms*rand1)+1
        delx = maxmove*(rand2 - 0.5)
        dely = maxmove*(rand3 - 0.5)
        delz = maxmove*(rand4 - 0.5)

        xo = x(atom)
        yo = y(atom)
        zo = z(atom)

        x(atom) = x(atom) + delx
        y(atom) = y(atom) + dely
        z(atom) = z(atom) + delz
        x(atom) = x(atom) - boxsize*anint(x(atom)/boxsize)
        y(atom) = y(atom) - boxsize*anint(y(atom)/boxsize)
        z(atom) = z(atom) - boxsize*anint(z(atom)/boxsize)
    enddo
end subroutine jitter



end module

program hw1
    use lattice
    use IFPORT
    implicit none
    integer :: nuc ! Number of unit cells
    double precision, parameter :: sigma = 3.4 ! For LJ, in Angstroms
    double precision, parameter :: eps = 0.0104 ! For LJ, in eV
    !double precision, parameter :: sigma = 1.0 ! For real untis, in Angstroms
    !double precision, parameter :: eps = 1.0 ! For real units, in eV
    double precision, parameter :: latcon = 5.260/sigma ! Lattice constant in Angstroms
    double precision :: cutoff = 8.5/sigma ! 2.5 in reduced units
    double precision :: maxmove = 0.75/sigma
    !double precision :: movecutoff = 3.3/sigma
    double precision :: movecutoff = 0.0/sigma
    ! Atomic coordinates
    double precision, dimension(:), allocatable :: x,y,z
    !double precision, dimension(:), allocatable :: fx,fy,fz
    integer :: natoms ! Number of atoms
    integer :: i,j,k,n ! iterables
    integer :: istat, iseed
    integer :: nhis ! Number of bins in gr histogram
    double precision :: temperature
    double precision :: kb = 8.617332478E-5
    double precision :: beta
    integer :: step = 0
    double precision :: epot, ekin, eold, enew, delE, emin, optmove, eold_atom, enew_atom
    double precision :: rho ! density info stuff
    double precision, dimension(:), allocatable :: e1, e1_old ! Individual energy for each atom
    double precision, dimension(:), allocatable :: vir_atom, vir_atom_old ! Individual energy for each atom
    double precision :: boxsize
    double precision, dimension(:), allocatable :: g ! gr histogram
    double precision :: randnum
    double precision :: xo,yo,zo
    double precision :: pres, pres_old, vir = 0.0, pres_correction, pkt
    integer :: atom, maxstep
    type(hutch_array) :: ha ! hutch array
    double precision, dimension(100) :: accept = 0.5
    double precision :: acceptance

    nuc = 7
    temperature = 240.0 *kb/eps
    beta = 1.0 / (kb*temperature)

    call generate_lat(nuc, latcon, x, y, z, natoms)
    write(*,*) "Number of atoms:", natoms
    boxsize = nuc*latcon
    write(*,*) "Initial box size (reduced,real)", boxsize, boxsize*sigma
    !call init_hutch(ha,x,y,z,natoms,boxsize)

    !allocate(fx(natoms))
    !allocate(fy(natoms))
    !allocate(fz(natoms))

    call lat_save(x,y,z,natoms,"lattice.xyz",sigma)

    !call jitter(x,y,z,ha,natoms,boxsize,cutoff,0.4)
    rho = natoms/boxsize**3
    write(*,*) "Inital density (reduced,real):", rho, rho/(sigma**3)
    rho = 0.84 ! Desired density (reduced) for input into volume rescaling
    call rescale_vol(x,y,z,natoms,boxsize,rho) ! VOLUME rescaling
    write(*,*) "Rescaled boxsize to (redu,real)", boxsize, boxsize*sigma
    write(*,*) "Rescaled density (reduce,real):", rho, rho/(sigma**3)

    call lat_save(x,y,z,natoms,"lattice-post-rescale-jitter.xyz",sigma)

    nhis = 250
    call gr(x,y,z,natoms,boxsize,g,nhis,sigma)
    open(20,file="gr_initial.txt",status='unknown')
    write(20,*) "r, gr"
    do i=1,nhis
        write(20,*) boxsize/(2*nhis)*i,g(i)
    enddo
    close(20)

    ! Need to equilibrate after changing volume/density
    ! MONTE CARLO
    step = 0
    iseed = 0
    call seed(iseed)
    !call etot(x,y,z,fx,fy,fz,ha,natoms,boxsize,cutoff,e1,eold,vir)
    call etot(x,y,z,natoms,boxsize,cutoff,e1,eold,vir,vir_atom)
    allocate(e1_old(natoms)); e1_old = e1
    allocate(vir_atom_old(natoms)); vir_atom_old = vir_atom
    write(*,*) "Initial energy:", eold
    write(*,*) "Equilibrating..."
    write(*,*) "step, energy, pressure"
    !write(*,*) "step, energy, maxmove"
    pres_correction = 16.0/3.0*3.14159*rho**2* (2.0/3.0*1/cutoff**9 - 1/cutoff**3)
    pkt = rho*kb*temperature
    pres = pkt + vir/boxsize**3 + pres_correction
    pres_old = pres
    do while(step .lt. 1000000)
        call randommove(x,y,z,ha,natoms,boxsize,movecutoff,maxmove,atom,xo,yo,zo)
        !write(*,*) "Moving atom", atom
        !write(*,*) "From", xo,yo,xo
        !write(*,*) "To  ", x(atom),y(atom),z(atom)
        !call etot_atom(x,y,z,fx,fy,fz,ha,atom,xo,yo,zo,natoms,boxsize,cutoff,e1,enew,vir)
        call etot_atom(x,y,z,atom,xo,yo,zo,natoms,boxsize,cutoff,e1,enew,vir,vir_atom)
        !write(*,*) "enew - atom: ", enew, delE
        !write(*,*) "vir-atom: ", vir
        !write(*,*) "sumFx-atom: ", sum(fx)
        !call etot(x,y,z,fx,fy,fz,ha,natoms,boxsize,cutoff,e1,enew,vir)
        !write(*,*) "sumFx-whole: ", sum(fx)
        !write(*,*) "enew - whole:", enew
        !write(*,*) "vir-whole:", vir
        !write(*,*) "step,enew,eold=", step,enew,eold
        delE = enew - eold
        if(delE < 0.0) then
            ! Accept the move
            !write(*,*) "Auto-accept!", delE
            eold = enew
            e1_old = e1
            vir_atom_old = vir_atom
            pres = pkt + vir/boxsize**3 + pres_correction
            pres_old = pres
            !write(*,*) "vir, pkT, correction", vir/boxsize**3, pkt, pres_correction
            accept(mod(step,100)) = 1
        else
            randnum = random(0)
            !write(*,*) "step,randnum=",step, randnum
            !call random(randnum)
            if(log(1.-randnum)<-delE*beta .and. delE .gt. 0.0) then
                ! Accept the move
                !write(*,*) "Accepted due to probability.", log(1.-randnum),-delE*beta
                !write(*,*) "eold,enew", eold, enew
                !write(*,*) "moved atom:", atom
                !write(*,*) x(atom),y(atom),z(atom)
                !write(*,*) xo,yo,zo
                !stop
                eold = enew
                e1_old = e1
                vir_atom_old = vir_atom
                pres = pkt + vir/boxsize**3 + pres_correction
                pres_old = pres
                !write(*,*) "vir, pkT, correction", vir/boxsize**3, pkt, pres_correction, pres_old
                accept(mod(step,100)) = 1
            else
                ! Reject the move
                e1 = e1_old
                vir_atom = vir_atom_old
                call moveatom(x,y,z,ha,atom,boxsize,xo,yo,zo)
                !write(*,*) "Rejected.", log(1.-randnum),-delE*beta
                !write(*,*) "R!", log(1.-randnum),-delE*beta, delE
                accept(mod(step,100)) = 0
            endif
        endif
        if(mod(step,100) .eq. 0) then
            acceptance = sum(accept)/100.0
            if(acceptance .gt. 0.6) then
                maxmove = maxmove*(1.0 + (acceptance-0.5)/5.0)
                !write(*,*) "Changing mm by", (1.0 + (acceptance-0.5)/5.0)
            else if(acceptance .lt. 0.4) then
                if(maxmove .gt. 1.0E-5) maxmove = maxmove*(1.0 - (0.5-acceptance)/5.0)
                !write(*,*) "Changing mm by", (1.0 - (0.5-acceptance)/5.0)
            endif
        endif
        write(*,*) step, eold*eps, pres_old*eps/sigma**3
        !write(*,*) step, eold, maxmove, acceptance
        step = step + 1
    enddo

    write(*,*) "Final energy:", eold

    call lat_save(x,y,z,natoms,"lattice_final.xyz",sigma)

    nhis = 250
    call gr(x,y,z,natoms,boxsize,g,nhis,sigma)
    open(20,file="gr_final.txt",status='unknown')
    write(20,*) "r, gr"
    do i=1,nhis
        write(20,*) boxsize/(2*nhis)*i,g(i)
    enddo
    close(20)


end program


