
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
    !write(*,*) "Box size = ",nuc, "unit cells =", latcon*nuc
end subroutine generate_lat


subroutine etot(x,y,z,fx,fy,fz,natoms,nuc,latcon,e1,etotal)
    ! Calculates the energy of the lattice
    ! Calculates the forces within the lattice
    implicit none
    ! Atomic coordinates
    double precision, dimension(:), allocatable, intent(in) :: x
    double precision, dimension(:), allocatable, intent(in) :: y
    double precision, dimension(:), allocatable, intent(in) :: z
    double precision, dimension(:), allocatable, intent(inout) :: fx
    double precision, dimension(:), allocatable, intent(inout) :: fy
    double precision, dimension(:), allocatable, intent(inout) :: fz
    double precision, intent(in) :: latcon
    integer, intent(in) :: natoms, nuc
    double precision, dimension(:), allocatable, intent(out) :: e1
    double precision, intent(out) :: etotal
    integer :: i,j,k,n ! iterables
    integer :: istat
    double precision :: cutoff
    double precision :: rx, ry, rz, r2, r, boxsize, lj_r, lj_rc, devlj_r, devlj_rc
    double precision :: xi, yi, zi, cutoff2, rx12r2, ry12r2, rz12r2

    if( .not. allocated(fx)) allocate(fx(natoms))
    if( .not. allocated(fy)) allocate(fy(natoms))
    if( .not. allocated(fz)) allocate(fz(natoms))
    fx = 0.0
    fy = 0.0
    fz = 0.0

    cutoff = 2.5
    cutoff2 = cutoff**2
    boxsize = nuc*latcon
    if(.not. allocated(e1)) allocate(e1(natoms))
    if(size(e1) .ne. natoms) then
        deallocate(e1)
        allocate(e1(natoms))
    endif
    e1 = 0; etotal = 0

    lj_rc = LJ(cutoff)
    devlj_rc = derivLJ(cutoff)
    do i=1, natoms
        xi = x(i)
        yi = y(i)
        zi = z(i)
        do j=1, natoms
            if(j > i) then
                rx = xi - x(j)
                ry = yi - y(j)
                rz = zi - z(j)
                rx = rx - boxsize*anint(rx/boxsize)
                ry = ry - boxsize*anint(ry/boxsize)
                rz = rz - boxsize*anint(rz/boxsize)
                r2 = rx**2+ry**2+rz**2
                r = sqrt(r2)
                if(r2 .le. cutoff2) then
                    lj_r = LJ(r)
                    devlj_r = derivLJ(r)
                    e1(i) = e1(i) + lj_r - lj_rc - devlj_rc*(r-cutoff)
                    !e1(i) = e1(i) + lj_r
                    fx(i) = fx(i) + rx*(devlj_rc-devlj_r/r)
                    fy(i) = fy(i) + ry*(devlj_rc-devlj_r/r)
                    fz(i) = fz(i) + rz*(devlj_rc-devlj_r/r)
                    fx(j) = fx(j) - rx*(devlj_rc-devlj_r/r)
                    fy(j) = fy(j) - ry*(devlj_rc-devlj_r/r)
                    fz(j) = fz(j) - rz*(devlj_rc-devlj_r/r)
                endif
            endif
        enddo
    enddo
end subroutine etot


subroutine vv_update(x,y,z,vx,vy,vz,fx,fy,fz,natoms,nuc,latcon,timestep,e1,etotal)
    ! Velocity-Verlet single step
    implicit none
    double precision, dimension(:), allocatable, intent(inout) :: x
    double precision, dimension(:), allocatable, intent(inout) :: y
    double precision, dimension(:), allocatable, intent(inout) :: z
    double precision, dimension(:), allocatable, intent(inout) :: vx
    double precision, dimension(:), allocatable, intent(inout) :: vy
    double precision, dimension(:), allocatable, intent(inout) :: vz
    double precision, dimension(:), allocatable, intent(inout) :: fx
    double precision, dimension(:), allocatable, intent(inout) :: fy
    double precision, dimension(:), allocatable, intent(inout) :: fz
    double precision, intent(in) :: latcon
    integer, intent(in) :: natoms, nuc
    double precision, intent(in) :: timestep
    double precision, dimension(:), allocatable, intent(out) :: e1
    double precision, intent(out) :: etotal
    integer :: i
    double precision :: boxsize

    boxsize = nuc*latcon

    do i=1, natoms
        ! Update velocities
        vx(i) = vx(i) + fx(i)*timestep/2.0
        vy(i) = vy(i) + fy(i)*timestep/2.0
        vz(i) = vz(i) + fz(i)*timestep/2.0
        ! Update positions
        x(i) = x(i) + vx(i)*timestep
        y(i) = y(i) + vy(i)*timestep
        z(i) = z(i) + vz(i)*timestep
        x(i) = x(i) - boxsize*anint(x(i)/boxsize)
        y(i) = y(i) - boxsize*anint(y(i)/boxsize)
        z(i) = z(i) - boxsize*anint(z(i)/boxsize)
    enddo
        ! Recompute accelerations
        call etot(x,y,z,fx,fy,fz,natoms,nuc,latcon,e1,etotal)
    do i=1, natoms
        ! Update velocities
        vx(i) = vx(i) + fx(i)*timestep/2.0
        vy(i) = vy(i) + fy(i)*timestep/2.0
        vz(i) = vz(i) + fz(i)*timestep/2.0
    enddo
end subroutine vv_update

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


function KD(i,j)
    ! Kroneckers delta function
    implicit none
    integer, intent(in) :: i,j
    double precision :: KD
    if(i .eq. j) then
        KD = 1
    else
        KD = 0
    endif
end function KD


subroutine generate_velocities(vx,vy,vz,natoms,temperature)
    ! Temperature should be sent in in reduced units, not real
    implicit none
    ! Atomic coordinates
    double precision, dimension(:), allocatable, intent(inout) :: vx
    double precision, dimension(:), allocatable, intent(inout) :: vy
    double precision, dimension(:), allocatable, intent(inout) :: vz
    integer, intent(in) :: natoms
    !double precision,  intent(in) :: temperature
    double precision,  intent(inout) :: temperature
    double precision :: temp
    integer :: i,j,k,n ! iterables
    double precision :: q1, q2, s2, v0
    integer :: seed = 3

    temp = temperature

    if( .not. allocated(vx)) allocate(vx(natoms))
    if( .not. allocated(vy)) allocate(vy(natoms))
    if( .not. allocated(vz)) allocate(vz(natoms))
    vx = 0.0
    vy = 0.0
    vz = 0.0

    v0 = sqrt(3.0*temperature)
    do i=1, natoms
        s2 = 2.0
        do while(s2 .ge. 1)
            q1 = 2.0*ran(seed)-1.0
            q2 = 2.0*ran(seed)-1.0
            s2 = (q1*q1)+(q2*q2)
        enddo
        vx(i) = 2.0*sqrt(1-s2)*q1*v0
        vy(i) = 2.0*sqrt(1-s2)*q2*v0
        vz(i) = (1.0-(2.0*s2))*v0
        !write(*,*) "DEBUG0: |V(i)| = ", vx(i)**2 + vy(i)**2 + vz(i)**2
    enddo
end subroutine generate_velocities


subroutine calc_temp(temperature,natoms,vx,vy,vz)
    implicit none
    double precision, intent(out) :: temperature
    integer, intent(in) :: natoms
    double precision, dimension(:), allocatable, intent(in) :: vx
    double precision, dimension(:), allocatable, intent(in) :: vy
    double precision, dimension(:), allocatable, intent(in) :: vz
    integer :: i
    temperature = 0.0
    do i=1, natoms
        temperature = temperature + ( vx(i)**2 + vy(i)**2 + vz(i)**2 )
    enddo
    temperature = temperature/(3.0*natoms)
end subroutine calc_temp

end module

program hw1
    use lattice
    implicit none
    integer :: nuc ! Number of unit cells
    double precision, parameter :: sigma = 3.4 ! For LJ, in Angstroms
    double precision, parameter :: eps = 0.0104 ! For LJ, in eV
    double precision, parameter :: latcon = 5.260/sigma ! Lattice constant in Angstroms
    ! Atomic coordinates
    double precision, dimension(:), allocatable :: x
    double precision, dimension(:), allocatable :: y
    double precision, dimension(:), allocatable :: z
    ! Atomic velocities
    double precision, dimension(:), allocatable :: vx
    double precision, dimension(:), allocatable :: vy
    double precision, dimension(:), allocatable :: vz
    ! Atomic forces
    double precision, dimension(:), allocatable :: fx
    double precision, dimension(:), allocatable :: fy
    double precision, dimension(:), allocatable :: fz
    integer :: natoms ! Number of atoms
    integer :: i,j,k,n ! iterables
    integer :: istat
    double precision :: temperature
    double precision :: kb = 8.617332478E-5
    double precision :: timestep
    double precision :: epot, ekin
    double precision, dimension(:), allocatable :: e1 ! Individual energy for each atom

    nuc = 7
    temperature = 300.0 *kb/eps
    timestep = 0.02

    call generate_lat(nuc, latcon, x, y, z, natoms)
    write(*,*) "Number of atoms:", natoms
    call lat_save(x,y,z,natoms,"lattice.xyz", sigma)

    call generate_velocities(vx,vy,vz,natoms,temperature)

    call calc_temp(temperature,natoms,vx,vy,vz)
    write(*,*) "Temperature:", temperature/kb*eps


    write(*,*) "step, time (ps), etotal, epot, Ecohesive, ekin, temperature"
    ! MD Simulation
    ! Calculate accelerations in etot, stores them in fx, fy, fz
    call etot(x,y,z,fx,fy,fz,natoms,nuc,latcon,e1,epot)
        ! Write to file the first 50 forces
        open(30, file=trim('forces.txt'), iostat=istat, status='unknown')
        write(30,*) "Real (not reduced) forces for first 50 atoms:"
        do i=1, 50
            write(30,*) fx(i)*eps/sigma, fy(i)*eps/sigma, fz(i)*eps/sigma
        enddo
        close(30)
    do i=1, ceiling(20.0/timestep)
        call vv_update(x,y,z,vx,vy,vz,fx,fy,fz,natoms,nuc,latcon,timestep,e1,epot)
        epot = sum(e1)*eps ! Convert to real units
        call calc_temp(temperature,natoms,vx,vy,vz)
        !ekin = 3.0/2.0*kb*temperature*eps ! Convert to real units
        ekin = 3.0/2.0*temperature*eps ! Convert to real units
        !write(*,*) 3.0/2.0*kb*temperature*eps, 3.0/2.0*temperature*eps**2
        write(*,"(I4,6F15.9)") i, i*timestep*eps*1000, (epot/natoms+ekin), epot, epot/natoms, ekin, temperature/kb*eps
    enddo
end program


