
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


subroutine etot(x,y,z,fx,fy,fz,natoms,boxsize,e1,etotal)
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
    integer, intent(in) :: natoms
    double precision, intent(in) :: boxsize
    double precision, dimension(:), allocatable, intent(out) :: e1
    double precision, intent(out) :: etotal
    integer :: i,j,k,n ! iterables
    integer :: istat
    double precision :: cutoff
    double precision :: rx, ry, rz, r2, r, lj_r, lj_rc, devlj_r, devlj_rc
    double precision :: xi, yi, zi, cutoff2, rx12r2, ry12r2, rz12r2

    if( .not. allocated(fx)) allocate(fx(natoms))
    if( .not. allocated(fy)) allocate(fy(natoms))
    if( .not. allocated(fz)) allocate(fz(natoms))
    fx = 0.0
    fy = 0.0
    fz = 0.0

    cutoff = 2.5
    cutoff2 = cutoff**2
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


subroutine vv_update(x,y,z,vx,vy,vz,fx,fy,fz,natoms,boxsize,timestep,e1,etotal)
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
    double precision, intent(in) :: boxsize
    integer, intent(in) :: natoms
    double precision, intent(in) :: timestep
    double precision, dimension(:), allocatable, intent(out) :: e1
    double precision, intent(out) :: etotal
    integer :: i

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
        call etot(x,y,z,fx,fy,fz,natoms,boxsize,e1,etotal)
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


subroutine msd(x,y,z, xi,yi,zi, natoms, boxsize, t, outvalue)
    implicit none
    double precision, dimension(:), allocatable, intent(inout) :: x
    double precision, dimension(:), allocatable, intent(inout) :: y
    double precision, dimension(:), allocatable, intent(inout) :: z
    double precision, dimension(:), allocatable, intent(inout) :: xi
    double precision, dimension(:), allocatable, intent(inout) :: yi
    double precision, dimension(:), allocatable, intent(inout) :: zi
    integer, intent(in) :: natoms
    double precision, intent(in)  :: boxsize
    double precision, intent(in) :: t ! Time elapsed
    double precision, intent(out) :: outvalue
    integer :: i,j,k
    double precision :: rx,ry,rz,r2,r2sum
    r2sum = 0.0
    do i=1, natoms
        rx = xi(i) - x(i)
        ry = yi(i) - y(i)
        rz = zi(i) - z(i)
        rx = rx - boxsize*anint(rx/boxsize)
        ry = ry - boxsize*anint(ry/boxsize)
        rz = rz - boxsize*anint(rz/boxsize)
        r2 = rx**2+ry**2+rz**2
        r2sum = r2sum + r2
    enddo
    outvalue = r2sum / float(natoms) !(float(natoms) * 6.0*t)
end subroutine msd


subroutine vel_autocorr(vx,vy,vz, vxi,vyi,vzi, natoms, outvalue)
    implicit none
    double precision, dimension(:), allocatable, intent(inout) :: vx
    double precision, dimension(:), allocatable, intent(inout) :: vy
    double precision, dimension(:), allocatable, intent(inout) :: vz
    double precision, dimension(:), allocatable, intent(inout) :: vxi
    double precision, dimension(:), allocatable, intent(inout) :: vyi
    double precision, dimension(:), allocatable, intent(inout) :: vzi
    integer, intent(in) :: natoms
    double precision, intent(out) :: outvalue
    integer :: i,j,k
    double precision :: rx,ry,rz,v2,v2sum
    v2sum = 0.0
    do i=1, natoms
        rx = vxi(i) * vx(i)
        ry = vyi(i) * vy(i)
        rz = vzi(i) * vz(i)
        v2 = rx + ry + rz
        v2sum = v2sum + v2
    enddo
    outvalue = v2sum/float(natoms)
end subroutine vel_autocorr


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
    ! Atomic coordinates, initial to save
    double precision, dimension(:), allocatable :: xi
    double precision, dimension(:), allocatable :: yi
    double precision, dimension(:), allocatable :: zi
    ! Atomic velocities
    double precision, dimension(:), allocatable :: vx
    double precision, dimension(:), allocatable :: vy
    double precision, dimension(:), allocatable :: vz
    ! Atomic velocities, initial to save
    double precision, dimension(:), allocatable :: vxi
    double precision, dimension(:), allocatable :: vyi
    double precision, dimension(:), allocatable :: vzi
    ! Atomic forces
    double precision, dimension(:), allocatable :: fx
    double precision, dimension(:), allocatable :: fy
    double precision, dimension(:), allocatable :: fz
    integer :: natoms ! Number of atoms
    integer :: i,j,k,n ! iterables
    integer :: istat
    integer :: nhis ! Number of bins in gr histogram
    double precision :: temperature
    double precision :: kb = 8.617332478E-5
    double precision :: timestep
    double precision :: epot, ekin
    double precision :: rho ! density info stuff
    double precision, dimension(:), allocatable :: e1 ! Individual energy for each atom
    double precision :: boxsize
    double precision, dimension(:), allocatable :: g ! gr histogram
    double precision :: meansqd ! Mean square displacement
    double precision :: diffusivity, cvv

    nuc = 10
    temperature = 133.0 *kb/eps
    timestep = 0.005

    call generate_lat(nuc, latcon, x, y, z, natoms)
    write(*,*) "Number of atoms:", natoms
    boxsize = nuc*latcon
    write(*,*) "Initial box size =", boxsize

    !call lat_save(x,y,z,natoms,"lattice.xyz",sigma)

    rho = 0.85 ! Desired density (reduced) for input into volume rescaling
    call rescale_vol(x,y,z,natoms,boxsize,rho) ! VOLUME rescaling
    write(*,*) "Rescaled box size to", boxsize

    allocate(xi(natoms))
    allocate(yi(natoms))
    allocate(zi(natoms))
    xi = x
    yi = y
    zi = z

    nhis = 250
    call gr(x,y,z,natoms,boxsize,g,nhis,sigma)
    open(20,file="gr_initial.txt",status='unknown')
    write(20,*) "r, gr"
    do i=1,nhis
        write(20,*) boxsize/nhis*i,g(i)
    enddo
    close(20)

    call generate_velocities(vx,vy,vz,natoms,temperature)
    call calc_temp(temperature,natoms,vx,vy,vz)
    write(*,*) "Initial temperature:", temperature/kb*eps

    ! Need to equilibrate the temp and velocities before setting initial values
    write(*,*) "Equilibrating..."
    call etot(x,y,z,fx,fy,fz,natoms,boxsize,e1,epot)
    do i=1, ceiling(10.0/timestep)
        call vv_update(x,y,z,vx,vy,vz,fx,fy,fz,natoms,boxsize,timestep,e1,epot)
    enddo
    allocate(vxi(natoms))
    allocate(vyi(natoms))
    allocate(vzi(natoms))
    vxi = vx
    vyi = vx
    vzi = vx

    ! Used to calculate gr after volume rescaling
    !call gr(x,y,z,natoms,boxsize,g,nhis,sigma)
    !open(20,file="gr_post_vol_rescale_equilibrate.txt",status='unknown')
    !write(20,*) "r, gr"
    !do i=1,nhis
    !    write(20,*) boxsize/nhis*i,g(i)
    !enddo
    !close(20)
    !stop

    write(*,*) "step, time (ps), etotal, epot, Ecohesive, ekin, temperature, MSD, cvv"
    ! MD Simulation
    ! Calculate accelerations in etot, stores them in fx, fy, fz
    call etot(x,y,z,fx,fy,fz,natoms,boxsize,e1,epot)
    do i=1, ceiling(100.0/timestep)
        call vv_update(x,y,z,vx,vy,vz,fx,fy,fz,natoms,boxsize,timestep,e1,epot)
        call msd(x,y,z,xi,yi,zi,natoms,boxsize,timestep*i,meansqd)
        diffusivity = meansqd/(2.0*timestep*i)
        call vel_autocorr(vx,vy,vz, vxi,vyi,vzi, natoms, cvv)
        epot = sum(e1)*eps ! Convert to real units
        call calc_temp(temperature,natoms,vx,vy,vz)
        ekin = 3.0/2.0*temperature*eps ! Convert to real units
        write(*,"(I5,8F15.9)") i, i*timestep*eps, (epot/natoms+ekin), epot, epot/natoms, ekin, temperature/kb*eps, meansqd, cvv
    enddo

    call calc_temp(temperature,natoms,vx,vy,vz)
    write(*,*) "Final temperature:", temperature/kb*eps

    nhis = 250
    call gr(x,y,z,natoms,boxsize,g,nhis,sigma)
    open(20,file="gr_final.txt",status='unknown')
    write(20,*) "r, gr"
    do i=1,nhis
        write(20,*) boxsize/nhis*i,g(i)
    enddo
    close(20)
end program


