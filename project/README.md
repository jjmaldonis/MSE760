

Goal of project: Make progress towards figuring out what types of mutations would be good for our Zr-Cu-Al genetic algorithm simulation.
* Should we just be moving a cluster of atoms?
  * If so, is one shape better than another? Chains? Spheres? Ellipsoids?
* Should we be moving atoms at random places in the model?
* Should we be swapping atomic numbers?
* Is one atom better than another? Ie should we keep the 50% of the atoms that have the lowest individual energy assosciated with them? Is there something else that might make one atom better than another?
* What about doing a temperature jump and quench using MD (fast)?

* Try a couple of these ideas and change the number of atoms that we are changing as well.
* I would do an EAM energy minimization for each "mutation idea", recording how many steps it took to minimize the model.
  * Hopefully??? better mutations will minimize the energy faster. Not sure if this makes complete sense yet. May need to revise... (ie what makes a mutation "good"? - probably not E minimization... this could be a project in and of itself I bet... but hey, maybe doing the CGM after each mutation will help? worth a shot I think.)
  * Maybe something like jumping to a totally different basin (not metabasin). Some of my mutations above had better do that. Maybe I could move to a next-door basin (thinking kMC-ART sorta thing here).
  * One person moved atoms in a random direction a random number of times (5-50) while separating unphysically close atoms between each step.


Pseudocode:

```
do X = 1, ??? where X is the number of atoms to be changed in the model

    ! Monte carlo loop
    do while (stop when energy change from step to step is small enough)
        change X number of atoms (optionally close to each other) (eg a cluster of size X)
        do a conjugate gradient minimization (CGM) (LAMMPS)
        calculate energy (LAMMPS)
        accept or reject based on metropolis
    enddo ! MC loop
    record final step # it took to equilibrate

enddo ! X loop
```

