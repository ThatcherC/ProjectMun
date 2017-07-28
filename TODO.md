TODO:
- [x] Figure out if epsilon2 is in a righthanded frame
    * run two simulations, one on each side of the impact region
    * observe calculated and actual eps2
    * ---------Results---------
    * Seems to be left handed: negative eps has a positive mun-int-peri angle
    * Depends on the naming though... - a negative eps corresponds to a negative (CW) rotation
    * For a negative eps, we need to **add** the deflection angle to the entry angle
- [ ] Calculate exit angle from Mun
    * add/subtract entry angle and deflection angle depending on eps2 sign
    * depends on handedness of eps2
    * ---------Results---------
    * When eps is negative, `exitAngle = entryAngle+deflectionAngle`
- [ ] Calculate exit speed and angle in Kerbin frame
- [ ] Find out if PyPlot has a non-value indicator for handling undefined behavior
- [ ] Look into PyPlot interactivity

- [x] refactor the git repo

- [x] Relook at autopilots
- [x] Provide good documentation for the autopilots
- [ ] Choose the best orbit selection method for the autopilot - maybe give initial phase angle, like in an early version? maybe that version still exists

- [ ] explore orbits
    * conventional free return
    * slingshot - identify by kerbin speed and angle after intercept
    * that free return orbit that goes way out and then infront of the moon
      * play with KSP antenna mechanics
