Plan

Time of do a definitive test of all functions:
1. [ ] Get intercept?
2  [ ] Find orbit, esp - found a big problem! still might want ot check ra calc
    [x] Outbound
    [ ] Inbound
3. [x] getMunSOI time
4. Check getMunRms
    Doesn't match game's munRM
    Might not actually calculate RM either...... if ra is 0 how is rm non-zero??

Before rotation velocities are spot on.
After rotation entry velocity is a perfect match to KSPs. Exit velocity, however, is not.
In fact it is not even consistent with the entry velocity (total velocity off by .8m/s)
Let's correct this by whatever method we used to match exit velocity to begin with

2592744 - KSP's Mun periapsis until normal conditions (1:8:3:0:0, 100K, 35K, 10K)


Book says to change rtml and rtmr to offset each relative-velocity vector by the amount ra
How are relative-velocity vectors offset by a scalar distance?

Entry point rotation works well on moving actual ra value to target
Perhaps a small variation in thetaFL would is needed to get a good periapsis???


Notes on the graphical solver -
1 Choose a thetaFL and a Mun intercept position (also an angle)
2. Vary one (thetaFL, say) until proper RM is achieved
3. check rk
4. Vary other (posAngle, ie) such that rk is achieved (Newton's method)
