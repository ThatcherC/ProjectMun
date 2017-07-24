In my [last post](/flyby-part-1-periapsis/), I came up with an equation relating periapsis height to orbital height, the tangential component of orbital velocity, and total orbital speed:


$$(v\_T^2 - \frac{2GM\_k}{r\_x})r\_p^2 + 2GM\_k r\_p -  v\_\theta^2 r\_x^2 = 0$$

(Quick note: since I'll be using two 'GM' terms, $GM\_k$ will refer to Kerbin's standard gravitational parameter and $GM\_m$ will refer to the Mun's.)

In this post, I'm going to use this to relate some other variables: the angle and velocity of the craft's Mun periapsis:

![vp, thetap picture]

So we've got what I'll call $v\_p$ and $\theta\_p$ and I'm trying to get $v\_T$, $v\_\theta$, and $r\_x$. Once that's all lined up, I'll put those expressions into the equation above and graph it. There should be a bunch of different solutions - a lot of different possible values $v\_p$ and $\theta\_p$ that give the right Kerbin periapsis. The graph should reflect all of these solutions.

####But First, a Few More Derivations

Here's the plan: find $v\_T$, $v\_\theta$, and $r\_x$ in terms of the variables $v\_{mr}$, $v\_{M\theta}$, and $\theta\_M$. Next, find each of those variables in terms of $v\_P$ and $\theta\_P$ and then combine everything so $v\_{mr}$, $v\_{M\theta}$, and $\theta\_M$ go away.

![picture of vectors from notebook with T vector]

Up first: $v\_T$ and $v\_\theta$. Since the components of $v\_T$ are $v\_\theta$ and $v\_r$ and orthogonal to each other,
$$v\_T^2 = v\_r^2 + v\_\theta^2$$
We've got to find $v\_\theta$ anyway, so that can be reused. Nice.

[[explanation]]

$$v\_\theta = v\_{MR} \cos \theta\_m - v\_{M\theta} \sin \theta\_m + T$$

[[check sign on cosine term]]
$$v\_r = v\_{MR} \sin \theta\_m + v\_{M\theta} \cos \theta\_m$$

Since I'll be using a computer to evaluate these expressions, there's no need to analytically square and sum $v\_r$ and $v\_\theta$. I'll just write functions to evaluate them numerically.

Great! So we've got $v\_T$ and $v\_{M\theta}$. Time for $r\_x$. A picture will help:

![rx picture from notebook]

This is a classic Law of Cosines problem. I actually messed it up the first time though - it's always important to double check and make sure that you get answers that make sense. [elaborate?]

Anyway, $r\_x$ comes out to $\sqrt{d\_m^2+r\_{SOI}^2+2 d\_m r\_{SOI} \sin \theta\_m}$.

####Phase 2

Now we've got to define $v\_{MR}$, $v\_{M\theta}$, and $\theta\_M$ in terms of $v\_p$ and $\theta\_p$ and we'll be all set!

$v\_{M\theta}$ is the easiest: applying conservation of angular momentum, we get

$$v\_{M\theta} = \frac{v\_{P} r\_P}{r\_{SOI}}$$

at the edge of the Mun's sphere of influence.

$v\_{MR}$ can be gotten just by applying the equation that keeps coming up:

$$v\_r ( r\_{SOI} ) = \sqrt{ 2( E\_P + \frac{GM\_m}{r\_{SOI}} - \frac{v\_P^2 r\_P^2}{r\_{SOI}^2} ) } $$

where $E\_P$ is the energy of the orbit at the Mun periapsis. That energy is equal to the sum of the kinetic and potential energies of the craft and is conserved throughout the orbit:

$$E = E\_P = \frac{v\_P^2}{2} - \frac{GM\_m}{r\_P}$$

That's not really going to simplify to anything nicer, so when I graph everything, I'll leave $E$ and $v\_r$ as functions of $v\_P$ and $\theta\_P$.

$\theta\_M$ will be the sum of $\theta\_P$ (the angle of the periapsis around the Mun) and a quantity I'll call the deflection angle $\theta\_D$.

![picture of theta\_P]()

A bit of background: if you've got an orbit with the given parameters eccentricity $\epsilon$ and ??major axis?? $r\_0$, an expression for the path of the orbit is

$$r(1-\epsilon\cos\theta) = r\_0$$

...

$$\theta\_D = \pi - \cos^{-1}(-\frac{v\_P^2 r\_P^2 - GM\_m r\_{SOI}}{v\_P^2 r\_P r\_{SOI} - GM\_m r\_{SOI}})$$

Again, this doesn't simplify so it'll be left as a function for a computer to deal with. To get $\theta\_M$, we'll just add $\theta\_P$ and $\theta\_D$. Or so I thought! I had a lot of problems with I tried to test all of this initially. The graph I got kinda made sense visually but was way off when I tested it in KSP. One of the problem I found was that I had defined $\theta\_D$ as just the $\cos^{-1}$ term, not $\pi$ minus that. One other thing that I didn't initally think of, though, was the fact that the Mun moves appreciably throughout this whole maneuver. $\theta\_M+\theta\_D$ gives the correct angle if we assume the Mun isn't moving in the time in which we're flying past it (which could be a reasonable assumption above a certain speed). I found that the Mun moves around .1745 radians or 10 degrees in the time between periapsis and the edge of the Mun's sphere of influence. Once I added that into the $\theta\_M$ expression, things worked out much better.

So really, $\theta\_M$ should be:
$$\theta\_M = \theta\_P + \theta\_D + \tilde .1745$$

$\theta\_M$ also happens to be function of both $v\_P$ and $\theta\_P$.

And now I've got everything I wanted at the beginning! Total speed $v\_T$, tangential velocity $v\_\theta$, and altitude $r\_x$ in the Kerbin frame all defined in terms of Mun periapsis speed $v\_P$ and angle $\theta\_P$. Connecting those are the variables $v\_{MR}$, $v\_{M\theta}$, and $\theta\_M$. It's time for a graph!

####Graphing

My favorite graphing program is [GNUPlot](). It's got a kind of language of its own that makes it very configurable and powerful. What I did for this post was basically to compile all functions we defined above into one file along with some commands to makes everything look nice. This file is available on my [GitHub page for my Mun project](https://github.com/ThatcherC/ProjectMun) in case you want to check it out. What is actually plotted is the big quadratic function from the [previous post](/flyby-part-1-periapsis/) and a plane at z=0. Here's what the file looks like:
![GNUPlot code screenshot]()

And here's what the graph looks like:
![Plot shot]()

Any point along that purple-orange border should give a 35 km Kerbin periapsis if done correctly. Choose an angle (in radains) and slide along the chart until you reach the transition: that speed is the speed you want to have at that angle with an altitude of 10 km above the Mun's surface. In the next post I'll confirm that this works, but play around with it if you want to! GNUPlot is a lot of fun once you get to know it and the concepts in this post should be very applicable to any kind of transfer burn between orbiting bodies - I think you could pretty easily use the equations above to figure out how to get to Duna or another planet if you wanted. Check back soon for a post with much more KSP and a lot fewer equations!
