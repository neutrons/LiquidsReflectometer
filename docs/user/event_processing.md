# Event processing
The BL4B instrument leverages the concept of weighted events for several aspects of
the reduction process. Following this approach, each event is treated separately and is
assigned a weigth $w$ to accound for various corrections. Summing events then becomes the
sum of the weights for all events.

## Loading events and dead time correction
A dead time correction is available for rates above around 2000 counts/sec. Both
paralyzing and non-paralyzing implementation are available. Paralyzing refers to a detector
that extends its dead time period when events occur while the detector is already unavailable
to process events, while non-paralyzing refers to a detector that always becomes available
after the dead time period [1].

The dead time correction to be multiplied by the measured detector counts is given by
the following for the paralyzing case:

$$
C_{par} = -{\cal Re}W_0(-R\tau/\Delta_{TOF}) \Delta_{TOF}/R
$$

where $R$ is the number of triggers per accelerator pulse within a time-of-flight bin $\Delta_{TOF}$.
The dead time for the current BL4B detector is $\tau=4.2$ $\mu s$. In the equation avove, ${\cal Re}W_0$ referes to the principal branch of the Lambert W function.

The following is used for the non-paralyzing case:

$$
C_{non-par} = 1/(1-R\tau/\Delta_{TOF})
$$

By default, we use a paralyzing dead time correction with $\Delta_{TOF}=100$ $\mu s$. These parameters can be changed.

The BL4B detector is a wire chamber with a detector readout that includes digitization of the
position of each event. For a number of reasons, like event pileup, it is possible for the
electronics to be unable to assign a coordinate to a particular trigger event. These events are
labelled as error events and stored along with the good events. While only good events are used
to compute reflectivity, error events are included in the $R$ value defined above. For clarity, we chose to define $R$ in terms of number of triggers as opposed to events.

Once the dead time correction as a function for time-of-flight is computed, each event
in the run being processed is assigned a weight according to the correction.

$w_i = C(t_i)$

where $t_i$ is the time-of-flight of event $i$. The value of $C$ is interpolated from the
computed dead time correction distribution.

[1] V. Bécares, J. Blázquez, Detector Dead Time Determination and OptimalCounting Rate for a Detector Near a Spallation Source ora Subcritical Multiplying System, Science and Technology of Nuclear Installations, 2012, 240693, https://doi.org/10.1155/2012/240693

### Correct for emission time
Since neutrons of different wavelength will spend different amount of time on average
within the moderator, a linear approximation is used by the data acquisition system to
account for emission time when phasing choppers.

The time of flight for each event $i$ is corrected by an small value given by

$\Delta t_i = -t_{off} +  \frac{h L}{m_n} A  t_i $

where $h$ is Planck's constant, $m_n$ is the mass of the neutron, and $L$ is the distance
between the moderator and the detector.

The $t_{off}$, $A$, and $L$ parameters are process variables that are stored in the
data file and can be changed in the data acquisition system.

### Gravity correction
The reflected angle of each neutron is corrected for the effect of gravity according to
reference Campbell et al [2]. This correction is done individually for each neutron event according to its wavelength.

[2] R.A. Campbell et al, Eur. Phys. J. Plus (2011) 126: 107. https://doi.org/10.1140/epjp/i2011-11107-8

### Event selection
Following the correction described above, we are left with a list of events, each having
a detector position ($p_x, p_y$) and a wavelength $\lambda$.
As necessary, regions of interests can be defined to identify events to include in the specular
reflectivity calculation, and which will be used to estimate and subtract background.
Event selection is performed before computing the reflectivity as described in the following sections.

### Q calculation
The reflectivity $R(q)$ is computed by computing the $q$ value for each even and histogramming
in a predefined binning of the user's choice. This approach is slightly different from the
traditional approach of binning events in TOF, and then converting the TOF axis to $q$.
The event-based approach allows us to bin directly into a $q$ binning of our choice and avoid
the need for a final rebinning.

The standard way of computing the reflected signal is simply to compute $q$ for each event $i$
using the following equation:

$q_{z, i} = \frac{4\pi}{\lambda_i}\sin(\theta - \delta_{g,i})$

where the $\delta_{g,i}$ refers to the angular offset caused by gravity.

Once $q$ is computed for each neutron, they can be histogrammed, taking into account the
weight assigned to each event:

$S(q_z) = \frac{1}{Q} \sum_{i \in q_z \pm \Delta{q_z}/2}  w_i$

where the sum is over all event falling in the $q_z$ bin or width $\Delta q_z$, and $w_i$ is the
weight if the $i^{th}$ event. At this point we have an unnormalized $S(q_z)$, which remains to be
corrected for the neutron flux. The value of $Q$ is the integrated proton charge for the

### Constant-Q binning
When using a divergent beam, or when measuring a warped sample, it may be beneficial to take
into accound where a neutron landed on the detector in order to recalculate its angle, and its
$q$ value.

In this case, the $q_{z, i}$ equation above becomes:

$q_{z, i} = \frac{4\pi}{\lambda_i}\sin(\theta + \delta_{f,i} - \delta_{g,i})$

where $\delta_{f,i}$ is the angular offset between where the specular peak appears on the
detector and where the neutron was detected:

$\delta_{f,i} = \mathrm{sgn}(\theta)\arctan(d(p_i-p_{spec})/L_{det})/2$

where $d$ is the size of a pixel, $p_i$ is the pixel where event $i$ was detected,
$p_{spec}$ is the pixel at the center of the peak distribution, $L_{det}$ is the distance
between the sample and the detector. Care should be taken to asign the correct sign to
the angle offset. For this reason, we add the sign the scattering angle $\mathrm{sgn}(\theta)$ on from of the
previous equation to account for when we reflect up or down.


## Normalization options
The scattering signal computed above needs to be normalized by the incoming flux in order
to produce $R(q_z)$. For the simplest case, we follow the same procedure as above for the
relevant direct beam run, and simply compute the $S_1(q_z)$ using the standard procedure above,
using the same $q_z$ binning,
and replacing $\theta$ by the value at which the reflected beam was measured. We are then effectively computing what the measured signal would be if all neutron from the beam would reflect with a probability of 1. We refer this distribution at $S_1(q_z)$.

The measured reflectivity then becomes

$$
    R(q_z) = S(q_z) / S_1(q_z)
$$

This approach is equivalent to predetermining the TOF binning that would be needed to produce
the $q_z$ binning we actually want, summing counts in TOF for both scattered and direct beam,
taking the ratio of the two, and finally converting TOF to $q_z$. The only difference is that we
don't bother with the TOF bins and assign events directly into the $q_z$ we know they will contribute to the denominator of for normalization.

### Normalization using weighted events
An alternative approach to the normalization described above is also implemented to BL4B.
It leverages the weighted event approach. Using this approach, we can simply histogram the direct
beam event in a wavelenth distribution. In such a histogram, each bin in wavelength will have
a flux

$$\phi(\lambda) = N_{\lambda} / Q / \Delta_{\lambda}$$

where $N_{\lambda}$ is the number of neutrons in the bin of center $\lambda$, $Q$ is the
integrated proton charge, and $\Delta(\lambda)$ is the wavelength bin width for the distribution.

Coming back to the calculation of the reflected signal above, we now can add a new weight for
each event according to the flux for its particular wavelength:

$$
w_i \rightarrow w_i / \phi(\lambda_i) q_{z,i} / \lambda_i
$$

where $\phi(\lambda)$ is interpolated from the distribution we measured above. The $q_z/\lambda$
term is the Jacobian to account for the transformation of wavelength to $q$.
With this new weigth, we can compute reflectivity directly from the $S(q_z)$ equation above:

$$
R(q_z) = \frac{1}{Q} \sum_{i \in q_z \pm \Delta{q_z}/2}  w_i / \phi(\lambda_i) q_{z,i} / \lambda_i
$$
