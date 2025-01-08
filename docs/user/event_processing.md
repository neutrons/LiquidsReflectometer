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

$\Delta t_i = -t_{off} + t_{mult} * t_i * h L / m_n$

where $h$ is Planck's constant, $m_n$ is the mass of the neutron, and $L$ is the distance
between the moderator and the detector.

The $t_{off}$, $t_{mult}$, and $L$ parameters are process variables that are stored in the
data file and can be changed in the data acquisition system.


### Gravity correction


### Q calculation

### Constant-Q binning

## Normalization options