This modification includes some options to enforce the initial tracking direction.  For example, one might want to take the first step in the normal direction from a surface vertex.  While FSL allows one to specify an initial fiber bundle to sample from (```--fibst```), along with a hidden parameter that selects an initial fiber closest to direction of a preferred initial direction (```--prefdirfile```), one cannot explicitly enfore the first step direction.  For a full view of all default FSL options, view the probtrackxOptions.h file.

Within the class ```Particle``` (particle.h), for a single streamline emitting from a single seed point, ```Particle``` performs the tracking.  Within ```Particle```, FSL defines a ```jump``` member method, that takes in two spherical angles ( azimuthal and polar), among other parameters, that moves the streamline one step forward.  **Here, I provide an edit to accept an initial azimuthal and polar angles to enforce the first step.**