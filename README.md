# CancerSimulationR
Interactive simulation activity for cancer progression and evolution.

# Navigation
[Installing CancerSimulationR](#install)  
[Updating CancerSimulationR](#update)  
[Running CancerSimulationR](#howtorun)  


## Getting Started
CancerSimulationR was developed and tested on a 64-bit computer with a 3.1-GHz processor running Mac OS X (Version 10.13.6) with 16 GB of RAM.

The instructions below demonstrate how to install this package directly from Github to get the latest release.

### Software and Package Prerequisites:
Install version 3.6.0 or later of R. Users can install R by downloading the appropriate R-x.y.z.tar.gz  file from http://www.r-project.org and following the system-specific instructions. CancerSimulationR was developed and tested on version 3.6.0 of R. As of this release, we recommend using version 3.6.0.

CancerSimulationR depends on the following R libraries: RColorBrewer (version 1.1.2), ggplot2 (version 3.2.1), reshape (version 0.8.8), and shiny (version 1.3.2). The versions provided are the R package versions for which this CancerSimulationR code has been tested.

In order to install a package from github, you will need the devtools package. You can install this package with the following commands:

```
install.packages("devtools")
library(devtools)
```

VirusSimulationR package depends on several packages, which can be installed using the below commands:

```
install.packages("RColorBrewer") 
install.packages("ggplot2") 
install.packages("reshape") 
install.packages("shiny") 
```

<a name="install"></a>
### Installing CancerSimulationR:

To currently get the CancerSimulationR R package up and working on your computer, once you have installed all package dependencies (see above):

1. Open R studio and load devtools using `library(devtools)`. If you don't have devtools you may have to install it with `install.packages("devtools")` and then use `library(devtools)`.
2. Type the following into R studio: `install_github(repo = "mesako/CancerSimulationR")`. 
3. This should start installing all library dependencies so it may take a bit to finish. Check that it finishes without ERROR messages, though it may print WARNINGS.

Typical installation time should take no more than 5 minutes for the most up-to-date CancerSimulationR package. However, total installation time will vary depending on the installation time of other required packages and the speed of your internet connection.

<a name="update"></a>
### Updating CancerSimulationR:

To quickly update your CancerSimulationR R package up and get the latest version from GitHub:

1. Open R studio and load devtools using `library(devtools)`.
2. Type the following into R studio: `install_github(repo = "mesako/CancerSimulationR")`.
3. Load VirusSimulationR using `library(CancerSimulationR)`.

If the above commands run without error, you should have the latest version of CancerSimulationR.

## Running CancerSimulationR
<a name="howtorun"></a>

# TO BE ADDED LATER



