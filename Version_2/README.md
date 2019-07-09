# AGCT: A Geometric Clustering Tool
This is the Code and Dataset tarball for AGCT_vesion2 and paper  "AGCT: A Geometric Clustering Tool for robustly unravelling the inner cluster structures of gene expressions" by by R. Nock, N. Polouliakh, K. Oka, F. Nielsen, K. Shibanai and H. Kitano

## Content
This directory contains the following parts:
* Quickstart pdf presentation explaining the main steps when using AGCT: AGCT_kickstart_vesion2.pdf  
* AGCT paper supplementary information, in the parent directory AGCT_paper_supplementary_information  
* Datasets - Data files of the Yeast Cell Cycle (YCC) and the Yeast Metabolic Cycle (YMC): Yeast_6178.txt (Spellman at al), Yeast_9335.txt (Tu et al), Yeast_3565.txt (periodic genes form Tu et al). GEO tag file: Yeast.txt 
  - A small dataset (Test_570.txt) for testing and trying out additional algorithms is also provided (e.g. to see that the manifold construction is robust). 
* Sentinels - The sentinels folder includes groups of sentinels (probes that are known to be significantly correlated) for the Yeast Metabolic Cycle  (Yeast_9335.txt) and the Yeast Cell Cycle (Yeast_6178.txt).
* Scenarii - An example scenario for the small test dataset is provided.
* AGCT2alpha - The Java source files

## Getting Started
Download the pre-compiled JAR file [AGCT2alpha.jar](https://github.com/agct2018/agct/releases/download/v2.0/AGCT2alpha.jar). To run AGCT, we strongly suggest using a 64-bit JVM with sufficient memory:
```
java -d64 -Xmx8000m -jar AGCT2alpha.jar
```
Please watch the [demo](https://youtu.be/NoCLfP3t1sM) and/or read the quickstart that explains the main steps when using AGCT.

**HINT**: The top right icon blinks when there is computation going on.

### Using sentinels
Sentinels are genes/probes a priori verified as biologically important and whose locations/interaction to other genes might bring new clues/insights into a biological discussion. A user can optionally load them when AGCT asks for an “a priori genelist" during the verification stage. Multiple files from a single directory, or an entire directory should be selected.

### Using a scenario
AGCT allows saving previously recorded computation to a file, called a *scenario file*. Loading a scenario file will restore the comptation state (algorithms and parameters used) from when the scenario was created. To load a scenario file, simply click on the folder icon of the scenario bar (“load and execute scenario”) and select the scenario file. The scenario will then run automatically until its end (in this case, the computation of the non-linear manifold).

Please note:
 * Unlike Version 1, the scenario in AGCT_Vesion2 is recorded automatically upon start. To create a scenario file, the user need only stop recording (third button on scenario bar) and then save the scenario. You must **compute at least one clustering** before saving a scenario file. Otherwise, it will result in an incomplete scenario file.
 * Please always **stop recording** (third button on scenario bar) before loading/saving any scenario files.
 * To limit the size of the file, clusterings are no longer saved in the scenario file. Please **recompute any clusterings** after replaying the scenario file.
To run the provided scenario file /Version_2/Scenarii/Test_570scen.txt, you must edit Line 11 to give the correct path. However, we instead suggest to manually run the small dataset and, if desired, recreate the scenario file on the user's personal computer (AGCT will then automatically generate the right path). To replicate the scenario, please check the file AGCT_kickstart_vesion2.pdf.

## Building
You must compile with Java Development Kit 1.8 https://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html 

Thank you.



