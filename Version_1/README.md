# A Geometric Clustering Tool (AGCT)
This is the Code and Dataset tarball for AGCT and paper  "A Geometric Clustering Tool (AGCT) to robustly unravel the inner cluster structures of time-series gene expressions" by by R. Nock, N. Polouliakh, K. Oka, F. Nielsen, K. Shibanai and H. Kitano

**The version does NOT work with Java 7 or newer.**

## Content
This directory contains the following parts:
* Quickstart pdf presentation explaining the main steps when using AGCT: AGCT_kickstart.pdf
* AGCT paper supplementary information, in the parent directory AGCT_paper_supplementary_information. This pdf provides background knowledge for the algorithms used and explains the features of AGCT.
* Datasets - Data files of the Yeast Cell Cycle (YCC) containing the profiles of the 6178 Genes: Spellman_Cdc28_6178genes.agct. GEO tag file: Yeast_Spellman_tags.txt
  - Smaller datasets (Spellman_Cdc28_1000genes.agct, Spellman_Cdc28_3000genes.agct) for testing and trying out additional algorithms are also provided (e.g. to see that the manifold construction is robust).
* Scenarii - An example scenario for the dataset Spellman_Cdc28_1000genes.agct is provided.
* Source - Directory containing all Java classes needed to compile and run AGCT

## Getting Started
Download the pre-compiled JAR file [agct.jar](https://github.com/agct2019/agct/releases/download/v1.0/agct.jar). To run AGCT, we strongly suggest using a 64-bit JVM with sufficient memory:
```
java -d64 -Xmx8000m -jar agct.jar
```
Please watch the demo on [using a scenario](https://youtu.be/tY-5TeBRq7Y), the demo on [loading datasets manually](https://youtu.be/cQgMaZ1fLlk) and/or read the quickstart that explains the main steps when using AGCT.

### Using a Scenario
AGCT allows saving previously recorded computation to a file, called a *scenario file*. Loading a scenario file will restore the comptation state (algorithms and parameters used) from when the scenario was created. To load a scenario file, simply click on the folder icon of the scenario bar (“load and execute scenario”) and select the scenario file. The scenario will then run automatically until its end.

Please note:
 - When creating a scenario file, please start the recording (third button on scenario bar) **before** loading any datasets.


To run the provided toy scenario file /Version_1/Scenarii/NICTA_AGCT_Spellman_CDC28_1000.txt, you must first edit line 1 to give the correct path. The scenario ends after the computation of the manifold.

### Using the clustering result visualization frame
The clustering result visualization frame shows the results for different clusterings. For it to show anything useful, you have to create multiple clusterings of the same kind. See the `-L` option in the help on clustering commands.

Some graphs take a while to compute. Please be patient.

## Building
To compile the Java files and generate the agct.jar file,
 1. go to the Source repository and run ./compile_all.sh (compiles all Java files), then
 2. run ./makejar.sh, which will generate the jar file and place it in the directory `Version_1`.

You must compile with Java Development Kit 1.6 https://docs.oracle.com/javase/7/docs/webnotes/install/.

## License

See the [LICENSE](../LICENSE) file for details

---

Thank you.