This is the Code and Dataset tarball for AGCT_vesion2 and paper  "AGCT: A Geometric Clustering Tool for robustly unravelling the inner cluster structures of gene expressions" by by R. Nock, N. Polouliakh, K. Oka, F. Nielsen, K. Shibanai and H. Kitano
Pre-compiled JAR File: https://github.com/agct2018/agct/releases/download/v2.0/AGCT2alpha.jar

Demo Movie (Youtube): https://youtu.be/NoCLfP3t1sM
This tarball contains the following parts: 
* Quickstart pdf presentation to use AGCT: AGCT_kickstart_vesion2.pdf  
* AGCT paper supplementary information, in the parent directory AGCT_paper_supplementary_information  
* Datasets folder   Example file for YCC data: Yeast_6178.txt (Spellman at al), Yeast_9335.txt (Tu et al), Yeast_3565.txt (periodic genes form Tu et al) and Test_570.txt (small dataset for testing).
GEO tag file: Yeast.txt

In addition, we provide datasets containing the first 570 genes of YMC data, if the reviewer wants to do some fast/faster checks, and see e.g., that the manifold construction is robust to the number of genes as well (some operations involving matrices will also be 4 to 100 times faster, which may also be good to try various additional algorithms).

* The sentinels folder includes groups of sentinels (probes that are known to be significantly correlated) for the Yeast Metabolic Cycle and the Yeast Cell Cycle. “Sentinels” here are genes/probes a priori verified as biologically important and whose locations/interaction to other genes might bring new clues/insights into a biological discussion. A user can optionally load them when AGCT asks for an “a priori genelist" during the verification stage. Multiple files from a single directory, or an entire directory should be specified.
** To run AGCT, we strongly suggest using a 64-bit JVM with sufficient memory, for example:  
java -d64 -Xmx8000m -jar AGCT2alpha.jar

* A scenario example of a small test set is provided. However, we suggest to directly run the small dataset and save the scenario file on the user's personal computer instead, to let AGCT automatically generate the right path.
** To run the provided scenario file, you must edit Line 11 in /Version_2/Scenarii/Test.txt to give a correct path.
** To test the provided scenario, simply click on the “Scenario” folder in AGCT (“load and execute scenario”) and select the scenario file in its folder. The scenario will then run automatically until its end (in this case, the computation of the non-linear manifold). 
Please note about scenarios in general:
 (*) The scenario in AGCT_Vesion2 is recorded automatically upon start. To create a scenario file, the user need only stop recording by clicking on the scenario button (third button on scenario bar) and then save the scenario. You must compute at least one clustering before saving a scenario file. Otherwise, it will result in an incomplete scenario file.
 (*) Please always stop the recording before loading/saving any scenario files.
 (*) Clusterings are not saved in the scenario file. Please recompute all clusterings after replaying the scenario file.

* To run AGCT for yourself and replicate the scenario, check file AGCT_kickstart_vesion2.pdf.
Hint: The top right icon blinks when there is computation going on.

* To compile the java files and generate the agct.jar file, (i) go in the Source repository and run ./compile_all.sh (compiles all java files), then (ii) run ./makejar.sh, which will generate the jar file and place it at the right place to run agct.jar as above.
** You must compile with Java Development Kit 1.8 https://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html 

Thank you.



