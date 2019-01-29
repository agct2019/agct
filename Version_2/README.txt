This is the Code and Dataset tarball for AGCT_vesion2 and paper  "AGCT: A Geometric Clustering Tool for robustly unravelling the inner cluster structures of gene expressions" by by R. Nock, N. Polouliakh, K. Oka, F. Nielsen, K. Shibanai and H. Kitano  
Pre-compiled JAR File: https://github.com/agct2018/agct/releases/download/v2.0/AGCT2alpha.jar

Demo Movie (Youtube): https://youtu.be/NoCLfP3t1sM
This tarball contains the following parts: 
* Quickstart pdf presentation to run agct.jar: AGCT_kickstart_vesion2.pdf  
* AGCT paper supplementary information, in an upper directory AGCT_paper_supplementary_information  
* Datasets folder   Example file for YCC data: Yeast_6178.txt (Spellman at al), Yeast_9335.txt (Tu et al), Yeast_3565.txt (periodic genes form Tu et al) and Test_570.txt (small dataset for testing).
GEO tag file: Yeast.txt

In addition, we provide datasets containing the first 570 genes of YMC data, if the reviewer wants to do some fast/faster checks, and see e.g., that the manifold construction is robust to the number of genes as well  (some operations involving matrices will also be 4 to 100 times faster, which may also be good to try various additional algorithms). The scenario in AGCT_Vesion2 is recorded automatically upon start, and the User should only stop it by click on the scenario button and saving the scenario.

* Sentinels folder includes three groups of sentinels for Yeast Metabolic Dataset (probes that are known to be significantly correlated, i.e., biologically meaningful for the data). A user can load or ignore loading them when AGCT shows pop up window for their loading. 
“Sentinels” here are genes/probes "a priori" verified as biologically important and whose locations/interaction to other genes might bring new clues/insights into a biological discussion.
** To run agct, we strongly suggest to use a 64 bits JVM with sufficient memory, for example:  
java -d64 -Xmx8000m -jar agct.jar  (change jar file name here please)

* Scenario Example of small test set scenario is provided, however instead of modifying path before running we suggest directly run small dataset to the User personal computer and AGCT automatically generate the path.
** You must edit Line 11 in /Version_2/Scenarii/Test.txt to give a correct path.
** To test the scenario, simply click on the “Scenario” folder in AGCT (“load and execute scenario”) and select the scenario file in its folder. The scenario will then run automatically until its end (in this case, the computation of the non-linear manifold). 

* To run AGCT for yourself and replicate the scenario, check file AGCT_kickstart_vesion2.pdf 

* To compile the java files and generate the agct.jar file, (i) go in the Source repository and run ./compile_all.sh (compiles all java files), then (ii) run ./makejar.sh, which will generate the jar file and place it at the right place to run agct.jar as above.
** You must compile with Java Development Kit 1.8 https://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html 

Thank you.



