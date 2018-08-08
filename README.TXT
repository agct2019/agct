This is the Code and Dataset tarball for AGCT and paper

"AGCT: A Geometric Clustering Tool for robustly unravelling the inner cluster structures of gene expressions"
by by R. Nock, N. Polouliakh, K. Oka, F. Nielsen, K. Shibanai and H. Kitano

This tarball contains 4 parts:

1* AGCT jarfile
  jar file of AGCT classes: agct.jar

* Quickstart pdf presentation to run agct.jar: AGCT_kickstart.pdf

* AGCT paper supplementary information, in directory AGCT_paper_supplementary_information

* Datasets folder
  Example file for YCC data, containing the profiles of the 6178 Genes: Spellman_Cdc28_6178genes.agct
  GEO tag file: Yeast_Spellman_tags.txt

  In addition, we provide datasets containing the first 1000 genes of YCC data, if the reviewer wants to do some fast/faster checks,
  and see e.g., that the manifold construction is robust to the number of genes as well 
  (some operations involving matrices will also be 4 to 100 times faster, which may also be good to try various additional algorithms)

* Example Scenario folder
  Scenario recovering the paper’s results for YCC: NICTA_AGCT_Spellman_CDC28_6178.txt
  Equivalent scenario on the 1000 genes subset: NICTA_AGCT_Spellman_CDC28_1000.txt

* Source file, containing all Java classes needed to compile and run AGCT

=========================================================
THIS SOFTWARE IS PROVIDED WITHOUT ANY WARRANTY WHATSOEVER
USE IT AT YOUR OWN RISKS
=========================================================

** Before running agct, CHANGE the path file of the dataset in the scenario file, which is the FIRST line of the Scenario file

It currently looks like:	Load|/Users/rnock/Submitted_AGCT_Jar/Datasets/Spellman_Cdc28_6178genes.agct
It SHOULD look like:		Load|The_Complete_path_to_access_this_subdirectory_on_your_computer/Submitted_AGCT_Jar/Datasets/Spellman_Cdc28_6178genes.agct

For example, if you have uncompressed the .zip file in /Users/JohnDoe/, then the line should read
Load|Users/JohnDoe/Submitted_AGCT_Jar/Datasets/Spellman_Cdc28_6178genes.agct

After having made this change, you are ready to go !

** To run agct, we strongly suggest to use a 64 bits JVM with sufficient memory, for example:

java -d64 -Xmx8000m -jar agct.jar

** To test the scenario, simply click on the “Scenario” folder in AGCT (“load and execute scenario”) and select the scenario file in its folder. The scenario will then run automatically until its end (in this case, the computation of the non-linear manifold).

** To run AGCT for yourself and replicate the scenario, check file AGCT_kickstart.pdf 

** To compile the java files and generate the agct.jar file, (i) go in
   the Source repository and run ./compile_all.sh (compiles all java
   files), then (ii) run ./makejar.sh, which will generate the jar
   file and place it at the right place to run agct.jar as above.
** You must compile with Java Development Kit 1.6 https://docs.oracle.com/javase/7/docs/webnotes/install/.
   The program does NOT work in Java 1.8 or above (untested in Java 1.7).

Thank you.
