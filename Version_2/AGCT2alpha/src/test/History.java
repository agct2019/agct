package test;

class History implements Debuggable {

    public static String[] DATES = {"05/05/2008",
            "05/09/2008",
            "05/10/2008",
            "05/11/2008",
            "05/26/2008",
            "07/04/2008",
            "07/07/2008",
            "07/09/2008",
            "07/25/2008",
            "07/26/2008",
            "07/27/2008",
            "08/05/2008",
            "08/08/2008",
            "08/11/2008",
            "04/06/2009",
            "04/08/2009",
            "04/13/2009",
            "04/15/2009",
            "04/28/2009",
            "04/29/2009",
            "05/03/2009",
            "06/01/2009",
            "06/15/2009",
            "15.12.2013"};

    public static String[] VERSION = {"1.0.0",
            "1.0.1",
            "1.0.2",
            "1.0.3",
            "1.0.4",
            "1.0.5",
            "1.0.6",
            "1.0.7",
            "1.0.8",
            "1.0.9",
            "1.0.10",
            "1.0.11",
            "1.0.12",
            "1.0.13",
            "1.0.14",
            "1.0.15",
            "1.0.16",
            "1.0.17",
            "1.0.18",
            "1.0.19",
            "1.0.20",
            "1.0.21",
            "1.0.22",
            "2.0"};

    public static String[] SUMMARY_OF_CHANGES = {"Bugs corrected; added search facility for clustering.",
            "Added -Best and -BestRatio to search facility for clusterings.",
            "Bugs corrected. Added Scenarii.",
            "Bugs corrected for clustering vs scenarii. Added features to Scenarii.",
            "Bugs corrected when there are <3 features.",
            "Added Wavelets with user-tunable number of time stamps.",
            "First feature selection procedure. Faster methods.",
            "Modified feature selection and selection pane. Bugs corrected.",
            "First prototype selection methods. Faster methods. Bugs corrected.",
            "Bugs corrected.",
            "Added saving possibility for prototype selection (keep).",
            "Faster. Bugs corrected.",
            "Saving features opened and improved.",
            "Bugs corrected in saving.",
            "Minor bugs corrected. Added bound for iterative clustering algorithms. Added Absolute Cosine Correlation. Added Iterative Bregman KMeans.",
            "Added Non Negative Matrix Factorization with L22 and KL divergences.",
            "Added Compartment(C:). Added hierarchy of F:/P:/C:.",
            "Added the P/F/C representation for clusters when All_Chi2 is performed",
            "Faster. Added Benjamini-Hochberg tests for tags. Save Delaunay / Chi2. Save Manifold in Scenario.",
            "Added consistency check for triangulation",
            "Added saving and loading of clusterings (KM)",
            "Added highlighting of genes, join, etc.",
            "Bypassed BugID:6264013"};

    public static String getLatestVersion() {
        return "Version " + VERSION[VERSION.length - 1];
    }

    public static String getLatestDate() {
        return DATES[DATES.length - 1];
    }
}