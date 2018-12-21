import java.io.*;
import java.util.*;
import javax.swing.*;
import java.text.*;

 class AGCTRessourceFileManager {

    public static String
	DEFAULT_PATH_DATA_FILES = "",
	DEFAULT_PATH_HIGHLIGHT_FILES = "",
	DEFAULT_PATH_CLUSTERING_FILES = "",
	DEFAULT_PATH_CHI2_FILES = "",
	DEFAULT_PATH_SCENARIO_FILES = "";

    public static String 
	COMMENT = "//",
	PATH_DATA_FILES = "@PATH_DATA_FILES",
	PATH_HIGHLIGHT_FILES = "@PATH_HIGHLIGHT_FILES",
	PATH_CLUSTERING_FILES = "@PATH_CLUSTERING_FILES",
	PATH_CHI2_FILES = "@PATH_CHI2_FILES",
	PATH_SCENARIO_FILES = "@PATH_SCENARIO_FILES";
	
    String Path_Data_Files,
	Path_Highlight_Files,
	Path_Clustering_Files,
	Path_Chi2_Files,
	Path_Scenario_Files;

    public void loadRessource(){
	load(new File("./Resources/AccessFiles.txt"));
    }

    public void toDefault(){
	Path_Data_Files = AGCTRessourceFileManager.DEFAULT_PATH_DATA_FILES;
	Path_Highlight_Files = AGCTRessourceFileManager.DEFAULT_PATH_HIGHLIGHT_FILES;
	Path_Clustering_Files = AGCTRessourceFileManager.DEFAULT_PATH_CLUSTERING_FILES;
	Path_Chi2_Files = AGCTRessourceFileManager.DEFAULT_PATH_CHI2_FILES;
	Path_Scenario_Files = AGCTRessourceFileManager.DEFAULT_PATH_SCENARIO_FILES;
    }

    public void load(File fl){
	FileReader e = null;
	BufferedReader br;
	String s, nf;
	StringTokenizer t;
	int nt;
	boolean zFound_PATH_DATA_FILES = false, zFound_PATH_HIGHLIGHT_FILES = false, zFound_PATH_CLUSTERING_FILES = false, zFound_PATH_CHI2_FILES = false, zFound_PATH_SCENARIO_FILES = false, zNoFile = false;

	try{
	    e = new FileReader(fl);
	}catch(FileNotFoundException ex){
	    toDefault();
	    zNoFile = true;
	}
	if (!zNoFile){
	    br = new BufferedReader(e);
	    try{
		while ( (s=br.readLine()) != null){
		    t = new StringTokenizer(s.replace('\t',';'), ";");
		    nt = t.countTokens();
		    if (nt>0) {
			s = t.nextToken();
			if ( (s.length()>1) && (!s.substring(0,AGCTRessourceFileManager.COMMENT.length()).equals(AGCTRessourceFileManager.COMMENT)) ){

			    if (s.equals(AGCTRessourceFileManager.PATH_DATA_FILES)){
				zFound_PATH_DATA_FILES = true;
				Path_Data_Files = t.nextToken();
			    }
			    if (s.equals(AGCTRessourceFileManager.PATH_HIGHLIGHT_FILES)){
				zFound_PATH_HIGHLIGHT_FILES = true;
				Path_Highlight_Files = t.nextToken();
			    }
			    if (s.equals(AGCTRessourceFileManager.PATH_CLUSTERING_FILES)){
				zFound_PATH_CLUSTERING_FILES = true;
				Path_Clustering_Files = t.nextToken();
			    }
			    if (s.equals(AGCTRessourceFileManager.PATH_CHI2_FILES)){
				zFound_PATH_CHI2_FILES = true;
				Path_Chi2_Files = t.nextToken();
			    }
			    if (s.equals(AGCTRessourceFileManager.PATH_SCENARIO_FILES)){
				zFound_PATH_SCENARIO_FILES = true;
				Path_Scenario_Files = t.nextToken();
			    }
			}
		    }
		}
		
		if (!zFound_PATH_DATA_FILES)
		    Path_Data_Files = AGCTRessourceFileManager.DEFAULT_PATH_DATA_FILES;
		if (!zFound_PATH_HIGHLIGHT_FILES)
		    Path_Highlight_Files = AGCTRessourceFileManager.DEFAULT_PATH_HIGHLIGHT_FILES;
		if (!zFound_PATH_CLUSTERING_FILES)
		    Path_Clustering_Files = AGCTRessourceFileManager.DEFAULT_PATH_CLUSTERING_FILES;
		if (!zFound_PATH_CHI2_FILES)
		    Path_Chi2_Files = AGCTRessourceFileManager.DEFAULT_PATH_CHI2_FILES;
		if (!zFound_PATH_SCENARIO_FILES)
		    Path_Scenario_Files = AGCTRessourceFileManager.DEFAULT_PATH_SCENARIO_FILES;

	    }catch(IOException ex){
		System.out.println(fl.getName() + " --- IOError when getting ressource !");
		return;
	    }
	}
    }
}


