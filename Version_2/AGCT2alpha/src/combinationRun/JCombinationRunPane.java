//package combinationRun;
//
//import java.util.ArrayList;
//import java.util.Collections;
//import java.util.HashMap;
//import java.util.Map.Entry;
//
//import gene.GeneList;
//
//import javax.swing.JFrame;
//
//public class JCombinationRunPane extends JFrame{
//	private GeneList geneList;
//	public void setGeneList(GeneList geneList){
//		this.geneList = geneList;
//	}
//	public JCombinationRunPane(String caption){
//		super(caption);
//	}
//	private ArrayList<State> stateList = new ArrayList<State>();
//	
//	public void calc(){
//		HashMap<FilterData, GeneList> map = geneList.calcFilteredGeneListIterative();
//		for(Entry<FilterData,GeneList> entry:map.entrySet()){
//			stateList.add(new State(entry.getKey(),entry.getValue()));
//		}
//		Collections.sort(stateList);
//		for(State state:stateList){
//			System.out.println(state.toString());
//		}
//	}
//	
//	
//	
//	private class State implements Comparable<State>{
//		private GeneList geneList;
//		private FilterData filterData;
//		State(FilterData filterData,GeneList geneList){
//			this.geneList = geneList;
//			this.filterData = filterData;
//		}
//		@Override
//		public int compareTo(State o) {
//			return -(int)Math.signum(geneList.getPercentageOfCoveringResponsiveGenes() - o.geneList.getPercentageOfCoveringResponsiveGenes());
//		}
//		
//		@Override
//		public String toString() {
//			// TODO Auto-generated method stub
//			return Integer.toBinaryString(filterData.getBitMask()) + " " + filterData.getPercentile()*100 + " " + geneList.getPercentageOfCoveringResponsiveGenes();
//		}
//	}
//}
