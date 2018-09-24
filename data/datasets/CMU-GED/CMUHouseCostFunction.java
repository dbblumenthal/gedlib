package util;

import java.text.DecimalFormat;
import java.util.Locale;

import algorithms.Constants;


public class CMUHouseCostFunction implements ICostFunction {

	private double nodeCosts;

	private double edgeCosts;
	
	private double alpha;
	
	public CMUHouseCostFunction(double nodeCosts, double edgeCosts, double alpha){
		this.nodeCosts = nodeCosts;
		this.edgeCosts = edgeCosts;
		this.alpha = alpha;
	}
	
	@Override
	public double getCosts(GraphComponent start, GraphComponent end) {
		// TODO Auto-generated method stub
		
		if( Constants.nodecostmatrix!=null && Constants.edgecostmatrix!=null){

			if( Constants.nodecostmatrix!=null && Constants.edgecostmatrix!=null){
				return precomputedcosts(start,end);


			}
		}
		/**
		 * node handling
		 */ 
		if (start.isNode() || end.isNode()) {
			double xStart;
			double yStart;
			String startType;
			double xEnd;
			double yEnd;
			String endType;
			// start is not empty
			if (!start.getComponentId().equals(Constants.EPS_ID)) {
				startType = (String) start.getTable().get("type");
				String xStartString = (String) start.getTable().get("x");
				xStart = Double.parseDouble(xStartString);
				String yStartString = (String) start.getTable().get("y");
				yStart = Double.parseDouble(yStartString);				
			} else {
				// insertion
				return this.alpha * this.nodeCosts;
			}
			// end is not empty
			if (!end.getComponentId().equals(Constants.EPS_ID)) {
				endType = (String) end.getTable().get("type");
				String xEndString = (String) end.getTable().get("x");
				xEnd = Double.parseDouble(xEndString);
				String yEndString = (String) end.getTable().get("y");
				yEnd = Double.parseDouble(yEndString);
			} else {
				// deletion
				return this.alpha * this.nodeCosts;
			}
			
			double distance = Math.sqrt(Math.pow((xEnd - xStart), 2.)
					+ Math.pow((yEnd - yStart), 2.));
			DecimalFormat decFormat = (DecimalFormat) DecimalFormat
					.getInstance(Locale.ENGLISH);
			decFormat.applyPattern("0.00000");
			String distanceString = decFormat.format(distance);
			distance = Double.parseDouble(distanceString);
			return this.alpha * 0;
				

		}/**
		 * edge handling
		 */ 
		else {
			double label1=0;
			double label2=0;;
			String stringlabel1;
			String stringlabel2;
			// start is not empty
			if (!start.getComponentId().equals(Constants.EPS_ID) && (start.isNode()==false)) {
				stringlabel1 = (String) start.getTable().get("dist");
				label1 = Double.parseDouble(stringlabel1);						
			} 
			// end is not empty
			if (!end.getComponentId().equals(Constants.EPS_ID)) {
				stringlabel2 = (String) end.getTable().get("dist");
				label2 = Double.parseDouble(stringlabel2);
			} 
			
			// insertion
			if( (start.getComponentId().equals(Constants.EPS_ID) ==true) ){
				return (1-this.alpha)*label2;
			}
			
			// deletion
			if( (end.getComponentId().equals(Constants.EPS_ID) ==true) ){
				return (1-this.alpha)*label1;
			}
			
			double distance = Math.abs(label1 - label2);//Math.sqrt(Math.pow((label1 - label2), 2.));
			DecimalFormat decFormat = (DecimalFormat) DecimalFormat
					.getInstance(Locale.ENGLISH);
			decFormat.applyPattern("0.00000");
			String distanceString = decFormat.format(distance);
			distance = Double.parseDouble(distanceString);
			return  (1-this.alpha) * distance;

		}
	}
	
	
	private double precomputedcosts(GraphComponent start, GraphComponent end) {
		// TODO Auto-generated method stub
		
		double[][] matrix=null;
		
		if (start.isNode() || end.isNode()) {
			matrix=Constants.nodecostmatrix;
		}else{
			matrix=Constants.edgecostmatrix;
		}
		
		int n1 = matrix.length;
		int n2 = matrix[0].length;
		int insertindexg1=n2-2;
		int insertindexg2=n1-2;
		int delindexg1=n2-1;
		int delindexg2=n1-1;


		if (start.getComponentId().equals(Constants.EPS_ID)) {
			if(end.belongtosourcegraph){
				return matrix[end.id][insertindexg1];
			}else{
				return matrix[insertindexg2][end.id];
			}
		}

		if (end.getComponentId().equals(Constants.EPS_ID)) {
			if(start.belongtosourcegraph){
				return matrix[start.id][delindexg1];
			}else{
				return matrix[delindexg2][start.id];
			}
		}
		
		if(start.belongtosourcegraph){
			return matrix[start.id][end.id];
		}else{
			return matrix[end.id][start.id];
		}




	}

	/**
	 * @return the cost of an edge operation
	 */
	public double getEdgeCosts() {
		return (1-alpha)*edgeCosts;
	}
	
	/**
	 * @return the cost of a node operation
	 */
	public double getNodeCosts() {
		return this.alpha * nodeCosts;
	}

}
