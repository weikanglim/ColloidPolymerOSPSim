package org.opensourcephysics.sip.CPM;

import java.util.List;

public class Utils {
	public static double [] toPrimitiveArray(List<Double> d){
		if(d == null) return null;
		
		double [] arr = new double[d.size()];
		
		for(int i = 0; i < arr.length; i++){
			arr[i] = d.get(i);
		}
		
		return arr;
	}
	
	public static void averageRuns(double [][][] runData, double [][] results){
		int runs = runData.length;
		int dataPoints = results.length;
		int numVariables = runData[0][0].length;
		
		for(int dataPoint = 0; dataPoint < dataPoints; dataPoint++){
			results[dataPoint][0] = runData[0][dataPoint][0];
		}
		
		int resultsIndex = 1;
		for(int var = 1; var < numVariables; var++){
			for(int dataPoint = 0; dataPoint < dataPoints; dataPoint++){
				for(int run = 0; run < runs; run++){
					results[dataPoint][resultsIndex] += runData[run][dataPoint][var];
					results[dataPoint][resultsIndex+1] += Math.pow(runData[run][dataPoint][var],2);
				}
				
				
				results[dataPoint][resultsIndex] /= runs;
				results[dataPoint][resultsIndex+1] /= runs;
				results[dataPoint][resultsIndex+1] = Math.sqrt(results[dataPoint][resultsIndex+1] - Math.pow(results[dataPoint][resultsIndex],2));
			}
			
			resultsIndex += 2;
		}
	}
}
