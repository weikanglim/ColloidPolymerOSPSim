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
}
