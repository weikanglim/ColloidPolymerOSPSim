package org.opensourcephysics.sip.CPM;

import java.util.Arrays;
import java.util.Iterator;
import java.lang.Iterable;
import java.util.NoSuchElementException;

/**
 * This class provides an implementation of a discretely normalized histogram.
 * 
 * It is created by reading in a dataset. The Histogram is strictly read-only after instantiating with dataset.
 *  
 * @author Wei Kang Lim
 *
 */
public class Histogram implements Iterator<Histogram.Point>, Iterable<Histogram.Point>{
	private double [] histogramY; // frequency of each bin
	private double [] histogramX; // x-value of each bin
	private int bins;
	private int iteratorCount = 0;
	private double binWidth;
	
	/**
	 * Represents a 2-D point in Cartesian space.
	 * @author Wei Kang Lim
	 *
	 */
	public static class Point{
		public double x; public double y;
		
		public Point(double x, double y){
			this.x = x; this.y = y;
		}
		
		public boolean equals(Histogram.Point p){
			return this.x == p.x && this.y == p.y; 
		}
	}
	
	/**
	 * Creates a (discretely normalzed) histogram from a random assortment of data.
	 * @param data The dataset.
	 * @param start The start range of the histogram.
	 * @param end The end range of the histogram.
	 * @param bins The number of bins in the histogram.
	 */
	public Histogram(double [] data, double start, double end, int bins){
		histogramY = new double[bins];
		histogramX = new double[bins];
		this.bins = bins;
		this.binWidth = (end - start) / (bins);
		Arrays.sort(data); // ascending order
		
		// Perform bin counting.
		int count = 0;
		int j = 0; // data pointer
		
		for(int i = 0; i < bins ; i++){
			count = 0;
			
			while(j < data.length && data[j] <= start + binWidth*(i+1)){
				count++;
				j++;
			}
			
			histogramY[i] = count;
			histogramX[i] = start + binWidth*i + binWidth/2; 
		}
		
		// Perform discrete normalization
		int totalCount = 0;
		for(int i = 0; i < bins ;i++){
			totalCount += histogramY[i];
		}
		
		for(int i = 0; i < bins; i++){
			histogramY[i] /= totalCount;
		}
	}
	
	/**
	 * Returns the x-value of the given bin.
	 * @param binIndex
	 * @return
	 */
	public double getX(int binIndex){
		return histogramX[binIndex];
	}
	
	/**
	 * Returns the y-value of the given bin.
	 * @param binIndex
	 * @return
	 */
	public double getY(int binIndex){
		return histogramY[binIndex];
	}
	
	/**
	 * Returns the total number of bins.
	 * @return
	 */
	public int getBins(){
		return this.bins; 
	}
	
	/**
	 * Returns the width of each bin.
	 * @return
	 */
	public double getBinWidth(){
		return this.binWidth;
	}
	
	/**
	 * Calculates the entropy of the histogram.
	 * @return The entropy of the distribution.
	 */
	public double calculateEntropy(){
		double entropy = 0;
		for(int i=0; i < this.getBins(); i++){
			entropy += this.getY(i) * Math.log(this.getY(i) / this.getBinWidth()); 
		}
		
		return entropy;
	}
	
	/**
	 * Returns a string representation of the histogram.
	 */
	public String toString(){
		String s = "";
		StringBuffer sb = new StringBuffer();
		for(Point p : this){
			sb.append(p.x);
			sb.append(",");
			sb.append(p.y);
			sb.append("\n");
		}
		
		s = sb.toString();		
		return s;
	}

	@Override
	public boolean hasNext() {
		if(iteratorCount >= bins){
			return false;
		}
		
		return true;
	}

	@Override
	public Histogram.Point next() {
		if(iteratorCount >= bins){
			throw new NoSuchElementException();
		}
		
		Point p = new Point(histogramX[iteratorCount], histogramY[iteratorCount]);
		iteratorCount++;
		
		return p;
	}

	@Override
	public Iterator<Point> iterator() {
		this.iteratorCount = 0;
		return this;
	}

	@Override
	public void remove() {
		throw new UnsupportedOperationException("Histogram is unmutable!");
	}
}
