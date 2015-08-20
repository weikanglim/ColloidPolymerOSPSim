package org.opensourcephysics.sip.CPM;

import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * This class provides an implementation of a histogram.
 * 
 *  
 * @author Wei Kang Lim
 *
 */
public class Histogram implements Iterator<Histogram.Point>, Iterable<Histogram.Point>{
	private double [] histogramY; // frequency of each bin
	private double [] normalizedY; // normalized frequency (probability of each bin)
	private double [] histogramX; // x-value of each bin
	private int bins;
	private int iteratorCount = 0;
	private double binWidth;
	private boolean normalized = false;
	
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
	 * Creates an empty histogram.
	 * @param start The start range of the histogram.
	 * @param end The end range of the histogram.
	 * @param bins The number of bins in the histogram.
	 */
	public Histogram(double start, double end, int bins){
		histogramY = new double[bins];
		histogramX = new double[bins];
		normalizedY = new double[bins];
		this.bins = bins;
		this.binWidth = (end - start) / (bins);
		
		for(int i=0; i < bins; i++){
			histogramY[i] = 0;
			histogramX[i] = start + binWidth*i + binWidth/2; 
		}
	}
	
	/**
	 * 
	 * Creates a histogram from a random assortment of data.
	 * @param data The dataset.
	 * @param start The start range of the histogram.
	 * @param end The end range of the histogram.
	 * @param bins The number of bins in the histogram.
	 */
	public Histogram(List<Double> data, double start, double end, int bins){
		this(start, end, bins);
		Collections.sort(data); // ascending order
		
		// Perform bin counting.
		int count = 0;
		int j = 0; // data pointer
		
		for(int i = 0; i < bins ; i++){
			count = 0;
			
			while(j < data.size() && data.get(j) <= endRangeOfBin(i)){
				count++;
				j++;
			}
			
			histogramY[i] = count;
		}
		
		normalize();
	}
	
	/**
	 * Creates a (discretely normalzed) histogram from a random assortment of data.
	 * @param data The dataset.
	 * @param start The start range of the histogram.
	 * @param end The end range of the histogram.
	 * @param bins The number of bins in the histogram.
	 */
	public Histogram(double [] data, double start, double end, int bins){
		this(start, end, bins);
		Arrays.sort(data); // ascending order
		
		// Perform bin counting.
		int count = 0;
		int j = 0; // data pointer
		
		for(int i = 0; i < bins ; i++){
			count = 0;
			
			while(j < data.length && data[j] <= endRangeOfBin(i)){
				count++;
				j++;
			}
			
			histogramY[i] = count;
			histogramX[i] = start + binWidth*i + binWidth/2; 
		}
		
		normalize();
	}
	
	public void add(double newValue){
		int low = 0;
		int high = this.bins;
		int mid;
		
		while(low <= high){
			mid = (low + high)/2;
			if(mid < 0 || mid >= this.bins) break;	
			
			if(insideBin(newValue, mid)){
				histogramY[mid]++;
				break;
			} else if(newValue <= startRangeOfBin(mid)){
				high = mid-1;
			} else{
				low = mid+1;
			}
		}
		
		normalized = false;
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
	 * Returns the y-value (frequency) of the given bin.
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
		if(!normalized) normalize();
		
		double entropy = 0;
		for(int i=0; i < this.getBins(); i++){
			if(normalizedY[i] != 0){
				entropy += normalizedY[i] * Math.log(normalizedY[i] / this.getBinWidth());
			}
		}
		
		return entropy;
	}
	
	/**
	 * Removes all frequency data from the histogram.
	 */
	public void clear(){
		this.normalized = false;
		this.iteratorCount = 0;

		for(int i=0; i < this.getBins(); i++){
			histogramY[i] = 0;
			normalizedY[i] = 0;
		}
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
	public boolean equals(Object o){
		if(o instanceof Histogram){
			Histogram h = (Histogram) o;
			if(h.getBins() != this.getBins() || h.getBinWidth() != this.getBinWidth()) return false;
			
			for(int i=0; i < this.getBins(); i++){
				if(this.getX(i) != h.getX(i) || this.getY(i) != h.getY(i)) return false;
			}
			
			return true;
		}
		
		return false;
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
	
	/**
	 * Performs discrete normalization of the histogram distribution. 
	 * 
	 * The normalized data is stored in normalizedY.
	 */
	private void normalize(){
		// Perform discrete normalization
		int totalCount = 0;
		System.arraycopy(histogramY, 0, normalizedY, 0, histogramY.length);
		
		for(int i = 0; i < bins ;i++){
			totalCount += normalizedY[i];
		}
		
		for(int i = 0; i < bins; i++){
			normalizedY[i] /= totalCount;
		}
		
		normalized = true;
	}
	
	private boolean insideBin(double value, int index){
		return (value > startRangeOfBin(index) && value <= endRangeOfBin(index));
	}
	
	private double startRangeOfBin(int index){
		return (histogramX[index] - binWidth / 2);
	}
	
	private double endRangeOfBin(int index){
		return (histogramX[index] + binWidth / 2);
	}
}
