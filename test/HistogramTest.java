package org.opensourcephysics.sip.CPM.test;

import static org.junit.Assert.*;

import org.junit.Test;
import org.opensourcephysics.sip.CPM.Histogram;

public class HistogramTest {
	double [] testData = { 4.7, 3.6, 2.4, 1.23};


	@Test
	public void testHistogram() {
		Histogram h = new Histogram(testData, 1, 5, 4);
	}
	
	@Test
	public void testGetBin(){
		Histogram h = new Histogram(testData, 1, 5, 4);
		assertEquals(h.getBins(), 4);
	}

	@Test
	public void testGetX() {
		Histogram h = new Histogram(testData, 1, 5, 4);
		assertTrue(h.getX(0) == 1.5);
		assertTrue(h.getX(1) == 2.5);
		assertTrue(h.getX(2) == 3.5);
		assertTrue(h.getX(3) == 4.5);
	}

	@Test
	public void testGetY() {
		Histogram h = new Histogram(testData, 1, 5, 4);
		for(int i = 0; i < h.getBins(); i++){
			assertTrue(h.getY(i) == 0.25);
		}	
	}
	
	@Test 
	public void testGetBinWidth(){
		Histogram h = new Histogram(testData, 1, 5, 4);
		assertTrue(h.getBinWidth() == 1);
	}
	
	@Test
	public void testHasNext(){
		Histogram h = new Histogram(testData, 1, 5, 4);
		assertTrue(h.hasNext());
	}
	
	@Test
	public void testNext(){
		Histogram h = new Histogram(testData, 1, 5, 4);
		assertTrue(h.next().equals(new Histogram.Point(1.5, 0.25)));
		assertTrue(h.next().equals(new Histogram.Point(2.5, 0.25)));
		assertTrue(h.next().equals(new Histogram.Point(3.5, 0.25)));
		assertTrue(h.next().equals(new Histogram.Point(4.5, 0.25)));
	}
	
	@Test
	public void testCalculateEntropy(){
		Histogram h = new Histogram(testData, 1, 5, 4);
		
		double expectedEntropy = Math.log(0.25); // 4 * 0.25 * log( 0.25 / 1)
		assertTrue((h.calculateEntropy() - expectedEntropy) < 0.001);
	}
	
	@Test
	public void testIterator_success(){
		Histogram h = new Histogram(testData, 1, 5, 4);
		
		int index = 0;
		for(Histogram.Point p : h){
			assertTrue(p.y == 0.25);
			switch(index){
			case 0: assertTrue(p.x == 1.5); break;
			case 1: assertTrue(p.x == 2.5); break;
			case 2: assertTrue(p.x == 3.5); break;
			case 3: assertTrue(p.x == 4.5); break;
			}
			
			index++;
		}
	}
}
