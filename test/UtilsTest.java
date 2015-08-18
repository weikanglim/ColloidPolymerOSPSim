package org.opensourcephysics.sip.CPM.test;

import static org.junit.Assert.*;

import java.util.ArrayList;

import org.junit.Test;
import org.opensourcephysics.sip.CPM.Utils;

public class UtilsTest {

	@Test
	public void testToPrimitiveArray() {
		ArrayList<Double> test = new ArrayList<Double>();
		test.add(1.2);
		test.add(5.0);
		test.add(6.0);
		test.add(7.7);
		
		double [] expected = {1.2, 5.0 , 6.0, 7.7};
		
		double [] results = Utils.toPrimitiveArray(test);
		
		assertTrue(expected.length == results.length);
		for(int i=0; i < results.length; i++){
			assertTrue(expected[i] == results[i]);
		}
	}

}
