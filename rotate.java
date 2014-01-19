	/**
	 * Attempts a trial rotation of a polymer. The rotation magnitude is specified by the instance variable, rotMagnitude.
	 * @param poly The polymer to be rotated.
	 */
	public void rotate(Polymer poly){
		double [] a = poly.getAxis();
		double [] v = new double[3];

		// generate randomly oriented vector v 
		for(int i = 0 ; i < v.length; i++){
			v[i] = Math.random() - 0.5;
		}

		// normalize new (randomly oriented) vector v 

		VectorMath.normalize(v);

		/* Alternative way to generate v:

		v[2] = 2.*(Math.random() - 0.5);
		double sintheta = Math.signum(v[2])*Math.sqrt(1.-v[2]*v[2]);
		double phi = 2.*(Math.random() - 0.5)*Math.PI;
		v[0] = sintheta*Math.cos(phi);
		v[1] = sintheta*Math.sin(phi);

		*/
				
		// addition of the old and new vector 
		// Note: rotTolerance, which replaces rotMagnitude, should be << 1 (e.g. 0.1)
		for(int i = 0; i < v.length; i++){
			a[i] = a[i] + rotTolerance*v[i];
		}

		// normalize result
		VectorMath.normalize(a);
		poly.setTransformAxis(a);
		
		poly.setAxis(a);  // added by Alan (needed to update vector a)

// No further changes by Alan from here on ...

		// Check for overlaps
		Stack<Nano> overlapNanos = new Stack<Nano>();
		int overlapCount = 0;

		// Check for intersections with nanoparticles
		for (int i = 0; i < nanos.length; i++) {
			if (poly.overlap(nanos[i]) 
					&& !poly.intersectPairs.contains(nanos[i])
						&& !nanos[i].intersectPairs.contains(poly)) { // Check for previous overlap
					overlapCount++;
					overlapNanos.push(nanos[i]);
			}
		}
		
		// Acceptance probability
		if (Math.random() < Math.exp(-Ep*overlapCount)) {		
		 // Update the intersecting pairs
			while(!overlapNanos.empty()){
				overlapNanos.peek().intersectPairs.add(poly);
				poly.intersectPairs.add(overlapNanos.pop());
				totalIntersectCount++;
			}
		
			// Remove particles that are no longer overlapping
			for (int j = 0; j < nanos.length; j++) {
				if (!poly.overlap(nanos[j]) &&
						poly.intersectPairs.remove(nanos[j]) && 
							nanos[j].intersectPairs.remove(poly)){ 
						totalIntersectCount--;
				}
			}
		} else {
			poly.setTransformAxis(poly.getAxis());
		}

	}
}
