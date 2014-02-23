package org.opensourcephysics.sip.CPM;

import java.util.Stack;


/**
 * NanoPolyMix is an abstraction for a binary mixture of colloids(nanoparticles)
 * and polymers model.
 * 
 * @author Wei Kang Lim, Alan Denton
 * @version 1.0 11-10-2013 * 
 */
public class CPM {
	public enum Vector {x, y, z};
	
	// constants
	public final double K = 0.015923;
	public final double a = 0.0802;
	public final double D = 1.842;
	// end constants

	public Polymer polymers[];
	public Nano nanos[];
	public int nN, nP, nx; // number of nanoparticles, number of polymers,
							// number of columns and rows
	// public double rho; // Ratio of polymer diameter to nanoparticle (colloid)
	// diameter, sigP/sigN
	public double init_U;
	public double sigmaX;
	public double sigmaY;
	public double q;
	public int trialShapeChangePerMcs;
	public int trialDisplacementPerMcs;
	public int moveToShapeRatio;
	public double lc; // lattice constant, defined as lc = d / sigN, where d is
						// defined below.
	public double Lx; // dimension of the boundary (along x)
	public double Ly; // dimension of the boundary (along y)
	public double Lz; // dimension of the boundary (along z)
	public double steps; // number of monte carlo steps
	public double tolerance;
	public double shapeTolerance;
	public double d; // the center-to-center distance between particles of same
						// species (hexagonal lattice),
						// or columns and rows (square lattice)
	public double totalIntersectCount;
	public double Ep; // Penetration Energy
	public double mcs;
//	public double rotTolerance;

	// end declaration


	/**
	 * initializes the model.
	 * 
	 * @param configuration
	 *            the initial lattice structure of the model
	 */
	public void initialize(String configuration) {
		// set-up variables
		totalIntersectCount = 0;
		mcs = 0;
		steps = 0;
		polymers = new Polymer[nP];
		nanos = new Nano[nN];
		d = lc; // distance between two nanoparticles
		Ep = 3/q;
		Nano.setTolerance(tolerance);
		Polymer.setTolerance(tolerance);
		Polymer.setQ(q);
		Nano.setDefault_r(0.5);
		Polymer.setDefault_U(init_U);

		// initialize positions
		if (configuration.toUpperCase().equals("SQUARE")) {
			setSqrPositions();
		} 
	}

	/**
	 * Places particles on a square lattice depending on the number of nano particles.
	 * If number of nano particles is 0, polymers are placed by themselves.
	 */
	public void setSqrPositions() {
		int ix, iy, iz;
		if (nN > 0) {
			double dnx = Math.cbrt(nN);
			nx = (int) dnx;
			if (dnx - nx > 0.00001) {
				nx++; // N is not a perfect cube
			}
			Lx = Ly = Lz = d * nx;
			Particle.setBoundaries(Lx, Ly, Lz);

			// place nano particles
			int i = 0;
			for (iy = 0; iy < nx; iy++) { // loops through particles in a column
				for (ix = 0; ix < nx; ix++) { // loops through particles in a
												// row
					for (iz = 0; iz < nx; iz++) {
						if (i < nanos.length) { // checks for remaining
												// particles
							nanos[i] = new Nano(ix * d, iy * d, iz * d);
							i++;
						}
					}
				}
			}

			// place polymer particles
			i = 0;
			for (iy = 0; iy < nx; iy++) {
				for (ix = 0; ix < nx; ix++) {
					for (iz = 0; iz < nx; iz++) {
						if (i < polymers.length) { // checks for remaining
													// particles
							polymers[i] = new Polymer(ix * d + d / 2, iy * d
									+ d / 2, iz * d + d / 2);
							i++;

							if ((ix + 1) * d + d / 2 > Lx
									&& (iy + 1) * d + d / 2 > Ly
									&& (iz + 1) * d + d / 2 > Lz) { // the next
																	// position
																	// is out of
																	// the
																	// lattice
								ix = iy = 0; // reset position counter,
												// overlap
												// poylmers
								iz = -1;
							}
						} else
							return;
					}
				}
			}
		} else { // no nanoparticles on the lattice
			double dnx = Math.cbrt(nP); // count lattice points based on number of polymers instead
			nx = (int) dnx;
			if (dnx - nx > 0.00001) {
				nx++; // N is not a perfect cube
			}
			Lx = Ly = Lz = d * nx;
			Particle.setBoundaries(Lx, Ly, Lz);

			// place polymers
			int i = 0;
			for (iy = 0; iy < nx; iy++) {
				for (ix = 0; ix < nx; ix++) {
					for (iz = 0; iz < nx; iz++) {
						if (i < polymers.length) { // checks for remaining
													// particles
							polymers[i] = new Polymer(ix * d + d / 2, iy * d
									+ d / 2, iz * d + d / 2);
							i++;

							if ((ix + 1) * d + d / 2 > Lx
									&& (iy + 1) * d + d / 2 > Ly
									&& (iz + 1) * d + d / 2 > Lz) { // the next
																	// position
																	// is out of
																	// the
																	// lattice
								ix = iy = 0; // reset position counter,
												// overlap
												// poylmers
								iz = -1;
							}
						} else
							return;
					}
				}
			}
		}
	}

	/**
	 * Does a monte carlo simulation step.
	 */
	public void step() {
		mcs++;
		trialMoves();
	}

	/**
	 * Attempts trial moves, rotations, and shape changes.
	 */
	public void trialMoves() {
		for(int x = 1; x <= trialDisplacementPerMcs; x++){
			for (int i = 0; i < nanos.length; ++i) {
				nanoTrialMove(nanos[i]);
			}
			// Polymer Trial Displacements
			for (int i = 0; i < polymers.length; ++i) {
				polyTrialMove(polymers[i]);
			}
		}
				
		// Polymer Trial Shape Changes
		for(int x = 1; x <= trialShapeChangePerMcs; x++){
			for (int i = 0; i < polymers.length; ++i) {
					shapeChange(polymers[i]);
			}
		}
	}

	/**
	 * Trial move for a polymer.
	 * 
	 * @param poly
	 *            A polymer object to be moved.
	 */
	public void polyTrialMove(Polymer poly) {
		double oldX = poly.getX();
		double oldY = poly.getY();
		double oldZ = poly.getZ();
		int overlapCount = 0;
		Stack<Nano> overlapNanos = new Stack<Nano>();
		poly.move();

		// Check for intersections with nanoparticles
		for (int i = 0; i < nanos.length; i++) {
			boolean overlap = poly.overlap(nanos[i]);
			boolean wasOverlapping = poly.intersectPairs.contains(nanos[i]) || 
					nanos[i].intersectPairs.contains(poly);
			if (overlap && !wasOverlapping) { // Check for gain in overlap
					overlapCount++;
					overlapNanos.push(nanos[i]);
			}
			else if(!overlap && wasOverlapping) { // Check for loss of previous overlap
				overlapCount--;
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
			poly.setX(oldX);
			poly.setY(oldY);
			poly.setZ(oldZ);
		}
	}

	/**
	 * Nanoparticle trial moves
	 * 
	 * @param nano
	 *            A nanoparticle to be moved.
	 */
	public void nanoTrialMove(Nano nano) {
		double oldX = nano.getX();
		double oldY = nano.getY();
		double oldZ = nano.getZ();
		int overlapCount = 0;
		Stack<Polymer> overlapPolymers = new Stack<Polymer>();
		nano.move();

		// Nano-nano intersections, reject immediately
		for (int j = 0; j < nanos.length; ++j) {
			if (!nano.equals(nanos[j]) && nano.overlap(nanos[j])) { // Nano-nano overlap
					nano.setX(oldX);
					nano.setY(oldY);
					nano.setZ(oldZ);
					return;
				}
		}

		// Count number of intersections
		for (int i = 0; i < polymers.length; i++) {
			boolean overlap = nano.overlap(polymers[i]);
			boolean wasOverlapping = polymers[i].intersectPairs.contains(nano) || 
					nano.intersectPairs.contains(polymers[i]);
			if (overlap && !wasOverlapping) { // Check for gain in overlap
					overlapCount++;
					overlapPolymers.push(polymers[i]);
			}
			else if(!overlap && wasOverlapping) { // Check for loss of previous overlap
				overlapCount--;
			}
		}
		
		// Acceptance probability
		if (Math.random() < Math.exp(-Ep*overlapCount)) { 
				while(!overlapPolymers.isEmpty()){
					overlapPolymers.peek().intersectPairs.add(nano);
					nano.intersectPairs.add(overlapPolymers.pop()); // Add
					totalIntersectCount++;
				}
				
				// Remove particles that are no longer overlapping
				for (int j = 0; j < polymers.length; j++) {
					if (!nano.overlap(polymers[j]) &&
							nano.intersectPairs.remove(polymers[j])
								&& polymers[j].intersectPairs.remove(nano)){
							totalIntersectCount--;
					}
				}
		}
		else 
		{
			nano.setX(oldX);
			nano.setY(oldY);
			nano.setZ(oldZ);
		}
	}

	/**
	 * Shape changes for a polymer.
	 * 
	 * @param poly
	 *            A polymer object to be attempted for a trial shape change.
	 */
	public void shapeChange(Polymer poly) {
		int overlapCount = 0;
		Stack<Nano> overlapNanos = new Stack<Nano>();
		
		double oldU = poly.getU();

		// trial shape changes
		double newU = oldU + shapeTolerance * 2. * (Math.random() - 0.5);
		
		// Reject changes for negative radii eigenvalue immediately
		if(newU < 0){
			return;
		}
		
		poly.setU(newU);
		
		// Check for intersections with nanoparticles
		for (int i = 0; i < nanos.length; i++) {
			boolean overlap = poly.overlap(nanos[i]);
			boolean wasOverlapping = poly.intersectPairs.contains(nanos[i]) || 
					nanos[i].intersectPairs.contains(poly);
			if (overlap && !wasOverlapping) { // Check for gain in overlap
					overlapCount++;
					overlapNanos.push(nanos[i]);
			}
			else if(!overlap && wasOverlapping) { // Check for loss of previous overlap
				overlapCount--;
			}
		}

		double p = (prob(newU)  ) / 
				   (prob(oldU)  )  * Math.exp(-Ep*overlapCount);
		
		// Acceptance probability
		if (p > 1 || Math.random() < p) {		
		 // Update the intersecting pairs
			while(!overlapNanos.empty()){
				overlapNanos.peek().intersectPairs.add(poly);
				poly.intersectPairs.add(overlapNanos.pop());
				totalIntersectCount++;
			}
		
			// Since shape change was accepted, update possible intersections that were removed as a result.
			for (int j = 0; j < nanos.length; j++) {
				if (!poly.overlap(nanos[j])) { // particles that are no longer
												// overlapping
					if (poly.intersectPairs.remove(nanos[j])
							&& nanos[j].intersectPairs.remove(poly)) // update the
																		// intersecting
																		// pairs and
																		// count
						totalIntersectCount--;
				}
			}
			
		} else {
			poly.setU(oldU);
		}
	}

	/**
	 * Returns the polymer shape probability for a given radius axis.
	 * 
	 * @param ei
	 *            The eigenvalue of the radius axis
	 * @return Polymer shape probability.
	 */
	public double prob(double u) {
		return 1 / (2 * K * u) * Math.exp(-u/a - D*D*a/u);
	}

	
	/**
	 * Calculates the polymer colloid size ratio.
	 * 
	 * @param v
	 *            enumeration representing the radius axis
	 * @return Polymer colloid size ratio.
	 */
	public double polymerColloidSizeRatio() {
		double total = 0;
		double average = 0;
		double ratio = 0;

		for (Polymer poly : polymers) {
			total += poly.getrX() + poly.getrY() + poly.getrZ();
		}

		average = Math.sqrt(total / polymers.length);
		ratio = average / Nano.getDefault_r();
		return ratio;
	}
	
}

