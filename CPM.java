package org.opensourcephysics.sip.CPM;

import java.util.Stack;

import org.opensourcephysics.numerics.PBC;
import org.opensourcephysics.numerics.VectorMath;

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
	public final double NX = 0.5;
	public final double NY = 2.5;
	public final double NZ = 4.0;

	public final double AX = 0.08065;
	public final double AY = 0.01813;
	public final double AZ = 0.006031;

	public final double DX = 1.096;
	public final double DY = 1.998;
	public final double DZ = 2.684;

	public final double KX = 0.094551;
	public final double KY = 0.0144146;
	public final double KZ = 0.0052767;
	// end constants

	public Polymer polymers[];
	public Nano nanos[];
	public int nN, nP, rows; // number of nanoparticles, number of polymers,
							// number of columns and rows
	public double init_eX;
	public double init_eY;
	public double init_eZ;
	public double sigmaX;
	public double sigmaY;
	public double nano_r;
	public double q; // Rp / Rc
	public int trialMovesPerMcs;
	public double lc; // lattice constant, defined as lc = d / sigN, where d is
						// defined below.
	public double Lx; // dimension of the boundary (along x)
	public double Ly; // dimension of the boundary (along y)
	public double Lz; // dimension of the boundary (along z)
	public double steps; // number of monte carlo steps
	public double volFraction;
	public double tolerance;
	public double shapeTolerance;
	public double d; // the center-to-center distance between particles of same
						// species (hexagonal lattice),
						// or columns and rows (square lattice)
	public double totalIntersectCount;
	public double Ep; // Penetration Energy
	public double mcs;
	public double rotTolerance;

	// end declaration


	/**
	 * initializes the model.
	 * 
	 * @param configuration
	 *            the initial lattice structure of the model
	 */
	public void initialize(String configuration, boolean penetrationEnergy) {
		// set-up variables
		totalIntersectCount = 0;
		mcs = 0;
		steps = 0;
		polymers = new Polymer[nP];
		nanos = new Nano[nN];
		if(penetrationEnergy){
			Ep = 3./q;
		}else{
			Ep = 0;
			System.out.println("Penetration energy turned off.");
		}
		Nano.setTolerance(tolerance);
		Polymer.setTolerance(tolerance);
		Polymer.setRotTolerance(rotTolerance);
		Polymer.setShapeTolerance(shapeTolerance);
		Polymer.setQ(q);
		Nano.setDefault_r(nano_r); // hardcoded default value
		Polymer.setDefault_eX(init_eX);
		Polymer.setDefault_eY(init_eY);
		Polymer.setDefault_eZ(init_eZ);

		// initialize positions
		if (configuration.toUpperCase().equals("SQUARE")) {
			setPositions();
		} 
		Particle.setBoundaries(Lx, Ly, Lz);
		volFraction = nN * (4 * Math.PI / 3) * Math.pow(nano_r / Lx , 3);
	}

	/**
	 * Initializes the configuration.
	 * A nanoparticle is placed at the center, and the polymers are put into an evenly divided lattice.
	 */
	public void setPositions() {
		int ix, iy, iz;
		double dnx = Math.cbrt(nP); // count lattice points based on number of polymers
		rows = (int) dnx;
		if (dnx - rows > 0.00001) {
			rows++; // N is not a perfect cube
		}

		// place nanoparticle at the center of the lattice box
		nanos[0] = new Nano(Lx/2f, Ly/2f, Lz/2f);

		
		// place polymers
		int i = 0;
		for (iy = 0; iy < rows; iy++) {
			for (ix = 0; ix < rows; ix++) {
				for (iz = 0; iz < rows; iz++) {
					if (i < polymers.length) { // checks for remaining
												// particles
						polymers[i] = new Polymer(ix * Lx/rows, iy * Ly/rows, iz * Lz/rows);
						i++;

						if ((ix + 1) * Lx/rows > Lx
								&& (iy + 1) * Ly/rows > Ly
								&& (iz + 1) * Lz/rows > Lz) { // the next
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
		
		for(int x = 1; x <= trialMovesPerMcs; x++){
			// Polymer Trial Displacements
			for (int i = 0; i < polymers.length; ++i) {
				int j = (int) (Math.random()*polymers.length);
				polyTrialMove(polymers[j]);
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
		double oldEX = poly.geteX();
		double oldEY = poly.geteY();
		double oldEZ = poly.geteZ();
		double [] initialOldAxis = poly.getOldAxis();
		double [] initialNewAxis = poly.getNewAxis();


		int overlapCount = 0;
		Stack<Nano> overlapNanos = new Stack<Nano>();
		/** Trial Displacement **/
		poly.move();
		
		/** Shape Changes **/
		if(Polymer.getShapeTolerance() > 0){
			poly.shapeChange();
		}
		double newEX = poly.geteX();
		double newEY = poly.geteY();
		double newEZ = poly.geteZ();
		
		/** Rotation **/
		if(Polymer.getRotTolerance() > 0){
			poly.rotate();
		}
		
		// Check for intersections with nanoparticles
		for (int i = 0; i < nanos.length; i++) {
			boolean overlap = poly.overlap(nanos[i]);
			boolean wasOverlapping = poly.intersectPairs.contains(nanos[i]) || nanos[i].intersectPairs.contains(poly);
			if (overlap && !wasOverlapping) { // Check for gain in overlap
					overlapCount++;
					overlapNanos.push(nanos[i]);
			}
			else if (!overlap && wasOverlapping) { // Check for loss of previous overlap
					overlapCount--;
			}
		}
		
		double p;
		if(Polymer.getShapeTolerance() == 0){
			p = Math.exp(-Ep*overlapCount);
		} else {
			p = (prob(newEX, Vector.x) * prob(newEY, Vector.y) * prob(newEZ, Vector.z) ) / 
				   (prob(oldEX, Vector.x) * prob(oldEY, Vector.y) * prob(oldEZ, Vector.z) )  * Math.exp(-Ep*overlapCount);
		}

		// Acceptance probability
		if (p > 1 || Math.random() < p) {		
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

			poly.setOldAxis(initialOldAxis);
			poly.setNewAxis(initialNewAxis);
			
			poly.seteX(oldEX);
			poly.seteY(oldEY);
			poly.seteZ(oldEZ);
		}
	}

	/**
	 * Performs a trial placement of nanoparticle at radius r.
	 * 
	 * @param r
	 *            The radial distance to perform the trial placement
	 * @return Probability of acceptance, e^(-delta_U)
	 */
	public double nanoTrialPlacement(double r) {
		
		// Generate random angles in spherical coordinates
		double phi =  2*Math.random()*Math.PI;
		double cosTheta = 2*Math.random() - 1; // (-1,1)
		double theta = Math.acos(cosTheta);
		
		// Calculate cartesian coordinates
		double x = r*Math.cos(phi)*Math.sin(theta);  
		double y = r*Math.sin(phi)*Math.sin(theta);
		double z = r*Math.cos(theta);
		
		// Translate to lattice coordinates
		x += Lx/2f;
		y += Ly/2f;
		z += Lz/2f;
		
		// Count number of intersections
		Nano nano = new Nano(x,y,z);
		int overlapCount = 0;
		for(Nano n : nanos){
			if(nano.overlap(n)){
				return 0;
			}
		}
		
		for (int i = 0; i < polymers.length; i++) {
			if(nano.overlap(polymers[i])){
				overlapCount++;
			}
		}
		
		// Count other nano-polymer intersections
		for(Nano n : nanos){
			for(Polymer p : polymers){
				if(n.overlap(p)){
					overlapCount++;
				}
			}
		}
		
		return Math.exp(-Ep*overlapCount); 
	}


	/**
	 * Returns the polymer shape probability for a given radius axis.
	 * 
	 * @param ei
	 *            The eigenvalue of the radius axis
	 * @param v
	 *            Enumeration representing the radius axis.
	 * @return Polymer shape probability.
	 */
	public double prob(double ei, Vector v) {
		switch (v) {
		case x:
			return (Math.pow(ei, -NX) * Math.pow(AX * DX, NX - 1) / (2 * KX))
					* Math.exp(-ei / AX - DX * DX * AX / ei);
		case y:
			return (Math.pow(ei, -NY) * Math.pow( (AY * DY), (NY - 1)) / (2 * KY))
					* Math.exp( (-ei / AY) - (DY * DY * AY / ei));
		case z:
			return (Math.pow(ei, -NZ) * Math.pow(AZ * DZ, NZ - 1) / (2 * KZ))
					* Math.exp(-ei / AZ - DZ * DZ * AZ / ei);
		}
		return 0;
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
			total += poly.geteX() + poly.geteY() + poly.geteZ();
		}

		average = Math.sqrt(total / polymers.length);
		ratio = average / Nano.getDefault_r();
		return ratio;
	}
}
