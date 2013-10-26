package org.opensourcephysics.sip.CPM;

import java.awt.Color;
import java.awt.Graphics;

import org.opensourcephysics.display.Drawable;
import org.opensourcephysics.display.DrawingPanel;

/**
 * NanoPolyMix is an abstraction for a binary mixture of colloids(nanoparticles)
 * and polymers model.
 * 
 * @author Wei Kang Lim, Alan Denton
 * @version 0.5 beta 7-1-2013
 * 
 */
public class CPM {
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
	public int nN, nP, nx; // number of nanoparticles, number of polymers,
							// number of columns and rows
	// public double rho; // Ratio of polymer diameter to nanoparticle (colloid)
	// diameter, sigP/sigN
	public double init_eX;
	public double init_eY;
	public double init_eZ;
	public double sigmaX;
	public double sigmaY;
	public double nano_r;
	public double q;
	public int moveToShapeRatio;
	public double lc; // lattice constant, defined as lc = d / sigN, where d is
						// defined below.
	public double Lx; // dimension of the boundary (along x)
	public double Ly; // dimension of the boundary (along y)
	public double Lz; // dimension of the boundary (along z)
	public double steps = 0; // number of monte carlo steps
	public double tolerance;
	public double shapeTolerance;
	public double d; // the center-to-center distance between particles of same
						// species (hexagonal lattice),
						// or columns and rows (square lattice)
	public double intersectCount = 0;
	public double Ep; // Penetration Energy
	public double mcs = 0;

	// end declaration

	public enum Vector {
		x, y, z
	};

	/**
	 * initializes the model.
	 * 
	 * @param configuration
	 *            the initial lattice structure of the model
	 */
	public void initialize(String configuration) {
		polymers = new Polymer[nP];
		nanos = new Nano[nN];
		d = lc; // distance between two nanoparticles
		
		Ep = 3/q;

		// Set a constant distribution
		Particle.setBoundaries(Lx, Ly, Lz);
		Nano.setTolerance(tolerance);
		Polymer.setTolerance(tolerance);
		Polymer.setQ(q);
		Nano.setDefault_r(nano_r);
		Polymer.setDefault_eX(init_eX);
		Polymer.setDefault_eY(init_eY);
		Polymer.setDefault_eZ(init_eZ);
		intersectCount = 0;

		if (configuration.toUpperCase().equals("SQUARE")) {
			setSqrPositions();
		} else {
			// setHexPositions();
		}
	}

	/**
	 * Places particles on a square lattice.
	 */
	public void setSqrPositions() {
		int ix, iy, iz;
		if(nN > 0){
			double dnx = Math.cbrt(nN);
			nx = (int) dnx;
			if (dnx - nx > 0.00001) {
				nx++; // N is not a perfect cube
			}
			Lx = Ly = Lz = d * nx;
	
			int i = 0;
			for (iy = 0; iy < nx; iy++) { // loops through particles in a column
				for (ix = 0; ix < nx; ix++) { // loops through particles in a
												// row
					for (iz = 0; iz < nx; iz++) {
						if (i < nanos.length) { // checks for remaining particles
							nanos[i] = new Nano(ix * d, iy * d, iz * d);
							i++;
						}
					}
				}
			}
	
			i = 0;
			for (iy = 0; iy < nx; iy++) {
				for (ix = 0; ix < nx; ix++) {
					for (iz = 0; iz < nx; iz++) {
						if (i < polymers.length) { // checks for remaining
													// particles
							polymers[i] = new Polymer(ix * d + d / 2, iy * d + d
									/ 2, iz * d + d / 2);
							System.out.println(polymers[i].getrX() + " " + polymers[i].getrY() + " " + polymers[i].getrZ());
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
		} else{
			double dnx = Math.cbrt(nP);
			nx = (int) dnx;
			if (dnx - nx > 0.00001) {
				nx++; // N is not a perfect cube
			}
			Lx = Ly = Lz = d * nx;
			int i = 0;
			
			for (iy = 0; iy < nx; iy++) {
				for (ix = 0; ix < nx; ix++) {
					for (iz = 0; iz < nx; iz++) {
						if (i < polymers.length) { // checks for remaining
													// particles
							polymers[i] = new Polymer(ix * d + d / 2, iy * d + d
									/ 2, iz * d + d / 2);
							System.out.println(polymers[i].getrX() + " " + polymers[i].getrY() + " " + polymers[i].getrZ());
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

	public void trialMoves() {
		// Nanoparticles Trial Moves
		for (int i = 0; i < nanos.length; ++i) {
			nanoTrialMove(nanos[i]);
		}

		// Polymer Trial Moves
		for (int i = 0; i < polymers.length; ++i) {
			polyTrialMove(polymers[i]);
			 shapeChange(polymers[i]);
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
		poly.move();

		for (int i = 0; i < nanos.length; i++) {
			if (poly.overlap(nanos[i])) { // Nano-polymer now overlaps
				if (!poly.intersectPairs.contains(nanos[i])
						&& !nanos[i].intersectPairs.contains(poly)) { // Check
																		// if
																		// the
																		// particles
																		// were
																		// overlapping
					if (Math.random() < Math.exp(-Ep)) { // Change in number of
															// CP penetrations =
															// +1, accepted
						intersectCount++;
						poly.intersectPairs.add(nanos[i]); // Add
						nanos[i].intersectPairs.add(poly);
					} else {
						poly.setX(oldX);
						poly.setY(oldY);
						poly.setZ(oldZ);
						return;
					}
				}
			}
		}

		// Since trial move was accepted, let's look at possible intersections
		// being removed.
		for (int j = 0; j < nanos.length; j++) {
			if (!poly.overlap(nanos[j])) { // particles that are no longer
											// overlapping
				if (poly.intersectPairs.remove(nanos[j])
						&& nanos[j].intersectPairs.remove(poly)) // update the
																	// intersecting
																	// pairs and
																	// count
					intersectCount--;
			}
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
		nano.move();

		for (int j = 0; j < nanos.length; ++j) {
			if (!nano.equals(nanos[j])) { // Iterate through all the other
											// particles and check for an
											// overlap
				if (nano.overlap(nanos[j])) { // Nano-nano overlap
					nano.setX(oldX);
					nano.setY(oldY);
					nano.setZ(oldZ);
					return;
				}
			}
		}

		for (int i = 0; i < polymers.length; i++) {
			if (nano.overlap(polymers[i])) { // Nano-polymer now overlaps
				if (!nano.intersectPairs.contains(polymers[i])
						&& !polymers[i].intersectPairs.contains(nano)) { // Check
																			// if
																			// the
																			// particles
																			// were
																			// overlapping
					if (Math.random() < Math.exp(-Ep)) { // Change in number of
															// CP penetrations =
															// +1, accepted
						intersectCount++;
						nano.intersectPairs.add(polymers[i]); // Add
						polymers[i].intersectPairs.add(nano);
					} else {
						nano.setX(oldX);
						nano.setY(oldY);
						nano.setZ(oldZ);
						return;
					}
				}
			}
		}

		// Since trial move was accepted, let's look at possible intersections
		// being removed.
		for (int j = 0; j < polymers.length; j++) {
			if (!nano.overlap(polymers[j])) { // particles that are no longer
												// overlapping
				if (nano.intersectPairs.remove(polymers[j])
						&& polymers[j].intersectPairs.remove(nano)) // update
																	// the
																	// intersecting
																	// pairs and
																	// count
					intersectCount--;
			}
		}
	}

	/**
	 * Shape changes for a polymer.
	 * 
	 * @param poly A polymer object to be attempted for a trial shape change.
	 */
	public void shapeChange(Polymer poly) {
//		double q = polymerColloidSizeRatio();
		double oldRX = poly.getrX();
		double oldRY = poly.getrY();
		double oldRZ = poly.getrZ();
		
		double oldEX = poly.geteX();
		double oldEY = poly.geteY();
		double oldEZ = poly.geteZ();

		// trial shape changes
		double newEX = oldEX + shapeTolerance * 2. * (Math.random() - 0.5);
		double newEY = oldEY + shapeTolerance * 2. * (Math.random() - 0.5);
		double newEZ = oldEZ + shapeTolerance * 2. * (Math.random() - 0.5);

		double P_x = prob(newEX, Vector.x) / prob(oldEX, Vector.x);
		double P_y = prob(newEY, Vector.y) / prob(oldEY, Vector.y);
		double P_z = prob(newEZ, Vector.z) / prob(oldEZ, Vector.z);

		boolean rx_changed = attemptEigenRadiusChange(P_x, poly,
				newEY, q, Vector.x);
		boolean ry_changed = attemptEigenRadiusChange(P_y, poly,
				newEY, q, Vector.y);
		boolean rz_changed = attemptEigenRadiusChange(P_z, poly,
				newEZ, q, Vector.z);
		
		if(rx_changed){
			if((poly.getX() + poly.getrX()) > Lx || (poly.getX() - poly.getrX()) < 0 ){
				rx_changed = false;
				System.out.println("reject");
				poly.seteX(oldEX, q);
				System.out.println(poly.getrX() - oldRX);
				return;
			}
			System.out.println("Accepted DR_X= " + (poly.getrX() - oldRX) + " with P " + P_x);
		}
		
		if(ry_changed){
			if((poly.getY() + poly.getrY()) > Lx || (poly.getY() - poly.getrY()) < 0 ){
				ry_changed = false;
				poly.seteY(oldEY, q);
				return;
			}
		}
		if(rz_changed){
			if((poly.getZ() + poly.getrZ()) > Lx || (poly.getZ() - poly.getrZ()) < 0 ){
				rz_changed = false;
				poly.seteZ(oldEZ, q);
				return;
			}
		}

		// Check for intersections with nanoparticles
		for (int i = 0; i < nanos.length; i++) {
			if (poly.overlap(nanos[i])) { // Nano-polymer now overlaps
				if (!poly.intersectPairs.contains(nanos[i])
						&& !nanos[i].intersectPairs.contains(poly)) { // Check
																		// for
																		// previous
																		// overlap
					if (Math.random() < Math.exp(-Ep)) {
						intersectCount++;
						poly.intersectPairs.add(nanos[i]); // Add
						nanos[i].intersectPairs.add(poly);
					} else {
						if (rx_changed)
							poly.seteX(oldEX, q);
						if (ry_changed)
							poly.seteY(oldEY, q);
						if (rz_changed)
							poly.seteZ(oldEZ, q);
						return;
					}
				}
			}
		}

		// Since trial move was accepted, let's look at possible intersections
		// being removed.
		for (int j = 0; j < nanos.length; j++) {
			if (!poly.overlap(nanos[j])) { // particles that are no longer
											// overlapping
				if (poly.intersectPairs.remove(nanos[j])
						&& nanos[j].intersectPairs.remove(poly)) // update the
																	// intersecting
																	// pairs and
																	// count
					intersectCount--;
			}
		}
	}

	/**
	 * Returns the polymer shape probability for a given radius axis. 
	 * @param ei The eigenvalue of the radius axis
	 * @param v Enumeration representing the radius axis.
	 * @return Polymer shape probability.
	 */
	public double prob(double ei, Vector v) {
		switch (v) {
		case x:
			return (Math.pow(ei, -NX) * Math.pow(AX * DX, NX - 1) / (2 * KX))
					* Math.exp(-ei / AX - DX * DX * AX / ei);
		case y:
			return (Math.pow(ei, -NY) * Math.pow(AY * DY, NY - 1) / (2 * KY))
					* Math.exp(-ei / AY - DY * DY * AY / ei);
		case z:
			return (Math.pow(ei, -NZ) * Math.pow(AZ * DZ, NZ - 1) / (2 * KZ))
					* Math.exp(-ei / AZ - DZ * DZ * AZ / ei);
		default:
			return 0;
		}
	}

	/**
	 * Attempts a change on one of the radii eigenvalue of a polymer based on the
	 * acceptance probability. The radius of the polymer changes accordingly.
	 * 
	 * @param p
	 *            Acceptance probability
	 * @param poly
	 *            Polymer in question
	 * @param newE
	 *            The new value of the radius eigenvalue
	 * @param v
	 *            enumeration representing the radius axis
	 * @param q  Colloid-polymer ratio
	 * @return true if the change succeeded, false otherwise
	 */
	public boolean attemptEigenRadiusChange(double p, Polymer poly, double newE,
			double q,Vector v) {
		if (p > 1 || Math.random() < p) {
			switch (v) {
			case x:
				poly.seteX(newE, q);	// Note that the polymer calculates its new radii whenever you change its eigenvalue
			case y:
				poly.seteY(newE, q);
			case z:
				poly.seteZ(newE, q);
			default:
				return true;
			}
		} else {
			return false;
		}
	}

	/**
	 * Calculates the polymer colloid size ratio.
	 * @param v enumeration representing the radius axis
	 * @return Polymer colloid size ratio.
	 */
	public double polymerColloidSizeRatio() {
		double total = 0;
		double average = 0;
		double ratio = 0;

		for (Polymer poly : polymers) {
				total += poly.geteX()*poly.geteX() + poly.geteY()*poly.geteY() + poly.geteZ()*poly.geteZ();
		}

		average = Math.sqrt( total / polymers.length );
		ratio = average / Nano.getDefault_r();
		return ratio;
	}
}
