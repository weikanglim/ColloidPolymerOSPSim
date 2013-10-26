package org.opensourcephysics.sip.CPM;
import java.awt.Color;

import org.opensourcephysics.controls.*;
import org.opensourcephysics.display3d.simple3d.ElementEllipsoid;
import org.opensourcephysics.display3d.simple3d.ElementSphere;
import org.opensourcephysics.frames.*;

/**
 * NanoPolyMixApp is a simulation framework for a binary mixture of colloids(nanoparticles) and polymers model.
 * @author Wei Kang Lim, Alan Denton
 * @version 0.5 beta 7-1-2013
 *
 */
public class CPMApp extends AbstractSimulation{
	 CPM np = new CPM();
	 Display3DFrame display3d = new Display3DFrame("3D Frame");
	 PlotFrame plotframe = new PlotFrame("Monte Carlo Steps", "Number of Intersections", "Number of Intersections");
	 int mcs;
	 double totalIntersections = 0;
	 ElementSphere nanoSphere[];
	 ElementEllipsoid polySphere[];
	 boolean added = false;
	 
	  /**
	   * Initializes the model.
	   */
	  public void initialize() {
		totalIntersections = 0;
		mcs = 0;
		np.nP = control.getInt("N Polymers");
	    np.nN = control.getInt("N Nano");
	    np.tolerance = control.getDouble("tolerance");
	    np.shapeTolerance = control.getDouble("Shape Tolerance");
	    np.q = control.getDouble("Polymer Colloid Ratio");
	    np.init_eX = control.getDouble("x");
	    np.init_eY = control.getDouble("y");
	    np.init_eZ = control.getDouble("z");
	    np.nano_r = control.getDouble("Nanoparticle radius");
	    np.lc = control.getDouble("Lattice constant");
	    np.Ep = control.getDouble("Penetration Energy");
	    String configuration = control.getString("initial configuration");
	    np.moveToShapeRatio = control.getInt("Trial Moves to Shape Changes Ratio");
	    np.initialize(configuration);
	    display3d.setPreferredMinMax(0, np.Lx, 0, np.Ly, 0, np.Lz);
	    display3d.setSquareAspect(true);

		if (!added) {
			nanoSphere = new ElementSphere[np.nN];
			polySphere = new ElementEllipsoid[np.nP];
			
			for (int i = 0; i < np.nN; i++) {
				nanoSphere[i] = new ElementSphere();
				display3d.addElement(nanoSphere[i]);
			}
			
			for (int i = 0; i < np.nP; i++) {
				polySphere[i] = new ElementEllipsoid();
				display3d.addElement(polySphere[i]);
			}
			
			added = true;
		}
	    
	    for(int i = 0; i < np.nN; i++){
	        nanoSphere[i].setSizeXYZ(2*np.nanos[i].getrX(), 2*np.nanos[i].getrY(), 2*np.nanos[i].getrZ());
	        nanoSphere[i].getStyle().setFillColor(Color.BLACK);
	        nanoSphere[i].setXYZ(np.nanos[i].getX(), np.nanos[i].getY(), np.nanos[i].getZ());
	    }
	    
	    for(int i = 0; i < np.nP; i++){
	        polySphere[i].setSizeXYZ(2*np.polymers[i].getrX(), 2*np.polymers[i].getrY(), 2*np.polymers[i].getrZ());
	        polySphere[i].getStyle().setFillColor(Color.RED);
	        polySphere[i].setXYZ(np.polymers[i].getX(), np.polymers[i].getY(), np.polymers[i].getZ());
	    }
	    
	    plotframe.append(0, np.mcs, np.intersectCount);
	  }

	  /**
	   * Does a simulation step.
	   */
	  public void doStep() {
		np.step();
		
		for (int i = 0; i < np.nN; i++) {
	        nanoSphere[i].setXYZ(np.nanos[i].getX(), np.nanos[i].getY(), np.nanos[i].getZ());
		}
		for (int i = 0; i < np.nP; i++) {
	        polySphere[i].setXYZ(np.polymers[i].getX(), np.polymers[i].getY(), np.polymers[i].getZ());
	        polySphere[i].setSizeXYZ(2*np.polymers[i].getrX(), 2*np.polymers[i].getrY(), 2*np.polymers[i].getrZ());
		}
		
		plotframe.append(0, np.mcs, np.intersectCount);
		totalIntersections += np.intersectCount;
		display3d.setMessage("Number of mcs steps: " + np.mcs);
	}

	  /**
	   * Resets the colloid polymer mixture model to its default state.
	   */
	  public void reset() {
	    enableStepsPerDisplay(true);
	    control.setValue("N Polymers", 27);
	    control.setValue("N Nano", 27);
	    control.setValue("tolerance", 0.1);
	    control.setValue("Shape Tolerance", 0.1);
	    control.setValue("Nanoparticle radius", 0.1);
	    control.setValue("x", 0.001);
	    control.setValue("y", 0.001);
	    control.setValue("z", 0.001);
	    control.setValue("Polymer Colloid Ratio", 5);
	    control.setValue("Lattice constant", 3);
	    control.setValue("initial configuration", "square");
	    control.setValue("Penetration Energy", 0);
	    control.setValue("Trial Moves to Shape Changes Ratio", 5);	    
	    initialize();
	  }
	  
	  public void stop(){
		  double averageIntersections = totalIntersections / np.mcs;
		  double volSpheres = np.nP*(4/3d*Math.PI*np.init_eX*np.init_eX*np.init_eX)*np.nN;
		  double volFract = volSpheres / (np.Lx*np.Ly*np.Lz);
		  control.println("Average no. of Intersections: " + averageIntersections);
		  control.println("Vol fraction: " + volFract);
	  }

	  /**
	   * Starts the Java application.
	   * @param args  command line parameters
	   */
	  public static void main(String[] args) { // set up animation control structure using this class
	    @SuppressWarnings("unused")
		SimulationControl control = SimulationControl.createApp(new CPMApp());
	  }
}
