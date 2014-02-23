package org.opensourcephysics.sip.CPM;

import java.awt.Color;
import java.text.DecimalFormat;

import org.opensourcephysics.controls.AbstractSimulation;
import org.opensourcephysics.controls.SimulationControl;
import org.opensourcephysics.display3d.simple3d.ElementEllipsoid;
import org.opensourcephysics.display3d.simple3d.ElementSphere;
import org.opensourcephysics.frames.Display3DFrame;
import org.opensourcephysics.frames.PlotFrame;

/**
 * NanoPolyMixApp is a simulation framework for a binary mixture of
 * colloids(nanoparticles) and polymers model.
 * 
 * @author Wei Kang Lim, Alan Denton
 * @version 1.0 11-10-2013
 * 
 */
public class CPMApp extends AbstractSimulation {
	public enum WriteModes{WRITE_NONE,WRITE_SHAPES,WRITE_ROTATIONS,WRITE_ALL};
	CPM np = new CPM();
	Display3DFrame display3d = new Display3DFrame("3D Frame");
	PlotFrame plotframe = new PlotFrame("Monte Carlo Steps",
			"Number of Intersections", "Number of Intersections");
	double totalIntersections = 0;
	double snapshotIntervals = 1;
	double polar;
	double azimuth;
	WriteModes writeMode;
	ElementSphere nanoSphere[];
	ElementEllipsoid polySphere[];
	boolean added = false;
	boolean penetrationEnergyToggle;
	DataFile dataFile;

	/**
	 * Initializes the model.
	 */
	public void initialize() {
		added = false;
		totalIntersections = 0;
		np.nP = control.getInt("N Polymers");
		np.nN = control.getInt("N Nano");
		np.tolerance = control.getDouble("tolerance");
		np.shapeTolerance = control.getDouble("Shape Tolerance");
		np.q = control.getDouble("Polymer Colloid Ratio");
		np.init_U = control.getDouble("u");
		np.lc = control.getDouble("Lattice constant");
		String configuration = control.getString("initial configuration");
		np.trialDisplacementPerMcs = control.getInt("Trial Disp Per Mcs");
		np.trialShapeChangePerMcs = control.getInt("Trial Shape Change Per Mcs");
		snapshotIntervals = control.getInt("Snapshot Interval");
		penetrationEnergyToggle =control.getBoolean("Penetration Energy");
		switch(control.getInt("Write Mode")){
		case 0: writeMode = WriteModes.WRITE_NONE; break;
		case 1: writeMode = WriteModes.WRITE_SHAPES; break;
		case 2: writeMode = WriteModes.WRITE_ROTATIONS; break;
		case 3: writeMode = WriteModes.WRITE_ALL; break;
		}
		
		if(!penetrationEnergyToggle){ 
			System.out.println("Penetration energy turned off.");// warning for user
			np.Ep = 0;
		}
		
		np.initialize(configuration);
		if(display3d != null) display3d.dispose(); // closes an old simulation frame is present
		display3d = new Display3DFrame("3D Frame");
		display3d.setPreferredMinMax(0, np.Lx, 0, np.Ly, 0, np.Lz);
		display3d.setSquareAspect(true);
		
		// add simple3d.Element particles to the arrays 
		if (!added) { // particles only allowed to be added once, this is to prevent clone particles when initialize is called twice by the simulation
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
		
		// Initialize visualization elements for nano particles
		for (int i = 0; i < np.nN; i++) {
			nanoSphere[i].setSizeXYZ(2 * np.nanos[i].getrX(),
					2 * np.nanos[i].getrY(), 2 * np.nanos[i].getrZ());
			nanoSphere[i].getStyle().setFillColor(Color.BLACK);
			nanoSphere[i].setXYZ(np.nanos[i].getX(), np.nanos[i].getY(),
					np.nanos[i].getZ());
		}

		// Initialize visualization elements for polymer particles
		for (int i = 0; i < np.nP; i++) {
			polySphere[i].setSizeXYZ(2 * np.polymers[i].getrX(),
					2 * np.polymers[i].getrY(), 2 * np.polymers[i].getrZ());
			polySphere[i].getStyle().setFillColor(Color.RED);
			polySphere[i].setXYZ(np.polymers[i].getX(), np.polymers[i].getY(),
					np.polymers[i].getZ());
		}

		plotframe.append(0, np.mcs, np.totalIntersectCount);
		}
	
	/**
	 * Does a simulation step.
	 */
	public void doStep() {
		// Initialize files for writing output data during the first step
		// This code is being ran here 
		if(np.mcs == 0 && writeMode != WriteModes.WRITE_NONE){
			DecimalFormat largeDecimal = new DecimalFormat("0.##E0");
			DecimalFormat threeDecimal = new DecimalFormat("#0.###");
			String configurations = "# Number of Polymers:" + np.nP +
					"\n# Number of Nanoparticles: "+np.nN +
					"\n# Move Tolerance: "+threeDecimal.format(np.tolerance)+
					"\n# Shape Change Tolerance: "+threeDecimal.format(np.shapeTolerance)+
					"\n# Nanoparticle Radius :"+threeDecimal.format(Nano.getDefault_r()) + 
					"\n# Polymer Colloid Ratio: "+threeDecimal.format(np.q)+
					"\n# Lattice Constant: " +threeDecimal.format(np.lc)+
					"\n# Trial Moves to Attempt Shape Change Ratio: "+np.moveToShapeRatio+ // !TODO
					"\n# Snapshot Interval: "+largeDecimal.format(this.snapshotIntervals)+
					"\n# Penetration Energy On: " + this.penetrationEnergyToggle
					;
			switch(writeMode){
			case WRITE_SHAPES:
				dataFile = new DataFile("u", configurations);
				break;
			case WRITE_ALL:
				dataFile = new DataFile("u", configurations);
				break;
			default:
				break;
			}
		}

		// logical step in the CPM class
		np.step();
		
		// update intersect count and mcs step
		plotframe.append(0, np.mcs, np.totalIntersectCount);
		totalIntersections += np.totalIntersectCount;
		display3d.setMessage("Number of mcs steps: " + np.mcs);

		// Visualization updates
		if(control.getBoolean("Visualization on")){
			// nanoparticle visualization updates
			for (int i = 0; i < np.nN; i++) {
				if(np.nanos[i].intersectPairs.isEmpty()){
					nanoSphere[i].getStyle().setFillColor(Color.BLACK);
				}else{
					nanoSphere[i].getStyle().setFillColor(Color.GREEN);
				}
				nanoSphere[i].setXYZ(np.nanos[i].getX(), np.nanos[i].getY(),
						np.nanos[i].getZ());
			}
			// polymer visualization updates
			for (int i = 0; i < np.nP; i++) {
				polySphere[i].setXYZ(np.polymers[i].getX(), np.polymers[i].getY(),
						np.polymers[i].getZ());
				polySphere[i].setSizeXYZ(2 * np.polymers[i].getrX(),
						2 * np.polymers[i].getrY(), 2 * np.polymers[i].getrZ());
			}
		}
		
		// writing of data
		if (writeMode != WriteModes.WRITE_NONE && np.mcs % snapshotIntervals == 0) {
			switch(writeMode){
			case WRITE_SHAPES:
				if(np.mcs > 50000){ // hardcoded 
					for(Polymer poly: np.polymers){
					dataFile.record(String.valueOf(poly.getU()));
					}
				}
				break;
			case WRITE_ALL:
					for(Polymer poly: np.polymers){
						if(np.mcs > 50000){ // hardcoded 
							dataFile.record(String.valueOf(poly.getU()));
						}
					}
			break;
			default:break;
			}
		}
		
		// write data onto harddisk for every 100 data values
		if(writeMode != WriteModes.WRITE_NONE && np.mcs % (100*snapshotIntervals) == 0){
				dataFile.write();
		}
	}

	/**
	 * Resets the colloid polymer mixture model to its default state.
	 */
	public void reset() {
		enableStepsPerDisplay(true);
		control.setValue("N Polymers", 27);
		control.setValue("N Nano", 27);
		control.setValue("tolerance", 0.1);
		control.setValue("Shape Tolerance", 0.001);
		control.setValue("u", 0.005);
		control.setValue("Polymer Colloid Ratio", 5);
		control.setValue("Lattice constant", 2);
		control.setValue("initial configuration", "square");
		control.setValue("Trial Disp Per Mcs", 1);
		control.setValue("Trial Shape Change Per Mcs", 1);
		control.setAdjustableValue("Visualization on", true);
		control.setValue("Snapshot Interval", 1000);
		control.setValue("Penetration Energy", true);
		control.setValue("Write Mode", 1);
		control.setAdjustableValue("Save", false);
		initialize();
	}

	public void stop() {
		double averageIntersections = totalIntersections / np.mcs;
		
		double volSpheres = np.nP
				* (4 / 3d * Math.PI * np.polymers[0].getrX()
						* np.polymers[0].getrY() * np.polymers[0].getrZ())
				* np.nN;
		double volFract = volSpheres / (np.Lx * np.Ly * np.Lz);
		control.println("Average no. of Intersections: " + averageIntersections);
		control.println("Expected average no. of intersections with Ep = 0: " + volFract);
		
		// close streams
		if(control.getBoolean("Save")){
			if(writeMode != WriteModes.WRITE_NONE){
					dataFile.close();
			}
		}
	}
	
	/**
	 * Starts the Java application.
	 * 
	 * @param args
	 *            command line parameters
	 */
	public static void main(String[] args) { // set up animation control
												// structure using this class
		@SuppressWarnings("unused")
		SimulationControl control = SimulationControl.createApp(new CPMApp());
	}
}
