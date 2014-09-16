package org.opensourcephysics.sip.CPM;

import java.awt.Color;
import java.text.DecimalFormat;

import org.opensourcephysics.controls.AbstractSimulation;
import org.opensourcephysics.controls.SimulationControl;
import org.opensourcephysics.display3d.simple3d.ElementEllipsoid;
import org.opensourcephysics.display3d.simple3d.ElementSphere;
import org.opensourcephysics.frames.Display3DFrame;
import org.opensourcephysics.frames.PlotFrame;
import org.opensourcephysics.numerics.PBC;

/**
 * NanoPolyMixApp is a simulation framework for a binary mixture of
 * colloids(nanoparticles) and polymers model.
 * 
 * @author Wei Kang Lim, Alan Denton
 * @version 1.0 11-10-2013
 * 
 */
public class CPMApp extends AbstractSimulation {
	public enum WriteModes{WRITE_NONE,WRITE_SHAPES,WRITE_ROTATIONS,WRITE_RADIAL,WRITE_POMF,WRITE_ALL;};
	CPM np = new CPM();
	Display3DFrame display3d = new Display3DFrame("3D Frame");
	PlotFrame plotframe = new PlotFrame("Monte Carlo Steps",
			"Number of Intersections", "Number of Intersections");
	double totalIntersections = 0;
	double snapshotIntervals = 1;
	double [] zAxis = {0,0,1};
	double [] xAxis = {1,0,0};
	RDF rdf;
	double polar;
	double azimuth;
	double volFraction;
	final double RADIAL_START = 0.9; 
	final double RADIAL_END = 2.6;
	double sumVolume;
	double volumeSnapshots;
	double sumDistribution;
	int conformations;
	int maxConformations;
	int dataPoints;
	int maxDataPoints;
	double steps;
	int i =0;
	WriteModes writeMode;
	ElementSphere nanoSphere[];
	ElementEllipsoid polySphere[];
	boolean added = false;
	boolean penetrationEnergyToggle;
	DataFile [] dataFiles;

	public CPMApp(){
		if(display3d != null) display3d.dispose();
		display3d = new Display3DFrame("3D Frame");

	}
	/**
	 * Initializes the model.
	 */
	public void initialize() {
		i++;
		sumDistribution = 0;
		totalIntersections = 0;
		dataPoints = 0;
		conformations = 0;
		np.nP = control.getInt("N Polymers");
		np.nN = 1;
		np.q = control.getDouble("Polymer colloid ratio");
		np.Lx = np.Ly = np.Lz = control.getDouble("Lattice length");
		np.nano_r = 0.5;
		np.init_eX = control.getDouble("x");
		np.init_eY = control.getDouble("y");
		np.init_eZ = control.getDouble("z");
		String configuration = control.getString("Initial configuration");
		np.tolerance = control.getDouble("Tolerance");
		np.rotTolerance = control.getDouble("Rotation tolerance");
		np.shapeTolerance = control.getDouble("Shape tolerance");
		np.trialMovesPerMcs = control.getInt("Trial moves per MCS");
		snapshotIntervals = control.getInt("Snapshot interval");
		maxConformations = control.getInt("Number of conformations");
		maxDataPoints = control.getInt("Number of datapoints");
		penetrationEnergyToggle =control.getBoolean("Penetration energy");
		steps = (RADIAL_END-RADIAL_START) / (maxDataPoints - 1);
		switch(control.getInt("Write Mode")){
		case 0: writeMode = WriteModes.WRITE_NONE; break;
		case 1: writeMode = WriteModes.WRITE_SHAPES; break;
		case 2: writeMode = WriteModes.WRITE_ROTATIONS; break;
		case 3: writeMode = WriteModes.WRITE_RADIAL; break;
		case 4: writeMode = WriteModes.WRITE_POMF; break;
		case 5: writeMode = WriteModes.WRITE_ALL; break;
		}
		
		np.initialize(configuration, penetrationEnergyToggle);
		
		// add simple3d.Element particles to the arrays 
		if (i == 2) { // particles only allowed to be added once, this is to prevent clone particles when initialize is called twice by the simulation
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
			
			// Initialize visualization elements for nano particles
			for (int i = 0; i < np.nN; i++) {
				nanoSphere[i].setSizeXYZ(2 * np.nanos[i].getrX(),
						2 * np.nanos[i].getrY(), 2 * np.nanos[i].getrZ());
				nanoSphere[i].getStyle().setFillColor(new Color(92, 146, 237, 100));  // light blue with half transparency
				nanoSphere[i].getStyle().setLineColor(new Color(92, 146, 237, 100));
				nanoSphere[i].setXYZ(np.nanos[i].getX(), np.nanos[i].getY(),
						np.nanos[i].getZ());
			}

			// Initialize visualization elements for polymer particles
			for (int i = 0; i < np.nP; i++) {
				polySphere[i].setSizeXYZ(2 * np.polymers[i].getrX(),
						2 * np.polymers[i].getrY(), 2 * np.polymers[i].getrZ());
				polySphere[i].getStyle().setFillColor(Color.RED);
				polySphere[i].getStyle().setLineColor(Color.RED);
				polySphere[i].setXYZ(np.polymers[i].getX(), np.polymers[i].getY(),
						np.polymers[i].getZ());
	            if(np.polymers[i].isRotated()){ 
	                    polySphere[i].setTransformation(np.polymers[i].getTransformation());
	            }
			}

			plotframe.append(0, np.mcs, np.totalIntersectCount);
			display3d.setPreferredMinMax(0, np.Lx, 0, np.Ly, 0, np.Lz);
			display3d.setSquareAspect(true);
		}		
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
			String configurations = "# Number of Polymers: " + np.nP +
					"\n# Number of Nanoparticles: "+np.nN +
					"\n# Move Tolerance: "+threeDecimal.format(np.tolerance)+
					"\n# Shape Change Tolerance: "+threeDecimal.format(np.shapeTolerance)+
					"\n# Nanoparticle Volume Fraction: "+threeDecimal.format(np.volFraction) + 
					"\n# Polymer Colloid Ratio: "+threeDecimal.format(np.q)+
					"\n# Lattice Constant: " +threeDecimal.format(np.lc)+
					"\n# Rotation Tolerance: "+threeDecimal.format(np.rotTolerance)+
					"\n# Trial Moves Per Mcs: "+np.trialMovesPerMcs+
					"\n# Snapshot Interval: "+largeDecimal.format(this.snapshotIntervals)+
					"\n# Number of Coformations Sampled: " + maxConformations +
					"\n# Number of dataPoints: " + maxDataPoints + 
					"\n# Penetration Energy On: " + this.penetrationEnergyToggle
					;
			switch(writeMode){
			case WRITE_SHAPES:
				dataFiles = new DataFile[3];
				dataFiles[0] = new DataFile("eX", configurations);
				dataFiles[1] = new DataFile("eY", configurations);
				dataFiles[2] = new DataFile("eZ", configurations);
				break;
			case WRITE_ROTATIONS:
				dataFiles = new DataFile[2];
				dataFiles[0] = new DataFile("polar", configurations);
				dataFiles[1] = new DataFile("azimuth", configurations);
				break;
			case WRITE_RADIAL:
				dataFiles = new DataFile[1]; 
				dataFiles[0] = new DataFile("radial", configurations);
				rdf = new RDF(np.nanos, np.Lx);
				break;
			case WRITE_POMF:
				dataFiles = new DataFile[1];
				dataFiles[0] = new DataFile("POMF", configurations);
				break;
			case WRITE_ALL:
				dataFiles = new DataFile[4];
				dataFiles[0] = new DataFile("eX", configurations);
				dataFiles[1] = new DataFile("eY", configurations);
				dataFiles[2] = new DataFile("eZ", configurations);
				dataFiles[3] = new DataFile("radial", configurations);
				rdf = new RDF(np.nanos, np.Lx);
				break;
			default:
				break;
			}
			
			// Write out the initialization data
			for(DataFile df : dataFiles){
				df.write();
			}
		}

		// logical step in the CPM class
		np.step();
		
		// update intersect count and mcs step
		display3d.setMessage("Number of mcs steps: " + np.mcs);

		// Visualization updates
		if(control.getBoolean("Visualization on")){
			// nanoparticle visualization updates
			for (int i = 0; i < np.nN; i++) {
//				if(np.nanos[i].intersectPairs.isEmpty()){
//					nanoSphere[i].getStyle().setFillColor(new Color(92, 146, 237, 100));  // light blue with half transparency
//					nanoSphere[i].getStyle().setLineColor(new Color(92, 146, 237, 100));
//				}else{
//					nanoSphere[i].getStyle().setFillColor(new Color(0, 255, 0, 100));  // light blue with half transparency
//					nanoSphere[i].getStyle().setLineColor(new Color(0, 255, 0, 100));
//				}
				nanoSphere[i].setXYZ(np.nanos[i].getX(), np.nanos[i].getY(),
						np.nanos[i].getZ());
			}
			// polymer visualization updates
			for (int i = 0; i < np.nP; i++) {
				polySphere[i].setXYZ(np.polymers[i].getX(), np.polymers[i].getY(),
						np.polymers[i].getZ());
				polySphere[i].setSizeXYZ(2 * np.polymers[i].getrX(),
						2 * np.polymers[i].getrY(), 2 * np.polymers[i].getrZ());
				if(np.polymers[i].isRotated()){ // this is done to save computation
					polySphere[i].setTransformation(np.polymers[i].getTransformation());
				}
			}
		}
		
		// perform insertion algorithm and record the data
		if (writeMode != WriteModes.WRITE_NONE && np.mcs >= 50000 && np.mcs % snapshotIntervals == 0) {
			if(conformations > maxConformations){
				double avgDistribution = sumDistribution / maxConformations;
				dataFiles[0].record(RADIAL_START+dataPoints*steps + " " + avgDistribution);
				dataPoints++;
				if(dataPoints > maxDataPoints){
					control.setAdjustableValue("Save", true);
					this.stopAnimation();
					return;
				}
				conformations = 0;
				sumDistribution = 0;
			}
			
			// record volume fraction
			for(Polymer poly : np.polymers){
				sumVolume += (4./3.) * Math.PI * poly.getrX() * poly.getrY() * poly.getrZ();
			}		
			volumeSnapshots++;
			
			double e_delU = np.nanoTrialPlacement(RADIAL_START+dataPoints*steps);
			sumDistribution += e_delU;
			plotframe.append(0, np.mcs, e_delU);
			conformations++;
		}		
		
		if(writeMode != WriteModes.WRITE_NONE && np.mcs >= 50000 && np.mcs % (snapshotIntervals*100) == 0){
			dataFiles[0].write();
		}
	}

	/**
	 * Resets the colloid polymer mixture model to its default state.
	 */
	public void reset() {
		enableStepsPerDisplay(true);
		control.setValue("N Polymers", 13);
		control.setValue("Polymer colloid ratio", 2);
		control.setValue("Lattice length", 3.74);
		control.setValue("x", 0.01);
		control.setValue("y", 0.01);
		control.setValue("z", 0.01);
		control.setValue("Tolerance", 0.1);
		control.setValue("Rotation tolerance", 0.1);
		control.setValue("Shape tolerance", 0.001);
		control.setValue("Initial configuration", "square");
		control.setValue("Trial moves per MCS", 1);
		control.setAdjustableValue("Visualization on", true);
		control.setValue("Snapshot interval", 1000);
		control.setValue("Number of datapoints", 50);
		control.setValue("Number of conformations", 100);
		control.setValue("Penetration energy", true);
		control.setValue("Write Mode", 4);
		control.setAdjustableValue("Save", false);
		initialize();
	}

	public void stop() {
		double averageIntersections = totalIntersections / np.mcs;
		if(control.getInt("N Polymers")> 0){
			
		double volFract = sumVolume / (volumeSnapshots* np.Lx * np.Ly * np.Lz);
		control.println("Volume fraction of nanoparticles: " + Math.PI/(6 * np.Lx * np.Ly * np.Lz));
		control.println("Volume fraction of polymers: " + volFract);
		control.println("Average no. of Intersections: " + averageIntersections);
		control.println("Expected average no. of intersections with Ep = 0: " + volFract);
		
		}

		// close streams
		if(control.getBoolean("Save")){
			if(writeMode != WriteModes.WRITE_NONE){
				control.setAdjustableValue("Save", false);
				if(writeMode == WriteModes.WRITE_RADIAL){ // Radial Distribution Function 
					dataFiles[0].record(rdf.distributionData());
					System.out.println(rdf.nrData());
					dataFiles[0].write();
				}
				
				if(writeMode == WriteModes.WRITE_ALL){ // Radial Distribution Function 
					dataFiles[3].record(rdf.distributionData());
					System.out.println(rdf.nrData());
					dataFiles[3].write();
				}
				
				for(DataFile df : dataFiles){
					df.record("#Number of data points: " + dataPoints);
					df.close();
				}
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
