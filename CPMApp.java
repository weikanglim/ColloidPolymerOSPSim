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
	double [] zAxis = {0,0,1};
	double [] xAxis = {1,0,0};
	double polar;
	double azimuth;
	WriteModes writeMode;
	ElementSphere nanoSphere[];
	ElementEllipsoid polySphere[];
	boolean added = false;
	boolean penetrationEnergyToggle;
	boolean visualizeOn = true;
	DataFile [] dataFiles;

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
		np.init_eX = control.getDouble("x");
		np.init_eY = control.getDouble("y");
		np.init_eZ = control.getDouble("z");
		np.nano_r = control.getDouble("Nanoparticle radius");
		np.lc = control.getDouble("Lattice constant");
		String configuration = control.getString("initial configuration");
		np.rotMagnitude = control.getDouble("Rotation magnitude");
		visualizeOn = control.getBoolean("Visualization");
		np.moveToShapeRatio = control
				.getInt("Trial Moves to Shape Changes Ratio");
		snapshotIntervals = control.getInt("Snapshot Interval");
		penetrationEnergyToggle =control.getBoolean("Penetration Energy");
		switch(control.getInt("Write Mode")){
		case 0: writeMode = WriteModes.WRITE_NONE; break;
		case 1: writeMode = WriteModes.WRITE_SHAPES; break;
		case 2: writeMode = WriteModes.WRITE_ROTATIONS; break;
		case 3: writeMode = WriteModes.WRITE_ALL; break;
		}
		if(!penetrationEnergyToggle){
			System.out.println("Executed!");
			np.Ep = 0;
		}
		np.initialize(configuration);
		if(display3d != null) display3d.dispose(); // closes an old simulation frame is present
		display3d = new Display3DFrame("3D Frame");
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

		for (int i = 0; i < np.nN; i++) {
			nanoSphere[i].setSizeXYZ(2 * np.nanos[i].getrX(),
					2 * np.nanos[i].getrY(), 2 * np.nanos[i].getrZ());
			nanoSphere[i].getStyle().setFillColor(Color.BLACK);
			nanoSphere[i].setXYZ(np.nanos[i].getX(), np.nanos[i].getY(),
					np.nanos[i].getZ());
		}

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
		// Initialize files for writing output data
		if(np.mcs == 0 && writeMode != WriteModes.WRITE_NONE){
			DecimalFormat largeDecimal = new DecimalFormat("0.##E0");
			DecimalFormat threeDecimal = new DecimalFormat("#0.###");
			String configurations = "# Number of Polymers:" + np.nP +
					"\n# Number of Nanoparticles: "+np.nN +
					"\n# Move Tolerance: "+threeDecimal.format(np.tolerance)+
					"\n# Shape Change Tolerance: "+threeDecimal.format(np.shapeTolerance)+
					"\n# Nanoparticle Radius :"+threeDecimal.format(np.nano_r) + 
					"\n# Polymer Colloid Ratio: "+threeDecimal.format(np.q)+
					"\n# Lattice Constant: " +threeDecimal.format(np.lc)+
					"\n# Rotation Tolerance: "+threeDecimal.format(np.rotMagnitude)+
					"\n# Trial Moves to Attempt Shape Change Ratio: "+np.moveToShapeRatio+
					"\n# Snapshot Interval: "+largeDecimal.format(this.snapshotIntervals)+
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
			case WRITE_ALL:
				dataFiles = new DataFile[5];
				dataFiles[0] = new DataFile("eX", configurations);
				dataFiles[1] = new DataFile("eY", configurations);
				dataFiles[2] = new DataFile("eZ", configurations);
				dataFiles[3] = new DataFile("polar", configurations);
				dataFiles[4] = new DataFile("azimuth", configurations);
				break;
			default:
				break;
			}
		}
		
		np.step();

		for (int i = 0; i < np.nN; i++) {
			nanoSphere[i].setXYZ(np.nanos[i].getX(), np.nanos[i].getY(),
					np.nanos[i].getZ());
		}
		for (int i = 0; i < np.nP; i++) {
			polySphere[i].setXYZ(np.polymers[i].getX(), np.polymers[i].getY(),
					np.polymers[i].getZ());
			polySphere[i].setSizeXYZ(2 * np.polymers[i].getrX(),
					2 * np.polymers[i].getrY(), 2 * np.polymers[i].getrZ());
			if(visualizeOn && np.polymers[i].isRotated()){
				polySphere[i].setTransformation(np.polymers[i].getTransformation());
			}
		}

		if (writeMode != WriteModes.WRITE_NONE && np.mcs % snapshotIntervals == 0) {
			switch(writeMode){
			case WRITE_SHAPES:
				if(np.mcs > 50000){ // hardcoded 
					for(Polymer poly: np.polymers){
					dataFiles[0].record(String.valueOf(poly.geteX()));
					dataFiles[1].record(String.valueOf(poly.geteY()));
					dataFiles[2].record(String.valueOf(poly.geteZ()));
					}
				}
				break;
			case WRITE_ROTATIONS:
			for (Polymer poly : np.polymers) {
				double[] ellipseAxis = poly.getTransformAxis();
				control.println("(" + ellipseAxis[0] + ", " + ellipseAxis[1] + ", " + ellipseAxis[2] + ")");
				polar = ellipseAxis[2];
				azimuth = Math.atan(ellipseAxis[1]/ ellipseAxis[0]);
				dataFiles[0].record(String.valueOf(polar));
				dataFiles[1].record(String.valueOf(azimuth));
			} break;
			case WRITE_ALL:
					for(Polymer poly: np.polymers){
						if(np.mcs > 50000){ // hardcoded 
							dataFiles[0].record(String.valueOf(poly.geteX()));
							dataFiles[1].record(String.valueOf(poly.geteY()));
							dataFiles[2].record(String.valueOf(poly.geteZ()));
						}
						double[] ellipseAxis = poly.getTransformAxis();
						control.println("(" + ellipseAxis[0] + ", " + ellipseAxis[1] + ", " + ellipseAxis[2] + ")");
						polar = ellipseAxis[2];
						azimuth = Math.atan(ellipseAxis[1]/ ellipseAxis[0]);
						dataFiles[3].record(String.valueOf(polar));
						dataFiles[4].record(String.valueOf(azimuth));
					}
			break;
			default:break;
			}
		}
		
		// write data out when there's 100 data values
		if(writeMode != WriteModes.WRITE_NONE && np.mcs % (100*snapshotIntervals) == 0){
			for(DataFile df : dataFiles){
				df.write();
			}
		}

		plotframe.append(0, np.mcs, np.totalIntersectCount);
		totalIntersections += np.totalIntersectCount;
		display3d.setMessage("Number of mcs steps: " + np.mcs);
	}

	/**
	 * Resets the colloid polymer mixture model to its default state.
	 */
	public void reset() {
		enableStepsPerDisplay(true);
		control.setValue("N Polymers", 64);
		control.setValue("N Nano", 0);
		control.setValue("tolerance", 0);
		control.setValue("Shape Tolerance", 0);
		control.setValue("Nanoparticle radius", 0.1);
		control.setValue("x", 0.005);
		control.setValue("y", 0.005);
		control.setValue("z", 0.005);
		control.setValue("Polymer Colloid Ratio", 5);
		control.setValue("Lattice constant", 10);
		control.setValue("initial configuration", "square");
		control.setValue("Rotation magnitude", 5);
		control.setValue("Trial Moves to Shape Changes Ratio", 1);
		control.setAdjustableValue("Visualization", false);
		control.setValue("Snapshot Interval", 1000);
		control.setValue("Penetration Energy", false);
		control.setValue("Write Mode", 2);
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
		control.println("Vol fraction: " + volFract);
		
		if(writeMode != WriteModes.WRITE_NONE){
			for(DataFile df : dataFiles){
				df.close();
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
