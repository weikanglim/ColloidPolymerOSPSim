package org.opensourcephysics.sip.CPM;

import java.awt.Color;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Date;

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
	CPM np = new CPM();
	Display3DFrame display3d = new Display3DFrame("3D Frame");
	PlotFrame plotframe = new PlotFrame("Monte Carlo Steps",
			"Number of Intersections", "Number of Intersections");
	double totalIntersections = 0;
	double snapshotIntervals = 1;
	ElementSphere nanoSphere[];
	ElementEllipsoid polySphere[];
	boolean added = false;
	BufferedWriter bw1;
	BufferedWriter bw2;
	BufferedWriter bw3;
	Date date = new Date();
	Path f1;
	Path f2;
	Path f3;
	Path dir;


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
		np.moveToShapeRatio = control
				.getInt("Trial Moves to Shape Changes Ratio");
		snapshotIntervals = control.getInt("Snapshot Interval");
		np.initialize(configuration);
		if(display3d != null) display3d.dispose();
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
		if(np.mcs == 0 && snapshotIntervals > 0){
			dir = Paths.get("data");
			f1 = Paths.get("data/" + date + " x .dat");
			f2 = Paths.get("data/" + date + " y .dat");
			f3 = Paths.get("data/" + date + " z .dat");
			
			try {
				if (Files.notExists(dir, LinkOption.values())) {
					Files.createDirectory(dir);
				}
				
				Charset ascii = Charset.forName("US-ASCII");
				StandardOpenOption append = StandardOpenOption.APPEND;
				if(f1.toFile().exists()){
					bw1 = Files.newBufferedWriter(f1, ascii, append);
				} else{
					bw1 = Files.newBufferedWriter(f1, ascii);
				}
				
				if(f2.toFile().exists()){
					bw2 = Files.newBufferedWriter(f2, ascii, append);
				} else{
					bw2 = Files.newBufferedWriter(f2, ascii);
				}
				
				if(f3.toFile().exists()){
					bw3 = Files.newBufferedWriter(f3, ascii, append);
				} else{
					bw3 = Files.newBufferedWriter(f3, ascii);
				}
			} catch (IOException e) {
				e.printStackTrace();
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
		}

		if (snapshotIntervals > 0 && np.mcs % snapshotIntervals == 0 && np.mcs >= 50000) {
			for (Polymer poly : np.polymers) {
				try {
					bw1.write(String.valueOf(poly.geteX()) + "\n");
					bw2.write(String.valueOf(poly.geteY()) + "\n");
					bw3.write(String.valueOf(poly.geteZ()) + "\n");
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
		
		if(snapshotIntervals > 0 && np.mcs % 100*snapshotIntervals == 0 && np.mcs >= 50000){
			try {
				bw1.flush();
				bw2.flush();
				bw3.flush();
			} catch (IOException e) {
				e.printStackTrace();
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
		control.setValue("tolerance", 0.1);
		control.setValue("Shape Tolerance", 0.001);
		control.setValue("Nanoparticle radius", 0.01);
		control.setValue("x", 0.005);
		control.setValue("y", 0.005);
		control.setValue("z", 0.005);
		control.setValue("Polymer Colloid Ratio", 5);
		control.setValue("Lattice constant", 10);
		control.setValue("initial configuration", "square");
		control.setValue("Trial Moves to Shape Changes Ratio", 1);
		control.setValue("Snapshot Interval", 0);
		initialize();
	}

	public void stop() {
		double averageIntersections = totalIntersections / np.mcs;
		
		double volSpheres = np.nP
				* (4 / 3d * Math.PI * np.polymers[0].geteX()
						* np.polymers[0].geteY() * np.polymers[0].geteZ())
				* np.nN;
		double volFract = volSpheres / (np.Lx * np.Ly * np.Lz);
		control.println("Average no. of Intersections: " + averageIntersections);
		control.println("Vol fraction: " + volFract);
		
		if(snapshotIntervals > 0){
			try {
				bw1.close();
				bw2.close();
				bw3.close();
			} catch (IOException e) {
				e.printStackTrace();
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
