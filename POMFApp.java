package org.opensourcephysics.sip.CPM;

import java.awt.Color;
import java.text.DecimalFormat;

import org.opensourcephysics.controls.AbstractSimulation;
import org.opensourcephysics.controls.SimulationControl;
import org.opensourcephysics.display.Dataset;
import org.opensourcephysics.display3d.simple3d.ElementEllipsoid;
import org.opensourcephysics.display3d.simple3d.ElementPoints;
import org.opensourcephysics.display3d.simple3d.ElementSphere;
import org.opensourcephysics.frames.Display3DFrame;
import org.opensourcephysics.frames.DisplayFrame;
import org.opensourcephysics.frames.PlotFrame;
import org.opensourcephysics.numerics.Matrix3DTransformation;
import org.opensourcephysics.numerics.PBC;

/**
 * NanoPolyMixApp is a simulation framework for a binary mixture of
 * colloids(nanoparticles) and polymers model.
 * 
 * @author Wei Kang Lim, Alan Denton
 * @version 1.0 11-10-2013
 * 
 */
public class POMFApp extends AbstractSimulation {
	final int MCS_WAIT_TO_EQUIL = 5000;
	public enum WriteModes{WRITE_NONE,WRITE_SHAPES,WRITE_ROTATIONS,WRITE_RADIAL,WRITE_POMF,WRITE_ALL;};
	CPM np = new CPM();
	Display3DFrame display3d = new Display3DFrame("3D Frame");
	Display3DFrame display3d2 = new Display3DFrame("3D Frame");
	PlotFrame plotframe = new PlotFrame("Monte Carlo Steps",
			"e^-U", "Acceptance Probability Plot");
	DisplayFrame resultsFrame = new DisplayFrame("Radial Distance",
			"Effective Potential", "Effective Potential Plot");
	Dataset data;
	String configuration;
	String insertionType;
	RDF rdf;

	double totalIntersections = 0;
	double snapshotIntervals = 1;
	double [] zAxis = {0,0,1};
	double [] xAxis = {1,0,0};
	double [][][] radialData; 
	int runs;
	int currentRun;
	double polar;
	double azimuth;
	double volFraction;
	double radialStart; 
	double radialEnd;
	double sumVolume;
	double volumeSnapshots;
	double sumDistribution;
	double sumSquaredDistribution;
	long timeStarted = 0;
	long timeElapsed = 0;
	int conformations;
	int maxConformations;
	int dataPoints;
	int maxDataPoints;
	double waitMCS;
	double totalMCS;
	double steps;
	double placementPosition;
	int i =0;
	WriteModes writeMode;
	ElementSphere nanoSphere[];
	ElementEllipsoid polySphere[];
	ElementSphere closest;
	ElementSphere sphere;
	ElementEllipsoid ellipse;
	boolean added = false;
	boolean penetrationEnergyToggle;
	boolean clearNano = false;
	boolean debug = false;
	String minuteInfo = "";
	DataFile [] dataFiles;

	public POMFApp(){
		if(display3d != null) display3d.dispose();
		display3d = new Display3DFrame("3D Frame");

		if(display3d2 != null) display3d2.dispose();
		display3d2 = new Display3DFrame("3D Frame");
	}
	
	public void clearCounters(){
		sumDistribution = 0;
		sumVolume = 0;
		volumeSnapshots = 0;
		totalIntersections = 0;
		dataPoints = 0;
		conformations = 0;
		clearNano =false;
	}
	
	public void getInput(){
		// CPM (Simulation logic) input
		np.insertionType =insertionType = control.getString("Insertion method");
		np.q = control.getDouble("Polymer colloid ratio");
		np.Lx = np.Ly = np.Lz = control.getDouble("Lattice length");
		configuration = control.getString("Initial configuration");
		np.tolerance = control.getDouble("Tolerance");
		np.trialMovesPerMcs = control.getInt("Trial moves per MCS");
		Polymer.setExact(control.getBoolean("Exact overlap"));
		
		if(control.getBoolean("Spherical polymers")){
			np.init_eY = np.init_eZ =np.init_eX = 1/18f;
			np.rotTolerance = 0;
			np.shapeTolerance = 0;
		} else{
			np.init_eX = control.getDouble("x");
			np.init_eY = control.getDouble("y");
			np.init_eZ = control.getDouble("z");
			np.rotTolerance = control.getDouble("Rotation tolerance");
			np.shapeTolerance = control.getDouble("Shape tolerance");
		}

		
		// Simulation control input
		snapshotIntervals = control.getInt("Snapshot interval");
		maxConformations = control.getInt("Number of conformations");
		maxDataPoints = control.getInt("Number of datapoints"); 
		penetrationEnergyToggle =control.getBoolean("Penetration energy");
		runs = control.getInt("Runs");
	}
	/**
	 * Initializes the model.
	 */
	public void initialize() {
		i++;
		waitMCS = MCS_WAIT_TO_EQUIL;
		currentRun = 0;
		np.nano_r = 0.5;
		clearCounters();
		getInput();
		if(np.q < 1){
			radialEnd = Math.min(1 + np.q + 0.1, 1+3*np.q) ;
		} else {
			radialEnd = Math.min(1 + np.q + 0.1, (np.Lx-1)/2) ; // 2*Rp+2*Rn
		}
		radialStart = 1;
		steps = (radialEnd-radialStart) / maxDataPoints; // calculate dr needed to iterate through from [radialEnd, radialStart]
		maxDataPoints++; // increase datapoint by 1 to account for 1 extra datapoint for run at r = 0
		System.out.println(radialStart + " " + radialEnd + " "  + " by " + steps );
		if(insertionType.equals("polymer")){
			totalMCS = runs* maxConformations * maxDataPoints;
		} else{
			totalMCS = runs* maxConformations * snapshotIntervals * maxDataPoints;			
		}
		radialData = new double[runs][maxDataPoints+1][3]; // radialData[i][0] = r, radialData[i][1] = e^[-U(r)], radialData[i][2] = uncertainty

		switch(control.getInt("Write Mode")){
		case 0: writeMode = WriteModes.WRITE_NONE; break;
		case 1: writeMode = WriteModes.WRITE_SHAPES; break;
		case 2: writeMode = WriteModes.WRITE_ROTATIONS; break;
		case 3: writeMode = WriteModes.WRITE_RADIAL; break;
		case 4: writeMode = WriteModes.WRITE_POMF; break;
		case 5: writeMode = WriteModes.WRITE_ALL; break;
		}
		
		np.initialize(configuration, penetrationEnergyToggle);
		placementPosition = radialStart;
		if(insertionType.equals("polymer")) np.placeNano2(placementPosition);
		plotframe.setMessage("r = " + placementPosition);

		// add simple3d.Element particles to the arrays 
		if (i == 2) { // particles only allowed to be added once, this is to prevent clone particles when initialize is called twice by the simulation
			nanoSphere = new ElementSphere[np.nN];
			polySphere = new ElementEllipsoid[np.nP];
			closest = new ElementSphere();
			closest.setSizeXYZ(0.1,0.1,0.1);
			closest.getStyle().setFillColor(Color.BLACK); 
			display3d.addElement(closest);


			for (int i = 0; i < np.nN; i++) {
				nanoSphere[i] = new ElementSphere();
				display3d.addElement(nanoSphere[i]);
			}

			for (int i = 0; i < np.nP; i++) {
				if(control.getBoolean("Spherical polymers")){
					polySphere[i] = new ElementSphere();
				}else{
					polySphere[i] = new ElementEllipsoid();
				}
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
//	            if(np.polymers[i].updateRotation()){
//	                    polySphere[i].setTransformation(np.polymers[i].getRotationTransformation());
//	            }
			}
			 

			display3d.setPreferredMinMax(0, np.Lx, 0, np.Ly, 0, np.Lz);
			display3d.setSquareAspect(true);
		}		
	}
	
	/**
	 * Does a simulation step.
	 */
	public void doStep() {
		if(np.mcs == 0){ // Initialization
			initializeDataFiles();
			timeStarted = System.nanoTime();
		}

		// logical step in the CPM class
		np.step();
		
		// Insertion algorithms
		// Polymer Insertion
		if(insertionType.equals("polymer")){
			double e_negU = np.polyTrialPlacement(Math.random()*np.Lx, Math.random()*np.Ly, Math.random()*np.Lz);
			
			if(Polymer.getShapeTolerance() == 0 || np.mcs >= 100){
				if(debug){
					if(np.polymers[0].debug){ // overlap occurred
						// Update our visualization for debugging purposes
						for (int i = 0; i < np.nN; i++) {
							nanoSphere[i].setXYZ(np.nanos[i].getX(), np.nanos[i].getY(),
									np.nanos[i].getZ());
						}
						
						for (int i = 0; i < np.nP; i++) {
							polySphere[i].setXYZ(np.polymers[i].getX(), np.polymers[i].getY(),
									np.polymers[i].getZ());
							polySphere[i].setSizeXYZ(2 * np.polymers[i].getrX(),
									2 * np.polymers[i].getrY(), 2 * np.polymers[i].getrZ());
							if(np.polymers[i].updateRotation()){ // this is done to save computation
								polySphere[i].setTransformation(np.polymers[i].getRotationTransformation());
							}
						}
	
						closest.setXYZ(np.polymers[0].closestPoint[0], 
								np.polymers[0].closestPoint[1],
								np.polymers[0].closestPoint[2]);
						if(ellipse == null || sphere == null){
							ellipse = new ElementEllipsoid();
							ellipse.getStyle().setFillColor(Color.RED);
							ellipse.getStyle().setLineColor(Color.RED);
	
							sphere = new ElementSphere();
							sphere.setRadius(np.nanos[0].getrX());
							sphere.getStyle().setFillColor(new Color(92, 146, 237, 100));  // light blue with half transparency
							sphere.getStyle().setLineColor(new Color(92, 146, 237, 100));
	
							display3d2.addElement(ellipse);
							display3d2.addElement(sphere);
							display3d2.setPreferredMinMax(-np.q, np.q, -np.q, np.q, -np.q, np.q);
							display3d2.setSquareAspect(true);
						}
						ellipse.setSizeXYZ(2 * np.polymers[0].getrX(),
								2 * np.polymers[0].getrY(), 2 * np.polymers[0].getrZ());
						ellipse.setXYZ(0, 0, 0);
						sphere.setXYZ(np.polymers[0].overlapSphere[0], np.polymers[0].overlapSphere[1], np.polymers[0].overlapSphere[2]);
						display3d2.render();
						display3d.render();
						System.out.println("Pause debug");
					} else {
						closest.setXYZ(0,0,0);
					}
				}
				
				
				sumDistribution += e_negU;
				plotframe.append(0, np.mcs, e_negU);
				conformations++;
			}
		}
		
		// Nanoparticle insertion
		if(insertionType.equals("nano") ){
			if(waitMCS > 0){
				waitMCS--;
			} else if( np.mcs % snapshotIntervals == 0){
				double e_negU = np.nanoTrialPlacement(placementPosition);
				sumDistribution += e_negU;
				plotframe.append(0, np.mcs, e_negU);
				conformations++;
			}
		}
		
		// Insertion algorithm snapshots
		if (writeMode != WriteModes.WRITE_NONE && conformations >= maxConformations) {		
			// Enough conformations for a data point, analyze distribution and record.
			double avgDistribution = sumDistribution / conformations;
			System.out.println(avgDistribution);
			radialData[currentRun][dataPoints][0] = placementPosition;
			radialData[currentRun][dataPoints][1] = avgDistribution; // temporarily place the variable, will be replaced with V_r calculations later on.
			dataPoints++;
			if(dataPoints >= maxDataPoints +1 ){ // all data points + U(inf)				
				System.out.println("Run " + currentRun);

				// Calculate V_r from averaged e^-U
				int U_zero_index = maxDataPoints; // last position
				double U_inf = 2*radialData[currentRun][U_zero_index][1]-1;
				double lnU_inf = Math.log(U_inf);
				System.out.println(lnU_inf + " U_inf: " + U_inf);
				for(int i =0; i < maxDataPoints; i++){
					double r = radialData[currentRun][i][0];
					double lnU_r = Math.log(radialData[currentRun][i][1]);
					System.out.println(lnU_r);
					double V_r = -lnU_r + lnU_inf;
					radialData[currentRun][i][1] = V_r;
					System.out.println(r + "\t" + V_r);
				}
				
				// Reset counters for next run
				clearCounters();
				if(insertionType.equals("polymer")){
					np.nN = 2;
					np.nanos = new Nano[2];
					np.nanos[0] = new Nano(np.Lx/2f, np.Ly/2f, np.Lz/2f);
				} else {
					np.nN = 1;
					np.nanos = new Nano[1];
					np.nanos[0] = new Nano(np.Lx/2f, np.Ly/2f, np.Lz/2f);
				}
				currentRun++;

				if(currentRun >= runs){
					double [] avgPotential = new double[maxDataPoints];
					double [] stdDevPotential = new double[maxDataPoints];
					for(int i = 0; i < maxDataPoints; i++){
						for(int j = 0; j < runs; j++){
							avgPotential[i] += radialData[j][i][1];
							stdDevPotential[i] += Math.pow(radialData[j][i][1],2);
						}
						
						
						avgPotential[i] /= runs;
						stdDevPotential[i] /= runs;
						stdDevPotential[i] = Math.sqrt(stdDevPotential[i] - Math.pow(avgPotential[i],2));
					}
					
					data = new Dataset();
					dataFiles[0].record("# r\tV_r\tStd-dev" );
					for(int i = 0; i < maxDataPoints; i++){
						dataFiles[0].record(radialData[0][i][0] + "\t" + avgPotential[i] + "\t" + stdDevPotential[i]);
						data.append(radialData[0][i][0], avgPotential[i], 0, stdDevPotential[i]);
					}

					resultsFrame.addDrawable(data);
					control.setAdjustableValue("Save", true);
					control.clearMessages();
					control.println("Runs completed.");
//					double avgPhiP = (sumVolume / volumeSnapshots) / (np.nP*np.Lx*np.Ly*np.Lz);
//					System.out.println(avgPhiP);
//					DecimalFormat threeDecimals = new DecimalFormat("#0.###");
//					String phiP = threeDecimals.format(avgPhiP);
//					dataFiles[0].record("Average phiP: " +  phiP);
					int elapsedMinutes = (int) Math.floor(timeElapsed/(1000*60));
					int elapsedSeconds = (int) Math.round(timeElapsed/1000) % 60;
					String formatTimeElapsed = (elapsedMinutes == 0) ? elapsedSeconds + "s ": elapsedMinutes + "m " + elapsedSeconds + "s"; 
					dataFiles[0].record("# Total simulation time: " + formatTimeElapsed); 
					dataFiles[0].record("# Number of runs: " + currentRun);
					this.stopAnimation();
					return;
				}
			}
			
			// reset counters for next data point.
			conformations = 0;
			sumDistribution = 0;			
			// set placement position to be 0 to calculate U at inf for last run, otherwise perform increment in radial distance from radialStart by step
			placementPosition = dataPoints == maxDataPoints ? 0 : radialStart+dataPoints*steps;
			
			if(placementPosition == 0 && !clearNano){ 
				if(insertionType.equals("polymer")){
					clearNano = true; // only perform once
					nanoSphere[1].setVisible(false);
					plotframe.clearDataAndRepaint();
					
					// Keep only one nanoparticle at the center
					np.nN = 1;
					np.nanos = new Nano[1];
					np.nanos[0] = new Nano(np.Lx/2f, np.Ly/2f, np.Lz/2f);
				} else {
					clearNano = true;
					nanoSphere[0].setVisible(false);
					plotframe.clearDataAndRepaint();
					np.nN = 0;
					np.nanos = new Nano[0];
					waitMCS = MCS_WAIT_TO_EQUIL;
				}
			} 

			// set new nanoparticle position
			if(insertionType.equals("polymer")) np.placeNano2(placementPosition);
			plotframe.setMessage("r = " + placementPosition);
		}		
		
//		if(writeMode != WriteModes.WRITE_NONE && np.mcs % (snapshotIntervals*10) == 0){
//			dataFiles[0].write();
//		}
		

		//!-- Visualization --! // START
		// update intersect count and mcs step
		display3d.setMessage("Number of mcs steps: " + np.mcs);

		// Visualization updates
		if(control.getBoolean("Visualization on")){
			// nanoparticle visualization updates
			for (int i = 0; i < np.nN; i++) {
				nanoSphere[i].setXYZ(np.nanos[i].getX(), np.nanos[i].getY(),
						np.nanos[i].getZ());
			}
			// polymer visualization updates
			for (int i = 0; i < np.nP; i++) {
				polySphere[i].setXYZ(np.polymers[i].getX(), np.polymers[i].getY(),
						np.polymers[i].getZ());
				polySphere[i].setSizeXYZ(2 * np.polymers[i].getrX(),
						2 * np.polymers[i].getrY(), 2 * np.polymers[i].getrZ());
				if(np.polymers[i].updateRotation()){ // this is done to save computation
					polySphere[i].setTransformation(np.polymers[i].getRotationTransformation());
				}
			}
		}
		// Simulation info
		timeElapsed = (System.nanoTime() - timeStarted)/1000000;
		if(timeElapsed % 1000 == 0){
			control.clearMessages();
			int elapsedMinutes = (int) Math.floor(timeElapsed/(1000*60));
			int elapsedSeconds = (int) Math.round(timeElapsed/1000) % 60;
			String formatTimeElapsed = (elapsedMinutes == 0) ? elapsedSeconds + "s ": elapsedMinutes + "m " + elapsedSeconds + "s"; 
			control.println("Time Elapsed: " + formatTimeElapsed);
			control.println(minuteInfo);
			
			if(timeElapsed % 60000 == 0){
				DecimalFormat largeDecimal = new DecimalFormat("0.##E0");
				double mcsPerMinute = 60000* np.mcs / timeElapsed;
				double timeRemain = (totalMCS - np.mcs) / mcsPerMinute;
				int minutes = (int) Math.floor(timeRemain);
				int seconds = (int) Math.round((timeRemain - minutes) * 60); 
				minuteInfo = "MCS per minute: " + largeDecimal.format(mcsPerMinute)
						+"\nTime remaining: " + minutes + "m " + seconds + "s";
			}
		}
		//!-- Visualization --! // END		
	}

	/**
	 * Resets the colloid polymer mixture model to its default state.
	 */
	public void reset() {
		enableStepsPerDisplay(true);
		control.setValue("Polymer colloid ratio", 0.776);
		control.setValue("Spherical polymers", false);
		control.setValue("Lattice length", Math.cbrt(Math.PI/6*1/0.01)); // Testing 2 + 2q for q < 1
		control.setValue("x", 0.01);
		control.setValue("y", 0.01);
		control.setValue("z", 0.01);
		control.setValue("Runs", 5);
		control.setValue("Tolerance", 0);
		control.setValue("Rotation tolerance", 0.1);
		control.setValue("Shape tolerance", 0.001);
		control.setValue("Insertion method", "polymer");
		control.setValue("Exact overlap", true);
		control.setValue("Initial configuration", "square");
		control.setValue("Trial moves per MCS", 1);
		control.setAdjustableValue("Visualization on", true);
		control.setValue("Snapshot interval", 1000);
		control.setValue("Number of datapoints", 8);
		control.setValue("Number of conformations", 2000000);
		control.setValue("Penetration energy", true);
		control.setValue("Write Mode", 4);
		control.setAdjustableValue("Save", false);
		i=0;
		initialize();
	}

	public void stop() {
		double averageIntersections = totalIntersections / np.mcs;

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
					df.close();
				}
			}
		}
	}
	
	public void initializeDataFiles(){
		DecimalFormat largeDecimal = new DecimalFormat("0.##E0");
		DecimalFormat threeDecimal = new DecimalFormat("#0.###");
		DecimalFormat fourDecimal = new DecimalFormat("#0.####");
		String phi_n;
		if(np.volFraction < 0.0001){
			phi_n = largeDecimal.format(np.volFraction);
		} else if(np.volFraction < 0.001){
			phi_n = fourDecimal.format(np.volFraction);
		} else {
			phi_n = threeDecimal.format(np.volFraction);
		}
		String configurations = "# Number of Polymers: " + np.nP +
				"\n# Number of Nanoparticles: "+np.nN +
				"\n# Move Tolerance: "+threeDecimal.format(np.tolerance)+
				"\n# Shape Change Tolerance: "+threeDecimal.format(np.shapeTolerance)+
				"\n# Insertion type: " + insertionType +
				"\n# Nanoparticle Volume Fraction: "+phi_n + 
				"\n# Polymer Colloid Ratio: "+threeDecimal.format(np.q)+
				"\n# Exact overlap: " + Polymer.getExact() + 
				"\n# Lattice Length: " +threeDecimal.format(np.Lx)+
				"\n# Rotation Tolerance: "+threeDecimal.format(np.rotTolerance)+
				"\n# Trial Moves Per Mcs: "+np.trialMovesPerMcs+
				"\n# Snapshot Interval: "+largeDecimal.format(this.snapshotIntervals)+
				"\n# Number of Coformations Sampled: " + maxConformations +
				"\n# Number of dataPoints: " + (maxDataPoints-1) + // 1 datapoint is used for e^(-U(r=0))
				"\n# Penetration Energy On: " + this.penetrationEnergyToggle
				;
		switch(writeMode){
		case WRITE_SHAPES:
			dataFiles = new DataFile[3];
			dataFiles[0] = new DataFile("eX", configurations, DataFile.FileIdentifier.SIZE_AND_FRACTION, true);
			dataFiles[1] = new DataFile("eY", configurations, DataFile.FileIdentifier.SIZE_AND_FRACTION, true);
			dataFiles[2] = new DataFile("eZ", configurations, DataFile.FileIdentifier.SIZE_AND_FRACTION, true);
			break;
		case WRITE_ROTATIONS:
			dataFiles = new DataFile[2];
			dataFiles[0] = new DataFile("polar", configurations, DataFile.FileIdentifier.SIZE_AND_FRACTION, true);
			dataFiles[1] = new DataFile("azimuth", configurations, DataFile.FileIdentifier.SIZE_AND_FRACTION, true);
			break;
		case WRITE_RADIAL:
			dataFiles = new DataFile[1]; 
			dataFiles[0] = new DataFile("radial", configurations, DataFile.FileIdentifier.SIZE_AND_FRACTION, false);
			rdf = new RDF(np.nanos, np.Lx);
			break;
		case WRITE_POMF:
			dataFiles = new DataFile[1];
			dataFiles[0] = new DataFile("POMF", configurations, DataFile.FileIdentifier.SIZE_AND_FRACTION, false);
			break;
		case WRITE_ALL:
			dataFiles = new DataFile[4];
			dataFiles[0] = new DataFile("eX", configurations, DataFile.FileIdentifier.SIZE_AND_FRACTION, true);
			dataFiles[1] = new DataFile("eY", configurations, DataFile.FileIdentifier.SIZE_AND_FRACTION, true);
			dataFiles[2] = new DataFile("eZ", configurations, DataFile.FileIdentifier.SIZE_AND_FRACTION, true);
			dataFiles[3] = new DataFile("radial", configurations, DataFile.FileIdentifier.SIZE_AND_FRACTION, true);
			rdf = new RDF(np.nanos, np.Lx);
			break;
		default:
			break;
		}
		
		if(writeMode != WriteModes.WRITE_NONE){
		// Write out the initialization data
			for(DataFile df : dataFiles){
				df.write();
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
		SimulationControl control = SimulationControl.createApp(new POMFApp());
	}
}
