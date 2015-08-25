package org.opensourcephysics.sip.CPM;

import java.awt.Color;
import java.text.DecimalFormat;

import org.opensourcephysics.controls.AbstractSimulation;
import org.opensourcephysics.controls.SimulationControl;
import org.opensourcephysics.display.Dataset;
import org.opensourcephysics.display3d.simple3d.ElementEllipsoid;
import org.opensourcephysics.display3d.simple3d.ElementSphere;
import org.opensourcephysics.frames.Display3DFrame;
import org.opensourcephysics.frames.DisplayFrame;
import org.opensourcephysics.frames.PlotFrame;

/**
 * NanoPolyMixApp is a simulation framework for a binary mixture of
 * colloids(nanoparticles) and polymers model.
 * 
 * @author Wei Kang Lim, Alan Denton
 * 
 */
public class POMFApp extends AbstractSimulation {
	final int MCS_WAIT_TO_EQUIL = 5000;
	final double LAMBDA1_END = 0.4;
	final double LAMBDA2_END = 0.08;
	final double LAMBDA3_END = 0.03;
	final int HISTOGRAM_BINS = 1000;
	
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
	double [] lnP_inf;
	Histogram lambda1 = new Histogram(0.0, LAMBDA1_END, HISTOGRAM_BINS);
	Histogram lambda2 = new Histogram(0.0, LAMBDA2_END, HISTOGRAM_BINS);
	Histogram lambda3 = new Histogram(0.0, LAMBDA3_END, HISTOGRAM_BINS);
	// radialData[i][0] = r, radialData[i][1] = v(r), radialData[i][2] = uncertainty of v(r)
	int runs;
	int currentRun;
	double polar;
	double azimuth;
	double volFraction;
	double radialStart; 
	double radialEnd;
	double sumVolume;
	double volumeSnapshots;
	double sumInteractionEnergy;
	double sumLn_Pr;
	double sumSquaredInteraction;
	double r_U_inf;
	long timeStarted = 0;
	long timeElapsed = 0;
	int conformations;
	int maxConformations;
	int dataPoints;
	int userDataPoints; // user specified data points
	int simDataPoints; // simulation required data points.
	int countShape;
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
	boolean bruteForce = false;
	boolean added = false;
	boolean penetrationEnergyToggle;
	boolean debug = false;
	String minuteInfo = "";
	DataFile [] dataFiles;

	public void clearCounters(){
		lambda1.clear();
		lambda2.clear();
		lambda3.clear();
		sumInteractionEnergy = 0;
		sumLn_Pr = 0;
		sumVolume = 0;
		volumeSnapshots = 0;
		totalIntersections = 0;
		dataPoints = 0;
		conformations = 0;
		countShape = 0;
	}
	
	public void getInput(){
		// CPM (Simulation logic) input
		np.insertionType =insertionType = control.getString("Insertion method");
		np.q = control.getDouble("Polymer colloid ratio");
		np.Lx = np.Ly = np.Lz = control.getDouble("Lattice length");
		
		if(control.getBoolean("Auto width")){
			np.Lz = np.Lx = 1 + 3*np.q;
			np.Ly = 2 + 4*np.q;
		}
		
		configuration = control.getString("Initial configuration");
		np.tolerance = control.getDouble("Tolerance");
		np.trialMovesPerMcs = control.getInt("Trial moves per MCS");
		np.C = control.getDouble("Penetration free parameter:");
		np.energyProfile = control.getBoolean("Energy profile");
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
		userDataPoints = control.getInt("Number of datapoints"); 
		simDataPoints = userDataPoints + 2; // Includes runs for U_inf and U_inf_shape
		penetrationEnergyToggle =control.getBoolean("Penetration energy");
		runs = control.getInt("Runs");
		bruteForce = control.getBoolean("U_inf bruteforce");
		switch(control.getInt("Write Mode")){
		case 0: writeMode = WriteModes.WRITE_NONE; break;
		case 1: writeMode = WriteModes.WRITE_SHAPES; break;
		case 2: writeMode = WriteModes.WRITE_ROTATIONS; break;
		case 3: writeMode = WriteModes.WRITE_RADIAL; break;
		case 4: writeMode = WriteModes.WRITE_POMF; break;
		case 5: writeMode = WriteModes.WRITE_ALL; break;
		}
	}
	/**
	 * Initializes the model.
	 */
	public void initialize() {
		if(display3d != null) display3d.dispose();
		if(display3d2 != null) display3d2.dispose();
		if(plotframe != null) plotframe.dispose();
		if(resultsFrame != null) resultsFrame.dispose();
		
		if(debug) display3d2 = new Display3DFrame("3D Frame");
		plotframe = new PlotFrame("Monte Carlo Steps",
				"e^-U", "Acceptance Probability Plot");
		resultsFrame = new DisplayFrame("Radial Distance",
				"Effective Potential", "Effective Potential Plot");
		display3d = new Display3DFrame("3D Frame");


		waitMCS = MCS_WAIT_TO_EQUIL;
		currentRun = 0;
		np.nano_r = 0.5;
		getInput();
		clearCounters();

		// Setup insertion parameters.
		radialStart = 1;
		// Place the second nanoparticle up to 1 + q, when V(r) -> 0, or the maximum box length it can go before PBC kicks in.
		if(np.q >= 0.5){
			radialEnd = Math.min(1 + np.q + 0.1, np.Ly/2 - 0.5 - np.q);  // Nanoparticle placed at center.
		} else {
			radialEnd = Math.min(1 + np.q + np.q/2, np.Ly - 1 - 3 * np.q ); // The nanoparticle is not placed at the center for q < 0.5
		}
				
		steps = (radialEnd-radialStart) / (userDataPoints-1); // calculate dr needed to iterate through from [radialEnd, radialStart]
		System.out.println(radialStart + " " + radialEnd + " "  + " by " + steps );
		if(insertionType.equals("polymer")){
			totalMCS = runs* maxConformations * simDataPoints; 
		} else{
			totalMCS = runs* maxConformations * snapshotIntervals * userDataPoints;			
		}
		radialData = new double[runs][simDataPoints][3];
		lnP_inf = new double[runs];
		
		np.initialize(configuration, penetrationEnergyToggle);
		placementPosition = radialStart;
		
		if(bruteForce) r_U_inf = 4 * ( 1 + np.q);
		if(insertionType.equals("polymer")) np.placeNano2(placementPosition);
		plotframe.setMessage("r = " + placementPosition);

		// Initialize Visualization
		// add simple3d.Element particles to the arrays 
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
			display3d.setPreferredMinMax(0, np.Lx, 0, np.Ly, 0, np.Lz);
			display3d.setSquareAspect(true);
		}		
	}
	
	/**
	 * Does a simulation step.
	 */
	@SuppressWarnings("deprecation")
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
				
				// Record polymer internal conformations
				if(e_negU < 1){ // delta U > 0, energy unfavourable
					 // Accept with boltzmann factor
					if(Math.random() < e_negU){
						// Accepted.
						// Record polymer shape.
						lambda1.add(np.polymers[0].geteX());
						lambda2.add(np.polymers[0].geteY());
						lambda3.add(np.polymers[0].geteZ());
					}
					
					// Rejected, do nothing.
				} else { 
					// Accepted.
					// Record polymer shape.
					lambda1.add(np.polymers[0].geteX());
					lambda2.add(np.polymers[0].geteY());
					lambda3.add(np.polymers[0].geteZ());
				}
				
				// Record polymer-nanoparticle interaction
				sumInteractionEnergy += e_negU;
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
				sumInteractionEnergy += e_negU;
				plotframe.append(0, np.mcs, e_negU);
				conformations++;
			}
		}
		
		// Insertion algorithm snapshots
		if (writeMode != WriteModes.WRITE_NONE && conformations >= maxConformations) {		
			// Enough conformations for a data point, analyze distribution and record.
			double avgInteractionEnergy = sumInteractionEnergy / conformations;
			double lnP1_r = lambda1.calculateEntropy(),
				   lnP2_r = lambda2.calculateEntropy(),
				   lnP3_r = lambda3.calculateEntropy();
			double lnP_r = lnP1_r + lnP2_r + lnP3_r;
			System.out.println(placementPosition + "\t" + avgInteractionEnergy + "\t" + lnP_r);
			radialData[currentRun][dataPoints][0] = placementPosition;
			radialData[currentRun][dataPoints][1] = avgInteractionEnergy; // temporarily place the variable, will be replaced with V_r calculations later on.
			radialData[currentRun][dataPoints][2] = lnP_r; // temporary variable storage, will be replaced with uncertainty of V_r when a single run finishes.
			
			// Increase counter for next datapoint.
			dataPoints++;
			
			// For non-bruteforce calculations, we just need lnP_1.
			// Hence, the last datapoint (r_inf) run can be skipped.
			if(dataPoints == simDataPoints - 1 && !bruteForce){ 
				dataPoints++; // Skip the last datapoint run.
			}

			// Enough datapoints, start calculating V_r
			if(dataPoints >= simDataPoints ){ // all data points + U(inf)				
				System.out.println("Run " + currentRun + " completed.");
				System.out.println("Results:");

				// Get the infinite separation data.
				// Last datapoint stores data from one nanoparticle insertions
				double U_inf = 2*radialData[currentRun][simDataPoints-2][1]-1;
				double U_inf_shape;
				if(bruteForce){
					U_inf_shape = radialData[currentRun][simDataPoints-1][2];
				} else {
					double lnP_one = radialData[currentRun][simDataPoints-2][2];
					U_inf_shape = 2*lnP_one - 8.44002894726;
				}

				System.out.println("r\tV(r)\tf(r)_poly-nano\tf(r)_shape");
				
				// Calculate V_r from averaged e^-U 
				for(int i =0; i < userDataPoints; i++){
					double r = radialData[currentRun][i][0];
					double U_r = radialData[currentRun][i][1];
					
					double internalFree = U_inf_shape - radialData[currentRun][i][2];
					double V_r = U_inf - U_r + internalFree;
					radialData[currentRun][i][1] = V_r;
					lnP_inf[currentRun] = U_inf_shape;
					System.out.println(r + "\t" + V_r  + "\t" + (U_inf - U_r) + "\t" + internalFree);
				}
				System.out.println("ln_P_inf:" + U_inf_shape);
				System.out.println("U_inf: " + U_inf);
								
				// Reset counters for next run
				clearCounters();
				
				// Reset simulation parameters
				if(insertionType.equals("polymer")){
					np.nN = 2;
					np.nanos = new Nano[2];
					if(bruteForce){
						np.Ly -= r_U_inf;
						
						// Redraw display
						display3d.setPreferredMinMax(0, np.Lx, 0, np.Ly, 0, np.Lz);
						display3d.setSquareAspect(true);
					}
					np.setPolyInsertionPositions();
					
					nanoSphere[1].setVisible(true);
				} else {
					np.nN = 1;
					np.nanos = new Nano[1];
					nanoSphere[1].setVisible(true);
					np.nanos[0] = new Nano(np.Lx/2f, np.Ly/2f, np.Lz/2f);
				}
				currentRun++;

				if(currentRun >= runs){ // All runs completed.
					double [] avgPotential = new double[userDataPoints];
					double [] stdDevPotential = new double[userDataPoints];
					for(int i = 0; i < userDataPoints; i++){
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
					for(int i = 0; i < userDataPoints; i++){
						dataFiles[0].record(radialData[0][i][0] + "\t" + avgPotential[i] + "\t" + stdDevPotential[i]);
						data.append(radialData[0][i][0], avgPotential[i], 0, stdDevPotential[i]);
					}
					
					double avg_lnP_inf = 0;
					double stdDevLnP_inf = 0;
					
					for(int i=0; i < runs; i++){
						avg_lnP_inf += lnP_inf[i];
						stdDevLnP_inf += Math.pow(lnP_inf[i], 2);
					}
					
					avg_lnP_inf /= runs;
					stdDevLnP_inf /= runs;
					stdDevLnP_inf = Math.sqrt(stdDevLnP_inf - Math.pow(avg_lnP_inf, 2));

					resultsFrame.addDrawable(data);
					control.setAdjustableValue("Save", true);
					control.clearMessages();
					control.println("Runs completed.");
					int elapsedMinutes = (int) Math.floor(timeElapsed/(1000*60));
					int elapsedSeconds = (int) Math.round(timeElapsed/1000) % 60;
					String formatTimeElapsed = (elapsedMinutes == 0) ? elapsedSeconds + "s ": elapsedMinutes + "m " + elapsedSeconds + "s";
					dataFiles[0].comment("lnP_inf: " + avg_lnP_inf + "\t" + stdDevLnP_inf);
					dataFiles[0].comment("Total simulation time: " + formatTimeElapsed); 
					dataFiles[0].comment("Number of runs: " + currentRun);
					this.stopAnimation();
					return;
				} else{
					System.out.println("Run " + currentRun + " start.");
					System.out.println("r\t<e^-(del_U)>\t<ln P_2(r)>");
				}
			}
			
			// reset counters for next data point.
			conformations = 0;
			sumInteractionEnergy = 0;			
			
			// Set placement position for next datapoint.
			if(dataPoints == simDataPoints - 2){ // U_inf run
				placementPosition = 0;
				
				if(insertionType.equals("polymer")){
					nanoSphere[1].setVisible(false);
					plotframe.clearDataAndRepaint();
					
					// Keep only one nanoparticle at the center
					np.nN = 1;
					np.nanos = new Nano[1];
					np.setPolyInsertionPositions();
				} else {
					nanoSphere[0].setVisible(false);
					plotframe.clearDataAndRepaint();
					np.nN = 0;
					np.nanos = new Nano[0];
					waitMCS = MCS_WAIT_TO_EQUIL;
				}
			} else if(dataPoints == simDataPoints - 1 && bruteForce){ // U_inf_shape run
				np.nN = 2;
				np.nanos = new Nano[2];
				placementPosition = r_U_inf;
				
				np.Ly += r_U_inf; // Increase box length
				np.setPolyInsertionPositions(); // Place nanoparticle at new center.
								
				// Redraw display
				nanoSphere[1].setVisible(true);				
				display3d.setPreferredMinMax(0, np.Lx, 0, np.Ly, 0, np.Lz);
				display3d.setSquareAspect(true);

			} else{
				placementPosition = radialStart+(dataPoints)*steps;
			}
			
			// set new nanoparticle position
			if(insertionType.equals("polymer"))  np.placeNano2(placementPosition);
			plotframe.setMessage("r = " + placementPosition);
		}		
		
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
		control.setValue("Polymer colloid ratio", 10);
		control.setValue("Spherical polymers", false);
		control.setValue("Lattice length", Math.cbrt(Math.PI/6*1/0.0002)); // Testing 2 + 2q for q < 1
		control.setValue("x", 0.01);
		control.setValue("y", 0.01);
		control.setValue("z", 0.01);
		control.setValue("Penetration free parameter:", 6.839);
		control.setValue("Runs", 5);
		control.setValue("Tolerance", 0);
		control.setValue("Rotation tolerance", 0.1);
		control.setValue("Shape tolerance", 0.001);
		control.setValue("Insertion method", "polymer");
		control.setValue("Auto width", true);
		control.setValue("Exact overlap", true);
		control.setValue("Energy profile", true);
		control.setValue("U_inf bruteforce", false);
		control.setValue("Trial moves per MCS", 1);
		control.setAdjustableValue("Visualization on", true);
		control.setValue("Snapshot interval", 1);
		control.setValue("Number of datapoints", 8);
		control.setValue("Number of conformations", 2000000);
		control.setValue("Penetration energy", true);
		control.setValue("Write Mode", 4);
		control.setAdjustableValue("Save", false);
		i=0;
		initialize();
	}

	public void stop() {
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
		String [] configurations = { "Number of Polymers: " + np.nP,
				"Number of Nanoparticles: "+np.nN,
				"Move Tolerance: "+threeDecimal.format(np.tolerance),
				"Shape Change Tolerance: "+threeDecimal.format(np.shapeTolerance),
				"Insertion type: " + insertionType,
				"Nanoparticle Volume Fraction: "+np.volFraction, 
				"Polymer Colloid Ratio: "+threeDecimal.format(np.q),
				"Exact overlap: " + Polymer.getExact(),
				"x length: " +largeDecimal.format(np.Lx),
				"y length: " +largeDecimal.format(np.Ly),
				"z length: " +largeDecimal.format(np.Lz),
				"Rotation Tolerance: "+threeDecimal.format(np.rotTolerance),
				"Trial Moves Per Mcs: "+np.trialMovesPerMcs,
				"Snapshot Interval: "+largeDecimal.format(this.snapshotIntervals),
				"Number of Coformations Sampled: " + maxConformations,
				"Number of dataPoints: " + userDataPoints,
				"U_inf bruteforce: " + bruteForce,
				"Penetration Energy: " + (np.energyProfile? np.C + "/r" : np.C + "/q" + "="  + np.step_Ep)}
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
