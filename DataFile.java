package org.opensourcephysics.sip.CPM;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.text.DecimalFormat;
import java.util.Calendar;

/**
 * <p>This class handles the creations of data files, and writing of data. All data are stored in a temporary buffer before being written to file.</p>
 * 
 * <p>The data added to the DataFile via methods {@link #record(String) record}, {@link #comment(String) comment} and {@link #commentLines(String[]) commentLines}
 * are stored in a data buffer. The data can be written to file by using the {@link #write() write} method.</p>
 * 
 * <p>By default, the data files are stored in "$CWD/data". <br></br>
 * For data files being stored for multiple runs, they are stored under "$CWD/data/$RUN/$RunNumber" for each multiple run.
 * For example, "data/phiN=0.5,q=0.125_18-6-2014 23:44/1" for the first run.
 * </p>
 *  
 *  
 * @author Wei Kang Lim
 *
 */
public class DataFile {
	private static Calendar calendar;
	private final static String baseDir = "data";
	private BufferedWriter bw;
	private Path filePath;
	public enum FileIdentifier { 
		SIZE,
		FRACTION,
		SIZE_AND_FRACTION
	};
	
	public DataFile(String name, String [] configurations, FileIdentifier type, boolean multiple){
		if(DataFile.calendar == null){
			calendar = Calendar.getInstance();
		}
		String phiN = "";
		String runNumber = "";
		String q = "";

		for(int i=0; i < configurations.length; i++){
			String [] tokens = configurations[i].split(": ");
			String key = tokens[0];
			String value = tokens[1];
			
			if(key.equals("Nanoparticle Volume Fraction")){
				phiN = value;
				double phi_N = Double.parseDouble(phiN);
				DecimalFormat threeDecimal = new DecimalFormat("#0.###");
				DecimalFormat largeDecimal = new DecimalFormat("0.##E0");
				DecimalFormat fourDecimal = new DecimalFormat("#0.####");

				if(phi_N < 0.0001){
					phiN = largeDecimal.format(phi_N);
				} else if(phi_N < 0.001){
					phiN = fourDecimal.format(phi_N);
				} else {
					phiN = threeDecimal.format(phi_N);
				}

			}
			
			if(key.equals("Polymer Colloid Ratio")){
				q = value;
			}
			
			if(key.equals("Run Number")){
				runNumber = value;
			}
			
		}
		
		String prefix;
		switch(type){
			case SIZE: 
				prefix = "phiN=" + phiN; break;
			case FRACTION:
				prefix = "q=" + q; break;
			case SIZE_AND_FRACTION:
				prefix = "phiN=" + phiN + ",q=" + q; break;
			default:
				prefix = "phiN=" + phiN; break;
		}

		String creationTime = calendar.get(Calendar.HOUR_OF_DAY) + "." + 
							  calendar.get(Calendar.MINUTE);
		String creationDate = calendar.get(Calendar.DAY_OF_MONTH) + "-" + 
				  (calendar.get(Calendar.MONTH) + 1) + "-" + // calendar.MONTH starts counting at 0
				  calendar.get(Calendar.YEAR) + " " + creationTime;
		Path dir;
		if(multiple){
			dir = Paths.get(DataFile.baseDir + "/" +  prefix + "_" + creationDate + "/" + runNumber);
		} else {
			dir = Paths.get(DataFile.baseDir + "/");
		}
		
		// create directory to store datafiles, categorized by date if it doesn't exist
		try {
			Files.createDirectories(dir);
		} catch (IOException e) {
			e.printStackTrace();
		}

		
		if(multiple){
			filePath = Paths.get(dir + "/" + prefix + "_" + name + ".dat");
		}else{
			filePath = Paths.get(dir + "/" + prefix + "_" + name + "." + creationDate +".dat");
		}
		
		try {
			bw = Files.newBufferedWriter(filePath, 
					Charset.forName("US-ASCII"), filePath.toFile().exists()? StandardOpenOption.APPEND : StandardOpenOption.CREATE);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		commentLines(configurations);
	}
	
	/**
	 * Adds the data to the file buffer.
	 * @param data Data to be added
	 */
	public void record(String data){
		try {
			bw.write(data);
			bw.newLine();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (NullPointerException e){
			System.out.println("record method called without calling constructor!");
		}
	}
	
	/**
	 * Adds a line of comment to the file buffer. 
	 * @param comment A line of comment.
	 */
	public void comment(String comment){
		try{
			bw.write("# " + comment);
			bw.newLine();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (NullPointerException e){
			System.out.println("commentLines method called without calling constructor!");
		}
	}
	
	/**
	 * Adds multiple lines of comments to the file buffer.
	 * @param comments An array of comments, with each being written on a new line.
	 */
	public void commentLines(String [] comments){
		try{
			for(int i = 0; i < comments.length; i++){
				bw.write("# " + comments[i]);
				bw.newLine();
			}
		} catch (IOException e) {
			e.printStackTrace();
		} catch (NullPointerException e){
			System.out.println("commentLines method called without calling constructor!");
		}
	}
	
	/**
	 * Writes the data to the file.
	 */
	public void write(){
		try{
			bw.flush();
		} catch (IOException e){
			e.printStackTrace();
			System.out.println("Error writing to file.");
		}
	}
	
	/**
	 * Closes the stream for writing.
	 */
	public void close(){
		try{
			bw.close();
		} catch(IOException e){
			e.printStackTrace();
			System.out.println("Problem closing the buffer stream.");
		}
	}
}
