package org.opensourcephysics.sip.CPM;

import org.opensourcephysics.numerics.Matrix3DTransformation;
import org.opensourcephysics.numerics.PBC;
import org.opensourcephysics.numerics.Root;
import org.opensourcephysics.numerics.VectorMath;

public class Polymer extends Particle{
	private static boolean exact;
	private static double tolerance;
	private static double shapeTolerance;
	private static double rotTolerance;
	private static double default_eX;
	private static double default_eY;
	private static double default_eZ;
	private static double q;
	private double eX;
	private double eY;
    private double eZ;
    private double oldAxis[] = {0,0,1};
    private double newAxis[] = {0,1,0}; // transformation axis
    private ESOverlapPolynomial overlapPolynomial;



	public Polymer(double x_, double y_, double z_, double tolerance_, double eX, double eY, double eZ, double q_){
		setX(x_);
		setY(y_);
		setZ(z_);
		seteX(eX);
		seteY(eY);
		seteZ(eZ);
		tolerance = tolerance_;
	}
	
	public Polymer(double x_, double y_, double z_){
		setX(x_);
		setY(y_);
		setZ(z_);
		seteX(default_eX);
		seteY(default_eY);
		seteZ(default_eZ);
	}
	
	/**
	 * Checks if a polymer overlaps a given particle.
	 * @param Particle A particle object.
	 * @return boolean true if an overlap occurs, false otherwise.
	 */
	public boolean overlap(Particle particle) {
		if(particle instanceof Polymer){
			return false;
		} else if(this.equals(particle)){
			return false;
		} else{
			Nano nano = (Nano) particle;
			
			// Check if polymer is spherical. 
			if(this.geteX() == this.geteY() && this.geteX() == this.geteZ()){
				return this.separation(nano) < Math.pow(this.getrX() + nano.getrX(), 2); // Calculate the pbc distance directly.
			} else { // Polymer is ellipsoidal				
				if(exact){
					/** Initial filter: **/
					// If a sphere is outside of a sphere formed by the longest radius of the ellispoid, reject it right away.
					double maxEllipsRadius = Math.max(Math.max(this.getrX(), this.getrY()), this.getrZ());
					if(this.separation(nano) >= Math.pow(maxEllipsRadius + nano.getrX(), 2)){ 
						return false;
					}
					
					double [] originAxis = {0,1,0};
					double [] point = {nano.getX(), nano.getY(), nano.getZ()};
					boolean normalAxis = true;
					for(int i = 0; i < 3; i++){
						normalAxis = normalAxis && (newAxis[i] == originAxis[i]);
					}	
					
					if(!normalAxis  && this.getRotTolerance() != 0){ // Rotations needed
						Matrix3DTransformation transformation =  Matrix3DTransformation.createAlignmentTransformation(originAxis,newAxis);
						transformation.setOrigin(this.getX(), this.getY(), this.getZ());
						point = transformation.direct(point);
					}

					/** Exact solution **/
					double [] ellips = {Math.pow(this.getrX(),2), Math.pow(this.getrY(),2), Math.pow(this.getrZ(),2)}; // ellipsoid radii
					// sphere distances from center of ellipsoid
					double [] sphere = {point[0] - this.getX(), point[1] - this.getY(), point[2] - this.getZ()}; 
					
					// Update the overlap polynomial with the new ellipsoid shape, and new distance of sphere.
					if(overlapPolynomial == null){ // If the function hasn't been instantiated, instantiate a new object
						overlapPolynomial = new ESOverlapPolynomial(ellips, sphere);
					} else { // Object instantiated, reuse it instead of creating new references for performance.
						overlapPolynomial.update(ellips, sphere);
					}
					
					double yRatio = ellips[1] / ellips[0]; //  ratio of the radii, lambda2/lambda1
					double zRatio = ellips[2] / ellips[0]; //  lambda3/lambda1
					// solve for root, which must be in between (0, maxEllipsRadius).
					// 0.001 is currently set as the acceptable computational tolerance.
					// x, y, z are the coordinates of the closest point
					
					double max = sphere[0] > 0 ? maxEllipsRadius : 0;
					double min = sphere[0] > 0 ? 0 : -maxEllipsRadius;
					double x = Root.bisection(overlapPolynomial, min , max, 0.0001); 
					double y = yRatio*x*sphere[1]/(sphere[0] + (yRatio-1)*x);
					double z = zRatio*x*sphere[2]/(sphere[0] + (zRatio-1)*x);
					double ellipsEquation = Math.pow(x/this.getrX(),2)+Math.pow(y/this.getrY(),2)+Math.pow(z/this.getrZ(),2);
					
					
					System.out.println("Ellipsoid radii: " + this.getrX() + " " + this.getrY() + " " + this.getrZ());
					System.out.println("Closest: " + x + " " + y + " " + z);
					System.out.println("Equation of ellipsoid: " + ellipsEquation);
					System.out.println("Nano center: " + sphere[0] + " " + sphere[1] + " " + sphere[2]);
					return Math.pow(x-sphere[0], 2) + Math.pow(y-sphere[1], 2) + Math.pow(z-sphere[2], 2) < Math.pow(nano.getrX(),2);
				} else{
					double [] originAxis = {0,1,0}; // The unrotated vector of the ellipsoid.
					
					/** The current way the overlap detection works:
					 *  All periodic images of the polymer are generated, and if a nanoparticle is found to be within the approximate exclusion shell
					 *  given by the sum of the two radii in each dimension, an overlap is detected.
					 */
					// Initialize the starting position of the polymer corresponding to the periodic image closest to (0,0,0)
					double [] startImage = {this.getX() > Particle.getLx()/2 ? this.getX() - Particle.getLx() : this.getX(), 
									   this.getY() > Particle.getLy()/2 ? this.getY() - Particle.getLy() : this.getY(), 
									   this.getZ() > Particle.getLz()/2 ? this.getZ() - Particle.getLz() : this.getZ()};
					double [] boundary = {Particle.getLx(), Particle.getLy(), Particle.getLz()};
					double [] dist = new double[3];
		
					boolean normalAxis = true;
					for(int i = 0; i < 3; i++){
						normalAxis = normalAxis && (newAxis[i] == originAxis[i]);
					}
					
					if(!normalAxis ){ // Rotations needed
						Matrix3DTransformation transformation =  Matrix3DTransformation.createAlignmentTransformation(originAxis,newAxis);
						for(int x = 0; x <= 1; x++){
							for(int y = 0; y <= 1; y++){
								for(int z = 0; z <= 1; z++){
									double polymer[] = {startImage[0] + boundary[0]*x, startImage[1] + boundary[1]*y, startImage[2] + boundary[2]*z};
									double [] point = {nano.getX(), nano.getY(), nano.getZ()};
									if(!normalAxis ){
										transformation.setOrigin(polymer[0], polymer[1], polymer[2]);
										point = transformation.direct(point);
									}
									dist[0] = Math.abs(point[0]-(startImage[0] + boundary[0]*x) );
									dist[1] = Math.abs(point[1]-(startImage[1] + boundary[1]*y) );
									dist[2] = Math.abs(point[2]-(startImage[2] + boundary[2]*z) );
									if(Math.pow(dist[0]/(this.getrX()+nano.getrX()),2) + Math.pow(dist[1]/(this.getrY()+nano.getrY()),2) + Math.pow(dist[2]/(this.getrZ()+nano.getrZ()),2) < 1){ // AD
									//if(Math.pow(dist[0]/this.getrX(),2) + Math.pow(dist[1]/this.getrY(), 2) + Math.pow(dist[2]/this.getrZ(),2) < 1){
										return true;
									}
								}
							}
						}
					} else {
						dist[0] = PBC.separation(nano.getX()-this.getX(), boundary[0]);
						dist[1] = PBC.separation(nano.getY()-this.getY(), boundary[1]);
						dist[2] = PBC.separation(nano.getZ()-this.getZ(), boundary[2]);
						if(Math.pow(dist[0]/(this.getrX()+nano.getrX()),2) + Math.pow(dist[1]/(this.getrY()+nano.getrY()), 2) + Math.pow(dist[2]/(this.getrZ()+nano.getrZ()),2) < 1){ // AD
								return true;
						}
					}
					return false;
				}
			}
		}
	}
	
	/** Performs a shape change
	 * 
	 */
	public void shapeChange(){
		double newEX = this.geteX() + 10 * shapeTolerance * 2. * (Math.random() - 0.5);
		double newEY = this.geteY() + 3 * shapeTolerance * 2. * (Math.random() - 0.5);
		double newEZ = this.geteZ() + shapeTolerance * 2. * (Math.random() - 0.5);
		
		// Reject changes for negative radii eigenvalue immediately
		if(newEY < 0 || newEY < 0 || newEZ < 0){
			return;
		}
		
		this.seteX(newEX);
		this.seteY(newEY);
		this.seteZ(newEZ);
	}
	
	/**
	 * Performs a rotation on the polymer.
	 */
	public void rotate(){
		this.setOldAxis(this.getNewAxis());
		double [] a = this.getOldAxis();
		double [] v = new double[3];
		
		// generate randomly oriented vector v 
		for(int i = 0 ; i < v.length; i++){
			v[i] = Math.random() - 0.5;
		}

		// normalize new (randomly oriented) vector v 
		VectorMath.normalize(v);
				
		// addition of the old and new vector 
		// Note: rotTolerance, which replaces rotMagnitude, should be << 1 (e.g. 0.1)
		for(int i = 0; i < v.length; i++){
			a[i] = a[i] + rotTolerance*v[i];
		}

		// normalize result
		VectorMath.normalize(a);
		this.setNewAxis(a);

	}
	
	/**
	 * Moves the polymer with the tolerance specified by Polymer.tolerance.
	 */
	public void move(){
		super.move(getTolerance());
	}
	
	public static void setTolerance(double tolerance) {
		Polymer.tolerance = tolerance;
	}
	
	public static double getDefault_eX() {
		return default_eX;
	}

	public static void setDefault_eX(double defaulteX_) {
		Polymer.default_eX = defaulteX_;
	}
	public static void setDefault_eY(double defaulteY_) {
		Polymer.default_eY = defaulteY_;
	}
	public static void setDefault_eZ(double defaulteZ_) {
		Polymer.default_eZ = defaulteZ_;
	}
	public void seteX(double eX_){
		eX = eX_;
		setrX(toRadius(eX));
	}
	
	public void seteY(double eY_){
		eY = eY_;
		setrY(toRadius(eY));
	}
	
	public void seteZ(double eZ_){
		eZ = eZ_;
		setrZ(toRadius(eZ));
	}
	public void setOldAxis(double [] axis){
		if(axis.length == 3){
			this.oldAxis = axis;
		}
	}

	public static void setQ(double q_){
		Polymer.q = q_;
	}
	public void setNewAxis(double newAxis[]) {
		this.newAxis = newAxis;
	}


  // -------------------------------------
  // Getter methods
  // -------------------------------------
	public Matrix3DTransformation getTransformation(){
		Matrix3DTransformation transform;
		if(!isRotated()){
			double [][] identity = {
					{1, 0, 0},
					{0, 1, 0},
					{0, 0, 1}
			};
			 transform = new Matrix3DTransformation(identity);
		} else{
			double [] originAxis = {0,0,1};
			transform = Matrix3DTransformation.createAlignmentTransformation(originAxis, newAxis);
		}
		return transform;
	}

	
  public double[] getNewAxis() {
	  	double [] returnAxis = new double[newAxis.length];
	  	for(int i =0; i < newAxis.length; i++){
	  		returnAxis[i] = newAxis[i];
	  	}
		return returnAxis;
	}


	public double geteX(){
		return eX;
	}
	
	public double geteY(){
		return eY;
	}
	
	public double geteZ(){
		return eZ;
	}
	public double[] getOldAxis(){
		double a[] = new double[oldAxis.length];		
		for(int i = 0; i < oldAxis.length; i ++){
			a[i] = oldAxis[i];
		}
		return a;
	}
	public static double getQ(){
		return Polymer.q;
	}
	public static double getDefault_eZ() {
		return default_eZ; // AD: eY --> eZ
	}
	public static double getTolerance() {
		return tolerance;
	}
	public static double getDefault_eY() {
		return default_eY;
	}
	
	/**
	 * Converts radius eigenvalue to radius.
	 * @param ei Eigenvalue
	 * @param v enumeration representing the radius axis
	 * @return The radius
	 */
	public double toRadius(double ei) {
		return (q/2.) * Math.sqrt(18. * ei);
	}
	
	public boolean isRotated(){
		for(int i = 0; i < oldAxis.length; i++){
			if(oldAxis[i] != newAxis[i]) return true;
		}
		return false;
	}
	
	public static double getRotTolerance() {
		return rotTolerance;
	}

	public static void setRotTolerance(double rotTolerance) {
		Polymer.rotTolerance = rotTolerance;
	}
	
	public static double getShapeTolerance() {
		return shapeTolerance;
	}

	public static void setShapeTolerance(double shapeTolerance) {
		Polymer.shapeTolerance = shapeTolerance;
	}
	
	public static boolean getExact(){
		return exact;
	}
	
	public static void setExact(boolean exact){
		Polymer.exact = exact;
	}

	public String toString(){
		return "oldAxis: " + oldAxis[0] + ", " + oldAxis[1] + ", " + oldAxis[2] + "\n" +
			   "newAxis: " + newAxis[0] + ", " + newAxis[1] + ", " + newAxis[2];
	}
	
	public String toPovray(){
		double [] oldAxis = this.getOldAxis();
		double [] newAxis = this.getNewAxis();
		String pString = String.format("sphere{\n"
				+ "\t<%.3f, %.3f, %.3f>, 1\n"
				+ "\tscale <%.3f, %.3f, %.3f>\n"
				+ "\ttexture {\n"
				+ "\t		pigment { color Red }\n"
				+ "\t}\n"
				+ "Reorient_Trans(<%.3f, %.3f, %.3f>, <%.3f, %.3f, %.3f>)"
				+ "}", 
				this.getX(), this.getY(), this.getZ(),
				this.getrX(), this.getrY(), this.getrZ(),
				oldAxis[0], oldAxis[1], oldAxis[2],
				newAxis[0], newAxis[1], newAxis[2]);
		return pString;
	}
}
