package org.opensourcephysics.sip.CPM;

import org.opensourcephysics.numerics.Matrix3DTransformation;
import org.opensourcephysics.numerics.VectorMath;

public class Polymer extends Particle{
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
	
	@Override
	public boolean overlap(Particle particle) {
		if(particle instanceof Polymer){
			return false;
		} else if(this.equals(particle)){
			return false;
		} else{
			Nano nano = (Nano) particle;
			double [] originAxis = {0,1,0}; // AD
			//double [] originAxis = {0,0,1};
			double [] start = {this.getX() - Particle.getLx(), this.getY() - Particle.getLy(), this.getZ() - Particle.getLz()};
			double [] boundary = {Particle.getLx(), Particle.getLy(), Particle.getLz()};
			double [] dist = new double[3];

			boolean normalAxis = true;
			for(int i = 0; i < 3; i++){
				normalAxis = normalAxis && (newAxis[i] == originAxis[i]);
			}
			if(!normalAxis ){
				Matrix3DTransformation transformation =  Matrix3DTransformation.createAlignmentTransformation(originAxis,newAxis);
				for(int x = 0; x <= 2; x++){
					for(int y = 0; y <= 2; y++){
						for(int z = 0; z <= 2; z++){
							double polymer[] = {start[0] + boundary[0]*x, start[1] + boundary[1]*y, start[2] + boundary[2]*z};
							double [] point = {nano.getX(), nano.getY(), nano.getZ()};
							if(!normalAxis ){
								transformation.setOrigin(polymer[0], polymer[1], polymer[2]);
								point = transformation.direct(point);
							}
							dist[0] = Math.abs(point[0]-(start[0] + boundary[0]*x) );
							dist[1] = Math.abs(point[1]-(start[1] + boundary[1]*y) );
							dist[2] = Math.abs(point[2]-(start[2] + boundary[2]*z) );
							if(Math.pow(dist[0]/(this.getrX()+nano.getrX()),2) + Math.pow(dist[1]/(this.getrY()+nano.getrY()),2) + Math.pow(dist[2]/(this.getrZ()+nano.getrZ()),2) < 1){ // AD
							//if(Math.pow(dist[0]/this.getrX(),2) + Math.pow(dist[1]/this.getrY(), 2) + Math.pow(dist[2]/this.getrZ(),2) < 1){
								return true;
							}
						}
					}
				}
			} else {
				double [] point = {nano.getX(), nano.getY(), nano.getZ()};
				for(int x = 0; x <= 2; x++){
					for(int y = 0; y <= 2; y++){
						for(int z = 0; z <= 2; z++){
							dist[0] = Math.abs(point[0]-(start[0] + boundary[0]*x) );
							dist[1] = Math.abs(point[1]-(start[1] + boundary[1]*y) );
							dist[2] = Math.abs(point[2]-(start[2] + boundary[2]*z) );
							if(Math.pow(dist[0]/(this.getrX()+nano.getrX()),2) + Math.pow(dist[1]/(this.getrY()+nano.getrY()), 2) + Math.pow(dist[2]/(this.getrZ()+nano.getrZ()),2) < 1){ // AD
							//if(Math.pow(dist[0]/this.getrX(),2) + Math.pow(dist[1]/this.getrY(), 2) + Math.pow(dist[2]/this.getrZ(),2) < 1){
								return true;
							}
						}
					}
				}
			}
			return false;	
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

	public String toString(){
		return "oldAxis: " + oldAxis[0] + ", " + oldAxis[1] + ", " + oldAxis[2] + "\n" +
			   "newAxis: " + newAxis[0] + ", " + newAxis[1] + ", " + newAxis[2];
	}
}
