package org.opensourcephysics.sip.CPM;

import org.opensourcephysics.numerics.Matrix3DTransformation;

public class Polymer extends Particle{
	private static double tolerance;
	private static double default_eX;
	private static double default_eY;
	private static double default_eZ;
	private static double q;
	private double eX;
	private double eY;
    private double eZ;
    private double axis[] = {0,0,1};
    private double transformAxis[];


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
			double dx = Math.abs(nano.getX()-this.getX());
			double dy = Math.abs(nano.getY()-this.getY());
			double dz = Math.abs(nano.getZ()-this.getZ());
			return Math.pow(dx/this.getrX(),2) + Math.pow(dy/this.getrY(), 2) + Math.pow(dz/this.getrZ(),2) < 1;
		}	
	}

	
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
	public void setAxis(double [] axis){
		if(axis.length == 3){
			this.axis = axis;
		}
	}

	public static void setQ(double q_){
		Polymer.q = q_;
	}

	public Matrix3DTransformation getTransformation(){
		return Matrix3DTransformation.createAlignmentTransformation(axis, transformAxis);
	}

	
  public double[] getTransformAxis() {
		return transformAxis;
	}

	public void setTransformAxis(double transformAxis[]) {
		this.transformAxis = transformAxis;
	}

	// -------------------------------------
  // Getter methods
  // -------------------------------------
	/**
	 * 
	 * @return
	 */
	public double getLX(){
		return Math.pow(this.getrX(), 2);
	}
	
	public double getLY(){
		return Math.pow(this.getrY(), 2);
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
	
	public double[] getAxis(){
		double a[] = new double[axis.length];		
		for(int i = 0; i < axis.length; i ++){
			a[i] = axis[i];
		}
		return a;
	}
	public static double getQ(){
		return Polymer.q;
	}
	public static double getDefault_eZ() {
		return default_eY;
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
		return q / 2 * Math.sqrt(18 * ei);
	}
}
