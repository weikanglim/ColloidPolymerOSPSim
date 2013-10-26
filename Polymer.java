package org.opensourcephysics.sip.CPM;

public class Polymer extends Particle{
	private static double tolerance;
	private static double default_eX;
	private static double default_eY;
	private static double default_eZ;
	private static double q;
	private double eX;
	private double eY;
    private double eZ;


	public Polymer(double x_, double y_, double z_, double tolerance_, double eX, double eY, double eZ, double q_){
		setX(x_);
		setY(y_);
		setZ(z_);
		seteX(eX, q_);
		seteY(eY, q_);
		seteZ(eZ, q_);
		tolerance = tolerance_;
	}
	
	public Polymer(double x_, double y_, double z_){
		setX(x_);
		setY(y_);
		setZ(z_);
		seteX(default_eX, q);
		seteY(default_eY, q);
		seteZ(default_eZ, q);
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
	
	public static double getTolerance() {
		return tolerance;
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
	
	public static double getDefault_eY() {
		return default_eY;
	}

	public static void setDefault_eY(double defaulteY_) {
		Polymer.default_eY = defaulteY_;
	}
	
	public static double getDefault_eZ() {
		return default_eY;
	}

	public static void setDefault_eZ(double defaulteZ_) {
		Polymer.default_eZ = defaulteZ_;
	}
	
	public double getLX(){
		return Math.pow(this.getrX(), 2);
	}
	
	public double getLY(){
		return Math.pow(this.getrY(), 2);
	}
	
	
	public void seteX(double eX_, double q_){
		eX = eX_;
		setrX(toRadius(eX_, q_));
	}
	
	public void seteY(double eY_, double q_){
		eY = eY_;
		setrY(toRadius(eY_, q_));
	}
	
	public void seteZ(double eZ_, double q_){
		eZ = eZ_;
		setrZ(toRadius(eZ_,q_));
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
	
	/**
	 * Converts radius eigenvalue to radius.
	 * @param ei Eigenvalue
	 * @param v enumeration representing the radius axis
	 * @return The radius
	 */
	public double toRadius(double ei, double q_) {
		return q_ / 2 * Math.sqrt(18 * ei);
	}
	
	public static void setQ(double q_){
		Polymer.q = q_;
	}
	
	public static double getQ(){
		return Polymer.q;
	}
}
