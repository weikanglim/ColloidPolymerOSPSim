package org.opensourcephysics.sip.CPM;


public class Polymer extends Particle{
	private static double tolerance;
	private static double default_U;
	private static double q;
	public double u;
    private double oldAxis[] = {0,0,1};
    private double newAxis[] = {0,1,0}; // transformation axis


	public Polymer(double x_, double y_, double z_, double tolerance_, double u, double q_){
		setX(x_);
		setY(y_);
		setZ(z_);
		setU(u);
		tolerance = tolerance_;
	}
	
	public Polymer(double x_, double y_, double z_){
		setX(x_);
		setY(y_);
		setZ(z_);
		setU(default_U);

	}
	
	@Override
	public boolean overlap(Particle particle) {
		if(particle instanceof Polymer){
			return false;
		} else if(this.equals(particle)){
			return false;
		} else{
			Nano nano = (Nano) particle;
			double [] start = {this.getX() - Particle.getLx(), this.getY() - Particle.getLy(), this.getZ() - Particle.getLz()};
			double [] boundary = {Particle.getLx(), Particle.getLy(), Particle.getLz()};
			double [] dist = new double[3];

			double [] point = {nano.getX(), nano.getY(), nano.getZ()};
			for(int x = 0; x <= 2; x++){
				for(int y = 0; y <= 2; y++){
					for(int z = 0; z <= 2; z++){
						dist[0] = Math.abs(point[0]-(start[0] + boundary[0]*x) );
						dist[1] = Math.abs(point[1]-(start[1] + boundary[1]*y) );
						dist[2] = Math.abs(point[2]-(start[2] + boundary[2]*z) );
						if(Math.pow(dist[0]/this.getrX(),2) + Math.pow(dist[1]/this.getrY(), 2) + Math.pow(dist[2]/this.getrZ(),2) < 1){
							return true;
						}
					}
				}
			}
			return false;	
		}
	}

	
	public void move(){
		super.move(getTolerance());
	}
	
	public static void setTolerance(double tolerance) {
		Polymer.tolerance = tolerance;
	}
	
	public static void setDefault_U(double defaultU_) {
		Polymer.default_U = defaultU_;
	}
	public static void setQ(double q_){
		Polymer.q = q_;
	}
	public void setU(double u){
		this.u = u;
		setrX(toRadius(u));
		setrY(toRadius(u));
		setrZ(toRadius(u));
	}

	// -------------------------------------
  // Getter methods
  // -------------------------------------	
	public double getU(){
		return u;
	}
	public static double getQ(){
		return Polymer.q;
	}
	public static double getDefault_U() {
		return default_U;
	}
	public static double getTolerance() {
		return tolerance;
	}	
	/**
	 * Converts radius eigenvalue to radius.
	 * @param ei Eigenvalue
	 * @param v enumeration representing the radius axis
	 * @return The radius
	 */
	public double toRadius(double u) {
		return q / 2 * Math.sqrt(6 * u);
	}
		
	public String toString(){
		return "oldAxis: " + oldAxis[0] + ", " + oldAxis[1] + ", " + oldAxis[2] + "\n" +
			   "newAxis: " + newAxis[0] + ", " + newAxis[1] + ", " + newAxis[2];
	}
}
