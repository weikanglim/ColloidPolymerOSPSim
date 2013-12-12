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
    private double transformAxis[] = {0,0,1}; // transformation axis


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
			double [] point = {nano.getX(), nano.getY(), nano.getZ()};
			double [] originAxis = {0,0,1};
			boolean normalAxis = true;
			for(int i = 0; i < 3; i++){
				normalAxis = normalAxis && (transformAxis[i] == originAxis[i]);
			}
			
			if(!normalAxis ){
				Matrix3DTransformation transformation =  Matrix3DTransformation.createAlignmentTransformation(originAxis,transformAxis);
				transformation.setOrigin(this.getX(), this.getY(), this.getZ());
				point = transformation.direct(point);
			}
			double [] start = {this.getX() - Particle.getLx(), this.getY() - Particle.getLy(), this.getZ() - Particle.getLz()};
			double [] boundary = {Particle.getLx(), Particle.getLy(), Particle.getLz()};
			double [] dist = new double[3];
			
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
		Matrix3DTransformation transform;
		if(!isRotated()){
			double [][] identity = {
					{1, 0, 0},
					{0, 1, 0},
					{0, 0, 1}
			};
			 transform = new Matrix3DTransformation(identity);
		} else{
//		double [] origin = {this.getX(), this.getY(), this.getZ()};
			transform = Matrix3DTransformation.createAlignmentTransformation(axis, transformAxis);
//			double [] matrix = new double[16];
//			transform.getFlatMatrix(matrix);
//			System.out.println(String.format("(%e, %e, %e)", matrix[0],matrix[1],matrix[2]));
//			System.out.println(String.format("(%e, %e, %e)", matrix[4],matrix[5],matrix[6]));
//			System.out.println(String.format("(%e, %e, %e)", matrix[8],matrix[9],matrix[10]));
//		transform.setOrigin(origin);
		}
		return transform;
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
	
	public boolean isRotated(){
		for(int i = 0; i < axis.length; i++){
			if(axis[i] != transformAxis[i]) return true;
		}
		return false;
	}
	
}
