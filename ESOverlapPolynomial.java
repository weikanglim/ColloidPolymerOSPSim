package org.opensourcephysics.sip.CPM;

/**
 * This class represents the Ellipsoid-Sphere overlap characteristic polynomial, 
 * which root gives one of the coordinates of the point on the ellipsoid, that is closest to the sphere.
 * The other coordinates can be found with this given coordinate, using the constraint of the point being on the ellipsoid (1D freedom).
 * <br>
 * @author Wei Kang Lim
 *
 */
public class ESOverlapPolynomial implements org.opensourcephysics.numerics.Function {
	private double l1,l2,l3; // eigenvalues
	private double x0,y0,z0; // coordinate of sphere
	
	@SuppressWarnings("unused")
	private ESOverlapPolynomial(){} // hide constructor
	
	
	/**
	 * Instantiates the function with the required variables.
	 * @param ellips double[] An array containing the ellipsoid radii squared.<br>
	 * 				 [0,1,2] =  The x,y,z ellipsoid radii.
	 * @param sphere double[] An array containing the sphere coordinates. <br>
	 * 				 [0,1,2] = The sphere coordinates, relative to the center of ellipsoid. <br>
	 */
	public ESOverlapPolynomial(double [] ellips, double [] sphere){
		update(ellips, sphere);
	}
	
	/**
	 * Updates the polynomial with new parameters.
	 * @param ellips double[] An array containing the ellipsoid radii squared.<br>
	 * 				 [0,1,2] =  The x,y,z ellipsoid radii.
	 * @param sphere double[] An array containing the sphere coordinates. <br>
	 * 				 [0,1,2] = The sphere coordinates, relative to the center of ellipsoid. <br>
	 */
	public void update(double [] ellips, double [] sphere){
		l1 = ellips[0]; l2=ellips[1]; l3=ellips[2];
		x0 = sphere[0]; y0 = sphere[1]; z0=sphere[2];
	}
	
	
	/**
	 * Evaluates the polynomial at a given point.
	 */
	public double evaluate(double x){
		double a = l2/(l1*l1) * y0*y0;
		double b = l2/l1 - 1;
		double c = l3/(l1*l1) * z0*z0;
		double d = l3/l1 - 1;
		double e = this.x0*this.x0/l1;
		double u = this.x0 / x;
		
		 return e/(u*u) + a/((u+b)*(u+b)) + c/((u+d)*(u+d)) - 1;
//		return e*(u+b)*(u+b)*(u+d)*(u+d) + a*u*u*(u+d)*(u+d) + c*u*u*(u+b)*(u+b) - u*u*(u+b)*(u+b)*(u+d)*(u+d);
	}
	
	public String toString(){
		double a = l2/(l1*l1) * y0*y0;
		double b = l2/l1 - 1;
		double c = l3/(l1*l1) * z0*z0;
		double d = l3/l1 - 1;
		
		String result = "x^2/$l1+$a/($x0/x+$b)^2+$c/($x0/x+$d)^2";
		result = result.replace("$l1", Double.toString(l1));
		result = result.replace("$x0", Double.toString(x0));
		result = result.replace("$a", Double.toString(a));
		result = result.replace("$b", Double.toString(b));
		result = result.replace("$c", Double.toString(c));
		result = result.replace("$d", Double.toString(d));
		
		return result;
	}
}