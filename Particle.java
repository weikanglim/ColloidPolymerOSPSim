package org.opensourcephysics.sip.CPM;
import java.util.ArrayList;

import org.opensourcephysics.numerics.PBC;
import org.opensourcephysics.sip.CPM.CPM.Vector;


public abstract class Particle {
	private double x;
	private double y;
	private double z;
	private double rX;
	private double rY;
	private double rZ;
	private static double Lx;
	private static double Ly;
	private static double Lz;
	public ArrayList<Particle> intersectPairs = new ArrayList<Particle>();

	
	public double separation(Particle particle){
		double dx = this.x-particle.x;
		double dy = this.y-particle.y;
		double dz = this.z-particle.z;
		return dx*dx+dy*dy+dz*dz;
	}
	
	public abstract boolean overlap(Particle particle);
	
	public void move(double tolerance){
		setX(PBC.position(getX() + tolerance*2.*(Math.random()-0.5), Lx));
		setY(PBC.position(getY() + tolerance*2.*(Math.random()-0.5), Ly));
		setZ(PBC.position(getZ() + tolerance*2.*(Math.random()-0.5), Lz));
	}
	
	public static void setBoundaries(double Lx_, double Ly_, double Lz_){
		Lx = Lx_;
		Ly = Ly_;
		Lz = Lz_;
	}

	
	public void setrX(double rX_){
		if(rX_ > 0) rX = rX_;
	}
	
	public void setrY(double rY_){
		if(rY_ > 0) rY = rY_;
	}
	
	public void setrZ(double rZ_){
		if(rZ_ > 0 ) rZ = rZ_;
	}
	
	public double getrX(){
		return rX;
	}
	
	public double getrY(){
		return rY;
	}
	
	public double getrZ(){
		return rZ;
	}
	
	public void setX(double x_){
		x = x_;
	}
	
	public void setY(double y_){
		y = y_;
	}
	
	public void setZ(double z_){
		z = z_;
	}
	
	public double getX(){
		return x;
	}
	
	public double getY(){
		return y;
	}
	
	public double getZ(){
		return z;
	}
}
