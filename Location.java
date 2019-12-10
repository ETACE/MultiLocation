import java.awt.Color;
import java.util.ArrayList;
import uchicago.src.sim.gui.Drawable;
import uchicago.src.sim.gui.SimGraphics;

/**
 * This class represents the locations where firms can decide to locate
 */

public class Location implements Cloneable, Drawable {

	
	/****************Cloning to avoid problems in deep copying************************/
	
	public Location clone()  {
        try {
			return (Location) super.clone();
		} catch (CloneNotSupportedException e) {
			
			e.printStackTrace();
			return null;
		}
		
    }
	
	/*Declaration of internal memory*/
	
	int locationID;
	double academicActivity;
	int numberFirmsInLocation;
	int coordinateX, coordinateY; 
	
	/*Container for local firms objects*/
	protected ArrayList<Firm> firmList = new ArrayList<Firm>();
	
	// List of ordered id
	ArrayList<Integer> orderedIDList = new ArrayList<Integer>();
	
	
	/*Constructor*/
	Location(int id, double a){
		
		locationID = id;
		academicActivity = a;
		coordinateX = Model.randomInt(0,Model.numLocations*Model.numFirms);
		coordinateY = locationID* Model.numFirms - (int)(Model.numFirms/2);
		numberFirmsInLocation =0;
		
	}

	
	/*These functions are used for drawing in the GUI mode*/
	int setXCoordinate(){
		
		
		return 0;
		
	}
	

	int setYCoordinate(){
	
		return 0;

	}

	@Override
	public void draw(SimGraphics g) {

		double value = 1.0;
		
		int red, green, blue;
		
		if(value<0.5){

			red  = 255;
			green=Math.min(255,(int) ((( value)/(0.5))*255));
			blue =0;
			
		}else{
			
			red  = 255;
			green =255;
			blue=(int) ((( value-0.5)/(1.0-0.5))*255);
			
			
		}
	
		
		Color color = new Color(red, green, blue);

	        g.drawHollowFastOval(color);	
		
	}

	@Override
	public int getX() {
		
		return ((int)(ModelGUI.dimX*0.5));
	}

	@Override
	public int getY() {
		
		return ((int)(ModelGUI.dimX*0.5)) ;
	}

}
