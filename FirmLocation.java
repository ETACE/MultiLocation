import java.awt.Color;
import java.util.ArrayList;

import uchicago.src.sim.gui.Drawable;
import uchicago.src.sim.gui.SimGraphics;

/**
 * This class holds the firm specific data related to that location
 */
public class FirmLocation implements Cloneable, Drawable {
    

	
	
	/****************Cloning   to avoid problems in deep copying************************/
	
	public FirmLocation clone()  {
        try {
			return (FirmLocation) super.clone();
		} catch (CloneNotSupportedException e) {
			
			e.printStackTrace();
			return null;
		}
		
        
    }
		
		int firmID;
		int locationID;
		double academicActivity;
		int numInnovators;
		int numImitators;
		
		int type;
		boolean toBeConsidered;
		
		protected ArrayList<FirmReference> competitorList;
		
		ArrayList<Integer> orderedIDList = new ArrayList<Integer>();

		
		/*Constructor No. 1*/
	FirmLocation(int id, int locID,  double a){
		
		firmID = id;
		locationID = locID;
		competitorList = new ArrayList<FirmReference>();
		academicActivity = a;
		
	}

	/*Constructor No 2*/
	FirmLocation(){
		
	}
	
	void countFirms(){
		
		numInnovators=0;
		numImitators =0;
		
		for(int i=0; i < competitorList.size();i++){
			
			if(competitorList.get(i).innovator)
				numInnovators++;
			else
				numImitators++;
			
		}
		
	}
	
	/*Methods for drawing*/

	@Override
	public void draw(SimGraphics G) {


		 G.drawFastRoundRect(Color.blue);
		
	}

	@Override
	public int getX() {
		
		return 2;
	}

	@Override
	public int getY() {
		return 5;
	}

}

