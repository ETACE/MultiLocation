import java.awt.Color;

import uchicago.src.sim.gui.DrawableEdge;
import uchicago.src.sim.gui.SimGraphics;
import uchicago.src.sim.network.DefaultDrawableEdge;
import uchicago.src.sim.network.DefaultEdge;
import uchicago.src.sim.network.Node;


public class FirmEdge extends DefaultEdge implements DrawableEdge{
	
	
	Color color;
	
	FirmEdge(){
		
		super();
		
	}
	
	FirmEdge(Node node1, Node node2){
		
		super(node1, node2);
	
		
	}
	
	
	public void draw(SimGraphics g, int fromX, int toX, int fromY, int toY) {

	  //    g.drawLink(color.RED, fromX, toX, fromY, toY, (float) super.strength);
	    
	  }

	
	
	
	

}
