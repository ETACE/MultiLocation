import uchicago.src.sim.gui.NetworkDrawable;
import uchicago.src.sim.network.DefaultDrawableNode;



public class FirmNode extends DefaultDrawableNode  {
	
	
	Firm aFirm;
	
	FirmNode(Firm firm ,NetworkDrawable drawable) {
	    super(drawable);
	
		aFirm = firm;
		
	
	}
	
	
	FirmNode(Firm firm  ){
			
		
			aFirm = firm;
			

			
		}
	
	

}
