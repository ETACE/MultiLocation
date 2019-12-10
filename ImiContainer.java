import java.util.ArrayList;

/*
 * This class is used to store the imitation scheme across the different scenario based batches of the MC in the location decision.
 * */

public class ImiContainer {

	//Imitation schedule
	ArrayList<ImitationScheme> listOfEvents;
	
	/*Constructor*/
	ImiContainer(){
	
		listOfEvents = new ArrayList<ImitationScheme>();
		
	}

}
