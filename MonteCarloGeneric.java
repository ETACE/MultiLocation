import java.util.ArrayList;

/**
*This class is a generic interface for different implementations of the Montecarlo simulator for computing the NPV. 
*(Note: eventually we have only used one, i.e. MonteCarloSimulatorSerial )
*/


public class MonteCarloGeneric {
	
	MonteCarloSimulatorSerial monteCarloSerial ;
	
	
	MonteCarloGeneric(Firm bFirm, ArrayList<FirmReference> compList){
	
		monteCarloSerial = new MonteCarloSimulatorSerial( bFirm,  compList);
	
	}
	
	/**Setup function
	 * */
	void setup(){
		
		monteCarloSerial.setup();
	}
	
	
	/**Run MC simulations for baseline case
	 * */
	double monteCarloBase(){
		
	
		return monteCarloSerial.monteCarloBase();
	
	}
	
	/**Run MC simulations for exit case
	 * */
	double monteCarloExit(int locationID){
		
		return monteCarloSerial.monteCarloExit(locationID);
	
	}
	
	
	/**Run MC simulations for entry case
	 * */
	double monteCarloEntry(int locationID){
	
		return monteCarloSerial.monteCarloEntry( locationID);
	
	}
	
	/**Run MC simulations for switching case
	 * */
	double monteCarloEntryExit(int idLocationExit, int idLocationEntry){
		
		return monteCarloSerial.monteCarloEntryExit( idLocationExit,  idLocationEntry);
	
	}
	
	/**clean up
	 * */
	void destructMontecarloSimulator(){
		
		 monteCarloSerial.destructMontecarloSimulator( );
	
	}
	
	
	
	void setEntryID(int id){
		
		 monteCarloSerial.entryID = id;
		
	}
	
	
	void setPendingEntry(boolean entry){
		
			 monteCarloSerial.pendingEntry = entry;
		
	}
	

}
