import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.concurrent.TimeUnit;
import uchicago.src.sim.analysis.DataRecorder;
import uchicago.src.sim.engine.IController;
import uchicago.src.sim.engine.SimInit;

/**
 * This class is used to run the model in Batch mode, where parameters that deviate from the defaults defined in Model.java are read from
 * command line as additional arguments. The simulation output is stored in SQLite databases.  Note that this is a sub class of Model.java
 *
 */


public class BachParallelModel extends SetGetModel {
	
	
	/*This defines the number of iterations per simulation run*/
	int numOfTimeSteps, numRuns ; 
	int firmToStore, entryExitID;
	int currentRun = 1;
	long startTime;

	
	 String currentPathXML;
	 protected IController control;
	
	/*This holds the data recorders for the firms*/
	 protected ArrayList<DataRecorder> firmRecorderList = new ArrayList<DataRecorder>();
	
	 protected ArrayList<String> parameterVaried = new ArrayList<String>();
	 protected ArrayList<String> parameterValuesVaried = new ArrayList<String>();
	
	/*Repast function for data recording*/
	DataRecorder recorderTimeSeries, recorderLastIteration, recorderEntry, recorderExit;
	
	// SQLite data base
	SQLiteJDBC dataBase;
	
	
	/***
	 * This is the setup method inherited from Model.java
	 */
	public void setup(){
		
		printDebug = false; /*Write data to text file*/
		printDebugLocation = false; /*Write data to text file*/
		
		startTime = System.currentTimeMillis();
		
		//original model has to be initialized first
		super.setup();
		
		System.out.println("buildmodel");

		setStoppingTime(numOfTimeSteps);

	}
	
	
	/***
	 * This is the buildModel method inherited from Model.java
	 */
	public void buildModel(){
		
		/*Execute buildModel of Model.java first*/
		super.buildModel();

		
		//Init data base
		dataBase = new SQLiteJDBC(currentPathXML);	
		dataBase.createFirmTable();
		dataBase.createLocationTable();
		dataBase.createEngryTable();
		dataBase.createExitTable();
		dataBase.createAggregateTable();
	
	}
	
	
	/***
	 * This is the method defining the rquired steps at the end of the simulation
	 */
	public void atEnd(){
		
		
		/*Close the data base connection*/
		dataBase.atEnd();
		long endTime   = System.currentTimeMillis();
		long totalTime = endTime - startTime;
	
		
		System.out.println("Simulation time: "+String.format("%d min, %d sec", 
			    TimeUnit.MILLISECONDS.toMinutes(totalTime),
			    TimeUnit.MILLISECONDS.toSeconds(totalTime) - 
			    TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(totalTime))
			));
		
	}
	
	
	
	/***
	 * This is the step() method inherited from Model.java
	 */
		public void step(){
	
		super.step();
		
	
		System.out.println(getTickCount());
	
		//Write data to data base:
		
		dataBase.insertAggregatedData((int) getTickCount(),this);
		
		for (int i=0; i< firmList.size(); i++){
	
			dataBase.insertFirm((int) getTickCount(),firmList.get(i));
			
		}
		
		
		for (int i=0; i< locationList.size(); i++){
	
			dataBase.insertLocation((int) getTickCount(),locationList.get(i));
			
		}
		

		for(int i=0; i < firmList.size();i++){
			
			for(int j=0; j < firmList.get(i).entryList.size(); j++){
		
				dataBase.insertEntry((int) getTickCount(), firmList.get(i).firmID, firmList.get(i).entryList.get(j));
				
			}
			
		}
		
		
		for(int i=0; i < firmList.size();i++){
			
			
			for(int j=0; j < firmList.get(i).exitList.size(); j++){
				
				
				dataBase.insertExit((int) getTickCount(), firmList.get(i).firmID,firmList.get(i).exitList.get(j));
				
			}
		
		}
		
		
		dataBase.commit();

	}

	
	
	
	/*Methods to write instances' variable values to the data recorder */
	
	
	/*Aggregated variables*/
	public double computeTotalOutput(){
		
		
		return totalOutput;
		
	}
	
	
	public double ComputePriceIndex(){
		
		return (double) (priceIndex);
		
	}
					
			
	public double computeQualityIndex(){
		
		return (double) (qualityIndex);
		
	}
			
	
	
	public double computeAverageFirmLocations(){
		
		return (double) (averageFirmLocations);
		
	}
	
	
	public int computeNumActiveFirms(){
		
		return (int) (numActiveFirms);
		
	}
	
		
	public double computeSdOutput(){
		
		return (double) (sdOutput);
		
	}
		
		
	
	public double computeSdQuality(){
		
		return (double) (sdQuality);
		
	}
		

	
	public double compueSdAverageLocation(){
		
		return (double) (sdAverageLocation);
		
	}
	

	public double computeAverageInnoProbability(){
		
		return (double) (averageInnoProbability);
		
	}
	
	
	public double computeAverageImiProbability(){
		
		return (double) (averageImiProbability);
		
	}


	/*Firm specific variables*/
	public double computeFirmID(){
		
		return firmList.get(firmToStore).firmID;
		
	}
	public double computeFirmOutput(){
		
		return firmList.get(firmToStore).equilQuantity;
		
	}
	public double computeFirmPrice(){
		
		return firmList.get(firmToStore).equilPrice;
		
	}
	public double computeFirmQuality(){
		
		return firmList.get(firmToStore).quality;
		
	}
	
	public double computeNumLocations(){
		
		return firmList.get(firmToStore).numLocationsActive;
		
	}

	public double computeFirmProfit(){
		
		return firmList.get(firmToStore).equilProfit;
		
	}
	
	
	
	public int computeEntryLocationID(){
		
		return firmList.get(firmToStore).entryList.get(entryExitID).locationID;
		
	}
	
	public int computeEntryNumCompetitors(){
		
		return firmList.get(firmToStore).entryList.get(entryExitID).numCompetitors;
		
	}
	
	


	public double computeEntryAverageQualityConcept(){
		
		return firmList.get(firmToStore).entryList.get(entryExitID).averageQualityConcept;
		
	}
	
	
	
	public int computeExitLocationID(){
		
		return firmList.get(firmToStore).exitList.get(entryExitID).locationID;
		
	}
	
	public int computeExitNumCompetitors(){
		
		return firmList.get(firmToStore).exitList.get(entryExitID).numCompetitors;
		
	}
	


	public double computeExitAverageQualityConcept(){
		
		return firmList.get(firmToStore).exitList.get(entryExitID).averageQualityConcept;
		
	}

	
	
	// this creates the main method to start the model from the BachCobweb class
	public static void main(String [] args){

		BachParallelModel m =new BachParallelModel();
		
	// read the input arguments (Number of iterations and deviating model parameters)
		
		
		Field[] fields = Model.class.getDeclaredFields();
		
		for(int i = 0; i < args.length;i++){
			
			if(i==0){
				
				
				m.numOfTimeSteps = Integer.parseInt(args[i]);
				
			}
			
			
			for(int j=0; j < fields.length; j++){

				if(fields[j].getName().equals(args[i])){
				
				
				try {
					
					if(fields[j].getType().getName()=="int"){
						fields[j].setInt(m,Integer.parseInt(args[i+1]));
					}else if(fields[j].getType().getName()=="double"){
						fields[j].setDouble(m,Double.parseDouble(args[i+1]));
					}else if(fields[j].getType().getName()=="boolean"){
						
						if(args[i+1].equals("true")){
							
							fields[j].setBoolean(m,true);
							
						}else{
							
							fields[j].setBoolean(m,false);
						}
					}
		
				} catch (IllegalArgumentException e) {
					e.printStackTrace();
				} catch (IllegalAccessException e) {
					e.printStackTrace();
				}
				}
				
			}
			
			if(i==(args.length-1)){
				
				m.currentPathXML = args[i];
						
			}
			
		}
		
	
		System.out.println(args[0]);
		
		System.out.println(args[1]);
		System.out.println(m.enteringCosts);
		
		if(Model.monteCarloMedian)
			System.out.println("Monte calo mode: median");
		else
			System.out.println("Monte calo mode: mean");
		
		SimInit init = new  SimInit();
		init.loadModel(m, null, true);
	

	}

}
