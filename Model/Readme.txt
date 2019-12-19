This java code is the implementation of the model described in the paper "R&D location strategies" authored by Luca Colombo, Herbert
Dawid and Philipp Harting. 

The model uses the Repast 3.1 libraries as a backbone, and Model.java is the top level java file of the implementation. The model can be executed in three different modes:

	1. Executing the model from Model.java runs the model without any further analysis of the simulation data. In this case, no data is collected and only the console output is displayed. This mode is useful for debugging, where one wants to run the model without the data processing classes (see below) that are not directly linked to the model.
	
	2. Executing the model in the GUI mode. If the model is launched from the ModelGUI class, then one can run a single run using a GUI. The GUI can thereby be used to change the parameters, and to start, pause and stop the simulation. Furthermore, during the simulations, selected variables are graphically displayed in form of time series, network graphics and other useful graphs to illustrate the spatial interaction of firms.
	
	3. Executing the model in batch mode. If the model is run from BatchParallelModel.java, the model is executed for a specified number of iterations and writes the output (all integer, double and boolean variables of different classes) into a SQLite database. If the java executable is launched from the command line, it accepts the number of iterations, the parameters, for which the defaul values are not used, as as well as the path to which the database should be written, as inputs from the command line.
	
One should note that the implementation features several java files with different roles. First, there are the three model files:
	-Model.java
	-ModelGUI.java
	-BatchParallelModel.java
	
In particular, the Model.java file is important as this one defines the structure of the model and controls which activities are executed from where. So, in order to understand the sequence of activities of the model, it is best to start from the Model.java file.	
	
Then we have the agent classes
	-Firm.java
	-Location.java
	
where one can see how the agent specific activities are implemented. 

Furthermore, there are some classes used to store certain information, especially used by the Firm agent:
	-FirmLocation.java
	-Firmrefernce.java
	-ImiContainer.java
	-ImitationScheme.java
	-InnoContainer.java
	-InnovationItem.java
	-InnovationScheme.java,

some classes that hold procedures executed by the firm, but thar are too complex to have them within the Firm class
	-InnovationProcess.java
	-MonteCarloGeneric.java
 	-MonteCarloSimulatorSerial.java
 
and finally some more technical java classes used to store and process data but are not directly linked to the actual functioning of the model.
 	-FirmEdge.java
	-FirmNode.java
	-FirmWriter.java
	-SetGetModel.java
	-SQLiteJDBC.java
	-SQLMontecarlo.java
	    
