import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

/**
 * This class is used to write certain information in txt files. Only used for debugging
 * */
public class FirmWriter {
	
	
	ArrayList<File> writerList = new ArrayList<File>();
	
	FirmWriter(int numFirms){
		
		
		for(int i=0; i < numFirms;i++){
		
			File file = new File("Firm"+(i+1)+".txt");
			writerList.add(file);
			
			
			if(writerList.get(i).exists()){
		        System.out.println("Clean");
		        writerList.get(i).delete();
		    }
		
		}
		
	}
	
	
	
	void write(int id, String string){
	

	try{
	    if(!writerList.get(id-1).exists()){
	        System.out.println("We had to make a new file.");
	        writerList.get(id-1).createNewFile();
	    }

	    FileWriter fileWriter = new FileWriter(writerList.get(id-1), true);
	    
	    BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);
	    bufferedWriter.write(string);
	    bufferedWriter.close();

	   
	} catch(IOException e) {
	    System.out.println("COULD NOT LOG!!");
	}
	
	
	
	}

}
