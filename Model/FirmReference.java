/**
 * This class is used by firms to hold information about competitors
 *
 */
public class FirmReference implements Cloneable {
    

	
	
	/****************Cloning   to avoid problems in deep copying************************/
	
	public FirmReference clone()  {
        try {
			return (FirmReference) super.clone();
		} catch (CloneNotSupportedException e) {
			
			e.printStackTrace();
			return null;
		}
		
        
	}
	
	
	int firmID;
	double quality;
	double qualityConcept;
	double output;
	boolean innovator = false ;
	double maxPendingQualityConcept;

	
	FirmReference(int id, double qual, double qualConcept, boolean inno, double maxPendingQualConc){
		
		firmID = id;
		quality = qual;
		
		qualityConcept= qualConcept;
		innovator = inno;
		maxPendingQualityConcept=  maxPendingQualConc;
		
		
	}
	


	FirmReference(int id, Object nu,  double qualConcept, boolean inno, double maxPendingQualConc){
	
	firmID = id;
	qualityConcept= qualConcept;
	innovator = inno;
	maxPendingQualityConcept=  maxPendingQualConc;
	}
	
}