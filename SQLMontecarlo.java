
import java.sql.*;

/**
 * 
 * This class is used to store the data from the MC location analysis in a data base. Only used for
 * debugging
 *
 */


public class SQLMontecarlo {

	
	
		
		String pathToDB;
		 Connection c;
		 Statement stmt = null;
		
	  public SQLMontecarlo(String path, int firmID)
	  {
		  
		  System.out.println("Start");
		  
		  pathToDB = "montecarlo_"+firmID+".db";
		  
		  
	    Connection c = null;
	    try {
	      Class.forName("org.sqlite.JDBC");
	      c = DriverManager.getConnection("jdbc:sqlite:"+pathToDB);
	    } catch ( Exception e ) {
	      System.err.println( e.getClass().getName() + ": " + e.getMessage() );
	      System.exit(0);
	    }
	   System.out.println("Opened database successfully");
	    
	    
	    try {
	    
	    	     stmt = null;
	    	 try {
				Class.forName("org.sqlite.JDBC");
			} catch (ClassNotFoundException e) {
			
				e.printStackTrace();
			}
	       
	  		stmt = c.createStatement();
	  		 String sql = "DROP TABLE IF EXISTS Switch"; 
	  	      stmt.executeUpdate(sql);
	  	      
	  	     sql = "DROP TABLE IF EXISTS Entry"; 
		      stmt.executeUpdate(sql);
		      
		      
		       sql = "DROP TABLE IF EXISTS Exit"; 
	  	      stmt.executeUpdate(sql);
	  	      
	  	     

		       sql = "DROP TABLE IF EXISTS Baseline"; 
	  	      stmt.executeUpdate(sql);
	  	      
	  	      
	  	  
		  	} catch (SQLException e1) {
		  		
		  		e1.printStackTrace();
		  	}
	    
	    
	   
	  }
	  
	  
	  
	  
	  
	  public SQLMontecarlo(String path)
	  {

		  System.out.println("Start");
		  pathToDB = path;
		  
		  
	    Connection c = null;
	    try {
	      Class.forName("org.sqlite.JDBC");
	      c = DriverManager.getConnection("jdbc:sqlite:"+pathToDB);
	    } catch ( Exception e ) {
	      System.err.println( e.getClass().getName() + ": " + e.getMessage() );
	      System.exit(0);
	    }
	   //System.out.println("Opened database successfully");
	   
	  }
	  
	  
	  
	  public void createEntryExitTable()
	  {
	 
	  
	    
	  
		  System.out.println("Start createEntryExitTable");
	    
	    
	    try {
	      Class.forName("org.sqlite.JDBC");
	      c = DriverManager.getConnection("jdbc:sqlite:"+pathToDB);
	     //System.out.println("Opened database successfully");
	      
	      
	     
	      

	      stmt = c.createStatement();
	      String sql = "CREATE TABLE EntryExit " +
	              "(_ITERATION_NO	INT ";
	      
	      
	      	
	      				sql = sql +","+"run"+"	INT";
						sql = sql +","+"t"+"	INT";
						sql = sql +","+"region_id_entry"+"	INT";
						sql = sql +","+"region_id_exit"+"	INT";
						sql = sql +","+"firm_id"+"	INT";
						
						sql = sql +","+"output"+"	REAL";
						sql = sql +","+"profit"+"	REAL";
						sql = sql +","+"npv"+"	REAL";
						sql = sql +","+"quality"+"	REAL";
						sql = sql +","+"concept_quality"+"	REAL";
						
						sql = sql +","+"max_pending_concept_quality"+"	REAL";

		
	      sql = sql+")";
	      stmt.executeUpdate(sql);
	     
	   
	    } catch ( Exception e ) {
	      System.err.println( e.getClass().getName() + ": " + e.getMessage() );

	    }
	   //System.out.println("Table Firm created successfully");
	  }
	  
	  
	  
	  
	  
	
	  
	  
	  public void createEntryTable()
	  {
	  
		  
		  System.out.println("Start createEntryTable");
	  
	    try {

	      
	      String sql = "CREATE TABLE Entry " +
	              "(_ITERATION_NO	INT ";
	      sql = sql +","+"run"+"	INT";
	      	sql = sql +","+"t"+"	INT";
			sql = sql +","+"region_id"+"	INT";
			sql = sql +","+"firm_id"+"	INT";
			
			sql = sql +","+"output"+"	REAL";
			sql = sql +","+"profit"+"	REAL";
			sql = sql +","+"npv"+"	REAL";
			sql = sql +","+"quality"+"	REAL";
			sql = sql +","+"concept_quality"+"	REAL";
			
			sql = sql +","+"max_pending_concept_quality"+"	REAL";

	      sql = sql+")";

	      stmt = c.createStatement();
	      
	      
	      stmt.executeUpdate(sql);
	     
	     
	    } catch ( Exception e ) {
	      System.err.println( e.getClass().getName() + ": " + e.getMessage() );
	      System.exit(0);
	    }
	   //System.out.println("Table Entry created successfully");
	  }
	  
	  
	  
	  
	  public void createExitTable()
	  {
	  
		  
		  System.out.println("Start createExitTable");
	  
	    try {
	     

	      stmt = c.createStatement();
	      String sql = "CREATE TABLE Exit " +
	              "(_ITERATION_NO	INT ";
	      
	      sql = sql +","+"run"+"	INT";   
  	sql = sql +","+"t"+"	INT";
	sql = sql +","+"region_id"+"	INT";
	sql = sql +","+"firm_id"+"	INT";
	
	sql = sql +","+"output"+"	REAL";
	sql = sql +","+"profit"+"	REAL";
	sql = sql +","+"npv"+"	REAL";
	sql = sql +","+"quality"+"	REAL";
	sql = sql +","+"concept_quality"+"	REAL";
	
	sql = sql +","+"max_pending_concept_quality"+"	REAL";
	      sql = sql+")";
	      
	      
	      stmt.executeUpdate(sql);
	     
	    
	    } catch ( Exception e ) {
	      System.err.println( e.getClass().getName() + ": " + e.getMessage() );
	      System.exit(0);
	    }
	   //System.out.println("Table Exit created successfully");
	  }
	  
	  
	  
	  
	  public void createBaselineTable()
	  {
	  
	  
		  System.out.println("Start createBaselineTable");
		  
	    try {

	      
	      String sql = "CREATE TABLE Baseline " +
	              "(_ITERATION_NO	INT ";
	      sql = sql +","+"run"+"	INT";
	      	sql = sql +","+"t"+"	INT";
			sql = sql +","+"firm_id"+"	INT";
			
			sql = sql +","+"output"+"	REAL";
			sql = sql +","+"profit"+"	REAL";
			sql = sql +","+"npv"+"	REAL";
			sql = sql +","+"quality"+"	REAL";
			sql = sql +","+"concept_quality"+"	REAL";
			
			sql = sql +","+"max_pending_concept_quality"+"	REAL";

	      sql = sql+")";

	      stmt = c.createStatement();
	      
	      
	      stmt.executeUpdate(sql);
	     
	     
	    } catch ( Exception e ) {
	      System.err.println( e.getClass().getName() + ": " + e.getMessage() );
	      System.exit(0);
	    }
	   //System.out.println("Table Entry created successfully");
	  }
	  
	 
	  
	  
	  public void insertEntryExit(int run, int iteration, int time, int regionIDEntry, int regionIDExitx, double npv , Firm aFirm  )
	  {
	    
	  
	    try {
	   
	      c.setAutoCommit(false);
	     //System.out.println("Opened database successfully");

	      stmt = c.createStatement();
	      
	      String sql = "INSERT INTO EntryExit (_ITERATION_NO,run, t, region_id_entry, region_id_exit, firm_id, output, profit, npv, quality,concept_quality,max_pending_concept_quality";
	      
	     
	  	
	  	sql = sql +") VALUES ("+iteration;
	  	
	  	sql = sql +","+run+","+time+","+regionIDEntry+","+regionIDExitx+","+aFirm.firmID+","+aFirm.equilQuantity+","+aFirm.equilProfit+","+npv+","+aFirm.quality+","+aFirm.qualityConcept+","+aFirm.maxPendingQualityConcept;   
	      
	  		
	  	
		sql = sql +");";
	      
	      
	     
	      stmt.executeUpdate(sql);

	      

	     
	      //c.commit();
	     
	    } catch ( Exception e ) {
	      System.err.println( e.getClass().getName() + ": " + e.getMessage() );
	      System.exit(0);
	    }
	   //System.out.println("Records created successfully");
	  }
	  
	  
	  
	  
	 
	  
	  
	  
	  public void insertEntry(int run,int iteration, int time,  int regionID, double npv , Firm aFirm  )
	  {
	  
	  
	    try {
	     
	      c.setAutoCommit(false);
	     //System.out.println("Opened database successfully");

	      stmt = c.createStatement();
	      
	      
	      String sql = "INSERT INTO Entry (_ITERATION_NO, run, t, region_id, firm_id, output, profit, npv, quality,concept_quality,max_pending_concept_quality";
	      
		     
		  	
		  	sql = sql +") VALUES ("+iteration;
		  	
		  	sql = sql +","+run +","+time+","+regionID+","+aFirm.firmID+","+aFirm.equilQuantity+","+aFirm.equilProfit+","+npv+","+aFirm.quality+","+aFirm.qualityConcept+","+aFirm.maxPendingQualityConcept;   
		      
	  	
		sql = sql +");";
	      stmt.executeUpdate(sql);

	      

	     
	      //c.commit();
	   
	    } catch ( Exception e ) {
	      System.err.println( e.getClass().getName() + ": " + e.getMessage() );
	      System.exit(0);
	    }
	   //System.out.println("Records created successfully");
	  }
	  
	  
	  

	  public void insertExit(int run, int iteration, int time,  int regionID, double npv , Firm aFirm )
	  {
	    //Connection c = null;
	  
	    try {
	      //Class.forName("org.sqlite.JDBC");
	      //c = DriverManager.getConnection("jdbc:sqlite:"+pathToDB+"/iters.db");
	      c.setAutoCommit(false);
	     //System.out.println("Opened database successfully");

	      stmt = c.createStatement();
	      
	      

	      String sql = "INSERT INTO Exit (_ITERATION_NO ,run, t, region_id, firm_id, output, profit, npv, quality,concept_quality,max_pending_concept_quality";
	      
		     
		  	
		  	sql = sql +") VALUES ("+iteration;
		  	
		  	sql = sql  +","+run+","+time+","+regionID+","+aFirm.firmID+","+aFirm.equilQuantity+","+aFirm.equilProfit+","+npv+","+aFirm.quality+","+aFirm.qualityConcept+","+aFirm.maxPendingQualityConcept;   
		      
	  	
	  	
		sql = sql +");";
	      stmt.executeUpdate(sql);

	      

	     
	      //c.commit();
	      //c.close();
	    } catch ( Exception e ) {
	      System.err.println( e.getClass().getName() + ": " + e.getMessage() );
	      System.exit(0);
	    }
	   //System.out.println("Records created successfully");
	  }
	  
	  
	  
	  
	  
	  
	  
	  
	  public void insertbaseline(int run ,int iteration, int time, double npv , Firm aFirm )
	  {
	    //Connection c = null;
	  
	    try {
	      //Class.forName("org.sqlite.JDBC");
	      //c = DriverManager.getConnection("jdbc:sqlite:"+pathToDB+"/iters.db");
	      c.setAutoCommit(false);
	     //System.out.println("Opened database successfully");

	      stmt = c.createStatement();
	      
	      

	      String sql = "INSERT INTO Baseline (_ITERATION_NO, run, t, firm_id, output, profit, npv, quality,concept_quality,max_pending_concept_quality";
	      
		     
		  	
		  	sql = sql +") VALUES ("+iteration;
		  	
		  	sql = sql  +","+run+","+time+","+aFirm.firmID+","+aFirm.equilQuantity+","+aFirm.equilProfit+","+npv+","+aFirm.quality+","+aFirm.qualityConcept+","+aFirm.maxPendingQualityConcept;   
		      
	  	
	  	
		sql = sql +");";
	      stmt.executeUpdate(sql);

	      

	     
	      //c.commit();
	      //c.close();
	    } catch ( Exception e ) {
	      System.err.println( e.getClass().getName() + ": " + e.getMessage() );
	      System.exit(0);
	    }
	   //System.out.println("Records created successfully");
	  }
	
	  

	  
	  void commit(){
		  
		  
	 try{
			c.commit();
			  
		  }catch ( Exception e ) {
		      System.err.println( e.getClass().getName() + ": " + e.getMessage() );
		      System.exit(0);
		    }
		  
		  
	  }
	  void atEnd(){
		  
		  
		  try{
			  
			  stmt.close();
			  c.close();
			  
		  }catch ( Exception e ) {
		      System.err.println( e.getClass().getName() + ": " + e.getMessage() );
		      System.exit(0);
		    }
		  
	  }
	  
	}
	
	
	
	
	

