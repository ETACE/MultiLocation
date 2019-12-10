import java.lang.reflect.Field;

/**
 * This classs is used in the BatchModel mode to write the data to SQLite data bases. We use this kind of data base in
 *  combination with R. This is an attempt to write a relatively generic method in a sense that for a single agent type there 
 *  is no need to adjust the code after introducing a new variable.
 * */
import java.sql.*;


public class SQLiteJDBC
{
	
	String pathToDB;
	 Connection c;
	 Statement stmt = null;
	
  public SQLiteJDBC(String path)
  {
	  
	  pathToDB = path;
	  
	  
    Connection c = null;
    try {
      Class.forName("org.sqlite.JDBC");
      c = DriverManager.getConnection("jdbc:sqlite:"+pathToDB+"/iters.db");
    } catch ( Exception e ) {
      System.err.println( e.getClass().getName() + ": " + e.getMessage() );
      System.exit(0);
    }
   //System.out.println("Opened database successfully");
    
    
    try {
    
    	     stmt = null;
    	 try {
			Class.forName("org.sqlite.JDBC");
		} catch (ClassNotFoundException e) {
			
			e.printStackTrace();
		}
       
  		stmt = c.createStatement();
  		 String sql = "DROP TABLE IF EXISTS Firm"; 
  	      stmt.executeUpdate(sql);
  	      
  	     sql = "DROP TABLE IF EXISTS Entry"; 
	      stmt.executeUpdate(sql);
	      
	      
	       sql = "DROP TABLE IF EXISTS Exit"; 
  	      stmt.executeUpdate(sql);
  	      
  	    sql = "DROP TABLE IF EXISTS AggregatedData"; 
	      stmt.executeUpdate(sql);
	      
	      sql = "DROP TABLE IF EXISTS Location"; 
	      stmt.executeUpdate(sql);
  	     
  	  
	  	} catch (SQLException e1) {
	  		
	  		e1.printStackTrace();
	  	}
    
    
   
  }

  public void createFirmTable()
  {

    try {
      Class.forName("org.sqlite.JDBC");
      c = DriverManager.getConnection("jdbc:sqlite:"+pathToDB+"/iters.db");
     //System.out.println("Opened database successfully");
      
      
     
      

      stmt = c.createStatement();
      String sql = "CREATE TABLE Firm " +
              "(_ITERATION_NO	INT ";
      
      
  	Field[] fields = Firm.class.getDeclaredFields();
	
	for(int i=0; i < fields.length; i++){

		if(fields[i].getType().getName().equals("int") ){
			sql = sql +","+fields[i].getName()+"	INT";
		}
		else if( fields[i].getType().getName().equals("double")){
			sql = sql +","+fields[i].getName()+"	REAL";
		}else if(fields[i].getType().getName().equals("boolean")){
			
			sql = sql +","+fields[i].getName()+"	INT"; 
			
		}
		

	}
      sql = sql+")";
      stmt.executeUpdate(sql);
     
   
    } catch ( Exception e ) {
      System.err.println( e.getClass().getName() + ": " + e.getMessage() );

    }
   //System.out.println("Table Firm created successfully");
  }
  
  
  
  
  
  public void createLocationTable()
  {

    try {
      Class.forName("org.sqlite.JDBC");
      c = DriverManager.getConnection("jdbc:sqlite:"+pathToDB+"/iters.db");
     //System.out.println("Opened database successfully");


      stmt = c.createStatement();
      String sql = "CREATE TABLE Location " +
              "(_ITERATION_NO	INT ";
      
      
  	Field[] fields = Location.class.getDeclaredFields();
	
	for(int i=0; i < fields.length; i++){

		if(fields[i].getType().getName().equals("int") ){
			sql = sql +","+fields[i].getName()+"	INT";
		}
		else if( fields[i].getType().getName().equals("double")){
			sql = sql +","+fields[i].getName()+"	REAL";
		}else if(fields[i].getType().getName().equals("boolean")){
			
			sql = sql +","+fields[i].getName()+"	INT"; 
			
		}
		

	}
	
	
	sql = sql +",firms   TEXT"; 
	
      sql = sql+")";
      stmt.executeUpdate(sql);
     
   
    } catch ( Exception e ) {
      System.err.println( e.getClass().getName() + ": " + e.getMessage() );

    }
   //System.out.println("Table Firm created successfully");
  }
  
  
  
  public void createEngryTable()
  {
  
  
    try {

      
      String sql = "CREATE TABLE Entry " +
              "(_ITERATION_NO	INT, "+
    		  "firmID 	INT " ;
      
      
  	Field[] fields = Firm.EntryCharacteristics.class.getDeclaredFields();
	
	for(int i=0; i < fields.length; i++){

		if(fields[i].getType().getName().equals("int") ){
			sql = sql +","+fields[i].getName()+"	INT";
		}
		else if( fields[i].getType().getName().equals("double")){
			sql = sql +","+fields[i].getName()+"	REAL";
		}else if(fields[i].getType().getName().equals("boolean")){
			
			sql = sql +","+fields[i].getName()+"	INT"; 
			
		}
		

	}
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
  
  
    try {
     

      stmt = c.createStatement();
      String sql = "CREATE TABLE Exit " +
              "(_ITERATION_NO	INT, "+
              "firmID 	INT" ;
      
      
  	Field[] fields = Firm.ExitCharacteristics.class.getDeclaredFields();
	
	for(int i=0; i < fields.length; i++){

		if(fields[i].getType().getName().equals("int") ){
			sql = sql +","+fields[i].getName()+"	INT";
		}
		else if( fields[i].getType().getName().equals("double")){
			sql = sql +","+fields[i].getName()+"	REAL";
		}else if(fields[i].getType().getName().equals("boolean")){
			
			sql = sql +","+fields[i].getName()+"	INT"; 
			
		}
		

	}
      sql = sql+")";
      
      
      stmt.executeUpdate(sql);
     
    
    } catch ( Exception e ) {
      System.err.println( e.getClass().getName() + ": " + e.getMessage() );
      System.exit(0);
    }
   //System.out.println("Table Exit created successfully");
  }
  
  
  
  public void createAggregateTable()
  {
  
  
    
  
   
    
    
    try {
     
      
      stmt = c.createStatement();
    
      String sql = "CREATE TABLE AggregatedData " +
              "(_ITERATION_NO	INT ";
      
      
  	Field[] fields = Model.class.getDeclaredFields();
	
	for(int i=0; i < fields.length; i++){

		if(fields[i].getType().getName().equals("int") ){
			sql = sql +","+fields[i].getName()+"	INT";
		}
		else if( fields[i].getType().getName().equals("double")){
			sql = sql +","+fields[i].getName()+"	REAL";
		}else if(fields[i].getType().getName().equals("long")){
			
			sql = sql +","+fields[i].getName()+"	TEXT"; 
			
		}else if(fields[i].getType().getName().equals("boolean")){
			
			sql = sql +","+fields[i].getName()+"	INT"; 
			
		}
		

	}
	
	
	
	
      sql = sql+")";
      stmt.executeUpdate(sql);
     
    
    } catch ( Exception e ) {
      System.err.println( e.getClass().getName() + ": " + e.getMessage() );

    }
   //System.out.println("Table Firm created successfully");
  }
  
  
  public void insertFirm(int iteration, Firm aFirm  )
  {
    
  
    try {
   
      c.setAutoCommit(false);
     //System.out.println("Opened database successfully");

      stmt = c.createStatement();
      
      String sql = "INSERT INTO Firm (_ITERATION_NO";
      
      Field[] fields = Firm.class.getDeclaredFields();
  	
  	for(int i=0; i < fields.length; i++){

  		if(fields[i].getType().getName().equals("int") || fields[i].getType().getName().equals("double") || fields[i].getType().getName().equals("boolean")){
  			sql = sql +","+fields[i].getName();
      
  		}
  	}
  	
  	sql = sql +") VALUES ("+iteration;
  	
  	for(int i=0; i < fields.length; i++){

  		if(fields[i].getType().getName().equals("int") || fields[i].getType().getName().equals("double")){
  			sql = sql +","+fields[i].get(aFirm);
      
  		}else if(fields[i].getType().getName().equals("boolean")){
  			
  			//System.out.println(fields[i].get(aFirm));
  			
  			if(fields[i].get(aFirm).equals(true)){
  				
  				sql = sql +", 1";
  				
  			}else{
  				
  				sql = sql +", 0";
  			}
  			
  			
  		}
  	}
  	
	sql = sql +");";
      
      
     
      stmt.executeUpdate(sql);

      

     
      //c.commit();
     
    } catch ( Exception e ) {
      System.err.println( e.getClass().getName() + ": " + e.getMessage() );
      System.exit(0);
    }
   //System.out.println("Records created successfully");
  }
  
  
  
  
  
  
  public void insertLocation(int iteration, Location aLocation  )
  {
    
  
    try {
   
      c.setAutoCommit(false);
     //System.out.println("Opened database successfully");

      stmt = c.createStatement();
      
      String sql = "INSERT INTO Location (_ITERATION_NO";
      
      Field[] fields = Location.class.getDeclaredFields();
  	
  	for(int i=0; i < fields.length; i++){

  		if(fields[i].getType().getName().equals("int") || fields[i].getType().getName().equals("double")|| fields[i].getType().getName().equals("boolean")){
  			sql = sql +","+fields[i].getName();
      
  		}
  	}
  	
  	sql = sql + ", firms";
  	
  	sql = sql +") VALUES ("+iteration;
  	
  	for(int i=0; i < fields.length; i++){

  		if(fields[i].getType().getName().equals("int") || fields[i].getType().getName().equals("double")){
  			sql = sql +","+fields[i].get(aLocation);
      
  		}else if(fields[i].getType().getName().equals("boolean")){
  			
  			if(fields[i].get(aLocation).equals(true)){
  				
  				sql = sql +", 1";
  				
  			}else{
  				
  				sql = sql +", 0";
  			}
  			
  			
  		}
  	}
  	
  	String str = ",'";
  	
  	for(int i=0; i < aLocation.firmList.size();i++){
  		
  		if(i==0)
  			str = str+Integer.toString(aLocation.firmList.get(i).firmID);
  		else
  			str = str+","+Integer.toString(aLocation.firmList.get(i).firmID);
  	}
  	
  	
	sql = sql +str+"');";
	
	
      
      
     
      stmt.executeUpdate(sql);

      

     
      //c.commit();
     
    } catch ( Exception e ) {
      System.err.println( e.getClass().getName() + ": " + e.getMessage() );
      System.exit(0);
    }
   //System.out.println("Records created successfully");
  }
  
  
  
  
  
  public void insertEntry(int iteration, int firmID, Firm.EntryCharacteristics entry  )
  {
  
  
    try {
     
      c.setAutoCommit(false);
     //System.out.println("Opened database successfully");

      stmt = c.createStatement();
      
      
      
      String sql = "INSERT INTO Entry (_ITERATION_NO,firmID";
      
      Field[] fields = Firm.EntryCharacteristics.class.getDeclaredFields();
  	
  	for(int i=0; i < fields.length; i++){

  		if(fields[i].getType().getName().equals("int") || fields[i].getType().getName().equals("double")|| fields[i].getType().getName().equals("boolean")){
  			sql = sql +","+fields[i].getName();
      
  		}
  	}
  	
  	sql = sql +") VALUES ("+iteration+","+firmID;
  	
  	for(int i=0; i < fields.length; i++){

  		if(fields[i].getType().getName().equals("int") || fields[i].getType().getName().equals("double")){
  			sql = sql +","+fields[i].get(entry);
      
  		}else if(fields[i].getType().getName().equals("boolean")){
  			
  			if(fields[i].get(entry).equals(true)){
  				
  				sql = sql +", 1";
  				
  			}else{
  				
  				sql = sql +", 0";
  			}
  			
  			
  		}
  	}
  	
	sql = sql +");";
      stmt.executeUpdate(sql);

      

     
      //c.commit();
   
    } catch ( Exception e ) {
      System.err.println( e.getClass().getName() + ": " + e.getMessage() );
      System.exit(0);
    }
   //System.out.println("Records created successfully");
  }
  
  
  

  public void insertExit(int iteration, int firmID, Firm.ExitCharacteristics exit )
  {
    //Connection c = null;
  
    try {
      //Class.forName("org.sqlite.JDBC");
      //c = DriverManager.getConnection("jdbc:sqlite:"+pathToDB+"/iters.db");
      c.setAutoCommit(false);
     //System.out.println("Opened database successfully");

      stmt = c.createStatement();
      
      
      String sql = "INSERT INTO Exit (_ITERATION_NO,firmID";
      
      Field[] fields = Firm.ExitCharacteristics.class.getDeclaredFields();
  	
  	for(int i=0; i < fields.length; i++){

  		if(fields[i].getType().getName().equals("int") || fields[i].getType().getName().equals("double")|| fields[i].getType().getName().equals("boolean")){
  			sql = sql +","+fields[i].getName();
      
  		}
  	}
  	
  	sql = sql +") VALUES ("+iteration+","+firmID;
  	
  	for(int i=0; i < fields.length; i++){

  		if(fields[i].getType().getName().equals("int") || fields[i].getType().getName().equals("double")){
  			sql = sql +","+fields[i].get(exit);
      
  		}else if(fields[i].getType().getName().equals("boolean")){
  			
  			if(fields[i].get(exit).equals(true)){
  				
  				sql = sql +", 1";
  				
  			}else{
  				
  				sql = sql +", 0";
  			}
  			
  			
  		}
  	}
  	
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
  
  
  
  public void insertAggregatedData(int iteration , Model m)
  {
  
  
    try {
     
      c.setAutoCommit(false);
     //System.out.println("Opened database successfully");

      stmt = c.createStatement();
  
      
      String sql = "INSERT INTO AggregatedData (_ITERATION_NO";
      
      Field[] fields = Model.class.getDeclaredFields();
  	
  	for(int i=0; i < fields.length; i++){

  		if(fields[i].getType().getName().equals("int") || fields[i].getType().getName().equals("double")|| fields[i].getType().getName().equals("long")|| fields[i].getType().getName().equals("boolean")){
  			sql = sql +","+fields[i].getName();
  		
  			
  			
  			
  			//System.out.println(fields[i].getType().getName());
			//System.out.println(fields[i].get(m));
      
  		}
  		
  	
  	}
  	
  	sql = sql +") VALUES ("+iteration;
  	
  	for(int i=0; i < fields.length; i++){

  		if(fields[i].getType().getName().equals("int") || fields[i].getType().getName().equals("double")|| fields[i].getType().getName().equals("long")){
  			sql = sql +","+fields[i].get(m);
  			
      
  		}else if(fields[i].getType().getName().equals("boolean")){
  			
  			//System.out.println(fields[i].get(m));
  			
  			if(fields[i].get(m).equals(true)){
  				
  				sql = sql +", 1";
  				
  			}else{
  				
  				sql = sql +", 0";
  			}
  			
  			
  		}
  	}
  	sql = sql +");";
    
      stmt.executeUpdate(sql);

      

     
      //c.commit();
   
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