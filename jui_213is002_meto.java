

import java.lang.*;

class Criteria{
    public double latency;
    public double energyConsumption;
    public Criteria(double x,double y){
        latency=x;
        energyConsumption=y;
    }
}

class PAIR{
    public double val;
    public int no;
    public PAIR(double x,int y){
        val=x;
        no=y;
    }
}


public class Main {
	static Criteria[] fogDevices;
	static Criteria[] tasks;
//	static Vector<Vector<Vector<double>>> decisionMatrix; //decision matrix A = T  U  F
//	static Vector<Vector<double>> weightMatrix;
//	static Vector<Vector<double>> performance;  
    static double[][][] decisionMatrixTask;
    static double[][][] decisionMatrixFN;
    static double[][] weightMatrixTask;
    static double[][] weightMatrixFN;
    static double[][] performanceTask;
    static double[][] performanceFN;
    static int noOfFN;
	static int noOfTask;
	
	public static void main(String[] args) {
        int n=3;
        int m=4;
        noOfTask=n;
        noOfFN=m;
        tasks=new Criteria[n];
	    fogDevices=new Criteria[m];
	    double x=10.5;
        double y=5.5;
        for(int i=0;i<n;i++){
            //double x=10.5;
            //double y=5.5;
            x=Math.abs(Math.random());
            y=Math.abs(Math.random());
            Criteria item = new Criteria(x ,y);
            tasks[i]=item;
        }
        for(int i=0;i<m;i++){
            x=Math.abs(Math.random());
            y=Math.abs(Math.random());
            Criteria item = new Criteria(x ,y);
            fogDevices[i]=item;
        }
		decisionMatrixTask=new double[n][m][2];
		decisionMatrixFN=new double[m][n][2];
		
		weightMatrixTask=new double[n][2];
		weightMatrixFN=new double[m][2];
		
		performanceTask=new double[n][m];
		performanceFN=new double[m][n];
		
		createDecisionMatrix();
		print1();
		
		CRITIC();
		print2();
		
		TOPSIS();
		print3();
		
		matching();
		
	}
	
	
	

	private static void createDecisionMatrix(){
		

		// creating decision matrix for Tasks of IOT devices
		// matrix dTask size = noOfFN x 2
		
		for(int k=0;k<noOfTask;k++){  
		    double[][] dTask=new double[noOfFN][2];  
			for(int i=0;i<noOfFN;i++){
				
				dTask[i][0]=fogDevices[i].latency+Math.abs(Math.random());
				dTask[i][1]=fogDevices[i].energyConsumption+Math.abs(Math.random());
			}
			decisionMatrixTask[k]=dTask;
		}
		

		// creating decision matrix for FN
		// matrix dFN size = noOfTask x 2
		for(int k=0;k<noOfFN;k++){
			
			double[][] dFN=new double[noOfTask][2];  
			for(int i=0;i<noOfTask;i++){
				
				dFN[i][0]=tasks[i].latency+Math.abs(Math.random());
				dFN[i][1]=tasks[i].energyConsumption+Math.abs(Math.random());
			}
			decisionMatrixFN[k]=dFN;
		}

	}
	
    private static double max(double x,double y){
	    if(x>y) return x;
	    else return y;
	}
	private static double min(double x,double y){
	    if(x<y) return x;
	    else return y;
	}



	private static void CRITIC(){

    	

    	for(int l=0;l<noOfTask;l++){
        	
        	int n=noOfFN;
        	int c=2;
            double[][] B=new double[n][2];
            B=decisionMatrixTask[l];
        	// Normalize the decision matrix of an agent ‘a’ as per Eq. (16).
        
        	double[] best=new double[c];
        	double[] worst=new double[c];
        	for(int j=0;j<c;j++){
            	best[j]=B[0][j];
            	worst[j]=B[0][j];
            	for(int i=0;i<n;i++){
                	best[j]=max(best[j],B[i][j]);
                	worst[j]=min(worst[j],B[i][j]);
            	}
        	}
        	
        	for(int j=0;j<c;j++){
   		 		for(int i=0;i<n;i++){
                	B[i][j]=(B[i][j]-worst[j])/(best[j]-worst[j]);  // Executing Eq 16
            	}
        	}
            
            

        	//Evaluate the standard deviation σk for each criterion in the normalized decision matrix
        	double[] SD=new double[c];
        	for(int j=0;j<c;j++){
            	
            	double[] v=new double[n];
            	for(int i=0;i<n;i++) v[i]=B[i][j];
            	SD[j]=calculateSD(v,n);
        	}
            

        	// constructing a criteria correlation matrix which is a symmetric matrix of size  c*c
        	
        
        	double[][] S=new double[c][c];
        	for(int i=0;i<c;i++){
        		
            	for(int j=0;j<c;j++){
            		double[] v1=new double[n];
            		double[] v2=new double[n];
                	for(int x=0;x<n;x++) v1[x]=(B[x][i]);
                	for(int x=0;x<n;x++) v2[x]=(B[x][j]);
                	S[i][j]=find_coefficient(v1,v2,n);
            	}
        	}
         

           // Determine each criterion weight wk calculated as per Eq. (18) and form the criteria weight vector Wa for agent ‘a’
        	
        	double[] w=new double[c];
        	float total_sum=0;
        	for(int i=0;i<c;i++){
            	float sum=0;
            	for(int j=0;j<c;j++)
                	sum+=(1-S[i][j]);
            	w[i]=SD[i]*sum;        // Executing Eq (17)
	     	}
	     	
	     
			for(int i=0;i<c;i++)
				total_sum+=w[i];

        	for(int i=0;i<c;i++)
            	w[i]=(w[i]/total_sum);   //  Executing Eq (18)
            	

            //Add Wa as next row in W.
   			weightMatrixTask[l]=w;     	        
    	}
    	
    	
// --------------------------------------------------------------------    	
    	
    	
    	for(int l=0;l<noOfFN;l++){
        	int n=noOfTask;
        	int c=2;
            double[][] B=new double[n][2];
            B=decisionMatrixFN[l];
        	// Normalize the decision matrix of an agent ‘a’ as per Eq. (16).
        	
        	double[] best=new double[c];
        	double[] worst=new double[c];
        	for(int j=0;j<c;j++){
            	best[j]=B[0][j];
            	worst[j]=B[0][j];
            	for(int i=0;i<n;i++){
                	best[j]=max(best[j],B[i][j]);
                	worst[j]=min(worst[j],B[i][j]);
            	}
        	}
        	for(int j=0;j<c;j++){
   		 		for(int i=0;i<n;i++){
                	B[i][j]=(B[i][j]-worst[j])/(best[j]-worst[j]);  // Executing Eq 16
            	}
        	}


        	//Evaluate the standard deviation σk for each criterion in the normalized decision matrix
        	double[] SD=new double[c];
        	for(int j=0;j<c;j++){
            	
            	double[] v=new double[n];
            	for(int i=0;i<n;i++) v[i]=B[i][j];
            	SD[j]=calculateSD(v,n);
        	}


        	// constructing a criteria correlation matrix which is a symmetric matrix of size  c*c
        	
        
        	double[][] S=new double[c][c];
        	for(int i=0;i<c;i++){
        		
            	for(int j=0;j<c;j++){
            		double[] v1=new double[n];
            		double[] v2=new double[n];
                	for(int x=0;x<n;x++) v1[x]=(B[x][i]);
                	for(int x=0;x<n;x++) v2[x]=(B[x][j]);
                	S[i][j]=find_coefficient(v1,v2,n);
            	}
        	}


           // Determine each criterion weight wk calculated as per Eq. (18) and form the criteria weight vector Wa for agent ‘a’
        	
        	double[] w=new double[c];
        	float total_sum=0;
        	for(int i=0;i<c;i++){
            	float sum=0;
            	for(int j=0;j<c;j++)
                	sum+=(1-S[i][j]);
            	w[i]=SD[i]*sum;        // Executing Eq (17)
	     	}
			for(int i=0;i<c;i++)
				total_sum+=w[i];

        	for(int i=0;i<c;i++)
            	w[i]=(w[i]/total_sum);   //  Executing Eq (18)

            //Add Wa as next row in W.
   			weightMatrixFN[l]=w;     	        
    	}
   	}

	private static double find_coefficient(double[] X, double[] Y, int n){
   		double sum_X = 0, sum_Y = 0, sum_XY = 0;
   		double squareSum_X = 0, squareSum_Y = 0;
   		for (int i = 0; i < n; i++){
      		sum_X = sum_X + X[i];
      		sum_Y = sum_Y + Y[i];
      		sum_XY = sum_XY + X[i] * Y[i];
      		squareSum_X = squareSum_X + X[i] * X[i];
      		squareSum_Y = squareSum_Y + Y[i] * Y[i];
   		}
   		double corr = (float)(n * sum_XY - sum_X * sum_Y) / Math.sqrt((n * squareSum_X - sum_X * sum_X) * (n * squareSum_Y - sum_Y * sum_Y));
   		return corr;
	}

	private static double calculateSD(double[] data,int n) {
  		double sum = 0.0, mean, standardDeviation = 0.0;
  		int i;

  		for(i = 0; i < n; ++i) {
    		sum += data[i];
  		}

  		mean = sum / n;

  		for(i = 0; i < n; ++i) {
    		standardDeviation += Math.pow(data[i] - mean, 2);
  		}

  		return Math.sqrt(standardDeviation / n);
	}

	private static void TOPSIS(){
    	

    	for(int l=0;l<noOfTask;l++){
       		int n=noOfFN;
        	int c=2;
            double[][] B=new double[n][2];
            B=decisionMatrixTask[l];
            
    
                
                
        	//Normalize the decision matrix using Eq. (19).
        	for(int j=0;j<c;j++){
            	double square_sum=0;
            	for(int i=0;i<n;i++) square_sum+=B[i][j]*B[i][j];
        	    square_sum=Math.sqrt(square_sum);
       	 		for(int i=0;i<n;i++) B[i][j]=(B[i][j]/square_sum); // Eq. (19).
        	}
        	
        	
        	

        	// Calculate the weighted normalized decision matrix based on Eq. (20)
        	for(int j=0;j<c;j++){
            	for(int i=0;i<n;i++) B[i][j]=B[i][j]*weightMatrixTask[l][j];  // Eq. (20)
        	}



           
                
                
        	double[] positive=new double[c];
        	double[] negative=new double[c];
        	// Evaluate the positive and negative ideal solutions for each criterion k as 
        	// minimum and maximum value of kth column weighted normalized matrix.
        	for(int j=0;j<c;j++){
            	negative[j]=positive[j]=B[0][j];
            	for(int i=0;i<n;i++){
                	negative[j]=max(negative[j],B[i][j]);
                	positive[j]=min(positive[j],B[i][j]);
            	}
        	}
            
           
                
        	double[] d_pos=new double[n];
        	double[] d_neg=new double[n];
        	double[] p=new double[n];
 			// Determine the distance of each alternative from the positive and negative ideal 
 			// solutions based on Eq.(21) and (22).
        	for(int i=0;i<n;i++){
            	for(int j=0;j<c;j++){
                	d_pos[i]+=(B[i][j]-positive[j])*(B[i][j]-positive[j]); 
                	d_neg[i]+=(B[i][j]-negative[j])*(B[i][j]-negative[j]);   
            	}
            	d_pos[i]=Math.sqrt(d_pos[i]);  // Eq.(21)
            	d_neg[i]=Math.sqrt(d_neg[i]);  // Eq.(22)
                                 
        	}
        	
        	
           
                
           
                
                
        	// Compute the performance score for each alternative following Eq. (23),
        	for(int i=0;i<n;i++){
        	 	p[i]=d_neg[i]/(d_neg[i]+d_pos[i]);     // Eq.(23)    
        	}
            
             System.out.println("p >>>>> "+" "+l);
            for(int i=0;i<n;i++)
                System.out.println(p[i]);
                
        	// ranked in decreasing order of their performance scores 
        	//Vector<Pair <float, Integer>> vv=new Vector<Pair <float, Integer>>(n);
        	PAIR[] vv=new PAIR[n];
        	for(int i=0;i<n;i++){
        		vv[i]=new PAIR(p[i],i);
        	}
        	//sort_reverse(vv,n);
        	for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    PAIR temp = new PAIR(0.0,0);
                    if (vv[j].val < vv[i].val) {
                        // Swapping
                        temp = vv[i];
                        vv[i] = vv[j];
                        vv[j] = temp;
                    }
                }
            }
            PAIR temp = new PAIR(0.0,0);
        	int i1=0,j1=n-1;
        	while(i1<=j1){
        	    temp = vv[i1];
                vv[i1] = vv[j1];
                vv[j1] = temp;
                i1++;
                j1--;
        	}
        	for(int i=0;i<n;i++)
        		p[i]=vv[i].no;
        	performanceTask[l]=p;
    	}
    	
    	
    	
//   ------------------------------------------------------------------------
        for(int l=0;l<noOfFN;l++){
       		int n=noOfTask;
        	int c=2;
            double[][] B=new double[n][2];
            B=decisionMatrixFN[l];

        	//Normalize the decision matrix using Eq. (19).
        	for(int j=0;j<c;j++){
            	double square_sum=0;
            	for(int i=0;i<n;i++) square_sum+=B[i][j]*B[i][j];
        	    square_sum=Math.sqrt(square_sum);
       	 		for(int i=0;i<n;i++) B[i][j]=(B[i][j]/square_sum); // Eq. (19).
        	}

        	// Calculate the weighted normalized decision matrix based on Eq. (20)
        	for(int j=0;j<c;j++){
            	for(int i=0;i<n;i++) B[i][j]=B[i][j]*weightMatrixFN[l][j];  // Eq. (20)
        	}


        	double[] positive=new double[c];
        	double[] negative=new double[c];
        	// Evaluate the positive and negative ideal solutions for each criterion k as 
        	// minimum and maximum value of kth column weighted normalized matrix.
        	for(int j=0;j<c;j++){
            	negative[j]=positive[j]=B[0][j];
            	for(int i=0;i<n;i++){
                	negative[j]=max(negative[j],B[i][j]);
                	positive[j]=min(positive[j],B[i][j]);
            	}
        	}

        	double[] d_pos=new double[n];
        	double[] d_neg=new double[n];
        	double[] p=new double[n];
 			// Determine the distance of each alternative from the positive and negative ideal 
 			// solutions based on Eq.(21) and (22).
        	for(int i=0;i<n;i++){
            	for(int j=0;j<c;j++){
                	d_pos[i]+=(B[i][j]-positive[j])*(B[i][j]-positive[j]); 
                	d_neg[i]+=(B[i][j]-negative[j])*(B[i][j]-negative[j]);   
            	}
            	d_pos[i]=Math.sqrt(d_pos[i]);  // Eq.(21)
            	d_neg[i]=Math.sqrt(d_neg[i]);  // Eq.(22)
                                 
        	}

        	// Compute the performance score for each alternative following Eq. (23),
        	for(int i=0;i<n;i++){
        	 	p[i]=d_neg[i]/(d_neg[i]+d_pos[i]);     // Eq.(23)    
        	}

        	// ranked in decreasing order of their performance scores 
        	//Vector<Pair <float, Integer>> vv=new Vector<Pair <float, Integer>>(n);
        	PAIR[] vv=new PAIR[n];
        	for(int i=0;i<n;i++){
        		vv[i]=new PAIR(p[i],i);
        	}
        	//sort_reverse(vv,n);
        	for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    PAIR temp = new PAIR(0.0,0);
                    if (vv[j].val < vv[i].val) {
                        // Swapping
                        temp = vv[i];
                        vv[i] = vv[j];
                        vv[j] = temp;
                    }
                }
            }
            PAIR temp = new PAIR(0.0,0);
        	int i1=0,j1=n-1;
        	while(i1<=j1){
        	    temp = vv[i1];
                vv[i1] = vv[j1];
                vv[j1] = temp;
                i1++;
                j1--;
        	}
        	for(int i=0;i<n;i++)
        		p[i]=vv[i].no;
        	performanceFN[l]=p;
    	}
    
    }


    public static void matching(){


    // Q -- Quota of each FN
   	Integer[] Q=new Integer[noOfFN];
   	for(int i=0;i<noOfFN;i++)
   	    Q[i]=1; // store some value

    // Assign  -- Assigned tasks in FN
    Integer[][] Assign=new Integer[noOfFN][noOfTask];
    	for(int i=0;i<noOfTask;i++){
    	    double[] tj=new double[noOfFN];
    	    tj=performanceTask[i];
        	int n=noOfFN;
        	for(int x=0;x<n;x++){
        		// fi∗= highest ranked FN in P(tj) to which tj has not proposed yet
        		// Send proposal to fi∗.
            	int fi=(int)tj[x];
            	if(Q[fi]>0){    // if Qi∗ > 0 then
                	Assign[fi][x]=i+100;
                	System.out.println("Task "+i+" is being assigned to Fog Device "+fi);
                	Q[fi]=Q[fi]-1;
                	break;
            	}
            	else{
                	// Reject the assignment request;
           	}
        }
    	
    	    
    	}
    	


	}
	
	
	
	private static void print1(){
	    System.out.println("Task decisionMatrix");
	    for(int i=0;i<noOfTask;i++){
	        System.out.println("Task"+i+" -----------");
	        for(int j=0;j<noOfFN;j++){
	            System.out.println(decisionMatrixTask[i][j][0]+" "+decisionMatrixTask[i][j][1]);
	        }
	        
	    }
	    
	    System.out.println("FN decisionMatrix");
	    for(int i=0;i<noOfFN;i++){
	        System.out.println("FN"+i+" -----------");
	        for(int j=0;j<noOfTask;j++){
	            System.out.println(decisionMatrixFN[i][j][0]+" "+decisionMatrixFN[i][j][1]);
	        }
	        
	    }
	}
	
	private static void print2(){
	    System.out.println("Task weightMatrixTask");
	    for(int i=0;i<noOfTask;i++){
	        System.out.println("Task"+i+" -----------");
	        System.out.println(weightMatrixTask[i][0]+" "+weightMatrixTask[i][1]);
	        
	        
	    }
	    
	    System.out.println("FN weightMatrixFN");
	    for(int i=0;i<noOfFN;i++){
	        System.out.println("FN"+i+" -----------");
	       
	        System.out.println(weightMatrixFN[i][0]+" "+weightMatrixFN[i][1]);
	        
	        
	    }
	}
	
	
	private static void print3(){
	    System.out.println("Task performanceTask");
	    for(int i=0;i<noOfTask;i++){
	        System.out.println("Task"+i+" -----------");
	        for(int j=0;j<noOfFN;j++)
	            System.out.println(performanceTask[i][j]);
	      
	        
	    }
	    
	    System.out.println("FN performanceFN");
	    for(int i=0;i<noOfFN;i++){
	        System.out.println("FN"+i+" -----------");
	       
	        for(int j=0;j<noOfTask;j++)
	            System.out.println(performanceFN[i][j]);
	        
	    }
	}

}
