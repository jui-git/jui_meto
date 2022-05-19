package org.fog.test.perfeval;

import java.util.ArrayList;
import java.util.Calendar;
import java.util.LinkedList;
import java.util.List;
import org.cloudbus.cloudsim.Host;
import org.cloudbus.cloudsim.Log;
import org.cloudbus.cloudsim.Pe;
import org.cloudbus.cloudsim.Storage;
import org.cloudbus.cloudsim.core.CloudSim;
import org.cloudbus.cloudsim.power.PowerHost;
import org.cloudbus.cloudsim.provisioners.RamProvisionerSimple;
import org.cloudbus.cloudsim.sdn.overbooking.BwProvisionerOverbooking;
import org.cloudbus.cloudsim.sdn.overbooking.PeProvisionerOverbooking;
import org.fog.application.AppEdge;
import org.fog.application.AppLoop;
import org.fog.application.Application;
import org.fog.application.selectivity.FractionalSelectivity;
import org.fog.entities.Actuator;
import org.fog.entities.FogBroker;
import org.fog.entities.FogDevice;
import org.fog.entities.FogDeviceCharacteristics;
import org.fog.entities.Sensor;
import org.fog.entities.Tuple;
import org.fog.placement.Controller;
import org.fog.placement.ModuleMapping;
import org.fog.placement.ModulePlacementEdgewards;
import org.fog.placement.ModulePlacementMapping;
import org.fog.policy.AppModuleAllocationPolicy;
import org.fog.scheduler.StreamOperatorScheduler;
import org.fog.utils.FogLinearPowerModel;
import org.fog.utils.FogUtils;
import org.fog.utils.TimeKeeper;
import org.fog.utils.distribution.DeterministicDistribution;


public class METO {
	static List<FogDevice> fogDevices = new ArrayList<FogDevice>();
	static List<Task> tasks = new ArrayList<Task>();
	static Vector<Vector<Vector<float>>> decisionMatrix; //decision matrix A = T  U  F
	static Vector<Vector<float>> weightMatrix;
	static Vector<Vector<float>> performance;  
	
	public static void main(String[] args) {

		Log.printLine("Starting METO...");

		try {
			

		} catch (Exception e) {
			e.printStackTrace();
			Log.printLine("Unwanted errors happen");
		}
	}
	
	/**
	 * Creates the fog devices in the physical topology of the simulation.
	 * @param userId
	 * @param appId
	 */
	private static void createFogDevices(int userId, String appId) {

		System.out.println("enter no. of fog devices and tasks");
		int f,t;

		FogDevice cloud = createFogDevice("cloud", 44800, 40000, 100, 10000, 0, 0.01, 16*103, 16*83.25);
		cloud.setParentId(-1);
		fogDevices.add(cloud);
		FogDevice proxy = createFogDevice("proxy-server", 2800, 4000, 10000, 10000, 1, 0.0, 107.339, 83.4333);
		proxy.setParentId(cloud.getId());
		proxy.setUplinkLatency(100); // latency of connection between proxy server and cloud is 100 ms
		fogDevices.add(proxy);
		for(int i=0;i<numOfAreas;i++){
			addArea(i+"", userId, appId, proxy.getId());
		
	}
		//array list of fog devices
		//each fog device has parameters latency, energy ,no. of vrus


		
	}
	private static void createIotDevices(int userId, String appId) {
		//array list of tasks
		//each task has input size , output size, computation cycles, deadline
		FogDevice camera = createFogDevice("m-"+id, 500, 1000, 10000, 10000, 3, 0, 87.53, 82.44);
		camera.setParentId(parentId);
		Sensor sensor = new Sensor("s-"+id, "CAMERA", userId, appId, new DeterministicDistribution(5)); // inter-transmission time of camera (sensor) follows a deterministic distribution
		sensors.add(sensor);
		Actuator ptz = new Actuator("ptz-"+id, userId, appId, "PTZ_CONTROL");
		actuators.add(ptz);
		sensor.setGatewayDeviceId(camera.getId());
		sensor.setLatency(1.0);  // latency of connection between camera (sensor) and the parent Smart Camera is 1 ms
		ptz.setGatewayDeviceId(camera.getId());
		ptz.setLatency(1.0);  // latency of connection between PTZ Control and the parent Smart Camera is 1 ms
		return camera;

		
	}


	private static void createDecisionMatrix(){
		int noOfFN=fogDevices.size();
		int noOfTask=tasks.size();

		// creating decision matrix for Tasks of IOT devices
		// matrix dTask size = noOfFN x 2
		
		for(int k=0;k<noOfTask;k++){  
			Vector<Vector<float>> dTask= new Vector<Vector<float>>(noOfFN);   
			for(int i=0;i<noOfFN;i++){
				Vector<float> row=new Vector<float>(2);
				row[0]=fogDevices[i].latency;
				row[1]=fogDevices[i].energyConsumption;
				dTask[i]=row;
			}
			decisionMatrix.add(dTask);
		}
		

		// creating decision matrix for FN
		// matrix dFN size = noOfTask x 2
		for(int k=0;k<noOfFN;k++){
			Vector<Vector<float>> dFN= new Vector<Vector<float>>(noOfTask);     
			for(int i=0;i<noOfTask;i++){
				Vector<float> row=new vector<float>(2);
				row[0]=tasks[i].energyConsumption;
				row[1]=tasks[i].deadLine;
				dFN[i]=row;
			}
			decisionMatrix.add(dFN);
		}

	}




	private static void CRITIC(){

    	int m=decisionMatrix.size();//m= np. of task nodes + no. of fog nodes

    	for(int l=0;l<m;l++){
        	Vector<Vector<float>> B=decisionMatrix[l]; //taking each agent 'a' belongs to A=T U F
        	int n=B.size();
        	int c=B[0].size();

        	// Normalize the decision matrix of an agent ‘a’ as per Eq. (16).
        	Vector<float> best=new Vector<float>(c);
        	Vector<float> worst=new Vector<float>(c);
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
                	B[i][j]=(B[i][j]-worst[j])/(best[j]-best[j]);  // Executing Eq 16
            	}
        	}


        	//Evaluate the standard deviation σk for each criterion in the normalized decision matrix
        	Vector<float> SD=new Vector<float>(c);
        	for(int j=0;j<c;j++){
            	Vector<float> v=new Vector<float>(n);
            	for(int i=0;i<n;i++) v[i]=B[i][j];
            	SD[j]=calculateSD(v);
        	}


        	// constructing a criteria correlation matrix which is a symmetric matrix of size  c*c
        	
        	Vector<Vector<float>> S=new Vector<Vector<float>>(c);
        	for(int i=0;i<c;i++){
        		S[i]=new Vector<float>(c);
            	for(int j=0;j<c;j++){
            		Vector<float> v1=new Vector<float>(c);
            		Vector<float> v2=new Vector<float>(c);
                	for(int x=0;x<n;x++) v1[x]=(B[x][i]);
                	for(int x=0;x<n;x++) v2[x]=(B[x][j]);
                	S[i][j]=find_coefficient(v1,v2);
            	}
        	}


           // Determine each criterion weight wk calculated as per Eq. (18) and form the criteria weight vector Wa for agent ‘a’
        	Vector<float> w=new Vector<float>(c);
        	float total_sum=0;
        	for(int i=0;i<c;i++){
            	float sum=0;
            	for(int j=0;j<n;j++)
                	sum+=(1-S[i][j]);
            	w[i]=SD[i]*sum;        // Executing Eq (17)
	     	}
			for(int i=0;i<c;i++)
				total_sum+=w[i];

        	for(int i=0;i<c;i++)
            	w[i]=(w[i]/total_sum);   //  Executing Eq (18)

            //Add Wa as next row in W.
   			weightMatrix.add(w);     	        
    	}
   	}

	private static float find_coefficient(Vector<float> X, Vector<float> Y, int n){
   		int sum_X = 0, sum_Y = 0, sum_XY = 0;
   		int squareSum_X = 0, squareSum_Y = 0;
   		for (int i = 0; i < n; i++){
      		sum_X = sum_X + X[i];
      		sum_Y = sum_Y + Y[i];
      		sum_XY = sum_XY + X[i] * Y[i];
      		squareSum_X = squareSum_X + X[i] * X[i];
      		squareSum_Y = squareSum_Y + Y[i] * Y[i];
   		}
   		float corr = (float)(n * sum_XY - sum_X * sum_Y) / sqrt((n * squareSum_X - sum_X * sum_X) * (n * squareSum_Y - sum_Y * sum_Y));
   		return corr;
	}

	private static float calculateSD(Vector<float> data) {
  		float sum = 0.0, mean, standardDeviation = 0.0;
  		int i;

  		for(i = 0; i < 10; ++i) {
    		sum += data[i];
  		}

  		mean = sum / 10;

  		for(i = 0; i < 10; ++i) {
    		standardDeviation += pow(data[i] - mean, 2);
  		}

  		return sqrt(standardDeviation / 10);
	}


	private static void TOPSIS(){
    	int m=decisionMatrix.size();

    	for(int l=0;l<m;l++){
       		Vector<Vector<float>> B=decisionMatrix[l]; //taking each agent 'a' belongs to A=T U F
        	int n=B.size();
        	int c=B[0].size();

        	//Normalize the decision matrix using Eq. (19).
        	for(int j=0;j<c;j++){
            	float square_sum=0;
            	for(int i=0;i<n;i++) square_sum+=B[i][j]*B[i][j];
        	    square_sum=sqrt(square_sum);
       	 		for(int i=0;i<n;i++) B[i][j]=(B[i][j]/square_sum); // Eq. (19).
        	}

        	// Calculate the weighted normalized decision matrix based on Eq. (20)
        	for(int j=0;j<c;j++){
            	for(int i=0;i<n;i++) B[i][j]=B[i][j]*weightMatrix[l][j];  // Eq. (20)
        	}

        	Vector<float> positive=new Vector<float>(c);
        	Vector<float> negetive=new Vector<float>(c);
        	// Evaluate the positive and negative ideal solutions for each criterion k as 
        	// minimum and maximum value of kth column weighted normalized matrix.
        	for(int j=0;j<c;j++){
            	negetive[j]=positive[j]=B[0][j];
            	for(int i=0;i<n;i++){
                	negetive[j]=max(negetive[j],B[i][j]);
                	positive[j]=min(positive[j],B[i][j]);
            	}
        	}

        	Vector<float> d_pos=new Vector<float>(n);
        	Vector<float> d_neg=new Vector<float>(n);
        	Vector<float> p=new Vector<float>(n);
 			// Determine the distance of each alternative from the positive and negative ideal 
 			// solutions based on Eq.(21) and (22).
        	for(int i=0;i<n;i++){
            	for(int j=0;j<c;j++){
                	d_pos[i]+=(B[i][j]-positive[j])*(B[i][j]-positive[j]); 
                	d_neg[i]+=(B[i][j]-negtive[j])*(B[i][j]-negtive[j]);   
            	}
            	d_pos[i]=sqrt(d_pos[i]);  // Eq.(21)
            	d_neg[i]=sqrt(d_neg[i]);  // Eq.(22)
                                 
        	}

        	// Compute the performance score for each alternative following Eq. (23),
        	for(int i=0;i<n;i++){
        	 	p[i]=d_neg[i]/(d_neg[i]+d_pos[i]);     // Eq.(23)    
        	}

        	// ranked in decreasing order of their performance scores 
        	sort(p.begin(),p.end());
        	reverse(p.begin(),p.end());
        
        	performance.add(p);
    	}
    }



    public static void matching(){


    	// Q -- Quota of each FN
   		 Vector<Integer> Q=new Vector<Integer>(no_of_FN);

    	// Assign  -- Assigned tasks in FN
    	Vector<Vector<Integer>> Assign=new Vector<Vector<Integer>>(no_of_FN);

    	for(vector<float> tj:performance){
        	int n=tj.size();
        	for(int x=0;x<n;x++){
        		// fi∗= highest ranked FN in P(tj) to which tj has not proposed yet
        		// Send proposal to fi∗.
            	int fi=tj[x];
            	if(Q[fi]>0){    // if Qi∗ > 0 then
                	Assign[fi].push_back(tj);
                	Q[fi]=Q[fi]-1;
            	}
            	else{
                	// Reject the assignment request;
           		}
        	}
    	}


	}


	private static Application createApplication(String appId, int userId){
		
		Application application = Application.createApplication(appId, userId);
		/*
		 * Adding modules (vertices) to the application model (directed graph)
		 */
		
		return application;
	}
}