package dynPgFWA;




import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Random;

import dynPgFWA.Info;
import dynPgFWA.testfunc;

public class dynPgFWA {
	 public static void main(String[] args) throws Exception
	  {
		  System.out.println("\nBegin fireworks algorithm optimization demo\n");
		  System.out.println("Goal is to find solution to Ackley's function for 10 variables");
		  System.out.println("Function has known min value = 0.0 at (0, 0, 0, 0, 0, 0, 0, 0, 0, 0)");
		  int fun_num = 15;
		  while(fun_num < 16)
		  {
		      int dim = 100; // often number weights to find in ML training
		      int n = 5;    // 烟花数量
		      int time = 0;
		      int maxEpochs = 10000 * dim;
		      System.out.println("\nSetting Ackley's dimension to " + dim);
		      System.out.println("Setting number fireworks to " + n);
		      System.out.println("Setting maxEpochs to " + maxEpochs);
	
		      System.out.println("\nBegin algorithm\n");
		      while(time < 30)
		      {
		    	  System.out.println(time);
		          double[] bestPosition = Solve(dim, n, maxEpochs, fun_num); // use fireworks algorithm
		          System.out.println("\nAlgorithm complete");
	
		          System.out.println("\nBest solution found: ");
		          for (int i = 0; i < dim; ++i)
		    	      System.out.print(bestPosition[i] + "   ");
		          System.out.println(" ");
	
		          testfunc tf3 = new testfunc();
		          double[] f3 = new double[1];
		          tf3.test_func(bestPosition, f3, dim, 1, fun_num);
		          double error = f3[0];
		          System.out.println("\nError of best solution found = " + error );
			      FileOutputStream fos = new FileOutputStream("e://改进FWA实验//data//dynPgFWA_fun"+fun_num+"_d"+dim+".txt",true);
			      OutputStreamWriter osw = new OutputStreamWriter(fos);
			      PrintWriter pw = new PrintWriter(osw);
			      pw.print(error + "\r\n");
			      pw.close();
		          ++time;
		      }
	
		      System.out.println("\nEnd fireworks algorithm optimization demo\n");
		      fun_num++;
		  }
	  }

	public static double[] Solve(int dim, int n, int maxEpochs,int fun_num) throws Exception
	  {
		int m = 150;  // total number of regular sparks for all fireworks
	    int mHat = 5;    // number of Gaussian sparks
	    double a = 0.04; // controls min sparks per firework
	    double b = 0.8;  // controls max sparks per firework
	    int A = 40;      // max amplitude
	    double minX = -100.0;
	    double maxX = 100.0;
	    int t = 0;
	    int number = 0;
	    testfunc tf0 = new testfunc();
	    double[] dimrange = new double[dim];
	    double[] pro = new double[dim];			//存放CF的爆炸半径
	    double[] delet =  new double[dim];
	    //double[] expro =  new double[dim];
	    int[] changed = new int[dim];
	    boolean better = false;//判断CF是否好于上一代的标识
	    double[] paraB = {1.3, 2.0, 1.3, 2, 1.5, 1.5, 1.5, 1.4, 1.5, 1.0, 1.2, 1.5, 1.5, 1.4, 1.6};
	    double[] paraLambda = {50,50,50,70,90,70,30,70,90,50,70,90,70,90,70};
	    double[] paraCR = {0.8,0.8,0.8,0.5,0.5,0.8,0.3,0.8,0.8,0.8,0.8,0.8,0.8,0.3,0.8};
	    double[] paraMinAP = {0.13,0.04,0.04,0.01,0.01,0.07,0.04,0.10,0.01,0.13,0.01,0.13,0.01,0.08,0.07};
	    
	    for(int i = 0; i < dim; i++)
	    {
	    	changed[i] = 0;
	    }
	    
	    for(int i = 0; i < dim; ++i)
	    {
	    	pro[i] = 1.0;
	    	delet[i] = 1;
	    	//expro[i] = 0;
	    }
	    
		int seed = (int) (Math.random()*100);    
	    Random rnd = new Random(seed); // seed = 3 gives representative demo
	    Info[] fireworks = new Info[n];//Info 存放位置和error的类   
	    for(int i = 0; i < n; i++)
	    {
	      fireworks[i] = new Info();
	      fireworks[i].position =  new double[dim];
	      for(int j = 0; j < dim; j++)
	    	  fireworks[i].position[j] = (maxX - minX) * rnd.nextDouble() + minX;
	    }
	    double[] x0 = new double[n * dim];
	    double[] f0 = new double[n];
	    for(int i = 0; i < n; ++i)
	    {
	    	for(int j = 0; j < dim; ++j)
	    	{
	    		x0[number] = fireworks[i].position[j];
	    		++number;
	    	}
	    }
	    tf0.test_func(x0, f0, dim, n, fun_num);
	    t += n;
	    for(int i = 0; i < n; ++i)
	    	fireworks[i].error = f0[i];
	    double[] bestPosition = new double[dim];
	    double besterr = fireworks[0].error;
	    int num = 0;
	    for(int i = 1 ; i < n; ++i)
	    {
	    	if(besterr > fireworks[i].error)
	    	{
	    		besterr = fireworks[i].error;
	    		num = i;
	    	}
	    }
	    for (int k = 0; k < dim; ++k)
	      bestPosition[k] = fireworks[num].position[k];
	    double bestError = fireworks[num].error;
	    
		ArrayList<ArrayList<Info>> sparksList = new ArrayList<ArrayList<Info>>();
		sparksList.clear();
		
		int epoch = 0;
		double Ca = 1;
		double range = 1;
	    while (t < maxEpochs)
	    {
	    	boolean isdyn = false;
	    	/*if (epoch % 10 == 0)
	    	{
		        FileOutputStream fos = new FileOutputStream("h://data4//dynCrFWA_fun"+fun_num+"_d"+dim+".txt",true);
		        OutputStreamWriter osw = new OutputStreamWriter(fos);
		        PrintWriter pw = new PrintWriter(osw);
		        pw.print(bestError + "\r\n");
		        pw.close();
	    	}*/
	    	
	    	int[] numberSparks = NumberSparks(fireworks, m, a, b);
	    	double[] amplitudes = Amplitudes(fireworks, A, t, maxEpochs, minX, maxX);
	    	
	    	
	        sparksList.clear();
	        if(epoch == 0)
	        	range = amplitudes[num];
	    	
	    	for (int i = 0; i < n; ++i)
	    	{
	    		double amp = amplitudes[i];
	    		int ns = numberSparks[i];
	    		Info[] spark = new Info[ns];
	    		for (int ks = 0; ks < ns; ++ks)
	    			spark[ks] = new Info();
	    		
	    		ArrayList<Info> temp = new ArrayList<Info>();
	    		temp.clear();
	    		if(epoch != 0 && i == 0){
	    			for(int k = 0; k < dim; ++k)
	    			{
	    				dimrange[k] = pro[k];
	    			}
	    		}
	    		for (int j = 0; j < ns; ++j)
	    		{
	    			spark[j].position =  new double[dim];
	    			for (int k = 0; k < dim; ++k)
	    				spark[j].position[k] = fireworks[i].position[k];
	    			//int z = (int)Math.round(dim * rnd.nextDouble());
	    			//int[] dimensions = PickDimensions(dim, z, epoch);
	    			//int[] dimChange = PickDimensions(dim);
	    			int[] dimChange = PickDimensions(dim);
	    			int[] dimChangeCF;
	    			if(better)
	    			{	dimChangeCF = PickDimensions(dim, changed, paraCR[fun_num-1]); }
	    			else
	    			{	dimChangeCF = PickDimensions(dim); }
	    			for (int ii = 0; ii < dim; ++ii)
	    			{
	    				if(i != 0)
	    				{
	    					if(dimChange[ii] == 1){
	    						double h = amp * (2 * rnd.nextDouble() - 1); // displacement hi = +1, lo = -1, (hi - lo) * r + lo
	    						spark[j].position[ii] += h; // displace from parent firework
	    						if (spark[j].position[ii] < minX || spark[j].position[ii] > maxX) // bring out-of-range values back in
	    						{
	    							spark[j].position[ii] = (maxX - minX) * rnd.nextDouble() + minX;
	    						}
	    					}
	    				}
	    				else
	    				{
	    					if(dimChangeCF[ii] == 1){
	    						double h = dimrange[ii] *( 2 * rnd.nextDouble() - 1); // displacement hi = +1, lo = -1, (hi - lo) * r + lo
			    	  
	    						
	    						//if((Math.random()*100) > 50)
	    						if((Math.random()*100) > paraLambda[fun_num - 1])
	    						{
	    							if(delet[ii] > 0)
	    							{
	    								if(h > 0)
	    									h = h * -1;
	    							}
	    							else if(delet[ii] < 0)
	    							{
	    								if(h < 0)
	    									h = h * -1;
	    							}
	    						}
	    						spark[j].position[ii] += h; // displace from parent firework
	    						if (spark[j].position[ii] < minX || spark[j].position[ii] > maxX) // bring out-of-range values back in
	    						{
	    							spark[j].position[ii] = (maxX - minX) * rnd.nextDouble() + minX;
	    						}
	    					}
	    				}
	    			}
	    		}
	    		double[] x1 = new double[ns * dim];
	    		double[] f1 = new double[ns];
	    		int number1 = 0;
	    		testfunc tf1 = new testfunc();
	    		for(int k = 0; k < ns; ++k)
	    		{
	    			for(int kk = 0; kk < dim; ++kk)
	    			{
	    				x1[number1] = spark[k].position[kk];
	    				++number1;
	    			}
	    		}
	    		tf1.test_func(x1, f1, dim, ns, fun_num);
	    		for(int k =0; k < ns; ++k)
	    		{
	    			spark[k].error = f1[k];
	    		}
	    		t += ns;
	    		for(int k = 0; k < ns; ++k)
	    		{
	    			 if (spark[k].error < bestError)
	    	          {
	    				 isdyn = true;
	    	            bestError = spark[k].error;
	    	            for (int kk = 0; kk < dim; ++kk)
		    	              bestPosition[kk] = spark[k].position[kk];
	    	          }
	    		}
	    		for(int k = 0; k < ns; ++k)
	    			temp.add(spark[k]);
	    		sparksList.add(temp);
	    	}
	    	AddGaussianSparks(t, fun_num, fireworks, sparksList, dim, mHat, epoch, minX, maxX, bestPosition, bestError, rnd);
	    	
	    	double[] bestSparkPos = new double[dim];
	        double bestSparkErr = 1.79E+308; 

	        double[] worstSparkPos = new double[dim];
	        double worstSparkErr = -1.79E+308;
	        
	        for (int i = 0; i < n; ++i) // numFireworks
	        {
	          for (int j = 0; j < sparksList.get(i).size(); ++j) // number sparks in each firework
	          {
	            if (sparksList.get(i).get(j).error < bestSparkErr)
	            {
	              bestSparkErr = sparksList.get(i).get(j).error;
	              for (int k = 0; k < sparksList.get(i).get(j).position.length; ++k)
	                bestSparkPos[k] = sparksList.get(i).get(j).position[k];
	            }
	            if (sparksList.get(i).get(j).error > worstSparkErr)
	            {
	              worstSparkErr = sparksList.get(i).get(j).error;
	              for (int k = 0; k < sparksList.get(i).get(j).position.length; ++k)
	                worstSparkPos[k] = sparksList.get(i).get(j).position[k];
	            }
	          } // each spark
	        } 
	        double[] exfireworkpos = new double[dim];
	        for (int k = 0; k < dim; ++k) // first new firework is best spark
	        {
	        	exfireworkpos[k] = fireworks[0].position[k];
	            fireworks[0].position[k] = bestSparkPos[k];
	        }
	          fireworks[0].error = bestSparkErr;
	          /*
	          for (int k = 0; k < dim; ++k) // second new firework is worst spark
	             fireworks[1].position[k] = worstSparkPos[k];
	          fireworks[1].error = worstSparkErr;
			*/
	          //selection part
	          
	          for (int i= 1; i < n; ++i) // n-1 random sparks 
	          {
	            int row = rnd.nextInt(n); // more likely to be a good spark
	            int cols = sparksList.get(row).size();
	            int col = rnd.nextInt(cols);
	            for (int k = 0; k < dim; ++k)
	              fireworks[i].position[k] = sparksList.get(row).get(col).position[k];
	            fireworks[i].error = sparksList.get(row).get(col).error;
	          }
	          
	          //new selection method
	          /*
	          int selectCount = 1;
	          int [] select = {0, 0, 0, 0, 0};
	          while(selectCount <5)
	          {
	        	  int row = rnd.nextInt(n);
	        	  if(select[row] >= 2)
	        	  {
	        		  continue;
	        	  }
	        	  else
	        	  {
	        		  int cols = sparksList.get(row).size();
	        		  int col = rnd.nextInt(cols);
	        		  for (int k = 0; k < dim; ++k)
	        		  {
	        			  fireworks[selectCount].position[k] = sparksList.get(row).get(col).position[k];
	        			  fireworks[selectCount].error = sparksList.get(row).get(col).error;
	        		  }
	        		  selectCount++;
	        		  select[row]++;
	        	  }
	          }*/
	          if(isdyn == true)
	        	  Ca = 1.2;
	          else
	        	  Ca = 0.9;
	          range = range * Ca;
	          int nozeronum = 0;
	          double total = 0;
	          for(int i = 0; i < dim; ++i)
	          {
	        	  delet[i] = fireworks[0].position[i] - exfireworkpos[i];
	        	  total += delet[i];
	        	  if(delet[i] != 0)
	        	  {
	        		  ++nozeronum;	  
	        	  	 changed[i] = 1;
    		  	  }
    		      else
    		      {
    			      changed[i] = 0;
    		      }
	          }
	          if(nozeronum > 1)
	          {
	        	  better = true;
	        	  for(int i = 0; i < dim; ++i)
	        	  {	pro[i] = paraB[fun_num - 1] * nozeronum * range * (delet[i] / total);}
	        	  //{	pro[i] = 2 * nozeronum * range * (delet[i] / total);}
	        	  
	        	  //double minAmpPara = 0.10;
	        	  
	        	  double minAmpPara = paraMinAP[fun_num-1];
	        	  double minAmp = (minAmpPara-minAmpPara*(double)epoch/(maxEpochs))*range;
	        	  
	        	  for(int i = 0; i < dim; i++)
	        	  {
	        		  
	        		  if(pro[i] < minAmp)
	        		  {pro[i] = minAmp; }
	        		 
	        	  }
	          }
	          else
	          {
	        	  better = false;
	        	  for(int i = 0; i < dim; ++i)
	        		  pro[i] = range;
	          }
	          epoch++;
	          // consider, instead selecting a row based on row-count . . 
	    }
	    return bestPosition;
	  }
	  
/*
	  private static int[] PickDimensions(int dim, int z, int seed)
	  {
	    // pick z random dimensions of a position
	    int[] result = new int[z];
	    int[] indices = new int[dim];
	    for (int i = 0; i < dim; ++i)
	      indices[i] = i;

	    Random rnd = new Random(seed); // shuffle indices
	    for (int i = 0; i < indices.length; ++i)
	    {
	      int ri = rnd.nextInt(indices.length - i) + i;
	      int tmp = indices[ri];
	      indices[ri] = indices[i];
	      indices[i] = tmp;
	    }
	    // copy first z indices to result
	    for (int i = 0; i < z; ++i)
	      result[i] = indices[i];

	    return result;
	  }*/
	
	private static int[] PickDimensions(int dim)
	{
	    // pick z random dimensions of a position
	    int[] result = new int[dim];
	    //Random rnd = new Random(seed);
	    for (int i = 0; i < dim; ++i)
	    {
	    	//double rndNumber = rnd.nextDouble();
	    	double rndNumber = Math.random();
	    	if(rndNumber > 0.5){
	    		result[i] = 0;
	    	}
	    	else
	    	{
	    		result[i] = 1;
	    	}
	    }
	    return result;
	  }
	
	private static int[] PickDimensions(int dim, int[] changed, double cr)
	{
		int[] result = new int[dim];
		for (int i = 0; i < dim; ++i)
	    {
	    	//double rndNumber = rnd.nextDouble();
	    	double rndNumber = Math.random();
	    	if(changed[i] == 1)
	    	{
	    		if(rndNumber > cr){
	    			result[i] = 1;
	    		}
	    		else
	    		{
	    			result[i] = 0;
	    		}
	    	}
	    	else
	    	{
	    		if(rndNumber > cr){
	    			result[i] = 0;
	    		}
	    		else
	    		{
	    			result[i] = 1;
	    		}
	    	}
	    }
		return result;
	}
	
	  private static double YMax(Info[] fireworks)
	  {
	    // largest (worst) error in any firework
	    double result = fireworks[0].error;
	    for (int i = 1; i < fireworks.length; ++i)
	      if (fireworks[i].error > result)
	        result = fireworks[i].error;
	    return result;
	  }

	  private static double YMin(Info[] fireworks)
	  {
	    // smallest (best) error in any firework
	    double result = fireworks[0].error;
	    for (int i = 1; i < fireworks.length; ++i)
	      if (fireworks[i].error < result)
	        result = fireworks[i].error;
	    return result;
	  }
	  
	  private static int[] NumberSparks(Info[] fireworks, int m, double a, double b)
	  {
	    // number sparks for each firework
	    int n = fireworks.length;
	    int minSparks = (int)Math.round(a * m); // if n=5, m=50, a=.04, -> 2
	    if (minSparks < 1) minSparks = 1;
	    int maxSparks = (int)Math.round(b * m); // if n=5, m=50, b=.8 -> 40
	    if (maxSparks > m - (n - 1) * minSparks)
	      maxSparks = m - (n - 1) * minSparks;

	    double yMax = YMax(fireworks);
	    double sumDeltas = 0.0; // sum diffs between yMax and each error
	    for (int i = 0; i < n; ++i)
	      sumDeltas += yMax - fireworks[i].error;

	    int[] numSparks = new int[n]; // the result
	    for (int i = 0; i < n; ++i)
	    {
	      numSparks[i] = (int)Math.round(m * (yMax - fireworks[i].error + 1.0E-10) / (sumDeltas + 1.0E-10));
	      if (numSparks[i] < minSparks)
	        numSparks[i] = minSparks;
	      else if (numSparks[i] > maxSparks)
	        numSparks[i] = maxSparks;
	    }
	    return numSparks;
	  }

	  private static double[] Amplitudes(Info[] fireworks, int A, int epoch, int maxEpochs, double minX, double maxX)
	  {
	    int n = fireworks.length;
	    double yMin = YMin(fireworks);
	    double sumDeltas = 0.0; // sum  diffs between yMin and each error
	    for (int i = 0; i < n; ++i)
	      sumDeltas += fireworks[i].error - yMin;

	    double[] result = new double[n]; // an amplitude for each firework
	    double minAmplitude = MinAmplitude(epoch, maxEpochs, minX, maxX);
	    for (int i = 0; i < n; ++i)
	    {
	      result[i] = A * (fireworks[i].error - yMin + 1.0E-10) / (sumDeltas + 1.0E-10);
	      if (result[i] < minAmplitude)
	        result[i] = minAmplitude;
	      if(i == 0)
	    	  result[i] = result[i]; 
	    }
	    return result;
	  }

	  private static double MinAmplitude(int epoch, int maxEpochs, double minX, double maxX)
	  {
	    // minimum amplitude for any firework at curr epoch 
	    double Ainit = (maxX - minX) * 0.02;
	    double Afinal = (maxX - minX) * 0.001;
	    return Ainit - (Ainit - Afinal) / maxEpochs * (Math.sqrt((2 * maxEpochs - epoch) * epoch));
	  }
	  /*
	  private static void AddGaussianSparks(int t, int fun_num, Info[] fireworks, ArrayList<ArrayList<Info>> sparksList, int dim, int mHat, int epoch, double minX, double maxX, double[] bestPosition,  double bestError, Random rnd)throws Exception
	  {
		  int n = fireworks.length;
		  for (int i = 0; i < mHat; i++)
			{
				double randflag = 0;						//随机数标识
				int selectGaussianID = rnd.nextInt(n);
				Info gSpark = new Info();
				gSpark.position = new double[dim];
				for (int j = 0; j < dim; j++)
				{
					//产生0到10的随机数，再除以10得到0到1的精确到十分位小数
					randflag = Math.random();
					if (randflag < 0.5)
					{
						int p = rnd.nextInt(n);
						int q = rnd.nextInt(n);
						int better = 0;
						//better =
						//	benchmarkEvaluation(newpopulation[p]) < benchmarkEvaluation(newpopulation[q]) ?
						//p : q;
						testfunc tf1 = new testfunc();
						double[] p_value = new double[1];
						double[] q_value = new double[1];
						tf1.test_func(fireworks[p].position, p_value, dim, 1, fun_num);
						tf1.test_func(fireworks[q].position, q_value, dim, 1, fun_num);
						if(p_value[0] < q_value[0])
						{
							better = p;
						}
						else
						{
							better = q;
						}
						double u1 = rnd.nextDouble(); // make Gaussian displacement u = 0, sd = 1
					    double u2 = rnd.nextDouble();
						double left = Math.cos(2.0 * Math.PI * u1);
					    double right = Math.sqrt(-2.0 * Math.log(u2));
					    double e = left * right;
						gSpark.position[j] = fireworks[selectGaussianID].position[j] + (fireworks[better].position[j] - fireworks[selectGaussianID].position[j])*e;
						if (gSpark.position[j] < minX || gSpark.position[j] > maxX) // bring out-of-range values back in
			        		gSpark.position[j] = (maxX - minX) * Math.random() + minX;
						
						//newpopulation[pointer][j] = newpopulation[selectGaussianID][j]
						//	+ (newpopulation[better][j] - newpopulation[selectGaussianID][j])
						//	*randGaussian(0, 1);
					}
				}
				testfunc tf2 = new testfunc();
				double[] x2 = new double[dim];
				double[] f2 = new double[1];
				for(int iii = 0; iii < dim; ++iii)
				{	x2[iii] = gSpark.position[iii];}
				tf2.test_func(x2, f2, dim, 1, fun_num);
				gSpark.error = f2[0];
				sparksList.get(i).add(gSpark);
				t++;
				if (gSpark.error < bestError)
			      {
			        bestError = gSpark.error;
			        for (int k = 0; k < dim; ++k)
			          bestPosition[k] = gSpark.position[k];
			      }

			}

	  }*/
	  
	  private static void AddGaussianSparks(int t, int fun_num, Info[] fireworks, ArrayList<ArrayList<Info>> sparksList, int dim, int mHat, int epoch, double minX, double maxX, double[] bestPosition,  double bestError, Random rnd) throws Exception
	  {
	    // generate mHat Gaussian sparks, add to sparksList, update bestPOsition[], bestError
	    int n = fireworks.length;
	    for (int g = 0; g < mHat; ++g)
	    {
	      Info gSpark = new Info();
	      gSpark.position = new double[dim];

	      int i = rnd.nextInt(n); // pick a random firework
	      
	      //System.out.println(i);
	      
	      for (int k = 0; k < dim; ++k) // spark position based on its parent firework
	        gSpark.position[k] = fireworks[i].position[k];
	      //int z = (int)Math.round(dim * rnd.nextDouble()); // number of random dimensions
	      int[] dimChange = PickDimensions(dim); // pick random dimensions
	      double u1 = rnd.nextDouble(); // make Gaussian displacement u = 0, sd = 1
	      double u2 = rnd.nextDouble();
	      double left = Math.cos(2.0 * Math.PI * u1);
	      double right = Math.sqrt(-2.0 * Math.log(u2));
	      double e = left * right; // mean = 0, sd = 1
	      for (int ii = 0; ii < dim; ++ii) // each of the randomly selected dimensions
	      {
	    	  if(dimChange[ii] == 1)
		        {
		        	gSpark.position[ii] = gSpark.position[ii] + (bestPosition[ii] - gSpark.position[ii]) * e;
		        	if (gSpark.position[ii] < minX || gSpark.position[ii] > maxX) // bring out-of-range values back in
		        		gSpark.position[ii] = (maxX - minX) * Math.random() + minX;
		        }
	      }
	      testfunc tf2 = new testfunc();
	      double[] x2 = new double[dim];
	      double[] f2 = new double[1];
	      for(int iii = 0; iii < dim; ++iii)
	    	  x2[iii] = gSpark.position[iii];
	      tf2.test_func(x2, f2, dim, 1, fun_num);
	      gSpark.error = f2[0];
	      ++t;

	      sparksList.get(i).add(gSpark);

	      if (gSpark.error < bestError)
	      {
	        bestError = gSpark.error;
	        for (int k = 0; k < dim; ++k)
	          bestPosition[k] = gSpark.position[k];
	      }
	    } // each Gaussian spark
	  }
}

class Info
{
  public double[] position;
  public double error;
}





