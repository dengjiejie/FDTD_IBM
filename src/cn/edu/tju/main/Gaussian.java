package cn.edu.tju.main;

public class Gaussian {

	/**
	 * @param args
	 */
	public void getGP(){
		double T = 0.646 * 1.0 / 10000;//Takemoto F0=10KHz
		int nn = 1;
		double E = 2.71828183;
		while(true){//dt(注意前面修改了别忘了修改)
			double t=0.000004;//*0.04 //
			double tempValue = -Math.pow((t* nn - T) * 1.0 / (0.29 * T), 2);
			double v = Math.pow(E, tempValue);
			if(nn > 2){
				if(v <= ParDefinition.gp[0])
					break;
			}
		
			ParDefinition.gp[nn - 1] = v;													
			nn++;
		}
		ParDefinition.nflag = nn - 1;
//		System.out.println("nn-----"+ParDefinition.nflag);
			for(int i = 0;i < ParDefinition.nflag;i++){
//						ParDefinition.gp[i] = Math.round(ParDefinition.gp[i] * 10000) * 1.0 / 10000;
			System.out.println(i+"\t"+ParDefinition.gp[i]);
		}
//		System.out.println("------------------");
	
	}

}
