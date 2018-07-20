package cn.edu.tju.main;


public class SetCondition {
	/**
	 * 设置初始条件和边界条件
	 */
	static double E = 2.71828183;
	public void initCondition(){
		for(int i = 0;i < ParDefinition.im;i ++){
			for(int j = 0;j < ParDefinition.jm; j ++){
				ParDefinition.u[i][j] = 0;
				ParDefinition.v[i][j] = 0;
				ParDefinition.u0[i][j] = 0;
				ParDefinition.v0[i][j] = 0;
				ParDefinition.fx[i][j] = 0;
				ParDefinition.fy[i][j] = 0;
				ParDefinition.px[i][j] = 0;
				ParDefinition.py[i][j] = 0;
				
				 //点高斯+Takemoto波源
				/*
				double c_air = 346.3;//在空气中的声速
				double rho_air = 1.17;//空气的密度
				double tpluse=rho_air*Math.pow(c_air, 2);
				double Q = ParDefinition.p[i][j];//* ParDefinition.dt * 1.0 / (ParDefinition.dx * ParDefinition.dy)
				ParDefinition.px[i][j] = 1.0 * Q / 2;
				ParDefinition.py[i][j] = 1.0 * Q / 2;
				*/
				
				//p值按照论文中所给，计算出来(自己设置的波源，见论文)
				double temp1 = (Math.pow((ParDefinition.x[i][j] - 0.005), 2) + Math.pow(ParDefinition.y[i][j], 2)) * 1.0 / Math.pow(0.002, 2);
				double temp2 = -Math.log(2) * temp1;
				ParDefinition.p[i][j] = Math.pow(E,temp2)*10;
				ParDefinition.px[i][j] = ParDefinition.p[i][j]/2.0;
				ParDefinition.py[i][j] = ParDefinition.p[i][j]/2.0;	
				
			}
		}		
	}
	
	public void setPeriodicalSource(int minX, int minY){
		int origin=1;// Use Futang's setting
		if (origin == 1){
			for(int i = 0;i < ParDefinition.im;i ++){
				for(int j = 0;j < ParDefinition.jm; j ++){					
					//p值按照论文中所给，计算出来(福棠论文设置的波源，见论文)
					double temp1 = (Math.pow((ParDefinition.x[i][j] - 0.005), 2) + Math.pow(ParDefinition.y[i][j], 2)) / Math.pow(0.002, 2);
					double temp2 = -Math.log(2) * temp1;
					ParDefinition.p0[i][j] += Math.pow(E,temp2)*10;
					ParDefinition.px0[i][j] += ParDefinition.p0[i][j]/2.0;
					ParDefinition.py0[i][j] += ParDefinition.p0[i][j]/2.0;		
				}	
			}
		} else {
			int Band=8; 
			for(int i = minX;i < minX+Band; i ++){ //sound source appears at the glottal part
				for(int j = minY-Band;j < minY+Band; j ++){
					double temp1 = (Math.pow(Math.abs(ParDefinition.x[i][j] - 0.005), 1.6) + Math.pow(Math.abs(ParDefinition.y[i][j]),2)) / Math.pow(0.002, 2);
					//double temp1 = (Math.pow(Math.abs(ParDefinition.x[i][j] - 0.005), 2)) / Math.pow(0.002, 2);
					double temp2 = -Math.log(2) * temp1;
					ParDefinition.p0[i][j] += Math.pow(E,temp2)*10;
					ParDefinition.px0[i][j] += ParDefinition.p0[i][j]/2.0;
					ParDefinition.py0[i][j] += ParDefinition.p0[i][j]/2.0;		
				}	
			}
		}
	}
	
	public void setPeriodicalSourceTest(int minX, int minY, double glottalGAF, int Adding){
		double temp2 = -Math.log(2)/ Math.pow(0.002, 2);
		int Band=7,Band2=5; 
		for(int i = minX;i < minX+Band; i ++){ //sound source appears at the glottal part
			for(int j = minY-Band2;j < minY+Band2; j ++){						

				double temp1 = temp2 * (Math.pow((ParDefinition.x[i][j] - 0.002), 2) + Math.pow(ParDefinition.y[i][j], 4));				
				if (Adding==1) {
					ParDefinition.p0[i][j] = ParDefinition.p0[i][j]/2 +Math.pow(E,temp1)*glottalGAF;
					ParDefinition.px0[i][j]= (ParDefinition.px0[i][j] +ParDefinition.p0[i][j])/2.0;
					ParDefinition.py0[i][j]= (ParDefinition.py0[i][j] +ParDefinition.p0[i][j])/2.0;	
				} else {
					ParDefinition.p0[i][j] = Math.pow(E,temp1)*glottalGAF;
					ParDefinition.px0[i][j] = ParDefinition.p0[i][j]/2.0;
					ParDefinition.py0[i][j] = ParDefinition.p0[i][j]/2.0;	
				}
				//if (i==minX && j==minY) System.out.println("temp1(x, y)= " + ParDefinition.x[i][j] +" "+ ParDefinition.y[i][j] +" "+ temp1 + "; P = " + ParDefinition.p0[i][j] );
			}
		}
	}

	public void boundCondition(){
		//设置左边界条件
		for(int j = 0;j < ParDefinition.jm;j ++)
			ParDefinition.u[0][j] = 0;
		
		//设置上下边界条件
		for(int i = 0;i < ParDefinition.im;i ++){
			ParDefinition.u[i][0] = 0;
			ParDefinition.u[i][ParDefinition.jm - 1] = 0;
		}
		
	}
}
