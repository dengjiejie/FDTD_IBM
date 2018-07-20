package cn.edu.tju.main;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class TimeOperation {
	
	double[][] s = new double[ParDefinition.im][ParDefinition.jm];
	static double E = 2.71828183;
	//--------------------------------- ������*/
	/*
	int npml = 8;
	double c_air = 1, rho_air = 1;
	double alpha = 0;
//	double flow_resistance = alpha * Math.pow(rho_air, 2) * Math.pow(c_air, 2); 
	double kappa = 1.0 / (rho_air * Math.pow(c_air, 2));
	int m = npml;
	int magicnum = 1;
	double alpha0 = Math.log(10) * kappa * c_air * 1.0 / ParDefinition.dx * magicnum;	
	double z = 0;
	*/
	
	//------------------------------------------------------------
	//--------------------------------- ������
	int npml = 8;
	double temperature = 25;
	double c_air = 331.45 * (1 + 0.0018 * temperature);//�ڿ����е�����
	double rho_air = 1.2929 * (1 - 0.0037 * temperature);//�������ܶ�
	//double c_air = 1;
	//double rho_air = 1;
	double alpha = Math.pow(10, -6); //�ڼ�����������Ϊ10�����η� һ������Ϊ0  Math.pow(10, -6)
	double alpha_wall=100;//wall�ϵ�������
	//alpha_air = 1.0e-6; %1.0e-8;
//	double flow_resistance = alpha * Math.pow(rho_air, 2) * Math.pow(c_air, 2); 
	double kappa = 1.0 / (rho_air * Math.pow(c_air, 2));
	int m = npml;
	int magicnum = 5;
	//PML_alpha_max
	double alpha0 = -Math.log(0.1) * kappa * c_air * 1.0 / ParDefinition.dx * magicnum;	
	//System.out.println("alpha0 ="+alpha0);
	
	//------------------------------------------------------------
	double[][] alpha_x = new double[ParDefinition.im][ParDefinition.jm];
	double[][] alpha_y = new double[ParDefinition.im][ParDefinition.jm];
	double[][] alpha_s_x = new double[ParDefinition.im][ParDefinition.jm];
	double[][] alpha_s_y = new double[ParDefinition.im][ParDefinition.jm];
	//û���õ�
	double[][] a_x = new double[ParDefinition.im][ParDefinition.jm];
	double[][] a_y = new double[ParDefinition.im][ParDefinition.jm];
	
	//���㸽����
	public void virForce(){
		
		int ii1,ii2,jj1,jj2;
		double xx1 = 100,xx2 = -100,yy1 = 100,yy2 =-100;

		//xx1�д�������������յ������������꣬xx2�д�������������յ������ҵ������
		//yy1�д�������������յ������ϵ�����꣬yy2�д�������������յ������µ�����꣬
		for(int i = 0;i < ParDefinition.nlp;i ++){
			if (ParDefinition.xl[i] < xx1)	xx1 = ParDefinition.xl[i];
			if (ParDefinition.xl[i] > xx2)	xx2 = ParDefinition.xl[i];
			if (ParDefinition.yl[i] < yy1)	yy1 = ParDefinition.yl[i];
			if (ParDefinition.yl[i] > yy2)	yy2 = ParDefinition.yl[i];
		}      
//		System.out.println(xx1+"  "+xx2+"  "+yy1+"  "+yy2+"  " + ParDefinition.dxy);
		//�ҳ�Բ�����С���α߽磬֮�󸽼����ļ���ֻ��Ҫ�ڸ����������⼴��
		ii1=0 + (int)((xx1 - ParDefinition.x[0][0]) / ParDefinition.dx) - 2;
		ii2=0 + (int)((xx2 - ParDefinition.x[0][0]) / ParDefinition.dx) + 3;
		jj1=0 + (int)((yy1 - ParDefinition.y[0][0]) / ParDefinition.dy) - 2;
		jj2=0 + (int)((yy2 - ParDefinition.y[0][0]) / ParDefinition.dy) + 3;
//		System.out.println(ii1+"  "+ii2+"  "+jj1+"  "+jj2+"  " + ParDefinition.dxy);
		
		for(int j = jj1;j <= jj2;j ++){
			for(int i = ii1;i <= ii2;i ++){
				double rmin = 100;
				double rdist;
				double vecx,vecy,recx,recy,temp;
				int iflag = 0;
				int nis,njs;
				double etx,ety,etr;
				double uuc,vvc,ppl,ppr,ppb,ppt;
				double termpreu,fxlp,termprev,fylp;
				
				//�ҳ��������(i,j)������������յ㣬����iflag��¼
				for(int ii = 0;ii < ParDefinition.nlp;ii ++){
					rdist = Math.sqrt(Math.pow((ParDefinition.xl[ii]-ParDefinition.x[i][j]), 2)+ Math.pow((ParDefinition.yl[ii]-ParDefinition.y[i][j]), 2));
					if(rdist < rmin) {
						rmin = rdist;
						iflag = ii;	
					}		
				}
				
				//���жϸ�������Ƿ������ڡ����˼�룺�ȸ���lag��������ڵ�������ߣ�֮��ɵõ����淨������Ϊ: vecy - vecx i
				//���ݣ��������յ�lp��ŷ����������ɵ�����������������������Ƿ����0�����ж��Ƿ�λ��Բ�����ڡ�
				 if(iflag == 0){
					 vecx = ParDefinition.xl[1] - ParDefinition.xl[0];
					 vecy = ParDefinition.yl[1] - ParDefinition.yl[0];
				        
				 }else if(iflag == (ParDefinition.nlp - 1)){
					 vecx = ParDefinition.xl[ParDefinition.nlp - 1] - ParDefinition.xl[ParDefinition.nlp - 2];
					 vecy = ParDefinition.yl[ParDefinition.nlp - 1] - ParDefinition.yl[ParDefinition.nlp - 2];
				        
				 }else{
					 vecx = ParDefinition.xl[iflag + 1] - ParDefinition.xl[iflag - 1];
					 vecy = ParDefinition.yl[iflag + 1] - ParDefinition.yl[iflag - 1];
				 }
	                    
	               
				recx = ParDefinition.x[i][j] - ParDefinition.xl[iflag];
				recy = ParDefinition.y[i][j] - ParDefinition.yl[iflag];
				temp =  vecy * recx - vecx * recy;
				
				
				//����������������ڣ���������֮��ľ���<dxy
				if((rmin <= ParDefinition.dxy) && (temp <= 0)){
					ParDefinition.just[i][j] = 1;
					//(nis,njs)��¼�������յ����������㣨�����㣩��λ��
					nis = (int) (0 + (ParDefinition.xl[iflag] - ParDefinition.x[0][0])/ParDefinition.dx);
					njs = (int) (0 + (ParDefinition.yl[iflag] - ParDefinition.y[0][0])/ParDefinition.dy);					
					
					//���������յ��ÿһ�����������Ϣ������������Χ�ĸ�������ֵ ˫���Բ�ֵ�õ���					
					//���pֵ�Ĳ�ֵϵ��
					etx = (ParDefinition.xl[iflag] - ParDefinition.x[nis][njs]) * 1.0 / ParDefinition.dx;
					ety = (ParDefinition.yl[iflag] - ParDefinition.y[nis][njs]) * 1.0 / ParDefinition.dy;
					
					//���ĵ㣨���������յ㣩u,vֵ���㡣
					uuc = (1.0 - etx) * (1.0 - ety) * ParDefinition.u0[nis][njs];
	                uuc = uuc + etx * (1.0 - ety) * ParDefinition.u0[nis + 1][njs];
	                uuc = uuc + (1.0 - etx) * ety * ParDefinition.u0[nis][njs + 1];
	                uuc = uuc + etx * ety * ParDefinition.u0[nis + 1][njs + 1];
	                
	                vvc = (1.0 - etx) * (1.0 - ety) * ParDefinition.v0[nis][njs];
	                vvc = vvc + etx * (1.0 - ety) * ParDefinition.v0[nis + 1][njs];
	                vvc = vvc + (1.0 - etx) * ety * ParDefinition.v0[nis][njs + 1];
	                vvc = vvc + etx * ety * ParDefinition.v0[nis + 1][njs + 1];							
					
					//�������pֵ
					ppl = (1.0 - etx) * (1.0 - ety) * ParDefinition.p0[nis - 1][njs];
					ppl = ppl + etx * (1.0 - ety) * ParDefinition.p0[nis][njs];
					ppl = ppl + (1.0 - etx) * ety * ParDefinition.p0[nis - 1][njs + 1];
					ppl = ppl + etx * ety * ParDefinition.p0[nis][njs + 1];
					
					//�Ҹ������pֵ
					ppr = (1.0 - etx) * (1.0 - ety) * ParDefinition.p0[nis + 1][njs];
					ppr = ppr + etx * (1.0 - ety) * ParDefinition.p0[nis + 2][njs];
					ppr = ppr + (1.0 - etx) * ety * ParDefinition.p0[nis + 1][njs + 1];
					ppr = ppr + etx * ety * ParDefinition.p0[nis + 2][njs + 1];
					
					//�¸������pֵ
					ppb = (1.0 - etx) * (1.0 - ety) * ParDefinition.p0[nis][njs - 1];
					ppb = ppb + etx * (1.0 - ety) * ParDefinition.p0[nis + 1][njs - 1];
					ppb = ppb + (1.0 - etx) * ety * ParDefinition.p0[nis][njs];
					ppb = ppb + etx * ety * ParDefinition.p0[nis + 1][njs];
					
					//�ϸ������pֵ
					ppt = (1.0 - etx) * (1.0 - ety) * ParDefinition.p0[nis][njs + 1];
					ppt = ppt + etx * (1.0 - ety) * ParDefinition.p0[nis + 1][njs + 1];
					ppt = ppt + (1.0 - etx) * ety * ParDefinition.p0[nis][njs + 2];
					ppt = ppt + etx * ety * ParDefinition.p0[nis + 1][njs + 2];
	
					//�����������յ㴦��x�����ϵ�ѹ��������������
					termpreu = (ppr - ppl) * 0.5 / ParDefinition.dx;
					//fxlp= (0 - uuc) * 1.0 / ParDefinition.dt + termpreu;
					fxlp= (0 - uuc) * 1.0 / ParDefinition.dt*rho_air+ termpreu;//tony_wang 2017.01.13
					
					//�����������յ㴦��y�����ϵ�ѹ��������������
	                termprev = (ppt - ppb) * 0.5 / ParDefinition.dy;
	                //fylp = (0 - vvc) * 1.0 / ParDefinition.dt + termprev;
	                fylp = (0 - vvc) * 1.0 / ParDefinition.dt*rho_air + termprev;//tony_wang 2017.01.13
	                
	                //��ֵ�õ���������ϵĸ�������ѹ���ݶ�ֵ
					etr = rmin * 1.0 / ParDefinition.dxy;
	                ParDefinition.fx[i][j] = (1 - etr) * fxlp;
	                ParDefinition.fy[i][j] = (1 - etr) * fylp;
	                ParDefinition.fxp[i][j] = (1 - etr) * termpreu;
	                ParDefinition.fyp[i][j] = (1 - etr) * termprev;
	                ParDefinition.tube[i][j] = 255;
	               
	                
	                // ����s�����ϵ�PML
					// alpha_x[i][j] = alpha_wall;
					// alpha_y[i][j] = alpha_wall;
					// alpha_s_x[i][j] = alpha;
					// alpha_s_y[i][j] = alpha;
	            
				}
				//���������������ⲿ��ֱ�����������ϵĸ��������ɡ�
				else if(temp >0&&rmin <=ParDefinition.dxy){//&&rmin <=2*ParDefinition.dx
					ParDefinition.just[i][j] = 2;
					//���Ĳ���󸽼���
					termpreu = (ParDefinition.p0[i + 1][j] - ParDefinition.p0[i - 1][j]) * 0.5 / ParDefinition.dx;
					ParDefinition.fx[i][j] = (0 - ParDefinition.u0[i][j]) * 1.0 / ParDefinition.dt*rho_air+ termpreu;
					ParDefinition.fxp[i][j] = termpreu;
					
					termprev = (ParDefinition.p0[i][j + 1] - ParDefinition.p0[i][j - 1]) * 0.5 / ParDefinition.dy;
					ParDefinition.fy[i][j] = (0 - ParDefinition.v0[i][j]) * 1.0 / ParDefinition.dt*rho_air	+ termprev;
					ParDefinition.fyp[i][j] = termprev;
					
					// ���������ϵ�PML
					 //alpha_x[i][j] = alpha_wall;
					 //alpha_y[i][j] = alpha_wall;
					 //alpha_s_x[i][j] = alpha;
					 //alpha_s_y[i][j] = alpha;
				}else{
					ParDefinition.just[i][j] = 3;
					//ParDefinition.fy[i][j]=0;
					//ParDefinition.fx[i][j]=0;
					 
				}//if end
				
			}//i
		}//j
		
	}//virForce()

	public void setPML(){	
		//System.out.println("alpha0 ="+alpha0);
		//L,R
		//double alpha = Math.pow(10, -6); 
		for(int i = 1;i < m;i ++){
			for(int j = m;j <= (ParDefinition.jm - m - 2);j ++){//���j�����������⣿�ѵ�����jm-1-m
				alpha_x[i][j] = alpha0 * Math.pow((m - i) * 1.0 / m, 2);
				alpha_x[ParDefinition.im - i][j] = alpha_x[i][j];
				alpha_y[i][j] = alpha;
				alpha_y[ParDefinition.im - i][j] = alpha_y[i][j];
			}
		}
		for(int i = 0;i < m;i ++){
			for(int j = m;j <= (ParDefinition.jm - m - 2);j ++){
				alpha_s_x[i][j] = alpha0 * Math.pow((m - i - 0.5) * 1.0 / m, 2) * rho_air / kappa;
				alpha_s_x[ParDefinition.im - 1 - i][j] = alpha_s_x[i][j];
				alpha_s_y[i][j] =alpha;
				alpha_s_y[ParDefinition.im - 1 - i][j] = alpha_s_y[i][j];			
			}
		}
		//T,B
		for(int j = 1;j < m;j ++){
			for(int i = m;i <= (ParDefinition.im - m - 2);i ++){
				alpha_y[i][j] = alpha0 * Math.pow((m - j) * 1.0 / m, 2);
				alpha_y[i][ParDefinition.jm - j] = alpha_y[i][j];
				alpha_x[i][j] = alpha;
				alpha_x[i][ParDefinition.jm - j] = alpha_x[i][j];
			}
		}
		for(int j = 0;j < m;j ++){
			for(int i = m;i <= (ParDefinition.im - m - 2);i ++){
				alpha_s_y[i][j] = alpha0 * Math.pow((m - j - 0.5) * 1.0 / m, 2) * rho_air / kappa;
				alpha_s_y[i][ParDefinition.jm - 1 - j] = alpha_s_y[i][j];
				alpha_s_x[i][j] = alpha;
				alpha_s_x[i][ParDefinition.jm - 1 - j] = alpha_s_x[i][j];			
			}
		}
		//�ĸ�Corner
		//C1
		for(int i = 1;i <= m;i ++){
			for(int j = 1;j <= m;j ++){
				alpha_x[i][j] = alpha0 * Math.pow((m + 1 - i) * 1.0 / (m + 1), 2);
				alpha_y[i][j] = alpha0 * Math.pow((m + 1 - j) * 1.0 / (m + 1), 2);
			}
		}
		for(int i = 0;i <= m;i ++){
			for(int j = 0;j <= m;j ++){
				alpha_s_x[i][j] = alpha0 * Math.pow((m + 1 - i - 0.5) * 1.0 / (m + 1), 2) * rho_air / kappa;
				alpha_s_y[i][j] = alpha0 * Math.pow((m + 1 - j - 0.5) * 1.0 / (m + 1), 2) * rho_air / kappa;
			}
		}
		//C2
		for(int i = 1;i <= m;i ++){
			for(int j = 1;j <= m;j ++){
				alpha_x[i][ParDefinition.jm - 1 - j] = alpha0 * Math.pow((m + 1 - i) * 1.0 / (m + 1), 2);
				alpha_y[i][ParDefinition.jm - 1 - j] = alpha0 * Math.pow((m + 1 - j) * 1.0 / (m + 1), 2);
			}
		}
		for(int i = 0;i <= m;i ++){
			for(int j = 0;j <= m;j ++){
				alpha_s_x[i][ParDefinition.jm - 1 - j] = alpha0 * Math.pow((m + 1 - i - 0.5) * 1.0 / (m + 1), 2) * rho_air / kappa;
				alpha_s_y[i][ParDefinition.jm - 1 - j] = alpha0 * Math.pow((m + 1 - j - 0.5) * 1.0 / (m + 1), 2) * rho_air / kappa;
			}
		}
		//C4
		for(int i = 1;i <= m;i ++){
			for(int j = 1;j <= m;j ++){
				alpha_x[ParDefinition.im - 1 - i][j] = alpha0 * Math.pow((m + 1 - i) * 1.0 / (m + 1), 2);
				alpha_y[ParDefinition.im - 1 - i][j] = alpha0 * Math.pow((m + 1 - j) * 1.0 / (m + 1), 2);
			}
		}
		for(int i = 0;i <= m;i ++){
			for(int j = 0;j <= m;j ++){
				alpha_s_x[ParDefinition.im - 1 - i][j] = alpha0 * Math.pow((m + 1 - i - 0.5) * 1.0 / (m + 1), 2) * rho_air / kappa;
				alpha_s_y[ParDefinition.im - 1 - i][j] = alpha0 * Math.pow((m + 1 - j - 0.5) * 1.0 / (m + 1), 2) * rho_air / kappa;
			}
		}
		//C3
		for(int i = 1;i <= m;i ++){
			for(int j = 1;j <= m;j ++){
				alpha_x[ParDefinition.im - 1- i][ParDefinition.jm - j] = alpha0 * Math.pow((m + 1 - i) * 1.0 / (m + 1), 2);
				alpha_y[ParDefinition.im - 1- i][ParDefinition.jm - j] = alpha0 * Math.pow((m + 1 - j) * 1.0 / (m + 1), 2);
			}
		}
		for(int i = 0;i <= m;i ++){
			for(int j = 0;j <= m;j ++){
				alpha_s_x[ParDefinition.im - 1 - i][ParDefinition.jm - 1 - j] = alpha0 * Math.pow((m + 1 - i - 0.5) * 1.0 / (m + 1), 2) * rho_air / kappa;
				alpha_s_y[ParDefinition.im - 1 - i][ParDefinition.jm - 1 - j] = alpha0 * Math.pow((m + 1 - j - 0.5) * 1.0 / (m + 1), 2) * rho_air / kappa;
			}
		}
		
		//for Computational domain
		for(int i = m; i < (ParDefinition.im - m);i ++){
			for(int j = m;j < (ParDefinition.jm - m);j ++){    
				  alpha_x[i][j] = alpha;
				  alpha_y[i][j] = alpha;
				  alpha_s_x[i][j] = alpha;
				  alpha_s_y[i][j] = alpha;
			}
			
		}
		
	}
	

	//���P U V
	public void solvePUV(){

		double[][] epx1 = new double[ParDefinition.im][ParDefinition.jm];
		double[][] epx2 = new double[ParDefinition.im][ParDefinition.jm];
		double[][] epy1 = new double[ParDefinition.im][ParDefinition.jm];
		double[][] epy2 = new double[ParDefinition.im][ParDefinition.jm];
		double[][] eux1 = new double[ParDefinition.im][ParDefinition.jm];
		double[][] eux2 = new double[ParDefinition.im][ParDefinition.jm];
		double[][] euy1 = new double[ParDefinition.im][ParDefinition.jm];
		double[][] euy2 = new double[ParDefinition.im][ParDefinition.jm];
		//setPML();
//		setImpedance();
		for(int i = 0; i < ParDefinition.im;i ++){
			for(int j = 0;j < ParDefinition.jm;j ++){
				double C1 = kappa * 1.0 / ParDefinition.dt - alpha_x[i][j] * 1.0 / 2;
				double C2 = kappa * 1.0 / ParDefinition.dt + alpha_x[i][j] * 1.0 / 2;
				epx1[i][j] = C1 * 1.0 / C2;
				epx2[i][j] = 1.0 / (C2 * ParDefinition.dx);
		            
				C1 = kappa * 1.0 / ParDefinition.dt - alpha_y[i][j] * 1.0 / 2;
				C2 = kappa * 1.0 / ParDefinition.dt + alpha_y[i][j] * 1.0 / 2;
				epy1[i][j] = C1 * 1.0 / C2;
				epy2[i][j] = 1.0 / (C2 * ParDefinition.dy);
				
				//2017.01.10 Tonywang
		        //C1 = rho_air * 1.0 / ParDefinition.dt - alpha_s_x[i][j]/2;
				C1 = rho_air * 1.0 / ParDefinition.dt - alpha_s_x[i][j]* 1.0 / 2;
				//C2 = rho_air * 1.0 / ParDefinition.dt + alpha_s_x[i][j]/2;
				C2 = rho_air * 1.0 / ParDefinition.dt + alpha_s_x[i][j]* 1.0 / 2 ;
		    	eux1[i][j] = C1 * 1.0 / C2;
		    	eux2[i][j] = 1.0 / (C2 * ParDefinition.dx);
		            
		    	
		    	C1 = rho_air * 1.0 / ParDefinition.dt - alpha_s_y[i][j] * 1.0/ 2 ;
		    	C2 = rho_air * 1.0 / ParDefinition.dt + alpha_s_y[i][j] * 1.0/ 2;
		    	
		    	//C1 = rho_air * 1.0 / ParDefinition.dt - alpha_s_y[i][j] * rho_air * 1.0  / (2 * kappa);
		    	//C2 = rho_air * 1.0 / ParDefinition.dt + alpha_s_y[i][j] * rho_air * 1.0  / (2 * kappa);
		     	euy1[i][j] = C1 * 1.0 / C2;
		     	euy2[i][j] = 1.0 / (C2 * ParDefinition.dy);	           		   
				
			}
		}
		
		
		//��p��ע���ڼ���ʱ����δ���Ǵ�������ĸ��߽��ϵ�puvֵ����һֱ��Ϊ���ǲ��䣬Ϊ0��֮��Ӧ�ÿ��Ǳ߽�����
		for(int i = 1;i < (ParDefinition.im - 1);i ++){
			for(int j = 1;j < (ParDefinition.jm - 1);j ++){
				
				ParDefinition.px[i][j] = epx1[i][j] * ParDefinition.px0[i][j] - epx2[i][j] * (ParDefinition.u0[i][j] - ParDefinition.u0[i - 1][j]);
				ParDefinition.py[i][j] = epy1[i][j] * ParDefinition.py0[i][j] - epy2[i][j] * (ParDefinition.v0[i][j] - ParDefinition.v0[i][j - 1]);
			
				ParDefinition.p[i][j] = ParDefinition.px[i][j] + ParDefinition.py[i][j];
				
				//if(ParDefinition.p[i][j]<10000)
				   //System.out.println("��ǰ�����"+"("+i+","+j+")"+"ѹ��ֵ"+ParDefinition.p[i][j]);
				//else
				  // break;
			}	
		}
		
	  /*
		
		if(TestProject.itime <= ParDefinition.nflag){
//			System.out.println("��ǰʱ�䲽��" + TestProject.itime+"     "+ParDefinition.nflag);
			double tempPulse = rho_air * Math.pow(c_air, 2) * ParDefinition.gp[TestProject.itime - 1] 
					* ParDefinition.dt * 1.0 / (ParDefinition.dx * ParDefinition.dy);
			System.out.println("----------------------"+tempPulse);
			//double tempPulse = rho_air * Math.pow(c_air, 2) * ParDefinition.gp[TestProject.itime - 1];
			ParDefinition.px[ParDefinition.si][ParDefinition.sj] = ParDefinition.px[ParDefinition.si][ParDefinition.sj] + tempPulse * 1.0 / 2;
			ParDefinition.py[ParDefinition.si][ParDefinition.sj] = ParDefinition.py[ParDefinition.si][ParDefinition.sj] + tempPulse * 1.0 / 2;
			ParDefinition.p[ParDefinition.si][ParDefinition.sj] = ParDefinition.px[ParDefinition.si][ParDefinition.sj] + ParDefinition.py[ParDefinition.si][ParDefinition.sj];
//			System.out.println("----nflag-"+ParDefinition.nflag);
		}  
		*/
		//��u,v	
		for(int i = 1;i < (ParDefinition.im - 1);i ++){
			for(int j = 1;j < (ParDefinition.jm - 1);j ++){
				double cc1 = rho_air * 1.0 / ParDefinition.dt + alpha_s_x[i][j] * 1.0 / 2;
				ParDefinition.u[i][j] = eux1[i][j] * ParDefinition.u0[i][j] - eux2[i][j] * (ParDefinition.p[i + 1][j] - ParDefinition.p[i][j]) 
						+ ParDefinition.fx[i][j] * 1.0 / cc1;
				double cc2 = rho_air * 1.0 / ParDefinition.dt + alpha_s_y[i][j] * 1.0 / 2;
				ParDefinition.v[i][j] = euy1[i][j] * ParDefinition.v0[i][j] - euy2[i][j] * (ParDefinition.p[i][j + 1] - ParDefinition.p[i][j]) 
						+ ParDefinition.fy[i][j] * 1.0 / cc2;				  
			}
		}
		
		

	}//solvePUV()
	
	//�ú������ڴ��ÿ������������ֵ �Լ���Ӧ��pֵ
	public void saveXYP(int fileCount){
		String fileName = "src/cn/edu/tju/main/values/puv_" + fileCount + ".dat";
		String fileName1= "src/cn/edu/tju/main/values/VTB_" + fileCount + ".dat";
		try {
			File outFile = new File(fileName);
			// if file doesn't exists, then create it
			if (!outFile.exists()) {
				outFile.createNewFile();
			}
			
			FileWriter fw = new FileWriter(outFile.getAbsoluteFile());
			BufferedWriter bw = new BufferedWriter(fw);

			for(int i = 0;i < ParDefinition.im;i ++){
				for(int j = 0;j < ParDefinition.jm;j ++){
					bw.write(ParDefinition.x[i][j] + "\t" +ParDefinition.y[i][j] + "\t" + ParDefinition.p[i][j]);
		            bw.newLine();//����
				}		
			}
	               
			bw.close();
			fw.close();
			
			//Record the VT shape
			
			File outFile1 = new File(fileName1);
			FileWriter fw1 = new FileWriter(outFile1.getAbsoluteFile());
			BufferedWriter bw1 = new BufferedWriter(fw1);
			for (int i=0; i<ParDefinition.nlp; i++ ){ //The boundary of VTs
				bw1.write(ParDefinition.xl[i] + "\t" +ParDefinition.yl[i]);
				bw1.newLine();//���� 
			}
	        bw1.close();
	        fw1.close();

    	  } catch (IOException e) {
    	   e.printStackTrace();
    	  }

	}
	
	//����������ѹ��ֵ�����У�����ox��oy��ʾ��������ǵ�����������
		public void saveValues_p(int ox,int oy,String fileName){

			try {
				File outFile = new File(fileName);
				// if file doesnt exists, then create it
				if (!outFile.exists()) {
					outFile.createNewFile();
				}
		
				FileWriter fw = new FileWriter(outFile.getAbsoluteFile(),true);
				BufferedWriter bw = new BufferedWriter(fw);			
				bw.write(ParDefinition.p[ox][oy] + "");
			    bw.newLine();//����
	    
				bw.close();
				fw.close();

	    	  } catch (IOException e) {
	    	   e.printStackTrace();
	    	  }   	 
		}
		
		//save the presure filtered by low pass filter p(n)=0.9*p(n-1)+x(n)
		public void saveValues_pNew(double p,String fileName){

			try {
				File outFile = new File(fileName);
				// if file doesnt exists, then create it
				if (!outFile.exists()) {
					outFile.createNewFile();
				}
		
				FileWriter fw = new FileWriter(outFile.getAbsoluteFile(),true);
				BufferedWriter bw = new BufferedWriter(fw);			
				bw.write(p + "");
			    bw.newLine();//����
	    
				bw.close();
				fw.close();

	    	  } catch (IOException e) {
	    	   e.printStackTrace();
	    	  }   	 
		}
	
	public void saveValues_u(int ox,int oy,String fileName){

		try {
			File outFile = new File(fileName);
			// if file doesnt exists, then create it
			if (!outFile.exists()) {
				outFile.createNewFile();
			}
	
			FileWriter fw = new FileWriter(outFile.getAbsoluteFile(),true);
			BufferedWriter bw = new BufferedWriter(fw);			
			bw.write(ParDefinition.u[ox][oy] + "");
		    bw.newLine();//����
    
			bw.close();
			fw.close();

    	  } catch (IOException e) {
    	   e.printStackTrace();
    	  }   	 
	}
	
	public void saveValues_v(int ox,int oy,String fileName){

		try {
			File outFile = new File(fileName);
			// if file doesnt exists, then create it
			if (!outFile.exists()) {
				outFile.createNewFile();
			}
	
			FileWriter fw = new FileWriter(outFile.getAbsoluteFile(),true);
			BufferedWriter bw = new BufferedWriter(fw);			
			bw.write(ParDefinition.v[ox][oy] + "");
		    bw.newLine();//����
    
			bw.close();
			fw.close();

    	  } catch (IOException e) {
    	   e.printStackTrace();
    	  }   	 
	}
	
	public void saveF(int ox,int oy,String fileName){

		try {
			File outFile = new File(fileName);
			// if file doesnt exists, then create it
			if (!outFile.exists()) {
				outFile.createNewFile();
			}
	
			FileWriter fw = new FileWriter(outFile.getAbsoluteFile(),true);
			BufferedWriter bw = new BufferedWriter(fw);			
			bw.write(ParDefinition.fx[ox][oy] + " " + ParDefinition.fy[ox][oy]);
		    bw.newLine();//����
    
			bw.close();
			fw.close();

    	  } catch (IOException e) {
    	   e.printStackTrace();
    	  }   	 
	}
	
	public void saveFlag(int fileCount){
		String fileName = "src/cn/edu/tju/main/values/flag" + fileCount + ".dat";
		try {
			File outFile = new File(fileName);
			// if file doesnt exists, then create it
			if (!outFile.exists()) {
				outFile.createNewFile();
			}
			
			FileWriter fw = new FileWriter(outFile.getAbsoluteFile());
			BufferedWriter bw = new BufferedWriter(fw);

			for(int i = 0;i < ParDefinition.im;i ++){
				for(int j = 0;j < ParDefinition.jm;j ++){
					bw.write(ParDefinition.x[i][j] + "\t\t" +ParDefinition.y[i][j] + "\t\t" + ParDefinition.just[i][j]+" ");
		            bw.newLine();//����
				}		
			}
	               
			bw.close();
			fw.close();

    	  } catch (IOException e) {
    	   e.printStackTrace();
    	  }   	 		
	}
	
	
}
