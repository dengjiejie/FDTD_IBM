package cn.edu.tju.main;

import java.io.IOException;
import java.lang.String;

public class TestProject {

	/**
	 * @param args
	 * @throws IOException
	 */
	static double E = 2.71828183;
	public static double[] dx = new double[ParDefinition.nlp]; // 存放拉格朗日点坐标的增量
	public static double[] dy = new double[ParDefinition.nlp];

	public TestProject() throws IOException {
		// TODO Auto-generated constructor stub
		// period of sound source
		double Tposs = 0.005, dT = ParDefinition.dt; // unit=s
		int SetPS = 0, OneGuassion = 0;
		long FS = 20000, P4signal = Math.round(1.0 / FS / dT), P400signal;
		long Nposs0 = Math.round(Tposs / dT), Nposs = 1; // computing steps;
															// time
															// step=0.0000002
		System.out.println("Tposs=" + Tposs + ";  Nposs=" + Nposs
				+ "; OneGaussion=" + OneGuassion);

		// ******Prepare the glottal air flow (Fant flow)
		int GAFn, GafFlag = 1;//GafFlag没有使用
		double Ap, a, b, t1, t2, t, pi = 3.1415926; //声波参数

		GAFn = Math.round(Nposs / 2);
		GAFn = 5000;
		if (GAFn > 15000)
			GAFn = 15000;
		Ap = 1;
		t2 = GAFn * dT;
		t1 = 2 * t2 / 3;
		a = pi / t1;
		b = 1 / (1 - Math.cos(pi * t2 / t1));
		for (int i = 0; i < GAFn; i++) {
			t = i * dT;
			ParDefinition.glottalAF[i] = 0;
			if (t < t1) {
				ParDefinition.glottalAF[i] = Ap * (1 - Math.cos(a * t));
			} else {
				ParDefinition.glottalAF[i] = Ap
						* (1 - b + b * Math.cos(a * (t - t1)));
			}
		}
		// ******END of Preparing the glottal air flow (Fant flow)

		// ************ define VT data file name root
		double MPS = 0.5;// the boundary moving speed= meter per second (MPS)
							// (max. MPS=0.33 for speech)
		int dataSet = 1, VTfn = 0, VTFN = 11, Ndiv = 1, VTnn = 0, aFPS = 5;//Ndiv没有使用
		double VTduration[] = { 0.05, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007,
				0.007, 0.007, 0.007, 0.007, 0.05, 0.03 }; // Unit:sec
		String VTfilenameRoot = "src/cn/edu/tju/main/";

		switch (dataSet) {
		case 1:
			VTfilenameRoot = VTfilenameRoot + "VTSCa-i/SCvocal_";
			VTFN = 12;
			aFPS = 2;// FPS=4000;
			break;
		case 2:
			VTfilenameRoot = VTfilenameRoot + "VTSCi-a/SCvocal_";
			VTFN = 12;
			aFPS = 2;// FPS=4000;
			break;
		case 3:
			VTfilenameRoot = VTfilenameRoot + "VTWSa-i/WSvocal_";
			// VTduration={0.05, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007,
			// 0.007, 0.007, 0.007, 0.05, 0.03}; //Unit:sec
			VTFN = 12;
			aFPS = 2;
			break;
		case 4:
			VTfilenameRoot = VTfilenameRoot + "VTuniform/Uniform_";
			VTduration[0] = 0.01; // Unit:sec
			VTduration[1] = 0.01; // Unit:sec
			VTduration[2] = 0.02; // Unit:sec
			VTduration[3] = 0.02; // Unit:sec
			VTFN = 4;
			aFPS = 2; // FPS=10000;
			break;
		case 5:
			//VTfilenameRoot = VTfilenameRoot + "VTuniform/ConvergeNN_";
			VTfilenameRoot = VTfilenameRoot + "VTuniform/convergeSTP0.0015_";
			if (MPS == 1.0) {
				VTduration[0] = 0.01; // Unit:sec
				VTduration[1] = 0.01; // Unit:sec
				VTduration[2] = 0.01; // Unit:sec
			} else if (MPS == 0.5) {
				VTduration[0] = 0.02; // Unit:sec
				VTduration[1] = 0.02; // Unit:sec
				VTduration[2] = 0.02; // Unit:sec
			} else if (MPS == 0.33) {
				VTduration[0] = 0.001; // Unit:sec
				VTduration[1] = 0.03; // Unit:sec
				VTduration[2] = 0.01; // Unit:sec
			}
			VTFN = 1;
			aFPS = 2; // FPS=10000;
			break;
		case 6:
			VTfilenameRoot = VTfilenameRoot + "VTuniform/DevergeN_";
			if (MPS == 1.0) {
				VTduration[0] = 0.02; // Unit:sec
				VTduration[1] = 0.01; // Unit:sec
				VTduration[2] = 0.02; // Unit:sec
			} else if (MPS == 1.5) {
				VTduration[0] = 0.0217; // Unit:sec
				VTduration[1] = 0.0066; // Unit:sec
				VTduration[2] = 0.0217; // Unit:sec
			} else if (MPS == 0.5) {
				VTduration[0] = 0.02; // Unit:sec
				VTduration[1] = 0.02; // Unit:sec
				VTduration[2] = 0.02; // Unit:sec
			}
			VTFN = 3;
			aFPS = 2; // FPS=10000;
			break;
		case 7:
			VTfilenameRoot = VTfilenameRoot + "VTuniform/DConvergeN_";
			if (MPS == 1.0) {
				VTduration[0] = 0.01; // Unit:sec
				VTduration[1] = 0.01; // Unit:sec
				VTduration[2] = 0.01; // Unit:sec
				VTduration[3] = 0.01; // Unit:sec
				VTFN = 4;
			} else if (MPS == 0.5) {
				VTduration[0] = 0.001; // Unit:sec
				VTduration[1] = 0.02; // Unit:sec
				VTduration[2] = 0.02; // Unit:sec
				VTduration[3] = 0.02; // Unit:sec
				VTduration[4] = 0.02; // Unit:sec
				VTFN = 5;
			}
			aFPS = 2; // FPS=10000;
			break;
		}
		String fn;
		// set the time for getting VT (vocal tract) data
		double TT = 0;
		long VTstep[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 5, 6, 7, 8, };
		for (int i = 0; i < VTFN; i++) {
			VTstep[i] = Math.round(VTduration[i] / dT); // step number of each
														// VTTT
			TT += VTduration[i];
			VTnn += VTstep[i];
			fn = VTfilenameRoot + i + "0.txt";
			System.out.println("MPS=" +MPS+ ";[" + (i - 1) + "," + i + "]=" + VTstep[i]
					+ ";\t dataFile: " + fn);
		}
		ParDefinition.tend = VTnn;
		VTstep[VTFN] = ParDefinition.tend + 1; // stop reading VT file
		VTnn = (int) VTstep[0];
		System.out.println("Total step =" + ParDefinition.tend
				+ ";\t Total Time=" + TT + "s;" +  "\t dt="  + ParDefinition.dt 
				+ "s;" + "\t dx=" + ParDefinition.dx + "m;" + "\t dy=" + ParDefinition.dy +"m"  );
		System.out
				.println("Calculation steps in the first petiod:" + VTstep[0]);

		// 第一步：建立拉格朗日点对象lp，读取拉格朗日点
		fn = VTfilenameRoot + "0.txt";
		System.out.println(fn);
		LagPoint lp = new LagPoint();
		System.out.println("----读取拉格朗日点----");
		lp.readLagFile(fn);// read the first data file
		lp.readLagFile01(fn);// read the first data file again as the next one
								// for calculating the increment
		// System.out.println("------读取完毕------");

		// 第二步：生成欧拉网格
		// 若之前未生成过欧拉网格，则进行网格生成以及初始化操作；否则，直接读取网格即可。本程序中，默认网格未生成。
		
		
		
		
		
		
		if (true) {
			// 欧拉网格，△x=△y=0.02
			EulaPoint ep = new EulaPoint();
			// System.out.println("*****网格生成中*****");
			ep.eulaGeneration();
			// System.out.println("*****网格已生成*****");

			// 设置初始条件
			SetCondition sc = new SetCondition();
			sc.initCondition();// 初始波源
			sc.boundCondition();
			System.out.println("*******初始化完成******");
		}

		int fileCount = 1;
		SetCondition sc = new SetCondition();

		// found the central end of VT
		int minX = (int) ((0 - ParDefinition.rl) / ParDefinition.dx);
		int minY = (int) ((0 - ParDefinition.rd) / ParDefinition.dy);
	
		
		System.out.println("VTend=" + minX + " " + minY + ";  V="
				+ ParDefinition.x[minX][minY] + " "
				+ ParDefinition.y[minX][minY]);

		for (int i = 0; i < ParDefinition.nlp; i++) { // Increment between two
														// VTs
			dx[i] = (ParDefinition.x0l[i] - ParDefinition.xl[i]) / VTstep[0];
			dy[i] = (ParDefinition.y0l[i] - ParDefinition.yl[i]) / VTstep[0];
		}

		// 第三步：主要操作即对每一时间步进行处理
		TimeOperation to = new TimeOperation();
		// 完美匹配层的设置
		to.setPML();
		
		
		
		
		// 下面这部分代码主要是为了保存上一个时间步计算之后的值，方便下次计算（因为fdtd是交替计算的）
		for (int itime = ParDefinition.tfrom; itime <= 250000; itime++) {

			// 对每一个欧拉网格点，将上一时间步所得到的uvp值赋值给u0、v0、p0
			for (int i = 0; i < ParDefinition.im; i++) {
				for (int j = 0; j < ParDefinition.jm; j++) {
					ParDefinition.u0[i][j] = ParDefinition.u[i][j];
					ParDefinition.v0[i][j] = ParDefinition.v[i][j];
					ParDefinition.p0[i][j] = ParDefinition.p[i][j];
					ParDefinition.px0[i][j] = ParDefinition.px[i][j];
					ParDefinition.py0[i][j] = ParDefinition.py[i][j];
					ParDefinition.fx[i][j] = 0;
					ParDefinition.fy[i][j] = 0;
				}
			}
			// 将附加力以及压力倒数初始化，just用来记录网格点是否在圆柱体内
			for (int j = 0; j < ParDefinition.jm; j++) {
				for (int i = 0; i < ParDefinition.im; i++) {
					ParDefinition.fx[i][j] = 0;
					ParDefinition.fy[i][j] = 0;
					ParDefinition.fxp[i][j] = 0;
					ParDefinition.fyp[i][j] = 0;
					ParDefinition.just[i][j] = 0;
					ParDefinition.tube[i][j] = 0;
				}
			}
			// 移动边界时处理的是元音/a/的变化，若做其他的元音，该处代码可以注释
			// 由于实验设置的问题目前，同时通过计算，几乎1-2个网格间距，声道变化特别小，
			// 20ms
			// set periodical sound source
			if (OneGuassion == 1) { // 1: sound source is one Gaussion, the same
									// as the initial one
				if (--Nposs == 0) {
					// Nposs=Nposs0+(int) ((Math.random()-0.5)*400);//add a
					// perturbation on the period
					sc.setPeriodicalSource(minX, minY);
					System.out.println(" One-Gaussion source at  " + itime
							+ " ------");
				}
			} else { // 0: Test the source with multiple Gausion
				int Adding = 1; /*
								 * Adding=0: reset the end part of VT by Fant
								 * flow; Adding=1: Add Fant flow on the current
								 * pressure
								 */
				if (--Nposs == 0) {//Nposs 可以通过每次的--Nposs来使Nposs为零时，加入Fant声波
					Nposs = Nposs0 + (int) ((Math.random() - 0.5) * 400);// add
																			// a
																			// perturbation
																			// on
																			// the
																			// period
					SetPS = GAFn;
					System.out.println("Fant flow source at " + itime
							+ "; Adding=" + Adding + "; GAFn=" + GAFn
							+ "; Nposs=" + Nposs);
				}
				if (SetPS > 1) {// Set the shape of sound source
					sc.setPeriodicalSourceTest(minX, minY,
							ParDefinition.glottalAF[GAFn - SetPS], Adding);
					SetPS--;
				}
			}
			if (itime > VTnn) {
				VTfn++;//进入12阶段中的下一阶段
				VTnn += VTstep[VTfn];//VTnn记录什么时候进入下一阶段
				fn = VTfilenameRoot + (VTfn - 1) + ".txt";//VTfn减1是因为VTfn提前加了1
				LagPoint lp1 = new LagPoint();
				System.out.println("----读取拉格朗日点: " + fn + "----");//为什么眼读取两次
				lp1.readLagFile(fn);
				fn = VTfilenameRoot + VTfn + ".txt";
				System.out.println("----读取拉格朗日点: " + fn + "----");
				lp1.readLagFile01(fn);
				System.out.println("Calculation steps of this petiod:"
						+ VTstep[VTfn]);
				
				
				
				for (int i = 0; i < ParDefinition.nlp; i++) { // Increment between two VTs
					dx[i] = (ParDefinition.x0l[i] - ParDefinition.xl[i])
							/ VTstep[VTfn];
					dy[i] = (ParDefinition.y0l[i] - ParDefinition.yl[i])
							/ VTstep[VTfn];
				}
			}
			for (int i = 0; i < ParDefinition.nlp; i++) { // Add the Increment
															// between two VTs
				ParDefinition.xl[i] += dx[i];
				ParDefinition.yl[i] += dy[i];
			}

			
			
			
			
			// 第四步：计算附加力，更新 p u v
			to.virForce();
			// FDTD求解puv
			to.solvePUV();
			
			
			
			
			

			double observationX6 = 0.19, observationY6 = 0;
			double observationX1 = 0.185, observationY1 = 0;
			double observationX0 = 0.005, observationY0 = 0;
			String fileName6 = "src/cn/edu/tju/main/values/observation("
					+ observationX6 + "," + observationY6 + ")p.dat";
			int ox6 = (int) ((observationX6 - ParDefinition.rl) / ParDefinition.dx);
			int oy6 = (int) ((observationY6 + ParDefinition.rt) / ParDefinition.dy);
			String fileName1 = "src/cn/edu/tju/main/values/observation("
					+ observationX1 + "," + observationY1 + ")p.dat";
			int ox1 = (int) ((observationX1 - ParDefinition.rl) / ParDefinition.dx);
			int oy1 = (int) ((observationY1 + ParDefinition.rt) / ParDefinition.dy);
			String fileName0 = "src/cn/edu/tju/main/values/observation("
					+ observationX0 + "," + observationY0 + ")p.dat";
			int ox0 = (int) ((observationX0 - ParDefinition.rl) / ParDefinition.dx);
			int oy0 = (int) ((observationY0 + ParDefinition.rt) / ParDefinition.dy);

			// Low pass filter p(n)=lps*p(n-1)+ParDefinition.p(n)
			double p66 = 0, p11 = 0, p00 = 0, lps = 0.08;
			p66 = lps * p66 + ParDefinition.p[ox6][oy6];
			p11 = lps * p11 + ParDefinition.p[ox1][oy1];
			p00 = lps * p00 + ParDefinition.p[ox0][oy0];

			// 第五步：将计算结果进行保存，也就是采样 ，每隔20个时间步，生成 一个输出文件，存放欧拉网格点三列值：x y p
			if ((itime % P4signal) == 0) {// FS=20kHz
				P400signal = P4signal * aFPS; // FPS=4000Hz
				if ((itime % P400signal) == 0) {
//					to.saveXYP(fileCount);
					fileCount++;
				}
				// 采样点
				to.saveValues_pNew(p66, fileName6);
				to.saveValues_pNew(p11, fileName1);
				to.saveValues_pNew(p00, fileName0);
				/*
				 * double observationX2 = 0.14,observationY2 = 0; String
				 * fileName2 = "src/cn/edu/tju/main/values/observation(" +
				 * observationX2 + "," + observationY2 + ")p.dat"; int
				 * ox2=(int)((observationX2 +0.02)/ParDefinition.dx); int
				 * oy2=(int)((observationY2 +0.04)/ParDefinition.dy);
				 * to.saveValues_p(ox2,oy2,fileName2);
				 * //to.saveValues_p(40,50,fileName2);
				 * 
				 * double observationX3 = 0.15,observationY3 = 0; String
				 * fileName3 = "src/cn/edu/tju/main/values/observation(" +
				 * observationX3 + "," + observationY3 + ")p.dat"; int
				 * ox3=(int)((observationX3 +0.02)/ParDefinition.dx); int
				 * oy3=(int)((observationY3 +0.04)/ParDefinition.dy);
				 * to.saveValues_p(ox3,oy3,fileName3);
				 * 
				 * double observationX4 = 0.17,observationY4 = 0; String
				 * fileName4 = "src/cn/edu/tju/main/values/observation(" +
				 * observationX4 + "," + observationY4 + ")p.dat"; int
				 * ox4=(int)((observationX4 +0.02)/ParDefinition.dx); int
				 * oy4=(int)((observationY4 +0.04)/ParDefinition.dy);
				 * to.saveValues_p(ox4,oy4,fileName4);
				 * 
				 * double observationX5 = 0.16,observationY5 = 0; String
				 * fileName5 = "src/cn/edu/tju/main/values/observation(" +
				 * observationX5 + "," + observationY5 + ")p.dat"; int
				 * ox5=(int)((observationX5 +0.02)/ParDefinition.dx); int
				 * oy5=(int)((observationY5 +0.04)/ParDefinition.dy);
				 * to.saveValues_p(ox5,oy5,fileName5);
				 */

				if ((itime % 5000) == 0) {
					System.out.println("------" + itime + "------");
				}
			}

		}// for结束，时间步处理结束
		System.out.println("---the END ::: Total time=" + TT + "s ------");
	}

	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		long startTime=System.currentTimeMillis();   //获取开始时间  
		 //测试的代码段  
		
		
		
		// 读取时间步长、初始时间步、结束时间步。
		ReadFile readInit = new ReadFile();
		String initString = readInit.readInitFile("src/cn/edu/tju/main/in.dat");
		String[] init = initString.split(" ");

		ParDefinition.dt = Double.parseDouble(init[0]);
		ParDefinition.tfrom = Integer.parseInt(init[1]);
		ParDefinition.tend = Integer.parseInt(init[2]);
		System.out.println("时间间隔：" + ParDefinition.dt + "  初始时间步："
				+ ParDefinition.tfrom + "  结束时间步：" + ParDefinition.tend);

		// 创建对象，调用默认的构造函数
		TestProject tp = new TestProject();
		
		
		
		long endTime=System.currentTimeMillis(); //获取结束时间  
		
		System.out.println("程序运行时间： "+formatDuring(endTime-startTime));  
	}
	
	
	public static String formatDuring(long mss) {  
	    long days = mss / (1000 * 60 * 60 * 24);  
	    long hours = (mss % (1000 * 60 * 60 * 24)) / (1000 * 60 * 60);  
	    long minutes = (mss % (1000 * 60 * 60)) / (1000 * 60);  
	    long seconds = (mss % (1000 * 60)) / 1000;  
	    return days + " days " + hours + " hours " + minutes + " minutes "  
	            + seconds + " seconds ";  
	}  
	
}
