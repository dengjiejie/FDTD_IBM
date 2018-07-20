package cn.edu.tju.main;


public class ParDefinition {
	/**
	 * 定义所要用到的各种参数
	 */
	public static double dt; 			//时间步长(如果修改dt记得修改高斯波源处的dt)
	public static int tfrom,tend;		//初始时间、结束时间
	public static double dx = 0.001;		//欧拉网格点间距
	public static double dy = 0.001;
	public static double dxy = Math.sqrt(Math.pow(dx, 2) + Math.pow(dy, 2));
	//public static double rl = -0.02,rr = 0.20,rd = -0.04,rt = 0.04;	//计算域大小：Fu Tang
	public static double rl = -0.01,rr = 0.24,rd = -0.04,rt = 0.04;	//计算域大小：Jianwu Dang
	public static int im = (int)((rr-rl)/dx)+1;		//x方向上欧拉网格点的个数
	public static int jm = (int)((rt-rd)/dy)+1;		//y方向上欧拉网格点的个数
	public static int nlp = 359;//293(a)		//圆柱的点数（即拉格朗日点）
	public static int si = (int)((rr-rl)/dx), sj = (int)(rt/dy);
	public static double nflag = 0;
	public static double[] gp = new double[500];
	public static double[] glottalAF = new double[15000]; // glottal air flow

	public static double[] xl=new double[nlp]; 			//存放拉格朗日点坐标
	public static double[] yl=new double[nlp]; 
	public static double[] x0l=new double[nlp]; 			//存放下一个时间的拉格朗日点坐标
	public static double[] y0l=new double[nlp]; 
	public static double[][] x = new double[im][jm];  	//存放欧拉点坐标
	public static double[][] y = new double[im][jm];  
	public static double[][] u = new double[im][jm];	//存放欧拉点的u、v、p
	public static double[][] v = new double[im][jm];
	public static double[][] p = new double[im][jm];
	public static double[][] px = new double[im][jm];
	public static double[][] py = new double[im][jm];
	public static double[][] u0 = new double[im][jm];	//存放上一时间步、欧拉点的u、v、p	
	public static double[][] v0 = new double[im][jm];
	public static double[][] p0 = new double[im][jm];
	public static double[][] px0 = new double[im][jm];
	public static double[][] py0 = new double[im][jm];
	public static double[][] fx = new double[im][jm];	//附加力
	public static double[][] fy = new double[im][jm];
	public static double[][] fxp = new double[im][jm];	//压力倒数
	public static double[][] fyp = new double[im][jm];
	public static int[][] just = new int[im][jm];		//标记网格点是否在圆柱体外，在体外为0，在体内为1
	public static int[][] tube = new int[im][jm];
	
/*	public static int npml = 8;
	public static double c_air = 1, rho_air = 1,alpha = 1;
	public static double flow_resistance = alpha * Math.pow(rho_air, 2) * Math.pow(c_air, 2); 
	public static double kappa = 1.0 / (rho_air * Math.pow(c_air, 2));*/
}
