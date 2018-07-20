package cn.edu.tju.main;


public class ParDefinition {
	/**
	 * ������Ҫ�õ��ĸ��ֲ���
	 */
	public static double dt; 			//ʱ�䲽��(����޸�dt�ǵ��޸ĸ�˹��Դ����dt)
	public static int tfrom,tend;		//��ʼʱ�䡢����ʱ��
	public static double dx = 0.001;		//ŷ���������
	public static double dy = 0.001;
	public static double dxy = Math.sqrt(Math.pow(dx, 2) + Math.pow(dy, 2));
	//public static double rl = -0.02,rr = 0.20,rd = -0.04,rt = 0.04;	//�������С��Fu Tang
	public static double rl = -0.01,rr = 0.24,rd = -0.04,rt = 0.04;	//�������С��Jianwu Dang
	public static int im = (int)((rr-rl)/dx)+1;		//x������ŷ�������ĸ���
	public static int jm = (int)((rt-rd)/dy)+1;		//y������ŷ�������ĸ���
	public static int nlp = 359;//293(a)		//Բ���ĵ��������������յ㣩
	public static int si = (int)((rr-rl)/dx), sj = (int)(rt/dy);
	public static double nflag = 0;
	public static double[] gp = new double[500];
	public static double[] glottalAF = new double[15000]; // glottal air flow

	public static double[] xl=new double[nlp]; 			//����������յ�����
	public static double[] yl=new double[nlp]; 
	public static double[] x0l=new double[nlp]; 			//�����һ��ʱ����������յ�����
	public static double[] y0l=new double[nlp]; 
	public static double[][] x = new double[im][jm];  	//���ŷ��������
	public static double[][] y = new double[im][jm];  
	public static double[][] u = new double[im][jm];	//���ŷ�����u��v��p
	public static double[][] v = new double[im][jm];
	public static double[][] p = new double[im][jm];
	public static double[][] px = new double[im][jm];
	public static double[][] py = new double[im][jm];
	public static double[][] u0 = new double[im][jm];	//�����һʱ�䲽��ŷ�����u��v��p	
	public static double[][] v0 = new double[im][jm];
	public static double[][] p0 = new double[im][jm];
	public static double[][] px0 = new double[im][jm];
	public static double[][] py0 = new double[im][jm];
	public static double[][] fx = new double[im][jm];	//������
	public static double[][] fy = new double[im][jm];
	public static double[][] fxp = new double[im][jm];	//ѹ������
	public static double[][] fyp = new double[im][jm];
	public static int[][] just = new int[im][jm];		//���������Ƿ���Բ�����⣬������Ϊ0��������Ϊ1
	public static int[][] tube = new int[im][jm];
	
/*	public static int npml = 8;
	public static double c_air = 1, rho_air = 1,alpha = 1;
	public static double flow_resistance = alpha * Math.pow(rho_air, 2) * Math.pow(c_air, 2); 
	public static double kappa = 1.0 / (rho_air * Math.pow(c_air, 2));*/
}
