package cn.edu.tju.main;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;


public class LagPoint {
	String tempString = null;
	static double PI = 3.14159265;
	double xxx = PI * 2.0 / 3;
	String fileName0 = "src/cn/edu/tju/main/test.txt";
	String fileName1 = "src/cn/edu/tju/main/lptest.dat";
	//String fileName = "src/cn/edu/tju/main/VT-1.dat";
	//String fileName = "src/cn/edu/tju/main/VT_O.txt";
	//String fileName = "src/cn/edu/tju/main/LH_i_vocaltract.txt";
	//String fileName = "src/cn/edu/tju/main/CR_i_vocaltract.txt";
	//String fileName = "src/cn/edu/tju/main/vocal_a.txt";
	double xlp[] = new double[183];
	double ylp[] = new double[183];
	int count = 182;
	double tempsize = 0.05;

	//生成拉格朗日点，将圆柱以360个点进行表示，存放于xl[i]、yl[i]中
	public void writeLagFile() throws IOException {

		
		File inFile = new File(fileName0);
		BufferedReader reader = null;
		String tempString  = null;
		File outFile = new File(fileName1);
		// if file doesnt exists, then create it
		if (!outFile.exists()) {
			outFile.createNewFile();
		}
			
		FileWriter fw = new FileWriter(outFile.getAbsoluteFile());
		BufferedWriter bw = new BufferedWriter(fw);
		
		
		double tempx = 0, tempy = 0, temp = 0, length, mod_y,height;
		int num_y;
		
		bw.write(temp + "\t" + temp);
		bw.newLine();//换行
		/*xlp[count] = temp;
		ylp[count] = temp;
		count ++;
	*/
		reader = new BufferedReader(new FileReader(inFile));
		
		// 一次读入一行，直到读入null为文件结束
		while ((tempString = reader.readLine()) != null) {
			String[] lp = tempString.split(" ");
			tempx = Double.parseDouble(lp[0]) * 100;
			tempy = Double.parseDouble(lp[1]) * 100;
			height = temp - tempy;
			length = Math.sqrt(Math.pow(tempy - temp, 2) + Math.pow(tempsize, 2));
			num_y = (int)(length / tempsize);
		//	System.out.println(num_y);
			for(int i = 0;i < num_y;i ++){
				temp = temp + tempsize;
				bw.write(tempx + "\t" + (-temp));
				bw.newLine();//换行
			/*	xlp[count] = temp;
				ylp[count] = temp;
				count ++;*/
			}
			
			mod_y = height - num_y * tempsize;
			if(mod_y != 0){
				temp = tempy;
				bw.write(tempx + "\t" + (-temp));
				bw.newLine();//换行
				/*xlp[count] = temp;
				ylp[count] = temp;
				count ++;*/
			}
			
        }
		reader.close();
		
		for(int i=85;i <169;i++){
			bw.write(xlp[168 - i] + "\t" + ylp[168 - i]);
			bw.newLine();//换行
		}
		bw.close();
		fw.close();

    	 
		
	}
	
	
	//读取拉格朗日点信息
	public String readLagFile(String fileName) {
		File inFile = new File(fileName);
		BufferedReader reader = null;
		String tempString  = null;
		int i = 0;
		try {
			reader = new BufferedReader(new FileReader(inFile));
			// 一次读入一行，直到读入null为文件结束
			while ((tempString = reader.readLine()) != null) {
				String[] lp = tempString.split("\t");
				ParDefinition.xl[i] = Double.parseDouble(lp[0]);// * 100
				ParDefinition.yl[i] = Double.parseDouble(lp[1]);//* 100
				i ++;
	            }
	            reader.close();
	        } catch (IOException e) {
	            e.printStackTrace();
	        } finally {
	            if (reader != null) {
	                try {
	                    reader.close();
	                } catch (IOException e1) {
	                }
	            }
	        }
		return tempString;
	    }
	
	public String readLagFile01(String fileName) {
		File inFile = new File(fileName);
		BufferedReader reader = null;
		String tempString  = null;
		int i = 0;
		try {
			reader = new BufferedReader(new FileReader(inFile));
			// 一次读入一行，直到读入null为文件结束
			while ((tempString = reader.readLine()) != null) {
				String[] lp = tempString.split("\t");
				ParDefinition.x0l[i] = Double.parseDouble(lp[0]);// * 100
				ParDefinition.y0l[i] = Double.parseDouble(lp[1]);//* 100
				i ++;
	            }
	            reader.close();
	        } catch (IOException e) {
	            e.printStackTrace();
	        } finally {
	            if (reader != null) {
	                try {
	                    reader.close();
	                } catch (IOException e1) {
	                }
	            }
	        }
		return tempString;
	    }
	
}
