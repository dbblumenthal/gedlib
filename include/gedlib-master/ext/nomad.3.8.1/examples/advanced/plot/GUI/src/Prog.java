import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Vector;

public class Prog
{
	public static  Fenetre 	f;
	public static boolean 	isSingle = false;	//single or bi-objective
	public static boolean 	isOk = false;		//ready to run
	public static boolean 	fin = false;		//end of a run
	public static boolean 	det = true;
	// Vector of points
	public static Vector<Integer> 	listeX = new Vector<Integer>();
	public static Vector<Float> 	listeY = new Vector<Float>();
	public static Vector<Float> 	listeXX = new Vector<Float>();
	public static Vector<Float> 	listeYP = new Vector<Float>();
	public static Vector<Float> 	listeXP = new Vector<Float>();
	
	public static Graph graph;
	
	public static void main(String[] args) throws IOException
	{
		// Cleaning old files
		File h = new File("../../init.txt");
		if ( h.exists() ) h.delete();
		h = new File("../../in.txt");
		if ( h.exists() ) h.delete();
		isOk = false;
		fin = false;

		while(true)
		{
			// Wait the c++ program
			while (!isOk)
			{
				h = new File("../../init.txt");
				if ( h.exists() ) 
				{
					BufferedReader lecteurAvecBuffer = null;
				    try
				    {
				    	lecteurAvecBuffer = new BufferedReader(new FileReader("../../init.txt"));
				    }
				    	catch(FileNotFoundException exc)
				    {
				    	System.out.println("Error");
				    }
				    isOk = true;
				    // Choose single or bi-objective
					String s = lecteurAvecBuffer.readLine();
//					System.out.println(s);
				    if(Integer.parseInt(s) == 1) isSingle = true;
				    else isSingle = false;
				    
				    lecteurAvecBuffer.close();
					h.delete();
				}
			}
			// Creating window
			f = new Fenetre();
			h = new File("../../out.txt");
			if ( !h.exists() ) h.createNewFile();
			while(!fin)
			{
				// Waiting the next points
				File g = new File("../../in.txt");
				if (g.exists())
				{
					BufferedReader lecteurAvecBuffer = null;
	
				    try
				    {
				    	lecteurAvecBuffer = new BufferedReader(new FileReader("../../in.txt"));
				    }
				    catch(FileNotFoundException exc)
				    {
				    	System.out.println("Error");
				    }
				    // Reading the points
				    if (isSingle)
				    {
				    	int x = Integer.parseInt(lecteurAvecBuffer.readLine());
				    	float y = Float.parseFloat(lecteurAvecBuffer.readLine());
				    	listeX.addElement(x);
				    	listeY.addElement(y);
				    	// Update the window
				    	f.setContentPane(graph = new Graph(x, y));
						
				    }
			    	else
			    	{
						
			    		int n = Integer.parseInt(lecteurAvecBuffer.readLine());
			    		for (int i = 0; i < n; ++i)
			    		{
			    			listeXX.addElement(Float.parseFloat(lecteurAvecBuffer.readLine()));
			    			listeY.addElement(Float.parseFloat(lecteurAvecBuffer.readLine()));
			    			Graph.det = true;
			    		}
			    		// Update the window
			    		f.setContentPane(graph = new Graph(listeXX, listeY));
			    		n = Integer.parseInt(lecteurAvecBuffer.readLine());
			    		det = true;
			    		for (int i = 0; i < n; ++i)
			    		{
			    			float f1 = Float.parseFloat(lecteurAvecBuffer.readLine());
			    			if (det)
			    			{
			    				listeXP.removeAllElements();
			    			}
			    			listeXP.addElement(f1);
			    			float f2 = Float.parseFloat(lecteurAvecBuffer.readLine());
			    			if (det)
			    			{
			    				listeYP.removeAllElements();
			    				det = false;
			    			}
			    			listeYP.addElement(f2);
			    		}
			    	}
				    lecteurAvecBuffer.close();

					File destination = new File("../../out.txt");
					g.renameTo(destination);
				}
				g = new File("../../stop.txt");
				if (g.exists())
				{
					fin = true;
					g.delete();
				}
			}
			fin = false;
			h = new File("../../init.txt");
			if ( h.exists() ) h.delete();
			h = new File("../../in.txt");
			if ( h.exists() ) h.delete();
			isOk = false;
		}
	}
}