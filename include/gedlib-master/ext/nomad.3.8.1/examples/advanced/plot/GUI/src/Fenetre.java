import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import javax.swing.JFrame;
 
public class Fenetre extends JFrame
{	
	private static final long serialVersionUID = 1L;
	public static Graph graph;
	public static int 	w,
						h;
	 
    public Fenetre() throws NumberFormatException, IOException
    {
    	// Window settings
        this.setTitle("NOMAD");
        w = 600;
        h = 400;
        // Check for special resolution
        File f = new File("../resolution.txt");
		if ( f.exists() ) 
		{
			BufferedReader lecteurAvecBuffer = null;
		    try
		    {
		    	lecteurAvecBuffer = new BufferedReader(new FileReader("../../resolution.txt"));
		     }
		     catch(FileNotFoundException exc)
		     {
		    	System.out.println("Error");
		     }
		    w = Integer.parseInt(lecteurAvecBuffer.readLine());
		    h = Integer.parseInt(lecteurAvecBuffer.readLine());
		    lecteurAvecBuffer.close();
		}    
        this.setSize(w, h);
        this.setLocation(0,0);
        this.setContentPane(graph = new Graph());
        this.setVisible(true);
    }
}
