import java.awt.Color;
import java.awt.Graphics;
import java.awt.Point;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.Vector;
import javax.imageio.ImageIO;
import javax.swing.JComponent;
 
public class Graph extends JComponent
{
	private static final long serialVersionUID = 1L;
	public static int 		couleur = -1;
	public static boolean 	det = true;
	public static float 	// Pitch
							_x,
							_y,
							// Coordinates
							_x0,	// First point
							_y0,
							_x1,	// Last Point
	 						_yf,
	 						// First and last y-axe coordinates
				  			_yh,
				  			_yb,
				  			// First and last x-axe coordinates
				  			_xh,
				  			_xb;
	public static int		w = Fenetre.w,
							h = Fenetre.h;
	// List of Zone
	public static Vector<Zone> listeZ = new Vector<Zone>();
	
	// First call : graph with no point
	public Graph()
	{
		super();
		_x0 = 0;
		_y0 = 0;
		listeZ.removeAllElements();
		Prog.listeX.removeAllElements();
		Prog.listeY.removeAllElements();
		Prog.listeXX.removeAllElements();
		Prog.listeXP.removeAllElements();
		Prog.listeYP.removeAllElements();
	}
	
	// Add a point to the window
	@SuppressWarnings("deprecation")
	public Graph(int bbe, float y)
	{
		super();
		_x0 = Prog.listeX.firstElement();
		_y0 = Prog.listeY.firstElement();
		_yf = Prog.listeY.lastElement();
		
		// Calculate the new coordinate system
		calcul();
		
		// Add a zone corresponding to the new point
		Zone z = new Zone(bbe, y);
		z.setLocation(new Point((75 + (int)((bbe) * 65/_x) *w/800),
							((int)(480 - ((y - _yb) * 45/_y) *h/600))));
		add(z);
		listeZ.add(z);
		
		// Repaint the window
		resize(w,h);
	}
	
	@SuppressWarnings("deprecation")
	public Graph(Vector<Float> x, Vector<Float>y)
	{
		super();
		

		// Searching the min and max
		float min = x.firstElement(), max = min;
		
		for (int i = 1; i < x.size(); ++i)
		{
			if (x.elementAt(i) < min) min = x.elementAt(i);
			if (x.elementAt(i) > max) max = x.elementAt(i);
		}
		_x0 = min;
		_x1 = max;
		min = y.firstElement();
		max = min;
		for (int i = 1; i < y.size(); ++i)
		{
			if (y.elementAt(i) < min) min = y.elementAt(i);
			if (y.elementAt(i) > max) max = y.elementAt(i);
		}
		_y0 = max;
		_yf = min;
		calcul();
		for (int i = 0; i < x.size(); ++i)
		{
			Zone z = new Zone(x.elementAt(i), y.elementAt(i));
			z.setLocation(new Point((75 + (int)((x.elementAt(i)) * 65/_x) *w/800),
					(480 - (int)((y.elementAt(i) - _yb) * 45/_y )*h/600)));
			add(z);
			if (det)
			{
				listeZ.removeAllElements();
				det = false;
			}
			listeZ.add(z);
		}
		
		// Repaint the window
		resize(w,h);
	}
	
	// Method witch calculate the coordinate system
	public void calcul()
	{
		int cpt,
			test;
		boolean sor = false;
		// x-axe coordinates
		if (Prog.isSingle)
		{
			_x = Prog.listeX.lastElement();
			cpt = 1;
			test = 10;
	
			while (_x > test)
			{
				test *= 10;
				++cpt;
			}
			// The x-axe pitch is a multiple of 2 and 5
			while (	(_x % 20) != 0 || (_x % 50) != 0)
				++_x;
			_x /= 10;
		}
		else
		{
			_x = Math.abs(_x0 - _x1);
			while (!sor)
			{
				cpt = 0;
				if (_x != 0)
				{
					if(_x < 1)
						while (_x < 1)
						{
							_x *= 10;
							--cpt;
						}
					else 
						while (_x >= 10)
						{
							_x /= 10;
							++cpt;
						}
					int z = (int) _x;
					float tmp = _x;
					_x = z + _x/10;
					if (tmp != _x)
						_x = z+((int)((_x-z)*10) +1)/10;
										
					_x *= Math.pow(10, cpt-1);
	/*				if (_xb + 10 * _x < _x1)
					{
						System.out.println(_x);
						_x = 1 + Math.abs(_x0 - _x1) + Math.abs(_x1 - (_xb + 10*_x));
						System.out.println(_x);
					}
					else
					{
						sor = true;
					}
	*/
					sor = true;
				}
				else _x = 1;
			}
			
			// Calculate the left and right x-axe coordinates
			_xb = _x0;
			if (_xb > 0)
			{
				int i = 1;
				while (_xb > _x *i)
				{
					++i;
				}
				_xb = _x *(i-1);
			}
			else if (_xb < 0)
			{
				int i = 1;
				while (_xb < -_x *i)
				{
					++i;
				}
				_xb = -_x *i;
			}
		}

		// y-axe coordinate
		_y = Math.abs(_y0 - _yf);
		cpt = 0;
		
		if (_y != 0)
		{
			if(_y < 1)
				while (_y < 1)
				{
					_y *= 10;
					--cpt;
				}
			else 
				while (_y > 10)
				{
					_y /= 10;
					++cpt;
				}
			
			_y = (int) _y +1;
			_y *= Math.pow(10, cpt-1);
			
		}
		else _y = 1;
		
		// Calculate the top and bottom y-axe coordinates
		_yb = _yf;
		if (_yb > 0)
		{
			int i = 1;
			while (_yb > _y *i)
			{
				++i;
			}
			_yb = _y *(i-1);
		}
		else if (_yb < 0)
		{
			int i = 1;
			while (_yb < -_y *i)
			{
				++i;
			}
			_yb = -_y *i;
		}
		// Update the window
		invalidate();
		validate();
		repaint();
	}
	
	// Round up x-axe numbers
	public static float arrondi(float x)
	{
		x *= 100;
		float res;
		float tab[] = new float [3];
		tab[0] = (int) x;
		tab[1] =  (int)(x + 1);
		tab[2] = (int) (x - 1);
		if (Math.abs(x-tab[0]) < Math.abs(x-tab[1]))
			if (Math.abs(x-tab[0]) < Math.abs(x-tab[2]))
				res = tab[0]/100;
			else res = tab[2]/100;
		else if (Math.abs(x-tab[2]) < Math.abs(x-tab[1]))
				res = tab[2]/100;
			else res = tab[1]/100;
		
		return res;
	}

	// Paint the window
	public void paintComponent(Graphics g)
    {
		g.setColor(Color.WHITE);
		g.fillRect(0, 0, 3000, 3000);
		g.setColor(Color.BLACK);
		
		// Update the resolution
		h = Prog.f.getHeight();
		w = Prog.f.getWidth();
		
		// Trace the axes
		g.drawLine(65, 500*h/600,
					780*w/800, 500*h/600);

		g.drawLine(80*700/800, 10*h/600,
					80*700/800, 678*h/800);


		// Put the axes names
		BufferedImage image;
		if (Prog.isSingle)
		{
			try
			{
				image = ImageIO.read(new File("../obj.jpg"));
				g.drawImage(image,10*w/800,200*h/600,this);
			} 
			catch (IOException e)
			{
				e.printStackTrace();
			}
			try
			{
				image = ImageIO.read(new File("../bbe.jpg"));
				g.drawImage(image,300*w/800,530*h/600,this);
			} 
			catch (IOException e)
			{
				e.printStackTrace();
			}
		}
		else
		{
			try
			{
				image = ImageIO.read(new File("../f2.jpg"));
				g.drawImage(image,10*w/800,250*h/600,this);
			} 
			catch (IOException e)
			{
				e.printStackTrace();
			}
			try
			{
				image = ImageIO.read(new File("../f1.jpg"));
				g.drawImage(image,400*w/800,530*h/600,this);
			} 
			catch (IOException e)
			{
				e.printStackTrace();
			}
		}
		
		// Trace the spears
		g.setColor(Color.BLACK);
		int xx[] = {785*w/800,775*w/800, 775*w/800 };
		int yy[] = {500*h/600,495*h/600, 505*h/600};
		g.fillPolygon(xx, yy, 3);
		int xxx[] = {80*700/800,75*700/800,
				85*700/800 };
		int yyy[] = {5*h/600,15*h/600, 15*h/600};
		g.fillPolygon(xxx, yyy, 3);
		
		// Trace the scalings
		if (Prog.isSingle)
			for (int i = 1; i < 11; ++i)
				g.drawLine(70+(65*w/800)*i/1, (500-3)*h/600,
							70+(65*w/800)*i/1, (500+3)*h/600);
			else
				for (int i = 0; i < 11; ++i)
					g.drawLine(75+(i)*65*w/800, (500-3)*h/600,
							75+(i)*65*w/800, (500+3)*h/600);
	
		for (int i = 1; i < 12; ++i)
			g.drawLine((80-3)*700/800, (530-i*45)*h/600,
						(80+3)*700/800, (530-i*45)*h/600);
		
		
		// Write scaling value
		if (Prog.isSingle)
		{
			for (int i = 1; i < 11; ++i)
				g.drawString(String.valueOf((int)_x*(i)), 60+(i*65)*w/800,
						525*h/600);
		}
		else 
		for (int i = 0; i < 11; ++i)
		{
			g.drawString(String.valueOf(arrondi(_xb + _x*i)), 65+(i*65)*w/800,
					525*h/600);
		}
		
		for (int i = 0; i < 11; ++i)
		{
			if (arrondi(_yb + _y*i) >= 1000 || arrondi(_yb + _y*i) <= -100)
				g.drawString(String.valueOf(arrondi(_yb + _y*i)), 25*700/800,
							(490-45*i)*h/600);
			else if (arrondi(_yb + _y*i) >= 100 || arrondi(_yb + _y*i) <= -10)
				g.drawString(String.valueOf(arrondi(_yb + _y*i)), 30*700/800,
							(490-45*i)*h/600);
			else if (arrondi(_yb + _y*i) >= 10 || arrondi(_yb + _y*i) < 0)
				g.drawString(String.valueOf(arrondi(_yb + _y*i)), 40*700/800,
							(490-45*i)*h/600);
			else
				g.drawString(String.valueOf(arrondi(_yb + _y*i)), 50*700/800,
							(490-45*i)*h/600);
		}
		
		// Place the Zones
		g.setColor(Color.GREEN);
		if (Prog.isSingle)
			for (int i=0; i < listeZ.size(); ++i)
	   		{
	   			listeZ.elementAt(i).setLocation(new Point
	   					((int)(67 + Prog.listeX.elementAt(i)*65*w/_x/800),
	   					(int)((483 - (Prog.listeY.elementAt(i)- _yb) * 45/_y))*h/600));
	   			add(listeZ.elementAt(i));
	   		}
		else
		{
			for (int i=0; i < listeZ.size(); ++i)
	   		{
	   			listeZ.elementAt(i).setLocation(new Point
	   					((int)(73 + (Prog.listeXX.elementAt(i) - _xb)*65*w/_x/800),
	   					(int)((485 - ((Prog.listeY.elementAt(i)- _yb) * 45/_y)))*h/600));
	   			add(listeZ.elementAt(i));
	   		}
		}
		g.setColor(Color.BLACK);

		// Paint the points
		if (Prog.isSingle)
			for (int i = 0; i < Prog.listeX.size() ; ++i)
			{
				if (couleur == i) 
	    		   g.setColor(Color.RED);
	    	   	g.fillRect((int)(70 + Prog.listeX.elementAt(i)*65*w/_x/800),
	    	   		(int)((484 - ((Prog.listeY.elementAt(i) - _yb) * 45/_y))*h/600),3, 3);
	    	   	if (g.getColor()==Color.RED) g.setColor(Color.BLACK);
			}
		else
		{
			for (int i = 0; i < Prog.listeXX.size() ; ++i)
			{
				if (couleur == i) 
	    		   g.setColor(Color.RED);
	    	   	g.fillRect((int)(75 + (((Prog.listeXX.elementAt(i)-_xb)*(65*w/700)/_x))
	    	   			* 700/800), 
	    	   		(int)((484 - ((Prog.listeY.elementAt(i) - _yb) * 45/_y))*h/600), 3, 3);
	    	   	if (g.getColor()==Color.RED) g.setColor(Color.BLACK);
			}
			g.setColor(Color.RED);
			for (int i = 0; i < Prog.listeXP.size() ; ++i)
			{  
	    	   	g.fillRect((int)(75+(((Prog.listeXP.elementAt(i)-_xb)*(65*w/700)/_x))
	    	   			* 700/800), 
	    	   		(int)((484 - ((Prog.listeYP.elementAt(i) - _yb) * 45/_y))*h/600), 3, 3);  	
			}
			g.setColor(Color.BLACK);
		}
		
		// Trace the stepping function
		if (Prog.isSingle)
		{
			for (int i = 0; i < Prog.listeX.size() -1 ; ++i)
			{
				g.setColor(Color.BLUE);
				g.drawLine((int)(71 + Prog.listeX.elementAt(i)*65*w/_x/800), 
	    			   (int)((485 - ((Prog.listeY.elementAt(i) - _yb) * 45/_y))*h/600), 
	    			   (int)(71 + Prog.listeX.elementAt(i+1)*65*w/_x/800),
	    			   (int)((485 - ((Prog.listeY.elementAt(i) - _yb) * 45/_y))*h/600));
			}
			for (int i = 0; i < Prog.listeX.size() -1 ; ++i)
			{
				g.drawLine((int)(71 + Prog.listeX.elementAt(i+1)*65*w/_x/800), 
	    			   (int)((485 - ((Prog.listeY.elementAt(i) -_yb) * 45/_y))*h/600), 
	    			   (int)(71 + Prog.listeX.elementAt(i+1)*65*w/_x/800),
	    			   (int)((485 - ((Prog.listeY.elementAt(i+1) -_yb) * 45/_y))*h/600));
			}
		}
    }
}
