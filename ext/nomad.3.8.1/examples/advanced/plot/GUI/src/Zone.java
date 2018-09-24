import java.awt.Dimension;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import javax.swing.JComponent;
import javax.swing.ToolTipManager;

public class Zone extends JComponent implements MouseListener
{
	private static final long serialVersionUID = 1L;
	public static int 	ID = 0;
	public int 			id;
	public float 		_x;
	public float 		_y;
	
	public Zone(float x, float y)
	{
		// Settings
		addMouseListener(this);
		setOpaque(false);
		setSize(new Dimension(7, 7));
		id = ID;
		++ID;
		// Coordinates
		_x = x;
		_y = y;
	}

	public void mouseClicked(MouseEvent e){}

	public void mouseEntered(MouseEvent e)
	{
		// Display info, paint the point in red and update the window
		ToolTipManager.sharedInstance().registerComponent(this);
		ToolTipManager.sharedInstance().setInitialDelay( 0 );
		ToolTipManager.sharedInstance().setDismissDelay(100000);
		if (Prog.isSingle)
			setToolTipText("<html> BBE : " + (int)_x + "<br/ Objective value : " + _y + "</html>");
		else
			setToolTipText("<html> f1 : " + _x + "<br/ f2 : " + _y + "</html>");
		Graph.couleur = id;
		Prog.graph.calcul();
	}

	public void mouseExited(MouseEvent e)
	{
		// Paint the point in black and update the window
		Graph.couleur = -1;
		Prog.graph.calcul();
	}

	public void mousePressed(MouseEvent e){}

	public void mouseReleased(MouseEvent e){}
}
