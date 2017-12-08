package teleworx.worldwind;

public class LineOfSightClearance {
	public LineOfSightClearance(double totalDist, double dist, double elevation, double fromAntennaHt, double toAntennaHt)
	{
		double clearance = 0;
		double losHeight = 0;
	    if (dist == 0) {
	        clearance = fromAntennaHt - elevation;
	        losHeight = fromAntennaHt;
	    }
	    else {
	        losHeight = (((toAntennaHt - fromAntennaHt)/totalDist)*dist + fromAntennaHt);
	        clearance = losHeight - elevation;
	    }
	    losHeight = Math.round(losHeight*100)/100.0;
	    clearance = Math.round(clearance*100)/100.0;
	    
	    this.losHeight = losHeight;
	    this.clearance = clearance;
	}
	
	public double losHeight;
	public double clearance;
}
