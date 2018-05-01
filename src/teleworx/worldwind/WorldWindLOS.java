package teleworx.worldwind;

import gov.nasa.worldwind.geom.*;

public class WorldWindLOS {
	public static final double _equatorialRadius = 6378137.0;
	public static final double _polarRadius = 6356752.3;


	public static double computeDistanceBetween(LatLon pathStartPoint, LatLon pathEndPoint) {
		double pathDist = LatLon.ellipsoidalDistance(pathStartPoint, pathEndPoint, _equatorialRadius, _polarRadius);
//		double pathDist = equatorialRadius * LatLon.linearDistance(pathStartPoint, pathEndPoint).getRadians();
		if (Double.isNaN(pathDist))
			pathDist = 0;
		
		return pathDist;
	}
	
	public LosResult calculateLOSResult(ElevationData[] elevationData, double fromAntProposedHt, double toAntProposedHt, LosParams losParams, boolean collectProfileData)
	{
		LosResult result = new LosResult(); 

		if (elevationData == null)
			return result;
		
		double fromAntennaHt = fromAntProposedHt + elevationData[0].elevation;
		double toAntennaHt = toAntProposedHt + elevationData[elevationData.length - 1].elevation;

		double dist = 0;
		double elevation = 0;
		double elevationWithClutter = 0;
		
		double LIGHTSPEED = 300000000.0;
		double lambda = LIGHTSPEED / (losParams.freqMHz * 1000000);	// = C/f (f = f * 1000000.0)
		
		double d1 = 0;
		double d2 = 0;
//		double clearance = 0;
		double clearanceF1 = 0;
		double clearanceF2 = 0;
		double D = 0;
		double ht1 = 0;
		double ht2 = 0;

		double AXIS = 6378137.0;
		double R = AXIS * losParams.kFactor;	// WGS84 Ellipsoid, Radius of Earth in m.
	
		LatLon pathStartPoint = elevationData[0].location;
		LatLon pathEndPoint = elevationData[elevationData.length - 1].location;
		double pathDist = WorldWindLOS.computeDistanceBetween(pathStartPoint, pathEndPoint);
//        System.out.println(pathDist);
		if (pathDist <= 0) {
			result.minClearance = 0;
			result.minClearanceF1 = 0;
			result.minClearanceF1_Dist = 0;
			result.minClearanceF2 = 0;
			result.minClearanceF1LocLat = pathStartPoint.latitude.getDegrees();
			result.minClearanceF1LocLon = pathStartPoint.longitude.getDegrees();

			return result;
		}

		double m = (toAntennaHt - fromAntennaHt) / pathDist;
		double b = fromAntennaHt;
		
		double midPoint = pathDist / 2;
		double rise = R - Math.sqrt(R * R - midPoint * midPoint);

		pathDist = Math.round(pathDist)/1000.0; // convert to km
		for (int i = 0; i < elevationData.length; i++) {
			LatLon thisPoint = elevationData[i].location;
			dist = WorldWindLOS.computeDistanceBetween(pathStartPoint, thisPoint);
			
			// do fresnel zone calc for point.
			d1 = WorldWindLOS.computeDistanceBetween(pathStartPoint, thisPoint);
			d2 = WorldWindLOS.computeDistanceBetween(pathEndPoint, thisPoint);

			// Get Curvature factor
			D = Math.abs(midPoint - dist);
			double curvature = rise - (R - Math.sqrt(R * R - D * D));
			curvature = Math.round(curvature*100)/100.0;
			double signalCenterAtHeight = (d1 * m) + b;		// center of signal is at this ht.
			signalCenterAtHeight = Math.round(signalCenterAtHeight*100)/100.0;
			
			elevationWithClutter = elevationData[i].elevation + losParams.clutterHeight + curvature;
			elevationWithClutter = Math.round(elevationWithClutter*100)/100.0;

			elevation = elevationData[i].elevation + curvature;
			elevation = Math.round(elevation*100)/100.0;
			
			double F1Up = 0;
			double F1Dn = 0;
			double F2Up = 0;
			double F2Dn = 0;
			double freqClearanceF1 = 0;
			double freqClearanceF2 = 0;
			if (losParams.includeF1) {
			    result.includeF1 = true;
				clearanceF1 = Math.sqrt((losParams.fresnel1 * lambda * d1 * d2) / (d1 + d2));
				clearanceF1 = Math.round(clearanceF1*100)/100.0;
				freqClearanceF1 = (signalCenterAtHeight - clearanceF1) - elevationWithClutter ;
				freqClearanceF1 = Math.round(freqClearanceF1*100)/100.0;
				if (result.minClearanceF1 == -1 || result.minClearanceF1 > freqClearanceF1) {
					result.minClearanceF1 = freqClearanceF1;
					if (losParams.maxHopLengthKm == -1 || losParams.maxHopLengthKm <= 0) {
						result.minClearanceF1LocLat = thisPoint.latitude.getDegrees();
						result.minClearanceF1LocLon = thisPoint.longitude.getDegrees();
					}
				}
				if (losParams.maxHopLengthKm != -1 && losParams.maxHopLengthKm > 0) {
					if (d1 <= losParams.maxHopLengthKm * 1000 && d2 <= losParams.maxHopLengthKm * 1000) {
						if (result.minClearanceF1_Dist == -1 || result.minClearanceF1_Dist > freqClearanceF1) {
							result.minClearanceF1_Dist = freqClearanceF1;
							result.minClearanceF1LocLat = thisPoint.latitude.getDegrees();
							result.minClearanceF1LocLon = thisPoint.longitude.getDegrees();
						}
					}
				}				
				F1Up = signalCenterAtHeight + clearanceF1;
				F1Dn = signalCenterAtHeight - clearanceF1;
			}
			if (losParams.includeF2) {
			    result.includeF2 = true;
				clearanceF2 = Math.sqrt((losParams.fresnel2 * 2 * lambda * d1 * d2) / (d1 + d2));
				clearanceF2 = Math.round(clearanceF2*100)/100.0;
				freqClearanceF2 = (signalCenterAtHeight - clearanceF2) - elevationWithClutter ;
				freqClearanceF2 = Math.round(freqClearanceF2*100)/100.0;
				if (result.minClearanceF2 == -1|| result.minClearanceF2 > freqClearanceF2)
		    		result.minClearanceF2 = freqClearanceF2;
				F2Up = signalCenterAtHeight + clearanceF2;
				F2Dn = signalCenterAtHeight - clearanceF2;
			}

			double antennaHeight = 0;
			if (i == 0 || i == elevationData.length - 1)
				antennaHeight = (i == 0 ? fromAntennaHt : toAntennaHt);
	    	
			dist = Math.round(dist)/1000.0; // convert to km
			LineOfSightClearance clearance = new LineOfSightClearance(pathDist, dist, elevationWithClutter, fromAntennaHt, toAntennaHt);
	    	if (result.minClearance == -1 || result.minClearance > clearance.clearance) {
				result.minClearance = clearance.clearance;
			}
/*	    	
			if(collectProfileData)
			{
				result.chartData.push([dist, elevation, elevationWithClutter, clearance.losHeight, F1Up, F1Dn, F2Up, F2Dn, curvature]);
				//vR.los.losProfile.chartDataLOS.addRow([dist, elevation, elevationWithClutter, clearance.losHeight, F1Up, F1Dn, F2Up, F2Dn, curvature]);
				result.elevationPathInfoData.push(
				{ 
					location : elevationData[i].location, 
					distance : dist, 
					elevation: elevation, 
					elevationWithClutter: elevationWithClutter, 
					antennaHeight: antennaHeight, 
					losHeight: clearance.losHeight, 
					clearance: clearance.clearance, 
					includeF1: losParams.includeF1,
					clearanceF1: freqClearanceF1,
	                includeF2: losParams.includeF2,
					curvature: curvature    
				});
			}
*/	
//	        System.out.println(thisPoint.latitude.getDegrees() + " " + thisPoint.longitude.getDegrees() + " " + clearance.clearance + " " + freqClearanceF1 + " " + " " + elevationWithClutter + " " + elevationData[i].elevation);
		}
		result.status = "DONE";
		return result;
	}
}
