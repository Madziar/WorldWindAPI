package teleworx.worldwind;

import gov.nasa.worldwind.*;
import gov.nasa.worldwind.avlist.AVKey;
import gov.nasa.worldwind.awt.WorldWindowGLCanvas;
import gov.nasa.worldwind.geom.*;
import gov.nasa.worldwind.globes.ElevationModel;
import gov.nasa.worldwind.globes.EllipsoidalGlobe;
//import gov.nasa.worldwind.terrain.BasicElevationModelFactory;
import gov.nasa.worldwind.util.Logging;
import gov.nasa.worldwindx.examples.*;

import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.util.*;

public class WorldWindAPI {
	public static void main(String[] args) throws InterruptedException {
		Configuration.setValue("gov.nasa.worldwind.avkey.OfflineMode", true);
		boolean prop2 = Configuration.getBooleanValue("gov.nasa.worldwind.avkey.OfflineMode");
//		Configuration.setValue("gov.nasa.worldwind.avkey.EarthElevationModelConfigFile", "src/config/Earth/WorldWind topo30.xml");
		Configuration.setValue("gov.nasa.worldwind.avkey.EarthElevationModelConfigFile", "src/config/Earth/WorldWind srtm15 plus topo15.xml");
		String prop3 = Configuration.getStringValue("gov.nasa.worldwind.avkey.EarthElevationModelConfigFile");

		String apiTask = args[0];//"LOS"/"LOCLOS"/"ELV"/"ELVRNG"/"ELVINT"/"LOSELV";
		String inputFileName = args[1];//"testInputFile.csv";
		double maxHopLengthKm = Double.parseDouble(args[2]);//30
		String[] freqMHzArray = args[3].split(";");
		double deltaM = Double.parseDouble(args[4]);//100
//		File tmpDir = new File(inputFileName);
//		if (tmpDir.exists())
		calculateLos(apiTask, inputFileName, maxHopLengthKm, freqMHzArray, deltaM);
//		losTest();
//		elevationTest();
	}

	public static double roundElevationM(double elevationM) {
		return Math.round(elevationM * 100) / 100.0;		
	}
	
	public static void calculateLos(String apiTask, String inputFileName, double maxHopLengthKm, String[] freqMHzArray, double deltaM) {
		//synopsis
		//1. open input data file and recreate output file
		//2. for every elevation pair (row from input file), run los and dump result into output file
		//input file format - pipe separated
		//row_index|lat1|lon1|lat2|lon2|treeline_m

		try {
	        String[] ifnSplit = inputFileName.split("\\.");
			String outputFileName = ifnSplit[0] + "_out.csv";
	        String line = "";
	        try (
	        	BufferedReader br = new BufferedReader(new FileReader(inputFileName)); 
	        	FileWriter fileWriter = new FileWriter(outputFileName)
	        ) {
//	    		WorldWindowGLCanvas worldWindCanvas = new WorldWindowGLCanvas();
//	    		worldWindCanvas.setModel(new BasicModel());
	        	ElevationModel model = EllipsoidalGlobe.makeElevationModel(AVKey.EARTH_ELEVATION_MODEL_CONFIG_FILE, "config/Earth/EarthElevations2.xml");
	        	
	    		WorldWindLOS los = new WorldWindLOS();
	        	int lineCount = 0;
	            while ((line = br.readLine()) != null) {
	            	try {
		            	lineCount++;
		                String[] rowData = line.split("\\|");
	//	                System.out.println(rowData[0] + " " + rowData[1] + " " + rowData[2] + " " + rowData[3] + " " + rowData[4]);\
		                int rowIndex = Integer.parseInt(rowData[0]);
		                double lat1 = Double.parseDouble(rowData[1]);
		                double lon1 = Double.parseDouble(rowData[2]);
		            	String row = "";
		            	if (apiTask.equals("LOSELV")) {
			                double lat2 = Double.parseDouble(rowData[3]);
			                double lon2 = Double.parseDouble(rowData[4]);
		                	long startTime = System.currentTimeMillis();
		            		ElevationData[] elevationData = getElevationData(model, los, lat1, lon1, lat2, lon2);
		                    long stopTime = System.currentTimeMillis();
		                    long elapsedTime = stopTime - startTime;
		            		if (lineCount <= 10)
		            			System.out.println("getElevationData " + elapsedTime);
	            			row = String.valueOf(rowIndex);
		            		for (int i = 0; i < elevationData.length; i++) {
			            		row += "|" + elevationData[i].location.latitude.getDegrees() + ";" + elevationData[i].location.longitude.getDegrees() 
			            				+ ";" + roundElevationM(elevationData[i].elevation);
		            		}	            		
		            	}
		            	else if (apiTask.equals("ELVINT")) {
			                double lat2 = Double.parseDouble(rowData[3]);
			                double lon2 = Double.parseDouble(rowData[4]);
		                	long startTime = System.currentTimeMillis();
		                	Position position = calculateIntermediateSiteLocation(/*worldWindCanvas*/model, lat1, lon1, lat2, lon2, 1000 * maxHopLengthKm, deltaM);
		                    long stopTime = System.currentTimeMillis();
		                    long elapsedTime = stopTime - startTime;
		            		if (lineCount <= 10)
		            			System.out.println("calculateIntermediateSiteLocation " + elapsedTime);            
		            		row = rowIndex + "|" + position.latitude.getDegrees() + "|" + position.longitude.getDegrees() + "|" + roundElevationM(position.elevation);
		                }
		                else if (apiTask.equals("ELVRNG")) {
			                double radiusKm = Double.parseDouble(rowData[3]);
		                	long startTime = System.currentTimeMillis();
		                	Position position = calculateSearchRingPositionForItem(/*worldWindCanvas*/model, lat1, lon1, 1000 * radiusKm, deltaM);
		                    long stopTime = System.currentTimeMillis();
		                    long elapsedTime = stopTime - startTime;
		            		if (lineCount <= 10)
		            			System.out.println("calculateSearchRingPositionForItem " + elapsedTime);            
		            		row = rowIndex + "|" + position.latitude.getDegrees() + "|" + position.longitude.getDegrees() + "|" + roundElevationM(position.elevation);
		                }
		                else if (apiTask.equals("ELV")) {
		                	double elevation = calculateElevationForItem(/*worldWindCanvas*/model, lat1, lon1);
		            		row = rowIndex + "|" + roundElevationM(elevation);
		                }
		                else {
			                double lat2 = Double.parseDouble(rowData[3]);
			                double lon2 = Double.parseDouble(rowData[4]);
		//	                double treeLineM = 0;//Double.parseDouble(rowData[5]);
			                int length = 1;
			                if (apiTask.equals("LOS"))
			                	length = freqMHzArray.length;
			                for (int i = 0; i < length; i++) {
			                	double freqMHz = Double.parseDouble(freqMHzArray[i]);
			                	long startTime = System.currentTimeMillis();
				            	LosResult losResult = calculateLosForItem(/*worldWindCanvas*/model, los, lat1, lon1, lat2, lon2, maxHopLengthKm, freqMHz);
			                    long stopTime = System.currentTimeMillis();
			                    long elapsedTime = stopTime - startTime;
			            		if (lineCount <= 10)
			            			System.out.println("calculateLosForItem " + elapsedTime);            
				            	if (apiTask.equals("LOS")) {
				            		if (i == 0)
				            			row = rowIndex + "|" + Math.max(0, -losResult.minClearanceF1);
				            		else
				            			row += ";" + Math.max(0, -losResult.minClearanceF1);
				            	}
				            	else if (apiTask.equals("LOCLOS")) {
				            		if (losResult.minClearanceF1LocLat == 0) {
				            			Position p1 = new Position(LatLon.fromDegrees(lat1, lon1), 0);
				            			Position p2 = new Position(LatLon.fromDegrees(lat2, lon2), 0);
				            	        Position p = Position.interpolate(0.5, p1, p2);
				            	        losResult.minClearanceF1LocLat = p.latitude.getDegrees();
				            	        losResult.minClearanceF1LocLon = p.longitude.getDegrees();
				            		}
				            		row = rowIndex + "|" + losResult.minClearanceF1LocLat + "|" + losResult.minClearanceF1LocLon;
				            	}
			                }
		                }
		            	
	            		fileWriter.append(row + "\n");
	            		if (lineCount <= 10)
	            			System.out.println(row);
	            	}
	            	catch (Exception e) {
	    	            e.printStackTrace();
	            	}
	            }
	        }
	        catch (IOException e) {
	            e.printStackTrace();
	        }
		}
		catch (Exception e) {
            e.printStackTrace();
		}
    }

	public static Position calculateIntermediateSiteLocation(/*WorldWindowGLCanvas worldWindCanvas*/ElevationModel model, double lat, double lon, double lat2, double lon2, double searchRadiusM, double deltaM) {
//       	ElevationModel model = worldWindCanvas.getModel().getGlobe().getElevationModel();
		
		LatLon center1 = LatLon.fromDegrees(lat, lon);
		LatLon center2 = LatLon.fromDegrees(lat2, lon2);
		Position p1 = new Position(center1, 0);
		Position p2 = new Position(center2, 0);
        Position p = Position.interpolate(0.5, p1, p2);
		
		LatLon center = new LatLon(p);
		Position maxPosition = new Position(center, -100000);
		LatLon currentLeftCenter = new LatLon(center);
		LatLon currentRightCenter = new LatLon(center);
		int count = 0;
		while (true) {
			double distM1 = WorldWindLOS.computeDistanceBetween(center1, currentRightCenter);
			double distM2 = WorldWindLOS.computeDistanceBetween(center2, currentLeftCenter);
			if (distM1 > searchRadiusM || distM2 > searchRadiusM)
				break;

			for (int i = 0; i < 2; i++) {
				if (i == 1 && currentLeftCenter.latitude.getDegrees() == currentRightCenter.latitude.getDegrees())
					break;
				LatLon currentCenter = (i == 0 ? currentLeftCenter : currentRightCenter);
				for (int j = 0; j < 2; j++) {
					int k = 0;
					if (j == 1)
						k = -1;
					while (true) {
						Angle azimuth = Angle.ZERO;
						if (j == 1)
							azimuth = Angle.POS180;
						
						Angle deltaRad = Angle.fromRadians(Math.abs(k) * deltaM / WorldWindLOS._equatorialRadius);
						LatLon point = LatLon.linearEndPosition(currentCenter, azimuth, deltaRad);
						
						distM1 = WorldWindLOS.computeDistanceBetween(center1, point);
						distM2 = WorldWindLOS.computeDistanceBetween(center2, point);
						if (distM1 > searchRadiusM || distM2 > searchRadiusM)
							break;
						
						double elevationM = model.getUnmappedLocalSourceElevation(point.getLatitude(), point.getLongitude());
						if (elevationM > maxPosition.elevation) {
							maxPosition = new Position(point, elevationM);
						}
						
						count++;
						
						if (j == 0)
							k++;
						else
							k--;
					}
				}
			}
			currentLeftCenter = LatLon.linearEndPosition(currentLeftCenter, Angle.NEG90, Angle.fromRadians(deltaM / WorldWindLOS._equatorialRadius));
			currentRightCenter = LatLon.linearEndPosition(currentRightCenter, Angle.POS90, Angle.fromRadians(deltaM / WorldWindLOS._equatorialRadius));
		}
		
//    	System.out.println("count = " + count);

		return maxPosition;
	}
	
	public static Position calculateSearchRingPositionForItem(/*WorldWindowGLCanvas worldWindCanvas*/ElevationModel model, double lat, double lon, double searchRadiusM, double deltaM) {
//       	ElevationModel model = worldWindCanvas.getModel().getGlobe().getElevationModel();
                
		Position maxPosition = new Position(LatLon.fromDegrees(lat, lon), -100000);
		LatLon center = LatLon.fromDegrees(lat,  lon);
		LatLon currentLeftCenter = LatLon.fromDegrees(lat, lon);
		LatLon currentRightCenter = LatLon.fromDegrees(lat, lon);
		int count = 0;
		while (true) {
			double distM = WorldWindLOS.computeDistanceBetween(center, currentLeftCenter);
			if (distM > searchRadiusM)
				break;

			for (int i = 0; i < 2; i++) {
				if (i == 1 && currentLeftCenter.latitude.getDegrees() == currentRightCenter.latitude.getDegrees())
					break;
				LatLon currentCenter = (i == 0 ? currentLeftCenter : currentRightCenter);
				for (int j = 0; j < 2; j++) {
					int k = 0;
					if (j == 1)
						k = -1;
					while (true) {
						Angle azimuth = Angle.ZERO;
						if (j == 1)
							azimuth = Angle.POS180;
						
						Angle deltaRad = Angle.fromRadians(Math.abs(k) * deltaM / WorldWindLOS._equatorialRadius);
						LatLon point = LatLon.linearEndPosition(currentCenter, azimuth, deltaRad);
						distM = WorldWindLOS.computeDistanceBetween(center, point);
						if (distM > searchRadiusM)
							break;
						
						double elevationM = model.getUnmappedLocalSourceElevation(point.getLatitude(), point.getLongitude());
						if (elevationM > maxPosition.elevation) {
							maxPosition = new Position(point, elevationM);
						}
						
						count++;
						
						if (j == 0)
							k++;
						else
							k--;
					}
				}
			}
			currentLeftCenter = LatLon.linearEndPosition(currentLeftCenter, Angle.NEG90, Angle.fromRadians(deltaM / WorldWindLOS._equatorialRadius));
			currentRightCenter = LatLon.linearEndPosition(currentRightCenter, Angle.POS90, Angle.fromRadians(deltaM / WorldWindLOS._equatorialRadius));
		}
		
//    	System.out.println("count = " + count);
		
		return maxPosition;
	}
	
	public static double calculateElevationForItem(/*WorldWindowGLCanvas worldWindCanvas*/ElevationModel model, double lat, double lon) {
//       	ElevationModel model = worldWindCanvas.getModel().getGlobe().getElevationModel();
       	LatLon latLon = LatLon.fromDegrees(lat, lon);
       	double elevation = model.getUnmappedLocalSourceElevation(latLon.getLatitude(), latLon.getLongitude());
       	if (elevation == -32768) //for corrupt data, assume sea level
       		elevation = 0;
		
		return elevation;
	}
	
	public static ElevationData[] getElevationData(ElevationModel model, WorldWindLOS los, double lat1, double lon1, double lat2, double lon2) {
		ElevationData[] elevationData = null;
		elevationData = new ElevationData[50];

		elevationData[0] = new ElevationData();
		elevationData[0].location = LatLon.fromDegrees(lat1, lon1);
		elevationData[elevationData.length - 1] = new ElevationData();
		elevationData[elevationData.length - 1].location = LatLon.fromDegrees(lat2, lon2);

		Position p1 = new Position(elevationData[0].location, 0);
		Position p2 = new Position(elevationData[elevationData.length - 1].location, 0);
		
//       	ElevationModel model = worldWindCanvas.getModel().getGlobe().getElevationModel();
        ArrayList<LatLon> latLons = new ArrayList<LatLon>();
		for (int j = 0; j < elevationData.length; j++) {
			double ratio = j / (double)(elevationData.length - 1);
	        Position p = Position.interpolate(ratio, p1, p2);
			elevationData[j] = new ElevationData();
			latLons.add(LatLon.fromDegrees(p.getLatitude().degrees, p.getLongitude().degrees));
			elevationData[j].location = latLons.get(latLons.size() - 1);
			elevationData[j].elevation = model.getUnmappedLocalSourceElevation(elevationData[j].location.getLatitude(), elevationData[j].location.getLongitude());
//			System.out.println(elevationData[j].location.getLatitude().degrees + " " + elevationData[j].location.getLongitude().degrees + " " + elevationData[j].elevation);
		}
		
		return elevationData;
	}
	
	public static LosResult calculateLosForItem(/*WorldWindowGLCanvas worldWindCanvas*/ElevationModel model, WorldWindLOS los, double lat1, double lon1, double lat2, double lon2, double maxHopLengthKm, double freqMHz) {
		double fromAntProposedHt = 0;
		double toAntProposedHt = 0;
		LosParams losParams = new LosParams();
		losParams.freqMHz = freqMHz;
//		losParams.clutterHeight = treeLineM;
		losParams.maxHopLengthKm = maxHopLengthKm;
		boolean collectProfileData = false;

		ElevationData[] elevationData = getElevationData(model, los, lat1, lon1, lat2, lon2);
		
		LosResult losResult = los.calculateLOSResult(elevationData, fromAntProposedHt, toAntProposedHt, losParams, collectProfileData);
//        System.out.println(losResult.minClearanceF1);
        
        return losResult;
	}
	
	public static void losTest() {
		ElevationData[] elevationData = null;
		double fromAntProposedHt = 0;
		double toAntProposedHt = 0;
		LosParams losParams = new LosParams();
		boolean collectProfileData = false;
		elevationData = new ElevationData[50];

		//chainFinder clearance=-903.12
		losParams.clutterHeight = 10;
		elevationData[0] = new ElevationData();
		elevationData[0].location = LatLon.fromDegrees(-12.61964, -75.21204);
		elevationData[elevationData.length - 1] = new ElevationData();
		elevationData[elevationData.length - 1].location = LatLon.fromDegrees(-12.66454, -75.32431);

		Position p1 = new Position(elevationData[0].location, 0);
		Position p2 = new Position(elevationData[elevationData.length - 1].location, 0);
		WorldWindowGLCanvas worldWindCanvas = new WorldWindowGLCanvas();
		worldWindCanvas.setModel(new BasicModel());
		
       	ElevationModel model = worldWindCanvas.getModel().getGlobe().getElevationModel();
        ArrayList<LatLon> latLons = new ArrayList<LatLon>();
		for (int j = 0; j < elevationData.length; j++) {
			double ratio = j / (double)(elevationData.length - 1);
	        Position p = Position.interpolate(ratio, p1, p2);
			elevationData[j] = new ElevationData();
			latLons.add(LatLon.fromDegrees(p.getLatitude().degrees, p.getLongitude().degrees));
			elevationData[j].location = latLons.get(latLons.size() - 1);
//			elevationData[j].elevation = model.getElevation(elevationData[j].location.getLatitude(), elevationData[j].location.getLongitude());
			elevationData[j].elevation = model.getUnmappedLocalSourceElevation(elevationData[j].location.getLatitude(), elevationData[j].location.getLongitude());
//			System.out.println(elevationData[j].location.getLatitude().degrees + " " + elevationData[j].location.getLongitude().degrees + " " + elevationData[j].elevation);
		}
//        double[] elevations = new double[elevationData.length];
//        Sector sector = Sector.fromDegrees(-13, -12, -76, -75);
//        double targetResolution = Angle.fromDegrees(1d).radians;
//        double availability = model.getLocalDataAvailability(sector, targetResolution);
//		System.out.println(availability);
//        double resolutionAchieved = model.getElevations(
//            sector, latLons, targetResolution, elevations
//        );
//        for (int j = 0; j < elevationData.length - 1; j++) {
//        	elevationData[j].elevation = elevations[j];
//			System.out.println(elevationData[j].location.getLatitude().degrees + " " + elevationData[j].location.getLongitude().degrees + " " + elevationData[j].elevation);
//        }
//		System.out.println(targetResolution + " " + resolutionAchieved);
		
		WorldWindLOS los = new WorldWindLOS();
		LosResult losResult = los.calculateLOSResult(elevationData, fromAntProposedHt, toAntProposedHt, losParams, collectProfileData);
        System.out.println(losResult.minClearanceF1);
	}
	
	public static void elevationTest() {
    	long startTime = System.currentTimeMillis();
    	
        ArrayList<LatLon> latlons = new ArrayList<LatLon>();

        int max = 10;
        for (int i = 0; i < max; i++) {
            latlons.add(LatLon.fromDegrees(45.50d + i * 1, -123.3d));
        }

        Position p1 = new Position(latlons.get(0), 0);
        Position p2 = new Position(latlons.get(max - 1), 0);
        Position p = Position.interpolate(0.5, p1, p2);
        System.out.println(p1.latitude + "," + p1.longitude);
        System.out.println(p2.latitude + "," + p2.longitude);
        System.out.println(p.latitude + "," + p.longitude);
        
		WorldWindowGLCanvas worldWindCanvas = new WorldWindowGLCanvas();
		worldWindCanvas.setModel(new BasicModel());
        ElevationModel model = worldWindCanvas.getModel().getGlobe().getElevationModel();
        StringBuffer sb = new StringBuffer();
        for (LatLon ll : latlons)
        {
            double e = model.getElevation(ll.getLatitude(), ll.getLongitude());
            sb.append("\n").append(e);
        }

        Logging.logger().info(sb.toString());
        
        long stopTime = System.currentTimeMillis();
        long elapsedTime = stopTime - startTime;
        System.out.println("elevation " + elapsedTime + "ms [" + latlons.size() + "]");            
	}
	
	public static int junitTest1() {
		return 1;
	}
}
