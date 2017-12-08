package teleworx.worldwind;

import gov.nasa.worldwind.*;
import gov.nasa.worldwind.awt.WorldWindowGLCanvas;
import gov.nasa.worldwind.geom.*;
import gov.nasa.worldwind.globes.ElevationModel;
import gov.nasa.worldwind.terrain.BasicElevationModelFactory;
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

		String inputFileName = args[0];//"testInputFile.csv";
//		File tmpDir = new File(inputFileName);
//		if (tmpDir.exists())	
		calculateLos(inputFileName);
//		losTest();
//		elevationTest();
	}

	public static void calculateLos(String inputFileName) {
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
	    		WorldWindowGLCanvas worldWindCanvas = new WorldWindowGLCanvas();
	    		worldWindCanvas.setModel(new BasicModel());
	    		WorldWindLOS los = new WorldWindLOS();
	        	
	            while ((line = br.readLine()) != null) {
	                String[] rowData = line.split("\\|");
//	                System.out.println(rowData[0] + " " + rowData[1] + " " + rowData[2] + " " + rowData[3] + " " + rowData[4]);\
	                int rowIndex = Integer.parseInt(rowData[0]);
	                double lat1 = Double.parseDouble(rowData[1]);
	                double lon1 = Double.parseDouble(rowData[2]);
	                double lat2 = Double.parseDouble(rowData[3]);
	                double lon2 = Double.parseDouble(rowData[4]);
	                double treeLineM = Double.parseDouble(rowData[5]);
	            	LosResult losResult = calculateLosForItem(worldWindCanvas, los, lat1, lon1, lat2, lon2, treeLineM);
	                fileWriter.append(rowIndex + "|" + losResult.minClearanceF1 + "\n");
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
	
	public static LosResult calculateLosForItem(WorldWindowGLCanvas worldWindCanvas, WorldWindLOS los, double lat1, double lon1, double lat2, double lon2, double treeLineM) {
		ElevationData[] elevationData = null;
		double fromAntProposedHt = 0;
		double toAntProposedHt = 0;
		LosParams losParams = new LosParams();
		boolean collectProfileData = false;
		elevationData = new ElevationData[50];

		losParams.clutterHeight = treeLineM;
		elevationData[0] = new ElevationData();
		elevationData[0].location = LatLon.fromDegrees(lat1, lon1);
		elevationData[elevationData.length - 1] = new ElevationData();
		elevationData[elevationData.length - 1].location = LatLon.fromDegrees(lat2, lon2);

		Position p1 = new Position(elevationData[0].location, 0);
		Position p2 = new Position(elevationData[elevationData.length - 1].location, 0);
		
       	ElevationModel model = worldWindCanvas.getModel().getGlobe().getElevationModel();
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
		LosResult losResult = los.calculateLOSResult(elevationData, fromAntProposedHt, toAntProposedHt, losParams, collectProfileData);
        System.out.println(losResult.minClearanceF1);
        
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
	
	
}
