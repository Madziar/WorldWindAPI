package teleworx.worldwind;

public class LosParams {
    public double minAntennaHt = 0;
    public double maxAntennaHt = 0; 
    public double antennaHeightStep = 2; 
    public double defaultAntennaHeight = 30; 
    public double clutterHeight = 0; 
    public boolean includeF1 = true; 
    public double fresnel1 = 1;
    public boolean includeF2 =  true; 
    public double fresnel2 = 0.6;
    public double kFactor = 1.33; 
    public double freqMHz = 900;
    public double numSamples = 0;
    
    public double maxHopLengthKm = -1;
}
