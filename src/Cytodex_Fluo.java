//  This plugin extract branches of spherois and compute lengths branching ...


import Skeletonize3D_.Skeletonize3D_;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
//import ij.Macro;
import ij.WindowManager;
import ij.gui.*;
import ij.io.FileSaver;
import ij.io.Opener;
import ij.measure.Calibration;
import ij.measure.Measurements;
import ij.measure.ResultsTable;
import ij.plugin.Duplicator;
import ij.plugin.PlugIn;
import ij.plugin.filter.Analyzer;
import ij.plugin.filter.BackgroundSubtracter;
import ij.plugin.filter.GaussianBlur;
import ij.plugin.filter.MaximumFinder;
import ij.plugin.filter.RankFilters;
import ij.plugin.frame.RoiManager;
import ij.process.AutoThresholder;
import ij.process.ImageConverter;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import ij.util.ArrayUtil;
import java.awt.Color;
import java.awt.Polygon;
import java.io.*;
import static java.lang.Math.sqrt;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import skeleton_analysis.AnalyzeSkeleton_;
import skeleton_analysis.Edge;
import skeleton_analysis.Point;
import skeleton_analysis.SkeletonResult;


public class Cytodex_Fluo implements PlugIn {
	
    private final boolean canceled =false;
    private static int nbNucleus;   // total number of nucleus
    private static double sphXcent, sphYcent, sphFerret;   // Centroid coord and ferret of spheroid
    private static final ResultsTable table = new ResultsTable();
    private static boolean isStack;
    private static Roi spheroidRoi;
    private static String imgOutDir;
    private static String fileNameWithOutExt;
    
// clear outside selection with true background
    public void clearOutside(ImagePlus img, Roi cropRoi) {
        ImageProcessor ip = img.getProcessor();
        ip.setColor(Color.BLACK);
        img.setActivated();
        if (isStack) {
            for(int z = 1; z < img.getNSlices(); z++) {
                ip.fillOutside(cropRoi);
            }
        }
        else {
            ip.fillOutside(cropRoi);
        }
        img.deleteRoi();
    }
    
// find centroid and ferret diameter of the spheroid    
    public void getCentroid(ImagePlus img) {
        
        if (isStack) img.setSlice(2);
        Analyzer measure = new Analyzer(img,Measurements.CENTROID + Measurements.FERET,table);
        measure.measure();
	sphXcent = table.getValue("X",0);
	sphYcent = table.getValue("Y",0);
	sphFerret = table.getValue("Feret",0);
	table.reset();
	img.getProcessor().setColor(Color.BLACK);
        if (isStack) {  // do not fill DAPI channel do after Band pass filter in findNucleus
            for (int s = 2; s <= img.getNSlices(); s++) {
                img.setSlice(s);
                img.getProcessor().fill(img.getRoi());
            }
        }
        else {
            img.getProcessor().fill(img.getRoi());
            img.updateImage();
        }
	img.deleteRoi();
    }

    

/**
     * Returns the location of pixels clockwise along circumference
     * using Bresenham's Circle Algorithm
     * keep point only if pixel is inside image
     */
    public static ArrayList<Point> BresenhamCircle(int xc,int yc,int r) {    
        ArrayList<Point> ret = new ArrayList<Point>();
        int x,y,p;
        x=0;
        y=r;
        ret.add(new Point(xc+x,yc-y,0));
        p=3-(2*r);
        for(x=0;x<=y;x++) {
            if (p<0) {
                p=(p+(4*x)+6);
            }
            else {
                y=y-1;
                p=p+((4*(x-y)+10));
            }
            ret.add(new Point(xc+x,yc-y,0));
            ret.add(new Point(xc-x,yc-y,0));
            ret.add(new Point(xc+x,yc+y,0));
            ret.add(new Point(xc-x,yc+y,0));
            ret.add(new Point(xc+y,yc-x,0));
            ret.add(new Point(xc-y,yc-x,0));
            ret.add(new Point(xc+y,yc+x,0));
            ret.add(new Point(xc-y,yc+x,0));
        }
        return ret;
}


// compute local thickness and branches diameter
    public static double[] localThickness (ImagePlus imgMask, ImagePlus imgSkel) {
        imgMask.show();
        Calibration cal = new Calibration();
        cal = imgMask.getCalibration();
        double vxWH = cal.pixelWidth;
// parameters for center, min, max and step radius 
        final double incStep = 50;              // every 50 microns
        double cx = sphXcent;
        double cy = sphYcent;
        double startRadius = (sphFerret/2) + 1 ;
        final int wdth = imgMask.getWidth();
	final int hght = imgMask.getHeight();
	final double dx, dy, dz, maxEndRadius, stepRadius;      
        double [] radii;
        ArrayList<Double> diameters = new ArrayList<Double>();
        double [] meanDiameters;

        WindowManager.setCurrentWindow(imgMask.getWindow());
        IJ.run("Local Thickness (complete process)", "threshold=255");
// wait for end of local thickness process
        while (WindowManager.getImage("Branches_Mask_LocThk") == null) {
            IJ.wait(100);
        }
        ImagePlus imgLocalThickness = WindowManager.getImage("Branches_Mask_LocThk");
        
        imgLocalThickness.setCalibration(cal);
        dx = (cx<=wdth/2) ? (cx-wdth) : cx;
        dy = (cy<=hght/2) ? (cy-hght) : cy;
        maxEndRadius = Math.sqrt(dx*dx + dy*dy);
        stepRadius = incStep;
// save image map
        FileSaver imgLoth = new FileSaver(imgLocalThickness);
        imgLoth.saveAsTiff(imgOutDir+fileNameWithOutExt+"_Map.tif");
// Calculate how many samples will be taken
        final int size = (int) ((maxEndRadius-startRadius)/stepRadius)+1;
        
// Create arrays for radii (in physical units) and intersection counts
        radii = new double[size];
        meanDiameters = new double[size];
        for (int i = 0; i < size; i++) {
                radii[i] = startRadius + i*stepRadius;
        }
 
        ArrayList<Point> ptCircle = new ArrayList<Point>();
       
// compute points on circle every step microns and store in array
        for (int r = 0; r < radii.length; r++) {
            diameters.clear();
            ptCircle = BresenhamCircle((int)cal.getRawX(cx), (int)cal.getRawY(cy), (int)(radii[r]/cal.pixelWidth));            
            Polygon polygon = new Polygon();
            for (int i = 0; i < ptCircle.size(); i++) {
                polygon.addPoint(ptCircle.get(i).x, ptCircle.get(i).y);
            }
            
            for (int i = 0; i < ptCircle.size(); i++) {
                // check if pixel inside image
                
                int x = ptCircle.get(i).x;
                int y = ptCircle.get(i).y;
                if (((x > 0) || (x < wdth)) && ((y > 0) || (y < hght))) {
                    //inside skeleton branch
                    if (imgSkel.getProcessor().getPixelValue(x, y) == 255) { 
                        // read pixel value in localthickness image
                        double pixelValue = imgLocalThickness.getProcessor().getPixelValue(x, y);
                        diameters.add(pixelValue); 
                     }      
                 }                  
            }
            // calculate the mean diameters
            double avg = 0;
            for (double sum:diameters) {
                avg += sum;
            }
            meanDiameters[r] = avg/diameters.size();
        }
        imgLocalThickness.changes = false;
        imgLocalThickness.close();
        imgSkel.changes = false;
        imgSkel.close();
        imgSkel.flush();
        imgMask.changes = false;
        imgMask.close();
        imgMask.flush();
        return(meanDiameters);    
    }
   
/* count nucleus take the coordinates of pixel in nucleus image
/* if grey value = 0 in branchs mask then nucleus ++
*/

    public Polygon findNucleus(ImagePlus imgNuc, ImagePlus imgBranchs, int spheroid) {
	nbNucleus = 0;
        int xCoor, yCoor;
        Duplicator imgDup = new Duplicator();
        ImagePlus imgNucleus = imgDup.run(imgNuc,1,1);
        ImagePlus colorNucleus = imgNuc.duplicate();
        ImageConverter imgConv = new ImageConverter(colorNucleus);
        imgConv.convertToRGB();
// run difference of Gaussians
        IJ.run(imgNucleus,"Difference of Gaussians", "  sigma1=10 sigma2=2 enhance slice");
        imgNucleus.getProcessor().setColor(Color.BLACK);
        imgNucleus.setRoi(spheroidRoi);
        imgNucleus.getProcessor().fill(imgNucleus.getRoi());
        imgNucleus.deleteRoi();
// find maxima
        ImageProcessor ipNucleus = imgNucleus.getProcessor();
        MaximumFinder findMax = new MaximumFinder();
        Polygon  nucleusPoly = findMax.getMaxima(ipNucleus, 10,false);        
        colorNucleus.setColor(Color.BLUE);
        for (int i = 0; i < nucleusPoly.npoints; i++) {
		xCoor = (int)nucleusPoly.xpoints[i];
		yCoor = (int)nucleusPoly.ypoints[i];
		if ((int)imgBranchs.getProcessor().getPixelValue(xCoor,yCoor) == 255) {
                    nbNucleus++;
                    OvalRoi ptNucleus = new OvalRoi(xCoor, yCoor,4,4);
                    colorNucleus.getProcessor().draw(ptNucleus);
                    colorNucleus.updateAndDraw();                  
                }
	}
        FileSaver imgSave = new FileSaver(colorNucleus);
        imgSave.saveAsTiff(imgOutDir+fileNameWithOutExt+"_Crop_"+spheroid+"_nucleus.tif");
        colorNucleus.close();
        imgNucleus.changes = false;
        imgNucleus.close();
        imgNucleus.flush();
        spheroidRoi = null;
        return nucleusPoly;
    }
    
// Substract background
    public void backgroundSubstract(ImagePlus img) {
        BackgroundSubtracter imgSubstract = new BackgroundSubtracter();
        if (isStack) {
            for (int s = 1;s <= img.getNSlices(); s++) { 
                img.setSlice(s);
                imgSubstract.rollingBallBackground(img.getProcessor(), 50, false, false, false, false, false);
            }
            img.updateAndRepaintWindow();
        }
        else imgSubstract.rollingBallBackground(img.getProcessor(), 50, false, false, false, false, false);
    }
    
    public double calculateDistance(ImagePlus img, Point point1, Point point2) {
        return Math.sqrt(  Math.pow( (point1.x - point2.x) * img.getCalibration().pixelWidth, 2) 
                          + Math.pow( (point1.y - point2.y) * img.getCalibration().pixelHeight, 2)
                          + Math.pow( (point1.z - point2.z) * img.getCalibration().pixelDepth, 2));
    }

    // Calculate lenght of branches after skeletonize
    public void analyzeSkel (ImagePlus img, int spheroid, BufferedWriter output) {
	int nbSkeleton;             // number of skeleton
        double totalLength = 0;     // total branch lenght/spheroid
        int totalBranches = 0;   // total number of branches/spheroid
        double euclideanDist = 0;   // euclidean distance
        int nbTubes = 0;            // number of tubes/spheroid = number of skeleton + number of junction
        double meanTortuosity;      // mean euclidean distance / branch length
        double sdTortuosite;        // sd euclidean distance / branch length
        int nbJunctions = 0;         // number of junctions / spheroid
        
        AnalyzeSkeleton_ analyzeSkeleton = new AnalyzeSkeleton_();
        AnalyzeSkeleton_.calculateShortestPath = true;
        analyzeSkeleton.setup("",img);
        SkeletonResult skeletonResults = analyzeSkeleton.run(AnalyzeSkeleton_.NONE,false,true,null,true,false);
        nbSkeleton = skeletonResults.getNumOfTrees();
        int[] branchNumbers = skeletonResults.getBranches();
        int[] junctionNumbers = skeletonResults.getJunctions();
        // save labelled skeletons
        ImageStack labSkel = analyzeSkeleton.getLabeledSkeletons();
        ImagePlus imgLab = new ImagePlus("Labelled skeleton",labSkel);
        IJ.run(imgLab,"Fire","");
        FileSaver imgLab_save = new FileSaver(imgLab);
        imgLab_save.saveAsTiff(imgOutDir+fileNameWithOutExt+"_labSkel.tif");
        for (int b = 0; b < branchNumbers.length; b++) { 
                totalBranches += branchNumbers[b];
                nbJunctions += junctionNumbers[b];
        }
        float[] tortuosity = new float[totalBranches];
        for (int i = 0; i < nbSkeleton; i++) {
            ArrayList<Edge> listEdges;
            listEdges = skeletonResults.getGraph()[i].getEdges();
            for (int e = 0; e < listEdges.size(); e++) {
                totalLength += listEdges.get(e).getLength();
                euclideanDist = calculateDistance(img, listEdges.get(e).getV1().getPoints().get(0), listEdges.get(e).getV2().getPoints().get(0));
                tortuosity[e] = (float) (euclideanDist/listEdges.get(e).getLength());
            }
        }
        ij.util.ArrayUtil stats = new ArrayUtil(tortuosity);
        meanTortuosity = stats.getMean();
        sdTortuosite = sqrt(stats.getVariance());
        nbTubes = skeletonResults.getNumOfTrees() + nbJunctions;
       
        try {
            // write data
            output.write(fileNameWithOutExt + "\t" + spheroid + "\t" + nbSkeleton + "\t" + totalBranches + "\t" + totalLength +
                    "\t" + nbJunctions + "\t" + nbTubes + "\t" + meanTortuosity + "\t" + 
                    sdTortuosite + "\t" + nbNucleus + "\n");
            output.flush();
        } catch (IOException ex) {
            Logger.getLogger(Cytodex_Fluo.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    
    
    
    public void run(String arg) {
        try {
            if (canceled) {
                IJ.showMessage(" Pluging canceled");
                return;
            }           
            String imageDir = IJ.getDirectory("Choose Directory Containing TIF Files...");
            if (imageDir == null) return;
            File inDir = new File(imageDir);
            String [] imageFile = inDir.list();
            if (imageFile == null) return;
// create directory to store images
            imgOutDir = imageDir+"Images/";
            File imgTmpDir = new File(imgOutDir);
            if (!imgTmpDir.isDirectory())
                imgTmpDir.mkdir();
// write headers for global results
            FileWriter fwAnalyze;
            fwAnalyze = new FileWriter(imageDir + "Analyze_skeleton_results.xls",false);
            BufferedWriter outputAnalyze = new BufferedWriter(fwAnalyze);
            outputAnalyze.write("Image\t#Spheroids\t#Skeletons\t#Branches\tTotal branch length\t#Junctions\t"
                    + "#Tubes\tMean Tortuosity\tSD Tortuosity\t#Cells\n");
            outputAnalyze.flush();
// write headers for Sholl diameters results
            FileWriter fwDiameter;
            fwDiameter = new FileWriter(imageDir + "Analyze_skeleton_diameters_results.xls",false);
            BufferedWriter outputDiameter = new BufferedWriter(fwDiameter);
            outputDiameter.write("Image\t#Spheroids\tMean Diameter(d0)\tMean Diameter(d+50)\tMean Diameter(d+100)"
                    + "\tMean Diameter(d+150)\tMean Diameter(d+200)\tMean Diameter(d+250)\tMean Diameter(d+300)"
                    +"\tMean Diameter(d+350)\tMean Diameter(d+400)\tMean Diameter(d+450)\tMean Diameter(d+500)\n");
            outputDiameter.flush();

            Duplicator imgDup = new Duplicator();
            for (int i = 0; i < imageFile.length; i++) {
                if (imageFile[i].endsWith(".tif")) {
                    String imagePath = imageDir + imageFile[i];
                    Opener imgOpener = new Opener();
                    ImagePlus imgOrg = imgOpener.openImage(imagePath);
                    if (imgOrg.getNSlices() > 1) {
                        isStack = true;
                    }
                    fileNameWithOutExt = imageFile[i].substring(0, imageFile[i].length() - 4);

//convert to 8 bytes
//                    if (isStack) {
//                        new StackConverter(imgOrg).convertToGray8();
//                    }
//                    else {
//                        new ImageConverter(imgOrg).convertToGray8();
//                    }
                    
                    imgOrg.show();
                    IJ.run("Enhance Contrast", "saturated=0.35");
                    if (RoiManager.getInstance() != null) RoiManager.getInstance().close();
                    RoiManager rm = new RoiManager();
                    new WaitForUserDialog("Select part(s) of image to analyze\nPress t to add selection to the ROI manager.").show();
                    
                    for (int r = 0; r < rm.getCount(); r++) {
                        WindowManager.setTempCurrentImage(imgOrg);
                        rm.select(r);
                        int [] roiType = new int[rm.getCount()];
                        roiType[r] = imgOrg.getRoi().getType();
                        IJ.run("Duplicate...", "title=Crop duplicate range=1-" + imgOrg.getNSlices());
                        ImagePlus imgCrop = WindowManager.getCurrentImage();
                        imgCrop.show();
//// substract background                    
//                        backgroundSubstract(imgCrop);

                        
// fill outside roi with background except if ROI is a rectangle
                        if (roiType[r] != 0) clearOutside(imgCrop, imgCrop.getRoi());
// save croped image
                        FileSaver imgCrop_save = new FileSaver(imgCrop);
                        if (isStack) 
                            imgCrop_save.saveAsTiffStack(imgOutDir+fileNameWithOutExt+"_Crop"+r+".tif");
                        else 
                            imgCrop_save.saveAsTiff(imgOutDir+fileNameWithOutExt+"_Crop"+r+".tif");
                        if (isStack) 
                            imgCrop.setSlice(2);
                        IJ.setTool("oval");
                        new WaitForUserDialog("Outline the spheroid").show(); 			// ask for remove spheroid
                        spheroidRoi = imgCrop.getRoi();
// remove spheroid			
                        getCentroid(imgCrop);


// create image for branch mask
                        ImagePlus imgBranchsMask = imgDup.run(imgCrop,imgCrop.getNSlices(),imgCrop.getNSlices());

// threshold image                        
                        ImageProcessor ipBranchsMask = imgBranchsMask.getProcessor();
                        RankFilters median = new RankFilters();
                        GaussianBlur blur = new GaussianBlur();
                        if (isStack) {    
                            for (int s = 1; s <= imgBranchsMask.getNSlices(); s++) {
                                imgBranchsMask.setSlice(s);median.rank(ipBranchsMask, 2, RankFilters.MEDIAN);
                                blur.blurGaussian(ipBranchsMask, 1.5, 1.5, 1);
                                ipBranchsMask.setAutoThreshold(AutoThresholder.Method.Triangle, true, 1);
                                double minThreshold = ipBranchsMask.getMinThreshold();
                                double maxThreshold = ipBranchsMask.getMaxThreshold();
                                ipBranchsMask.setThreshold(minThreshold,maxThreshold,1);
                                imgBranchsMask.updateAndDraw();
                                WindowManager.setTempCurrentImage(imgBranchsMask);
                                IJ.run("Convert to Mask");
                            }   
                        }
                        else {
                            median.rank(ipBranchsMask, 2, RankFilters.MEDIAN);
                            blur.blurGaussian(ipBranchsMask, 1.5, 1.5, 1);
                            ipBranchsMask.setAutoThreshold(AutoThresholder.Method.Triangle, true, 1);
                            double minThreshold = ipBranchsMask.getMinThreshold();
                            double maxThreshold = ipBranchsMask.getMaxThreshold();
                            ipBranchsMask.setThreshold(minThreshold,maxThreshold,1);
                            imgBranchsMask.updateAndDraw();
                            WindowManager.setTempCurrentImage(imgBranchsMask);
                            IJ.run("Convert to Mask");
                        }
                        // Check if no branches
                        ImageStatistics stats = ImageStatistics.getStatistics(ipBranchsMask, ImageStatistics.MIN_MAX,imgBranchsMask.getCalibration());
                        
                        if (stats.max == 0) { // no branches
// write skeleton data with zero
                            outputAnalyze.write(fileNameWithOutExt + "\t" + (r+1) + "\t0\t0\t0\t0\t0\t0\t0\t0\n");
                            outputAnalyze.flush();                           
// write data in diameter file with zero
                            outputDiameter.write(fileNameWithOutExt + "\t" + (r+1) + "\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n");
                            outputDiameter.flush();
                            imgBranchsMask.changes = false;
                            imgBranchsMask.close();
                            imgBranchsMask.flush();
                        }
                        else {                        
                            imgBranchsMask.setTitle("Branches_Mask");
                            imgBranchsMask.setSlice(1);
                            imgBranchsMask.show();
                            new WaitForUserDialog("Correct the skeleton mask with paint tools").show();
    // Save branches mask
                            FileSaver imgMask_save = new FileSaver(imgBranchsMask);
                            if (imgBranchsMask.getNSlices() > 1) imgMask_save.saveAsTiffStack(imgOutDir+fileNameWithOutExt+"_Crop"+r+"_Mask.tif");
                            else imgMask_save.saveAsTiff(imgOutDir+fileNameWithOutExt+"_Crop"+r+"_Mask.tif");
    // find nucleus/spheroid                        
                            if (isStack) {
                                Polygon nucleusCoord = new Polygon();
                                nucleusCoord = findNucleus(imgCrop, imgBranchsMask, r);
                            }
                            imgCrop.close();
                            imgCrop.flush();
    // skeletonize
                            ImagePlus imgSkel = imgDup.run(imgBranchsMask,1,imgBranchsMask.getNSlices());
                            imgSkel.setTitle(fileNameWithOutExt+"_Crop"+ r + "_Skel.tif");
                            Skeletonize3D_ skeleton = new Skeletonize3D_();
                            //BinaryProcessor skeleton = new BinaryProcessor(imgSkel.getProcessor().convertToByteProcessor());
                            if (isStack) {
                                for (int s = 1;s <= imgSkel.getNSlices(); s++) {
                                    imgSkel.setSlice(s);
                                    skeleton.setup("",imgSkel);
                                    skeleton.run(imgSkel.getProcessor());
                                    //skeleton.skeletonize();
                                    imgSkel.updateAndDraw();
                                }
                            }
                            else {
                                skeleton.setup("",imgSkel);
                                skeleton.run(imgSkel.getProcessor());
                                //skeleton.skeletonize();
                                imgSkel.updateAndDraw();
                            }
                            imgSkel.show();
                            new WaitForUserDialog(" Correct skeleton with paint tools").show();
    // Save corrected skeleton image                        
                            FileSaver imgSkel_save = new FileSaver(imgSkel);
                            if (imgSkel.getNSlices() > 1) imgSkel_save.saveAsTiffStack(imgOutDir+fileNameWithOutExt+"_Crop"+r+"_Skel.tif");
                            else imgSkel_save.saveAsTiff(imgOutDir+fileNameWithOutExt+"_Crop"+r+"_Skel.tif");

                            analyzeSkel(imgSkel,r+1,outputAnalyze);  

    // compute mean branch diameter
                            double[] meanDiameter = localThickness(imgBranchsMask, imgSkel);

                                                 
    // write data in diameter file
                            outputDiameter.write(fileNameWithOutExt + "\t" + (r+1) + "\t");
                            for (int d = 0; d < meanDiameter.length; d++) {
                                outputDiameter.write(meanDiameter[d] + "\t");
                            }
                            outputDiameter.write("\n");
                        }
                    }
                    rm.close();                    
                    imgOrg.close();
                    imgOrg.flush();
                }
            }
           outputAnalyze.close();
//           outputCells.close();
           outputDiameter.close();
            IJ.showStatus("End of process");
        } catch (IOException ex) {
            Logger.getLogger(Cytodex_Fluo.class.getName()).log(Level.SEVERE, null, ex);
        }
    } 
    
}
