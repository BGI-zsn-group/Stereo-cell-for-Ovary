/*
 * QuPath Groovy script for exporting fixed 200×200 ROI-centered TIFF crops.
 * Applies viewer display transforms and manual gamma correction.
 */

import static qupath.lib.gui.scripting.QPEx.*
import qupath.lib.gui.dialogs.Dialogs
import qupath.lib.gui.viewer.QuPathViewer
import qupath.lib.regions.RegionRequest
import qupath.lib.roi.ROIs
import qupath.lib.awt.common.BufferedImageTools
import java.awt.image.BufferedImage
import java.awt.Graphics2D
import javax.imageio.ImageIO
import java.io.File
import qupath.lib.display.ImageDisplay
import java.awt.image.ByteLookupTable
import java.awt.image.LookupOp

// Check active viewer/image

def viewer = getCurrentViewer()
if (viewer == null) {
    Dialogs.showErrorMessage("Error", "No active image viewer. Please open an image first.")
    return
}

def imageData = getCurrentImageData()
if (imageData == null) {
    Dialogs.showErrorMessage("Error", "No image data loaded.")
    return
}

def server = imageData.getServer()
def hierarchy = imageData.getHierarchy()
int totalChannels = server.getMetadata().getChannels().size()

// User input

def outputDir = Dialogs.promptForDirectory(new File(System.getProperty("user.home")))
if (outputDir == null) {
    println "Warning: canceled because no output directory was selected."
    return
}
if (!outputDir.exists()) outputDir.mkdirs()
def outPath = outputDir.getAbsolutePath()

def prefix = Dialogs.showInputDialog("Filename prefix", "Example: Cell", "Cell")
prefix = prefix == null ? "" : prefix.trim()

def zText = Dialogs.showInputDialog("Z-slice index", "Use 0 for a single Z slice", "0")
if (zText == null) return
int zIndex
try {
    zIndex = Integer.parseInt(zText.trim())
    if (zIndex < 0 || zIndex >= server.nZSlices()) {
        Dialogs.showErrorMessage("Error", "Z-slice must be in the range 0 to ${server.nZSlices()-1}.")
        return
    }
} catch (Exception e) {
    Dialogs.showErrorMessage("Error", "Z-slice must be an integer.")
    return
}

boolean useA = Dialogs.showYesNoDialog("Name pattern", "Match A1/A2/... ?\nChoose No to match 1/2/...")

int size = 200
int halfSize = size / 2
int imgW = server.getWidth()
int imgH = server.getHeight()

def pattern = useA ? ~/^A\d+$/ : ~/^\d+$/
def anns = hierarchy.getAnnotationObjects()
        .findAll { it.getName() != null && it.getName().matches(pattern) }
        .sort { a, b ->
            def numA = a.getName().replaceAll("\\D", "").toInteger()
            def numB = b.getName().replaceAll("\\D", "").toInteger()
            numA <=> numB
        }

if (anns.isEmpty()) {
    Dialogs.showErrorMessage("Error", "No matching annotations found. Pattern: ${useA ? 'A1/A2/...' : '1/2/...'}")
    return
}
println "Info: found ${anns.size()} matching annotations."

// Replace ROI with centered fixed-size square

anns.each { ann ->
    def roi = ann.getROI()
    if (roi == null) return

    double cx = roi.getCentroidX()
    double cy = roi.getCentroidY()
    double x = Math.round(cx - halfSize)
    double y = Math.round(cy - halfSize)

    x = x < 0 ? 0 : (x + size > imgW ? imgW - size : x)
    y = y < 0 ? 0 : (y + size > imgH ? imgH - size : y)

    def plane = roi.getImagePlane()
    def newRoi = ROIs.createRectangleROI(x, y, size, size, plane)
    ann.setROI(newRoi)
}
hierarchy.fireHierarchyChangedEvent(this)
println "Info: all ROIs were replaced with 200×200 pixel squares."

int success = 0
anns.each { ann ->
    def name = ann.getName() ?: "Unknown"
    def roi = ann.getROI()
    if (roi == null) return

    int x = (int)Math.round(roi.getBoundsX())
    int y = (int)Math.round(roi.getBoundsY())
    int w = size
    int h = size

    def req = RegionRequest.createInstance(
        server.getPath(), 1.0, x, y, w, h, zIndex, 0
    )

    try {
        BufferedImage img = server.readBufferedImage(req)

        // Convert to ARGB before display transforms
        BufferedImage imgStandard = new BufferedImage(img.getWidth(), img.getHeight(), BufferedImage.TYPE_INT_ARGB)
        Graphics2D g2d = imgStandard.createGraphics()
        g2d.drawImage(img, 0, 0, null)
        g2d.dispose()

        def imageDisplay = viewer.getImageDisplay()
        BufferedImage adjusted = imageDisplay.applyTransforms(imgStandard, null)

        // Apply gamma manually and preserve alpha if present
        double gamma = viewer.getGamma()
        if (gamma != 1.0) {
            int numBands = adjusted.getRaster().getNumBands()
            byte[][] luts = new byte[numBands][256]
            for (int b = 0; b < Math.min(3, numBands); b++) {
                for (int i = 0; i < 256; i++) {
                    double val = Math.pow(i / 255.0, gamma)
                    luts[b][i] = (byte) (val * 255.0 + 0.5)
                }
            }
            if (numBands == 4) {
                for (int i = 0; i < 256; i++) {
                    luts[3][i] = (byte) i
                }
            }
            def lookupTable = new ByteLookupTable(0, luts)
            def lookupOp = new LookupOp(lookupTable, null)
            adjusted = lookupOp.filter(adjusted, null)
        }

        String fileName = (prefix.isEmpty() ? "" : prefix + "_") + name + ".tif"
        File file = new File(outPath, fileName)
        ImageIO.write(adjusted, "TIFF", file)

        success++
        println "Exported: ${file.getAbsolutePath()}"
    } catch (Exception e) {
        println "Export failed [${name}]: ${e.getMessage()}"
    }
}

// Final summary

def msg = "Export finished.\nTotal annotations: ${anns.size()} | Success: ${success}\nPath: ${outPath}"
println msg
Dialogs.showInfoNotification("Done", msg)
