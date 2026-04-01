import static qupath.lib.gui.scripting.QPEx.*
import qupath.fx.dialogs.Dialogs
import qupath.fx.dialogs.FileChoosers
import java.io.File

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

// User input
def homeDir = new File(System.getProperty("user.home"))
def outputDir = FileChoosers.promptForDirectory("Choose output directory", homeDir)
if (outputDir == null) {
    println "Warning: canceled because no output directory was selected."
    return
}
if (!outputDir.exists()) outputDir.mkdirs()

def prefix = Dialogs.showInputDialog("Filename prefix", "Example: A04227A4", "line_lengths")
if (prefix == null) {
    println "Warning: canceled because no filename prefix was provided."
    return
}
prefix = prefix.trim()

boolean useA = Dialogs.showYesNoDialog(
    "Name pattern",
    "Match A1/A2/... ?\nChoose No to match 1/2/..."
)

// Match annotation names like A1/A2/... or 1/2/...
def pattern = useA ? ~/^A\d+$/ : ~/^\d+$/

// Get line annotations only, and keep only matched names
def annotations = getAnnotationObjects()
def lineAnnotations = annotations.findAll { ann ->
    def roi = ann.getROI()
    roi != null &&
    roi.getRoiName() == "Line" &&
    ann.getName() != null &&
    ann.getName().matches(pattern)
}.sort { a, b ->
    def numA = a.getName().replaceAll("\\D", "").toInteger()
    def numB = b.getName().replaceAll("\\D", "").toInteger()
    numA <=> numB
}

if (lineAnnotations.isEmpty()) {
    Dialogs.showErrorMessage(
        "Error",
        "No matching line annotations found. Pattern: ${useA ? 'A1/A2/...' : '1/2/...'}"
    )
    return
}

println "Info: found ${lineAnnotations.size()} matching line annotations."

// CSV escaping helper
def escapeCsv = { String s ->
    s = s == null ? "" : s
    if (s.contains(",") || s.contains("\"") || s.contains("\n")) {
        return "\"" + s.replace("\"", "\"\"") + "\""
    }
    return s
}

// Output file
String fileName = (prefix.isEmpty() ? "" : prefix + "_") + "line_lengths.csv"
File outputFile = new File(outputDir, fileName)

// Build CSV
def header = "Name,Length_px\n"
def csvRows = lineAnnotations.collect { ann ->
    def name = ann.getName() ?: ""
    def length = ann.getROI().getLength()
    return "${escapeCsv(name)},${length}"
}

// Write CSV
outputFile.text = header + csvRows.join("\n")

def msg = "Export finished.\nTotal line annotations: ${lineAnnotations.size()}\nPath: ${outputFile.getAbsolutePath()}"
println msg
print "Done: " + msg