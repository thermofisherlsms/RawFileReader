from __future__ import print_function

import clr

clr.AddReference('ThermoFisher.CommonCore.Data')
clr.AddReference('ThermoFisher.CommonCore.RawFileReader')
clr.AddReference('ThermoFisher.CommonCore.BackgroundSubtraction')
clr.AddReference('ThermoFisher.CommonCore.MassPrecisionEstimator')

from System import *
from System.Collections.Generic import *

from ThermoFisher.CommonCore.Data import ToleranceUnits
from ThermoFisher.CommonCore.Data import Extensions
from ThermoFisher.CommonCore.Data.Business import ChromatogramSignal, ChromatogramTraceSettings, DataUnits, Device, GenericDataTypes, SampleType, Scan, TraceType
from ThermoFisher.CommonCore.Data.FilterEnums import IonizationModeType, MSOrderType
from ThermoFisher.CommonCore.Data.Interfaces import IChromatogramSettings, IScanEventBase, IScanFilter, RawFileClassification
from ThermoFisher.CommonCore.MassPrecisionEstimator import PrecisionEstimate
from ThermoFisher.CommonCore.RawFileReader import RawFileReaderAdapter

import sys, os


if len(sys.argv) > 1:
    filename = sys.argv[1]
else:
    print('No RAW file specified!')
    quit()

# Check to see if the specified RAW file exists
if not os.path.isfile(filename):
    print('The file doesn\'t exist in the specified location - {}'.format(filename))
    quit()

# Create the IRawDataPlus object for accessing the RAW file
rawFile = RawFileReaderAdapter.FileFactory(filename)

if not rawFile.IsOpen or rawFile.IsError:
    print('Unable to access the RAW file using the RawFileReader class!')
    quit()

# Check for any errors in the RAW file
if rawFile.IsError:
    print('Error opening ({}) - {}'.format(rawFile.FileError, filename))
    quit()

# Check if the RAW file is being acquired
if rawFile.InAcquisition:
    print('RAW file still being acquired - {}'.format(filename))
    quit()

# Get the number of instruments (controllers) present in the RAW file
# and set the selected instrument to the MS instrument, first instance
# of it
print('The RAW file has data from {} instruments'.format(rawFile.InstrumentCount))

rawFile.SelectInstrument(Device.MS, 1)

# Get the first and last scan from the RAW file
firstScanNumber = rawFile.RunHeaderEx.FirstSpectrum
lastScanNumber = rawFile.RunHeaderEx.LastSpectrum

# Get the start and end time from the RAW file
startTime = rawFile.RunHeaderEx.StartTime
endTime = rawFile.RunHeaderEx.EndTime

# Print some OS and other information
print('System Information:')
print('   OS Version: {}'.format(Environment.OSVersion))
print('   64 bit OS: {}'.format(Environment.Is64BitOperatingSystem))
print('   Computer: {}'.format(Environment.MachineName))
print('   # Cores: {}'.format(Environment.ProcessorCount))
print('   Date: {}'.format(DateTime.Now))
print()

# Get some information from the header portions of the RAW file and
# display that information.  The information is general information
# pertaining to the RAW file.
print('General File Information:')
print('   RAW file: {}'.format(rawFile.FileName))
print('   RAW file version: {}'.format(rawFile.FileHeader.Revision))
print('   Creation date: {}'.format(rawFile.FileHeader.CreationDate))
print('   Operator: {}'.format(rawFile.FileHeader.WhoCreatedId))
print('   Number of instruments: {}'.format(rawFile.InstrumentCount))
print('   Description: {}'.format(rawFile.FileHeader.FileDescription))
print('   Instrument model: {}'.format(rawFile.GetInstrumentData().Model))
print('   Instrument name: {}'.format(rawFile.GetInstrumentData().Name))
print('   Serial number: {}'.format(rawFile.GetInstrumentData().SerialNumber))
print('   Software version: {}'.format(rawFile.GetInstrumentData().SoftwareVersion))
print('   Firmware version: {}'.format(rawFile.GetInstrumentData().HardwareVersion))
print('   Units: {}'.format(Enum.GetName(DataUnits, rawFile.GetInstrumentData().Units)))
print('   Mass resolution: {:.3f}'.format(rawFile.RunHeaderEx.MassResolution))
print('   Number of scans: {}'.format(rawFile.RunHeaderEx.SpectraCount))
print('   Scan range: {} - {}'.format(firstScanNumber, lastScanNumber))
print('   Time range: {:.2f} - {:.2f}'.format(startTime, endTime))
print('   Mass range: {:.4f} - {:.4f}'.format(rawFile.RunHeaderEx.LowMass, rawFile.RunHeaderEx.HighMass))
print()

# Get information related to the sample that was processed
print('Sample Information:')
print('   Sample name: {}'.format(rawFile.SampleInformation.SampleName))
print('   Sample id: {}'.format(rawFile.SampleInformation.SampleId))
print('   Sample type: {}'.format(Enum.GetName(SampleType, rawFile.SampleInformation.SampleType)))
print('   Sample comment: {}'.format(rawFile.SampleInformation.Comment))
print('   Sample vial: {}'.format(rawFile.SampleInformation.Vial))
print('   Sample volume: {}'.format(rawFile.SampleInformation.SampleVolume))
print('   Sample injection volume: {}'.format(rawFile.SampleInformation.InjectionVolume))
print('   Sample row number: {}'.format(rawFile.SampleInformation.RowNumber))
print('   Sample dilution factor: {}'.format(rawFile.SampleInformation.DilutionFactor))
print()

# Read the first instrument method (most likely for the MS portion of
# the instrument).  NOTE: This method reads the instrument methods
# from the RAW file but the underlying code uses some Microsoft code
# that hasn't been ported to Linux or MacOS.  Therefore this method
# won't work on those platforms therefore the check for Windows.
if 'Windows' in str(Environment.OSVersion):
    deviceNames = rawFile.GetAllInstrumentNamesFromInstrumentMethod()
    for device in deviceNames:
        print('Instrument method: {}'.format(device))
    print()


def ListTrailerExtraFields(rawFile):
    '''Reads and reports the trailer extra data fields present in the RAW
    file.

    Args:
        rawFile (IRawDataPlus): the RAW file.
    '''

    # Get the Trailer Extra data fields present in the RAW file
    trailerFields = rawFile.GetTrailerExtraHeaderInformation()

    # Display each value
    i = 0
    print('Trailer Extra Data Information:')

    for field in trailerFields:
        print('   Field {} = {} storing data of type {}'.format(
            i, field.Label, Enum.GetName(GenericDataTypes, field.DataType)))
        i += 1

    print()


# Display all of the trailer extra data fields present in the RAW file
ListTrailerExtraFields(rawFile)

# Get the number of filters present in the RAW file
numberFilters = rawFile.GetFilters().Count

# Get the scan filter for the first and last spectrum in the RAW file
firstFilter = IScanFilter(rawFile.GetFilterForScanNumber(firstScanNumber))
lastFilter = IScanFilter(rawFile.GetFilterForScanNumber(lastScanNumber))

print('Filter Information:')
print('   Scan filter (first scan): {}'.format(firstFilter.ToString()))
print('   Scan filter (last scan): {}'.format(lastFilter.ToString()))
print('   Total number of filters: {}'.format(numberFilters))
print()


def GetChromatogram(rawFile, startScan, endScan, outputData):
    '''Reads the base peak chromatogram for the RAW file.

    Args:
        rawFile (IRawDataPlus): the RAW file being read.
        startScan (int): start scan for the chromatogram.
        endScan (int): end scan for the chromatogram.
        outputData (bool): the output data flag.
    '''

    # Define the settings for getting the Base Peak chromatogram
    settings = ChromatogramTraceSettings(TraceType.BasePeak)

    # Get the chromatogram from the RAW file.
    data = rawFile.GetChromatogramData([settings], startScan, endScan)

    # Split the data into the chromatograms
    trace = ChromatogramSignal.FromChromatogramData(data)

    if trace[0].Length:
        # Print the chromatogram data (time, intensity values)
        print('Base Peak chromatogram ({} points)'.format(trace[0].Length))

    if outputData:
        for i in range(trace[0].Length):
            print('  {} - {:.3f}, {:.0f}'.format(i, trace[0].Times[i], trace[0].Intensities[i]))

    print()


# Get the BasePeak chromatogram for the MS data
GetChromatogram(rawFile, firstScanNumber, lastScanNumber, True)


def ReadScanInformation(rawFile, firstScanNumber, lastScanNumber, outputData):
    '''Reads the general scan information for each scan in the RAW file
    using the scan filter object and also the trailer extra data
    section for that same scan.

    Args:
        rawFile (IRawDataPlus): the RAW file being read.
        firstScanNumber (int): the first scan in the RAW file.
        lastScanNumber (int): the last scan in the RAW file.
        outputData (bool): the output data flag.
    '''

    # Read each scan in the RAW File
    for scan in range(firstScanNumber, lastScanNumber):
        # Get the retention time for this scan number.  This is one of
        # two comparable functions that will convert between retention
        # time and scan number.
        time = rawFile.RetentionTimeFromScanNumber(scan)

        # Get the scan filter for this scan number
        scanFilter = IScanFilter(rawFile.GetFilterForScanNumber(scan))

        # Get the scan event for this scan number
        scanEvent = IScanEventBase(rawFile.GetScanEventForScanNumber(scan))

        # Get the ionizationMode, MS2 precursor mass, collision
        # energy, and isolation width for each scan
        if scanFilter.MSOrder == MSOrderType.Ms2:
            # Get the reaction information for the first precursor
            reaction = scanEvent.GetReaction(0)

            precursorMass = reaction.PrecursorMass
            collisionEnergy = reaction.CollisionEnergy
            isolationWidth = reaction.IsolationWidth
            monoisotopicMass = 0.0
            masterScan = 0
            ionizationMode = scanFilter.IonizationMode
            order = scanFilter.MSOrder

            # Get the trailer extra data for this scan and then look
            # for the monoisotopic m/z value in the trailer extra data
            # list
            trailerData = rawFile.GetTrailerExtraInformation(scan)

            for i in range(trailerData.Length):
                if trailerData.Labels[i] == 'Monoisotopic M/Z:':
                    monoisotopicMass = float(trailerData.Values[i])
                elif trailerData.Labels[i] in ('Master Scan Number:', 'Master Scan Number', 'Master Index:'):
                    masterScan = int(trailerData.Values[i])

            if outputData:
                print(
                    '''Scan number {} @ time {:.2f} - Master scan = {}, Ionization mode={},\
                    MS Order={}, Precursor mass={:.4f}, Monoisotopic Mass = {:.4f},\
                    Collision energy={:.2f}, Isolation width={:.2f}'''.format(
                        scan, time, masterScan, Enum.GetName(IonizationModeType, ionizationMode),
                        Enum.GetName(MSOrderType, order), precursorMass, monoisotopicMass,
                        collisionEnergy, isolationWidth))

        elif scanFilter.MSOrder == MSOrderType.Ms:
            scanDependents = rawFile.GetScanDependents(scan, 5)

            print(
                'Scan number {} @ time {:.2f} - Instrument type={}, Number dependent scans={}'.format(
                    scan, time, Enum.GetName(RawFileClassification, scanDependents.RawFileInstrumentType),
                    scanDependents.ScanDependentDetailArray.Length))


# Read the scan information for each scan in the RAW file
ReadScanInformation(rawFile, firstScanNumber, lastScanNumber, True)


def GetSpectrum(rawFile, scanNumber, scanFilter, outputData):
    '''Gets the spectrum from the RAW file.

    Args:
        rawFile (IRawDataPlus): the RAW file being read.
        scanNumber (int): the scan number being read.
        scanFilter (str): the scan filter for that scan.
        outputData (bool): the output data flag.
    '''

    # Check for a valid scan filter
    if not scanFilter:
        return

    # Get the scan statistics from the RAW file for this scan number
    scanStatistics = rawFile.GetScanStatsForScanNumber(scanNumber)

    # Check to see if the scan has centroid data or profile data.  Depending upon the
    # type of data, different methods will be used to read the data.  While the ReadAllSpectra
    # method demonstrates reading the data using the Scan.FromFile method, generating the
    # Scan object takes more time and memory to do, so that method isn't optimum.
    if scanStatistics.IsCentroidScan:
        # Get the centroid (label) data from the RAW file for this
        # scan
        centroidStream = rawFile.GetCentroidStream(scanNumber, False)

        print('Spectrum (centroid/label) {} - {} points'.format(scanNumber, centroidStream.Length))

        # Print the spectral data (mass, intensity, charge values).
        # Not all of the information in the high resolution centroid
        # (label data) object is reported in this example.  Please
        # check the documentation for more information about what is
        # available in high resolution centroid (label) data.
        if outputData:
            for i in range(centroidStream.Length):
                print('  {} - {:.4f}, {:.0f}, {:.0f}'.format(
                    i, centroidStream.Masses[i], centroidStream.Intensities[i], centroidStream.Charges[i]))
            print()

    else:
        # Get the segmented (low res and profile) scan data
        segmentedScan = rawFile.GetSegmentedScanFromScanNumber(scanNumber, scanStatistics)

        print('Spectrum (normal data) {} - {} points'.format(scanNumber, segmentedScan.Positions.Length))

        # Print the spectral data (mass, intensity values)
        if outputData:
            for i in range(segmentedScan.Positions.Length):
                print('  {} - {:.4f}, {:.0f}'.format(
                    i, segmentedScan.Positions[i], segmentedScan.Intensities[i]))
            print()


# Get a spectrum from the RAW file.
GetSpectrum(rawFile, firstScanNumber, firstFilter.ToString(), False)


def GetAverageSpectrum(rawFile, firstScanNumber, lastScanNumber, outputData):
    '''Gets the average spectrum from the RAW file.

    Args:
        rawFile (IRawDataPlus): the RAW file being read.
        firstScanNumber (int): the first scan to consider for the averaged spectrum.
        lastScanNumber (int): the last scan to consider for the averaged spectrum.
        outputData (bool): the output data flag.
    '''

    # Create the mass options object that will be used when averaging
    # the scans
    options = Extensions.DefaultMassOptions(rawFile)

    options.ToleranceUnits = ToleranceUnits.ppm
    options.Tolerance = 5.0

    # Get the scan filter for the first scan.  This scan filter will be used to located
    # scans within the given scan range of the same type
    scanFilter = IScanFilter(rawFile.GetFilterForScanNumber(firstScanNumber))

    # Get the average mass spectrum for the provided scan range. In addition to getting the
    # average scan using a scan range, the library also provides a similar method that takes
    # a time range.
    averageScan = Extensions.AverageScansInScanRange(
        rawFile, firstScanNumber, lastScanNumber, scanFilter, options)

    if averageScan.HasCentroidStream:
        print('Average spectrum ({} points)'.format(averageScan.CentroidScan.Length))

        # Print the spectral data (mass, intensity values)
        if outputData:
            for i in range(averageScan.CentroidScan.Length):
                print('  {:.4f} {:.0f}'.format(
                    averageScan.CentroidScan.Masses[i], averageScan.CentroidScan.Intensities[i]))

    # This example uses a different method to get the same average spectrum that was calculated in the
    # previous portion of this method.  Instead of passing the start and end scan, a list of scans will
    # be passed to the GetAveragedMassSpectrum function.
    scans = List[int]([1, 6, 7, 9, 11, 12, 14])

    averageScan = Extensions.AverageScans(rawFile, scans, options)

    if averageScan.HasCentroidStream:
        print('Average spectrum ({} points)', averageScan.CentroidScan.Length)

        # Print the spectral data (mass, intensity values)
        if outputData:
            for i in range(averageScan.CentroidScan.Length):
                print('  {:.4f} {:.0f}'.format(
                    averageScan.CentroidScan.Masses[i], averageScan.CentroidScan.Intensities[i]))

    print()


# Get a average spectrum from the RAW file for the first 15 scans in the file.
GetAverageSpectrum(rawFile, 1, 15, False)


def ReadAllSpectra(rawFile, firstScanNumber, lastScanNumber, outputData):
    '''Read all spectra in the RAW file.

    Args:
        rawFile (IRawDataPlus): the raw file.
        firstScanNumber (int): the first scan number.
        lastScanNumber (int): the last scan number.
        outputData (bool): the output data flag.
    '''

    for scanNumber in range(firstScanNumber, lastScanNumber):
        try:
            # Get the scan filter for the spectrum
            scanFilter = IScanFilter(rawFile.GetFilterForScanNumber(firstScanNumber))

            if not scanFilter.ToString():
                continue

            # Get the scan from the RAW file.  This method uses the Scan.FromFile method which returns a
            # Scan object that contains both the segmented and centroid (label) data from an FTMS scan
            # or just the segmented data in non-FTMS scans.  The GetSpectrum method demonstrates an
            # alternative method for reading scans.
            scan = Scan.FromFile(rawFile, scanNumber)

            # If that scan contains FTMS data then Centroid stream
            # will be populated so check to see if it is present.
            if scan.HasCentroidStream:
                labelSize = scan.CentroidScan.Length
            else:
                labelSize = 0

            # For non-FTMS data, the preferred data will be populated
            if scan.PreferredMasses is not None:
                dataSize = scan.PreferredMasses.Length
            else:
                dataSize = 0

            if outputData:
                print('Spectrum {} - {}: normal {}, label {} points'.format(
                    scanNumber, scanFilter.ToString(), dataSize, labelSize))

        except Exception as ex:
            print('Error reading spectrum {} - {}'.format(scanNumber, str(ex)))


# Read each spectrum
ReadAllSpectra(rawFile, firstScanNumber, lastScanNumber, True)


def CalculateMassPrecision(rawFile, scanNumber):
    '''Calculates the mass precision for a spectrum.

    Args:
        rawFile (IRawDataPlus): the RAW file being read.
        scanNumber (int): the scan to process.
    '''

    # Get the scan from the RAW file
    scan = Scan.FromFile(rawFile, scanNumber)

    # Get the scan event and from the scan event get the analyzer type for this scan
    scanEvent = IScanEventBase(rawFile.GetScanEventForScanNumber(scanNumber))

    # Get the trailer extra data to get the ion time for this file
    logEntry = rawFile.GetTrailerExtraInformation(scanNumber)

    trailerHeadings = List[str]()
    trailerValues = List[str]()
    for i in range(logEntry.Length):
        trailerHeadings.Add(logEntry.Labels[i])
        trailerValues.Add(logEntry.Values[i])

    # Create the mass precision estimate object
    precisionEstimate = PrecisionEstimate()

    # Get the ion time from the trailer extra data values
    ionTime = precisionEstimate.GetIonTime(scanEvent.MassAnalyzer, scan, trailerHeadings, trailerValues)

    # Calculate the mass precision for the scan
    listResults = precisionEstimate.GetMassPrecisionEstimate(
        scan, scanEvent.MassAnalyzer, ionTime, rawFile.RunHeader.MassResolution)

    # Output the mass precision results
    if listResults.Count:
        print('Mass Precision Results:')

        for result in listResults:
            print('Mass {:.5f}, mmu = {:.3f}, ppm = {:.2f}'.format(
                result.Mass, result.MassAccuracyInMmu, result.MassAccuracyInPpm))


# Calculate the mass precision for a spectrum
CalculateMassPrecision(rawFile, 1)

# Close (dispose) the RAW file
print()
print('Closing {}'.format(filename))

rawFile.Dispose()
