from .tables import *

class UserSpectrum(ResonanceData):
    __doc__="""Class to handle an arbitrary two-column CSV format neutron 
    flux spectrum."""

    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)

    def import_spectrum(self,csv_filename):
        """Import arbitrary neutron flux spectrum as a function of energy in 
        CSV format."""
        self.csv_filename = csv_filename
        with open(self.csv_filename, mode='r', encoding='utf-8-sig') as file:
            csvFile = csv.reader(file)
            next(csvFile, None)  # Skip header line
            En = []
            dndE = []
            for lines in csvFile:
                En.append(float(lines[0]))
                dndE.append(float(lines[1])) 
            return(np.array(En), np.array(dndE))

