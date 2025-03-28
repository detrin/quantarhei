    # -*- coding: utf-8 -*-
"""
    Quantarhei User Level Classes and Objects
    =========================================

    In Quantarhei, classes are losely grouped into three categories. First,
    there is agroup of classes, which represent basic concepts of quantum
    mechanics, provide access to important implementations of spectroscopic
    simulations and dynamics of open quantum systems, and classes which allow
    basic management of the simulation environment and numerical results.
    These classed are called **user level classes**, and they are all
    accessible in highest namespace level of the Quantarhei package.
    If you import Quantarhei like this:

    >>> import quantarhei as qr

    you can access user level classes through the qr. prefix, e.g.


    >>> manager = qr.Manager()
    >>> print(manager.version)
    0.0.63

    The list of user level classes is provided below. Tue latest and most
    uptodate information can be obtained by viewing the source code of the 
    root `__init__.py` file of the packages. All classes imported there are
    considered user level classes.
    
    
    Other Class Levels
    ------------------
    
    In this documentation we recognize two more groups (or levels) of classes.
    More specialized classes, which normal user does not need as often as the
    user level classes are called **advanced level classes**. These use the
    second level name space. For instance the class `SystemBathInteraction`
    is relatively rarely used directly. It is therefore *hidden* in the name
    space `qm` (as quantum mechanics) of the package. This class can be 
    instantiated e.g. like this
    
    >>> import quantarhei as qr
    >>> sbi = qr.qm.SystemBathInteraction()
    
    Advanced level classes are still intendend for relatively frequent use
    by the user. However, in order to reduced the *apparent* complexity of
    basic usage of Quantarhei, advanced level classes are documented in their
    respective sub-packages, one level deeper than user level classes. Complete
    documentation of advanced level classes is available in the Advanced Level
    Classes section of this documentation.
    
    Everything else in Quantarhei package goes under the banner of
    **expert level classes**. This includes all classes and objects used 
    internally in Quantarhei. We make every effort to document also this part
    of the package as completely as possible, but it is the last item on the
    list, so to say. The user is welcome to learn and use the expert level
    classes, but our aim is to structure Quantarhei in such a way, that this
    is not necessary. More on expert level classes in the section in 
    Quantarhei internals.

    User Level Objects and Convenience Functions
    ============================================
    
    Besides classes, Quantarhei also defines some user level objects and
    convenience functions. They are listed here under several categories
  
    Numeric types
    -------------

    .. toctree::
        :maxdepth: 2
        
        functions/numtypes        
    
    Convenience Functions
    ---------------------
    
    .. toctree::
        :maxdepth: 2
        
        functions/convenience
        
    
    Logging Functions and Loglevels
    -------------------------------

    .. toctree::
        :maxdepth: 2
        
        functions/logging

   
    .. 
        Builders
        --------
        
        Mode .......... represents a harmonic vibrational mode of a molecule
        Molecule ...... represents a molecule
        Aggregate ..... represents an aggregate of molecules
        PDBFile ....... reader and writter of structures from PDB format
        Disorder ...... class managing static disorder of molecular transition
                        energies
         
        Core classes
        ------------
        
        TimeAxis ......... linear axis of real values representing discrete time
        FrequencyAxis .... linear axis of real values representing discrete
                           frequency axis
        DFunction ........ discrete function
        
        
        Various managers
        ----------------
        
        Manager ............ the main behind-the-scenes manager of the package
        energy_units ....... energy units manager for use with the "with" construct
        frequency_units .... frequency units manager for use with 
                             the "with" construct
        eigenbasis_of ...... manager of the basis transformations to be used with 
                             the "with" construct
        set_current_units .. function to set current units globally
        
        ... to be continued


"""


###############################################################################
#
#
#            Imports of high level classes and functions 
#
#
###############################################################################

#
# Fix used numerical types
#
#import numpy
from quantarhei.core.managers import Manager
m = Manager()

REAL = m.get_real_type() #numpy.float64
COMPLEX = m.get_complex_type() #numpy.complex128

LOG_URGENT = 1
LOG_REPORT = 3
LOG_INFO = 5
LOG_DETAIL = 7
LOG_QUICK = 9


#
# Non-linear response signals
#
signal_REPH = "rephasing_2D_signal"
signal_NONR = "nonrephasing_2D_signal"
signal_TOTL = "total_2D_signal"
signal_DC = "double_coherence_signal"

TWOD_SIGNALS = dict(signal_REPH = "rephasing_2D_signal",
                    signal_NONR = "nonrephasing_2D_signal",
                    signal_TOTL = "total_2D_signal",
                    signal_DC = "double_coherence_signal")

#
# Parts of the complex data/signal
#
part_REAL = "real_part"
part_IMAGINARY = "imaginary_part"
part_COMPLEX = "complex"
part_ABS = "absolute_value"
part_PHASE = "phase"

SIGNAL_PARTS = dict(part_REAL = "real_part",
                    part_IMAGINARY= "imaginary_part",
                    part_COMPLEX = "complex",
                    part_ABS = "absolute_value",
                    part_PHASE = "phase")

DATA_PARTS = SIGNAL_PARTS

#
# Liouville pathway types
#
ptype_R1g = "pathway_type_R1g"
ptype_R2g = "pathway_type_R2g"
ptype_R3g = "pathway_type_R3g"
ptype_R4g = "pathway_type_R4g"
ptype_R1f = "pathway_type_R1f*"
ptype_R2f = "pathway_type_R2f*"
ptype_R3f = "pathway_type_R3f"
ptype_R4f = "pathway_type_R4f"

PATHWAY_TYPES = dict(ptype_R1g = "pathway_type_R1g",
                     ptype_R2g = "pathway_type_R2g",
                     ptype_R3g = "pathway_type_R3g",
                     ptype_R4g = "pathway_type_R4g",
                     ptype_R1f = "pathway_type_R1f*",
                     ptype_R2f = "pathway_type_R2f*",
                     ptype_R3f = "pathway_type_R3f",
                     ptype_R4f = "pathway_type_R4f")

LIOUVILLE_PATHWAY_TYPES = PATHWAY_TYPES

#
# Builders
#
from quantarhei.builders.modes import Mode
from quantarhei.builders.molecules import Molecule
from quantarhei.builders.molecule_test import TestMolecule
from quantarhei.builders.aggregates import Aggregate
from quantarhei.builders.aggregate_test import TestAggregate
from quantarhei.builders.pdb import PDBFile
from quantarhei.builders.disorder import Disorder

#
# Core classes
#
from quantarhei.core.time import TimeAxis
from quantarhei.core.frequency import FrequencyAxis
from quantarhei.core.valueaxis import ValueAxis
from quantarhei.core.dfunction import DFunction
#from quantarhei.core.saveable import Saveable

from quantarhei.core.saveable import Saveable
from quantarhei.core.parcel import Parcel

#
# Various managers
#
from quantarhei.core.managers import energy_units
from quantarhei.core.managers import frequency_units
from quantarhei.core.managers import length_units
from quantarhei.core.managers import eigenbasis_of
from quantarhei.core.managers import set_current_units

#
# Parallelization
#
from quantarhei.core.parallel import distributed_configuration
from quantarhei.core.parallel import start_parallel_region
from quantarhei.core.parallel import close_parallel_region
from quantarhei.core.parallel import parallel_function
from quantarhei.core.parallel import block_distributed_range
from quantarhei.core.parallel import block_distributed_list
from quantarhei.core.parallel import block_distributed_array
from quantarhei.core.parallel import collect_block_distributed_data
from quantarhei.core.parallel import asynchronous_range

###############################################################################
#                            SPECTROSCOPY
###############################################################################

#
# Linear absorption 
#
from quantarhei.spectroscopy.abs2 import AbsSpectrum
from quantarhei.spectroscopy.abscontainer import AbsSpectrumContainer
from quantarhei.spectroscopy.abscalculator import AbsSpectrumCalculator
from quantarhei.spectroscopy.mockabscalculator import MockAbsSpectrumCalculator
#
# Fluorescence
#
from quantarhei.spectroscopy.fluorescence import FluorSpectrum
from quantarhei.spectroscopy.fluorescence import FluorSpectrumContainer
from quantarhei.spectroscopy.fluorescence import FluorSpectrumCalculator
#
# Linear dichroism
#
from quantarhei.spectroscopy.linear_dichroism import LinDichSpectrum
from quantarhei.spectroscopy.linear_dichroism import LinDichSpectrumContainer
from quantarhei.spectroscopy.linear_dichroism import LinDichSpectrumCalculator
#
# Circular dichroism
#
from quantarhei.spectroscopy.circular_dichroism import CircDichSpectrum
from quantarhei.spectroscopy.circular_dichroism import CircDichSpectrumContainer
from quantarhei.spectroscopy.circular_dichroism import CircDichSpectrumCalculator
#
# Fourier transform Two-Dimensional Spectra
#
from quantarhei.spectroscopy.twod2 import TwoDResponse 
from quantarhei.spectroscopy.twodcontainer import TwoDResponseContainer, TwoDSpectrumContainer
from quantarhei.spectroscopy.twod import TwoDSpectrum
from quantarhei.spectroscopy.twodcalculator import TwoDResponseCalculator
from quantarhei.spectroscopy.mocktwodcalculator import MockTwoDResponseCalculator
from quantarhei.spectroscopy.responses import ResponseFunction
from quantarhei.spectroscopy.responses import LiouvillePathway

#
# Pump-probe spectrum
#
from quantarhei.spectroscopy.pumpprobe import PumpProbeSpectrum
from quantarhei.spectroscopy.pumpprobe import PumpProbeSpectrumContainer
from quantarhei.spectroscopy.pumpprobe import PumpProbeSpectrumCalculator
from quantarhei.spectroscopy.pumpprobe import MockPumpProbeSpectrumCalculator

from quantarhei.spectroscopy.pathwayanalyzer import LiouvillePathwayAnalyzer

from quantarhei.spectroscopy.labsetup import LabSetup
from quantarhei.spectroscopy.labsetup import LabField

from quantarhei.spectroscopy.dsfeynman import DSFeynmanDiagram
from quantarhei.spectroscopy.dsfeynman import R1g_Diagram
from quantarhei.spectroscopy.dsfeynman import R2g_Diagram
from quantarhei.spectroscopy.dsfeynman import R3g_Diagram
from quantarhei.spectroscopy.dsfeynman import R4g_Diagram
from quantarhei.spectroscopy.dsfeynman import R1f_Diagram
from quantarhei.spectroscopy.dsfeynman import R2f_Diagram


###############################################################################
#                           QUANTUM MECHANICS
###############################################################################


#
# State vectors
#
from quantarhei.qm import StateVector
from quantarhei.qm import OQSStateVector
#
# Operators
#
from quantarhei.qm import DensityMatrix
from quantarhei.qm import ReducedDensityMatrix
from quantarhei.qm import BasisReferenceOperator
from quantarhei.qm import Hamiltonian
from quantarhei.qm import Liouvillian
from quantarhei.qm import TransitionDipoleMoment
from quantarhei.qm import UnityOperator

#
# Propagators
#
from quantarhei.qm.propagators.poppropagator import PopulationPropagator
from quantarhei.qm.propagators.svpropagator import StateVectorPropagator
from quantarhei.qm import OQSStateVectorPropagator
from quantarhei.qm import ReducedDensityMatrixPropagator

#
# Evolutions (time-dependent operators)
#
from quantarhei.qm.propagators.statevectorevolution import StateVectorEvolution
from quantarhei.qm import OQSStateVectorEvolution
from quantarhei.qm import DensityMatrixEvolution
from quantarhei.qm import ReducedDensityMatrixEvolution

#
# Evolution operators
#
from quantarhei.qm.liouvillespace.evolutionsuperoperator import EvolutionSuperOperator



#
# System-bath interaction
#
from quantarhei.qm.corfunctions import CorrelationFunction
from quantarhei.qm.corfunctions import LineshapeFunction
from quantarhei.qm.corfunctions import SpectralDensity
from quantarhei.qm.corfunctions.correlationfunctions import oscillator_scalled_CorrelationFunction

#
# LINESHAPE FUNCTIONS
#
from quantarhei.qm.corfunctions import FunctionStorage
from quantarhei.qm.corfunctions import FastFunctionStorage


from quantarhei.qm.liouvillespace.heom import KTHierarchy
from quantarhei.qm.liouvillespace.heom import KTHierarchyPropagator

from quantarhei.symbolic.cumulant import evaluate_cumulant

###############################################################################
# Convenience functions
###############################################################################
#from quantarhei.core.saveable import load
#from quantarhei.core.saveable import read_info

from quantarhei.core.parcel import save_parcel
from quantarhei.core.parcel import load_parcel
from quantarhei.core.parcel import check_parcel

from quantarhei.core.units import convert
from quantarhei.core.units import in_current_units

from quantarhei.utils.vectors import normalize2
from quantarhei.utils.vectors import norm 

from quantarhei.utils.logging import printlog
from quantarhei.utils.logging import loglevels2bool
from quantarhei.utils.logging import log_urgent
from quantarhei.utils.logging import log_report
from quantarhei.utils.logging import log_info
from quantarhei.utils.logging import log_detail
from quantarhei.utils.logging import log_quick
from quantarhei.utils.logging import log_to_file
from quantarhei.utils.logging import init_logging

from quantarhei.utils.logging import tprint

from quantarhei.utils.timing import timeit
from quantarhei.utils.timing import untimeit
from quantarhei.utils.timing import finished_in
from quantarhei.utils.timing import done_in

from quantarhei.utils.paver import execute_paver

from quantarhei.wizard.input.input import Input


def exit(msg=None):
    """Exit to the level above the script with SystemExit exception
    
    """
    import sys
    if msg is not None:
        printlog("\n(SystemExit) Message: "+msg+"\n", loglevel=0)
    sys.exit()


def stop(msg=None):
    """Stop execution and leave to level above
    
    """

    exit("Execution stopped")
    
    
    
def show_plot(block=True):
    """Shows current plot
    
    This function is used to avoid explicit import of matplotlib
    """
    import matplotlib.pyplot as plt
    plt.show(block=block)



def savefig(fname):
    """Saves current plot to a file
    
    This function is used to avoid explicit import of matplotlib
    """
    import matplotlib.pyplot as plt
    plt.savefig(fname)    


def assert_version(check, vno):
    """Throws an exception if the condition is not satisfied
    
    """
    from packaging import version
    
    def ext():
        exit("Version requirement not satisfied.")

    if check == ">=":
        if not (version.parse(Manager().version) >= version.parse(vno)):
            ext()
            
    elif check == "==":
         if not (version.parse(Manager().version) == version.parse(vno)):
            ext()
            
    elif check == "<=":
        if not (version.parse(Manager().version) <= version.parse(vno)):
            ext()

    elif check == ">":
        if not (version.parse(Manager().version) == version.parse(vno)):
            ext()
            
    elif check == "<=":
        if not (version.parse(Manager().version) >= version.parse(vno)):
            ext()
            
    else:
        raise Exception("Unknown comparison operator `"+check+"`")
        
    