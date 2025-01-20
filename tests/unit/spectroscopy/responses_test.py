# -*- coding: utf-8 -*-
import unittest

import tempfile
import os

import numpy

from nose.tools import assert_raises

import quantarhei as qr

from quantarhei.utils.vectors import X 

import matplotlib.pyplot as plt


"""
*******************************************************************************


    Tests of the quantarhei.spectroscopy.responses package


*******************************************************************************
"""


class TestResponses(unittest.TestCase):
    """Tests for the response package
    
    
    """
    
    def setUp(self,verbose=False):
        pass

    """
    def test_MockTwoD(self):

        E0 = 12000.0
        E1 = E0 + 100.0
        #
        # Create a system and calculate 2D spectrum
        #
        with qr.energy_units("1/cm"):
            m1 = qr.Molecule([0.0, E0])
            m2 = qr.Molecule([0.0, E1])
            
            m1.set_dipole(0,1,[2.0, 0.0, 0.0])
            m2.set_dipole(0,1,[0.5, 0.5, 0.0])
            
            m1.position = [0.0, 0.0, 0.0]
            m2.position = [0.0, 0.0, 2.0]
            
            
        agg = qr.Aggregate(molecules=[m1, m2])
        agg.set_coupling_by_dipole_dipole()
        
        with qr.energy_units("1/cm"):
            m1.set_transition_width((0,1), 200)
            m2.set_transition_width((0,1), 200)
            
        
        agg2 = agg.deepcopy()
        agg.build(mult=1)
        
        H = agg.get_Hamiltonian()
        
        gammas = [1.0/200.0]
        with qr.eigenbasis_of(H):
            
            K = qr.qm.ProjectionOperator(1, 2, dim=H.dim)
            
        sbi = qr.qm.SystemBathInteraction(sys_operators=[K], rates=gammas,
                                          system=agg)
        
        LL = qr.qm.LindbladForm(H, sbi)
        
        t2axis = qr.TimeAxis(0.0, 200, 5.0)
        t1axis = qr.TimeAxis(0.0, 100, 10.0)
        t3axis = qr.TimeAxis(0.0, 100, 10.0)
        
        eUt = qr.EvolutionSuperOperator(time=t2axis, ham=H, relt=LL)
        eUt.set_dense_dt(10)
        eUt.calculate()

        msc = qr.MockTwoDSpectrumCalculator(t1axis, t2axis, t3axis)
        msc.bootstrap(rwa=qr.convert((E0+E1)/2.0,"1/cm","int"), 
                      all_positive=False, shape="Gaussian")
    
        agg2.build(mult=2)
        agg2.diagonalize()

        lab = qr.LabSetup()
        lab.set_polarizations(pulse_polarizations=[X,X,X],
                              detection_polarization=X)   
        
        cont = qr.TwoDSpectrumContainer(t2axis=t2axis)
        
        for t2 in t2axis.data:
            twod = msc.calculate_one_system(t2, agg2, H, eUt, lab) 
        
            cont.set_spectrum(twod)        
        
        cont.set_data_flag(qr.signal_TOTL)
        sp = cont.get_spectrum(800.0)
        
        with qr.energy_units("1/cm"):
            sp.plot(show=True)
            sp.savefig("twod.png")
#            
#            cont.make_movie("movie.mp4")
            
        pevol = cont.get_point_evolution(qr.convert(12600,"1/cm","int"),
                                         qr.convert(11500,"1/cm","int"),
                                         t2axis)
        
        plt.plot(t2axis.data, pevol)
        plt.show()
        
        
        pars = cont.global_fit_exponential()
        
        print("Optimum: ", pars.x)
        
        raise Exception()
        
    """
        
    def test_Responses(self):
        """Testing basic functions of the TwoDSpectrumBase class
        
        """
        t1 = qr.TimeAxis(0.0, 1000, 1.0)
        t3 = qr.TimeAxis(0.0, 1000, 1.0)
        
        resp = qr.ResponseFunction("R2g")
        
        def efce(t2,t1s,t3s,et13):
            
            N1 = t1s.shape[0]
            N3 = t3s.shape[0]
            ret = numpy.zeros((N1,N3), dtype=complex)
            ret[:,:] = 1.0*et13[:,:]
            
            return ret
        
        
        resp.set_evaluation_function(efce)
        
        d1 = [1.0, 0.0, 0.0]
        
        resp.set_dipoles(d1)
        
        rwa = 1.0
        mx = resp.calculate_matrix(None, None, 10, t1.data, t3.data, rwa)

        
        #with assert_raises(Exception):      
        #    twodB.data = data2


        #numpy.testing.assert_equal(twodB.data, data2)
        
        #with tempfile.TemporaryDirectory() as tdir:
        
        #    numpy.testing.assert_equal(twodB.data, twodB2.data)


    
    def test_LouvillePathway(self):
        """Testing basic functions of the TwoDSpectrumCalculator class
        
        """
        t1 = qr.TimeAxis(0.0, 1000, 1.0)
        t3 = qr.TimeAxis(0.0, 1000, 1.0)
        
        t2 = qr.TimeAxis(30, 10, 10.0)
        
        #twod_calc = qr.TwoDResponseCalculator(t1, t2, t3)
        
        
            
    
