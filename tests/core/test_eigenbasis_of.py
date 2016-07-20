# -*- coding: utf-8 -*-

import unittest

"""
*******************************************************************************


    Tests of the quantarhei.eigenbasis_of class


*******************************************************************************
"""        
        
import numpy
import itertools

from quantarhei import eigenbasis_of
from .test_BasisManaged import BasisManagedObject
    
        

class TestEigenbasisOf(unittest.TestCase):
    
    def setUp(self):

        hval = numpy.array([[0.1, 1.0],
                            [1.0, 0.0]])
        dval = numpy.array([[0.1, 3.0],
                            [3.0, 0.0]])

        self.H = BasisManagedObject(hval,"H")
        self.B = BasisManagedObject(numpy.diag(numpy.ones(2)),"B")
        self.B.data[1,1] = 0.1
        self.P = BasisManagedObject(hval.copy(),"P")
        self.P.data[0,0] = -0.2
        self.D = BasisManagedObject(dval,"D")

        self.h_copy = self.H.data.copy()
        self.b_copy = self.B.data.copy()
        self.p_copy = self.P.data.copy() 
        self.d_copy = self.D.data.copy()
        
        
    def test_of_single_context_single_var(self):
        """Testing single operator in a single basis management context
        
        
        """
        x,Sh = numpy.linalg.eigh(self.H.data)
        Sh1 = numpy.linalg.inv(Sh)
        h_copy_tr = numpy.dot(Sh1,numpy.dot(self.h_copy,Sh))
        
        with eigenbasis_of(self.H):
            self.assertTrue(numpy.allclose(self.H.data,h_copy_tr))
        
        self.assertTrue(numpy.allclose(self.H.data,self.h_copy))
        
        
        
    def test_of_many_contexts_single_var(self):
        """Testing single operator in nested basis management contexts
        
        
        """
        x,Sh = numpy.linalg.eigh(self.H.data)
        Sh1 = numpy.linalg.inv(Sh)
        h_copy_tr1 = numpy.dot(Sh1,numpy.dot(self.h_copy,Sh))
        b_copy_tr1 = numpy.dot(Sh1,numpy.dot(self.b_copy,Sh))
        
        x,Sb = numpy.linalg.eigh(b_copy_tr1)
        Sb1 = numpy.linalg.inv(Sb)
        h_copy_tr2 = numpy.dot(Sb1,numpy.dot(h_copy_tr1,Sb))   
        
        with eigenbasis_of(self.H):
            # in context 1
            #print("In 1",self.H.manager.get_current_basis())
            self.assertTrue(numpy.allclose(self.H.data,h_copy_tr1))
            with eigenbasis_of(self.B):   
                # in context 2
                #print("In 2",self.B.manager.get_current_basis())
                self.assertTrue(numpy.allclose(self.H.data,h_copy_tr2))
                
            # in context 1
            #print("Out 2 in 1",self.H.manager.get_current_basis())
            self.assertTrue(numpy.allclose(self.H.data,h_copy_tr1))
        
        # outside contexts
        #print("Out 1",self.H.manager.get_current_basis())
        self.assertTrue(numpy.allclose(self.H.data,self.h_copy))


    def test_of_many_contexts_many_vars_1(self):
        """Testing several operators in nested basis management contexts
        
        
        """
        x,Sh = numpy.linalg.eigh(self.H.data)
        Sh1 = numpy.linalg.inv(Sh)
        h_copy_tr1 = numpy.dot(Sh1,numpy.dot(self.h_copy,Sh))
        b_copy_tr1 = numpy.dot(Sh1,numpy.dot(self.b_copy,Sh))
        
        x,Sb = numpy.linalg.eigh(b_copy_tr1)
        Sb1 = numpy.linalg.inv(Sb)
        h_copy_tr2 = numpy.dot(Sb1,numpy.dot(h_copy_tr1,Sb))   
        b_copy_tr2 = numpy.dot(Sb1,numpy.dot(b_copy_tr1,Sb))
        
        with eigenbasis_of(self.H):
            # in context 1
            self.assertTrue(numpy.allclose(self.H.data,h_copy_tr1))
            self.assertTrue(numpy.allclose(self.B.data,b_copy_tr1))
            with eigenbasis_of(self.B):   
                # in context 2
                #print("In 2",self.B.manager.get_current_basis())
                self.assertTrue(numpy.allclose(self.H.data,h_copy_tr2))
                self.assertTrue(numpy.allclose(self.B.data,b_copy_tr2))
            # in context 1
            #print("Out 2 in 1",self.H.manager.get_current_basis())
            self.assertTrue(numpy.allclose(self.H.data,h_copy_tr1))
            self.assertTrue(numpy.allclose(self.B.data,b_copy_tr1))
        
        # outside contexts
        #print("Out 1",self.H.manager.get_current_basis())
        self.assertTrue(numpy.allclose(self.H.data,self.h_copy))
        self.assertTrue(numpy.allclose(self.B.data,self.b_copy))

        
    def test_of_many_contexts_many_vars_2(self):
        """Testing deeply nested basis management contexts (takes a while)
        
        
        """
        x,Sh = numpy.linalg.eigh(self.H.data)
        Sh1 = numpy.linalg.inv(Sh)
        h_copy_tr1 = numpy.dot(Sh1,numpy.dot(self.h_copy,Sh))
        b_copy_tr1 = numpy.dot(Sh1,numpy.dot(self.b_copy,Sh))
        p_copy_tr1 = numpy.dot(Sh1,numpy.dot(self.p_copy,Sh))
        d_copy_tr1 = numpy.dot(Sh1,numpy.dot(self.d_copy,Sh))
       
        x,Sb = numpy.linalg.eigh(b_copy_tr1)
        Sb1 = numpy.linalg.inv(Sb)
        h_copy_tr2 = numpy.dot(Sb1,numpy.dot(h_copy_tr1,Sb))   
        b_copy_tr2 = numpy.dot(Sb1,numpy.dot(b_copy_tr1,Sb))
        p_copy_tr2 = numpy.dot(Sb1,numpy.dot(p_copy_tr1,Sb))
        d_copy_tr2 = numpy.dot(Sb1,numpy.dot(d_copy_tr1,Sb))
        
        x,Sp = numpy.linalg.eigh(p_copy_tr2)
        Sp1 = numpy.linalg.inv(Sp)
        h_copy_tr3 = numpy.dot(Sp1,numpy.dot(h_copy_tr2,Sp))   
        b_copy_tr3 = numpy.dot(Sp1,numpy.dot(b_copy_tr2,Sp))
        p_copy_tr3 = numpy.dot(Sp1,numpy.dot(p_copy_tr2,Sp))
        d_copy_tr3 = numpy.dot(Sp1,numpy.dot(d_copy_tr2,Sp))

        x,Sd = numpy.linalg.eigh(d_copy_tr3)
        Sd1 = numpy.linalg.inv(Sd)
        h_copy_tr4 = numpy.dot(Sd1,numpy.dot(h_copy_tr3,Sd))   
        b_copy_tr4 = numpy.dot(Sd1,numpy.dot(b_copy_tr3,Sd))
        p_copy_tr4 = numpy.dot(Sd1,numpy.dot(p_copy_tr3,Sd))
        d_copy_tr4 = numpy.dot(Sd1,numpy.dot(d_copy_tr3,Sd))



#        knock_bunch = [[[True,True,True,True],
#                        [True,True,True,True],
#                        [True,True,True,True],
#                        [True,True,True,True]],
#                       [[True,False,True,True],
#                        [True,True,False,True],
#                        [True,False,True,True],
#                        [True,True,False,True]]]

        count = 0
        for k in generate(only=12):
            knocks = [k[0:4],k[4:8],k[8:12],k[12:16]]
            count += 1
            #print(count)

            with eigenbasis_of(self.H):
                # in context 1
                if knocks[0][0]:
                    self.assertTrue(numpy.allclose(self.H.data,h_copy_tr1))
                if knocks[0][1]:
                    self.assertTrue(numpy.allclose(self.B.data,b_copy_tr1))
                if knocks[0][2]:
                    self.assertTrue(numpy.allclose(self.P.data,p_copy_tr1))
                if knocks[0][3]:
                    self.assertTrue(numpy.allclose(self.D.data,d_copy_tr1))
                
                with eigenbasis_of(self.B):   
                    # in context 2
                    #print("In 2",self.B.manager.get_current_basis())
                    if knocks[1][0]:
                        self.assertTrue(numpy.allclose(self.H.data,h_copy_tr2))
                    if knocks[1][1]:
                        self.assertTrue(numpy.allclose(self.B.data,b_copy_tr2))
                    if knocks[1][2]:
                        self.assertTrue(numpy.allclose(self.P.data,p_copy_tr2))
                    if knocks[1][3]:
                        self.assertTrue(numpy.allclose(self.D.data,d_copy_tr2))
               
                    with eigenbasis_of(self.P):
                        if knocks[2][0]:
                            self.assertTrue(numpy.allclose(self.H.data,
                                                           h_copy_tr3))
                        if knocks[2][1]:                        
                            self.assertTrue(numpy.allclose(self.B.data,
                                                           b_copy_tr3))
                        if knocks[2][2]:
                            self.assertTrue(numpy.allclose(self.P.data,
                                                           p_copy_tr3))
                        if knocks[2][3]:
                            self.assertTrue(numpy.allclose(self.D.data,
                                                           d_copy_tr3))
                    
                        with eigenbasis_of(self.D):
                            if knocks[3][0]:
                                self.assertTrue(numpy.allclose(self.H.data,
                                                               h_copy_tr4))
                            if knocks[3][1]:
                                self.assertTrue(numpy.allclose(self.B.data,
                                                               b_copy_tr4))
                            if knocks[3][2]:
                                self.assertTrue(numpy.allclose(self.P.data,
                                                               p_copy_tr4))
                            if knocks[3][3]:
                                self.assertTrue(numpy.allclose(self.D.data,
                                                               d_copy_tr4))
                    
                
                # in context 1
                #print("Out 2 in 1",self.H.manager.get_current_basis())
                self.assertTrue(numpy.allclose(self.H.data,h_copy_tr1))
                self.assertTrue(numpy.allclose(self.B.data,b_copy_tr1))
                self.assertTrue(numpy.allclose(self.P.data,p_copy_tr1))
                self.assertTrue(numpy.allclose(self.D.data,d_copy_tr1))
            
                with eigenbasis_of(self.P):
                    self.assertTrue(numpy.allclose(self.P.data,p_copy_tr3))
                    
                with eigenbasis_of(self.B):   
                    self.assertTrue(numpy.allclose(self.D.data,d_copy_tr2))
        
                with eigenbasis_of(self.D): 
                    #print(self.B.data)
                    #print(b_copy_tr4)
                    #self.assertTrue(numpy.allclose(self.B.data,b_copy_tr4))
                    with eigenbasis_of(self.H):
                        self.assertTrue(numpy.allclose(self.D.data,d_copy_tr1))
                
            # outside contexts
            #print("Out 1",self.H.manager.get_current_basis())
            self.assertTrue(numpy.allclose(self.H.data,self.h_copy))
            self.assertTrue(numpy.allclose(self.B.data,self.b_copy))
            self.assertTrue(numpy.allclose(self.P.data,self.p_copy))
            self.assertTrue(numpy.allclose(self.D.data,self.d_copy))
            
        


def generate(only=-1):
    N = 16
    if only == -1:
        for i in range(N+1):
            l1 = [0]*(N-i)
            l2 = [1]*i
            l = l1 + l2
            for j in perm_unique(l):
                yield j
    else:
        l1 = [0]*(N-only)
        l2 = [1]*only
        l = l1 + l2
        for j in perm_unique(l):
            yield j
            
            
class unique_element:
    def __init__(self,value,occurrences):
        self.value = value
        self.occurrences = occurrences

def perm_unique(elements):
    eset=set(elements)
    listunique = [unique_element(i,elements.count(i)) for i in eset]
    u=len(elements)
    return perm_unique_helper(listunique,[0]*u,u-1)

def perm_unique_helper(listunique,result_list,d):
    if d < 0:
        yield tuple(result_list)
    else:
        for i in listunique:
            if i.occurrences > 0:
                result_list[d]=i.value
                i.occurrences-=1
                for g in  perm_unique_helper(listunique,result_list,d-1):
                    yield g
                i.occurrences+=1




           
            
