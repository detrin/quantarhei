# -*- coding: utf-8 -*-
from sympy.physics.quantum import Operator, Dagger
from sympy.physics.quantum.qexpr import QExpr
from sympy import I, conjugate
from sympy import S
from sympy import Function, Wild, Mul, Pow
from sympy import sympify

import numpy


class CumulantException(Exception):
    pass


"""

    Cumulant expression class


"""
class CumulantExpr(QExpr):


    def _eval_simplify(self, ratio, measure, rational, inverse, doit):
        return self._evaluate_second_order_rules()
        #return A._getExpr()


    def getOrder(self,n):
        if n == self._calculate_order(self):
            return self
        else:
            add = self.args[0]
            A = sympify(0)
            for aa in add.args:
                if self._calculate_order(aa) == 2 :
                    A = A + aa

            return CumulantExpr(A)


    def evaluate(self,large=False):
        expr = self.expand()
        expr = expr.getOrder(2)
        #print(expr)
        C = expr.simplify()
        #print(C)
        if large:
            C = C._make_positive(large)
            C = C._expand_in_large(large)
        expr = C._getExpr()
        return expr.simplify()

    """
    ***************************************************************************
    Helper routines of Cumulant Expr
    ***************************************************************************
    """





    def _evaluate_second_order_rules(self):
        """
        Evaluates second order terms in terms of the line shape function

        """
        a = Wild('a')
        b = Wild('b')
        w = Wild('w')
        t1 = Wild('t1')
        t2 = Wild('t2')

        """
        Combinations of hh

        hh_plus * hh_plus --->
        hh_minus * hh_minus --->
        hh_minus**2 --->
        hh_plus * hh_minus --->
        hh_minus * hh_plus --->
        """
        A = self.replace(w*hh_plus(a,t1)*hh_plus(b,t2), \
                         w*(gg(a,b,t1)-gg(a,b,t1-t2)+conjugate(gg(a,b,t2))))
        A = A.replace(w*hh_plus(a,t1)**2,w*(gg(a,a,t1)+conjugate(gg(a,a,t1))))
        A = A.replace(w*hh_minus(a,t1)*hh_minus(b,t2), \
                         w*(conjugate(gg(b,a,t1)) \
                           -conjugate(gg(b,a,t1-t2))+gg(b,a,t2)))
        A = A.replace(w*hh_minus(a,t1)**2,w*(gg(a,a,t1)+conjugate(gg(a,a,t1))))
        A = A.replace(w*hh_plus(a,t1)*hh_minus(b,t2), \
                         w*(-gg(a,b,t1) \
                           +gg(a,b,t1+t2)-gg(a,b,t2)))
        A = A.replace(w*hh_minus(a,t1)*hh_plus(b,t2), \
                         w*(conjugate(-gg(b,a,t1) \
                           +gg(b,a,t1+t2)-gg(b,a,t2))))

        """
        Replacement rules for ggs

        First daggered ggs
        (
        and that the normal ones

        gg_plus ---> gg
        \hat{g}^{(-)}_{ab}(t) ---> g^{*}_{ab}(t)

        """
        A = A.replace(w*Dagger(gg_plus(a,b,t1)),w*conjugate(gg(a,b,t1)))
        A = A.replace(w*Dagger(gg_minus(a,b,t1)),w*gg(a,b,t1))
        A = A.replace(w*gg_plus(a,b,t1),w*gg(a,b,t1))
        A = A.replace(w*gg_minus(a,b,t1),w*conjugate(gg(a,b,t1)))

        """
        Replacement rules for dVs and their combinations with hh
        """
        A = A.replace(w*dV(a,t1)*dV(b,t2),w*g2(a,b,t1-t2))
        A = A.replace(w*dV(a,t1)**2,w*g2(a,a,0))

        A = A.replace(w*dV(a,t1)*hh_plus(b,t2),w*(-g1(a,b,t1-t2)+g1(a,b,t1)))
        A = A.replace(w*dV(a,t1)*hh_minus(b,t2),w*(g1(a,b,t1+t2)-g1(a,b,t1)))
        A = A.replace(w*hh_plus(a,t1)*dV(b,t2), \
                      w*(g1(a,b,t1-t2)+conjugate(g1(b,a,t2))))
        A = A.replace(w*hh_minus(a,t1)*dV(b,t2),
                      w*conjugate(g1(b,a,t1+t2)-g1(b,a,t2)))
        #A = A.replace(w*hh_plus(a,t1)*dV(b,t2), \
        #              w*(g1(a,b,t1-t2)-conjugate(g1(a,b,t2))))
        #A = A.replace(w*hh_minus(a,t1)*dV(b,t2),
        #              w*conjugate(g1(a,b,t1+t2)-g1(a,b,t2)))

        return A

    def _make_positive(self,arg):
        a = Wild('a')
        b = Wild('b')
        w = Wild('w')
        t1 = Wild('t')
        A = self.replace(w*gg(a,b,-arg+t1),w*conjugate(gg(b,a,arg-t1)))
        A = A.replace(w*gg(a,b,-arg),w*conjugate(gg(b,a,arg)))
        A = A.replace(w*g1(a,b,-arg+t1),-w*conjugate(g1(b,a,arg-t1)))
        A = A.replace(w*g1(a,b,-arg),-w*conjugate(g1(b,a,arg)))
        A = A.replace(w*g2(a,b,-arg+t1),w*conjugate(g2(b,a,arg-t1)))
        A = A.replace(w*g2(a,b,-arg),w*conjugate(g2(b,a,arg)))
        return A

    def _expand_in_large(self,arg):
        a = Wild('a')
        b = Wild('b')
        w = Wild('w')
        t1 = Wild('t1')
        A = self.replace(w*gg(a,b,arg+t1),w*(gg(b,a,arg)+(dd(a,b)-I*lam(a,b))*t1))
        A = A.replace(w*g1(a,b,arg+t1),w*(dd(a,b)-I*lam(a,b)))
        A = A.replace(w*g2(a,b,arg+t1),0)
        return A

    def _leading_index(self,arg):
        a = Wild('a')
        w = Wild('w')
        t1 = Wild('t1')
        A = self.replace(w*gg(a,arg,t1),w*gg(arg,a,t1))
        A = A.replace(w*g1(a,arg,t1),w*g1(arg,a,t1))
        A = A.replace(w*g2(a,arg,t1),w*g2(arg,a,t1))
        A = A.replace(w*dd(a,arg),w*dd(arg,a))
        A = A.replace(w*lam(a,arg),w*lam(arg,a))
        return A

    def _eliminate_off_diagonal(self):
        return self

    def _getExpr(self):
        return self.args[0]

    def _calculate_order(self,expr):
        """
        Calculates a perturbation order of the expression expr

        """
        content = expr
        order = 0

        if  content.func is Mul:
            for aa in content.args:
                order += self._calculate_order(aa)
            return order

        elif content.func is Pow:
            sorder = self._calculate_order(content.args[0])
            return sorder*content.args[1]

        elif content.func in (hh_plus,hh_minus,gg_plus,gg_minus,dV):
            return content.order()

        elif content.func is Dagger:
            return self._calculate_order(content.args[0])
        else:
            return 0


"""
    Special Operators

"""
class Ugde(Operator):
    nargs = 2
    def _eval_rewrite_as_gg(self,a,t):
        return (1-I*hh_plus(a,t)-gg_plus(a,a,t))

class Uedg(Operator):
    nargs = 2
    def _eval_rewrite_as_gg(self,a,t):
        return (1+I*hh_plus(a,t)-Dagger(gg_plus(a,a,t)))

class Uged(Operator):
    nargs = 2
    def _eval_rewrite_as_gg(self,a,t):
        return (1+I*hh_minus(a,t)-gg_minus(a,a,t))

class Uegd(Operator):
    nargs = 2
    def _eval_rewrite_as_gg(self,a,t):
        return (1-I*hh_minus(a,t)-Dagger(gg_minus(a,a,t)))

class ExpdV(Operator):
    nargs = 3
    def _eval_rewrite_as_gg(self,a,t,x):
        return (1+x*dV(a,t)+x**2*dV(a,t)*dV(a,t))

class hh_plus(Operator):

    def order(self):
        return 1

class hh_minus(Operator):

    def order(self):
        return 1

class gg_plus(Operator):

    def order(self):
        return 2

class gg_minus(Operator):

    def order(self):
        return 2

class dV(Operator):

    def order(self):
        return 1

"""

    Lineshape function and related stuff

"""

class gg(Function):
    nargs = (1,2,3)

    @classmethod
    def eval(cls, a, b, x):
        if x.is_Number:
            if x is S.Zero:
                return S.Zero
        #if not (a is b):
        #    return S.Zero
        
""" First derivative of gg """
class g21(Function):
    nargs = (1,2,3)

    @classmethod
    def eval(cls,a,b,x):
        if x.is_Number:
            if x is S.Zero:
                return S.Zero
""" First derivative of gg """
class g1(Function):
    nargs = (1,2,3)

    @classmethod
    def eval(cls,a,b,x):
        if x.is_Number:
            if x is S.Zero:
                return S.Zero





""" Second derivative of gg """
class g2(Function):
    nargs = (1,2,3)


""" dephasing rate """
class dd(Function):
    nargs = (1,2)

    def _eval_is_real(self):
        return True

""" Reorganization energy """
class lam(Function):
    nargs = (1,2)

    def _eval_is_real(self):
        return True

"""High level cumulant evaluation function """
def evaluate_cumulant(cuml, positive_times = [], leading_index=None,
                      lang = None, arrays=None):
    """
    
    """

    from .lang import python_code
    from .lang import fortran_code

    A = cuml.rewrite(gg)
    expr = CumulantExpr(A)
    expr = expr.evaluate()
    
    for tt in positive_times:
        expr = CumulantExpr(expr)._make_positive(tt)    
        
    #a = leading_index[0]
    if leading_index is not None:
        D = expr._leading_index(leading_index)
        expr = D._getExpr()
        
    if lang is None:
        ss = expr
    elif lang == "Fortran":
        ss = fortran_code(expr.__str__())
    elif lang == "Python":
        ss = python_code(expr.__str__(),arrays=arrays)
    else:
        raise Exception("Unknown language")
    
    return ss


# def evaluate_cumulant(cuml, positive_times = [], leading_index=None,
#                       lang = None, arrays=None):
#     """
    
#     """

#     from quantarhei.symbolic.cumulant import CumulantExpr
#     from quantarhei.symbolic.lang import python_code
#     from quantarhei.symbolic.lang import fortran_code

#     A = cuml.rewrite(gg)
#     expr = CumulantExpr(A)
#     expr = expr.evaluate()
    
#     for tt in positive_times:
#         expr = CumulantExpr(expr)._make_positive(tt)    
        
#     #a = leading_index[0]
#     if leading_index is not None:
#         D = expr._leading_index(leading_index)
#         expr = D._getExpr()
        
#     if lang is None:
#         ss = expr
#     elif lang == "Fortran":
#         ss = fortran_code(expr.__str__())
#     elif lang == "Python":
#         ss = python_code(expr.__str__(),arrays=arrays)
#     else:
#         raise Exception("Unknown language")
    
#     return ss

import sympy as sp

def transform_to_Python_vec(expr, target_name="gf.gg", conj_func="numpy.conj"):
    """ Transforms the sympy cumulant evaluation into a matrix Python code
    
    Transform gg(int1, int2, time_vars) into target_name_t1_t2_t3["int1int2"]
    Convert conjugate(gg(...)) into conj_func(target_name_t1_t2_t3["int1int2"]).
    
    Parameters:
    -----------
    - expr: The SymPy expression to transform
    - target_name: The base string to replace `gg` with (default: "gf.gg")
    - conj_func: The function name for complex conjugation (default: "numpy.conj")
    
    """
    import re
    
    # Define conjugate function dynamically
    conjugate = sp.conjugate  # Allows flexible handling of conjugation

    # Handle conjugate(gg(...))
    if isinstance(expr, sp.Basic) and expr.func == conjugate:
        inner_expr = expr.args[0]  # Extract the function inside conjugate
        transformed_inner = transform_to_Python_vec(inner_expr, target_name, conj_func)  # Recursively transform
        return sp.Symbol(f"{conj_func}({transformed_inner})")  # Wrap with custom conjugate function

    # Handle gg(int1, int2, times)
    if isinstance(expr, sp.Basic) and expr.func == gg:
        args = expr.args  # Extract function arguments
        
        if len(args) != 3:  # Ensure it has exactly three arguments (two indices + one time argument)
            raise ValueError(f"Expected gg(int1, int2, times), but got {expr}")

        # First two arguments are integer indices (convert to strings)
        int_indices = [str(args[0]), str(args[1])]
        indices_str = f'"{"".join(int_indices)}"'  # Format as "int1int2"

        # Last argument is time variables (single or sum)
        time_vars = []
        time_arg = args[2]

        if isinstance(time_arg, sp.Add):  # If it's a sum, extract terms
            terms = [str(term) for term in time_arg.args]
        else:
            terms = [str(time_arg)]
        
        # Convert negative time variables by replacing '-' with 'm' without an underscore
        time_vars = [term.replace('-', 'm') if term.startswith('-') else term for term in terms]
        
        # Sort time variables based on the integer within their names
        time_vars.sort(key=lambda var: int(re.search(r'\d+', var).group()))
        
        # Construct the formatted time string
        time_str = "_".join(time_vars) if time_vars else "0"

        # Return the transformed expression
        return sp.Symbol(f'{target_name}_{time_str}[{indices_str}]')

    # Apply transformation recursively to sub-expressions
    return expr.replace(lambda subexpr: isinstance(subexpr, sp.Basic) and subexpr.func in {gg, conjugate}, 
                        lambda subexpr: transform_to_Python_vec(subexpr, target_name, conj_func))


from sympy import Add, Symbol, Function
#from quantarhei.symbolic.cumulant import gg

def transform_to_einsum(expr, participation_matrix=None, index=0):
    """ Transforms the cumulant result into a Numpy einsum form

    """

    einsum = Function("np.einsum")
    
    result_terms = []

    if participation_matrix is not None:
        if isinstance(participation_matrix, str):
            MM = participation_matrix
            use_exciton = True
        else:
            raise Exception("String expected here")
    else:
        use_exciton = False

    for term in Add.make_args(expr):
        coeff = 1
        func = term

        if isinstance(term, Mul):
            coeff, func = term.as_coeff_Mul()

        if isinstance(func, gg) and func.func.__name__ == "gg":
            a, b, t = func.args
            t_str = str(t).replace(" ", "")

            
            if use_exciton:
                gg_part = Symbol(f'gg[:,"{t_str}"]')
                M_part = Symbol(f'{MM}[{a},{b},:]')
                einsum_string = Symbol('"i,ij"')  # key fix: Symbol with quotes
                einsum_term = einsum(einsum_string, M_part, gg_part)
            else:
                einsum_term = Symbol(f'gg[{index},"{t_str}"]')
                
            transformed = coeff * einsum_term
            result_terms.append(transformed)
        else:
            result_terms.append(term)

    return Add(*result_terms)


class GFInitiator:
    """A helper class to precalculate the values of line shape functions 
    
    
    Parameters
    ----------

    t1s : array
        Array of t1 times

    t3s : array
        Arrax of t3 times

    gg : function
        Function that returns the lineshape function value at arbitrary time 

    """

    def __init__(self, t1s, t3s, gdict):

        self.t2 = -1.0
        self.t3 = -1.0
        self.t1s = t1s
        self.t3s = t3s

        #self.gg_t1 = dict()
        #self.gg_t2 = dict()
        #self.gg_t3 = dict()
        #self.gg_t1_t2 = dict()
        #self.gg_t1_t3 = dict() 
        #self.gg_t2_t3 = dict() 
        #self.gg_t1_t2_t3 = dict()
        
        self.set_lineshape_functions(gdict)


    def set_lineshape_functions(self, gdict):
        """Sets the lineshape functions.
        
        Parameters:
        gdict (dict): A dictionary containing lineshape functions.

        Raises:
        TypeError: If gdict is not a dictionary.
        """
        if not isinstance(gdict, dict):
            raise TypeError(f"Expected a dictionary, but got {type(gdict).__name__}")

        self.gg = gdict
        
        
    def set_t2(self, t2):
        """Sets the waiting time t2 and precalculates values of the response functions
        
        """

        # for every t2 we precalculate all necessary combinations of g(t) arguments
        if t2 != self.t2:
            self.t2 = t2

            self.gg_t1 = dict()
            self.gg_t2 = dict()
            self.gg_t3 = dict()
            self.gg_t1_t2 = dict()
            self.gg_t2_t3 = dict()
            self.gg_t1_t3 = dict()
            self.gg_t1_t2_t3 = dict()

            for key in self.gg:
                self.gg_t1[key] = self.eval_gg_t2_tn(key, 0.0, self.t1s, self.t3s, zero=3)
                self.gg_t2[key] = self.gg[key](t2)
                self.gg_t3[key] = self.eval_gg_t2_tn(key, 0.0, self.t1s, self.t3s, zero=1)
                self.gg_t1_t2[key] = self.eval_gg_t2_tn(key, t2, self.t1s, self.t3s, zero=3)
                self.gg_t2_t3[key] = self.eval_gg_t2_tn(key, t2, self.t1s, self.t3s, zero=1)
                self.gg_t1_t3[key] = self.eval_gg_t1_t2_t3(key, self.t1s, 0.0, self.t3s)
                self.gg_t1_t2_t3[key] = self.eval_gg_t1_t2_t3(key, self.t1s, t2, self.t3s)
            


    def eval_gg_t2_tn(self, key, t2, t1, t3, zero=0):
        """Lineshape function g(t1 + t2) or g(t2 + t3)

        It covers all the combinations g(t1), g(t3), g(t1 + t2) and g(t2 + t3)

        Parameters
        ----------

        t2 : float
            Time t2

        t1 : float array
            Array of t1 times


        """
        if zero == 0:
            raise Exception("Parameter zero has to be set to 1 or 3. Current value is "+str(zero))
        N1 = t1.shape[0]
        N3 = t3.shape[0]
        gg = self.gg[key]
        #ret = numpy.zeros((N1, N3), dtype=complex)
        if zero == 3:
            #for ii in range(N3):
            #    ret[:, ii] = gg(t2 + t1[:])
            A = gg(t2 + t1[:])
            ret = numpy.broadcast_to(A[:,numpy.newaxis],(N1, N3))
        elif zero == 1:
            #for ii in range(N1):
            #    ret[ii, :] = gg(t2 + t3[:])
            A = gg(t2 + t3[:])
            ret = numpy.broadcast_to(A[numpy.newaxis,:],(N1, N3))

        return ret


    def eval_gg_t1_t2_t3(self, key, t1, t2, t3):
        """Lineshape function g(t1 + t2 + t3)
        
        
        """
        N1 = t1.shape[0]
        N3 = t3.shape[0]
        gg = self.gg[key]
        #ret = numpy.zeros((N1, N3), dtype=complex)
        #for ii in range(N1):
        #    ret[ii,:] = gg(t2 + t1[ii] + t3)
            
        tt = numpy.broadcast_to(t1[:,numpy.newaxis],(N1,N3)) \
            +numpy.broadcast_to(t3[numpy.newaxis,:],(N1,N3))

        ret = gg(t2 + tt)

        return ret

