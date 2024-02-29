"""
PGine: Software for PRS (Polygenic Risk Score) calculation in plants.
Functions for Cartesian Genetic Programming (CGP) algorithm implemented using HAL-CGP.

Author: Martin Hurta -  Faculty of Information Technology, Brno University of Technology, Czechia
Version: 1.0.1
Last update: 2024-02-29
"""
import cgp


class Identity(cgp.OperatorNode):
    """
    Class implementing simple identity function y = x_0.
    """
    _arity = 1
    _def_output = "x_0"


class Sin(cgp.OperatorNode):
    """
    Class implementing sine function y = sin(x_0).
    """
    _arity = 1
    _def_output = "np.sin(x_0)"


class Cos(cgp.OperatorNode):
    """
    Class implementing cosine function y = cos(x_0).
    """
    _arity = 1
    _def_output = "np.cos(x_0)"


class ProtectedSqrt(cgp.OperatorNode):
    """
    Class implementing protected variant of square root function.
    Negative value x_0 is first transformed into positive.
    Function is thus equivalent to y = abs(sqrt(x_0)).
    """
    _arity = 1
    _def_output = "np.sqrt(np.abs(x_0))"


class Pow(cgp.OperatorNode):
    """
    Class implementing power of two function y = pow(x_0, 2).
    """
    _arity = 1
    _def_output = "np.power(x_0, 2)"


class Abs(cgp.OperatorNode):
    """
    Class implementing absolute value function y = abs(x_0).
    """
    _arity = 1
    _def_output = "np.abs(x_0)"


class ProtectedDiv(cgp.OperatorNode):
    """
    Class protected division function. If the number is to be divided by 0, result will be constant 0.
    Function is thus equal to y = (x_0 / x_1 if x_1 != 0 else 0)
    """
    _arity = 2
    _def_output = "np.divide(x_0, x_1, out=np.zeros_like(x_0, dtype=np.float64), where=x_1!=0)"


class Minimum(cgp.OperatorNode):
    """
    Class implementing minimum function y = min(x_0, x_1).
    """
    _arity = 2
    _def_output = "np.minimum(x_0, x_1)"


class Maximum(cgp.OperatorNode):
    """
    Class implementing maximum function y = max(x_0, x_1).
    """
    _arity = 2
    _def_output = "np.maximum(x_0, x_1)"
