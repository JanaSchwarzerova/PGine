import cgp


class Identity(cgp.OperatorNode):
    _arity = 1
    _def_output = "x_0"


class Sin(cgp.OperatorNode):
    _arity = 1
    _def_output = "np.sin(x_0)"


class Cos(cgp.OperatorNode):
    _arity = 1
    _def_output = "np.cos(x_0)"


class ProtectedSqrt(cgp.OperatorNode):
    _arity = 1
    _def_output = "np.sqrt(np.abs(x_0))"


class Pow(cgp.OperatorNode):
    _arity = 1
    _def_output = "np.power(x_0, 2)"


class Abs(cgp.OperatorNode):
    _arity = 1
    _def_output = "np.abs(x_0)"


class ProtectedDiv(cgp.OperatorNode):
    _arity = 2
    _def_output = "np.divide(x_0, x_1, out=np.zeros_like(x_0, dtype=np.float64), where=x_1!=0)"


class Minimum(cgp.OperatorNode):
    _arity = 2
    _def_output = "np.minimum(x_0, x_1)"


class Maximum(cgp.OperatorNode):
    _arity = 2
    _def_output = "np.maximum(x_0, x_1)"
