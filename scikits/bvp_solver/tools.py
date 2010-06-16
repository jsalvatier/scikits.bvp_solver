import numpy

def has_nans(arg):
    return not (numpy.isfinite(arg).all())

def preparg(arg):
    """
    prepares an argument the user passed for use
    """
    if not (arg is None):
        if not (isinstance(arg, numpy.ndarray )):
            try:
                # if the argument was not a numpy array, try to cast it to one
                arg = numpy.array(arg, dtype = numpy.float)

            except ValueError:
                raise ValueError("argument not castable to array with dtype float: " + str(arg))

        if has_nans(arg):
            raise ValueError("argument has NaNs, Infs or -Infs")
    return arg

def farg(arg):
    """
    prepares an array argument for passing to fortran code by filtering out None's
    """
    if arg is None:
        return numpy.array([])
    else:
        return arg

def showdiff(arg1, arg2):
    """
    show the non-equality of two values
    """
    return " <"+ str(arg1) + " != " + str(arg2) + "> "

def argShapeTest(arg, shape, argName = "", recommendation = ""):
    """
    Test whether an argument is a numpy array and whether it has the right shape, otherwise raise an error
    """
    if not isinstance(arg, numpy.ndarray):
        raise ValueError(argName + " is not a numpy array but a " + str(type(arg)))

    if has_nans(arg):
        raise ValueError(argName + " has NaNs, Infs or -Infs")

    if arg.shape != shape:
        raise ValueError(argName + " is not the right shape. " + recommendation +"\n"+ showdiff(arg.shape, shape))

def fromf(arg):
    if arg is not None:
        return arg.copy()
    else:
        return arg