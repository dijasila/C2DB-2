"""Functionality to decorate plotting functionality in updated defaults."""


# Customized rc parameters for the c2db_2.0 paper
rcp = {}
# Fonts
rcp['font.family'] = 'serif'
rcp['font.serif'] = 'Times Roman'
rcp['font.monospace'] = 'Computer Modern Typewriter'
rcp['text.usetex'] = True
# Font sizes
rcp['font.size'] = 10
rcp['axes.labelsize'] = 10
rcp['legend.fontsize'] = 10
rcp['xtick.labelsize'] = 9
rcp['ytick.labelsize'] = 9


def plotter(rcp=rcp):
    """Decorator pattern for matplotlib plotters."""
    def decorator(plotter):
        return rcParamsDecorator(rcp, plotter)
    return decorator


class rcParamsDecorator:
    """Decorator class for matplotlib plotting methods."""
    def __init__(self, rcp, plotter):
        """Initiate rcParamsDecorator.

        Parameters
        ----------
        rcp : dict
            values to be updated from default rcParams
        plotter : function
            plotting method to decorate with custom rcParams
        """
        self._rcp = rcp
        self._plotter = plotter

    def __call__(self, *args, **kwargs):
        """Make plotter call with updated rcParams."""
        from matplotlib import rcParams
        rcParams.update(self._rcp)
        return self._plotter(*args, **kwargs)
