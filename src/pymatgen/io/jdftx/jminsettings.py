"""Store generic minimization settings read from a JDFTx out file.

This module contains the JMinSettings class for storing generic minimization
and mutants for storing specific minimization settings read from a JDFTx out
file.
"""

from dataclasses import dataclass


@dataclass
class JMinSettings:
    """Store generic minimization settings read from a JDFTx out file.

    Store generic minimization settings read from a JDFTx out file.
    """

    dirupdatescheme: str = None
    linminmethod: str = None
    niterations: int = None
    history: int = None
    knormthreshold: float = None
    energydiffthreshold: float = None
    nenergydiff: int = None
    alphatstart: float = None
    alphatmin: float = None
    updateteststepsize: bool = None
    alphatreducefactor: float = None
    alphatincreasefactor: float = None
    nalphaadjustmax: int = None
    wolfeenergy: float = None
    wolfegradient: float = None
    fdtest: bool = None
    maxthreshold: bool = None
    start_flag: str = None

    def __init__(
        self,
        dirupdatescheme: str = None,
        linminmethod: str = None,
        niterations: int = None,
        history: int = None,
        knormthreshold: float = None,
        energydiffthreshold: float = None,
        nenergydiff: int = None,
        alphatstart: float = None,
        alphatmin: float = None,
        updateteststepsize: bool = None,
        alphatreducefactor: float = None,
        alphatincreasefactor: float = None,
        nalphaadjustmax: int = None,
        wolfeenergy: float = None,
        wolfegradient: float = None,
        fdtest: bool = None,
        maxthreshold: bool = None,
    ) -> None:
        """Initialize a generic JMInSettings class.

        Parameters
        ----------
        dirupdatescheme : str
            The direction update scheme used in the minimization.
        linminmethod : str
            The line minimization method used in the minimization.
        niterations : int
            The number of iterations used in the minimization.
        history : int
            The number of previous steps used in the minimization.
        knormthreshold : float
            The threshold for the norm of the gradient.
        energydiffthreshold : float
            The threshold for the energy difference.
        nenergydiff : int
            The number of energy differences.
        alphatstart : float
            The starting step size.
        alphatmin : float
            The minimum step size.
        updateteststepsize : bool
            Whether to update the step size.
        alphatreducefactor : float
            The factor by which to reduce the step size.
        alphatincreasefactor : float
            The factor by which to increase the step size.
        nalphaadjustmax : int
            The maximum number of step size adjustments.
        wolfeenergy : float
            The energy Wolfe condition.
        wolfegradient : float
            The gradient Wolfe condition.
        fdtest : bool
            Whether to use finite difference testing.
        maxthreshold : bool
            Whether to use the maximum threshold.
        """
        # pre-commit was not a fan of the _assign_type method
        self.dirupdatescheme = None if dirupdatescheme is None else str(dirupdatescheme)
        self.linminmethod = None if linminmethod is None else str(linminmethod)
        self.niterations = None if niterations is None else int(niterations)
        self.history = None if history is None else int(history)
        self.knormthreshold = None if knormthreshold is None else float(knormthreshold)
        self.energydiffthreshold = (
            None if energydiffthreshold is None else float(energydiffthreshold)
        )
        self.nenergydiff = None if nenergydiff is None else int(nenergydiff)
        self.alphatstart = None if alphatstart is None else float(alphatstart)
        self.alphatmin = None if alphatmin is None else float(alphatmin)
        self.updateteststepsize = (
            None if updateteststepsize is None else bool(updateteststepsize)
        )
        self.alphatreducefactor = (
            None if alphatreducefactor is None else float(alphatreducefactor)
        )
        self.alphatincreasefactor = (
            None if alphatincreasefactor is None else float(alphatincreasefactor)
        )
        self.nalphaadjustmax = None if nalphaadjustmax is None else int(nalphaadjustmax)
        self.wolfeenergy = None if wolfeenergy is None else float(wolfeenergy)
        self.wolfegradient = None if wolfegradient is None else float(wolfegradient)
        self.fdtest = None if fdtest is None else bool(fdtest)
        self.maxthreshold = None if maxthreshold is None else bool(maxthreshold)


@dataclass
class JMinSettingsElectronic(JMinSettings):
    """JMInSettings mutant for electronic minimization settings.

    A class for storing electronic minimization settings read from a
    JDFTx out file.
    """

    start_flag: str = "electronic-minimize"

    def __init__(
        self,
        dirupdatescheme: str = None,
        linminmethod: str = None,
        niterations: int = None,
        history: int = None,
        knormthreshold: float = None,
        energydiffthreshold: float = None,
        nenergydiff: int = None,
        alphatstart: float = None,
        alphatmin: float = None,
        updateteststepsize: bool = None,
        alphatreducefactor: float = None,
        alphatincreasefactor: float = None,
        nalphaadjustmax: int = None,
        wolfeenergy: float = None,
        wolfegradient: float = None,
        fdtest: bool = None,
        maxthreshold: bool = None,
    ) -> None:
        super().__init__(
            dirupdatescheme=dirupdatescheme,
            linminmethod=linminmethod,
            niterations=niterations,
            history=history,
            knormthreshold=knormthreshold,
            energydiffthreshold=energydiffthreshold,
            nenergydiff=nenergydiff,
            alphatstart=alphatstart,
            alphatmin=alphatmin,
            updateteststepsize=updateteststepsize,
            alphatreducefactor=alphatreducefactor,
            alphatincreasefactor=alphatincreasefactor,
            nalphaadjustmax=nalphaadjustmax,
            wolfeenergy=wolfeenergy,
            wolfegradient=wolfegradient,
            fdtest=fdtest,
            maxthreshold=maxthreshold,
        )


@dataclass
class JMinSettingsFluid(JMinSettings):
    """JMInSettings mutant for fluid minimization settings.

    A class for storing fluid minimization settings read from a
    JDFTx out file.
    """

    start_flag: str = "fluid-minimize"

    def __init__(
        self,
        dirupdatescheme: str = None,
        linminmethod: str = None,
        niterations: int = None,
        history: int = None,
        knormthreshold: float = None,
        energydiffthreshold: float = None,
        nenergydiff: int = None,
        alphatstart: float = None,
        alphatmin: float = None,
        updateteststepsize: bool = None,
        alphatreducefactor: float = None,
        alphatincreasefactor: float = None,
        nalphaadjustmax: int = None,
        wolfeenergy: float = None,
        wolfegradient: float = None,
        fdtest: bool = None,
        maxthreshold: bool = None,
    ) -> None:
        super().__init__(
            dirupdatescheme=dirupdatescheme,
            linminmethod=linminmethod,
            niterations=niterations,
            history=history,
            knormthreshold=knormthreshold,
            energydiffthreshold=energydiffthreshold,
            nenergydiff=nenergydiff,
            alphatstart=alphatstart,
            alphatmin=alphatmin,
            updateteststepsize=updateteststepsize,
            alphatreducefactor=alphatreducefactor,
            alphatincreasefactor=alphatincreasefactor,
            nalphaadjustmax=nalphaadjustmax,
            wolfeenergy=wolfeenergy,
            wolfegradient=wolfegradient,
            fdtest=fdtest,
            maxthreshold=maxthreshold,
        )


@dataclass
class JMinSettingsLattice(JMinSettings):
    """JMInSettings mutant for lattice minimization settings.

    A class for storing lattice minimization settings read from a
    JDFTx out file.
    """

    start_flag: str = "lattice-minimize"

    def __init__(
        self,
        dirupdatescheme: str = None,
        linminmethod: str = None,
        niterations: int = None,
        history: int = None,
        knormthreshold: float = None,
        energydiffthreshold: float = None,
        nenergydiff: int = None,
        alphatstart: float = None,
        alphatmin: float = None,
        updateteststepsize: bool = None,
        alphatreducefactor: float = None,
        alphatincreasefactor: float = None,
        nalphaadjustmax: int = None,
        wolfeenergy: float = None,
        wolfegradient: float = None,
        fdtest: bool = None,
        maxthreshold: bool = None,
    ) -> None:
        super().__init__(
            dirupdatescheme=dirupdatescheme,
            linminmethod=linminmethod,
            niterations=niterations,
            history=history,
            knormthreshold=knormthreshold,
            energydiffthreshold=energydiffthreshold,
            nenergydiff=nenergydiff,
            alphatstart=alphatstart,
            alphatmin=alphatmin,
            updateteststepsize=updateteststepsize,
            alphatreducefactor=alphatreducefactor,
            alphatincreasefactor=alphatincreasefactor,
            nalphaadjustmax=nalphaadjustmax,
            wolfeenergy=wolfeenergy,
            wolfegradient=wolfegradient,
            fdtest=fdtest,
            maxthreshold=maxthreshold,
        )


@dataclass
class JMinSettingsIonic(JMinSettings):
    """JMInSettings mutant for ionic minimization settings.

    A class for storing ionic minimization settings read from a
    JDFTx out file.
    """

    start_flag: str = "ionic-minimize"

    def __init__(
        self,
        dirupdatescheme: str = None,
        linminmethod: str = None,
        niterations: int = None,
        history: int = None,
        knormthreshold: float = None,
        energydiffthreshold: float = None,
        nenergydiff: int = None,
        alphatstart: float = None,
        alphatmin: float = None,
        updateteststepsize: bool = None,
        alphatreducefactor: float = None,
        alphatincreasefactor: float = None,
        nalphaadjustmax: int = None,
        wolfeenergy: float = None,
        wolfegradient: float = None,
        fdtest: bool = None,
        maxthreshold: bool = None,
    ) -> None:
        super().__init__(
            dirupdatescheme=dirupdatescheme,
            linminmethod=linminmethod,
            niterations=niterations,
            history=history,
            knormthreshold=knormthreshold,
            energydiffthreshold=energydiffthreshold,
            nenergydiff=nenergydiff,
            alphatstart=alphatstart,
            alphatmin=alphatmin,
            updateteststepsize=updateteststepsize,
            alphatreducefactor=alphatreducefactor,
            alphatincreasefactor=alphatincreasefactor,
            nalphaadjustmax=nalphaadjustmax,
            wolfeenergy=wolfeenergy,
            wolfegradient=wolfegradient,
            fdtest=fdtest,
            maxthreshold=maxthreshold,
        )
