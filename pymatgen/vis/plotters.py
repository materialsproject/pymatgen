import collections

from pymatgen.util.plotting import pretty_plot

import palettable


class SpectrumPlotter(object):
    """
    Class for plotting Spectrum objects and subclasses. Note that the interface
    is extremely flexible given that there are many different ways in which
    people want to view spectra. The typical usage is::

        # Initializes plotter with some optional args. Defaults are usually
        # fine,
        plotter = SpectrumPlotter()

        # Adds a DOS (A kind of spectra) with a label.
        plotter.add_spectrum("Total DOS", dos)

        # Alternatively, you can add a dict of DOSs. This is the typical
        # form returned by CompleteDos.get_spd/element/others_dos().
        plotter.add_spectra({"dos1": dos1, "dos2": dos2})

    Args:
        xshift (float): A shift that is applied to the x values. This is
            commonly used to shift to an arbitrary zero. E.g., zeroing at the
            Fermi energy in DOS, or at the absorption edge in XAS spectra. The
            same xshift is applied to all spectra.
        yshift (float): A shift that is applied to the y values. This is
            commonly used to displace spectra for easier visualization.
            Successive spectra are applied successive shifts.
    """

    def __init__(self, xshift=0.0, yshift=0.0):
        self.xshift = xshift
        self.yshift = yshift
        # self.stack = stack
        # self.sigma = sigma
        self._spectra = collections.OrderedDict()

    def add_spectrum(self, label, spectrum):
        """
        Adds a Spectrum for plotting.

        Args:
            label:
                label for the Spectrum. Must be unique.
            spectrum:
                Dos object
        """
        self._spectra[label] = spectrum

    def add_spectra(self, spectra_dict, key_sort_func=None):
        """
        Add a dictionary of doses, with an optional sorting function for the
        keys.

        Args:
            dos_dict: dict of {label: Dos}
            key_sort_func: function used to sort the dos_dict keys.
        """
        if key_sort_func:
            keys = sorted(spectra_dict.keys(), key=key_sort_func)
        else:
            keys = spectra_dict.keys()
        for label in keys:
            self.add_spectra(label, spectra_dict[label])

    def get_plot(self, xlim=None, ylim=None):
        """
        Get a matplotlib plot showing the DOS.

        Args:
            xlim: Specifies the x-axis limits. Set to None for automatic
                determination.
            ylim: Specifies the y-axis limits.
        """

        ncolors = max(3, len(self._spectra))
        ncolors = min(9, ncolors)

        colors = palettable.colorbrewer.qualitative.Set1_9.mpl_colors

        plt = pretty_plot(12, 8)

        i = 0
        for key, sp in self._spectra.items():
            plt.plot(sp.x, sp.y + self.yshift * i, color=colors[i % ncolors],
                     label=str(key), linewidth=3)
            plt.xlabel(sp.XLABEL)
            plt.ylabel(sp.YLABEL)
            i += 1

        if xlim:
            plt.xlim(xlim)
        if ylim:
            plt.ylim(ylim)

        plt.legend()
        leg = plt.gca().get_legend()
        ltext = leg.get_texts()  # all the text.Text instance in the legend
        plt.setp(ltext, fontsize=30)
        plt.tight_layout()
        return plt

    def save_plot(self, filename, img_format="eps", **kwargs):
        """
        Save matplotlib plot to a file.

        Args:
            filename: Filename to write to.
            img_format: Image format to use. Defaults to EPS.
        """
        plt = self.get_plot(**kwargs)
        plt.savefig(filename, format=img_format)

    def show(self, **kwargs):
        """
        Show the plot using matplotlib.
        """
        plt = self.get_plot(**kwargs)
        plt.show()
