"""
Collection of helper functions
that are convenient in jupyter
notebooks
"""
import pandas as pd
import os
import logging
import matplotlib.pyplot as plt

PD_DEFAULT_ROWS = 6

pd.set_option('display.max_rows', PD_DEFAULT_ROWS)
plt.rcParams.update({'figure.max_open_warning': 0})

# filter numba performance warnings
import warnings
from numba import NumbaWarning
warnings.filterwarnings("ignore", category=NumbaWarning)


def print_dim(adata):
    """
    Prints an info log about the current dimensions of adata.
    """
    dims = adata.shape
    print("Current dimensions of adata: {} cells x {} genes.".format(dims[0], dims[1]))


def display(obj, n=PD_DEFAULT_ROWS, ncol=10, *args, **kwargs):
    """
    extension of the display function, allows to
    specify the number of rows/cols for pandas.
    """
    from IPython.display import display as ipy_display
    with pd.option_context('display.max_rows', n, 'display.max_columns', ncol):
        ipy_display(obj)


def setwd(max_levels=7):
    """
    set the current working directory
    to the first parent directory that
    contains the ".git" directory.

    Don't run if ran from within nextflow.

    Args:
        max_levels: If no .git directory is found after going `max_levels` up,
            an error is raised.

    """
    if "NXF_HOME" in os.environ:
        print("Working directory did not change because calling from nextflow. ")
        return

    cnt = 0
    while ".git" not in os.listdir(os.getcwd()):
        if cnt > max_levels:
            raise FileNotFoundError(".git not found in the top {}"
                                    "directories".format(max_levels))
        os.chdir("..")
        cnt += 1

    if cnt:
        print("Changed to directory {}".format("/".join([".."] * cnt)))

    print("Working directory is {}".format(os.path.abspath(os.getcwd())))


def fix_logging(sc_settings):
    """
    Workaround for rstudio/reticulate#386.

    Get all logger instances and add a custom PrintHandler,
    as writing to Stdout doesn't seem to work.

    Args:
      sc_settings: reference to the sc.settings instance
    """
    class PrintHandler(logging.Handler):
        def emit(self, record):
            print(self.format(record))

    for name, logger in logging.Logger.manager.loggerDict.items():
        try:
            logger.handlers = []
            logger.addHandler(PrintHandler())
        except AttributeError:
            pass

    sc_settings._root_logger.handlers = []
    sc_settings._root_logger.addHandler(PrintHandler())

    # more verbose logging (e.g. tell how many cells filtered out in each step)
    sc_settings.verbosity = 3
