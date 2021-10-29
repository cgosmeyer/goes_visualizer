
"""
Utility functions for visualizers.

Authors:

    C.M. Gosmeyer, B. Tan (2021)
"""

import datetime

def round_dt2minute(dt):
    """
    Parameters
    ----------
    dt : DateTime
        The datetime to be rounded to nearest minute.
    """
    return dt - datetime.timedelta(seconds=dt.second)

def round_dt2nearest(dt, interval):
    """
    Parameters
    ----------
    dt : DateTime
        The datetime to be rounded to nearest minute.
    interval : int
        Either 10 or 15, depending on mode.
    """
    dt += datetime.timedelta(minutes=interval/2)
    dt -= datetime.timedelta(minutes=dt.minute % interval,
                             seconds=dt.second,
                             microseconds=dt.microsecond)
    return dt