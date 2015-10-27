"""
author: Tan Wei Hou (https://github.com/tanweihou)

Summary:

Update schemes for Total Variation, Mean Curvature, Gaussian Curvature and Bernstein filters.
Please refer to https://github.com/YuanhaoGong/CurvatureFilter for mathematical details behind
these filters. Regular user should not use or modify them, use cf_filter instead.
"""

from __future__ import division
import numpy as np

def update_tv(inputimg, rowbegin, colbegin):
    """
    Update scheme for Total Variation filter.
    Helper function, regular user should not use or modify it.
    """
    inputimg_ = 5 * np.copy(inputimg[rowbegin:-1:2, colbegin:-1:2])
    inputimg_ij = inputimg[rowbegin:-1:2, colbegin:-1:2]

    d1 = (inputimg[rowbegin - 1:-2:2, colbegin:-1:2] + inputimg[rowbegin + 1::2, colbegin:-1:2]) \
        + inputimg[rowbegin:-1:2, colbegin + 1::2] \
        + inputimg[rowbegin - 1:-2:2, colbegin + 1::2] \
        + inputimg[rowbegin + 1::2, colbegin + 1::2] \
        - inputimg_

    d2 = (inputimg[rowbegin:-1:2, colbegin - 1:-2:2] + inputimg[rowbegin:-1:2, colbegin + 1::2]) \
        + inputimg[rowbegin - 1:-2:2, colbegin:-1:2] \
        + inputimg[rowbegin - 1:-2:2, colbegin - 1:-2:2] \
        + inputimg[rowbegin - 1:-2:2, colbegin + 1::2] \
        - inputimg_

    d3 = (inputimg[rowbegin - 1:-2:2, colbegin:-1:2] + inputimg[rowbegin + 1::2, colbegin:-1:2]) \
        + inputimg[rowbegin:-1:2, colbegin - 1:-2:2] \
        + inputimg[rowbegin - 1:-2:2, colbegin - 1:-2:2] \
        + inputimg[rowbegin + 1::2, colbegin - 1:-2:2] \
        - inputimg_

    d4 = (inputimg[rowbegin:-1:2, colbegin - 1:-2:2] + inputimg[rowbegin:-1:2, colbegin + 1::2]) \
        + inputimg[rowbegin + 1::2, colbegin:-1:2] \
        + inputimg[rowbegin + 1::2, colbegin - 1:-2:2] \
        + inputimg[rowbegin + 1::2, colbegin + 1::2] \
        - inputimg_

    d5 = inputimg[rowbegin - 1:-2:2, colbegin + 1::2] + inputimg[rowbegin:-1:2, colbegin + 1::2] \
        + inputimg[rowbegin + 1::2, colbegin:-1:2] + inputimg[rowbegin + 1::2, colbegin + 1::2] \
        + inputimg[rowbegin + 1::2, colbegin - 1:-2:2] \
        - inputimg_

    d6 = inputimg[rowbegin - 1:-2:2, colbegin + 1::2] + inputimg[rowbegin:-1:2, colbegin + 1::2] \
        + inputimg[rowbegin - 1:-2:2, colbegin:-1:2] + inputimg[rowbegin + 1::2, colbegin + 1::2] \
        + inputimg[rowbegin - 1:-2:2, colbegin - 1:-2:2] \
        - inputimg_

    d7 = inputimg[rowbegin - 1:-2:2, colbegin + 1::2] + inputimg[rowbegin:-1:2, colbegin - 1:-2:2] \
        + inputimg[rowbegin - 1:-2:2, colbegin:-1:2] + inputimg[rowbegin:-1:2, colbegin - 1:-2:2] \
        + inputimg[rowbegin + 1::2, colbegin - 1:-2:2] \
        - inputimg_

    d8 = inputimg[rowbegin + 1::2, colbegin + 1::2] + inputimg[rowbegin + 1::2, colbegin - 1:-2:2] \
        + inputimg[rowbegin + 1::2, colbegin:-1:2] + inputimg[rowbegin:-1:2, colbegin - 1:-2:2] \
        + inputimg[rowbegin - 1:-2:2, colbegin - 1:-2:2] \
        - inputimg_


    d = d1 * (np.abs(d1) <= np.abs(d2)) + d2 * (np.abs(d2) < np.abs(d1))
    d = d * (np.abs(d) <= np.abs(d3)) + d3 * (np.abs(d3) < np.abs(d))
    d = d * (np.abs(d) <= np.abs(d4)) + d4 * (np.abs(d4) < np.abs(d))
    d = d * (np.abs(d) <= np.abs(d4)) + d5 * (np.abs(d4) < np.abs(d))
    d = d * (np.abs(d) <= np.abs(d4)) + d6 * (np.abs(d4) < np.abs(d))
    d = d * (np.abs(d) <= np.abs(d4)) + d7 * (np.abs(d4) < np.abs(d))
    d = d * (np.abs(d) <= np.abs(d4)) + d8 * (np.abs(d4) < np.abs(d))

    d /= 5

    inputimg_ij[...] +=d

def update_mc(inputimg, rowbegin, colbegin):
    """
    Update scheme for Mean Curvature filter.
    Helper function, regular user should not use or modify it.
    """
    inputimg_ = 8 * np.copy(inputimg[rowbegin:-1:2, colbegin:-1:2])
    inputimg_ij = inputimg[rowbegin:-1:2, colbegin:-1:2]

    d1 = 2.5 * (inputimg[rowbegin - 1:-2:2, colbegin:-1:2] + inputimg[rowbegin + 1::2, colbegin:-1:2]) \
        + 5.0 * inputimg[rowbegin:-1:2, colbegin + 1::2] \
        - inputimg[rowbegin - 1:-2:2, colbegin + 1::2] \
        - inputimg[rowbegin + 1::2, colbegin + 1::2] \
        - inputimg_

    d2 = 2.5 * (inputimg[rowbegin:-1:2, colbegin - 1:-2:2] + inputimg[rowbegin:-1:2, colbegin + 1::2]) \
        + 5.0 * inputimg[rowbegin - 1:-2:2, colbegin:-1:2] \
        - inputimg[rowbegin - 1:-2:2, colbegin - 1:-2:2] \
        - inputimg[rowbegin - 1:-2:2, colbegin + 1::2] \
        - inputimg_

    d3 = 2.5 * (inputimg[rowbegin - 1:-2:2, colbegin:-1:2] + inputimg[rowbegin + 1::2, colbegin:-1:2]) \
        + 5.0 * inputimg[rowbegin:-1:2, colbegin - 1:-2:2] \
        - inputimg[rowbegin - 1:-2:2, colbegin - 1:-2:2] \
        - inputimg[rowbegin + 1::2, colbegin - 1:-2:2] \
        - inputimg_

    d4 = 2.5 * (inputimg[rowbegin:-1:2, colbegin - 1:-2:2] + inputimg[rowbegin:-1:2, colbegin + 1::2]) \
        + 5.0 * inputimg[rowbegin + 1::2, colbegin:-1:2] \
        - inputimg[rowbegin + 1::2, colbegin - 1:-2:2] \
        - inputimg[rowbegin + 1::2, colbegin + 1::2] \
        - inputimg_

    d = d1 * (np.abs(d1) <= np.abs(d2)) + d2 * (np.abs(d2) < np.abs(d1))
    d = d * (np.abs(d) <= np.abs(d3)) + d3 * (np.abs(d3) < np.abs(d))
    d = d * (np.abs(d) <= np.abs(d4)) + d4 * (np.abs(d4) < np.abs(d))

    d /= 8

    inputimg_ij[...] +=d

def update_gc(inputimg, rowbegin, colbegin):
    """
    Update scheme for Gaussian Curvature filter.
    Helper function, regular user should not use or modify it.
    """
    inputimg_ij = inputimg[rowbegin:-1:2, colbegin:-1:2]

    d1 = (inputimg[rowbegin - 1:-2:2, colbegin:-1:2] + inputimg[rowbegin + 1::2, colbegin:-1:2])/2.0 \
        - inputimg[rowbegin:-1:2, colbegin:-1:2]
    d2 = (inputimg[rowbegin:-1:2, colbegin - 1:-2:2] + inputimg[rowbegin:-1:2, colbegin + 1::2])/2.0 \
        - inputimg[rowbegin:-1:2, colbegin:-1:2]
    d3 = (inputimg[rowbegin - 1:-2:2, colbegin - 1:-2:2] + inputimg[rowbegin + 1::2, colbegin + 1::2])/2.0 \
        - inputimg[rowbegin:-1:2, colbegin:-1:2]
    d4 = (inputimg[rowbegin - 1:-2:2, colbegin + 1::2] + inputimg[rowbegin + 1::2, colbegin - 1:-2:2])/2.0 \
        - inputimg[rowbegin:-1:2, colbegin:-1:2]

    d5 = (inputimg[rowbegin - 1:-2:2, colbegin:-1:2] + inputimg[rowbegin:-1:2, colbegin - 1:-2:2] \
        + inputimg[rowbegin - 1:-2:2, colbegin - 1:-2:2])/3.0 \
        - inputimg[rowbegin:-1:2, colbegin:-1:2]
    d6 = (inputimg[rowbegin - 1:-2:2, colbegin:-1:2] + inputimg[rowbegin:-1:2, colbegin + 1::2] \
        + inputimg[rowbegin - 1:-2:2, colbegin + 1::2])/3.0 \
        - inputimg[rowbegin:-1:2, colbegin:-1:2]
    d7 = (inputimg[rowbegin:-1:2, colbegin - 1:-2:2] + inputimg[rowbegin + 1::2, colbegin:-1:2] \
        + inputimg[rowbegin + 1::2, colbegin - 1:-2:2])/3.0 \
        - inputimg[rowbegin:-1:2, colbegin:-1:2]
    d8 = (inputimg[rowbegin:-1:2, colbegin + 1::2] + inputimg[rowbegin + 1::2, colbegin:-1:2] \
        + inputimg[rowbegin + 1::2, colbegin + 1::2])/3.0 \
        - inputimg[rowbegin:-1:2, colbegin:-1:2]

    d = d1 * (np.abs(d1) <= np.abs(d2)) + d2 * (np.abs(d2) < np.abs(d1))
    d = d * (np.abs(d) <= np.abs(d3)) + d3 * (np.abs(d3) < np.abs(d))
    d = d * (np.abs(d) <= np.abs(d4)) + d4 * (np.abs(d4) < np.abs(d))
    d = d * (np.abs(d) <= np.abs(d5)) + d5 * (np.abs(d5) < np.abs(d))
    d = d * (np.abs(d) <= np.abs(d6)) + d6 * (np.abs(d6) < np.abs(d))
    d = d * (np.abs(d) <= np.abs(d7)) + d7 * (np.abs(d7) < np.abs(d))
    d = d * (np.abs(d) <= np.abs(d8)) + d8 * (np.abs(d8) < np.abs(d))

    inputimg_ij[...] +=d

def update_bernstein(inputimg, rowbegin, colbegin):
    """
    Update scheme for Bernstein filter.
    Helper function, regular user should not use or modify it.
    """
    inputimg_ij = inputimg[rowbegin:-1:2, colbegin:-1:2]

    d1 = (inputimg[rowbegin - 1:-2:2, colbegin:-1:2] + inputimg[rowbegin + 1::2, colbegin:-1:2])/2.0 \
        - inputimg[rowbegin:-1:2, colbegin:-1:2]
    d2 = (inputimg[rowbegin:-1:2, colbegin - 1:-2:2] + inputimg[rowbegin:-1:2, colbegin + 1::2])/2.0 \
        - inputimg[rowbegin:-1:2, colbegin:-1:2]

    d = d1 * (np.abs(d1) <= np.abs(d2)) + d2 * (np.abs(d2) < np.abs(d1))

    inputimg_ij[...] +=d
