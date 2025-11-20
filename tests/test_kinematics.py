import pytest
import unittest
import numpy as np
import pandas as pd

import westcott
gw = westcott.Westcott()

class KinematicsTests(unittest.TestCase):

    __doc__="""Unit tests for methods belonging to the `Kinematics` class of 
    the `westcott_gfactors.py` module."""

    bad_args = [2, 3.14, "x"]

    def test_disiplay_thermal_properties_returns_NoneType(self):
        self.assertIsNone(gw.display_thermal_properties())

    def test_disiplay_thermal_properties_raises_TypeError_passing_args(self):
        for arg in KinematicsTests.bad_args:
            with self.assertRaises(TypeError):
                gw.display_thermal_properties(arg)

    def test_disiplay_constants_returns_NoneType(self):
        self.assertIsNone(gw.display_constants())

    def test_disiplay_constants_raises_TypeError_passing_args(self):
        for arg in KinematicsTests.bad_args:
            with self.assertRaises(TypeError):
                gw.display_constants(arg)

    
        
    
    
