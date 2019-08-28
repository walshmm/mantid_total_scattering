import numpy as np
from scipy.constants import Avogadro

# -------------------------------------------------------------------------
# Volume in Beam


class Shape(object):
    def __init__(self):
        self.shape = None

    def getShape(self):
        return self.shape


class Cylinder(Shape):
    def __init__(self):
        self.shape = 'Cylinder'

    def volume(self, Radius=None, Height=None, **kwargs):
        return np.pi * Height * Radius * Radius


class Sphere(Shape):
    def __init__(self):
        self.shape = 'Sphere'

    def volume(self, Radius=None, **kwargs):
        return (4. / 3.) * np.pi * Radius * Radius * Radius


class GeometryFactory(object):

    @staticmethod
    def factory(Geometry):
        factory = {"Cylinder": Cylinder(),
                   "Sphere": Sphere()}
        return factory[Geometry["Shape"]]


def get_number_atoms(PackingFraction, MassDensity, MolecularMass, Geometry=None):
    # setup the geometry of the sample
    if Geometry is None:
        Geometry = dict()
    if "Shape" not in Geometry:
        Geometry["Shape"] = 'Cylinder'

    # get sample volume in container
    space = GeometryFactory.factory(Geometry)
    volume_in_beam = space.volume(**Geometry)

    number_density = PackingFraction * MassDensity / \
        MolecularMass * Avogadro  # atoms/cm^3
    natoms = number_density * volume_in_beam  # atoms
    return natoms
    