#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy
from scipy.optimize import fsolve

'''
Provides a SU2 mesh of a cylinder of r=1 m in a circular domain
of r=40 m. This mesh can be used to numerically solve the Euler
equations as described in Hirsch, C., Numerical Computations of
Internal and External Flows, Second Edition, section 11.5 (p583)
'''

# number of nodes in the radial direction
nr = 33
# number of points in the circumferential direction
nt = 129
# cylinder radius
r0 = 1.
# first cell size in the radial direction
dr0 = 0.05
# farfield radius
r1 = 40.

if __name__ == "__main__":

    # partial sum of the n first element of a series with r>1
    def PartialSum(a0, r, n):
        return a0 * (r ** n - 1.) / (r - 1.)

    # determine geometrical ratio of cell growth in the radial direction
    f = lambda g: dr0 * (g ** (nr - 1.) - 1.) / (g - 1.) - (r1 - r0)
    g = fsolve(f, 1.2)[0]
    print("Geometrical ratio: %s" %(g))

    # define radius as a list
    radii = [1.0 + PartialSum(dr0, g, i) for i in range(nr)]

    # define azimuthal direction
    theta = numpy.linspace(0, 2 * numpy.pi, num=nt, endpoint=False)

    # compute a grid
    rr, tt = numpy.meshgrid(radii, theta, indexing="ij")
    # transform to cartesian
    xx, yy = rr * numpy.cos(tt), rr * numpy.sin(tt)

    # define the elements (ordering such that SU2 accepts it)
    nElems = (nr - 1) * nt
    elems = []
    for j in range(nt - 1):
        for i in range(nr - 1):
            firstNode = str(int((i + j * nr)))
            secondNode = str(int((i + j * nr + 1)))
            fourthNode = str(int((i + j * nr + nr)))
            thirdNode = str(int((i + j * nr + nr + 1)))
            elemID = str(int(i + j * (nr - 1)))
            elemDef = [
                "9 ",
                firstNode,
                " ",
                secondNode,
                " ",
                thirdNode,
                " ",
                fourthNode,
                " ",
                elemID
            ]
            elems.append(elemDef)

    # need to close the last element  (ordering such that SU2 accepts it)
    for i in range(nr - 1):
        j = nt - 1
        firstNode = str(int((i + j * nr)))
        secondNode = str(int((i + j * nr + 1)))
        fourthNode = str(int((i + 0 * nr)))
        thirdNode = str(int((i + 0 * nr + 1)))
        elemID = str(int(i + j * (nr - 1)))
        elemDef = [
            "9 ",
            firstNode,
            " ",
            secondNode,
            " ",
            thirdNode,
            " ",
            fourthNode,
            " ",
            elemID
        ]
        elems.append(elemDef)

    # define nodes
    nNodes = nr * nt
    nodes = []
    for j in range(nt):
        for i in range(nr):
            x = str(xx[i, j])
            y = str(yy[i, j])
            nodeID = str(int((i + j * nr)))
            nodeDef = [x, " ", y, " ", nodeID]
            nodes.append(nodeDef)

    # define cylinder skin bc
    cylinder = []
    for j in range(0, (nt - 1) * nr, nr):
        firstNodeID = str(j)
        secondNodeID = str(j + nr)
        faceType = '3'
        faceDef = [faceType, " ", firstNodeID, " ", secondNodeID]
        cylinder.append(faceDef)
    # close last face
    firstNodeID = str(j + nr)
    secondNodeID = str(0)
    faceType = '3'
    faceDef = [faceType, " ", firstNodeID, " ", secondNodeID]
    cylinder.append(faceDef)

    # define cylinder skin bc
    farfield = []
    for j in range(nr - 1, (nt - 1) * nr, nr):
        firstNodeID = str(j)
        secondNodeID = str(j + nr)
        faceType = '3'
        faceDef = [faceType, " ", firstNodeID, " ", secondNodeID]
        farfield.append(faceDef)
    # close last face
    firstNodeID = str(j + nr)
    secondNodeID = str(nr - 1)
    faceType = '3'
    faceDef = [faceType, " ", firstNodeID, " ", secondNodeID]
    farfield.append(faceDef)

    # write to file
    meshFile = open("cylinder.su2", "w")
    meshFile.write("NDIME= 2\n")
    meshFile.write("NELEM= %d\n"%(nElems))

    # write elements
    for element in elems:
        for e in element:
            meshFile.write(e)
        meshFile.write("\n")

    meshFile.write("NPOIN= %d\n"%(nNodes))

    # write nodes
    for node in nodes:
        for n in node:
            meshFile.write(n)
        meshFile.write("\n")

    meshFile.write("NMARK= 2\n")
    meshFile.write("MARKER_TAG= cylinder\n")
    meshFile.write("MARKER_ELEMS= %d\n"%(len(cylinder)))

    # write nodes
    for face in cylinder:
        for f in face:
            meshFile.write(f)
        meshFile.write("\n")

    meshFile.write("MARKER_TAG= farfield\n")
    meshFile.write("MARKER_ELEMS= %d\n"%(len(farfield)))

    # write nodes
    for face in farfield:
        for f in face:
            meshFile.write(f)
        meshFile.write("\n")

    meshFile.close()
