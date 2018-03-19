#!/usr/bin/python3
#
#
# CMIP.Grid - Grid details for CMIP calculation
#
# import CMIP.Grid
#
# gr = CMIP.Grid()
#
# METHODS
#
# value=gr.readGrid(readGrid)
#
# value=gr.perfill(perfill)
#
# array=gr.int(array)
#
# array=gr.cen(array)
#
# array=gr.dim(array)
#
# txt=gr.asParamString(g) # g true for outer grid
#
# array=gr.origin();
#
# array=gr.size();
#
# value=gr.ngrid();
#
# index=gr.mapindex(i,j,k)
#
# ijk=gr.getijk(index)
#
# xyz=gr.getxyz(index)
#

class Grid(object):
    def __init__(self, pbfocus=0):
        self.pbfocus=pbfocus
        self.readgrid= 0
        self.perfill = 0.8
        self.int = [0.5, 0.5, 0.5]
        self.dim = [64, 64, 64]
        self.cen =  [0.0, 0.0, 0.0]
        self.needUpdate=False

    def __str__ (self):
        if self.pbfocus:
            g='0'
        else:
            g=''
        out = " {}{}={}\n".format("READGRID", g, self.readgrid)
        if self.readgrid == 0:
            out = out + " {}{}={},{}{}={},{}{}={}\n".format('CENX', g, self.cen[0], 'CENY', g, self.cen[1], 'CENZ', g, self.cen[2])
            out = out + " {}{}={},{}{}={},{}{}={}\n".format('DIMX', g, self.dim[0], 'DIMY', g, self.dim[1], 'DIMZ', g, self.dim[2])
            out = out + " {}{}={},{}{}={},{}{}={}\n".format('INTX', g, self.int[0], 'INTY', g, self.int[1], 'INTZ', g, self.int[2])
        elif self.readgrid == 3:
            out = out + " {}{}={},{}{}={},{}{}={}\n".format('CENX', g, self.cen[0], 'CENY', g, self.cen[1], 'CENZ', g, self.cen[2])
            out = out + " {}{}={},{}{}={},{}{}={}\n".format('DIMX', g, self.dim[0], 'DIMY', g, self.dim[1], 'DIMZ', g, self.dim[2])
            out = out + " {}{}={},{}{}={},{}{}={}\n".format('INTX', g, self.int[0], 'INTY', g, self.int[1], 'INTZ', g, self.int[2])
        elif self.readgrid == 2:
            out = out + " {}{}={}\n".format('PERFILL', g, self.perfill)
            out = out + " {}{}={},{}{}={},{}{}={}\n".format('INTX', g, self.int[0], 'INTY', g, self.int[1], 'INTZ', g, self.int[2])
        elif self.readgrid > 3:
            out = out + " {}{}={}\n".format('PERFILL', g, self.perfill)
            out = out + " {}{}={},{}{}={},{}{}={}\n".format('INTX', g, self.int[0], 'INTY', g, self.int[1], 'INTZ', g, self.int[2])
        return out

    def setreadGrid(self, r=''):
        if r:
            self.readgrid = r
            self.needUpdate = True
        return self.readgrid

    def setperfill(self, p=''):
        if p:
            self.perfill = p
            self.needUpdate = 1 
        return self.perfill

    def ngrid(self):
        return self._getval('ngrid')

    def size (self):
        return self._getval('size')

    def origin(self):
        return self._getval('origin')

    def mapindex(self, i, j, k):
        dimv = self.dim
        return (k-1) * (dimv[1] + 1) * (dimv[0] + 1) + (j-1) * (dimv[0] + 1) + i;

    def mapijk(self, idgrid):
        dimv = self.dim
        kg = 1 + int(idgrid / (dimv[0] + 1) / (dimv[1] + 1))
        j = idgrid-(kg-1) * (dimv[0] + 1) * (dimv[1] + 1)
        jg = 1 + int(j / (dimv[0] + 1));
        ig = j-(jg-1) * (dimv[0] + 1);
# afegir per evitar problemes de paret en i i j
        if ig == 0:
            jg = jg-1
            ig = dimv[0] + 1
        if jg == 0: 
            kg = kg-1
            jg = dimv[1] + 1
        return [ig, jg, kg]

    def xyz2grid (self, x, y, z):
        [x1, y1, z1] = self.origin
        [xs, ys, zs] = self.int
        return [(x-x1) / xs + 1, (z-z1) / zs + 1, (z-z1) / zs + 1]

    def grid2xyz (self, i, j, k):
        [x1, y1, z1] = self.origin
        [xs, ys, zs] = self.int
        return [(x1 + (i-1) * xs, y1 + (j-1) * ys, z1 + (k-1) * zs)]
	
    def update(self):
        if self.needUpdate:
            self._calcNgrid()
            self._calcOrigin()
            self._calcSize()
            self.needUpdate = False
        return self

####Private
    
    def _getval(self, val):
        self._update();
        return self.val

    def _calcNgrid(self):
        self.ngrid = (self.dim[0] + 1) * (self.dim[1] + 1) * (self.dim[2] + 1)
        return self

    def _calcOrigin(self):
        self.origin = [
            self.cen[0] - 0.5 * self.int[0] * self.dim[0],
            self.cen[1] - 0.5 * self.int[1] * self.dim[1],
            self.cen[2] - 0.5 * self.int[2] * self.dim[2]
        ]
        return self

    def _calcSize (self):
        self.size = [
            self.int[0] * self.dim[0],
            self.int[1] * self.dim[1],
            self.int[2] * self.dim[2]
        ]
        return self
