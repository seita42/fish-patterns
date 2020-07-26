import gzip
import colorsys

import numpy as np
import numexpr as ne
from scipy.spatial import cKDTree
from skimage import measure
from vapory import (Camera, Scene, LightSource, Background,
                    Sphere, Isosurface, Box, Texture,
                    Pigment, Finish, ContainedBy, Function)


def save_mb_obj(out_file, metaballs, vertices, vcolors=[], vnormals=[]):
    """Save OBJ (+ Metaballs) file

    Save OBJ (+ Metaballs) file

    Args:
        out_file:
        metaballs:
        vertices:
        vcolors:
        vnormals:

    Returns:
        bool:

    """
    if out_file.endswith(".gz"):
        f_out = gzip.open(out_file, 'wt')
    else:
        f_out = open(out_file, 'w')

    f_out.write("####\n")
    f_out.write("#\n")
    f_out.write("# Metaballs: {}\n".format(len(metaballs.mbs)))
    f_out.write("#\n")
    f_out.write("# mth {}\n".format(str(metaballs.mth)))
    for mb in metaballs.mbs:
        mbstr = "# " + str(mb) + "\n"
        f_out.write(mbstr)
    f_out.write("#\n")

    f_out.write("####\n")
    f_out.write("#\n")
    f_out.write("# Vertices: {}\n".format(len(vertices)))
    f_out.write("#\n")
    f_out.write("####\n")

    for vi, v in enumerate(vertices):
        vstr = "v {} {} {}".format(v[0], v[1], v[2])
        if len(vcolors) > 0:
            vc = vcolors[vi]
            vstr += " {} {} {}".format(vc[0], vc[1], vc[2])
        vstr += "\n"
        f_out.write(vstr)

    f_out.write("# {} vertices\n\n".format(len(vertices)))

    if len(vnormals) > 0:
        for vn in vnormals:
            vnstr = "vn {} {} {}\n".format(vn[0], vn[1], vn[2])
            f_out.write(vnstr)

        f_out.write("# {} normals\n\n".format(len(vnormals)))

    f_out.write("# End of File")
    f_out.close()

    return True


def load_mb_obj(in_file):
    """Load OBJ (+ Metaballs) file

    Load OBJ (+ Metaballs) file

    Args:
        in_file:

    Returns:
        metaballs (Metaballs):
        vertices:
        vcolors:
        vnormals:
        faces: (0-based)

    """

    mth = 0
    metaballs = []
    vertices = []
    vcolors = []
    vnormals = []
    faces = []

    if in_file.endswith(".gz"):
        f_in = gzip.open(in_file, 'rt')
    else:
        f_in = open(in_file, 'r')

    for line in f_in:
        vals = line.split()

        if len(vals) == 0:
            continue

        if vals[0] == "#":
            if (len(vals) > 2):
                if vals[1] == "mth":
                    mth = float(vals[2])

                if vals[1] == "mb":
                    mb = Metaball(float(vals[2]),
                                  float(vals[3]),
                                  float(vals[4]),
                                  float(vals[5]),
                                  float(vals[6]))
                    metaballs.append(mb)

        if vals[0] == "v":
            v = [float(x) for x in vals[1:4]]
            vertices.append(v)

            if len(vals) == 7:
                vc = [float(x) for x in vals[4:7]]
                vcolors.append(vc)

        if vals[0] == "vn":
            vn = [float(x) for x in vals[1:4]]
            vnormals.append(vn)

        if vals[0] == "f":
            fvi = []

            for f in vals[1:]:
                w = f.split("/")
                fvi.append(int(w[0]) - 1)

            faces.append(fvi)

    f_in.close()

    print("load "+in_file+": {:d} cells".format(len(vertices)))

    return (Metaballs(metaballs, mth),
            np.array(vertices),
            np.array(vcolors),
            np.array(vnormals),
            np.array(faces))


def save_rd_prms(out_file,
                 vs,
                 mbs,
                 A, B, C, D, E, F, G,
                 synUmax, synVmax, ucmax,
                 dt,
                 Du, Dv,
                 RR):
    """Save rd_prms file

    Save RDprms file

    Args:
        out_file:
        vs:
        mbs:
        A, B, C, D, E, F, G:
        synUmax, synVmax, ucmax:
        dt:
        Du, Dv:
        RR:

    Returns:
        bool:

    """

    np.savez_compressed(out_file,
                        vs=vs,
                        mbs=mbs,
                        A=A, B=B, C=C, D=D, E=E, F=F, G=G,
                        synUmax=synUmax, synVmax=synVmax, ucmax=ucmax,
                        dt=dt, Du=Du, Dv=Dv, RR=RR)

    return True


def save_rd_uv(out_file, ucs, vcs):
    """Save rd_uv file

    Save rd_uv file

    Args:
        out_file:
        ucs, vcs:

    Returns:
        bool:

    """

    np.savez_compressed(out_file, ucs=ucs, vcs=vcs)

    return True


def load_rd_prms(in_file):
    """Load rd_prms file

    Load rd_prms file

    Args:
        in_file:

    Returns:
        vs:
        mbs:
        A, B, C, D, E, F, G:
        synUmax, synVmax, ucmax:
        dt:
        Du:
        Dv:
        RR:

    """

    prms = np.load(in_file)
    return (prms['vs'],
            prms['mbs'].item(),
            prms['A'],
            prms['B'],
            prms['C'],
            prms['D'],
            prms['E'],
            prms['F'],
            prms['G'],
            prms['synUmax'],
            prms['synVmax'],
            prms['ucmax'],
            prms['dt'],
            prms['Du'],
            prms['Dv'],
            prms['RR'])


def load_rd_uv(in_file):
    """Load rd_uv file

    Load rd_uv file

    Args:
        in_file:

    Returns:
        ucs, vcs:

    """
    uv_data = np.load(in_file)
    return uv_data['ucs'], uv_data['vcs']


def load_rd_mb(fnbase, time_point=2000):
    """Load rd_mb file

    Load rd_mb file

    Args:
        fnbase:

    Returns:
        vs, ucs, A, C:

    """
    rd_prms = load_rd_prms(fnbase+"_prms.npz")
    vs, mbs = rd_prms[0:2]
    A, B, C, D, E, F, G, = rd_prms[2:9]
    synUmax, synVmax, ucmax, dt, Du, Dv, RR = rd_prms[9:]
    ucs, vcs = load_rd_uv(fnbase+"_{:05}.npz".format(time_point))
    ucs = ucs / ucmax
    ucs[ucs > 1.0] = 1.0

    return vs, ucs, A, C


class Metaball():
    """Metaball object class

    Metaball object class

    Attributes:
        x, y, z:
        s:
        a:

    """

    def __init__(self, x=0, y=0, z=0, s=1.0, a=1.0):
        self.x = x
        self.y = y
        self.z = z
        self.s = s
        self.a = a

    def __str__(self):
        return "mb {} {} {} {} {}".format(self.x,
                                          self.y,
                                          self.z,
                                          self.s,
                                          self.a)


class Metaballs():
    """Metaballs class (a set of Metaball object)

    Metaballs class (a set of Metaball object)

    Attributes:
        mbs (Metaball list):
        mth:

    """

    def __init__(self, mbs=[], mth=0.65):
        self.mbs = []
        self.mth = mth

        if len(mbs) > 0:
            self.mbs = mbs
            self.update()

    def append(self, metaball):
        """Append metaball

        Append a metaball

        Args:
            metaball:

        """

        self.mbs.append(metaball)

    def update(self):
        """Update metaballs

        Update metaballs

        """

        self._pre_calc_mb()

    def _pre_calc_mb(self):
        self.mx = np.array([mb.x for mb in self.mbs])
        self.my = np.array([mb.y for mb in self.mbs])
        self.mz = np.array([mb.z for mb in self.mbs])
        self.ms = np.array([mb.s for mb in self.mbs])
        self.ma = np.array([mb.a for mb in self.mbs])

        self.mxT = np.array([mb.x for mb in self.mbs])[:, np.newaxis]
        self.myT = np.array([mb.y for mb in self.mbs])[:, np.newaxis]
        self.mzT = np.array([mb.z for mb in self.mbs])[:, np.newaxis]
        self.msT = np.array([mb.s for mb in self.mbs])[:, np.newaxis]
        self.maT = np.array([mb.a for mb in self.mbs])[:, np.newaxis]

        xmin, xmax, ymin, ymax, zmin, zmax = self.get_min_max()
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.zmin = zmin
        self.zmax = zmax

        self.marching_cubes()

    def set_threshold(self, mth):
        self.mth = mth

    def to_povray_func(self):
        """Function string for Pov-Ray Isosurface

        Function string for Pov-Ray Isosurface
        func_str:
        mth - sum(mmi * exp(-((x-mxi)**2 + (y-ymi)**2 + (z-mzi)**2)/(2*msi**2))

        Returns:
            func_str:

        """

        func_str = "{}".format(self.mth)

        for mb in self.mbs:
            func_str += (" - {} "
                         "* exp(-(pow(x-{},2) + pow(y-{},2) + pow(z-{},2))"
                         "/ (2 * pow({},2)))"
                         .format(mb.a, mb.x, mb.y, mb.z, mb.s))

        return func_str

    def f(self, x, y, z):
        """Metaball value at (x, y, z)

        Metaball value at (x, y, z)
        metaball_value:
        mth - sum(mmi * exp(-((x-mxi)**2 + (y-ymi)**2 + (z-mzi)**2)/(2*msi**2))

        Args:
            x, y, z:

        Returns:
            metaball_value:

        """

        me = ne.evaluate("maT "
                         "* exp(-((x-mxT)**2 + (y-myT)**2 + (z-mzT)**2)"
                         "/ (2*msT**2))",
                         local_dict={'x': x,
                                     'y': y,
                                     'z': z,
                                     'mxT': self.mxT,
                                     'myT': self.myT,
                                     'mzT': self.mzT,
                                     'msT': self.msT,
                                     'maT': self.maT})

        mesum = np.sum(me, axis=0)
        return self.mth - mesum

    def cr(self, x, y, z, k):
        """Contribution ratios at (x, y, z)

        An array of contribution ratio of each metaball
        to the metaball_value at (x, y, z)

        Args:
            x, y, z:
            k:

        Returns:
            cr:

        """

        if k == 0:
            cr_eq = "maT*exp(-((x-mxT)**2+(y-myT)**2+(z-mzT)**2)/(2*msT**2))"
        elif k == 1:
            cr_eq = "1/sqrt((x-mxT)**2+(y-myT)**2+(z-mzT)**2)"
        elif k == 1.5:
            cr_eq = "1/sqrt((x-mxT)**2+(y-myT)**2+(z-mzT)**2)**1.5"
        elif k == 2:
            cr_eq = "1/((x-mxT)**2+(y-myT)**2+(z-mzT)**2)"

        me = ne.evaluate(cr_eq,
                         local_dict={'x': x,
                                     'y': y,
                                     'z': z,
                                     'mxT': self.mxT,
                                     'myT': self.myT,
                                     'mzT': self.mzT,
                                     'msT': self.msT,
                                     'maT': self.maT})

        mesum = np.sum(me, axis=0)
        cr = me/mesum
        return cr

    def grad_f(self, x, y, z):
        """Gradient of metaball_value at (x, y, z)

        Gradient of metaball_value at (x, y, z)

        Args:
            x, y, z:

        Returns:
            (dfdx, dfdy, dfdz):

        """

        str_dfdx = ("maT * (x-mxT)"
                    "* exp(-((x-mxT)**2+(y-myT)**2+(z-mzT)**2)/(2*msT**2))"
                    "/ (msT**2)")
        dfdx = ne.evaluate(str_dfdx,
                           local_dict={'x': x,
                                       'y': y,
                                       'z': z,
                                       'mxT': self.mxT,
                                       'myT': self.myT,
                                       'mzT': self.mzT,
                                       'msT': self.msT,
                                       'maT': self.maT})
        str_dfdy = ("maT * (y-myT)"
                    "* exp(-((x-mxT)**2+(y-myT)**2+(z-mzT)**2)/(2*msT**2))"
                    "/ (msT**2)")
        dfdy = ne.evaluate(str_dfdy,
                           local_dict={'x': x,
                                       'y': y,
                                       'z': z,
                                       'mxT': self.mxT,
                                       'myT': self.myT,
                                       'mzT': self.mzT,
                                       'msT': self.msT,
                                       'maT': self.maT})
        str_dfdz = ("maT * (z-mzT)"
                    "* exp(-((x-mxT)**2+(y-myT)**2+(z-mzT)**2)/(2*msT**2))"
                    "/ (msT**2)")
        dfdz = ne.evaluate(str_dfdz,
                           local_dict={'x': x,
                                       'y': y,
                                       'z': z,
                                       'mxT': self.mxT,
                                       'myT': self.myT,
                                       'mzT': self.mzT,
                                       'msT': self.msT,
                                       'maT': self.maT})

        dfdx_sum = np.sum(dfdx, axis=0)
        dfdy_sum = np.sum(dfdy, axis=0)
        dfdz_sum = np.sum(dfdz, axis=0)

        return dfdx_sum, dfdy_sum, dfdz_sum

    def f_v(self, v):
        """Metaball value at v

        Metaball value at v
        metaball_value:
        mth - sum(mmi*exp(-((x-mxi)**2+(y-ymi)**2+(z-mzi)**2)/(2*msi**2))

        Args:
            v:

        Returns:
            metaball_value:

        """

        return self.f(v[:, 0], v[:, 1], v[:, 2])

    def cr_v(self, v, k):
        """Contribution ratios at v

        An array of contribution ratio of each metaball
        to the metaball_value at v

        Args:
            v:
            k:

        Returns:
            cr:

        """

        return self.cr(v[:, 0], v[:, 1], v[:, 2], k).T

    def grad_f_v(self, v):
        """Gradient of metaball_value at v

        Gradient of metaball_value at v

        Args:
            v:

        Returns:
            grad:

        """

        dfdx, dfdy, dfdz = self.grad_f(v[:, 0], v[:, 1], v[:, 2])
        return np.c_[dfdx, dfdy, dfdz]

    def normal_f_v(self, v):
        """Normal vector at v

        Normal vector at v

        Args:
            v:

        Returns:
            normal_vec:

        """

        gv = self.grad_f_v(v)
        gvn = np.linalg.norm(gv, axis=1)

        return (gv / gvn[:, np.newaxis])

    def get_min_max(self):
        """Loose bounding box for metaballs

        Loose bounding box for metaballs

        Returns:
            xmin, xmax, ymin, ymax, zmin, zmax:

        """

        mr = np.sqrt(2 * np.log(1/self.mth)) * self.ms
        mr[:] = np.max(mr)

        mxmin = self.mx - mr
        mxmax = self.mx + mr
        mymin = self.my - mr
        mymax = self.my + mr
        mzmin = self.mz - mr
        mzmax = self.mz + mr

        mb_xmin_idx = np.argmin(mxmin[self.ma > 0])
        mb_xmax_idx = np.argmax(mxmax[self.ma > 0])
        mb_ymin_idx = np.argmin(mymin[self.ma > 0])
        mb_ymax_idx = np.argmax(mymax[self.ma > 0])
        mb_zmin_idx = np.argmin(mzmin[self.ma > 0])
        mb_zmax_idx = np.argmax(mzmax[self.ma > 0])

        xmin0 = self.mx[mb_xmin_idx] - mr[mb_xmin_idx]
        xmax0 = self.mx[mb_xmax_idx] + mr[mb_xmax_idx]
        ymin0 = self.my[mb_ymin_idx] - mr[mb_ymin_idx]
        ymax0 = self.my[mb_ymax_idx] + mr[mb_ymax_idx]
        zmin0 = self.mz[mb_zmin_idx] - mr[mb_zmin_idx]
        zmax0 = self.mz[mb_zmax_idx] + mr[mb_zmax_idx]

        xmin = xmin0 - (xmax0 - xmin0) * 0.25
        xmax = xmax0 + (xmax0 - xmin0) * 0.25
        ymin = ymin0 - (ymax0 - ymin0) * 0.25
        ymax = ymax0 + (ymax0 - ymin0) * 0.25
        zmin = zmin0 - (zmax0 - zmin0) * 0.25
        zmax = zmax0 + (zmax0 - zmin0) * 0.25

        return xmin, xmax, ymin, ymax, zmin, zmax

    def set_min_max(self, xmin, xmax, ymin, ymax, zmin, zmax):
        """Set loose bounding box for metaballs

        Set loose bounding box for metaballs

        Args:
            xmin, xmax, ymin, ymax, zmin, zmax:

        """

        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.zmin = zmin
        self.zmax = zmax

    def get_mb_meshgrid(self, spc=0.02):
        """Meshgrid of loose bounding box for metaballs

        Meshgrid of loose bounding box for metaballs

        Args:
            spc:

        Returns:
            mb_meshgrid:
            xyz_spc:

        """

        kx = np.int(np.round((self.xmax - self.xmin) / spc))
        ky = np.int(np.round((self.ymax - self.ymin) / spc))
        kz = np.int(np.round((self.zmax - self.zmin) / spc))

        X, Xspc = np.linspace(self.xmin, self.xmax, kx+1, retstep=True)
        Y, Yspc = np.linspace(self.ymin, self.ymax, ky+1, retstep=True)
        Z, Zspc = np.linspace(self.zmin, self.zmax, kz+1, retstep=True)

        xyz_spc = (Xspc, Yspc, Zspc)

        XX, YY, ZZ = np.meshgrid(X, Y, Z, indexing='ij')
        XYZs = np.c_[(XX.ravel(), YY.ravel(), ZZ.ravel())]

        mb_meshgrid = self.f_v(XYZs).reshape(kx+1, ky+1, kz+1)

        return mb_meshgrid, xyz_spc

    def marching_cubes(self, spc=0.02):
        """Set initial vertices on metaballs using marching cubes

        Set initial vertices on metaballs using marching cubes

        Args:
            spc:

        """

        mb_meshgrid, xyz_spc = self.get_mb_meshgrid(spc)

        verts, faces, normals, values = measure.marching_cubes(
            mb_meshgrid,
            level=0.0,
            spacing=xyz_spc,
            gradient_direction='ascent',
            step_size=1)

        verts += np.c_[self.xmin, self.ymin, self.zmin]

        self.verts = verts
        self.faces = faces
        self.normals = normals
        self.values = values
        self.sa = measure.mesh_surface_area(verts, faces)

    def scatter_vertices(self, n=20000):
        """Scatter n vertices on metaballs

        Scatter n vertices on metaballs

        Args:
            n:

        Returns:
            vs:

        """

        vs = self.verts[np.random.choice(len(self.verts), n, replace=False)]

        return vs

    def effective_radius(self, n):
        """Effective radius

        Effective radius for repulsive force

        Args:
            n:

        Returns:
            er:

        """

        er2 = 5.0 * self.sa / n
        er = np.sqrt(er2)

        return er

    def stick_to_mb(self, vs0, max_ite=10, eps=1.0e-10):
        """Move a point onto the surface of metaball

        Move a point onto the surface of metaball

        Args:
            vs0:
            max_ite:
            eps:

        Returns:
            vs:

        """

        n = len(vs0)

        normal_vs0 = self.normal_f_v(vs0)

        ts = np.zeros(n)
        delta_ts = np.ones(n)

        vs = vs0

        for i in range(max_ite):
            delta_ts = 1.0 * self.f_v(vs) / np.einsum('ij,ij->i',
                                                      self.grad_f_v(vs),
                                                      normal_vs0)
            ts -= delta_ts

            vs = vs0 + ts[:, np.newaxis] * normal_vs0

            if (np.abs(delta_ts).max() < eps):
                break

        return vs


def get_nei_force(vs, mbs):
    """Get neighbor points and repulsive forces

    Get neighbor points and repulsive forces

    Args:
        vs: vertices
        mbs: Metaballs

    Returns:
        nei:
        dnei:
        fvnei:
        fv:
        dnei6mean:
        dnei6std:

    """

    n = len(vs)
    er = mbs.effective_radius(n)
    f_er = er * 1.0

    # repulsive force = alpha * exp(- ((dnei)**2) / (2 * sigma**2))
    alpha = 1.0
    sigma = 0.3 * f_er

    # maximum number of neighbor points to be obtained
    maxnei = 50
    # threads used by cKDTree query
    q_threads = 6

    nei = np.empty((n, maxnei+1), dtype=int)
    dnei = np.empty((n, maxnei+1), dtype=float)
    fnei = np.empty((n, maxnei), dtype=float)
    vnei = np.empty((n, maxnei, 3), dtype=float)
    fvnei = np.empty((n, maxnei, 3), dtype=float)
    fvneisum = np.empty((n, 3), dtype=float)
    Mv = np.empty((n, 3, 3), dtype=float)
    fv = np.empty((n, 3), dtype=float)

    tree = cKDTree(vs)

    # get neighbor points
    dnei0, nei0 = tree.query(vs, k=maxnei+1, n_jobs=q_threads)

    # remove first hit (self)
    dnei = dnei0[:, 1:]
    nei = nei0[:, 1:]

    # vector (neighbor -> self)
    vnei = np.einsum('ij,k->ikj', vs, np.ones(maxnei)) - vs[nei]

    # repulsive force from neighbor points
    # fnei = alpha * np.exp((-dnei**2)/(2*sigma**2))
    fnei[dnei < f_er] = alpha * np.exp((-dnei[dnei < f_er]**2)/(2*sigma**2))
    fnei[dnei >= f_er] = 0

    # repulsive force (vector)
    fvnei = np.einsum('ij,ijk->ijk', fnei/np.linalg.norm(vnei, axis=2), vnei)

    # sum of repulsive force (vector)
    fvneisum = fvnei.sum(axis=1)

    # projection matrix to the tangent plane of metaball at v
    Mv = (np.einsum('ij,k->kji', np.eye(3), np.ones(n))
          - np.einsum('i,ij,ik->ijk',
                      np.linalg.norm(mbs.grad_f_v(vs), axis=1)**2,
                      mbs.grad_f_v(vs),
                      mbs.grad_f_v(vs)))

    # tangential component of sum of repulsive forces
    fv = np.einsum('ijk,ik->ij', Mv, fvneisum)

    # mean/std distance to the nearby 6 points (used for end condition)
    dnei6mean = dnei[:, 0:6].mean()
    dnei6std = dnei[:, 0:6].std()

    return nei, dnei, fvnei, fv, dnei6mean, dnei6std


def render_povray(vs, rotx=0, roty=0, rotz=0,
                  width=400, height=400, angle=14, antialiasing=0.001):
    """Render vertices using Pov-Ray (Vapory)

    Render vertices using Pov-Ray (Vapory)

    Args:
        vs: vertices
        rotx, roty, rotz: rotation angle
        width, height:
        angle: camera angle

    Returns:
        rendered_scene:

    """

    rot1 = [rotx, 0, 0]
    rot2 = [0, roty, 0]
    rot3 = [0, 0, rotz]

    camera = Camera('location', [0, 0, -25],
                    'look_at', [0, 0, 0],
                    'angle', angle,
                    'right x*image_width/image_height')

    light = LightSource([-3, 2, -6], 'color', [1.0, 1.0, 1.0], 'parallel')
    light2 = LightSource([2, -2, -6], 'color', [0.6, 0.6, 0.6], 'parallel')
    background = Background('color', [1, 1, 1])

    spheres = [Sphere(v, 0.05,
                      Finish('ambient', 0.2, 'diffuse', 0.8, 'phong', 1.0),
                      Texture(Pigment('color', [1.0, 1.0, 1.0])),
                      'rotate', rot1,
                      'rotate', rot2,
                      'rotate', rot3) for v in vs]

    objects = [light, light2, background] + spheres

    scene = Scene(camera, objects=objects)

    return scene.render('ipython',
                        width=width, height=height,
                        antialiasing=antialiasing)


def render_povray_mb(mbs, rotx=0, roty=0, rotz=0,
                     width=400, height=400, angle=14):
    """Render metaballs using Pov-Ray (Vapory)

    Render metaballs using Pov-Ray (Vapory)

    Args:
        mbs: Metaballs
        width, height:

    Returns:
        rendered_scene:

    """

    rot1 = [rotx, 0, 0]
    rot2 = [0, roty, 0]
    rot3 = [0, 0, rotz]

    camera = Camera('location', [0, 0, -25],
                    'look_at', [0, 0, 0],
                    'angle', angle,
                    'right x*image_width/image_height')

    light = LightSource([-3, 2, -6], 'color', [1.0, 1.0, 1.0], 'parallel')
    # light2 = LightSource([2, -2, -6], 'color', [0.6, 0.6, 0.6], 'parallel')
    background = Background('color', [1, 1, 1, 1])

    mbs_function = mbs.to_povray_func()

    isosurface = Isosurface(Function(mbs_function),
                            ContainedBy(Box(-5, 5)),
                            'max_gradient', 1.8,
                            Pigment('color', [1.0, 0.15, 0.3]),
                            Finish('phong', 0.7,
                                   'specular', 0.2,
                                   'diffuse', 0.9,
                                   'ambient', 0.1),
                            'rotate', rot1,
                            'rotate', rot2,
                            'rotate', rot3,
                            'translate', [0, 0, 0],
                            'no_shadow')

    objects = [light, background] + [isosurface]

    scene = Scene(camera, objects=objects)

    return scene.render('ipython', width=width, height=height)


def make_mb_obj(file, mbs, mva=4e-4,
                do_ite=True, dt=1.0e-2, max_ite=5000, cv_end=0.08):
    """Make metaball_OBJ

    Make metaball_OBJ (consisting of evenly distributed vertices on metaball)

    Args:
        file: output file
        mbs: Metaballs
        mva: mean area per vertex
        do_ite: if True, iteration will be performed
        dt: iteration time step
        max_ite: maximum number of iterations
        cv_end: coefficient of variation threshold for end condition

    Returns:
        vs: vertices

    """

    # scatter initial vertices and calculate effective radius
    n = np.int(np.round(mbs.sa / mva))
    print("scatter vertices: ", n)
    vs = mbs.scatter_vertices(n)

    # iteration using repulsive force
    if (do_ite):
        for i in range(max_ite):
            nei, dnei, fvnei, fv, dnei6mean, dnei6std = get_nei_force(vs, mbs)
            print("[ite: {}] dnei6mean={}, dnei6std={}, CV={}"
                  .format(i, dnei6mean, dnei6std, dnei6std/dnei6mean))
            fv_norm = np.linalg.norm(fv, axis=1)
            print("fv_norm_mean={}, fv_norm_max={}"
                  .format(np.mean(fv_norm), np.max(fv_norm)))

            if dnei6std/dnei6mean < cv_end:
                break

            vs += fv * dt
            vs = mbs.stick_to_mb(vs)

    save_mb_obj(file, mbs, vs)

    return vs


def rdmb(vs, mbs, out_file_base, prms, max_ite=2000):
    """Reaction-diffusion on metaball

    Reaction-diffusion on metaball

    Args:
        vs: vertices
        mbs: Metaballs
        mb_obj_file: mb_obj_file
        out_file_base: output file basename
        prms: rd parameters A..G
        max_ite: maximum number of iterations

    """

    nei, dnei, fvnei, fv, dnei6mean, dnei6std = get_nei_force(vs, mbs)

    dneiN = dnei / dnei6mean

    n = len(vs)
    max_nei = len(nei[0])

    er = mbs.effective_radius(n)

    A, B, C, D, E, F, G = prms

    synUmax = 0.23
    synVmax = 0.50
    ucmax = 6.0

    dt = 1.0e-2

    Du = 0.5
    Dv = 10.0

    # RR = 30
    RR = 80

    u = np.random.rand(n) * ucmax
    v = np.random.rand(n) * ucmax

    rea_u = np.zeros(n)
    rea_v = np.zeros(n)
    syn_u = np.zeros(n)
    syn_v = np.zeros(n)

    Ru = Du / (dneiN**2)
    Rv = Dv / (dneiN**2)

    Ru[dnei > er] = 0
    Rv[dnei > er] = 0

    save_rd_prms(out_file_base+"_prms.npz",
                 vs, mbs,
                 A, B, C, D, E, F, G,
                 synUmax, synVmax, ucmax, dt, Du, Dv, RR)
    save_rd_uv(out_file_base+"_{:05}.npz".format(0), u, v)

    for ite in range(max_ite):
        syn_u = A * u - B * v + C
        syn_v = E * u - F

        syn_u[syn_u < 0] = 0
        syn_u[syn_u > synUmax] = synUmax
        syn_v[syn_v < 0] = 0
        syn_v[syn_v > synVmax] = synVmax

        rea_u = syn_u - D * u
        rea_v = syn_v - G * v

        uu = (Ru * (u[nei] - np.einsum('i,j->ij', u,
                                       np.ones(max_nei)))).sum(axis=1)
        vv = (Rv * (v[nei] - np.einsum('i,j->ij', v,
                                       np.ones(max_nei)))).sum(axis=1)

        u += (RR * rea_u + uu) * dt
        v += (RR * rea_v + vv) * dt

        if ((ite+1) % 500) == 0:
            print("[ite: {}]".format(ite+1))

    fname = out_file_base + "_{:05}.npz".format(max_ite)
    save_rd_uv(fname, u, v)

    return True


def rdmb_uni_AC(mb_obj_file, out_file_base, A=0.08, C=0.15, max_ite=2000):
    """Reaction-diffusion on metaball (uniform A, C)

    Reaction-diffusion on metaball (uniform A, C)

    Args:
        mb_obj_file: metaball OBJ file
        out_file_base: output file basename
        A: parameter A
        C: parameter C
        max_ite: maximum number of iterations

    """

    mbs, vs, vc, vn, fs = load_mb_obj(mb_obj_file)

    A = A
    B = 0.08
    C = C
    D = 0.03
    E = 0.10
    F = 0.12
    G = 0.06

    prms = (A, B, C, D, E, F, G)

    rdmb(vs, mbs, out_file_base, prms, max_ite)

    return True


def rdmb_grad_AC(mb_obj_file,
                 out_file_base,
                 px0=-1.0, px1=1.0,
                 pa0=0.08, pa1=0.08,
                 pc0=0.03, pc1=0.27,
                 max_ite=2000):
    """Reaction-diffusion on metaball (gradual A, C)

    Reaction-diffusion on metaball (gradual A, C)

    Args:
        mb_obj_file: metaball OBJ file
        out_file_base: output file basename
        px0, px1: gradation start/end edges
        pa0, pa1: A values at edges
        pc0, pc1: C values at edges
        max_ite: maximum number of iterations

    """

    mbs, vs, vc, vn, fs = load_mb_obj(mb_obj_file)
    n = len(vs)

    # make A gradient
    # A = pa * x + pb
    pa = (pa1 - pa0)/(px1 - px0)
    pb = pa0 - pa * px0

    A = np.ones(n) * 0.15
    A[:] = vs[:, 0] * pa + pb      # A = pa * x + pb
    A[vs[:, 0] <= px0] = pa0
    A[vs[:, 0] >= px1] = pa1

    B = 0.08

    # make C gradient
    # C = pa * x + pb
    pa = (pc1 - pc0)/(px1 - px0)
    pb = pc0 - pa * px0

    C = np.ones(n) * 0.15
    C[:] = vs[:, 0] * pa + pb      # C = pa * x + pb
    C[vs[:, 0] <= px0] = pc0
    C[vs[:, 0] >= px1] = pc1

    D = 0.03
    E = 0.10
    F = 0.12
    G = 0.06

    prms = (A, B, C, D, E, F, G)

    rdmb(vs, mbs, out_file_base, prms, max_ite)

    return True


def rdmb_blend_AC(mb_obj_file,
                  out_file_base,
                  pas, pcs, ccrs,
                  max_ite=2000):
    """Reaction-diffusion on metaball (blended A, C)

    Reaction-diffusion on metaball (blended A, C)
    params A, C for each vertex are calculated based on
    the contribution ratio of each metaball to the vertex

    Args:
        mb_obj_file: metaball OBJ file
        out_file_base: output file basename
        pas, pcs: np array of params A and C for each metaball
        ccrs: coefficient for contribution ratio (1/0)
        max_ite: maximum number of iterations

    """

    mbs, vs, vc, vn, fs = load_mb_obj(mb_obj_file)

    # crs: corrected contribution ratios
    crs = mbs.cr_v(vs, 1.5) * ccrs
    # normalization
    crs = crs / np.c_[np.linalg.norm(crs, ord=1, axis=1)]

    # make A gradient
    A = np.ravel(np.dot(crs, np.c_[pas]))

    B = 0.08

    # make C gradient
    C = np.ravel(np.dot(crs, np.c_[pcs]))

    D = 0.03
    E = 0.10
    F = 0.12
    G = 0.06

    prms = (A, B, C, D, E, F, G)

    rdmb(vs, mbs, out_file_base, prms, max_ite)

    return True


def rdmb_povray_save(out_file,
                     vs,
                     ucs, vcs,
                     width=800, height=600,
                     rotx=0, roty=0, rotz=0,
                     angle=14):
    """Render and save RD results using Pov-Ray

    Render and save RD results using Pov-Ray

    Args:
        out_file: output file
        vs: vertices
        ucs, vcs: u/v conc.
        width, height: width and height of output image
        rotx, roty, rotz: rotation angle
        angle: camera angle

    """

    ucmax = 6.0
    ucs = ucs / ucmax
    ucs[ucs > 1.0] = 1.0
    # ucs = ucs / np.max(ucs)

    rot1 = [rotx, 0, 0]
    rot2 = [0, roty, 0]
    rot3 = [0, 0, rotz]

    camera = Camera('location', [0, 0, -25],
                    'look_at', [0, 0, 0],
                    'angle', angle,
                    'right x*image_width/image_height')

    light = LightSource([-3, 2, -6], 'color', [1.0, 1.0, 1.0], 'parallel')
    light2 = LightSource([2, -2, -6], 'color', [0.6, 0.6, 0.6], 'parallel')
    background = Background('color', [1, 1, 1, 1])

    spheres = [Sphere(v, 0.02,
                      Finish('ambient', 0.2, 'diffuse', 0.8, 'phong', 1.0),
                      Texture(Pigment('color',
                                      [0.3+uc*0.7, 0.2+uc*0.8, 0.2+uc*0.8])),
                      'rotate', rot1,
                      'rotate', rot2,
                      'rotate', rot3) for v, uc in zip(vs, ucs)]

    objects = [light, light2, background] + spheres

    scene = Scene(camera, objects=objects)
    scene.render(out_file, width=width, height=height,
                 output_alpha=True, antialiasing=0.001,
                 tempfile=out_file+"__temp__.pov")


def rdmb_povray(file_base,
                time_point=2000,
                width=800, height=600,
                angle=14):
    """Load RD results and Render/Save image using Pov-Ray

    Load RD results and Render/Save image using Pov-Ray

    Args:
        file_base:
        time_point:
        width, height:
        angle:

    Returns:
        file_png

    """

    file_prms = file_base + "_prms.npz"
    vs, *_ = load_rd_prms(file_prms)
    file_uv = file_base + "_{:05}.npz".format(time_point)

    file_png = file_base + "_{:05}.png".format(time_point)
    ucs, vcs = load_rd_uv(file_uv)

    rdmb_povray_save(file_png,
                     vs,
                     ucs, vcs,
                     width=width, height=height,
                     rotx=0, roty=0, rotz=0,
                     angle=angle)

    return file_png


def rdmb_povray_save_q(out_file,
                       vs,
                       ucs, vcs,
                       width=800, height=600,
                       rotx=0, roty=0, rotz=0,
                       angle=14):
    """Render and save RD results using Pov-Ray (for quantification)

    Render and save RD results using Pov-Ray (for quantification)

    Args:
        out_file: output file
        vs: vertices
        ucs: u conc.
        vcs: v conc.
        rotx, roty, rotz: rotation angle
        width, height: width and height of output image
        angle: camera angle

    """

    ucmax = 6.0
    ucs = ucs / ucmax
    ucs[ucs > 1.0] = 1.0
    # ucs = ucs / np.max(ucs)

    rot1 = [rotx, 0, 0]
    rot2 = [0, roty, 0]
    rot3 = [0, 0, rotz]

    camera = Camera('location', [0, 0, -25],
                    'look_at', [0, 0, 0],
                    'angle', angle,
                    'right x*image_width/image_height')

    light = LightSource([0, 0, -10],
                        'color', [1.0, 1.0, 1.0], 'parallel', 'shadowless')
    light1 = LightSource([-10, 0, 0],
                         'color', [0.5, 0.5, 0.5], 'parallel', 'shadowless')
    light2 = LightSource([10, 0, 0],
                         'color', [0.5, 0.5, 0.5], 'parallel', 'shadowless')
    light3 = LightSource([0, -10, 0],
                         'color', [0.5, 0.5, 0.5], 'parallel', 'shadowless')
    light4 = LightSource([0, 10, 0],
                         'color', [0.5, 0.5, 0.5], 'parallel', 'shadowless')

    background = Background('color', [1, 1, 1, 1])

    spheres = [Sphere(v, 0.02,
                      Finish('ambient', 1.0),
                      Texture(Pigment('color',
                                      [0.3+uc*0.7, 0.2+uc*0.8, 0.2+uc*0.8])),
                      'rotate', rot1,
                      'rotate', rot2,
                      'rotate', rot3) for v, uc in zip(vs, ucs)]

    objects = [light, light1, light2, light3, light4, background] + spheres

    scene = Scene(camera, objects=objects)
    scene.render(out_file, width=width, height=height,
                 output_alpha=True, antialiasing=0.001,
                 tempfile=out_file+"__temp__.pov")


def rdmb_povray_q(file_base,
                  time_point=2000,
                  width=800, height=600,
                  angle=14):
    """Load RD results and Render/Save image using Pov-Ray (for quantification)

    Load RD results and Render/Save image using Pov-Ray (for quantification)

    Args:
        file_base:
        time_point:
        width, height:
        angle:

    Returns:
        file_png

    """

    file_prms = file_base + "_prms.npz"
    vs, *_ = load_rd_prms(file_prms)
    file_uv = file_base + "_{:05}.npz".format(time_point)

    file_png = file_base + "_{:05}.png".format(time_point)
    ucs, vcs = load_rd_uv(file_uv)

    rdmb_povray_save_q(file_png,
                       vs,
                       ucs, vcs,
                       width=width, height=height,
                       rotx=0, roty=0, rotz=0,
                       angle=angle)

    return file_png


def rdmb_povray_color(file_base,
                      time_point=2000,
                      width=800, height=600,
                      rotx=0, roty=0, rotz=0,
                      angle=14,
                      mode="C"):
    """Render and save RD results using Pov-Ray (color)

    Render and save RD results using Pov-Ray
    with color indicating parameter values

    Args:
        file_base:
        time_point:
        width, height:
        rotx, roty, rotz:
        angle:
        mode:
        
    """

    vs, ucs, As, Cs = load_rd_mb(file_base)
    
    file_png = file_base + "_color_{:05}.png".format(time_point)
    
    tempfile = file_png[:-4] + "__temp__" + ".pov"

    camera = Camera('location', [0, 0, -25],
                    'look_at', [0, 0, 0],
                    'angle', angle,
                    'right x*image_width/image_height')
    
    light = LightSource([-3, 2, -6],
                        'color', [1.0, 1.0, 1.0], 'parallel')
    light2 = LightSource([2, -2, -6],
                         'color', [0.2, 0.2, 0.2], 'parallel')
    background = Background('color', [1, 1, 1, 1])
    
    spheres = []
    spheres +=  sph(vs, ucs, As, Cs,
                    0, 0, 0,
                    rotx=rotx, roty=roty, rotz=rotz,
                    mode=mode)
    
    objects = [light, light2, background] + spheres
    
    scene = Scene(camera, objects=objects)
    
    scene.render(file_png,
                 width=width, height=height,
                 tempfile=tempfile,
                 output_alpha=True, antialiasing=0.001)

    return file_png


def sph(vs, ucs, A, C, x, y, z, rotx=0, roty=0, rotz=0, mode="C"):
    """Colored spheres for rendering RD results

    Colored spheres for rendering RD results

    Args:
        vs: vertices
        ucs: u conc.
        A, C: RD parameters
        x, y, z: translation
        rotx, roty, rotz: rotation angle
        mode: color mode for parameters

    """
    sph = []

    if mode == "A":
        nA = (A-0.07)/(0.12-0.07)
        nC = (C+0.1)/(0.31+0.1)
    elif mode == "C":
        nA = (A-0.07)/(0.12-0.07)
        nC = C/0.307
    elif mode == "AC":
        nA = (A-0.07)/(0.12-0.07)
        nC = (C+0.1)/(0.31+0.1)
    else:
        nA = (A-0.07)/(0.12-0.07)
        nC = C/0.307

    if (type(nA) is np.float64):
        nA = np.full(len(vs), nA)
    if (type(nC) is np.float64):
        nC = np.full(len(vs), nC)

    for v, uc, a, c in zip(vs, ucs, nA, nC):
        if mode == "A":
            H0 = a
            L0 = 1/(1+np.exp((-2.8*((a-0.52)*2.4))))
            # R0, G0, B0 = colorsys.hls_to_rgb(0.02+H0*0.35, 0.5-L0*0.4, 1.0)
            R0, G0, B0 = colorsys.hls_to_rgb(1.0-H0*0.40, 0.5-L0*0.4, 1.0)

        elif mode == "C":
            H0 = c
            L0 = 1/(1+np.exp((-2.8*((c-0.52)*2.4))))
            R0, G0, B0 = colorsys.hls_to_rgb(0.02+H0*0.35, 0.5-L0*0.4, 1.0)

        elif mode == "AC":
            R0 = a*1.0
            # G0 = max(0.8-(max(a+c, 0)), 0)
            G0 = 0.0
            B0 = c*1.0

        else:
            R0 = 0.3
            G0 = 0.2
            B0 = 0.2

        R1 = 1.0 - R0
        G1 = 1.0 - G0
        B1 = 1.0 - B0

        sph.append(Sphere(v, 0.022,
                          Texture(Pigment('color',
                                          [R0+uc*R1, G0+uc*G1, B0+uc*B1]),
                                  Finish('phong', 0.7,
                                         'specular', 0.2,
                                         'diffuse', 0.9,
                                         'ambient', 0.1)),
                          'rotate', [rotx, 0, 0],
                          'rotate', [0, roty, 0],
                          'rotate', [0, 0, rotz],
                          'translate', [x, y, z],
                          'no_shadow'))

    return sph
