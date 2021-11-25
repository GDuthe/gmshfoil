from airfoils import Airfoil
from scipy import interpolate
import numpy as np
import gmsh
import ntpath

class GMSHFoil:
    def __init__(self, foil_dat_file, naca_foil = '4812', mesh_name = None):
        """
        A thin wrapper class that creates a 2D airfoil mesh and sets
        names for the surfaces of the boundaries.
        """
        self.mesh_name = mesh_name
        if foil_dat_file is not None:
            self.airfoil_coords = np.loadtxt(foil_dat_file, skiprows=1)

            if self.mesh_name is None:
                self.mesh_name = ntpath.basename(foil_dat_file)[:-4]
        else:
            foil = Airfoil.NACA4(naca_foil)
            x_points = np.concatenate([np.linspace(0.0, 0.01, 100),np.linspace(0.01, 1.0, 200)[1:]])
            x_points_reversed = list(reversed(x_points))[:-1]
            self.airfoil_coords = np.array([[*x_points_reversed, *x_points],
                                            [*foil.y_upper(x_points_reversed), *foil.y_lower(x_points)]]).T

            if self.mesh_name is None:
                self.mesh_name = 'naca_' + naca_foil
        self.gmsh = gmsh
        
    def create_2d_unstructured_foil_mesh(self,
                                         use_base_coords = False,
                                         npoints_airfoil = 1500,
                                         h_0 = 0.2,
                                         R_b = 100,
                                         h_extrude = 1.,
                                         refine_wake_len = 0.25,
                                         h_w = 0.01):
        """
        Create the 2D unstructured mesh for simulating the airfoil.
        
        Arguments:
          use_base_coords  : (False) if False a spline fitting should be used to obtain coordinates
          npoints_airfoil   : (1500) number of points on the airfoil surface, can only be used if not using base coords
          h_0               : (0.2)  h-refinement size for overall sizing
          R_b               : the radius of the outer boundary.
          h_extrude         : (1.) openfoam needs an extrusion height. The geometry plane is extruded 
                              to a unit height in order to create cells.
          refine_wake_len   : (0.25) the length of the wake refinement line
          h_w               : (0.01) sizing of the refinement wake line
          
        """
        
        _gmsh = self.gmsh
        _gmsh.initialize()
        _gmsh.model.add(self.mesh_name)

        # points for the airfoil:
        if use_base_coords:
            coo = self.airfoil_coords
        else:
            t = np.linspace(0.0, 1.0, npoints_airfoil)
            foil_spline = interpolate.splprep([self.airfoil_coords[:, 0], self.airfoil_coords[:, 1]], s=0.0000002, k=3)
            coo = np.array(interpolate.splev(t, foil_spline[0], der=0)).transpose()

        gmsh_airf_points = [_gmsh.model.occ.addPoint(coo[k,0], coo[k,1], 0, 0.0, k) for k in range(len(coo))]
        n_airf_points = len(gmsh_airf_points)

        # create lines for airfoil:
        foil_curve = []
        start_end_lines_foil = []

        for k in range(n_airf_points):
            k_start, k_end = k, k + 1
            if k_end > max(gmsh_airf_points):
                k_end = min(gmsh_airf_points)

            _gmsh.model.occ.addLine(k_start, k_end, k)
            start_end_lines_foil.append((k_start, k_end))
            foil_curve.append(k)

        # create circular boundary for the boundary:
        boundary = _gmsh.model.occ.addCircle(0.5, 0.0, 0.0, R_b)

        # creating a circular surface with a foil-shaped hole:
        _gmsh.model.occ.addCurveLoop([boundary], 1)
        _gmsh.model.occ.addCurveLoop(foil_curve, 2)
        _gmsh.model.occ.addPlaneSurface([1,2], 0)

        # add refinement line in the airfoil wake
        p1, p2 = _gmsh.model.occ.addPoint(1.05, 0, 0, h_w), gmsh.model.occ.addPoint(1.05 + refine_wake_len, 0, 0, h_w)
        refine_line = _gmsh.model.occ.addLine(p1, p2)
        _gmsh.model.occ.synchronize()
        _gmsh.model.mesh.embed(1, [refine_line], 2, 0)

        # extrude the surface and create farfield surface
        e = _gmsh.model.occ.extrude([(2,0)], 0, 0, h_extrude, numElements=[1], recombine=True)
        farfield_surface =_gmsh.model.occ.addSurfaceLoop([0,e[0][1]])
        _gmsh.model.occ.synchronize()

        # definition of physical groups:
        front_surface = _gmsh.model.addPhysicalGroup(2,[e[0][1]])
        _gmsh.model.setPhysicalName(2, front_surface, 'FRONT')
        back_surface = _gmsh.model.addPhysicalGroup(2, [0])
        _gmsh.model.setPhysicalName(2, back_surface, 'BACK')
        airfoil_surface = _gmsh.model.addPhysicalGroup(2, list(range(2, n_airf_points + 2)))
        _gmsh.model.setPhysicalName(2, airfoil_surface, 'AIRFOIL')
        farfield_surface = _gmsh.model.addPhysicalGroup(2, [farfield_surface])
        _gmsh.model.setPhysicalName(2, farfield_surface, 'FARFIELD')
        internal_volume = _gmsh.model.addPhysicalGroup(3, [e[1][1]])
        _gmsh.model.setPhysicalName(3, internal_volume, "internalField")

        _gmsh.option.setNumber("Mesh.Algorithm", 6)
        _gmsh.option.setNumber("Mesh.CharacteristicLengthFactor", h_0)
        _gmsh.option.setNumber('Mesh.MshFileVersion' , 2.10)
        _gmsh.model.mesh.generate(3)
        _gmsh.model.mesh.optimize("Netgen")
        _gmsh.write("%s.msh"%self.mesh_name)
        
    def view(self):
        self.gmsh.fltk.run()

if __name__ == '__main__':
    gmsh_foil = GMSHFoil(foil_dat_file=None, naca_foil = '0012', mesh_name = None)

    # creates the mesh
    gmsh_foil.create_2d_unstructured_foil_mesh()

    
