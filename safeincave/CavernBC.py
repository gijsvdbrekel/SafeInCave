# Copyright 2025 The safeincave community.
#
# This file is part of safeincave.
#
# Licensed under the GNU GENERAL PUBLIC LICENSE, Version 3 (the "License"); you may not
# use this file except in compliance with the License.  You may obtain a copy
# of the License at
#
#     https://spdx.org/licenses/GPL-3.0-or-later.html
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  See the
# License for the specific language governing permissions and limitations under
# the License.
from abc import ABC, abstractmethod
import CoolProp.CoolProp as CP
from safeincave import CavernThermodynamics
import numpy as np


class CavernVolumeComputer:
    def __init__(self, grid, boundary_name: str, internal_point: list[float,float,float] = None, sym_scale: int = 1):
        self.grid = grid
        self.boundary_name = boundary_name
        self.sym_scale = self.validate_sym_scale(sym_scale)
        self.internal_point = internal_point
        cavern_data = self.__extract_cavern_surface_from_grid(self.boundary_name)
        self.coords_wall, self.conn_wall, self.ids_wall = cavern_data
        if internal_point is None:
            self.internal_point = self.__surface_centroid(self.coords_wall, self.conn_wall)
        self.gather_cells_for_wall_vertices()

    def validate_sym_scale(self, sym_scale):
        '''
        Validate that sym_scale is a positive integer and that the mesh represents a fraction 1/sym_scale of the full cavern.
        It can be either 1 (full cavern), 2 (half cavern), or 4 (quarter cavern).
        '''
        if sym_scale not in [1, 2, 4]:
            raise ValueError(f"sym_scale must be 1, 2, or 4. Got {sym_scale}.")
        return sym_scale


    def calculate_normals(self):
        normals = np.zeros((len(self.conn_wall), 3))
        for i, tri in enumerate(self.conn_wall):
            p0, p1, p2 = self.coords_wall[tri]
            normal = np.cross(p1 - p0, p2 - p0)
            normals[i] = normal / np.linalg.norm(normal)
        return normals
    

    def gather_cells_for_wall_vertices(self):
        tdim = self.grid.mesh.topology.dim
        self.grid.mesh.topology.create_connectivity(0, tdim)
        v2c = self.grid.mesh.topology.connectivity(0, tdim)
        self.wall_cells = np.empty(len(self.ids_wall), dtype=np.int32)
        for k, v in enumerate(self.ids_wall):
            links = v2c.links(v)
            if len(links) == 0:
                raise RuntimeError(f"Vertex {v} has no incident cell.")
            self.wall_cells[k] = links[0]


    def compute(self, u=None):
        if u is None:
            coords = self.coords_wall
        else:
            disp_wall = u.eval(self.coords_wall, self.wall_cells)
            coords = self.coords_wall + disp_wall
        volume = 0.0
        for tri in self.conn_wall:
            v0 = coords[tri][0] - self.internal_point
            v1 = coords[tri][1] - self.internal_point
            v2 = coords[tri][2] - self.internal_point
            volume += np.dot(v0, np.cross(v1, v2)) / 6.0
        return abs(volume) * self.sym_scale

    def __extract_cavern_surface_from_grid(self, boundary_name: str):
        """
        Return (coords_wall, tris_local, wall_ids) for the named boundary.

        The returned triangles are oriented so that their normals point outward
        from the 3D domain (i.e. away from the adjacent tetrahedral cell).

        Parameters
        ----------
        grid : GridHandlerGMSH-like
            Object exposing:
            - grid.mesh
            - grid.get_boundary_tag(boundary_name)
            - grid.boundaries or grid.facet_tags
        boundary_name : str
            Name of the boundary to extract.

        Returns
        -------
        coords_wall : (n_wall, 3) ndarray
            Coordinates of wall vertices, in the same order as wall_ids.
        tris_local : (n_tris, 3) ndarray
            Triangle connectivity indexing into coords_wall, with consistent
            outward orientation.
        wall_ids : (n_wall,) ndarray
            Global vertex ids of the wall vertices.
        """
        mesh = self.grid.mesh
        tdim = mesh.topology.dim
        fdim = tdim - 1

        if tdim != 3:
            raise RuntimeError("This function currently assumes a 3D tetrahedral mesh.")

        tag = self.grid.get_boundary_tag(boundary_name)

        mt = getattr(self.grid, "boundaries", None) or getattr(self.grid, "facet_tags", None)
        if mt is None:
            raise RuntimeError("Grid does not expose facet tags as 'boundaries' or 'facet_tags'.")

        # Facets carrying the requested boundary tag
        facets = mt.indices[mt.values == tag]
        if len(facets) == 0:
            raise RuntimeError(f"No facets found for boundary '{boundary_name}' (tag={tag}).")

        # Required connectivities
        mesh.topology.create_connectivity(fdim, 0)     # facet -> vertices
        mesh.topology.create_connectivity(fdim, tdim)  # facet -> cell(s)
        mesh.topology.create_connectivity(tdim, 0)     # cell  -> vertices

        f2v = mesh.topology.connectivity(fdim, 0)
        f2c = mesh.topology.connectivity(fdim, tdim)
        c2v = mesh.topology.connectivity(tdim, 0)

        coords_all = mesh.geometry.x

        tri_global = []
        wall_set = set()

        for f in facets:
            facet_vertices = np.array(f2v.links(f), dtype=np.int64)
            if len(facet_vertices) != 3:
                # Expected for tetrahedral boundary facets
                continue

            attached_cells = f2c.links(f)
            if len(attached_cells) != 1:
                raise RuntimeError(
                    f"Boundary facet {f} should have exactly one adjacent cell, got {len(attached_cells)}."
                )

            cell = attached_cells[0]
            cell_vertices = np.array(c2v.links(cell), dtype=np.int64)

            # The tetra vertex opposite to this facet
            opposite = np.setdiff1d(cell_vertices, facet_vertices)
            if len(opposite) != 1:
                raise RuntimeError(
                    f"Could not identify unique opposite vertex for facet {f} in cell {cell}."
                )
            opposite_vertex = opposite[0]

            v0, v1, v2 = facet_vertices
            p0 = coords_all[v0]
            p1 = coords_all[v1]
            p2 = coords_all[v2]
            p_opp = coords_all[opposite_vertex]

            # Normal from current ordering
            n = np.cross(p1 - p0, p2 - p0)

            # Vector from facet toward the cell interior (toward the opposite tetra vertex)
            to_interior = p_opp - p0

            # If normal points toward the interior, flip triangle orientation
            if np.dot(n, to_interior) > 0.0:
                facet_vertices = np.array([v0, v2, v1], dtype=np.int64)

            tri_global.append(facet_vertices)
            wall_set.update(facet_vertices.tolist())

        if not tri_global:
            raise RuntimeError(f"No triangular facets found for boundary '{boundary_name}' (tag={tag}).")

        wall_ids = np.array(sorted(wall_set), dtype=np.int64)
        gid2lid = {gid: i for i, gid in enumerate(wall_ids)}

        tris_local = np.array(
            [[gid2lid[v] for v in tri] for tri in tri_global],
            dtype=np.int32
        )
        coords_wall = coords_all[wall_ids]
        return coords_wall, tris_local, wall_ids

    
    def __surface_centroid(self, coordinates, triangles): 
        """Calculate the centroid of a surface defined by triangles."""
        total_area = 0
        weighted_sum = np.zeros(3)
        for tri in triangles:
            p0, p1, p2 = coordinates[tri]
            center = (p0 + p1 + p2) / 3.0
            area = np.linalg.norm(np.cross(p1 - p0, p2 - p0)) / 2.0
            weighted_sum += center * area
            total_area += area
        return weighted_sum / total_area







class Cavern(ABC):
    def __init__(self,
                 cavern_name: str,
                 fluid: str,
                 h_conv: float = None):
        self.cavern_name = cavern_name
        self.check_fluid(fluid)
        self.fluid = fluid
        self.h_conv = h_conv
        self.type = None

    def check_fluid(self, fluid: str) -> None:
        try:
            CP.AbstractState("HEOS", fluid)
        except ValueError:
            raise ValueError(f"Fluid '{fluid}' not recognized by CoolProp.")
        
    @abstractmethod
    def update_cavern(self, t: float) -> None:
        pass

    @abstractmethod
    def record_data(self, t: float) -> None:
        pass


class CavernHandler:
    def __init__(self):
        self.caverns_T = []
        self.caverns_PT = []
        self.caverns_MFlux = []

    def add_cavern(self, cavern: Cavern) -> None:
        if cavern.type == "Cavern_T":
            self.caverns_T.append(cavern)
        elif cavern.type == "Cavern_PT":
            self.caverns_PT.append(cavern)
        elif cavern.type == "Cavern_MFlux":
            self.caverns_MFlux.append(cavern)
        else:
            raise ValueError(f"Cavern type {cavern.type} not supported")
        
    def calculate_volumes(self, u: any) -> None:
        for cavern in self.caverns_MFlux:
            cavern.calculate_volume(u)

    def update_caverns(self, 
                       t: float,
                       dt: float = None,
                       u: any = None,
                       T: any = None) -> None:
        for cavern in self.caverns_T:
            cavern.update_cavern(t)

        for cavern in self.caverns_PT:
            cavern.update_cavern(t)

    def record_cavern_data(self, t: float) -> None:
        for cavern in self.caverns_T:
            cavern.record_data(t)

        for cavern in self.caverns_PT:
            cavern.record_data(t)

        for cavern in self.caverns_MFlux:
            cavern.record_data(t)

        


class Cavern_T(Cavern):
    def __init__(self, cavern_name: str, T_values: list, time_values: list, h_conv: float = None):
        super().__init__(cavern_name, None, h_conv)
        self.type = "Cavern_T"
        self.T_values = T_values
        self.time_values = time_values
        self.T = self.T_values[0]

        # Initialize histories
        self.T_hist = [self.T]
        self.t_hist = [self.time_values[0]]

    def update_cavern(self, t: float) -> None:
        self.T = np.interp(t, self.time_values, self.T_values)
        if self.T <= 0.0:
            raise ValueError(f"T must be > 0, got {self.T}")
    
    def record_data(self, t: float) -> None:
        self.T_hist.append(self.T)
        self.t_hist.append(t)



class Cavern_PT(Cavern):
    def __init__(self, 
                 cavern_name: str,
                 fluid: str,
                 P_values: list,        # Gauge pressure values (Pa)
                 T_values: list,        # Temperature values (K)
                 time_values: list,
                 ref_pos: float = 0.0,
                 direction: int = 0,
                 g: float = -9.81,
                 h_conv: float = None,
                 P_atm: float = 101325.0 # Atmospheric pressure in Pa
                ):
        super().__init__(cavern_name, fluid, h_conv)
        self.type = "Cavern_PT"
        self.P_values = P_values
        self.T_values = T_values
        self.time_values = time_values
        self.P_atm = P_atm
        self.ref_pos = ref_pos
        self.direction = direction
        self.gravity = g
        self.__initialize()
        
        # Initialize histories
        self.T_hist = [self.T]
        self.P_hist = [self.P]
        self.t_hist = [self.time_values[0]]
        self.density_hist = [self.density]

    def __initialize(self) -> None:
        self.AS = CP.AbstractState("HEOS", self.fluid)
        self.P = self.P_values[0]
        self.T = self.T_values[0]
        self.AS.update(CP.PT_INPUTS, self.P, self.T)
        self.density = self.AS.rhomass()

    def update_cavern(self, t: float) -> None:
        self.P = np.interp(t, self.time_values, self.P_values)
        self.T = np.interp(t, self.time_values, self.T_values)
        if self.P <= 0.0:
            raise ValueError(f"P must be > 0, got {self.P}")
        if self.T <= 0.0:
            raise ValueError(f"T must be > 0, got {self.T}")
        self.AS.update(CP.PT_INPUTS, self.P + self.P_atm, self.T)
        self.density = self.AS.rhomass()

    def record_data(self, t: float) -> None:
        self.P_hist.append(self.P)
        self.T_hist.append(self.T)
        self.t_hist.append(t)
        self.density_hist.append(self.density)



class Cavern_MassFlux(Cavern):
    def __init__(self,
                 *,
                 grid: any,
                 cavern_name: str,
                 fluid: str,
                 P_init: float,         # Initial gauge pressure (Pa)
                 T_init: float,         # Initial fluid temperature (K)
                 T_in: float,           # Temperature of injected fluid (K)
                 Mflux_values: list,
                 time_values: list,
                 ref_pos: float = 0.0,
                 direction: int = 0,
                 g: float = -9.81,
                 h_conv: float = None,
                 P_atm: float = 101325.0 # Atmospheric pressure in Pa
                 ):
        super().__init__(cavern_name, fluid, h_conv)
        self.type = "Cavern_MFlux"
        self.Mflux_values = Mflux_values
        self.time_values = time_values
        self.Mflux = self.Mflux_values[0]
        self.ref_pos = ref_pos
        self.direction = direction
        self.gravity = g
        self.P_atm = P_atm
        self.P_init = P_init
        self.T_init = T_init
        self.T_in = T_in

        # Initial absolut pressure and temperature
        self.P0 = P_init + self.P_atm
        self.T0 = T_init

        # Initialize cavern volume computer
        self.cvc = CavernVolumeComputer(grid, boundary_name=self.cavern_name)

        # Initialize thermodynamic model
        self.model = CavernThermodynamics(self.fluid)

        # Calculate initial density and mass
        AS = CP.AbstractState("HEOS", self.fluid)
        AS.update(CP.PT_INPUTS, self.P0, self.T0)
        self.density = AS.rhomass()

        # Compute initial mass
        self.M0 = self.density * self.V0

        # Initialize histories
        self.V_hist = []
        self.M_hist = []
        self.P_hist = []
        self.T_hist = []
        self.density_hist = []
        self.t_hist = []

    def calculate_volume(self, u: any) -> None:
        self.V1 = self.cvc.compute(u)

    def update_cavern(self, t: float, dt: float, u: any = None) -> None:
        self.Mflux = np.interp(t, self.time_values, self.Mflux_values)
        self.M1 = self.M0 + self.Mflux * dt
        dm = self.M1 - self.M0
        self.P, self.T, self.density = self.model.solve(
            dm = dm,
            Q_in = 0.0,
            T_in = self.T_in,
            P0 = self.P0,
            T0 = self.T0,
            V0 = self.V0,
            V1 = self.V1
        )

    def record_data(self, t: float) -> None:
        self.V_hist.append(self.V1)
        self.M_hist.append(self.M1)
        self.P_hist.append(self.P)
        self.T_hist.append(self.T)
        self.t_hist.append(t)
        self.density_hist.append(self.density)

        # Update for next step
        self.P0 = self.P
        self.T0 = self.T
        self.V0 = self.V1



    