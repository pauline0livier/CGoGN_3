/*******************************************************************************
 * CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
 * Copyright (C), IGG Group, ICube, University of Strasbourg, France            *
 *                                                                              *
 * This library is free software; you can redistribute it and/or modify it      *
 * under the terms of the GNU Lesser General Public License as published by the *
 * Free Software Foundation; either version 2.1 of the License, or (at your     *
 * option) any later version.                                                   *
 *                                                                              *
 * This library is distributed in the hope that it will be useful, but WITHOUT  *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
 * for more details.                                                            *
 *                                                                              *
 * You should have received a copy of the GNU Lesser General Public License     *
 * along with this library; if not, write to the Free Software Foundation,      *
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
 *                                                                              *
 * Web site: http://cgogn.unistra.fr/                                           *
 * Contact information: cgogn@unistra.fr                                        *
 *                                                                              *
 *******************************************************************************/

#include <cgogn/modeling/algos/incidenceGraph_to_hex.h>

#include <cgogn/core/types/cell_marker.h>

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/traversals/face.h>
#include <cgogn/core/functions/traversals/global.h>
#include <cgogn/core/functions/traversals/halfedge.h>
#include <cgogn/core/functions/traversals/vertex.h>

#include <cgogn/core/functions/mesh_ops/edge.h>
#include <cgogn/core/functions/mesh_ops/face.h>
#include <cgogn/core/functions/mesh_ops/vertex.h>
#include <cgogn/core/functions/mesh_ops/volume.h>

#include <cgogn/core/functions/mesh_info.h>

#include <cgogn/core/types/cmap/cmap_ops.h>
#include <cgogn/core/types/cmap/dart_marker.h>
#include <cgogn/core/types/mesh_views/cell_cache.h>

#include <cgogn/io/surface/surface_import.h>

#include <cgogn/geometry/algos/centroid.h>
#include <cgogn/geometry/algos/length.h>
#include <cgogn/geometry/algos/normal.h>
#include <cgogn/geometry/algos/picking.h>
#include <cgogn/modeling/algos/subdivision.h>

#include <cgogn/geometry/functions/angle.h>
#include <cgogn/geometry/functions/fitting.h>
#include <cgogn/geometry/functions/intersection.h>
#include <cgogn/geometry/functions/projection.h>

#include <fstream>
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <numeric>

namespace cgogn
{

namespace modeling
{

std::tuple<IG_GAttributes, IG_M2Attributes, IG_M3Attributes> incidenceGraph_to_hex(IncidenceGraph& ig, CMap2& m2,
																				   CMap3& m3)
{
	IncidenceGraphData igData;
	IG_GAttributes igAttribs;
	IG_M2Attributes m2Attribs;
	IG_M3Attributes m3Attribs;
	std::cout << " start" << std::endl;

	bool okay = get_incidenceGraph_data(ig, igData);
	std::cout << uint32(igData.intersections.size()) << " intersections" << std::endl;
	std::cout << uint32(igData.branches.size()) << " branches" << std::endl;
	std::cout << uint32(igData.efjunctures.size()) << " efjunctures" << std::endl;
	std::cout << uint32(igData.ffjunctures.size()) << " ffjunctures" << std::endl;
	std::cout << uint32(igData.leaflets.size()) << " leaflets" << std::endl;

	if (!okay)
		std::cout << "error incidenceGraph_to_hex: get_incidenceGraph_data" << std::endl;
	else
		std::cout << "incidenceGraph_to_hex (/): got incidenceGraph data" << std::endl;

	if (okay)
		okay = add_incidenceGraph_attributes(ig, igAttribs);
	if (!okay)
		std::cout << "error incidenceGraph_to_hex: add_incidenceGraph_attributes" << std::endl;
	else
		std::cout << "incidenceGraph_to_hex (/): added incidenceGraph attributes" << std::endl;

	if (okay)
		okay = compute_faces_geometry(ig, igData, igAttribs);
	if (!okay)
		std::cout << "error incidenceGraph_to_hex: compute_faces_geometry" << std::endl;
	else
		std::cout << "incidenceGraph_to_hex (/): computed incidence graph face geometries" << std::endl;

	if (okay)
		okay = add_cmap2_attributes(m2, m2Attribs);
	if (!okay)
		std::cout << "error incidenceGraph_to_hex: add_cmap2_attributes" << std::endl;
	else
		std::cout << "incidenceGraph_to_hex (/): added cmap2 attributes" << std::endl;

	if (okay)
		okay = build_contact_surfaces(ig, igAttribs, igData, m2, m2Attribs);
	if (!okay)
		std::cout << "error incidenceGraph_to_hex: build_contact_surfaces" << std::endl;
	else
		std::cout << "incidenceGraph_to_hex (/): contact surfaces built" << std::endl;

	if (okay)
		okay = create_intersection_frames(ig, igAttribs, igData, m2, m2Attribs);
	if (!okay)
		std::cout << "error incidenceGraph_to_hex: create_intersection_frames" << std::endl;
	else
		std::cout << "incidenceGraph_to_hex (/): create_intersection_frames completed" << std::endl;

	if (okay)
		okay = propagate_frames(ig, igAttribs, igData, m2);
	if (!okay)
		std::cout << "error incidenceGraph_to_hex: propagate_frames" << std::endl;
	else
		std::cout << "incidenceGraph_to_hex (/): propagate_frames completed" << std::endl;

	if (okay)
		okay = set_contact_surfaces_geometry(ig, igAttribs, m2, m2Attribs);
	if (!okay)
		std::cout << "error incidenceGraph_to_hex: set_contact_surfaces_geometry" << std::endl;
	else
		std::cout << "incidenceGraph_to_hex (/): set_contact_surfaces_geometry completed" << std::endl;

	if (okay)
		okay = build_volumes(ig, igAttribs, igData, m2, m2Attribs, m3);
	if (!okay)
		std::cout << "error incidenceGraph_to_hex: build_volumes" << std::endl;
	else
		std::cout << "incidenceGraph_to_hex (/): build_volumes completed" << std::endl;

	if (okay)
		okay = sew_sections_igh(m2, m2Attribs, m3);
	if (!okay)
		std::cout << "error incidenceGraph_to_hex: sew_sections" << std::endl;
	else
		std::cout << "incidenceGraph_to_hex (/): sew_sections completed" << std::endl;

	if (okay)
	{
		add_cmap3_attributes_igh(m3, m3Attribs);
		okay = set_volumes_geometry_igh(ig, igAttribs, igData, m2, m2Attribs, m3, m3Attribs);
	}
	if (!okay)
		std::cout << "error graph_to_hex: set_volumes_geometry" << std::endl;
	else
		std::cout << "graph_to_hex (/): set_volumes_geometry completed" << std::endl;

	std::cout << "CC: " << nb_cells<CMap3::CC>(m3) << std::endl;
	std::cout << "surface " << nb_cells<CMap2::Vertex>(m2) << ":" << nb_cells<CMap2::Edge>(m2) << ":"
			  << nb_cells<CMap2::Face>(m2) << std::endl;

	return {igAttribs, m2Attribs, m3Attribs};
}

/*****************************************************************************/
/* data preparation                                                          */
/*****************************************************************************/

bool add_incidenceGraph_attributes(IncidenceGraph& ig, IG_GAttributes& igAttribs)
{
	igAttribs.vertex_position = get_attribute<Vec3, IncidenceGraph::Vertex>(ig, "position");
	if (!igAttribs.vertex_position)
	{
		std::cout << "The incidence graph has no vertex position attribute" << std::endl;
		return false;
	}
	igAttribs.vertex_radius = get_attribute<Scalar, IncidenceGraph::Vertex>(ig, "radius");
	if (!igAttribs.vertex_radius)
	{
		std::cout << "The incidence graph has no vertex radius attribute" << std::endl;
		return false;
	}

	igAttribs.face_normal = add_attribute<Vec3, IncidenceGraph::Face>(ig, "face_normal");

	igAttribs.vertex_contact_surface = add_attribute<Dart, IncidenceGraph::Vertex>(ig, "vertex_contact_surface");
	igAttribs.halfedge_contact_surface_face =
		add_attribute<std::pair<Dart, Dart>, IncidenceGraph::Edge>(ig, "halfedge_contact_surface_face");
	igAttribs.halfedge_frame = add_attribute<std::pair<Mat3, Mat3>, IncidenceGraph::Edge>(ig, "halfedge_frame");

	igAttribs.face_edge_dart = add_attribute<std::vector<Dart>, IncidenceGraph::Face>(ig, "face_edge_dart");
	igAttribs.vertex_boundary_edge_dart =
		add_attribute<std::vector<Dart>, IncidenceGraph::Vertex>(ig, "vertex_boundary_edge_dart");

	return true;
}

bool add_cmap2_attributes(CMap2& m2, IG_M2Attributes& m2Attribs)
{
	m2Attribs.vertex_position = add_attribute<Vec3, CMap2::Vertex>(m2, "position");
	m2Attribs.dual_vertex_graph_branch =
		add_attribute<std::pair<IncidenceGraph::Vertex, IncidenceGraph::Edge>, CMap2::Vertex>(m2, "graph_branch");
	m2Attribs.volume_center = add_attribute<Vec3, CMap2::Volume>(m2, "center");
	m2Attribs.volume_igvertex = add_attribute<IncidenceGraph::Vertex, CMap2::Volume>(m2, "igvertex");
	m2Attribs.edge_mid = add_attribute<Vec3, CMap2::Edge>(m2, "edge_mid");
	m2Attribs.halfedge_volume_connection = add_attribute<Dart, CMap2::HalfEdge>(m2, "volume_connection");
	// m2Attribs.ortho_scaffold = add_attribute<CMap2*, CMap2::Volume>(m2, "ortho_scaffold");

	return true;
}

bool add_cmap3_attributes_igh(CMap3& m3, IG_M3Attributes& m3Attribs)
{
	m3Attribs.vertex_position = cgogn::add_attribute<Vec3, CMap3::Vertex>(m3, "position");

	return true;
}

bool compute_faces_geometry(const IncidenceGraph& ig, const IncidenceGraphData& incidenceGraph_data,
							IG_GAttributes& igAttribs)
{
	geometry::compute_normal<IncidenceGraph::Face>(ig, igAttribs.vertex_position.get(), igAttribs.face_normal.get());
	return true;
}

/*****************************************************************************/
/* utils                                                                     */
/*****************************************************************************/

std::vector<IncidenceGraph::Vertex> get_branch_vertices(const IncidenceGraph& ig, const Branch& b)
{
	using Vertex = IncidenceGraph::Vertex;
	using Edge = IncidenceGraph::Edge;

	const std::pair<Vertex, Edge>& branch_start = b.first;
	const std::pair<Vertex, Edge>& branch_end = b.second;

	Edge current_edge = branch_start.second;
	const std::pair<Vertex, Vertex>& inc_verts = (*ig.edge_incident_vertices_)[current_edge.index_];
	Vertex next_vertex = branch_start.first == inc_verts.first ? inc_verts.second : inc_verts.first;

	std::vector<Vertex> branch_vertices = {branch_start.first};
	while (current_edge != branch_end.second)
	{
		branch_vertices.push_back(next_vertex);
		const std::vector<Edge>& inc_edges = (*ig.vertex_incident_edges_)[next_vertex.index_];
		current_edge = inc_edges[0] == current_edge ? inc_edges[1] : inc_edges[0];
		const std::pair<Vertex, Vertex>& inc_verts = (*ig.edge_incident_vertices_)[current_edge.index_];
		next_vertex = inc_verts.first == next_vertex ? inc_verts.second : inc_verts.first;
	}
	branch_vertices.push_back(next_vertex);
	return branch_vertices;
}

// Get all edges belonging to the same leaflet as e that are also incident to v
std::vector<IncidenceGraph::Edge> get_incident_leaflet_edges(const IncidenceGraph& ig, IncidenceGraph::Vertex v,
															 IncidenceGraph::Edge e)
{
	std::vector<IncidenceGraph::Edge> edges = {e};
	CellMarker<IncidenceGraph, IncidenceGraph::Edge> edge_marker(ig);
	edge_marker.mark(e);
	for (uint32 i = 0; i < edges.size(); ++i)
	{
		foreach_adjacent_edge_through_face(ig, edges[i], [&](IncidenceGraph::Edge ie) -> bool {
			if (!edge_marker.is_marked(ie))
			{
				edge_marker.mark(ie);
				const std::pair<IncidenceGraph::Vertex, IncidenceGraph::Vertex>& inc_verts =
					(*ig.edge_incident_vertices_)[ie.index_];
				if (inc_verts.first == v || inc_verts.second == v)
					edges.push_back(ie);
			}
			return true;
		});
	}
	return edges;
}

IncidenceGraph::Edge get_shared_edge(const IncidenceGraph& ig, IncidenceGraph::Vertex v0, IncidenceGraph::Vertex v1)
{
	const std::vector<IncidenceGraph::Edge>& inc_e0 = (*ig.vertex_incident_edges_)[v0.index_];
	const std::vector<IncidenceGraph::Edge>& inc_e1 = (*ig.vertex_incident_edges_)[v1.index_];
	for (IncidenceGraph::Edge e0 : inc_e0)
	{
		for (IncidenceGraph::Edge e1 : inc_e1)
		{
			if (e0 == e1)
				return e0;
		}
	}
	return IncidenceGraph::Edge();
}

Dart add_chunk(CMap3& m3)
{
	std::vector<Dart> D = {add_prism(static_cast<CMap2&>(m3), 4).dart, add_prism(static_cast<CMap2&>(m3), 4).dart,
						   add_prism(static_cast<CMap2&>(m3), 4).dart, add_prism(static_cast<CMap2&>(m3), 4).dart};

	sew_volumes_igh(m3, phi2(m3, D[0]), phi<-1, 2>(m3, D[1]));
	sew_volumes_igh(m3, phi2(m3, D[1]), phi<-1, 2>(m3, D[2]));
	sew_volumes_igh(m3, phi2(m3, D[2]), phi<-1, 2>(m3, D[3]));
	sew_volumes_igh(m3, phi2(m3, D[3]), phi<-1, 2>(m3, D[0]));

	return D[0];
}

Dart add_plate(CMap3& m3)
{
	Dart d0 = add_prism(static_cast<CMap2&>(m3), 4).dart;
	Dart d1 = add_prism(static_cast<CMap2&>(m3), 4).dart;

	sew_volumes_igh(m3, d0, d1);

	return d0;
}

void index_volume_cells_igh(CMap2& m, CMap2::Volume vol)
{
	foreach_dart_of_orbit(m, vol, [&](Dart d) -> bool {
		if (is_indexed<CMap2::HalfEdge>(m))
		{
			if (index_of(m, CMap2::HalfEdge(d)) == INVALID_INDEX)
				set_index(m, CMap2::HalfEdge(d), new_index<CMap2::HalfEdge>(m));
		}
		if (is_indexed<CMap2::Vertex>(m))
		{
			if (index_of(m, CMap2::Vertex(d)) == INVALID_INDEX)
				set_index(m, CMap2::Vertex(d), new_index<CMap2::Vertex>(m));
		}
		if (is_indexed<CMap2::Edge>(m))
		{
			if (index_of(m, CMap2::Edge(d)) == INVALID_INDEX)
				set_index(m, CMap2::Edge(d), new_index<CMap2::Edge>(m));
		}
		if (is_indexed<CMap2::Face>(m))
		{
			if (index_of(m, CMap2::Face(d)) == INVALID_INDEX)
				set_index(m, CMap2::Face(d), new_index<CMap2::Face>(m));
		}
		return true;
	});
	if (is_indexed<CMap2::Volume>(m))
	{
		if (index_of(m, vol) == INVALID_INDEX)
			set_index(m, vol, new_index<CMap2::Volume>(m));
	}
}

void sew_volumes_igh(CMap3& m, Dart d0, Dart d1)
{
	cgogn_message_assert(codegree(m, CMap3::Face(d0)) == codegree(m, CMap3::Face(d1)),
						 "The faces to sew do not have the same codegree");
	Dart it0 = d0;
	Dart it1 = d1;
	do
	{
		cgogn_message_assert(phi3(m, it0) == it0 && phi3(m, it1) == it1, "The faces to sew are already sewn");
		phi3_sew(m, it0, it1);
		it0 = phi1(m, it0);
		it1 = phi_1(m, it1);
	} while (it0 != d0);
}

Mat3 rmf_step(const Vec3& x0, const Vec3& x1, const Mat3& U0, const Vec3& t1)
{
	// U0 = (r0, s0, t0)
	Vec3 r0 = U0.col(0);
	Vec3 t0 = U0.col(2);

	// compute reflexion of r0
	Vec3 v0 = x1 - x0;
	Scalar c0 = v0.dot(v0);
	Vec3 rl = r0 - (2 / c0) * (v0.dot(r0)) * v0;
	Vec3 tl = t0 - (2 / c0) * (v0.dot(t0)) * v0;
	Vec3 v1 = t1 - tl;
	Scalar c1 = v1.dot(v1);
	Vec3 r1 = rl - (2 / c1) * (v1.dot(rl)) * v1;
	Vec3 s1 = t1.cross(r1);

	Mat3 U1;
	U1.col(0) = r1;
	U1.col(1) = s1;
	U1.col(2) = t1;
	return U1;
}

Vec3 slerp_igh(Vec3 A, Vec3 B, Scalar alpha, bool in)
{
	Scalar phi = geometry::angle(A, B) - (in ? 0 : 2 * M_PI);
	Scalar s0 = std::sin(phi * (1 - alpha));
	Scalar s1 = std::sin(phi * alpha);
	Scalar s2 = std::sin(phi);
	Vec3 sl = A * (s0 / s2);
	sl += B * (s1 / s2);
	return sl;
}

Scalar angle_on_sphere_igh(Vec3 A, Vec3 B, Vec3 C)
{
	Vec3 sB = slerp_igh(A, B, 0.01, true);
	Vec3 sC = slerp_igh(A, C, 0.01, true);
	Vec3 AB = sB - A;
	Vec3 AC = sC - A;
	return geometry::angle(AB, AC);
}

Scalar edge_max_angle_igh(CMap2& m2, CMap2::Edge e, IG_M2Attributes& m2Attribs)
{
	Dart ed0 = e.dart;
	Dart ed1 = phi2(m2, ed0);

	Vec3 A = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(ed0));
	Vec3 B = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(ed1));
	Vec3 C = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi_1(m2, ed0)));

	Scalar a0 = angle_on_sphere_igh(A, B, C);
	C = B;
	B = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi1(m2, ed1)));
	a0 += angle_on_sphere_igh(A, B, C);

	A = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(ed1));
	B = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(ed0));
	C = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi_1(m2, ed1)));

	Scalar a1 = angle_on_sphere_igh(A, B, C);
	C = B;
	B = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi1(m2, ed0)));
	a1 += angle_on_sphere_igh(A, B, C);

	return std::max(a0, a1);
}

Scalar min_cut_angle_igh(CMap2& m2, CMap2::Vertex v0, CMap2::Vertex v1, IG_M2Attributes& m2Attribs)
{
	Vec3 A = value<Vec3>(m2, m2Attribs.vertex_position, v0);
	Vec3 B = value<Vec3>(m2, m2Attribs.vertex_position, v1);
	Vec3 C = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi_1(m2, v0.dart)));
	Scalar a0 = angle_on_sphere_igh(A, B, C);

	C = B;
	B = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi1(m2, v1.dart)));
	Scalar a1 = angle_on_sphere_igh(A, B, C);

	A = value<Vec3>(m2, m2Attribs.vertex_position, v1);
	B = value<Vec3>(m2, m2Attribs.vertex_position, v0);
	C = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi_1(m2, v1.dart)));

	Scalar a2 = angle_on_sphere_igh(A, B, C);
	C = B;
	B = value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi1(m2, v0.dart)));
	Scalar a3 = angle_on_sphere_igh(A, B, C);

	return std::min({a0, a1, a2, a3});
}

Vec3 spherical_barycenter_igh(std::vector<Vec3>& points, uint32 iterations)
{
	uint32 nb_pts = uint32(points.size());
	std::vector<Vec3> pts0(nb_pts);
	std::vector<Vec3> pts1(nb_pts);

	for (uint32 i = 0; i < nb_pts; ++i)
		pts0[i] = points[i];

	std::vector<Vec3> readvec = pts0;
	std::vector<Vec3> writevec = pts1;
	std::vector<Vec3> tempvec;

	for (uint32 it = 0; it < iterations; ++it)
	{
		for (uint32 i = 0; i < nb_pts; ++i)
		{
			writevec[i] = slerp_igh(readvec[i], readvec[(i + 1) % nb_pts], Scalar(0.5), true);
		}

		tempvec = writevec;
		writevec = readvec;
		readvec = tempvec;
	}

	Vec3 bary = {0, 0, 0};
	for (uint32 i = 0; i < nb_pts; ++i)
	{
		bary += readvec[i];
	}
	bary /= nb_pts;

	Vec3 normal = {0, 0, 0};
	for (uint32 i = 0; i < nb_pts; ++i)
	{
		normal += points[i].cross(points[(i + 1) % nb_pts]);
	}
	normal.normalize();

	if (normal.dot(bary) < 0)
		bary *= -1;

	return bary;
}

bool dijkstra_topo_igh(CMap2& m2, CMap2::Vertex v0, std::shared_ptr<CMap2::Attribute<Dart>> previous,
					   std::shared_ptr<CMap2::Attribute<uint32>> dist)
{
	foreach_incident_vertex(m2, CMap2::Volume(v0.dart), [&](CMap2::Vertex v) -> bool {
		value<Dart>(m2, previous, v) = Dart();
		value<uint32>(m2, dist, v) = UINT_MAX;
		return true;
	});

	DartMarker visited(m2);

	std::vector<CMap2::Vertex> vertices = {v0};
	value<uint32>(m2, dist, v0) = 0;
	CMap2::Vertex v_act;
	uint32 dist_act;
	while (vertices.size())
	{
		v_act = vertices[vertices.size() - 1];
		vertices.pop_back();
		dist_act = value<uint32>(m2, dist, v_act) + 1;

		std::vector<CMap2::Vertex> neighbors;
		foreach_dart_of_orbit(m2, v_act, [&](Dart d) -> bool {
			if (!visited.is_marked(d))
			{
				Dart d2 = phi2(m2, d);
				visited.mark(d);
				visited.mark(d2);

				uint32 dist_2 = value<uint32>(m2, dist, CMap2::Vertex(d2));
				if (dist_act < dist_2)
				{
					value<uint32>(m2, dist, CMap2::Vertex(d2)) = dist_act;
					value<Dart>(m2, previous, CMap2::Vertex(d2)) = d;
					neighbors.push_back(CMap2::Vertex(d2));
				}
			}

			return true;
		});

		vertices.insert(vertices.begin(), neighbors.begin(), neighbors.end());
	}

	return false;
}

Dart convex_hull_around_vertex(const IncidenceGraph& ig, IncidenceGraph::Vertex v, CMap2& m2,
							   IG_M2Attributes& m2Attribs, std::vector<Vec3>& Ppos,
							   std::vector<IncidenceGraph::Edge>& Pev)
{
	std::vector<uint32> Pid;
	Pid.reserve(Ppos.size());

	uint32 i = 0;

	for (IncidenceGraph::Edge e : Pev)
	{
		uint32 vertex_id = new_index<CMap2::Vertex>(m2);
		Pid.push_back(vertex_id);
		(*m2Attribs.vertex_position)[vertex_id] = Ppos[i++];
		(*m2Attribs.dual_vertex_graph_branch)[vertex_id] = {v, e};
	}

	std::vector<uint32> faces_nb_vertices;
	std::vector<uint32> faces_vertex_indices;

	std::vector<uint32> indices;
	for (uint32 i = 0; i < Ppos.size() - 2; ++i)
	{
		for (uint32 j = i + 1; j < Ppos.size() - 1; ++j)
		{
			for (uint32 k = j + 1; k < Ppos.size(); ++k)
			{
				Vec3 t0 = Ppos[j] - Ppos[i];
				Vec3 t1 = Ppos[k] - Ppos[i];
				Vec3 n = t0.cross(t1);
				int sign = 0;

				for (uint32 m = 0; m < Ppos.size(); ++m)
				{
					if (m == i || m == j || m == k)
						continue;

					Vec3 vec = Ppos[m] - Ppos[i];
					Scalar d = vec.dot(n);

					if (!sign)
						sign = (d < 0 ? -1 : 1);
					else
					{
						if (sign != (d < 0 ? -1 : 1))
						{
							sign = 0;
							break;
						}
					}
				}

				if (sign != 0)
				{
					if (sign == 1)
					{
						indices.push_back(Pid[j]);
						indices.push_back(Pid[i]);
						indices.push_back(Pid[k]);
					}
					else
					{
						indices.push_back(Pid[i]);
						indices.push_back(Pid[j]);
						indices.push_back(Pid[k]);
					}

					faces_nb_vertices.push_back(3);
					faces_vertex_indices.insert(faces_vertex_indices.end(), indices.begin(), indices.end());
					indices.clear();
				}
			}
		}
	}

	auto darts_per_vertex = add_attribute<std::vector<Dart>, CMap2::Vertex>(m2, "__darts_per_vertex");

	uint32 faces_vertex_index = 0u;
	std::vector<uint32> vertices_buffer;
	vertices_buffer.reserve(16u);

	Dart volume_dart;
	std::vector<Dart> all_darts;
	for (std::size_t i = 0u, end = faces_nb_vertices.size(); i < end; ++i)
	{
		uint32 nbv = faces_nb_vertices[i];

		vertices_buffer.clear();

		for (uint32 j = 0u; j < nbv; ++j)
		{
			uint32 idx = faces_vertex_indices[faces_vertex_index++];
			vertices_buffer.push_back(idx);
		}

		if (nbv > 2u)
		{
			CMap1::Face f = add_face(static_cast<CMap1&>(m2), nbv, false);
			Dart d = f.dart;
			for (uint32 j = 0u; j < nbv; ++j)
			{
				all_darts.push_back(d);
				const uint32 vertex_index = vertices_buffer[j];
				set_index<CMap2::Vertex>(m2, d, vertex_index);
				(*darts_per_vertex)[vertex_index].push_back(d);
				d = phi1(m2, d);
			}
		}
	}

	for (Dart d : all_darts)
	{
		if (phi2(m2, d) == d)
		{
			uint32 vertex_index = index_of(m2, CMap2::Vertex(d));

			const std::vector<Dart>& next_vertex_darts =
				value<std::vector<Dart>>(m2, darts_per_vertex, CMap2::Vertex(phi1(m2, d)));
			bool phi2_found = false;

			for (auto it = next_vertex_darts.begin(); it != next_vertex_darts.end() && !phi2_found; ++it)
			{
				if (index_of(m2, CMap2::Vertex(phi1(m2, *it))) == vertex_index)
				{
					if (phi2(m2, *it) == *it)
					{
						phi2_sew(m2, d, *it);
						phi2_found = true;
					}
				}
			}
		}
	}

	remove_attribute<CMap2::Vertex>(m2, darts_per_vertex);
	return all_darts[0];
}

void dualize_volume(CMap2& m, CMap2::Volume vol, IG_M2Attributes& m2Attribs, const IncidenceGraph& ig,
					IG_GAttributes& igAttribs)
{
	// set the new phi1
	foreach_dart_of_orbit(m, vol, [&](Dart d) -> bool {
		Dart dd = phi2(m, phi_1(m, d));
		(*(m.phi1_))[d.index] = dd;
		return true;
	});

	// set the new phi_1
	foreach_dart_of_orbit(m, vol, [&](Dart d) -> bool {
		Dart next = phi1(m, d);
		(*(m.phi_1_))[next.index] = d;
		return true;
	});

	DartMarkerStore<CMap2> face_marker(m);
	const IncidenceGraph::Vertex gv = value<IncidenceGraph::Vertex>(m, m2Attribs.volume_igvertex, vol);
	const Vec3& center = value<Vec3>(ig, igAttribs.vertex_position, gv);
	foreach_dart_of_orbit(m, vol, [&](Dart d) -> bool {
		if (!face_marker.is_marked(d))
		{
			CMap2::Face f(d);
			foreach_dart_of_orbit(m, f, [&](Dart d) -> bool {
				face_marker.mark(d);
				return true;
			});
			if (is_indexed<CMap2::Face>(m))
				set_index(m, f, new_index<CMap2::Face>(m)); // give a new index to the face

			// darts of the new face orbit are darts of the old dual vertex
			// get the outgoing graph branch from the old dual vertex
			std::pair<IncidenceGraph::Vertex, IncidenceGraph::Edge>& branch =
				value<std::pair<IncidenceGraph::Vertex, IncidenceGraph::Edge>>(m, m2Attribs.dual_vertex_graph_branch,
																			   CMap2::Vertex(d));
			if (degree(ig, branch.second) == 0) // dangling edge
			{
				std::pair<IncidenceGraph::Vertex, IncidenceGraph::Vertex>& inc_verts =
					(*ig.edge_incident_vertices_)[branch.second.index_];
				if (inc_verts.first == branch.first)
					value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, branch.second).first = d;
				else
					value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, branch.second).second = d;
			}
			else // leaflet
			{
				std::vector<IncidenceGraph::Edge> leaf_edges =
					get_incident_leaflet_edges(ig, branch.first, branch.second);
				for (IncidenceGraph::Edge le : leaf_edges)
				{
					std::pair<IncidenceGraph::Vertex, IncidenceGraph::Vertex>& inc_verts =
						(*ig.edge_incident_vertices_)[le.index_];
					if (inc_verts.first == branch.first)
						value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, le).first = d;
					else
						value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, le).second = d;
				}
			}

			return true;
		}
		return true;
	});

	DartMarkerStore<CMap2> vertex_marker(m);
	foreach_dart_of_orbit(m, vol, [&](Dart d) -> bool {
		if (!vertex_marker.is_marked(d))
		{
			CMap2::Vertex v(d);

			std::vector<Vec3> points;
			foreach_dart_of_orbit(m, v, [&](Dart d) -> bool {
				vertex_marker.mark(d);
				points.push_back(value<Vec3>(m, m2Attribs.vertex_position, CMap2::Vertex(d)) - center);
				return true;
			});
			Vec3 b;
			if (points.size() == 2)
				b = slerp_igh(points[0], points[1], 0.5, true) + center;
			else
				b = spherical_barycenter_igh(points, 10) + center;

			set_index(m, v, new_index<CMap2::Vertex>(m));	  // give a new index to the vertex
			value<Vec3>(m, m2Attribs.vertex_position, v) = b; // set the position to the computed position
			return true;
		}
		return true;
	});
}

Dart remesh_igh(CMap2& m2, CMap2::Volume vol, IG_M2Attributes& m2Attribs)
{
	std::cout << "remeshing contact surface " << std::endl;

	Dart vol_dart = vol.dart;

	std::vector<CMap2::Vertex> valence_sup4;
	std::vector<CMap2::Vertex> valence_3;

	foreach_incident_vertex(m2, vol, [&](CMap2::Vertex v) -> bool {
		uint32 valence = degree(m2, v);
		if (valence == 3)
			valence_3.push_back(v);
		if (valence > 4)
			valence_sup4.push_back(v);

		return true;
	});

	if (valence_sup4.size())
	{
		auto edge_angle_max = add_attribute<Scalar, CMap2::Edge>(m2, "angle_max");

		std::vector<CMap2::Edge> edges_n_n;
		std::vector<CMap2::Edge> edges_n_4;
		std::vector<CMap2::Edge> candidate_edges;

		do
		{
			edges_n_n.clear();
			edges_n_4.clear();
			foreach_incident_edge(m2, vol, [&](CMap2::Edge e) -> bool {
				auto vertices = incident_vertices(m2, e);

				uint32 deg_0 = degree(m2, vertices[0]);
				uint32 deg_1 = degree(m2, vertices[1]);
				uint32 deg_min = deg_0 < deg_1 ? deg_0 : deg_1;

				if (deg_min > 4)
				{
					edges_n_n.push_back(e);
				}
				else if (deg_min == 4 && deg_0 + deg_1 > 8)
				{
					edges_n_4.push_back(e);
				}

				return true;
			});

			candidate_edges = edges_n_n.size() ? edges_n_n : edges_n_4;
			if (!candidate_edges.size())
				break;
			for (CMap2::Edge e : candidate_edges)
			{
				value<Scalar>(m2, edge_angle_max, e) = edge_max_angle_igh(m2, e, m2Attribs);
			}

			std::sort(candidate_edges.begin(), candidate_edges.end(), [&](CMap2::Edge e0, CMap2::Edge e1) -> bool {
				return value<Scalar>(m2, edge_angle_max, e0) < value<Scalar>(m2, edge_angle_max, e1);
			});

			CMap2::Edge prime_edge = candidate_edges[0];
			auto neigh_vertices = incident_vertices(m2, prime_edge);

			vol_dart = phi_1(m2, prime_edge.dart);
			merge_incident_faces(m2, prime_edge, true);

		} while (candidate_edges.size() > 1);
		remove_attribute<CMap2::Edge>(m2, edge_angle_max);
	}

	std::vector<std::pair<std::pair<CMap2::Vertex, CMap2::Vertex>, Scalar>> candidate_vertices_pairs;

	valence_3.clear();
	foreach_incident_vertex(m2, vol, [&](CMap2::Vertex v) -> bool {
		if (degree(m2, v) == 3)
			valence_3.push_back(v);
		return true;
	});
	while (valence_3.size())
	{
		candidate_vertices_pairs.clear();
		std::vector<CMap2::Vertex> verts_3;
		verts_3.reserve(10);
		foreach_incident_face(m2, vol, [&](CMap2::Face f) -> bool {
			verts_3.clear();
			foreach_incident_vertex(m2, f, [&](CMap2::Vertex v) -> bool {
				if (degree(m2, v) == 3)
					verts_3.push_back(v);

				return true;
			});

			if (verts_3.size() >= 2)
				for (uint32 i = 0; i < verts_3.size() - 1; ++i)
				{
					for (uint32 j = i + 1; j < verts_3.size(); ++j)
					{
						candidate_vertices_pairs.push_back(
							{{verts_3[i], verts_3[j]}, min_cut_angle_igh(m2, verts_3[i], verts_3[j], m2Attribs)});
					}
				}
			return true;
		});

		if (candidate_vertices_pairs.size())
		{
			std::sort(candidate_vertices_pairs.begin(), candidate_vertices_pairs.end(),
					  [&](auto pair0, auto pair1) -> bool { return pair0.second < pair1.second; });

			auto pair = candidate_vertices_pairs[0].first;

			cut_face(m2, pair.first, pair.second);
		}
		else
		{
			uint32 max_path_length = 0;
			std::shared_ptr<CMap2::Attribute<Dart>> max_previous;
			std::shared_ptr<CMap2::Attribute<uint32>> max_dist;
			std::pair<CMap2::Vertex, CMap2::Vertex> max_shortest_path;
			for (uint32 i = 0; i < valence_3.size(); ++i)
			{
				CMap2::Vertex v = valence_3[i];
				std::shared_ptr<CMap2::Attribute<Dart>> previous =
					add_attribute<Dart, CMap2::Vertex>(m2, "previous" + std::to_string(v.dart.index));
				std::shared_ptr<CMap2::Attribute<uint32>> dist =
					add_attribute<uint32, CMap2::Vertex>(m2, "dist" + std::to_string(v.dart.index));

				dijkstra_topo_igh(m2, v, previous, dist);

				uint32 curr_min = UINT32_MAX;
				CMap2::Vertex curr_min_vert;
				for (CMap2::Vertex v2 : valence_3)
				{
					if (index_of(m2, v) == index_of(m2, v2))
						continue;

					uint32 new_min = value<uint32>(m2, dist, CMap2::Vertex(v2));
					if (new_min < curr_min)
					{

						curr_min = new_min; // value<uint32>(m2, dist, CMap2::Vertex(v2)); ????
						curr_min_vert = v2;
					}
				}

				if (curr_min > max_path_length)
				{
					if (max_previous)
						remove_attribute<CMap2::Vertex>(m2, max_previous);
					if (max_dist)
						remove_attribute<CMap2::Vertex>(m2, max_dist);

					max_previous = previous;
					max_dist = dist;
					max_shortest_path = {v, curr_min_vert};
					max_path_length = curr_min;
				}
			}

			// building path
			std::vector<Dart> path;
			path.reserve(max_path_length);
			Dart dp = value<Dart>(m2, max_previous, max_shortest_path.second);
			do
			{
				path.push_back(dp);
				dp = value<Dart>(m2, max_previous, CMap2::Vertex(dp));
			} while (dp.index != INVALID_INDEX);

			if (!(path.size() % 2))
			{
				std::vector<CMap2::Edge> first_edges;
				first_edges.push_back(CMap2::Edge(phi_1(m2, path.front())));
				first_edges.push_back(CMap2::Edge(phi2(m2, phi1(m2, path.front()))));
				first_edges.push_back(CMap2::Edge(phi1(m2, path.back())));
				first_edges.push_back(CMap2::Edge(phi2(m2, phi_1(m2, path.back()))));

				std::vector<std::pair<Scalar, uint32>> angles;
				angles.push_back({edge_max_angle_igh(m2, first_edges[0], m2Attribs), 0});
				angles.push_back({edge_max_angle_igh(m2, first_edges[1], m2Attribs), 1});
				angles.push_back({edge_max_angle_igh(m2, first_edges[2], m2Attribs), 2});
				angles.push_back({edge_max_angle_igh(m2, first_edges[3], m2Attribs), 3});

				std::sort(angles.begin(), angles.end(),
						  [&](std::pair<Scalar, uint32> p0, std::pair<Scalar, uint32> p1) -> bool {
							  return p0.first < p1.first;
						  });

				CMap2::Vertex v0, v1;
				CMap2::Edge e0;
				switch (angles[0].second)
				{
				case 0:
					v0 = CMap2::Vertex(phi1(m2, path.front()));
					v1 = CMap2::Vertex(phi_1(m2, path.front()));
					e0 = CMap2::Edge(v1.dart);
					path.erase(path.begin());
					break;
				case 1:
					v0 = CMap2::Vertex(phi2(m2, path.front()));
					v1 = CMap2::Vertex(phi<2, 1, 1>(m2, path.front()));
					e0 = CMap2::Edge(phi<2, 1>(m2, path.front()));
					path.erase(path.begin());
					break;

				case 2:
					v0 = CMap2::Vertex(path.back());
					v1 = CMap2::Vertex(phi<1, 1>(m2, path.back()));
					e0 = CMap2::Edge(phi1(m2, path.back()));
					path.pop_back();
					break;

				case 3:
					v0 = CMap2::Vertex(phi_1(m2, phi2(m2, path.back())));
					v1 = CMap2::Vertex(phi<2, 1>(m2, path.back()));
					e0 = CMap2::Edge(v0.dart);
					path.pop_back();
					break;
				default:
					break;
				}

				vol_dart = cut_face(m2, v0, v1, true).dart;
				merge_incident_faces(m2, e0, true);
			}

			for (uint32 i = 0; i < path.size(); ++i)
			{
				if (!(i % 2))
					cut_face(m2, CMap2::Vertex(path[i]), CMap2::Vertex(phi1(m2, path[i])), true);
				else
				{
					vol_dart = phi1(m2, path[i]);
					merge_incident_faces(m2, CMap2::Edge(path[i]), true);
				}
			}
		}

		valence_3.clear();
		foreach_incident_vertex(m2, vol, [&](CMap2::Vertex v) -> bool {
			if (degree(m2, v) == 3)
				valence_3.push_back(v);
			return true;
		});
	}

	// foreach_incident_vertex(m2, vol, [&](CMap2::Vertex v) -> bool {
	// 	// Vec3 p = value<Vec3>(m2, m2Attribs.vertex_position, v);
	// 	std::cout << degree(m2, v) << std::endl;
	// 	return true;
	// });

	return vol_dart;
}

/*****************************************************************************/
/* contact surfaces generation                                               */
/*****************************************************************************/

bool build_contact_surfaces(const IncidenceGraph& ig, IG_GAttributes& igAttribs,
							IncidenceGraphData& incidenceGraph_data, CMap2& m2, IG_M2Attributes& m2Attribs)
{
	bool success = true;

	foreach_cell(ig, [&](IncidenceGraph::Vertex v) -> bool {
		std::pair<uint32, uint32> info = pseudo_degree(ig, v); // { nb_leaflets, nb_edges }

		if (info.first == 0 && info.second == 1) // branch extremities
			build_contact_surface_1(ig, igAttribs, m2, m2Attribs, v);
		else if (info.first + info.second == 2) // efjunctures, ffjunctures or joint (eejunctures)
			build_contact_surface_2(ig, igAttribs, m2, m2Attribs, v, info.first);
		else if (info.first + info.second > 2) // intersections
			build_contact_surface_n(ig, igAttribs, m2, m2Attribs, v, info.first, info.second);

		return true;
	});

	return success;
}

void build_contact_surface_1(const IncidenceGraph& ig, IG_GAttributes& igAttribs, CMap2& m2, IG_M2Attributes& m2Attribs,
							 IncidenceGraph::Vertex v)
{
	Dart d = add_face(m2, 4, true).dart;

	value<Dart>(ig, igAttribs.vertex_contact_surface, v) = d;
	value<IncidenceGraph::Vertex>(m2, m2Attribs.volume_igvertex, CMap2::Volume(d)) = v;

	IncidenceGraph::Edge e = (*ig.vertex_incident_edges_)[v.index_][0];
	const std::pair<IncidenceGraph::Vertex, IncidenceGraph::Vertex>& inc_verts =
		(*ig.edge_incident_vertices_)[e.index_];
	if (inc_verts.first == v)
		value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, e).first = d;
	else
		value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, e).second = d;
}

void build_contact_surface_2(const IncidenceGraph& ig, IG_GAttributes& igAttribs, CMap2& m2, IG_M2Attributes& m2Attribs,
							 IncidenceGraph::Vertex v, uint32 nb_leaflets)
{
	Dart d0 = add_face(static_cast<CMap1&>(m2), 4, false).dart;
	Dart d1 = add_face(static_cast<CMap1&>(m2), 4, false).dart;
	phi2_sew(m2, d0, d1);
	phi2_sew(m2, phi1(m2, d0), phi_1(m2, d1));
	phi2_sew(m2, phi<1, 1>(m2, d0), phi<1, 1>(m2, d1));
	phi2_sew(m2, phi_1(m2, d0), phi1(m2, d1));
	d1 = phi<1, 2>(m2, d0);

	index_volume_cells_igh(m2, CMap2::Volume(d0));

	value<Dart>(ig, igAttribs.vertex_contact_surface, v) = d0;
	value<IncidenceGraph::Vertex>(m2, m2Attribs.volume_igvertex, CMap2::Volume(d0)) = v;

	if (nb_leaflets == 0) // eejuncture
	{
		const std::vector<IncidenceGraph::Edge>& inc_edges = (*ig.vertex_incident_edges_)[v.index_];

		// first edge is bound to d0
		const std::pair<IncidenceGraph::Vertex, IncidenceGraph::Vertex>& inc_verts_0 =
			(*ig.edge_incident_vertices_)[inc_edges[0].index_];
		if (v == inc_verts_0.first)
			value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, inc_edges[0]).first = d0;
		else
			value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, inc_edges[0]).second = d0;

		// second edge is bound to d1
		const std::pair<IncidenceGraph::Vertex, IncidenceGraph::Vertex>& inc_verts_1 =
			(*ig.edge_incident_vertices_)[inc_edges[1].index_];
		if (v == inc_verts_1.first)
			value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, inc_edges[1]).first = d1;
		else
			value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, inc_edges[1]).second = d1;
	}
	else if (nb_leaflets == 1) // efjuncture
	{
		foreach_incident_edge(ig, v, [&](IncidenceGraph::Edge e) -> bool {
			const std::pair<IncidenceGraph::Vertex, IncidenceGraph::Vertex>& inc_verts =
				(*ig.edge_incident_vertices_)[e.index_];

			// dangling edge is bound to d0, leaflet edges are bound to d1
			Dart d = degree(ig, e) == 0 ? d0 : d1;
			if (inc_verts.first == v)
				value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, e).first = d;
			else
				value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, e).second = d;

			return true;
		});
	}
	else // nb_leaflets == 2 // ffjuncture
	{
		Dart d = d0; // first leaflet edges are bound to d0
		CellMarker<IncidenceGraph, IncidenceGraph::Edge> edge_marker(ig);
		foreach_incident_edge(ig, v, [&](IncidenceGraph::Edge ie) -> bool {
			if (edge_marker.is_marked(ie))
				return true;
			std::vector<IncidenceGraph::Edge> inc_leaflet_edges = get_incident_leaflet_edges(ig, v, ie);
			for (IncidenceGraph::Edge ile : inc_leaflet_edges)
			{
				edge_marker.mark(ile);
				const std::pair<IncidenceGraph::Vertex, IncidenceGraph::Vertex>& inc_verts =
					(*ig.edge_incident_vertices_)[ile.index_];
				if (inc_verts.first == v)
					value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, ile).first = d;
				else
					value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, ile).second = d;
			}
			d = d1; // second leaflet edges are bound to d1
			return true;
		});
	}
}

void build_contact_surface_orange(const IncidenceGraph& ig, IG_GAttributes& igAttribs, CMap2& m2,
								  IG_M2Attributes& m2Attribs, IncidenceGraph::Vertex v, std::vector<Vec3>& Ppos,
								  std::vector<IncidenceGraph::Edge>& Pev)
{
	uint32 nbf = Ppos.size();

	// create the topology of the surface
	std::vector<Dart> faces;
	faces.reserve(nbf);
	for (uint32 i = 0; i < nbf; ++i)
		faces.push_back(add_face(static_cast<CMap1&>(m2), 4, false).dart);
	for (uint32 i = 0; i < nbf; ++i)
	{
		Dart d = faces[i];
		Dart e = faces[(i + 1) % nbf];
		phi2_sew(m2, d, phi1(m2, e));
		phi2_sew(m2, phi_1(m2, d), phi<1, 1>(m2, e));
	}

	index_volume_cells_igh(m2, CMap2::Volume(faces[0]));

	value<Dart>(ig, igAttribs.vertex_contact_surface, v) = faces[0];

	const Vec3& center = value<Vec3>(ig, igAttribs.vertex_position, v);

	// get the best fitting plane normal and build a local frame based on this normal and first branch
	std::pair<Vec3, Scalar> plane = geometry::plane_fitting(Ppos);
	Vec3 L1 = (Ppos[0] - center).normalized();
	geometry::project_on_plane(L1, plane.first, 0);
	Vec3 L3 = plane.first;
	Vec3 L2 = L3.cross(L1);
	Mat3 L;
	L << L1[0], L1[1], L1[2], L2[0], L2[1], L2[2], L3[0], L3[1], L3[2];

	// sort the incident branches on the plane in CW order
	std::vector<uint32> permutation(nbf);
	std::iota(permutation.begin(), permutation.end(), 0);
	std::sort(permutation.begin(), permutation.end(), [&](uint32 i, uint32 j) -> bool {
		Vec3 proji = (Ppos[i] - center).normalized();
		geometry::project_on_plane(proji, plane.first, 0);
		proji = L * proji;
		Scalar anglei = atan2(proji[1], proji[0]);

		Vec3 projj = (Ppos[j] - center).normalized();
		geometry::project_on_plane(projj, plane.first, 0);
		projj = L * projj;
		Scalar anglej = atan2(projj[1], projj[0]);

		if (anglei >= 0)
		{
			if (anglej >= 0)
				return anglei < anglej;
			else
				return true;
		}
		else
		{
			if (anglej >= 0)
				return false;
			else
				return anglei < anglej;
		}
	});

	// apply the permutation to branches point & edge
	std::vector<Vec3> sorted_Ppos(nbf);
	std::vector<IncidenceGraph::Edge> sorted_Pev(nbf);
	std::transform(permutation.begin(), permutation.end(), sorted_Ppos.begin(), [&](uint32 i) { return Ppos[i]; });
	std::transform(permutation.begin(), permutation.end(), sorted_Pev.begin(), [&](uint32 i) { return Pev[i]; });

	// put the geometry on the surface mesh vertices
	Vec3 Q1 = center + plane.first;
	Vec3 Q2 = center - plane.first;
	for (uint32 i = 0; i < nbf; ++i)
	{
		IncidenceGraph::Edge e = sorted_Pev[i];
		if (degree(ig, e) == 0) // dangling edge
		{
			const std::pair<IncidenceGraph::Vertex, IncidenceGraph::Vertex>& inc_verts =
				(*ig.edge_incident_vertices_)[e.index_];
			if (inc_verts.first == v)
				value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, e).first = faces[i];
			else
				value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, e).second = faces[i];
		}
		else // leaflet
		{
			std::vector<IncidenceGraph::Edge> leaf_edges = get_incident_leaflet_edges(ig, v, e);
			for (IncidenceGraph::Edge le : leaf_edges)
			{
				const std::pair<IncidenceGraph::Vertex, IncidenceGraph::Vertex>& inc_verts =
					(*ig.edge_incident_vertices_)[le.index_];
				if (inc_verts.first == v)
					value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, le).first = faces[i];
				else
					value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, le).second = faces[i];
			}
		}

		value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(faces[i])) =
			center + (sorted_Ppos[(i + 1) % nbf] - sorted_Ppos[i]).normalized().cross(plane.first);
	}
	value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi1(m2, faces[0]))) = Q1;
	value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi_1(m2, faces[0]))) = Q2;

	value<IncidenceGraph::Vertex>(m2, m2Attribs.volume_igvertex, CMap2::Volume(faces[0])) = v;
}

void build_contact_surface_n(const IncidenceGraph& ig, IG_GAttributes& igAttribs, CMap2& m2, IG_M2Attributes& m2Attribs,
							 IncidenceGraph::Vertex v, uint32 nb_leaflets, uint32 nb_edges)
{
	using Vertex = IncidenceGraph::Vertex;
	using Edge = IncidenceGraph::Edge;

	uint32 nbf = nb_leaflets + nb_edges;

	const Vec3& center = value<Vec3>(ig, igAttribs.vertex_position, v);
	std::vector<Vec3> Ppos;
	Ppos.reserve(nbf);

	std::vector<IncidenceGraph::Edge> Pev;
	Pev.reserve(nbf);

	CellMarker<IncidenceGraph, IncidenceGraph::Edge> edge_marker(ig);
	foreach_incident_edge(ig, v, [&](IncidenceGraph::Edge e) -> bool {
		// filter out edges of already managed leaflets
		if (edge_marker.is_marked(e))
			return true;

		Vec3 p = {0, 0, 0};
		if (degree(ig, e) == 0) // dangling edge
		{
			edge_marker.mark(e);
			const std::pair<IncidenceGraph::Vertex, IncidenceGraph::Vertex>& inc_verts =
				(*ig.edge_incident_vertices_)[e.index_];
			p = value<Vec3>(ig, igAttribs.vertex_position, inc_verts.first == v ? inc_verts.second : inc_verts.first);
		}
		else // leaflet
		{
			std::vector<IncidenceGraph::Edge> leaf_edges = get_incident_leaflet_edges(ig, v, e);
			for (IncidenceGraph::Edge le : leaf_edges)
			{
				edge_marker.mark(le);
				const std::pair<IncidenceGraph::Vertex, IncidenceGraph::Vertex>& inc_verts =
					(*ig.edge_incident_vertices_)[le.index_];
				p += value<Vec3>(ig, igAttribs.vertex_position,
								 inc_verts.first == v ? inc_verts.second : inc_verts.first);
			}
			p /= leaf_edges.size();
		}

		geometry::project_on_sphere(p, center, 1);
		Ppos.push_back(p);
		Pev.push_back(e);
		return true;
	});

	std::pair<Vec3, Scalar> plane = geometry::plane_fitting(Ppos);
	bool planar = true;
	for (const Vec3& p : Ppos)
	{
		Scalar dist = geometry::distance_plane_point(plane.first, plane.second, p);
		if (dist > 0.05)
		{
			planar = false;
			break;
		}
	}
	if (planar)
	{
		build_contact_surface_orange(ig, igAttribs, m2, m2Attribs, v, Ppos, Pev);
		return;
	}

	// compute the n points on the sphere
	// generate Delaunay mesh from the n points
	// store the graph branch on their respective delaunay vertex (m2Attribs.dual_vertex_graph_branch)
	// modify connectivity until all vertices are valence 4
	// call dualize_volume

	Dart vol_dart = convex_hull_around_vertex(ig, v, m2, m2Attribs, Ppos, Pev);

	index_volume_cells_igh(m2, CMap2::Volume(vol_dart));
	value<IncidenceGraph::Vertex>(m2, m2Attribs.volume_igvertex, CMap2::Volume(vol_dart)) = v;

	vol_dart = remesh_igh(m2, CMap2::Volume(vol_dart), m2Attribs);
	dualize_volume(m2, CMap2::Volume(vol_dart), m2Attribs, ig, igAttribs);

	value<Dart>(ig, igAttribs.vertex_contact_surface, v) = vol_dart;

	return;
}

/*****************************************************************************/
/* frames initialization & propagation                                       */
/*****************************************************************************/

bool create_intersection_frames(const IncidenceGraph& ig, IG_GAttributes& igAttribs, const IncidenceGraphData& igData,
								CMap2& m2, IG_M2Attributes m2Attribs)
{
	bool res = true;
	for (IncidenceGraph::Vertex v : igData.intersections)
		res = create_intersection_frame_n(ig, igAttribs, m2, m2Attribs, v);
	for (IncidenceGraph::Vertex v : igData.efjunctures)
		res = create_ef_frame(ig, igAttribs, v);
	for (IncidenceGraph::Vertex v : igData.ffjunctures)
		res = create_ff_frame(ig, igAttribs, v);
	return res;
}

bool create_intersection_frame_n(const IncidenceGraph& ig, IG_GAttributes& igAttribs, CMap2& m2,
								 IG_M2Attributes& m2Attribs, IncidenceGraph::Vertex v)
{
	const Vec3& center = value<Vec3>(ig, igAttribs.vertex_position, v);

	CellMarker<IncidenceGraph, IncidenceGraph::Edge> edge_marker(ig);
	foreach_incident_edge(ig, v, [&](IncidenceGraph::Edge e) -> bool {
		// filter out edges of already managed leaflets
		if (edge_marker.is_marked(e))
			return true;

		std::pair<Dart, Dart>& inc_d = value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, e);

		if (degree(ig, e) == 0) // dangling edge
		{
			edge_marker.mark(e);
			const std::pair<IncidenceGraph::Vertex, IncidenceGraph::Vertex>& inc_verts =
				(*ig.edge_incident_vertices_)[e.index_];

			const Vec3& p =
				value<Vec3>(ig, igAttribs.vertex_position, inc_verts.first == v ? inc_verts.second : inc_verts.first);

			Dart d0 = inc_verts.first == v ? inc_d.first : inc_d.second;
			Dart d1 = phi<1, 1>(m2, d0);

			Vec3 T = (p - center).normalized();
			Vec3 diag = (value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(d1)) -
						 value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(d0)))
							.normalized();
			Vec3 R = diag.cross(T).normalized();
			Vec3 S = T.cross(R).normalized();

			Mat3& f = inc_verts.first == v ? value<std::pair<Mat3, Mat3>>(ig, igAttribs.halfedge_frame, e).first
										   : value<std::pair<Mat3, Mat3>>(ig, igAttribs.halfedge_frame, e).second;
			f.col(0) = R;
			f.col(1) = S;
			f.col(2) = T;
		}
		else // leaflet
		{
			std::vector<IncidenceGraph::Edge> leaf_edges = get_incident_leaflet_edges(ig, v, e);

			Vec3 leaflet_dir{0, 0, 0};
			for (IncidenceGraph::Edge le : leaf_edges)
			{
				edge_marker.mark(le);
				const std::pair<IncidenceGraph::Vertex, IncidenceGraph::Vertex>& inc_verts =
					(*ig.edge_incident_vertices_)[le.index_];
				IncidenceGraph::Vertex iv = inc_verts.first == v ? inc_verts.second : inc_verts.first;
				leaflet_dir += (value<Vec3>(ig, igAttribs.vertex_position, iv) - center).normalized();
			}
			leaflet_dir.normalize();

			const std::pair<IncidenceGraph::Vertex, IncidenceGraph::Vertex>& inc_verts =
				(*ig.edge_incident_vertices_)[e.index_];
			Dart d0 = inc_verts.first == v ? inc_d.first : inc_d.second;
			Dart d1 = phi<1, 1>(m2, d0);

			Vec3 T = leaflet_dir;
			Vec3 diag = (value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(d1)) -
						 value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(d0)))
							.normalized();
			Vec3 R = diag.cross(T).normalized();
			Vec3 S = T.cross(R).normalized();

			for (IncidenceGraph::Edge el : leaf_edges)
			{
				edge_marker.mark(e);
				const std::pair<IncidenceGraph::Vertex, IncidenceGraph::Vertex>& inc_verts =
					(*ig.edge_incident_vertices_)[el.index_];
				Mat3& f = inc_verts.first == v ? value<std::pair<Mat3, Mat3>>(ig, igAttribs.halfedge_frame, e).first
											   : value<std::pair<Mat3, Mat3>>(ig, igAttribs.halfedge_frame, e).second;
				f.col(0) = R;
				f.col(1) = S;
				f.col(2) = T;
			}
		}

		return true;
	});

	return true;
}

bool create_ef_frame(const IncidenceGraph& ig, IG_GAttributes& igAttribs, IncidenceGraph::Vertex v)
{
	using Vertex = IncidenceGraph::Vertex;
	using Edge = IncidenceGraph::Edge;
	using Face = IncidenceGraph::Face;

	const Vec3& center = value<Vec3>(ig, igAttribs.vertex_position, v);

	Vec3 leaflet_dir{0, 0, 0};
	Vec3 leaflet_normal{0, 0, 0};
	Edge branch_edge;
	foreach_incident_edge(ig, v, [&](Edge e) -> bool {
		uint32 d = degree(ig, e);
		if (d > 0)
		{
			const std::pair<Vertex, Vertex>& inc_verts = (*ig.edge_incident_vertices_)[e.index_];
			Vertex iv = inc_verts.first == v ? inc_verts.second : inc_verts.first;
			leaflet_dir += (value<Vec3>(ig, igAttribs.vertex_position, iv) - center).normalized();
			const std::vector<Face>& inc_faces = (*ig.edge_incident_faces_)[e.index_];
			for (Face iface : inc_faces)
				leaflet_normal += value<Vec3>(ig, igAttribs.face_normal, iface);
		}
		else if (d == 0)
			branch_edge = e;
		return true;
	});
	leaflet_dir.normalize();
	leaflet_normal.normalize();

	Vec3 T = leaflet_dir;
	Vec3 R = leaflet_normal;
	Vec3 S = T.cross(R).normalized();
	R = S.cross(T).normalized();

	foreach_incident_edge(ig, v, [&](Edge e) -> bool {
		const std::pair<Vertex, Vertex>& inc_verts = (*ig.edge_incident_vertices_)[e.index_];
		Mat3& f = inc_verts.first == v ? value<std::pair<Mat3, Mat3>>(ig, igAttribs.halfedge_frame, e).first
									   : value<std::pair<Mat3, Mat3>>(ig, igAttribs.halfedge_frame, e).second;
		if (degree(ig, e) == 0) // dangling edge
		{
			f.col(0) = R;
			f.col(1) = -S;
			f.col(2) = -T;
		}
		else // leaflet edge
		{
			f.col(0) = R;
			f.col(1) = S;
			f.col(2) = T;
		}
		return true;
	});

	return true;
}

bool create_ff_frame(const IncidenceGraph& ig, IG_GAttributes& igAttribs, IncidenceGraph::Vertex v)
{
	using Vertex = IncidenceGraph::Vertex;
	using Edge = IncidenceGraph::Edge;
	using Face = IncidenceGraph::Face;

	const Vec3& center = value<Vec3>(ig, igAttribs.vertex_position, v);

	std::vector<Edge> leaflet_edges1;
	Vec3 leaflet_dir1, leaflet_dir2;

	std::vector<Edge> leaflet_edges2;
	Vec3 leaflet_normal1, leaflet_normal2;

	bool first = true;

	CellMarker<IncidenceGraph, Edge> edge_marker(ig);
	foreach_incident_edge(ig, v, [&](Edge e) -> bool {
		// filters out edges of already managed leaflets
		if (edge_marker.is_marked(e))
			return true;

		std::vector<Edge> leaf_edges = get_incident_leaflet_edges(ig, v, e);
		Vec3 leaflet_dir{0, 0, 0};
		Vec3 leaflet_normal{0, 0, 0};
		for (Edge le : leaf_edges)
		{
			edge_marker.mark(le);

			const std::pair<Vertex, Vertex>& inc_verts = (*ig.edge_incident_vertices_)[le.index_];
			Vertex iv = inc_verts.first == v ? inc_verts.second : inc_verts.first;
			leaflet_dir += (value<Vec3>(ig, igAttribs.vertex_position, iv) - center).normalized();
			const std::vector<Face>& inc_faces = (*ig.edge_incident_faces_)[le.index_];
			for (Face iface : inc_faces)
				leaflet_normal += value<Vec3>(ig, igAttribs.face_normal, iface);
		}
		leaflet_dir.normalize();
		leaflet_normal.normalize();

		if (first)
		{
			leaflet_edges1 = leaf_edges;
			leaflet_dir1 = leaflet_dir;
			leaflet_normal1 = leaflet_normal;
			first = false;
		}
		else
		{
			leaflet_edges2 = leaf_edges;
			leaflet_dir2 = leaflet_dir;
			leaflet_normal2 = leaflet_normal;
		}

		return true;
	});

	Vec3 meandir = (leaflet_dir1 - leaflet_dir2) / 2.0;
	meandir.normalize();
	if (leaflet_normal2.dot(leaflet_normal1) < 0)
		leaflet_normal2 *= -1.0;
	Vec3 meannormal = (leaflet_normal1 + leaflet_normal2) / 2.0;
	meannormal.normalize();

	Vec3 T = meandir;
	Vec3 R = meannormal;
	Vec3 S = T.cross(R).normalized();
	R = S.cross(T).normalized();

	for (Edge e : leaflet_edges1)
	{
		const std::pair<Vertex, Vertex>& inc_verts = (*ig.edge_incident_vertices_)[e.index_];
		Mat3& f = inc_verts.first == v ? value<std::pair<Mat3, Mat3>>(ig, igAttribs.halfedge_frame, e).first
									   : value<std::pair<Mat3, Mat3>>(ig, igAttribs.halfedge_frame, e).second;
		f.col(0) = R;
		f.col(1) = S;
		f.col(2) = T;
	}
	for (Edge e : leaflet_edges2)
	{
		const std::pair<Vertex, Vertex>& inc_verts = (*ig.edge_incident_vertices_)[e.index_];
		Mat3& f = inc_verts.first == v ? value<std::pair<Mat3, Mat3>>(ig, igAttribs.halfedge_frame, e).first
									   : value<std::pair<Mat3, Mat3>>(ig, igAttribs.halfedge_frame, e).second;
		f.col(0) = R;
		f.col(1) = -S;
		f.col(2) = -T;
	}

	return true;
}

bool create_extremity_frame(const IncidenceGraph& ig, IG_GAttributes& igAttribs, IncidenceGraph::Vertex v)
{
	const Vec3& center = value<Vec3>(ig, igAttribs.vertex_position, v);
	IncidenceGraph::Edge e = (*ig.vertex_incident_edges_)[v.index_][0];
	const std::pair<IncidenceGraph::Vertex, IncidenceGraph::Vertex>& inc_verts =
		(*ig.edge_incident_vertices_)[e.index_];

	Vec3 T =
		(value<Vec3>(ig, igAttribs.vertex_position, inc_verts.first == v ? inc_verts.second : inc_verts.first) - center)
			.normalized();
	Vec3 temp = Vec3(T[1], -T[0], T[2]);

	Vec3 R = temp.cross(T).normalized();
	Vec3 S = T.cross(R).normalized();

	Mat3& f = inc_verts.first == v ? value<std::pair<Mat3, Mat3>>(ig, igAttribs.halfedge_frame, e).first
								   : value<std::pair<Mat3, Mat3>>(ig, igAttribs.halfedge_frame, e).second;
	f.col(0) = R;
	f.col(1) = S;
	f.col(2) = T;

	return true;
}

bool propagate_frames(const IncidenceGraph& ig, IG_GAttributes& igAttribs, const IncidenceGraphData& igData, CMap2& m2)
{
	using Vertex = IncidenceGraph::Vertex;
	using Edge = IncidenceGraph::Edge;

	for (const Branch& branch : igData.branches)
	{
		std::vector<Vertex> branch_vertices = get_branch_vertices(ig, branch);
		uint32 d0 = degree(ig, branch_vertices.front());
		uint32 d1 = degree(ig, branch_vertices.back());

		if (d0 > 1 && d1 == 1)
			propagate_frame_n_1(ig, igAttribs, branch_vertices);
		else if (d0 == 1 && d1 > 1)
		{
			std::reverse(branch_vertices.begin(), branch_vertices.end());
			propagate_frame_n_1(ig, igAttribs, branch_vertices);
		}
		else if (d0 > 1 && d1 > 1)
			propagate_frame_n_n(ig, igAttribs, m2, branch_vertices);
		else
		{
			// if(degree(ig, branch_vertices[0]) == 1) erreur
			// create_extremity_frame(ig, igAttribs, branch_vertices[0]);
		}
	}

	return true;
}

void propagate_frame_n_1(const IncidenceGraph& ig, IG_GAttributes& igAttribs,
						 std::vector<IncidenceGraph::Vertex>& branch_vertices)
{
	std::vector<Vec3> tangents;
	for (uint32 i = 1; i < branch_vertices.size(); ++i)
	{
		const Vec3& p0 = value<Vec3>(ig, igAttribs.vertex_position, branch_vertices[i - 1]);
		const Vec3& p1 = value<Vec3>(ig, igAttribs.vertex_position, branch_vertices[i]);
		Vec3 t = (p1 - p0).normalized();
		if (i < branch_vertices.size() - 1)
			t += (value<Vec3>(ig, igAttribs.vertex_position, branch_vertices[i + 1]) - p1).normalized();
		t.normalize();
		tangents.push_back(t);
	}

	for (uint32 i = 1; i < branch_vertices.size(); ++i)
	{
		IncidenceGraph::Vertex v0 = branch_vertices[i - 1];
		IncidenceGraph::Vertex v1 = branch_vertices[i];

		const Vec3& p0 = value<Vec3>(ig, igAttribs.vertex_position, v0);
		const Vec3& p1 = value<Vec3>(ig, igAttribs.vertex_position, v1);

		IncidenceGraph::Edge e0 = get_shared_edge(ig, v0, v1);
		std::pair<IncidenceGraph::Vertex, IncidenceGraph::Vertex> inc_verts = (*ig.edge_incident_vertices_)[e0.index_];

		std::pair<Mat3, Mat3>& frames = value<std::pair<Mat3, Mat3>>(ig, igAttribs.halfedge_frame, e0);
		Mat3& U0 = inc_verts.first == v0 ? frames.first : frames.second;
		Mat3 U1 = rmf_step(p0, p1, U0, tangents[i - 1]);
		Mat3& U_1 = inc_verts.first == v1 ? frames.first : frames.second;
		U_1.col(0) = U1.col(0);
		U_1.col(1) = -U1.col(1);
		U_1.col(2) = -U1.col(2);

		if (i < branch_vertices.size() - 1)
		{
			const std::vector<IncidenceGraph::Edge>& inc_edges = (*ig.vertex_incident_edges_)[v1.index_];
			IncidenceGraph::Edge e1 = inc_edges[inc_edges[0] == e0 ? 1 : 0];
			inc_verts = (*ig.edge_incident_vertices_)[e1.index_];
			if (inc_verts.first == v1)
				value<std::pair<Mat3, Mat3>>(ig, igAttribs.halfedge_frame, e1).first = U1;
			else
				value<std::pair<Mat3, Mat3>>(ig, igAttribs.halfedge_frame, e1).second = U1;
		}
	}
}

bool propagate_frame_n_n(const IncidenceGraph& ig, IG_GAttributes& igAttribs, CMap2& m2,
						 std::vector<IncidenceGraph::Vertex>& branch_vertices)
{
	Mat3 lastU;
	std::vector<IncidenceGraph::Vertex> inc_verts;

	std::vector<Vec3> tangents;
	tangents.push_back(Vec3());
	for (uint32 i = 1; i < branch_vertices.size(); ++i)
	{
		const Vec3& p0 = value<Vec3>(ig, igAttribs.vertex_position, branch_vertices[i - 1]);
		const Vec3& p1 = value<Vec3>(ig, igAttribs.vertex_position, branch_vertices[i]);
		Vec3 t = (p1 - p0).normalized();
		if (i < branch_vertices.size() - 1)
			t += (value<Vec3>(ig, igAttribs.vertex_position, branch_vertices[i + 1]) - p1).normalized();
		t.normalize();
		tangents.push_back(t);
	}

	/// propagating to the second to last vertex

	for (uint32 i = 1; i < branch_vertices.size() - 1; ++i)
	{
		IncidenceGraph::Vertex v0 = branch_vertices[i - 1];
		IncidenceGraph::Vertex v1 = branch_vertices[i];
		IncidenceGraph::Edge e0 = get_shared_edge(ig, v0, v1);

		inc_verts = incident_vertices(ig, e0);
		uint32 vid0 = v0 == inc_verts[0] ? 0 : 1;
		uint32 vid1 = v1 == inc_verts[0] ? 0 : 1;

		const Vec3& p0 = value<Vec3>(ig, igAttribs.vertex_position, branch_vertices[i - 1]);
		const Vec3& p1 = value<Vec3>(ig, igAttribs.vertex_position, branch_vertices[i]);

		std::pair<Mat3, Mat3>& frames = value<std::pair<Mat3, Mat3>>(ig, igAttribs.halfedge_frame, e0);
		Mat3& U0 = vid0 == 0 ? frames.first : frames.second;
		lastU = U0;
		Mat3 U1 = rmf_step(p0, p1, U0, tangents[i]);
		Mat3& U_1 = vid1 == 0 ? frames.first : frames.second;
		U_1.col(0) = U1.col(0);
		U_1.col(1) = -U1.col(1);
		U_1.col(2) = -U1.col(2);

		if (i < branch_vertices.size() - 1)
		{
			std::vector<IncidenceGraph::Edge> inc_edges = incident_edges(ig, v1);
			IncidenceGraph::Edge e1 = inc_edges[e0 == inc_edges[0] ? 1 : 0];
			inc_verts = incident_vertices(ig, e1);
			vid1 = v1 == inc_verts[0] ? 0 : 1;
			if (vid1 == 0)
				value<std::pair<Mat3, Mat3>>(ig, igAttribs.halfedge_frame, e1).first = U1;
			else
				value<std::pair<Mat3, Mat3>>(ig, igAttribs.halfedge_frame, e1).second = U1;
		}
	}

	uint8 nb_vertices = branch_vertices.size();
	Mat3 endU = rmf_step(value<Vec3>(ig, igAttribs.vertex_position, branch_vertices[nb_vertices - 2]),
						 value<Vec3>(ig, igAttribs.vertex_position, branch_vertices[nb_vertices - 1]), lastU,
						 tangents[nb_vertices - 1]);

	endU.col(1) = -endU.col(1);
	endU.col(2) = -endU.col(2);

	IncidenceGraph::Edge e0 = get_shared_edge(ig, branch_vertices[nb_vertices - 2], branch_vertices[nb_vertices - 1]);
	inc_verts = incident_vertices(ig, e0);
	uint32 vidEnd = branch_vertices[nb_vertices - 1] == inc_verts[0] ? 0 : 1;
	Mat3 UE;
	if (vidEnd == 0)
		UE = value<std::pair<Mat3, Mat3>>(ig, igAttribs.halfedge_frame, e0).first;
	else
		UE = value<std::pair<Mat3, Mat3>>(ig, igAttribs.halfedge_frame, e0).second;

	Vec3 X = (endU.col(0) + endU.col(1)).normalized();
	Vec3 RE = UE.col(0), SE = UE.col(1);
	bool A = (RE.dot(X) >= 0);
	bool B = (SE.dot(X) >= 0);

	uint32 nb_shifts = 0;
	if (!A && B)
		nb_shifts = 1;
	else if (!A && !B)
		nb_shifts = 2;
	else if (A && !B)
		nb_shifts = 3;

	if (nb_shifts > 0)
	{
		// shift_frame(UE, nb_shifts);
		for (uint32 i = 0; i < nb_shifts; ++i)
		{
			Vec3 R = UE.col(1);
			Vec3 S = -UE.col(0);
			UE.col(0) = R;
			UE.col(1) = S;
		}

		if (vidEnd == 0)
			value<std::pair<Mat3, Mat3>>(ig, igAttribs.halfedge_frame, e0).first = UE;
		else
			value<std::pair<Mat3, Mat3>>(ig, igAttribs.halfedge_frame, e0).second = UE;

		Dart csface;
		if (vidEnd == 0)
			csface = value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, e0).first;
		else
			csface = value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, e0).second;

		for (uint32 i = 0; i < nb_shifts; ++i)
			csface = phi1(m2, csface);

		if (vidEnd == 0)
			value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, e0).first = csface;
		else
			value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, e0).second = csface;
	}

	if (nb_vertices > 2)
	{
		Scalar cos0 = UE.col(0).dot(endU.col(0));
		Scalar cos1 = UE.col(1).dot(endU.col(0));
		Scalar angle = std::acos(cos0) * (cos1 > 0 ? 1 : -1);
		Scalar angle_step = -angle / Scalar(nb_vertices - 1);

		uint32 step = 0;
		for (uint32 i = 1; i < branch_vertices.size() - 2; ++i)
		{
			++step;
			IncidenceGraph::Vertex v0 = branch_vertices[i - 1];
			IncidenceGraph::Vertex v1 = branch_vertices[i];
			IncidenceGraph::Vertex v2 = branch_vertices[i + 1];
			IncidenceGraph::Edge e0 = get_shared_edge(ig, v0, v1);
			IncidenceGraph::Edge e1 = get_shared_edge(ig, v1, v2);

			const std::pair<IncidenceGraph::Vertex, IncidenceGraph::Vertex>& inc_verts0 =
				(*ig.edge_incident_vertices_)[e0.index_];
			const std::pair<IncidenceGraph::Vertex, IncidenceGraph::Vertex>& inc_verts1 =
				(*ig.edge_incident_vertices_)[e1.index_];

			uint32 vid10 = v1 == inc_verts0.first ? 0 : 1;
			uint32 vid11 = v1 == inc_verts1.first ? 0 : 1;

			Mat3 U, U_;
			if (vid11 == 0)
				U = value<std::pair<Mat3, Mat3>>(ig, igAttribs.halfedge_frame, e1).first;
			else
				U = value<std::pair<Mat3, Mat3>>(ig, igAttribs.halfedge_frame, e1).second;

			Eigen::AngleAxisd rot(-angle_step * step, U.col(2));
			U.col(0) = rot * U.col(0);
			U.col(1) = U.col(2).cross(U.col(0));

			U_.col(0) = U.col(0);
			U_.col(1) = -U.col(1);
			U_.col(2) = -U.col(2);

			if (vid10 == 0)
				value<std::pair<Mat3, Mat3>>(ig, igAttribs.halfedge_frame, e0).first = U_;
			else
				value<std::pair<Mat3, Mat3>>(ig, igAttribs.halfedge_frame, e0).second = U_;

			if (vid11 == 0)
				value<std::pair<Mat3, Mat3>>(ig, igAttribs.halfedge_frame, e1).first = U;
			else
				value<std::pair<Mat3, Mat3>>(ig, igAttribs.halfedge_frame, e1).second = U;
		}
	}

	return true;
}

/*****************************************************************************/
/* contact surfaces geometry                                                 */
/*****************************************************************************/

bool set_contact_surfaces_geometry(const IncidenceGraph& ig, IG_GAttributes& igAttribs, CMap2& m2,
								   IG_M2Attributes& m2Attribs)
{
	foreach_cell(ig, [&](IncidenceGraph::Vertex v) -> bool {
		std::pair<uint32, uint32> info = pseudo_degree(ig, v); // { nb_leaflets, nb_edges }

		// inner (0,0) & boundary (1,0) leaflet vertices do not have a contact surface
		if (info.first < 2 && info.second == 0)
			return true;

		// TODO: get radius from attribute
		Scalar radius = value<Scalar>(ig, igAttribs.vertex_radius, v);
		const Vec3& center = value<Vec3>(ig, igAttribs.vertex_position, v);

		// set the center position of the contact surface (CMap2 volume attribute)
		CMap2::Volume contact_surface(value<Dart>(ig, igAttribs.vertex_contact_surface, v));
		value<Vec3>(m2, m2Attribs.volume_center, contact_surface) = center;

		// extremities (0,1), joints (eejunctures) (0,2), efjunctures (1,1) & ffjunctures (2,0)
		// -> simple contact surfaces with 4 vertices
		// -> get the computed frame from the first incident edge (they all have the info)
		if (info.first + info.second <= 2)
		{
			IncidenceGraph::Edge e = (*ig.vertex_incident_edges_)[v.index_][0];
			const std::pair<IncidenceGraph::Vertex, IncidenceGraph::Vertex>& inc_verts =
				(*ig.edge_incident_vertices_)[e.index_];

			Dart csf = inc_verts.first == v
						   ? value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, e).first
						   : value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, e).second;
			const Mat3& frame = inc_verts.first == v
									? value<std::pair<Mat3, Mat3>>(ig, igAttribs.halfedge_frame, e).first
									: value<std::pair<Mat3, Mat3>>(ig, igAttribs.halfedge_frame, e).second;

			value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(csf)) = center - frame.col(1) * radius;
			value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi1(m2, csf))) = center + frame.col(0) * radius;
			value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi<1, 1>(m2, csf))) =
				center + frame.col(1) * radius;
			value<Vec3>(m2, m2Attribs.vertex_position, CMap2::Vertex(phi_1(m2, csf))) = center - frame.col(0) * radius;
		}
		// intersections
		// -> contact surface vertices already have positions
		// -> just project them on the sphere with the right radius
		else if (info.first + info.second > 2)
		{
			foreach_incident_vertex(m2, contact_surface, [&](CMap2::Vertex v2) -> bool {
				Vec3 pos = value<Vec3>(m2, m2Attribs.vertex_position, v2);
				geometry::project_on_sphere(pos, center, radius);
				value<Vec3>(m2, m2Attribs.vertex_position, v2) = pos;
				return true;
			});
		}

		foreach_incident_edge(m2, contact_surface, [&](CMap2::Edge e) -> bool {
			std::vector<CMap2::Vertex> vertices = incident_vertices(m2, e);
			Vec3 mid = 0.5 * (value<Vec3>(m2, m2Attribs.vertex_position, vertices[0]) +
							  value<Vec3>(m2, m2Attribs.vertex_position, vertices[1]));
			geometry::project_on_sphere(mid, center, radius);
			value<Vec3>(m2, m2Attribs.edge_mid, e) = mid;
			return true;
		});

		return true;
	});
	return true;
}

/*****************************************************************************/
/* volume mesh generation                                                    */
/*****************************************************************************/

bool build_volumes(IncidenceGraph& ig, IG_GAttributes& igAttribs, IncidenceGraphData& incidenceGraph_data, CMap2& m2,
				   IG_M2Attributes& m2Attribs, CMap3& m3)
{
	// insert_ortho_chunks(ig, igAttribs, m2, m2Attribs, m3);

	bool success = true;
	for (auto& leaflet : incidenceGraph_data.leaflets)
	{
		success = build_leaflet(ig, igAttribs, m2, m2Attribs, m3, leaflet);
		if (!success)
			break;
	}

	if (!success)
		return false;

	foreach_cell(ig, [&](IncidenceGraph::Edge e) -> bool {
		if (degree(ig, e) == 0)
			success = build_branch_section(ig, igAttribs, m2, m2Attribs, m3, e);
		return success;
	});

	return success;
}

bool build_branch_section(IncidenceGraph& ig, IG_GAttributes& igAttribs, CMap2& m2, IG_M2Attributes& m2Attribs,
						  CMap3& m3, IncidenceGraph::Edge e)
{
	Dart m2f0 = value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, e).first;
	Dart m2f1 = value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, e).second;

	std::vector<Dart> F0 = {m2f0, phi1(m2, m2f0), phi<1, 1>(m2, m2f0), phi_1(m2, m2f0)};
	std::vector<Dart> F1 = {m2f1, phi1(m2, m2f1), phi<1, 1>(m2, m2f1), phi_1(m2, m2f1)};

	Dart m3d = add_chunk(m3);
	std::vector<Dart> D0 = {m3d, phi<2, 3, 2, 1>(m3, m3d), phi<2, 3, 2, 1, 2, 3, 2, 1>(m3, m3d),
							phi<-1, 2, 3, 2>(m3, m3d)};
	std::vector<Dart> D1 = {phi<2, 1, 1, 2>(m3, D0[0]), phi<2, 1, 1, 2>(m3, D0[1]), phi<2, 1, 1, 2>(m3, D0[2]),
							phi<2, 1, 1, 2>(m3, D0[3])};

	value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F0[0])) = phi1(m3, D0[0]);
	value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F0[1])) = phi1(m3, D0[1]);
	value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F0[2])) = phi1(m3, D0[2]);
	value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F0[3])) = phi1(m3, D0[3]);

	value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F1[0])) = phi<1, 1>(m3, D1[1]);
	value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F1[1])) = phi<1, 1>(m3, D1[0]);
	value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F1[2])) = phi<1, 1>(m3, D1[3]);
	value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F1[3])) = phi<1, 1>(m3, D1[2]);

	// value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_volume_connection, e).first = D0[0];
	// value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_volume_connection, e).second = phi1(m3, D1[0]);

	return true;
}

bool build_leaflet(IncidenceGraph& ig, IG_GAttributes& igAttribs, CMap2& m2, IG_M2Attributes& m2Attribs, CMap3& m3,
				   const std::vector<IncidenceGraph::Face>& leaflet)
{
	CellMarker<IncidenceGraph, IncidenceGraph::Edge> edge_marker(ig);

	std::vector<IncidenceGraph::Edge> inside_edges;
	std::vector<IncidenceGraph::Edge> boundary_edges;
	std::vector<uint8> boundary_edges_dir;

	// build plate elements for each face
	for (IncidenceGraph::Face f : leaflet)
	{
		Dart d = add_plate(m3);

		const std::vector<IncidenceGraph::Edge>& inc_edges = (*ig.face_incident_edges_)[f.index_];
		const std::vector<uint8>& inc_edges_dir = (*ig.face_incident_edges_dir_)[f.index_];

		value<std::vector<Dart>>(ig, igAttribs.face_edge_dart, f).resize(4);
		for (uint32 i = 0; i < 4; ++i)
		{
			IncidenceGraph::Edge e = inc_edges[i];
			if (!edge_marker.is_marked(e))
			{
				edge_marker.mark(e);
				if (degree(ig, e) == 1)
				{
					boundary_edges.push_back(e);
					boundary_edges_dir.push_back(inc_edges_dir[i]);
				}
				else
					inside_edges.push_back(e);
			}

			// attach dart to face edge
			value<std::vector<Dart>>(ig, igAttribs.face_edge_dart, f)[i] = d;

			d = phi_1(m3, d);
		}
	}

	// all edges of the leaflet are marked with edge_marker now

	// connect inside plate elements together
	for (IncidenceGraph::Edge e : inside_edges)
	{
		const std::vector<IncidenceGraph::Face>& inc_faces = (*ig.edge_incident_faces_)[e.index_];
		uint32 eid0 = get_incident_edge_id(ig, inc_faces[0], e);
		uint32 eid1 = get_incident_edge_id(ig, inc_faces[1], e);

		Dart edge_dart0 = value<std::vector<Dart>>(ig, igAttribs.face_edge_dart, inc_faces[0])[eid0];
		Dart edge_dart1 = value<std::vector<Dart>>(ig, igAttribs.face_edge_dart, inc_faces[1])[eid1];

		sew_volumes_igh(m3, phi2(m3, edge_dart0), phi2(m3, edge_dart1));
		sew_volumes_igh(m3, phi<3, 2>(m3, edge_dart0), phi<3, 2>(m3, edge_dart1));
	}

	// build plate elements for each boundary edge and connect to incident interior volumes
	for (IncidenceGraph::Edge e : boundary_edges)
	{
		Dart d = add_plate(m3);

		IncidenceGraph::Face f = (*ig.edge_incident_faces_)[e.index_][0];
		uint32 eid = get_incident_edge_id(ig, f, e);
		Dart edge_dart = value<std::vector<Dart>>(ig, igAttribs.face_edge_dart, f)[eid];

		sew_volumes_igh(m3, phi2(m3, d), phi2(m3, edge_dart));
		sew_volumes_igh(m3, phi<3, 2>(m3, d), phi<3, 2>(m3, edge_dart));
	}

	// connect corner boundary elements together (1,0)
	// start from one boundary edge
	IncidenceGraph::Edge be = boundary_edges[0];
	std::pair<IncidenceGraph::Vertex, IncidenceGraph::Vertex> be_inc_verts = (*ig.edge_incident_vertices_)[be.index_];
	IncidenceGraph::Vertex bv = boundary_edges_dir[0] == 0 ? be_inc_verts.second : be_inc_verts.first;
	std::vector<IncidenceGraph::Edge> inc_leaflet_edges = get_incident_leaflet_edges(ig, bv, be);
	IncidenceGraph::Edge be_next;
	for (IncidenceGraph::Edge le : inc_leaflet_edges)
	{
		if (le != be && degree(ig, le) == 1)
			be_next = le;
	}
	do
	{
		std::pair<uint32, uint32> info = pseudo_degree(ig, bv); // { nb_leaflets, nb_edges }

		if (info.first == 1 && info.second == 0) // corners (1,0)
		{
			IncidenceGraph::Face f0 = (*ig.edge_incident_faces_)[be.index_][0];
			uint32 eid0 = get_incident_edge_id(ig, f0, be);
			// uint8 edge_dir_0 = (*ig.face_incident_edges_dir_)[f0][eid0];
			Dart edge_dart_0 = value<std::vector<Dart>>(ig, igAttribs.face_edge_dart, f0)[eid0];

			IncidenceGraph::Face f1 = (*ig.edge_incident_faces_)[be_next.index_][0];
			uint32 eid1 = get_incident_edge_id(ig, f1, be_next);
			// uint8 edge_dir_1 = (*ig.face_incident_edges_dir_)[f1][eid1];
			Dart edge_dart_1 = value<std::vector<Dart>>(ig, igAttribs.face_edge_dart, f1)[eid1];

			sew_volumes_igh(m3, phi<2, 3, -1, 2>(m3, edge_dart_0), phi<2, 3, 1, 2>(m3, edge_dart_1));
			sew_volumes_igh(m3, phi<3, 2, 3, 1, 2>(m3, edge_dart_0), phi<3, 2, 3, -1, 2>(m3, edge_dart_1));
		}
		else if (info.first == 1 && info.second == 1) // efjunctures (1,1)
		{
			Dart m2f = bv == be_inc_verts.first
						   ? value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, be).first
						   : value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, be).second;
			std::vector<Dart> F0 = {m2f, phi1(m2, m2f), phi<1, 1>(m2, m2f), phi_1(m2, m2f)};

			IncidenceGraph::Face f0 = (*ig.edge_incident_faces_)[be.index_][0];
			uint32 eid0 = get_incident_edge_id(ig, f0, be);
			// uint8 edge_dir_0 = (*ig.face_incident_edges_dir_)[f0][eid0];
			Dart edge_dart_0 = value<std::vector<Dart>>(ig, igAttribs.face_edge_dart, f0)[eid0];

			IncidenceGraph::Face f1 = (*ig.edge_incident_faces_)[be_next.index_][0];
			uint32 eid1 = get_incident_edge_id(ig, f1, be_next);
			// uint8 edge_dir_1 = (*ig.face_incident_edges_dir_)[f1][eid1];
			Dart edge_dart_1 = value<std::vector<Dart>>(ig, igAttribs.face_edge_dart, f1)[eid1];

			std::vector<Dart> D0 = {phi<2, 3, -1, 2>(m3, edge_dart_0), phi<2, 3, 1, 2, 1>(m3, edge_dart_1),
									phi<3, 2, 3, -1, 2>(m3, edge_dart_1), phi<3, 2, 3, 1, 2, 1>(m3, edge_dart_0)};

			value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F0[0])) = phi1(m3, D0[0]);
			value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F0[1])) = phi1(m3, D0[1]);
			value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F0[2])) = phi1(m3, D0[2]);
			value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F0[3])) = phi1(m3, D0[3]);
		}
		else if (info.first == 2 && info.second == 0) // ffjunctures (2,0)
		{
			Dart m2f = bv == be_inc_verts.first
						   ? value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, be).first
						   : value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, be).second;
			std::vector<Dart> F0 = {m2f, phi1(m2, m2f), phi<1, 1>(m2, m2f), phi_1(m2, m2f)};

			const Mat3& frame = bv == be_inc_verts.first
									? value<std::pair<Mat3, Mat3>>(ig, igAttribs.halfedge_frame, be).first
									: value<std::pair<Mat3, Mat3>>(ig, igAttribs.halfedge_frame, be).second;

			Vec3 leaflet_normal{0, 0, 0};
			for (IncidenceGraph::Edge le : inc_leaflet_edges)
			{
				const std::vector<IncidenceGraph::Face>& inc_faces = (*ig.edge_incident_faces_)[le.index_];
				for (IncidenceGraph::Face iface : inc_faces)
				{
					const Vec3& n = value<Vec3>(ig, igAttribs.face_normal, iface);
					leaflet_normal += n;
				}
			}
			leaflet_normal.normalize();

			// TODO: check other cases
			if (frame.col(0).dot(leaflet_normal) < 0)
				std::transform(F0.begin(), F0.end(), F0.begin(), [&](Dart d) { return phi<1, 1>(m2, d); });

			IncidenceGraph::Face f0 = (*ig.edge_incident_faces_)[be.index_][0];
			uint32 eid0 = get_incident_edge_id(ig, f0, be);
			// uint8 edge_dir_0 = (*ig.face_incident_edges_dir_)[f0][eid0];
			Dart edge_dart_0 = value<std::vector<Dart>>(ig, igAttribs.face_edge_dart, f0)[eid0];

			IncidenceGraph::Face f1 = (*ig.edge_incident_faces_)[be_next.index_][0];
			uint32 eid1 = get_incident_edge_id(ig, f1, be_next);
			// uint8 edge_dir_1 = (*ig.face_incident_edges_dir_)[f1][eid1];
			Dart edge_dart_1 = value<std::vector<Dart>>(ig, igAttribs.face_edge_dart, f1)[eid1];

			std::vector<Dart> D0 = {phi<2, 3, -1, 2>(m3, edge_dart_0), phi<2, 3, 1, 2, 1>(m3, edge_dart_1),
									phi<3, 2, 3, -1, 2>(m3, edge_dart_1), phi<3, 2, 3, 1, 2, 1>(m3, edge_dart_0)};

			value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F0[0])) = phi1(m3, D0[0]);
			value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F0[1])) = phi1(m3, D0[1]);
			value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F0[2])) = phi1(m3, D0[2]);
			value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F0[3])) = phi1(m3, D0[3]);
		}
		else // intersections
		{
			Dart m2f = bv == be_inc_verts.first
						   ? value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, be).first
						   : value<std::pair<Dart, Dart>>(ig, igAttribs.halfedge_contact_surface_face, be).second;
			std::vector<Dart> F0 = {m2f, phi1(m2, m2f), phi<1, 1>(m2, m2f), phi_1(m2, m2f)};

			const Mat3& frame = bv == be_inc_verts.first
									? value<std::pair<Mat3, Mat3>>(ig, igAttribs.halfedge_frame, be).first
									: value<std::pair<Mat3, Mat3>>(ig, igAttribs.halfedge_frame, be).second;

			Vec3 leaflet_normal{0, 0, 0};
			for (IncidenceGraph::Edge le : inc_leaflet_edges)
			{
				const std::vector<IncidenceGraph::Face>& inc_faces = (*ig.edge_incident_faces_)[le.index_];
				for (IncidenceGraph::Face iface : inc_faces)
				{
					const Vec3& n = value<Vec3>(ig, igAttribs.face_normal, iface);
					leaflet_normal += n;
				}
			}
			leaflet_normal.normalize();

			// TODO: check other cases
			if (frame.col(0).dot(leaflet_normal) < 0)
				std::transform(F0.begin(), F0.end(), F0.begin(), [&](Dart d) { return phi<1, 1>(m2, d); });

			IncidenceGraph::Face f0 = (*ig.edge_incident_faces_)[be.index_][0];
			uint32 eid0 = get_incident_edge_id(ig, f0, be);
			// uint8 edge_dir_0 = (*ig.face_incident_edges_dir_)[f0][eid0];
			Dart edge_dart_0 = value<std::vector<Dart>>(ig, igAttribs.face_edge_dart, f0)[eid0];

			IncidenceGraph::Face f1 = (*ig.edge_incident_faces_)[be_next.index_][0];
			uint32 eid1 = get_incident_edge_id(ig, f1, be_next);
			// uint8 edge_dir_1 = (*ig.face_incident_edges_dir_)[f1][eid1];
			Dart edge_dart_1 = value<std::vector<Dart>>(ig, igAttribs.face_edge_dart, f1)[eid1];

			std::vector<Dart> D0 = {phi<2, 3, -1, 2>(m3, edge_dart_0), phi<2, 3, 1, 2, 1>(m3, edge_dart_1),
									phi<3, 2, 3, -1, 2>(m3, edge_dart_1), phi<3, 2, 3, 1, 2, 1>(m3, edge_dart_0)};

			value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F0[0])) = phi1(m3, D0[0]);
			value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F0[1])) = phi1(m3, D0[1]);
			value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F0[2])) = phi1(m3, D0[2]);
			value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(F0[3])) = phi1(m3, D0[3]);
		}

		be = be_next;
		be_inc_verts = (*ig.edge_incident_vertices_)[be.index_];
		bv = bv == be_inc_verts.first ? be_inc_verts.second : be_inc_verts.first;
		inc_leaflet_edges = get_incident_leaflet_edges(ig, bv, be);
		for (IncidenceGraph::Edge le : inc_leaflet_edges)
		{
			if (le != be && degree(ig, le) == 1)
				be_next = le;
		}
	} while (be != boundary_edges[0]);

	return true;
}

bool sew_sections_igh(CMap2& m2, IG_M2Attributes& m2Attribs, CMap3& m3)
{
	parallel_foreach_cell(m2, [&](CMap2::Edge e) -> bool {
		if (is_incident_to_boundary(m2, e))
			return true;

		std::vector<CMap2::HalfEdge> halfedges = incident_halfedges(m2, e);
		Dart d1 = value<Dart>(m2, m2Attribs.halfedge_volume_connection, halfedges[0]);
		Dart d2 = value<Dart>(m2, m2Attribs.halfedge_volume_connection, halfedges[1]);

		sew_volumes_igh(m3, d1, phi1(m3, d2));
		return true;
	});

	std::cout << "holes " << close(m3, false) << std::endl;
	return true;
}

bool set_volumes_geometry_igh(IncidenceGraph& ig, IG_GAttributes& igAttribs, IncidenceGraphData& incidenceGraph_data,
							  CMap2& m2, IG_M2Attributes& m2Attribs, CMap3& m3, IG_M3Attributes& m3Attribs)
{
	parallel_foreach_cell(m2, [&](CMap2::Volume w) -> bool {
		Dart d = value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(w.dart));
		Dart m3d = phi_1(m3, d);
		value<Vec3>(m3, m3Attribs.vertex_position, CMap3::Vertex(m3d)) = value<Vec3>(m2, m2Attribs.volume_center, w);
		return true;
	});

	parallel_foreach_cell(m2, [&](CMap2::Edge e) -> bool {
		Dart d = value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(e.dart));
		Dart m3d = phi1(m3, d);
		value<Vec3>(m3, m3Attribs.vertex_position, CMap3::Vertex(m3d)) = value<Vec3>(m2, m2Attribs.edge_mid, e);
		return true;
	});

	CellMarker<CMap2, CMap2::Vertex> vertex_marker_m2(m2);
	for (Dart m2d = m2.begin(), end = m2.end(); m2d != end; m2d = m2.next(m2d))
	{
		if (!is_boundary(m2, m2d))
		{
			CMap2::Vertex m2v(phi1(m2, m2d));
			if (!vertex_marker_m2.is_marked(m2v))
			{
				vertex_marker_m2.mark(m2v);
				Dart m3d = value<Dart>(m2, m2Attribs.halfedge_volume_connection, CMap2::HalfEdge(m2d));
				value<Vec3>(m3, m3Attribs.vertex_position, CMap3::Vertex(m3d)) =
					value<Vec3>(m2, m2Attribs.vertex_position, m2v);
			}
		}
	}

	CellMarker<IncidenceGraph, IncidenceGraph::Vertex> vertex_marker_ig(ig);

	for (auto& leaflet : incidenceGraph_data.leaflets)
	{
		std::vector<IncidenceGraph::Vertex> inside_vertices;
		IncidenceGraph::Edge boundary_edge;
		IncidenceGraph::Edge boundary_edge_dir;

		// get leaflet inside vertices & one boundary edge
		for (IncidenceGraph::Face f : leaflet)
		{
			const std::vector<IncidenceGraph::Edge>& inc_edges = (*ig.face_incident_edges_)[f.index_];
			const std::vector<uint8>& inc_edges_dir = (*ig.face_incident_edges_dir_)[f.index_];

			for (uint32 i = 0; i < inc_edges.size(); ++i)
			{
				IncidenceGraph::Edge e = inc_edges[i];
				uint32 d = degree(ig, e);
				if (d == 1)
				{
					boundary_edge = e;
					boundary_edge_dir = inc_edges_dir[i];
				}
				const std::pair<IncidenceGraph::Vertex, IncidenceGraph::Vertex>& inc_verts =
					(*ig.edge_incident_vertices_)[e.index_];
				if (!vertex_marker_ig.is_marked(inc_verts.first))
				{
					vertex_marker_ig.mark(inc_verts.first);
					std::pair<uint32, uint32> info = pseudo_degree(ig, inc_verts.first); // { nb_leaflets, nb_edges }
					if (info.first == 0 && info.second == 0)
						inside_vertices.push_back(inc_verts.first);
				}
				if (!vertex_marker_ig.is_marked(inc_verts.second))
				{
					vertex_marker_ig.mark(inc_verts.second);
					std::pair<uint32, uint32> info = pseudo_degree(ig, inc_verts.second); // { nb_leaflets, nb_edges }
					if (info.first == 0 && info.second == 0)
						inside_vertices.push_back(inc_verts.second);
				}
			}
		}

		std::cout << "leaflet has " << inside_vertices.size() << " inside vertices" << std::endl;

		for (IncidenceGraph::Vertex iv : inside_vertices)
		{
			const Vec3& ivpos = value<Vec3>(ig, igAttribs.vertex_position, iv);
			Scalar radius = value<Scalar>(ig, igAttribs.vertex_radius, iv);

			Vec3 v_normal{0, 0, 0};
			foreach_incident_face(ig, iv, [&](IncidenceGraph::Face f) -> bool {
				const Vec3& n = value<Vec3>(ig, igAttribs.face_normal, f);
				v_normal += n;
				return true;
			});
			v_normal.normalize();

			std::cout << " pos -> " << ivpos[0] << " " << ivpos[1] << " " << ivpos[2] << std::endl;
			std::cout << " rad -> " << radius << std::endl;
			std::cout << " nor -> " << v_normal[0] << " " << v_normal[1] << " " << v_normal[2] << std::endl;

			IncidenceGraph::Edge e = (*ig.vertex_incident_edges_)[iv.index_][0];
			IncidenceGraph::Face f = (*ig.edge_incident_faces_)[e.index_][0];
			uint32 eid = get_incident_edge_id(ig, f, e);
			Dart edge_dart = value<std::vector<Dart>>(ig, igAttribs.face_edge_dart, f)[eid];
			uint8 edge_dir = (*ig.face_incident_edges_dir_)[f.index_][eid];
			edge_dart = edge_dir == 0 ? edge_dart : phi1(m3, edge_dart);

			value<Vec3>(m3, m3Attribs.vertex_position, CMap3::Vertex(edge_dart)) = ivpos;
			value<Vec3>(m3, m3Attribs.vertex_position, CMap3::Vertex(phi<2, 1, 1>(m3, edge_dart))) =
				ivpos + radius * v_normal;
			value<Vec3>(m3, m3Attribs.vertex_position, CMap3::Vertex(phi<3, 2, -1>(m3, edge_dart))) =
				ivpos - radius * v_normal;
		}

		// set boundary geometry
		// start from one boundary edge
		IncidenceGraph::Edge be = boundary_edge;
		std::pair<IncidenceGraph::Vertex, IncidenceGraph::Vertex> be_inc_verts =
			(*ig.edge_incident_vertices_)[be.index_];
		IncidenceGraph::Vertex bv = boundary_edge_dir == 0 ? be_inc_verts.second : be_inc_verts.first;
		std::vector<IncidenceGraph::Edge> inc_leaflet_edges = get_incident_leaflet_edges(ig, bv, be);
		IncidenceGraph::Edge be_next;
		for (IncidenceGraph::Edge le : inc_leaflet_edges)
		{
			if (le != be && degree(ig, le) == 1)
				be_next = le;
		}
		do
		{
			const Vec3& bvpos = value<Vec3>(ig, igAttribs.vertex_position, bv);
			Scalar radius = value<Scalar>(ig, igAttribs.vertex_radius, bv);

			std::pair<uint32, uint32> info = pseudo_degree(ig, bv); // { nb_leaflets, nb_edges }

			Vec3 leaflet_dir{0, 0, 0};
			Vec3 leaflet_normal{0, 0, 0};
			for (IncidenceGraph::Edge le : inc_leaflet_edges)
			{
				const std::pair<IncidenceGraph::Vertex, IncidenceGraph::Vertex>& inc_verts =
					(*ig.edge_incident_vertices_)[le.index_];
				IncidenceGraph::Vertex iv = inc_verts.first == bv ? inc_verts.second : inc_verts.first;
				leaflet_dir += (value<Vec3>(ig, igAttribs.vertex_position, iv) - bvpos).normalized();
				const std::vector<IncidenceGraph::Face>& inc_faces = (*ig.edge_incident_faces_)[le.index_];
				for (IncidenceGraph::Face iface : inc_faces)
				{
					const Vec3& n = value<Vec3>(ig, igAttribs.face_normal, iface);
					leaflet_normal += n;
				}
			}
			leaflet_dir.normalize();
			leaflet_normal.normalize();

			if (info.first == 1 && info.second == 0) // corners (1,0)
			{
				IncidenceGraph::Face f = (*ig.edge_incident_faces_)[be.index_][0];
				uint32 eid = get_incident_edge_id(ig, f, be);
				// uint8 edge_dir = (*ig.face_incident_edges_dir_)[f][eid];
				Dart edge_dart = value<std::vector<Dart>>(ig, igAttribs.face_edge_dart, f)[eid];

				value<Vec3>(m3, m3Attribs.vertex_position, CMap3::Vertex(edge_dart)) = bvpos;
				value<Vec3>(m3, m3Attribs.vertex_position, CMap3::Vertex(phi<2, 1, 1>(m3, edge_dart))) =
					bvpos + radius * leaflet_normal;
				value<Vec3>(m3, m3Attribs.vertex_position, CMap3::Vertex(phi<3, 2, -1>(m3, edge_dart))) =
					bvpos - radius * leaflet_normal;

				value<Vec3>(m3, m3Attribs.vertex_position, CMap3::Vertex(phi<2, 3, 2, 1, 1>(m3, edge_dart))) =
					bvpos - radius * leaflet_dir;
				value<Vec3>(m3, m3Attribs.vertex_position, CMap3::Vertex(phi<2, 3, 2, 1, 2, -1>(m3, edge_dart))) =
					bvpos - radius * leaflet_dir + radius * leaflet_normal;
				value<Vec3>(m3, m3Attribs.vertex_position, CMap3::Vertex(phi<2, 3, 2, 1, 3, 2, 1, 1>(m3, edge_dart))) =
					bvpos - radius * leaflet_dir - radius * leaflet_normal;
			}

			be = be_next;
			be_inc_verts = (*ig.edge_incident_vertices_)[be.index_];
			bv = bv == be_inc_verts.first ? be_inc_verts.second : be_inc_verts.first;
			inc_leaflet_edges = get_incident_leaflet_edges(ig, bv, be);
			for (IncidenceGraph::Edge le : inc_leaflet_edges)
			{
				if (le != be && degree(ig, le) == 1)
					be_next = le;
			}
		} while (be != boundary_edge);
	}

	return true;
}

} // namespace modeling

} // namespace cgogn
