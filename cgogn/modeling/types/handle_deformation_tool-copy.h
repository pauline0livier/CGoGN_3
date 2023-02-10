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

#ifndef CGOGN_MODELING_HANDLE_DEFORMATION_TOOL_H_
#define CGOGN_MODELING_HANDLE_DEFORMATION_TOOL_H_

#include <cgogn/modeling/types/space_deformation_tool.h>

namespace cgogn
{

namespace modeling
{

template <typename MESH>
class HandleDeformationTool : public SpaceDeformationTool<MESH>
{
	using Graph = cgogn::IncidenceGraph;

public:
	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;
	using MeshVertex = typename mesh_traits<MESH>::Vertex;
	using MeshFace = typename mesh_traits<MESH>::Face;

	using Vec2 = geometry::Vec2;
	using Vec3 = geometry::Vec3;
	using Scalar = geometry::Scalar;

	Graph* control_handle_;

	std::shared_ptr<Graph::Attribute<Vec3>> control_handle_vertex_position_;
	std::shared_ptr<boost::synapse::connection> handle_attribute_update_connection_;

	HandleDeformationTool() : SpaceDeformationTool<MESH>(), control_handle_vertex_position_(nullptr)
	{
	}

	~HandleDeformationTool()
	{
	}

	void create_space_tool(Graph* g, Graph::Attribute<Vec3>* vertex_position, Graph::Attribute<Scalar>* vertex_radius,
						   const Vec3& center1, const Vec3& center2)
	{
		control_handle_ = g;
		handle_vertex_ = cgogn::modeling::create_handle(*g, vertex_position, vertex_radius, center1, center2);

		control_handle_vertex_position_ = cgogn::get_attribute<Vec3, Graph::Vertex>(*g, "position");

		/*std::shared_ptr<Graph::Attribute<uint32>> vertex_index =
			cgogn::add_attribute<uint32, Graph::Vertex>(*control_handle_, "vertex_index");*/
		// cgogn::modeling::set_graph_attribute_vertex_index(*control_handle_, vertex_index.get());
	}

	void set_influence_cage_handle(MESH* m, CMap2::Attribute<Vec3>* vertex_position,
								   CMap2::Attribute<Vec3>* local_vertex_position, const Vec3& bb_min,
								   const Vec3& bb_max, const Vec3& handle_position, Vec3& normal)
	{
		this->influence_cage_ = m;
		handle_normal_ = normal;
		handle_position_ = handle_position;

		// plane ax+by+cz+d = 0, with (a,b,c) = handle_normal and (x,y,z) = bb_min

		// bbmin belongs to the plane
		float d_handle_min = -(handle_normal_.dot(bb_min));

		// find alpha such that center + alpha*handle_normal belongs to plane
		float alpha_min = -d_handle_min - (handle_normal_.dot(handle_position));

		Eigen::Vector3d center_min_plane = handle_position + alpha_min * handle_normal_;

		Vec3 u = bb_min - center_min_plane;
		u.normalize();

		Vec3 v = handle_normal_.cross(u);
		v.normalize();

		local_frame_.row(0) = u;
		local_frame_.row(1) = v;
		local_frame_.row(2) = handle_normal_;
		local_frame_inverse_ = local_frame_.inverse();

		const Vec3 local_bb_min = local_frame_ * (bb_min - handle_position_);
		const Vec3 local_bb_max = local_frame_ * (bb_max - handle_position_);

		const double radius = local_bb_min.norm();

		if (local_bb_min[2] > local_bb_max[2])
		{
			cgogn::modeling::create_handle_box(*m, vertex_position, local_vertex_position, handle_position, radius,
											   local_frame_inverse_, local_bb_min[2], local_bb_max[2]);
		}
		else
		{
			cgogn::modeling::create_handle_box(*m, vertex_position, local_vertex_position, handle_position, radius,
											   local_frame_inverse_, local_bb_max[2], local_bb_min[2]);
		}

		this->influence_cage_vertex_position_ = cgogn::get_attribute<Vec3, MeshVertex>(*m, "position");

		std::shared_ptr<Attribute<uint32>> vertex_index =
			cgogn::add_attribute<uint32, MeshVertex>(*this->influence_cage_, "vertex_index");
		cgogn::modeling::set_attribute_vertex_index(*this->influence_cage_, vertex_index.get());

		std::shared_ptr<Attribute<bool>> marked_vertices =
			cgogn::add_attribute<bool, MeshVertex>(*this->influence_cage_, "marked_vertices");
		cgogn::modeling::set_attribute_marked_vertices(*this->influence_cage_, marked_vertices.get());
	}

	void update_influence_cage_position()
	{
		handle_position_ = value<Vec3>(*control_handle_, control_handle_vertex_position_, handle_vertex_);

		std::shared_ptr<Attribute<Vec3>> local_position =
			cgogn::get_attribute<Vec3, MeshVertex>(*(this->influence_cage_), "local_position");

		foreach_cell(*(this->influence_cage_), [&](MeshVertex v) -> bool {
			Vec3 local_point = value<Vec3>(*(this->influence_cage_), local_position, v);

			value<Vec3>(*(this->influence_cage_), this->influence_cage_vertex_position_, v) =
				local_frame_inverse_ * local_point + handle_position_;

			return true;
		});
	}

	void set_up_attenuation(MESH& object, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{

		uint32 nbv_object = nb_cells<MeshVertex>(object);

		this->attenuation_.resize(nbv_object);
		this->attenuation_.setZero();

		compute_attenuation_cage(object);
	}

private:
	Graph::Vertex handle_vertex_;
	Vec3 handle_position_; // also frame_origin
	Vec3 handle_normal_;
	Eigen::Matrix3d local_frame_;
	Eigen::Matrix3d local_frame_inverse_;

	void compute_attenuation_cage(MESH& object)
	{
		std::shared_ptr<Attribute<uint32>> cage_face_indices =
			add_attribute<uint32, MeshFace>(*(this->influence_cage_), "face_indices");

		cgogn::modeling::set_attribute_face_index(*(this->influence_cage_), cage_face_indices.get());

		std::shared_ptr<Attribute<uint32>> object_vertex_index =
			get_attribute<uint32, MeshVertex>(object, "vertex_index");

		std::shared_ptr<Attribute<Vec3>> object_position = get_attribute<Vec3, MeshVertex>(object, "position");

		std::shared_ptr<Attribute<Scalar>> cage_face_indices =
			add_attribute<Scalar, MeshVertex>(m, "vertex_geodesic_distance");

		uint32 nbf_cage = 2 * nb_cells<MeshFace>(*(this->influence_cage_));
		uint32 nbv_cage = nb_cells<MeshVertex>(*(this->influence_cage_));

		float h = 0.0f;
		float max_dist = 0.0f;
		std::vector<Vec2> attenuation_points;
		this->influence_area_->foreach_cell([&](MeshVertex v) {
			uint32 surface_point_idx = value<uint32>(object, object_vertex_index, v);

			Vec3 surface_point = value<Vec3>(object, object_position, v);

			Vec3 point_to_handle = (surface_point - handle_position_);
			// double i_dist = point_to_handle.squaredNorm();
			double i_dist = point_to_handle.norm();
			// this->cage_influence_distance(surface_point_idx, nbf_cage, nbv_cage);

			// std::cout << i_dist << std::endl;
			if (i_dist > max_dist)
			{
				max_dist = i_dist;
			}

			// this->attenuation_(surface_point_idx) = (float)sin(0.5*M_PI * (i_dist ));
			this->attenuation_(surface_point_idx) = i_dist;
		});

		this->influence_area_->foreach_cell([&](MeshVertex v) {
			uint32 surface_point_idx = value<uint32>(object, object_vertex_index, v);

			this->attenuation_(surface_point_idx) = 1.0 - (this->attenuation_(surface_point_idx) / max_dist);
		});
	}

	void geodesic_distance(MESH& m, const Attribute<Vec3>* vertex_position,
						    Attribute<Scalar>* vertex_geodesic_distance)
	{

		std::shared_ptr<Attribute<uint32>> vertex_index = cgogn::get_attribute<uint32, Vertex>(m, "vertex_index");

		uint32 nb_vertices = nb_cells<Vertex>(m);

		Eigen::SparseMatrix<Scalar, Eigen::ColMajor> Lc =
			cgogn::geometry::cotan_operator_matrix(m, vertex_index.get(), vertex_position);

		auto vertex_area = add_attribute<Scalar, Vertex>(m, "__vertex_area");
		cgogn::geometry::compute_area<Vertex>(m, vertex_position, vertex_area.get());

		Eigen::VectorXd A(nb_vertices);
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			uint32 vidx = value<uint32>(m, vertex_index, v);
			A(vidx) = value<Scalar>(m, vertex_area, v);
			return true;
		});

		Eigen::VectorXd u0(nb_vertices);
		u0.setZero();
		/*source_vertices->foreach_cell([&](Vertex v) {
			uint32 vidx = value<uint32>(m, vertex_index, v);
			u0(vidx) = 1.0;
		});*/

		Scalar h = cgogn::geometry::mean_edge_length(m, vertex_position);
		Scalar t = h * h;

		Eigen::SparseMatrix<Scalar, Eigen::ColMajor> Am(A.asDiagonal());
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar, Eigen::ColMajor>> heat_solver(Am - t * Lc);
		Eigen::VectorXd u = heat_solver.solve(u0);

		auto vertex_heat = get_or_add_attribute<Scalar, Vertex>(m, "__vertex_heat");
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			uint32 vidx = value<uint32>(m, vertex_index, v);
			value<Scalar>(m, vertex_heat, v) = u(vidx);
			return true;
		});

		auto face_heat_gradient = get_or_add_attribute<Vec3, Face>(m, "__face_heat_gradient");
		parallel_foreach_cell(m, [&](Face f) -> bool {
			Vec3 g(0, 0, 0);
			Vec3 n = cgogn::geometry::normal(m, f, vertex_position);
			Scalar a = cgogn::geometry::area(m, f, vertex_position);
			std::vector<Vertex> vertices = incident_vertices(m, f);
			for (uint32 i = 0; i < vertices.size(); ++i)
			{
				Vec3 e = value<Vec3>(m, vertex_position, vertices[(i + 2) % vertices.size()]) -
						 value<Vec3>(m, vertex_position, vertices[(i + 1) % vertices.size()]);
				g += value<Scalar>(m, vertex_heat, vertices[i]) * n.cross(e);
			}
			g /= 2 * a;
			value<Vec3>(m, face_heat_gradient, f) = -1.0 * g.normalized();
			return true;
		});

		auto vertex_heat_gradient_div = get_or_add_attribute<Scalar, Vertex>(m, "__vertex_heat_gradient_div");
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			Scalar d = vertex_gradient_divergence(m, v, face_heat_gradient.get(), vertex_position);
			value<Scalar>(m, vertex_heat_gradient_div, v) = d;
			return true;
		});

		Eigen::VectorXd b(nb_vertices);
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			uint32 vidx = value<uint32>(m, vertex_index, v);
			b(vidx) = value<Scalar>(m, vertex_heat_gradient_div, v);
			return true;
		});

		Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar, Eigen::ColMajor>> poisson_solver(Lc);
		Eigen::VectorXd dist = poisson_solver.solve(b);

		Scalar min = dist.minCoeff();
		parallel_foreach_cell(m, [&](Vertex v) -> bool {
			uint32 vidx = value<uint32>(m, vertex_index, v);
			value<Scalar>(m, vertex_geodesic_distance, v) = dist(vidx) - min;
			return true;
		});

		//remove_attribute<Vertex>(m, vertex_index);
		remove_attribute<Vertex>(m, vertex_area);
		// remove_attribute<Vertex>(m, vertex_heat);
		// remove_attribute<Face>(m, face_heat_gradient);
		// remove_attribute<Vertex>(m, vertex_heat_gradient_div);

		mesh_provider_->emit_attribute_changed(m, vertex_geodesic_distance);
	}
};

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_HANDLE_DEFORMATION_TOOL_H_
