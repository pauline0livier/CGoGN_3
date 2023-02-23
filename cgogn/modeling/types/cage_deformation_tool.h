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

#ifndef CGOGN_MODELING_CAGE_DEFORMATION_TOOL_H_
#define CGOGN_MODELING_CAGE_DEFORMATION_TOOL_H_

#include <cgogn/modeling/types/space_deformation_tool.h>

namespace cgogn
{

namespace modeling
{

template <typename MESH>

/**
 * @Class CageDeformationTool
 * Represents the local cage deformation tool
 */
class CageDeformationTool : public SpaceDeformationTool<MESH>
{
	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;

	using Vec2 = geometry::Vec2;
	using Vec3 = geometry::Vec3;

	struct Triangle
	{
		std::vector<CMap2::Vertex> vertices;
		std::vector<Vec3> positions;
		std::vector<uint32_t> indices;
		std::vector<uint32_t> virtual_indices;
		Vec3 normal;
		std::pair<Vec3, Vec3> edges;
	};

	struct Local_direction_control_planes
	{
		double d_min;
		double d_max;
		double d_gap;

		std::pair<Triangle, Triangle> triangles_d_min;
		std::pair<Triangle, Triangle> triangles_d_max;

		Vec3 direction;
		Vec3 shift_after_d_max;
		Vec3 shift_before_d_min;
	};

public:
	MESH* control_cage_;
	std::shared_ptr<Attribute<Vec3>> control_cage_vertex_position_;
	std::shared_ptr<Attribute<Vec3>> influence_cage_vertex_position_;

	Eigen::VectorXd attenuation_;

	std::shared_ptr<boost::synapse::connection> cage_attribute_update_connection_;

	std::vector<Triangle> cage_triangles_;

	MESH* influence_cage_;

	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> control_cage_coords_;
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> control_cage_normal_coords_;

	CageDeformationTool() : SpaceDeformationTool<MESH>(), m_hFactor(-1.0f), control_cage_vertex_position_(nullptr)
	{
	}

	~CageDeformationTool()
	{
	}

	void create_space_tool(MESH* m, CMap2::Attribute<Vec3>* vertex_position, const Vec3& bb_min, const Vec3& bb_max,
						   const Vec3& center, const Vec3& normal)
	{
		control_cage_ = m;

		/*cgogn::modeling::create_cage_box(*m, vertex_position, bb_min, bb_max, center, normal);*/ // not working well
																								   // so far
		cgogn::modeling::create_bounding_box(*m, vertex_position, bb_min, bb_max);

		control_cage_vertex_position_ = cgogn::get_attribute<Vec3, Vertex>(*m, "position");

		control_cage_bb_min_ = bb_min;
		control_cage_bb_max_ = bb_max;

		control_cage_vertex_index_ = cgogn::add_attribute<uint32, Vertex>(*control_cage_, "vertex_index");
		cgogn::modeling::set_attribute_vertex_index(*control_cage_, control_cage_vertex_index_.get());

		init_triangles();

		init_control_cage_plane();
	}

	void set_influence_cage(MESH& object, const CMap2::Attribute<Vec3>* object_vertex_position, MESH* m,
							CMap2::Attribute<Vec3>* vertex_position)
	{
		influence_cage_ = m;

		/*std::pair<Vec3, Vec3> influence_area_borders = cgogn::modeling::get_border_values_in_set(object,
		 * object_vertex_position, this->influence_area_); */

		// first simplification 3x control cage
		std::tuple<Vec3, Vec3, Vec3> res_extended_bounding_box =
			cgogn::modeling::get_extended_bounding_box(control_cage_bb_min_, control_cage_bb_max_, 3.0);

		influence_cage_bb_min_ = std::get<0>(res_extended_bounding_box);
		influence_cage_bb_max_ = std::get<1>(res_extended_bounding_box);

		cgogn::modeling::create_bounding_box(*m, vertex_position, influence_cage_bb_min_, influence_cage_bb_max_);

		influence_cage_vertex_position_ = cgogn::get_attribute<Vec3, Vertex>(*m, "position");

		std::shared_ptr<Attribute<bool>> marked_vertices =
			cgogn::add_attribute<bool, Vertex>(*influence_cage_, "marked_vertices");
		cgogn::modeling::set_attribute_marked_vertices(*influence_cage_, marked_vertices.get());

		// update influence area
		foreach_cell(object, [&](Vertex v) -> bool {
			const Vec3& surface_point = value<Vec3>(object, object_vertex_position, v);

			bool inside_cage = check_mvc_on_point_inside_influence_cage(surface_point);

			if (inside_cage)
			{
				this->influence_area_->select(v);
			}

			return true;
		});
	}

	void bind_mvc(MESH& object, const std::shared_ptr<Attribute<Vec3>>& object_vertex_position)
	{

		std::shared_ptr<Attribute<uint32>> object_vertex_index = get_attribute<uint32, Vertex>(object, "vertex_index");

		uint32 nbv_object = nb_cells<Vertex>(object);
		uint32 nbv_cage = nb_cells<Vertex>(*control_cage_);

		control_cage_coords_.resize(nbv_object, nbv_cage);
		control_cage_coords_.setZero();

		this->influence_area_->foreach_cell([&](Vertex v) {

			const Vec3& surface_point = value<Vec3>(object, object_vertex_position, v);

			uint32 surface_point_index = value<uint32>(object, object_vertex_index, v);

			const double d_x = -(surface_point.dot(local_x_direction_control_planes_.direction)),

						 d_y = -(surface_point.dot(local_y_direction_control_planes_.direction)),

						 d_z = -(surface_point.dot(local_z_direction_control_planes_.direction));

			const bool valid_x_dir = (d_x <= local_x_direction_control_planes_.d_max &&
									  d_x >= local_x_direction_control_planes_.d_min),

					   valid_y_dir = (d_y <= local_y_direction_control_planes_.d_max &&
									  d_y >= local_y_direction_control_planes_.d_min),

					   valid_z_dir = (d_z <= local_z_direction_control_planes_.d_max &&
									  d_z >= local_z_direction_control_planes_.d_min);

			std::cout << "before big loop" << std::endl; 

			if (valid_x_dir && valid_y_dir && valid_z_dir)
			{
				std::cout << "inside cage " << std::endl; 
				compute_mvc_on_point_inside_cage(surface_point, surface_point_index);
			}
			else
			{
				std::vector<Triangle> virtual_cube_triangles;
				if (valid_x_dir)
				{
					if (valid_y_dir)
					{
						if (d_z > local_z_direction_control_planes_.d_max)
						{
							// fakeCube[5]
							std::cout << "fake cube 5" << std::endl; 
							// Beware, need to inverse the normal of the triangle1 and triangle2

							virtual_cube_triangles =
								get_virtual_cube_triangles(local_z_direction_control_planes_.triangles_d_max,
														   local_z_direction_control_planes_.shift_after_d_max);
						}
						else
						{
							// fakeCube[4]
							std::cout << "fake cube 4" << std::endl; 
							virtual_cube_triangles =
								get_virtual_cube_triangles(local_z_direction_control_planes_.triangles_d_min,
														   local_z_direction_control_planes_.shift_before_d_min);
						}
					}
					else if (valid_z_dir)
					{
						if (d_y > local_y_direction_control_planes_.d_max)
						{
							// fakeCube[3]
							std::cout << "fake cube 3" << std::endl; 
							virtual_cube_triangles =
								get_virtual_cube_triangles(local_y_direction_control_planes_.triangles_d_max,
														   local_y_direction_control_planes_.shift_after_d_max);
						}
						else
						{
							// fakeCube[2]
							std::cout << "fake cube 2" << std::endl; 
							virtual_cube_triangles =
								get_virtual_cube_triangles(local_y_direction_control_planes_.triangles_d_min,
														   local_y_direction_control_planes_.shift_before_d_min);
						}
					}
					else
					{
						if (d_y > local_y_direction_control_planes_.d_max)
						{

							if (d_z > local_z_direction_control_planes_.d_max)
							{
								// fakeCubeDiag[11]
								std::cout << "fake diag cube 11" << std::endl; 
								// find intersection edge between y and z triangle
								std::vector<Vec3> intersect_positions =
									find_intersection_positions_face(local_y_direction_control_planes_.triangles_d_max,
																	local_z_direction_control_planes_.triangles_d_max);

								std::pair<Triangle, Triangle> new_face = get_face_from_intersecting_edge(
									intersect_positions, local_y_direction_control_planes_.shift_after_d_max);

								virtual_cube_triangles = get_virtual_cube_triangles(
									new_face, local_z_direction_control_planes_.shift_after_d_max);
							}
							else
							{
								// fakeCubeDiag[9]
								std::cout << "fake diag cube 9" << std::endl; 
								std::vector<Vec3> intersect_positions =
									find_intersection_positions_face(local_y_direction_control_planes_.triangles_d_max,
																	local_z_direction_control_planes_.triangles_d_min);

								std::pair<Triangle, Triangle> new_face = get_face_from_intersecting_edge(
									intersect_positions, local_y_direction_control_planes_.shift_after_d_max);

								virtual_cube_triangles = get_virtual_cube_triangles(
									new_face, local_z_direction_control_planes_.shift_before_d_min);
							}
						}
						else
						{

							if (d_z > local_z_direction_control_planes_.d_max)
							{
								// fakeCubeDiag[10]
								std::cout << "fake diag cube 10" << std::endl; 
								std::vector<Vec3> intersect_positions =
									find_intersection_positions_face(local_y_direction_control_planes_.triangles_d_min,
																	local_z_direction_control_planes_.triangles_d_max);

								std::pair<Triangle, Triangle> new_face = get_face_from_intersecting_edge(
									intersect_positions, local_y_direction_control_planes_.shift_before_d_min);

								virtual_cube_triangles = get_virtual_cube_triangles(
									new_face, local_z_direction_control_planes_.shift_after_d_max);
							}
							else
							{
								// fakeCubeDiag[8]
								std::cout << "fake diag cube 8" << std::endl; 
								std::vector<Vec3> intersect_positions =
									find_intersection_positions_face(local_y_direction_control_planes_.triangles_d_min,
																	local_z_direction_control_planes_.triangles_d_min);

								std::pair<Triangle, Triangle> new_face = get_face_from_intersecting_edge(
									intersect_positions, local_y_direction_control_planes_.shift_before_d_min);

								virtual_cube_triangles = get_virtual_cube_triangles(
									new_face, local_z_direction_control_planes_.shift_before_d_min);
							}
						}
					}
				}
				else if (valid_y_dir)
				{
					if (valid_z_dir)
					{
						if (d_x > local_x_direction_control_planes_.d_max)
						{
							// fakeCube[1]
							std::cout << "fake cube 1" << std::endl; 
							virtual_cube_triangles =
								get_virtual_cube_triangles(local_x_direction_control_planes_.triangles_d_max,
														   local_x_direction_control_planes_.shift_after_d_max);
						}
						else
						{
							// fakeCube[0]
							std::cout << "fake cube 0" << std::endl; 
							virtual_cube_triangles =
								get_virtual_cube_triangles(local_x_direction_control_planes_.triangles_d_min,
														   local_x_direction_control_planes_.shift_before_d_min);
						}
					}
					else
					{
						if (d_x > local_x_direction_control_planes_.d_max)
						{

							if (d_z > local_z_direction_control_planes_.d_max)
							{
								// fakeCubeDiag[3]
								std::cout << "fake diag cube 3" << std::endl; 
								std::vector<Vec3> intersect_positions =
									find_intersection_positions_face(local_x_direction_control_planes_.triangles_d_max,
																	local_z_direction_control_planes_.triangles_d_max);

								std::pair<Triangle, Triangle> new_face = get_face_from_intersecting_edge(
									intersect_positions, local_x_direction_control_planes_.shift_after_d_max);

								virtual_cube_triangles = get_virtual_cube_triangles(
									new_face, local_z_direction_control_planes_.shift_after_d_max);
							}
							else
							{
								// fakeCubeDiag[1]
								std::cout << "fake diag cube 1" << std::endl; 
								std::vector<Vec3> intersect_positions =
									find_intersection_positions_face(local_x_direction_control_planes_.triangles_d_max,
																	local_z_direction_control_planes_.triangles_d_min);

								std::pair<Triangle, Triangle> new_face = get_face_from_intersecting_edge(
									intersect_positions, local_x_direction_control_planes_.shift_after_d_max);

								virtual_cube_triangles = get_virtual_cube_triangles(
									new_face, local_z_direction_control_planes_.shift_before_d_min);
							}
						}
						else
						{

							if (d_z > local_z_direction_control_planes_.d_max)
							{
								// fakeCubeDiag[2]
								std::cout << "fake diag cube 2" << std::endl; 
								std::vector<Vec3> intersect_positions =
									find_intersection_positions_face(local_x_direction_control_planes_.triangles_d_min,
																	local_z_direction_control_planes_.triangles_d_max);

								std::pair<Triangle, Triangle> new_face = get_face_from_intersecting_edge(
									intersect_positions, local_x_direction_control_planes_.shift_before_d_min);

								std::vector<Triangle> virtual_cube_triangles = get_virtual_cube_triangles(
									new_face, local_z_direction_control_planes_.shift_after_d_max);
							}
							else
							{
								// fakeCubeDiag[0]
								std::cout << "fake diag cube 0" << std::endl; 
								std::vector<Vec3> intersect_positions =
									find_intersection_positions_face(local_x_direction_control_planes_.triangles_d_min,
																	local_z_direction_control_planes_.triangles_d_min);

								std::pair<Triangle, Triangle> new_face = get_face_from_intersecting_edge(
									intersect_positions, local_x_direction_control_planes_.shift_before_d_min);

								virtual_cube_triangles = get_virtual_cube_triangles(
									new_face, local_z_direction_control_planes_.shift_before_d_min);
							}
						}
					}
				}
				else if (valid_z_dir)
				{
					if (d_y > local_y_direction_control_planes_.d_max)
					{
						if (d_x > local_x_direction_control_planes_.d_max)
						{
							// fakeCubeDiag[7]
							std::cout << "fake diag cube 7" << std::endl; 
							std::vector<Vec3> intersect_positions =
								find_intersection_positions_face(local_y_direction_control_planes_.triangles_d_max,
																local_x_direction_control_planes_.triangles_d_max);

							std::pair<Triangle, Triangle> new_face = get_face_from_intersecting_edge(
								intersect_positions, local_y_direction_control_planes_.shift_after_d_max);

							virtual_cube_triangles = get_virtual_cube_triangles(
								new_face, local_x_direction_control_planes_.shift_after_d_max);
						}
						else
						{
							// fakeCubeDiag[5]
							std::cout << "fake diag cube 5" << std::endl; 
							std::vector<Vec3> intersect_positions =
								find_intersection_positions_face(local_y_direction_control_planes_.triangles_d_max,
																local_x_direction_control_planes_.triangles_d_min);

							std::pair<Triangle, Triangle> new_face = get_face_from_intersecting_edge(
								intersect_positions, local_y_direction_control_planes_.shift_after_d_max);

							virtual_cube_triangles = get_virtual_cube_triangles(
								new_face, local_x_direction_control_planes_.shift_before_d_min);
						}
					}
					else
					{

						if (d_x > local_x_direction_control_planes_.d_max)
						{
							// fakeCubeDiag[6]
							std::cout << "fake diag cube 6" << std::endl; 
							std::vector<Vec3> intersect_positions =
								find_intersection_positions_face(local_y_direction_control_planes_.triangles_d_min,
																local_x_direction_control_planes_.triangles_d_max);

							std::pair<Triangle, Triangle> new_face = get_face_from_intersecting_edge(
								intersect_positions, local_y_direction_control_planes_.shift_before_d_min);

							virtual_cube_triangles = get_virtual_cube_triangles(
								new_face, local_x_direction_control_planes_.shift_after_d_max);
						}
						else
						{
							// fakeCubeDiag[4]
							std::cout << "fake diag cube 4" << std::endl; 
							std::vector<Vec3> intersect_positions =
								find_intersection_positions_face(local_y_direction_control_planes_.triangles_d_min,
																local_x_direction_control_planes_.triangles_d_min);

							std::pair<Triangle, Triangle> new_face = get_face_from_intersecting_edge(
								intersect_positions, local_y_direction_control_planes_.shift_before_d_min);

							virtual_cube_triangles = get_virtual_cube_triangles(
								new_face, local_x_direction_control_planes_.shift_before_d_min);
						}
					}
				}
				else
				{
					if (d_x > local_x_direction_control_planes_.d_max)
					{
						if (d_y > local_y_direction_control_planes_.d_max)
						{
							if (d_z > local_z_direction_control_planes_.d_max)
							{
								// fakeCubeOut[7]
								std::cout << "fake out cube 7" << std::endl; 
								Vec3 intersection_position =
									find_intersection_position(local_x_direction_control_planes_.triangles_d_max,
															 local_y_direction_control_planes_.triangles_d_max,
															 local_z_direction_control_planes_.triangles_d_max);

								std::pair<Triangle, Triangle> new_face = get_face_from_intersecting_vertex(
									intersection_position, local_x_direction_control_planes_.shift_after_d_max,
									local_y_direction_control_planes_.shift_after_d_max);

								virtual_cube_triangles = get_virtual_cube_triangles(
									new_face, local_z_direction_control_planes_.shift_after_d_max);
							}
							else
							{
								// fakeCubeOut[5]
								std::cout << "fake out cube 5" << std::endl; 
								Vec3 intersection_position =
									find_intersection_position(local_x_direction_control_planes_.triangles_d_max,
															 local_y_direction_control_planes_.triangles_d_max,
															 local_z_direction_control_planes_.triangles_d_min);

								std::pair<Triangle, Triangle> new_face = get_face_from_intersecting_vertex(
									intersection_position, local_x_direction_control_planes_.shift_after_d_max,
									local_y_direction_control_planes_.shift_after_d_max);

								virtual_cube_triangles = get_virtual_cube_triangles(
									new_face, local_z_direction_control_planes_.shift_before_d_min);
							}
						}
						else
						{
							if (d_z > local_z_direction_control_planes_.d_max)
							{
								// fakeCubeOut[3]
								std::cout << "fake out cube 3" << std::endl; 
								Vec3 intersection_position =
									find_intersection_position(local_x_direction_control_planes_.triangles_d_max,
															 local_y_direction_control_planes_.triangles_d_min,
															 local_z_direction_control_planes_.triangles_d_max);

								std::pair<Triangle, Triangle> new_face = get_face_from_intersecting_vertex(
									intersection_position, local_x_direction_control_planes_.shift_after_d_max,
									local_y_direction_control_planes_.shift_before_d_min);

								virtual_cube_triangles = get_virtual_cube_triangles(
									new_face, local_z_direction_control_planes_.shift_after_d_max);
							}
							else
							{
								// fakeCubeOut[1]
								std::cout << "fake out cube 1" << std::endl; 
								Vec3 intersection_position =
									find_intersection_position(local_x_direction_control_planes_.triangles_d_max,
															 local_y_direction_control_planes_.triangles_d_min,
															 local_z_direction_control_planes_.triangles_d_min);

								std::pair<Triangle, Triangle> new_face = get_face_from_intersecting_vertex(
									intersection_position, local_x_direction_control_planes_.shift_after_d_max,
									local_y_direction_control_planes_.shift_before_d_min);

								virtual_cube_triangles = get_virtual_cube_triangles(
									new_face, local_z_direction_control_planes_.shift_before_d_min);
							}
						}
					}
					else
					{
						if (d_y > local_y_direction_control_planes_.d_max)
						{
							if (d_z > local_z_direction_control_planes_.d_max)
							{
								// fakeCubeOut[6]
								std::cout << "fake out cube 6" << std::endl; 
								Vec3 intersection_position =
									find_intersection_position(local_x_direction_control_planes_.triangles_d_min,
															 local_y_direction_control_planes_.triangles_d_max,
															 local_z_direction_control_planes_.triangles_d_max);

								std::pair<Triangle, Triangle> new_face = get_face_from_intersecting_vertex(
									intersection_position, local_x_direction_control_planes_.shift_before_d_min,
									local_y_direction_control_planes_.shift_after_d_max);

								virtual_cube_triangles = get_virtual_cube_triangles(
									new_face, local_z_direction_control_planes_.shift_after_d_max);
							}
							else
							{
								// fakeCubeOut[4]
								std::cout << "fake out cube 4" << std::endl; 
								Vec3 intersection_position =
									find_intersection_position(local_x_direction_control_planes_.triangles_d_min,
															 local_y_direction_control_planes_.triangles_d_max,
															 local_z_direction_control_planes_.triangles_d_min);

								std::pair<Triangle, Triangle> new_face = get_face_from_intersecting_vertex(
									intersection_position, local_x_direction_control_planes_.shift_before_d_min,
									local_y_direction_control_planes_.shift_after_d_max);

								virtual_cube_triangles = get_virtual_cube_triangles(
									new_face, local_z_direction_control_planes_.shift_before_d_min);
							}
						}
						else
						{
							if (d_z > local_z_direction_control_planes_.d_max)
							{
								// fakeCubeOut[2]
								std::cout << "fake out cube 2" << std::endl; 
								Vec3 intersection_position =
									find_intersection_position(local_x_direction_control_planes_.triangles_d_min,
															 local_y_direction_control_planes_.triangles_d_min,
															 local_z_direction_control_planes_.triangles_d_max);

								std::pair<Triangle, Triangle> new_face = get_face_from_intersecting_vertex(
									intersection_position, local_x_direction_control_planes_.shift_before_d_min,
									local_y_direction_control_planes_.shift_before_d_min);

								virtual_cube_triangles = get_virtual_cube_triangles(
									new_face, local_z_direction_control_planes_.shift_after_d_max);
							}
							else
							{
								// fakeCubeOut[0]
								std::cout << "fake out cube 0" << std::endl; 
								Vec3 intersection_position =
									find_intersection_position(local_x_direction_control_planes_.triangles_d_min,
															 local_y_direction_control_planes_.triangles_d_min,
															 local_z_direction_control_planes_.triangles_d_min);

								std::pair<Triangle, Triangle> new_face = get_face_from_intersecting_vertex(
									intersection_position, local_x_direction_control_planes_.shift_before_d_min,
									local_y_direction_control_planes_.shift_before_d_min);

								virtual_cube_triangles = get_virtual_cube_triangles(
									new_face, local_z_direction_control_planes_.shift_before_d_min);
							}
						}
					}
				}

				compute_mvc_on_point_outside_cage(surface_point, surface_point_index, virtual_cube_triangles);
			}
		});
	}

	void set_center_control_cage(Vec3& center)
	{
		control_cage_center_ = center;
	}

	void update_influence_cage_position()
	{
		foreach_cell(*control_cage_, [&](Vertex v) -> bool {
			const Vec3& cage_point = value<Vec3>(*control_cage_, control_cage_vertex_position_, v);

			value<Vec3>(*(this->influence_cage_), this->influence_cage_vertex_position_, v) =
				((cage_point - control_cage_center_) * 1.5) + control_cage_center_;

			return true;
		});
	}

	void update_deformation_object(MESH& object, const std::shared_ptr<Attribute<Vec3>>& object_vertex_position)
	{
		std::shared_ptr<Attribute<uint32>> object_vertex_indices =
			get_attribute<uint32, Vertex>(object, "vertex_index");

		this->influence_area_->foreach_cell([&](Vertex v) -> bool {
			uint32 object_vertex_index = value<uint32>(object, object_vertex_indices, v);

			Vec3 new_pos_ = {0.0, 0.0, 0.0};

			foreach_cell(*control_cage_, [&](Vertex cv) -> bool {
				const Vec3& cage_point = value<Vec3>(*control_cage_, control_cage_vertex_position_, cv);

				uint32 cage_point_index = value<uint32>(*control_cage_, control_cage_vertex_index_, cv);

				new_pos_ += control_cage_coords_(object_vertex_index, cage_point_index) * cage_point;

				return true;
			});

			if ((new_pos_[0] != 0.0) && (new_pos_[1] != 0.0) && (new_pos_[2] != 0.0))
			{
				value<Vec3>(object, object_vertex_position, v) = new_pos_;
			}
			else
			{
				std::cout << "empty " << std::endl;
			}

			return true;
		});
	}

	/*void set_up_attenuation(MESH& object, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
	{

		std::shared_ptr<Attribute<uint32>> object_vertex_index =
			get_attribute<uint32, Vertex>(object, "vertex_index");

		uint32 nbv_object = nb_cells<Vertex>(object);

		control_area_validity_.resize(nbv_object);
		control_area_validity_.setZero();

		foreach_cell(object, [&](Vertex v) -> bool {
			const Vec3& surface_point = value<Vec3>(object, vertex_position, v);

			//bool inside_cage = local_mvc_pt_control_area(surface_point);

			if (inside_cage)
			{
				uint32 surface_point_index = value<uint32>(object, object_vertex_index, v);
				control_area_validity_(surface_point_index) = 1.0f;
			}

			return true;
		});

		this->attenuation_.resize(nbv_object);
		this->attenuation_.setZero();

		//compute_attenuation_cage(object);

	}*/

private:
	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> global_matrix_;

	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> global_matrix_comp_;

	float m_hFactor;

	Vec3 control_cage_center_;
	Vec3 control_cage_bb_min_;
	Vec3 control_cage_bb_max_;

	Vec3 influence_cage_bb_min_;
	Vec3 influence_cage_bb_max_;

	Local_direction_control_planes local_x_direction_control_planes_;
	Local_direction_control_planes local_y_direction_control_planes_;
	Local_direction_control_planes local_z_direction_control_planes_;

	Eigen::VectorXd control_area_validity_;

	Eigen::Matrix<Vec2, Eigen::Dynamic, Eigen::Dynamic> normal_weights_;

	std::shared_ptr<Attribute<uint32>> control_cage_vertex_index_;

	std::shared_ptr<Attribute<bool>> control_cage_marked_vertices_;

	void init_triangles()
	{

		foreach_cell(*control_cage_, [&](Face fc) -> bool {
			std::vector<CMap2::Vertex> face_vertices_ = incident_vertices(*control_cage_, fc);

			std::vector<CMap2::Vertex> triangle1_vertex(3);
			triangle1_vertex[0] = face_vertices_[1];
			triangle1_vertex[1] = face_vertices_[3];
			triangle1_vertex[2] = face_vertices_[0];

			std::vector<CMap2::Vertex> triangle2_vertex(3);
			triangle2_vertex[0] = face_vertices_[1];
			triangle2_vertex[1] = face_vertices_[2];
			triangle2_vertex[2] = face_vertices_[3];

			std::vector<Vec3> triangle1_position(3);
			std::vector<uint32_t> triangle1_index(3);

			std::vector<Vec3> triangle2_position(3);
			std::vector<uint32_t> triangle2_index(3);
			for (std::size_t i = 0; i < triangle1_vertex.size(); i++)
			{
				triangle1_position[i] = value<Vec3>(*control_cage_, control_cage_vertex_position_, triangle1_vertex[i]);

				triangle1_index[i] = value<uint32_t>(*control_cage_, control_cage_vertex_index_, triangle1_vertex[i]);

				triangle2_position[i] = value<Vec3>(*control_cage_, control_cage_vertex_position_, triangle2_vertex[i]);

				triangle2_index[i] = value<uint32_t>(*control_cage_, control_cage_vertex_index_, triangle2_vertex[i]);
			}

			const Vec3 t1_normal =
				(cgogn::geometry::normal(triangle1_position[0], triangle1_position[1], triangle1_position[2]))
					.normalized();

			const Vec3 t2_normal =
				(cgogn::geometry::normal(triangle2_position[0], triangle2_position[1], triangle2_position[2]))
					.normalized();

			Triangle triangle1;
			triangle1.vertices = triangle1_vertex;
			triangle1.positions = triangle1_position;
			triangle1.indices = triangle1_index;
			triangle1.normal = t1_normal;
			cage_triangles_.push_back(triangle1);

			Triangle triangle2;
			triangle2.vertices = triangle2_vertex;
			triangle2.positions = triangle2_position;
			triangle2.indices = triangle2_index;
			triangle2.normal = t2_normal;
			cage_triangles_.push_back(triangle2);

			return true;
		});
	}

	void init_control_cage_plane()
	{
		// TODO set for local frame
		const Vec3 x_dir = {1.0, 0.0, 0.0}, y_dir = {0.0, 1.0, 0.0}, z_dir = {0.0, 0.0, 1.0};

		double d_x_min = 1000.0, d_x_max = -1000.0, d_y_min = 1000.0, d_y_max = -1000.0, d_z_min = 1000.0,
			   d_z_max = -1000.0;

		std::pair<Triangle, Triangle> face_x_min, face_x_max, face_y_min, face_y_max, face_z_min, face_z_max;

		for (std::size_t i = 0; i < 6; i++)
		{ // for each face
			const Triangle triangle1 = cage_triangles_[2 * i];
			const Triangle triangle2 = cage_triangles_[2 * i + 1];

			const Vec3 normal = triangle1.normal;

			const Vec3 triangle_point = triangle1.positions[0];

			if (normal.dot(x_dir) != 0.0)
			{
				const double d = -(triangle_point.dot(x_dir));
				if (d < d_x_min)
				{
					d_x_min = d;
					face_x_min = std::make_pair(triangle1, triangle2);
				}

				if (d > d_x_max)
				{
					d_x_max = d;
					face_x_max = std::make_pair(triangle1, triangle2);
				}
			}
			else if (normal.dot(y_dir) != 0.0)
			{
				const double d = -(triangle_point.dot(y_dir));
				if (d < d_y_min)
				{
					d_y_min = d;
					face_y_min = std::make_pair(triangle1, triangle2);
				}

				if (d > d_y_max)
				{
					d_y_max = d;
					face_y_max = std::make_pair(triangle1, triangle2);
				}
			}
			else
			{
				const double d = -(triangle_point.dot(z_dir));
				if (d < d_z_min)
				{
					d_z_min = d;
					face_z_min = std::make_pair(triangle1, triangle2);
				}

				if (d > d_z_max)
				{
					d_z_max = d;
					face_z_max = std::make_pair(triangle1, triangle2);
				}
			}
		}

		local_x_direction_control_planes_.d_min = d_x_min;
		local_x_direction_control_planes_.d_max = d_x_max;
		local_x_direction_control_planes_.d_gap = d_x_max - d_x_min;
		local_x_direction_control_planes_.direction = x_dir;
		local_x_direction_control_planes_.triangles_d_min = face_x_min;
		local_x_direction_control_planes_.triangles_d_max = face_x_max;
		local_x_direction_control_planes_.shift_after_d_max = {-(2.0 * d_x_max - d_x_min), 0.0, 0.0}; // (d_max + d_gap)
		local_x_direction_control_planes_.shift_before_d_min = {-(2.0 * d_x_min - d_x_max), 0.0,
																0.0}; // (d_min - d_gap)

		local_y_direction_control_planes_.d_min = d_y_min;
		local_y_direction_control_planes_.d_max = d_y_max;
		local_y_direction_control_planes_.d_gap = d_y_max - d_y_min;
		local_y_direction_control_planes_.direction = y_dir;
		local_y_direction_control_planes_.triangles_d_min = face_y_min;
		local_y_direction_control_planes_.triangles_d_max = face_y_max;
		local_y_direction_control_planes_.shift_after_d_max = {0.0, -(2.0 * d_y_max - d_y_min), 0.0}; // (d_max + d_gap)
		local_y_direction_control_planes_.shift_before_d_min = {0.0, -(2.0 * d_y_min - d_y_max),
																0.0}; // (d_min - d_gap)

		local_z_direction_control_planes_.d_min = d_z_min;
		local_z_direction_control_planes_.d_max = d_z_max;
		local_z_direction_control_planes_.d_gap = d_z_max - d_z_min;
		local_z_direction_control_planes_.direction = z_dir;
		local_z_direction_control_planes_.triangles_d_min = face_z_min;
		local_z_direction_control_planes_.triangles_d_max = face_z_max;
		local_z_direction_control_planes_.shift_after_d_max = {0.0, 0.0, -(2.0 * d_z_max - d_z_min)}; // (d_max + d_gap)
		local_z_direction_control_planes_.shift_before_d_min = {0.0, 0.0,
																-(2.0 * d_z_min - d_z_max)}; // (d_min - d_gap)
	}

	std::vector<Triangle> get_virtual_cube_triangles(const std::pair<Triangle, Triangle> face, const Vec3& shift_vector)
	{
		std::vector<Triangle> virtual_cube_triangles;

		Triangle triangle1 = face.first;
		triangle1.virtual_indices = {1, 3, 0};
		Triangle triangle2 = face.second;
		triangle2.virtual_indices = {1, 2, 3};
		virtual_cube_triangles.push_back(triangle1);
		virtual_cube_triangles.push_back(triangle2);

		std::vector<Vec3> new_positions;

		for (std::size_t i = 0; i < 3; i++)
		{
			const Vec3 shifted_position = triangle1.positions[i] - shift_vector;
			new_positions.push_back(shifted_position);
		}

		for (std::size_t i = 0; i < 3; i++)
		{
			const Vec3 shifted_position = triangle2.positions[i] - shift_vector;
			new_positions.push_back(shifted_position);
		}

		Triangle triangle3, triangle4, triangle5, triangle6, triangle7, triangle8, triangle9, triangle10, triangle11,
			triangle12;

		triangle3.positions = {new_positions.begin(), new_positions.end() - 3};
		triangle3.virtual_indices = {5, 7, 4};
		// triangle3.normal = -1.0 * triangle1.normal;

		triangle4.positions = {new_positions.begin() + 3, new_positions.end()};
		triangle4.virtual_indices = {5, 6, 7};
		// triangle4.normal = -1.0 * triangle2.normal;

		triangle5.positions = {triangle1.positions[0], triangle3.positions[2], triangle1.positions[2]};
		triangle5.virtual_indices = {triangle1.virtual_indices[0], triangle3.virtual_indices[2],
									 triangle1.virtual_indices[2]};

		triangle6.positions = {triangle1.positions[0], triangle3.positions[0], triangle3.positions[2]};
		triangle6.virtual_indices = {triangle1.virtual_indices[0], triangle3.virtual_indices[0],
									 triangle3.virtual_indices[2]};

		triangle7.positions = {triangle1.positions[2], triangle3.positions[1], triangle3.positions[2]};
		triangle7.virtual_indices = {triangle1.virtual_indices[2], triangle3.virtual_indices[1],
									 triangle3.virtual_indices[2]};

		triangle8.positions = {triangle1.positions[2], triangle1.positions[1], triangle3.positions[1]};
		triangle8.virtual_indices = {triangle1.virtual_indices[2], triangle1.virtual_indices[1],
									 triangle3.virtual_indices[1]};

		triangle9.positions = {triangle1.positions[0], triangle4.positions[1], triangle3.positions[0]};
		triangle9.virtual_indices = {triangle1.virtual_indices[0], triangle4.virtual_indices[1],
									 triangle3.virtual_indices[0]};

		triangle10.positions = {triangle1.positions[0], triangle2.positions[1], triangle4.positions[1]};
		triangle10.virtual_indices = {triangle1.virtual_indices[0], triangle2.virtual_indices[1],
									  triangle4.virtual_indices[1]};

		triangle11.positions = {triangle4.positions[1], triangle1.positions[1], triangle3.positions[1]};
		triangle11.virtual_indices = {triangle4.virtual_indices[1], triangle1.virtual_indices[1],
									  triangle3.virtual_indices[1]};

		triangle12.positions = {triangle4.positions[1], triangle2.positions[1], triangle1.positions[1]};
		triangle12.virtual_indices = {triangle4.virtual_indices[1], triangle2.virtual_indices[1],
									  triangle1.virtual_indices[1]};

		virtual_cube_triangles.push_back(triangle3);
		virtual_cube_triangles.push_back(triangle4);
		virtual_cube_triangles.push_back(triangle5);
		virtual_cube_triangles.push_back(triangle6);
		virtual_cube_triangles.push_back(triangle7);
		virtual_cube_triangles.push_back(triangle8);
		virtual_cube_triangles.push_back(triangle9);
		virtual_cube_triangles.push_back(triangle10);
		virtual_cube_triangles.push_back(triangle11);
		virtual_cube_triangles.push_back(triangle12);

		return virtual_cube_triangles;
	}

	bool compute_mvc_on_point_inside_cage(const Vec3& surface_point, const uint32& surface_point_index)
	{
		uint32 nbv_cage = nb_cells<Vertex>(*control_cage_);

		std::shared_ptr<Attribute<uint32>> cage_vertex_index =
			get_attribute<uint32, Vertex>(*control_cage_, "vertex_index");

		double epsilon = 0.00000001;
		double sumWeights = 0.0;

		Eigen::VectorXd w_control_cage_coords_;

		w_control_cage_coords_.resize(nbv_cage);
		w_control_cage_coords_.setZero();

		std::vector<double> d(nbv_cage);
		std::vector<Vec3> u(nbv_cage);

		parallel_foreach_cell(*control_cage_, [&](Vertex v) -> bool {
			const Vec3& cage_point = value<Vec3>(*control_cage_, control_cage_vertex_position_, v);

			uint32 cage_point_index = value<uint32>(*control_cage_, control_cage_vertex_index_, v);

			d[cage_point_index] = (surface_point - cage_point).norm();
			if (d[cage_point_index] < epsilon)
			{
				control_cage_coords_(surface_point_index, cage_point_index) = 1.0;
				return true;
			}

			u[cage_point_index] = (cage_point - surface_point) / d[cage_point_index];

			return true;
		});

		double l[3];
		double theta[3];
		double w[3];
		double c[3];
		double s[3];

		for (std::size_t t = 0; t < cage_triangles_.size(); t++)
		{

			std::vector<uint32> triangle_index = cage_triangles_[t].indices;

			for (std::size_t i = 0; i < 3; i++)
			{
				l[i] = (u[triangle_index[(i + 1) % 3]] - u[triangle_index[(i + 2) % 3]]).norm();
			}

			for (std::size_t i = 0; i < 3; i++)
			{
				theta[i] = 2.0 * asin(l[i] / 2.0);
			}

			double h = (theta[0] + theta[1] + theta[2]) / 2.0;
			if (M_PI - h < epsilon)
			{
				for (std::size_t i = 0; i < 3; i++)
				{
					w[i] = sin(theta[i]) * l[(i + 2) % 3] * l[(i + 1) % 3];
				}

				sumWeights = w[0] + w[1] + w[2];
				w_control_cage_coords_[triangle_index[0]] = w[0] / sumWeights;
				w_control_cage_coords_[triangle_index[1]] = w[1] / sumWeights;
				w_control_cage_coords_[triangle_index[2]] = w[2] / sumWeights;

				return true;
			}

			for (std::size_t i = 0; i < 3; i++)
			{
				c[i] = (2.0 * sin(h) * sin(h - theta[i])) / (sin(theta[(i + 1) % 3]) * sin(theta[(i + 2) % 3])) - 1.0;
			}

			double sign_Basis_u0u1u2 = 1;
			Vec3 crossVec = u[triangle_index[0]].cross(u[triangle_index[1]]);
			if (crossVec.dot(u[triangle_index[2]]) < 0.0)
			{
				sign_Basis_u0u1u2 = -1;
			}

			for (std::size_t i = 0; i < 3; i++)
			{
				s[i] = sign_Basis_u0u1u2 * sqrt(std::max<double>(0.0, 1.0 - c[i] * c[i]));
			}

			if (fabs(s[0]) < epsilon || fabs(s[1]) < epsilon || fabs(s[2]) < epsilon)
			{
				continue; // eta is on the same plane, outside t  ->  ignore triangle t :
			}

			for (std::size_t i = 0; i < 3; ++i)
			{
				w[i] = (theta[i] - c[(i + 1) % 3] * theta[(i + 2) % 3] - c[(i + 2) % 3] * theta[(i + 1) % 3]) /
					   (2.0 * d[triangle_index[i]] * sin(theta[(i + 1) % 3]) * s[(i + 2) % 3]);
			}

			sumWeights += (w[0] + w[1] + w[2]);
			w_control_cage_coords_[triangle_index[0]] += w[0];
			w_control_cage_coords_[triangle_index[1]] += w[1];
			w_control_cage_coords_[triangle_index[2]] += w[2];
		}

		parallel_foreach_cell(*control_cage_, [&](Vertex v) -> bool {
			uint32 cage_point_index = value<uint32>(*control_cage_, control_cage_vertex_index_, v);

			control_cage_coords_(surface_point_index, cage_point_index) =
				w_control_cage_coords_[cage_point_index] / sumWeights;

			return true;
		});

		return false;
	}

	bool compute_mvc_on_point_outside_cage(const Vec3& surface_point, const uint32& surface_point_index,
										   const std::vector<Triangle> virtual_cube_triangles)
	{
		uint32 nbv_cage = 8; // nb_cells<Vertex>(*control_cage_);

		double epsilon = 0.00000001;
		double sumWeights = 0.0;

		Eigen::VectorXd w_control_cage_coords_;

		w_control_cage_coords_.resize(nbv_cage);
		w_control_cage_coords_.setZero();

		std::vector<double> d(nbv_cage);
		std::vector<Vec3> u(nbv_cage);

		const Triangle triangle0 = virtual_cube_triangles[0];
		const Triangle triangle1 = virtual_cube_triangles[1];
		const Triangle triangle2 = virtual_cube_triangles[2];
		const Triangle triangle3 = virtual_cube_triangles[3];

		// TODO: find better solution
		std::vector<std::pair<Vec3, uint32_t>> virtual_cage_positions(8);
		virtual_cage_positions.push_back(std::make_pair(triangle0.positions[2], triangle0.virtual_indices[2]));
		virtual_cage_positions.push_back(std::make_pair(triangle1.positions[0], triangle1.virtual_indices[0]));
		virtual_cage_positions.push_back(std::make_pair(triangle1.positions[1], triangle1.virtual_indices[1]));
		virtual_cage_positions.push_back(std::make_pair(triangle1.positions[2], triangle1.virtual_indices[2]));
		virtual_cage_positions.push_back(std::make_pair(triangle2.positions[2], triangle2.virtual_indices[2]));
		virtual_cage_positions.push_back(std::make_pair(triangle3.positions[0], triangle3.virtual_indices[0]));
		virtual_cage_positions.push_back(std::make_pair(triangle3.positions[1], triangle3.virtual_indices[1]));
		virtual_cage_positions.push_back(std::make_pair(triangle3.positions[2], triangle3.virtual_indices[2]));

		for (std::size_t p = 0; p < 8; p++)
		{
			const Vec3& cage_point = virtual_cage_positions[p].first;

			uint32 cage_point_index = virtual_cage_positions[p].second;

			d[cage_point_index] = (surface_point - cage_point).norm();
			if (d[cage_point_index] < epsilon)
			{
				control_cage_coords_(surface_point_index, cage_point_index) = 1.0;
				return true;
			}

			u[cage_point_index] = (cage_point - surface_point) / d[cage_point_index];
		}

		double l[3];
		double theta[3];
		double w[3];
		double c[3];
		double s[3];

		for (std::size_t t = 0; t < virtual_cube_triangles.size(); t++)
		{

			std::vector<uint32> triangle_index = virtual_cube_triangles[t].virtual_indices;

			for (std::size_t i = 0; i < 3; i++)
			{
				l[i] = (u[triangle_index[(i + 1) % 3]] - u[triangle_index[(i + 2) % 3]]).norm();
			}

			for (std::size_t i = 0; i < 3; i++)
			{
				theta[i] = 2.0 * asin(l[i] / 2.0);
			}

			double h = (theta[0] + theta[1] + theta[2]) / 2.0;
			if (M_PI - h < epsilon)
			{
				for (std::size_t i = 0; i < 3; i++)
				{
					w[i] = sin(theta[i]) * l[(i + 2) % 3] * l[(i + 1) % 3];
				}

				sumWeights = w[0] + w[1] + w[2];
				w_control_cage_coords_[triangle_index[0]] = w[0] / sumWeights;
				w_control_cage_coords_[triangle_index[1]] = w[1] / sumWeights;
				w_control_cage_coords_[triangle_index[2]] = w[2] / sumWeights;

				return true;
			}

			for (std::size_t i = 0; i < 3; i++)
			{
				c[i] = (2.0 * sin(h) * sin(h - theta[i])) / (sin(theta[(i + 1) % 3]) * sin(theta[(i + 2) % 3])) - 1.0;
			}

			double sign_Basis_u0u1u2 = 1;
			Vec3 crossVec = u[triangle_index[0]].cross(u[triangle_index[1]]);
			if (crossVec.dot(u[triangle_index[2]]) < 0.0)
			{
				sign_Basis_u0u1u2 = -1;
			}

			for (std::size_t i = 0; i < 3; i++)
			{
				s[i] = sign_Basis_u0u1u2 * sqrt(std::max<double>(0.0, 1.0 - c[i] * c[i]));
			}

			if (fabs(s[0]) < epsilon || fabs(s[1]) < epsilon || fabs(s[2]) < epsilon)
			{
				continue; // eta is on the same plane, outside t  ->  ignore triangle t :
			}

			for (std::size_t i = 0; i < 3; ++i)
			{
				w[i] = (theta[i] - c[(i + 1) % 3] * theta[(i + 2) % 3] - c[(i + 2) % 3] * theta[(i + 1) % 3]) /
					   (2.0 * d[triangle_index[i]] * sin(theta[(i + 1) % 3]) * s[(i + 2) % 3]);
			}

			sumWeights += (w[0] + w[1] + w[2]);
			w_control_cage_coords_[triangle_index[0]] += w[0];
			w_control_cage_coords_[triangle_index[1]] += w[1];
			w_control_cage_coords_[triangle_index[2]] += w[2];
		}

		for (std::size_t p = 0; p < 8; p++)
		{
			uint32 cage_point_index = virtual_cage_positions[p].second;
			control_cage_coords_(surface_point_index, cage_point_index) =
				w_control_cage_coords_[cage_point_index] / sumWeights;
		}

		return false;
	}

	void bind_mvc_customed_vertices(const Vec3 surface_point, uint32 surface_point_index,
									std::vector<Vec3> position_vertices)
	{
	}

	bool check_mvc_on_point_inside_influence_cage(Vec3 point)
	{
		std::shared_ptr<Attribute<bool>> cage_vertex_marked =
			get_attribute<bool, Vertex>(*influence_cage_, "marked_vertices");

		DartMarker dm(*influence_cage_);

		bool checked = true;
		for (Dart d = influence_cage_->begin(), end = influence_cage_->end(); d != end; d = influence_cage_->next(d))
		{
			Vertex cage_vertex = CMap2::Vertex(d);
			bool vc_marked = value<bool>(*influence_cage_, cage_vertex_marked, cage_vertex);

			if (!dm.is_marked(d) && !vc_marked)
			{

				const Vec3& cage_point = value<Vec3>(*influence_cage_, influence_cage_vertex_position_, cage_vertex);

				float mvc_value =
					compute_mvc(point, d, *influence_cage_, cage_point, influence_cage_vertex_position_.get());

				dm.mark(d);

				value<bool>(*influence_cage_, cage_vertex_marked, cage_vertex) = true;

				if (mvc_value < 0)
				{
					checked = false;
					break;
				}
			}
		}

		parallel_foreach_cell(*influence_cage_, [&](Vertex vc) -> bool {
			value<bool>(*influence_cage_, cage_vertex_marked, vc) = false;
			return true;
		});

		return checked;
	}

	std::vector<Vec3> find_intersection_positions_face(const std::pair<Triangle, Triangle>& face1,
															   const std::pair<Triangle, Triangle>& face2)
	{

		std::vector<Vec3> intersect_positions; 
		double epsilon = 0.00000001;

		std::vector<Vec3> positions_face_1;
		positions_face_1.push_back(face1.first.positions[2]);
		positions_face_1.push_back(face1.first.positions[0]);
		positions_face_1.push_back(face1.second.positions[1]);
		positions_face_1.push_back(face1.first.positions[1]);

		std::vector<Vec3> positions_face_2;
		positions_face_2.push_back(face2.first.positions[2]);
		positions_face_2.push_back(face2.first.positions[0]);
		positions_face_2.push_back(face2.second.positions[1]);
		positions_face_2.push_back(face2.first.positions[1]);


		// TODO: find alternative to double loops
		for (std::size_t i = 0; i < 4; i++)
		{
			Vec3 target_position_face1 = positions_face_1[i];
			for (std::size_t j = 0; j < 4; j++)
			{
				Vec3 target_position_face2 = positions_face_2[j];

				const double distance_p1_p2 = (target_position_face1 - target_position_face2).norm();
				if (distance_p1_p2 < epsilon)
				{
					intersect_positions.push_back(target_position_face1);
					continue;
				}
			}
		}

		if (!intersect_positions.size() == 2){
			std::cout << "Warning, size not 2: " << intersect_positions.size() << std::endl; 
		}

		std::cout << "size intersect positions " << intersect_positions.size() << std::endl; 

		return intersect_positions;
	}

	Vec3 find_intersection_position(const std::pair<Triangle, Triangle>& face1,
										   const std::pair<Triangle, Triangle>& face2,
										   const std::pair<Triangle, Triangle>& face3)
	{
		Vec3 intersection_position;

		double epsilon = 0.00000001;

		std::vector<Vec3> position_face_1;
		position_face_1.push_back(face1.first.positions[2]);
		position_face_1.push_back(face1.first.positions[0]);
		position_face_1.push_back(face1.second.positions[1]);
		position_face_1.push_back(face1.first.positions[1]);

		std::vector<Vec3> position_face_2;
		position_face_2.push_back(face2.first.positions[2]);
		position_face_2.push_back(face2.first.positions[0]);
		position_face_2.push_back(face2.second.positions[1]);
		position_face_2.push_back(face2.first.positions[1]);

		std::vector<Vec3> position_face_3;
		position_face_3.push_back(face3.first.positions[2]);
		position_face_3.push_back(face3.first.positions[0]);
		position_face_3.push_back(face3.second.positions[1]);
		position_face_3.push_back(face3.first.positions[1]);

		// TODO: find alternative to double loops
		for (std::size_t i = 0; i < 4; i++)
		{
			Vec3 target_position_face1 = position_face_1[i];
			for (std::size_t j = 0; j < 4; j++)
			{
				Vec3 target_position_face2 = position_face_2[j];
				const double distance_p1_p2 = (target_position_face1 - target_position_face2).norm();
				if (distance_p1_p2 < epsilon)
				{
					for (std::size_t k = 0; k < 4; k++)
					{
						Vec3 target_position_face3 = position_face_3[k];
						const double distance_p1_p3 = (target_position_face1 - target_position_face3).norm();
						if (distance_p1_p3 < epsilon)
						{
							return target_position_face1; 
						}
					}
				}
			}
		}

		return intersection_position;
	}

	std::pair<Triangle, Triangle> get_face_from_intersecting_edge(const std::vector<Vec3>& intersect_positions,
																  const Vec3& shiftVector)
	{
		// create new face composed of the intersecting edge and this edge shifted by y
		std::cout << "first " << intersect_positions[0] << std::endl; 
		std::cout << "second " << intersect_positions[1] << std::endl; 
		const Vec3 face_position0 = intersect_positions[0];
		const Vec3 face_position1 = intersect_positions[1];

		std::cout << "test here " << std::endl; 

		const Vec3 face_position2 = face_position1 + shiftVector;
		const Vec3 face_position3 = face_position0 + shiftVector;

		std::cout << "valid here " << std::endl; 

		Triangle local_triangle1, local_triangle2;
		local_triangle1.positions = {face_position1, face_position3, face_position0};
		local_triangle2.positions = {face_position1, face_position2, face_position3};

		std::cout << "ok end function " << std::endl; 

		return std::make_pair(local_triangle1, local_triangle2);
	}

	std::pair<Triangle, Triangle> get_face_from_intersecting_vertex(const Vec3& intersection_position,
																	const Vec3& shiftVector1, const Vec3& shiftVector2)
	{
		// create new face composed of the intersection vertex and shifted in two directions
		const Vec3 face_position0 = intersection_position;

		const Vec3 face_position1 = face_position0 + shiftVector1;

		const Vec3 face_position2 = (face_position1 + shiftVector1) + shiftVector2;
		const Vec3 face_position3 = face_position0 + shiftVector2;

		Triangle local_triangle1, local_triangle2;
		local_triangle1.positions = {face_position1, face_position3, face_position0};
		local_triangle2.positions = {face_position1, face_position2, face_position3};

		return std::make_pair(local_triangle1, local_triangle2);
	}
};

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_CAGE_DEFORMATION_TOOL_H_

/*void compute_attenuation_cage(MESH& object)
{

std::shared_ptr<Attribute<uint32>> cage_face_indices =
	add_attribute<uint32, Face>(*(this->influence_cage_), "face_indices");

cgogn::modeling::set_attribute_face_index(*(this->influence_cage_),
											cage_face_indices.get());

std::shared_ptr<Attribute<uint32>> object_vertex_index =
	get_attribute<uint32, Vertex>(object, "vertex_index");

std::shared_ptr<Attribute<uint32>> cage_vertex_index =
	get_attribute<uint32, Vertex>(*control_cage_, "vertex_index");

std::shared_ptr<Attribute<uint32>> i_cage_vertex_index =
	get_attribute<uint32, Vertex>(*(this->influence_cage_), "vertex_index");

uint32 nbf_cage = 2 * nb_cells<Face>(*(this->influence_cage_));
uint32 nbv_cage = nb_cells<Vertex>(*(this->influence_cage_));

// first loop to find h
//
float h = 0.0f;
float max_dist = 0.0f;
std::vector<Vec2> attenuation_points;
this->influence_area_->foreach_cell([&](Vertex v) {
	uint32 surface_point_index = value<uint32>(object, object_vertex_index, v);

	float i_dist = this->cage_influence_distance(surface_point_index, nbf_cage, nbv_cage);

	this->attenuation_(surface_point_index) = (float)sin(0.5*M_PI * (i_dist ));
	/*if (control_area_validity_(surface_point_index) == 1.0f)
	{

		if (i_dist > h)
		{
			h = i_dist;
		}

		this->attenuation_(surface_const double point_index) = 1.0f;
	}
	else
	{

		if (i_dist > max_dist){
			max_dist = i_dist;
		}
		attenuation_points.push_back({surface_point_index, i_dist});
	}*/
//});

/*for (unsigned int i = 0; i < attenuation_points.size(); i++)
{
	//this->attenuation_(attenuation_points[i][0]) = 0.5f * ((float)sin(M_PI * ((attenuation_points[i][1] / h) -
0.5f))) + 0.5f; this->attenuation_(attenuation_points[i][0]) = (float)sin(0.5*M_PI * (attenuation_points[i][1] /
h));
}*/
//}

/*bool local_mvc_pt_control_area(Vec3 pt)
{
	std::shared_ptr<Attribute<uint32>> i_vertex_index =
		get_attribute<uint32, Vertex>(*control_cage_, "vertex_index");

	std::shared_ptr<Attribute<bool>> cage_vertex_marked =
		get_attribute<bool, Vertex>(*control_cage_, "marked_vertices");

	DartMarker dm(*control_cage_);

	bool checked = true;
	for (Dart d = control_cage_->begin(), end = control_cage_->end(); d != end; d = control_cage_->next(d))
	{
		Vertex cage_vertex = CMap2::Vertex(d);
		bool vc_marked = value<bool>(*control_cage_, cage_vertex_marked, cage_vertex);

		if (!dm.is_marked(d) && !vc_marked)
		{

			const Vec3& cage_point = value<Vec3>(*control_cage_, control_cage_vertex_position_, cage_vertex);

			float mvc_value = cgogn::modeling::compute_mvc(pt, d, *control_cage_, cage_point,
control_cage_vertex_position_.get());

			dm.mark(d);

			value<bool>(*control_cage_, cage_vertex_marked, cage_vertex) = true;

			if (mvc_value < 0)
			{
				checked = false;
				break;
			}
		}
	}

	parallel_foreach_cell(*control_cage_, [&](Vertex vc) -> bool {
		value<bool>(*control_cage_, cage_vertex_marked, vc) = false;
		return true;
	});

	return checked;
}*/
