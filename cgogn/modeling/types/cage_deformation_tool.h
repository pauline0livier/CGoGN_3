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

#include <cgogn/core/types/cmap/cmap2.h>
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
		Vec3 normal;
		std::pair<Vec3, Vec3> edges;
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

			const double d_x = -(surface_point.dot(std::get<0>(local_x_direction_control_planes_))),

						 d_y = -(surface_point.dot(std::get<0>(local_y_direction_control_planes_))),

						 d_z = -(surface_point.dot(std::get<0>(local_z_direction_control_planes_)));

			const bool valid_x_dir = (d_x <= std::get<2>(local_x_direction_control_planes_) &&
									  d_x >= std::get<1>(local_x_direction_control_planes_)),

					   valid_y_dir = (d_y <= std::get<2>(local_y_direction_control_planes_) &&
									  d_y >= std::get<1>(local_y_direction_control_planes_)),

					   valid_z_dir = (d_z <= std::get<2>(local_z_direction_control_planes_) &&
									  d_z >= std::get<1>(local_z_direction_control_planes_));

			if (valid_x_dir)
			{
				if (valid_y_dir && valid_z_dir)
				{
					// inside cube
					compute_mvc_on_point_inside_cage(surface_point, surface_point_index);
				}
				else if (valid_y_dir)
				{
					if (d_z > std::get<2>(local_z_direction_control_planes_))
					{
						// fakeCube[5]

						// Beware, need to inverse the normal of the triangle1 and triangle2
						const Triangle local_triangle1 = control_cage_face_d_z_max_.first;
						const Triangle local_triangle2 = control_cage_face_d_z_max_.second;

						// find plane of "fake" cube
						const double gap_z_plane = std::get<2>(local_z_direction_control_planes_) -
												   std::get<1>(local_z_direction_control_planes_);

						const double new_plane_d = std::get<2>(local_z_direction_control_planes_) + gap_z_plane;

						const Vec3 shift = {0.0, 0.0, -new_plane_d};

						std::vector<Triangle> virtual_cube_triangles =
							get_virtual_cube_triangles(local_triangle1, local_triangle2, shift);

						// compute_mvc_on_point_outside_cage(surface_point, surface_point_index,
						// virtual_cube_triangles);
					}
					else
					{
						// fakeCube[4]
						const Triangle local_triangle1 = control_cage_face_d_z_min_.first;
						const Triangle local_triangle2 = control_cage_face_d_z_min_.second;

						const double gap_z_plane = std::get<2>(local_z_direction_control_planes_) -
												   std::get<1>(local_z_direction_control_planes_);

						const double new_plane_d = std::get<1>(local_z_direction_control_planes_) - gap_z_plane;

						const Vec3 shift = {0.0, 0.0, -new_plane_d};

						std::vector<Triangle> virtual_cube_triangles =
							get_virtual_cube_triangles(local_triangle1, local_triangle2, shift);
					}
				}
				else if (valid_z_dir)
				{
					if (d_y > std::get<2>(local_y_direction_control_planes_))
					{
						// fakeCube[3]
						const Triangle local_triangle1 = control_cage_face_d_y_max_.first;
						const Triangle local_triangle2 = control_cage_face_d_y_max_.second;

						const double gap_y_plane = std::get<2>(local_y_direction_control_planes_) -
												   std::get<1>(local_y_direction_control_planes_);

						const double new_plane_d = std::get<2>(local_y_direction_control_planes_) + gap_y_plane;

						const Vec3 shift = {0.0, -new_plane_d, 0.0};

						std::vector<Triangle> virtual_cube_triangles =
							get_virtual_cube_triangles(local_triangle1, local_triangle2, shift);
					}
					else
					{
						// fakeCube[2]
						const Triangle local_triangle1 = control_cage_face_d_y_min_.first;
						const Triangle local_triangle2 = control_cage_face_d_y_min_.second;

						const double gap_y_plane = std::get<2>(local_y_direction_control_planes_) -
												   std::get<1>(local_y_direction_control_planes_);

						const double new_plane_d = std::get<1>(local_y_direction_control_planes_) - gap_y_plane;

						const Vec3 shift = {0.0, -new_plane_d, 0.0};

						std::vector<Triangle> virtual_cube_triangles =
							get_virtual_cube_triangles(local_triangle1, local_triangle2, shift);
					}
				}
				else
				{
					if (d_y > std::get<2>(local_y_direction_control_planes_))
					{
						const Triangle y_triangle1 = control_cage_face_d_y_max_.first;
						const Triangle y_triangle2 = control_cage_face_d_y_max_.second;

						const double gap_y_plane = std::get<2>(local_y_direction_control_planes_) -
												   std::get<1>(local_y_direction_control_planes_);

						const double new_plane_d_y = std::get<2>(local_y_direction_control_planes_) + gap_y_plane;

						const Vec3 shift_y = {0.0, -new_plane_d_y, 0.0};

						if (d_z > std::get<2>(local_z_direction_control_planes_))
						{
							// fakeCubeDiag[11]
							const Triangle z_triangle1 = control_cage_face_d_z_max_.first;
							const Triangle z_triangle2 = control_cage_face_d_z_max_.second;

							// find plane of "fake" cube
							const double gap_z_plane = std::get<2>(local_z_direction_control_planes_) -
													   std::get<1>(local_z_direction_control_planes_);

							const double new_plane_d_z = std::get<2>(local_z_direction_control_planes_) + gap_z_plane;

							const Vec3 shift_z = {0.0, 0.0, -new_plane_d_z};

							// find intersection edge between y and z triangle
							const std::vector<CMap2::Vertex> vertex_face_y;
							vertex_face_y.push_back(y_triangle1.vertices[2]);
							vertex_face_y.push_back(y_triangle1.vertices[0]);
							vertex_face_y.push_back(y_triangle2.vertices[1]);
							vertex_face_y.push_back(y_triangle1.vertices[1]);

							const std::vector<CMap2::Vertex> vertex_face_z;
							vertex_face_z.push_back(z_triangle1.vertices[2]);
							vertex_face_z.push_back(z_triangle1.vertices[0]);
							vertex_face_z.push_back(z_triangle2.vertices[1]);
							vertex_face_z.push_back(z_triangle1.vertices[1]);

							std::vector<CMap::Vertex> intersect_vertices;
							// TODO: find alternative to double loops
							for (std::size_t i = 0; i < 4; i++)
							{
								CMap2::Vertex target_y_vertex = vertex_face_y[i];
								for (std::size_t j = 0; j < 4; j++)
								{
									CMap2::Vertex target_z_vertex = vertex_face_z[j];

									if (target_y_vertex == target_z_vertex)
									{
										intersect_vertices.push_back(target_y_vertex);
									}
								}
							}

							// create new face composed of the intersecting edge and this edge shifted by y
							const Vec3 face_position0 =
								value<Vec3>(*control_cage_, control_cage_vertex_position_, intersect_vertices[0]);
							const Vec3 face_position1 =
								value<Vec3>(*control_cage_, control_cage_vertex_position_, intersect_vertices[1]);

							const Vec3 face_position2 = face_position1 + shift_y;
							const Vec3 face_position3 = face_position0 + shift_y;

							const Triangle local_triangle1, triangle2;
							local_triangle1.positions = {face_position1, face_position3, face_position0};
							local_triangle2.positions = {face_position1, face_position2, face_position3};

							std::vector<Triangle> virtual_cube_triangles =
								get_virtual_cube_triangles(local_triangle1, local_triangle2, shift_z);
						}
						else
						{
							// fakeCubeDiag[9]
							const Triangle z_triangle1 = control_cage_face_d_z_min_.first;
							const Triangle z_triangle2 = control_cage_face_d_z_min_.second;

							const double gap_z_plane = std::get<2>(local_z_direction_control_planes_) -
													   std::get<1>(local_z_direction_control_planes_);

							const double new_plane_d_z = std::get<1>(local_z_direction_control_planes_) - gap_z_plane;

							const Vec3 shift_z = {0.0, 0.0, -new_plane_d_z};

							// find intersection edge between y and z triangle
							const std::vector<CMap2::Vertex> vertex_face_y;
							vertex_face_y.push_back(y_triangle1.vertices[2]);
							vertex_face_y.push_back(y_triangle1.vertices[0]);
							vertex_face_y.push_back(y_triangle2.vertices[1]);
							vertex_face_y.push_back(y_triangle1.vertices[1]);

							const std::vector<CMap2::Vertex> vertex_face_z;
							vertex_face_z.push_back(z_triangle1.vertices[2]);
							vertex_face_z.push_back(z_triangle1.vertices[0]);
							vertex_face_z.push_back(z_triangle2.vertices[1]);
							vertex_face_z.push_back(z_triangle1.vertices[1]);

							std::vector<CMap::Vertex> intersect_vertices;
							// TODO: find alternative to double loops
							for (std::size_t i = 0; i < 4; i++)
							{
								CMap2::Vertex target_y_vertex = vertex_face_y[i];
								for (std::size_t j = 0; j < 4; j++)
								{
									CMap2::Vertex target_z_vertex = vertex_face_z[j];

									if (target_y_vertex == target_z_vertex)
									{
										intersect_vertices.push_back(target_y_vertex);
									}
								}
							}

							// create new face composed of the intersecting edge and this edge shifted by y
							const Vec3 face_position0 =
								value<Vec3>(*control_cage_, control_cage_vertex_position_, intersect_vertices[0]);
							const Vec3 face_position1 =
								value<Vec3>(*control_cage_, control_cage_vertex_position_, intersect_vertices[1]);

							const Vec3 face_position2 = face_position1 + shift_y;
							const Vec3 face_position3 = face_position0 + shift_y;

							const Triangle local_triangle1, triangle2;
							local_triangle1.positions = {face_position1, face_position3, face_position0};
							local_triangle2.positions = {face_position1, face_position2, face_position3};

							std::vector<Triangle> virtual_cube_triangles =
								get_virtual_cube_triangles(local_triangle1, local_triangle2, shift_z);
						}
					}
					else
					{
						const Triangle y_triangle1 = control_cage_face_d_y_min_.first;
						const Triangle y_triangle2 = control_cage_face_d_y_min_.second;

						const double gap_y_plane = std::get<2>(local_y_direction_control_planes_) -
												   std::get<1>(local_y_direction_control_planes_);

						const double new_plane_d_y = std::get<1>(local_y_direction_control_planes_) - gap_y_plane;

						const Vec3 shift_y = {0.0, -new_plane_d_y, 0.0};
						if (d_z > std::get<2>(local_z_direction_control_planes_))
						{
							// fakeCubeDiag[10]
							const Triangle z_triangle1 = control_cage_face_d_z_max_.first;
							const Triangle z_triangle2 = control_cage_face_d_z_max_.second;

							// find plane of "fake" cube
							const double gap_z_plane = std::get<2>(local_z_direction_control_planes_) -
													   std::get<1>(local_z_direction_control_planes_);

							const double new_plane_d_z = std::get<2>(local_z_direction_control_planes_) + gap_z_plane;

							const Vec3 shift_z = {0.0, 0.0, -new_plane_d_z};

							// find intersection edge between y and z triangle
							const std::vector<CMap2::Vertex> vertex_face_y;
							vertex_face_y.push_back(y_triangle1.vertices[2]);
							vertex_face_y.push_back(y_triangle1.vertices[0]);
							vertex_face_y.push_back(y_triangle2.vertices[1]);
							vertex_face_y.push_back(y_triangle1.vertices[1]);

							const std::vector<CMap2::Vertex> vertex_face_z;
							vertex_face_z.push_back(z_triangle1.vertices[2]);
							vertex_face_z.push_back(z_triangle1.vertices[0]);
							vertex_face_z.push_back(z_triangle2.vertices[1]);
							vertex_face_z.push_back(z_triangle1.vertices[1]);

							std::vector<CMap::Vertex> intersect_vertices;
							// TODO: find alternative to double loops
							for (std::size_t i = 0; i < 4; i++)
							{
								CMap2::Vertex target_y_vertex = vertex_face_y[i];
								for (std::size_t j = 0; j < 4; j++)
								{
									CMap2::Vertex target_z_vertex = vertex_face_z[j];

									if (target_y_vertex == target_z_vertex)
									{
										intersect_vertices.push_back(target_y_vertex);
									}
								}
							}

							// create new face composed of the intersecting edge and this edge shifted by y
							const Vec3 face_position0 =
								value<Vec3>(*control_cage_, control_cage_vertex_position_, intersect_vertices[0]);
							const Vec3 face_position1 =
								value<Vec3>(*control_cage_, control_cage_vertex_position_, intersect_vertices[1]);

							const Vec3 face_position2 = face_position1 + shift_y;
							const Vec3 face_position3 = face_position0 + shift_y;

							const Triangle local_triangle1, triangle2;
							local_triangle1.positions = {face_position1, face_position3, face_position0};
							local_triangle2.positions = {face_position1, face_position2, face_position3};

							std::vector<Triangle> virtual_cube_triangles =
								get_virtual_cube_triangles(local_triangle1, local_triangle2, shift_z);
						}
						else
						{
							// fakeCubeDiag[8]
							const Triangle z_triangle1 = control_cage_face_d_z_min_.first;
							const Triangle z_triangle2 = control_cage_face_d_z_min_.second;

							const double gap_z_plane = std::get<2>(local_z_direction_control_planes_) -
													   std::get<1>(local_z_direction_control_planes_);

							const double new_plane_d_z = std::get<1>(local_z_direction_control_planes_) - gap_z_plane;

							const Vec3 shift_z = {0.0, 0.0, -new_plane_d_z};

							// find intersection edge between y and z triangle
							const std::vector<CMap2::Vertex> vertex_face_y;
							vertex_face_y.push_back(y_triangle1.vertices[2]);
							vertex_face_y.push_back(y_triangle1.vertices[0]);
							vertex_face_y.push_back(y_triangle2.vertices[1]);
							vertex_face_y.push_back(y_triangle1.vertices[1]);

							const std::vector<CMap2::Vertex> vertex_face_z;
							vertex_face_z.push_back(z_triangle1.vertices[2]);
							vertex_face_z.push_back(z_triangle1.vertices[0]);
							vertex_face_z.push_back(z_triangle2.vertices[1]);
							vertex_face_z.push_back(z_triangle1.vertices[1]);

							std::vector<CMap::Vertex> intersect_vertices;
							// TODO: find alternative to double loops
							for (std::size_t i = 0; i < 4; i++)
							{
								CMap2::Vertex target_y_vertex = vertex_face_y[i];
								for (std::size_t j = 0; j < 4; j++)
								{
									CMap2::Vertex target_z_vertex = vertex_face_z[j];

									if (target_y_vertex == target_z_vertex)
									{
										intersect_vertices.push_back(target_y_vertex);
									}
								}
							}

							// create new face composed of the intersecting edge and this edge shifted by y
							const Vec3 face_position0 =
								value<Vec3>(*control_cage_, control_cage_vertex_position_, intersect_vertices[0]);
							const Vec3 face_position1 =
								value<Vec3>(*control_cage_, control_cage_vertex_position_, intersect_vertices[1]);

							const Vec3 face_position2 = face_position1 + shift_y;
							const Vec3 face_position3 = face_position0 + shift_y;

							const Triangle local_triangle1, triangle2;
							local_triangle1.positions = {face_position1, face_position3, face_position0};
							local_triangle2.positions = {face_position1, face_position2, face_position3};

							std::vector<Triangle> virtual_cube_triangles =
								get_virtual_cube_triangles(local_triangle1, local_triangle2, shift_z);
						}
					}
				}
			}
			else if (valid_y_dir)
			{
				if (valid_z_dir)
				{
					if (d_x > std::get<2>(local_x_direction_control_planes_))
					{
						// fakeCube[1]
						const Triangle local_triangle1 = control_cage_face_d_x_max_.first;
						const Triangle local_triangle2 = control_cage_face_d_x_max_.second;

						const double gap_x_plane = std::get<2>(local_x_direction_control_planes_) -
												   std::get<1>(local_x_direction_control_planes_);

						const double new_plane_d = std::get<2>(local_x_direction_control_planes_) + gap_x_plane;

						const Vec3 shift = {-new_plane_d, 0.0, 0.0};

						std::vector<Triangle> virtual_cube_triangles =
							get_virtual_cube_triangles(local_triangle1, local_triangle2, shift);
					}
					else
					{
						// fakeCube[0]
						const Triangle local_triangle1 = control_cage_face_d_x_min_.first;
						const Triangle local_triangle2 = control_cage_face_d_x_min_.second;

						const double gap_x_plane = std::get<2>(local_x_direction_control_planes_) -
												   std::get<1>(local_x_direction_control_planes_);

						const double new_plane_d = std::get<1>(local_x_direction_control_planes_) - gap_x_plane;

						const Vec3 shift = {-new_plane_d, 0.0, 0.0};

						std::vector<Triangle> virtual_cube_triangles =
							get_virtual_cube_triangles(local_triangle1, local_triangle2, shift);
					}
				}
				else
				{
					if (d_x > std::get<2>(local_x_direction_control_planes_))
					{
						const Triangle x_triangle1 = control_cage_face_d_x_max_.first;
						const Triangle x_triangle2 = control_cage_face_d_x_max_.second;

						const double gap_x_plane = std::get<2>(local_x_direction_control_planes_) -
												   std::get<1>(local_x_direction_control_planes_);

						const double new_plane_d_x = std::get<2>(local_x_direction_control_planes_) + gap_x_plane;

						const Vec3 shift_x = {-new_plane_d_y, 0.0, 0.0};

						if (d_z > std::get<2>(local_z_direction_control_planes_))
						{
							// fakeCubeDiag[3]
							const Triangle z_triangle1 = control_cage_face_d_z_max_.first;
							const Triangle z_triangle2 = control_cage_face_d_z_max_.second;

							// find plane of "fake" cube
							const double gap_z_plane = std::get<2>(local_z_direction_control_planes_) -
													   std::get<1>(local_z_direction_control_planes_);

							const double new_plane_d_z = std::get<2>(local_z_direction_control_planes_) + gap_z_plane;

							const Vec3 shift_z = {0.0, 0.0, -new_plane_d_z};

							// find intersection edge between y and z triangle
							const std::vector<CMap2::Vertex> vertex_face_x;
							vertex_face_x.push_back(y_triangle1.vertices[2]);
							vertex_face_x.push_back(y_triangle1.vertices[0]);
							vertex_face_x.push_back(y_triangle2.vertices[1]);
							vertex_face_x.push_back(y_triangle1.vertices[1]);

							const std::vector<CMap2::Vertex> vertex_face_z;
							vertex_face_z.push_back(z_triangle1.vertices[2]);
							vertex_face_z.push_back(z_triangle1.vertices[0]);
							vertex_face_z.push_back(z_triangle2.vertices[1]);
							vertex_face_z.push_back(z_triangle1.vertices[1]);

							std::vector<CMap::Vertex> intersect_vertices;
							// TODO: find alternative to double loops
							for (std::size_t i = 0; i < 4; i++)
							{
								CMap2::Vertex target_x_vertex = vertex_face_x[i];
								for (std::size_t j = 0; j < 4; j++)
								{
									CMap2::Vertex target_z_vertex = vertex_face_z[j];

									if (target_x_vertex == target_z_vertex)
									{
										intersect_vertices.push_back(target_x_vertex);
									}
								}
							}

							// create new face composed of the intersecting edge and this edge shifted by y
							const Vec3 face_position0 =
								value<Vec3>(*control_cage_, control_cage_vertex_position_, intersect_vertices[0]);
							const Vec3 face_position1 =
								value<Vec3>(*control_cage_, control_cage_vertex_position_, intersect_vertices[1]);

							const Vec3 face_position2 = face_position1 + shift_x;
							const Vec3 face_position3 = face_position0 + shift_x;

							const Triangle local_triangle1, triangle2;
							local_triangle1.positions = {face_position1, face_position3, face_position0};
							local_triangle2.positions = {face_position1, face_position2, face_position3};

							std::vector<Triangle> virtual_cube_triangles =
								get_virtual_cube_triangles(local_triangle1, local_triangle2, shift_z);
						}
						else
						{
							// fakeCubeDiag[1]
							const Triangle z_triangle1 = control_cage_face_d_z_min_.first;
							const Triangle z_triangle2 = control_cage_face_d_z_min_.second;

							const double gap_z_plane = std::get<2>(local_z_direction_control_planes_) -
													   std::get<1>(local_z_direction_control_planes_);

							const double new_plane_d_z = std::get<1>(local_z_direction_control_planes_) - gap_z_plane;

							const Vec3 shift_z = {0.0, 0.0, -new_plane_d_z};

							// find intersection edge between y and z triangle
							const std::vector<CMap2::Vertex> vertex_face_x;
							vertex_face_x.push_back(y_triangle1.vertices[2]);
							vertex_face_x.push_back(y_triangle1.vertices[0]);
							vertex_face_x.push_back(y_triangle2.vertices[1]);
							vertex_face_x.push_back(y_triangle1.vertices[1]);

							const std::vector<CMap2::Vertex> vertex_face_z;
							vertex_face_z.push_back(z_triangle1.vertices[2]);
							vertex_face_z.push_back(z_triangle1.vertices[0]);
							vertex_face_z.push_back(z_triangle2.vertices[1]);
							vertex_face_z.push_back(z_triangle1.vertices[1]);

							std::vector<CMap::Vertex> intersect_vertices;
							// TODO: find alternative to double loops
							for (std::size_t i = 0; i < 4; i++)
							{
								CMap2::Vertex target_x_vertex = vertex_face_y[i];
								for (std::size_t j = 0; j < 4; j++)
								{
									CMap2::Vertex target_z_vertex = vertex_face_z[j];

									if (target_x_vertex == target_z_vertex)
									{
										intersect_vertices.push_back(target_x_vertex);
									}
								}
							}

							// create new face composed of the intersecting edge and this edge shifted by y
							const Vec3 face_position0 =
								value<Vec3>(*control_cage_, control_cage_vertex_position_, intersect_vertices[0]);
							const Vec3 face_position1 =
								value<Vec3>(*control_cage_, control_cage_vertex_position_, intersect_vertices[1]);

							const Vec3 face_position2 = face_position1 + shift_x;
							const Vec3 face_position3 = face_position0 + shift_x;

							const Triangle local_triangle1, triangle2;
							local_triangle1.positions = {face_position1, face_position3, face_position0};
							local_triangle2.positions = {face_position1, face_position2, face_position3};

							std::vector<Triangle> virtual_cube_triangles =
								get_virtual_cube_triangles(local_triangle1, local_triangle2, shift_z);
						}
					}
					else
					{
						const Triangle x_triangle1 = control_cage_face_d_x_min_.first;
						const Triangle x_triangle2 = control_cage_face_d_x_min_.second;

						const double gap_x_plane = std::get<2>(local_x_direction_control_planes_) -
												   std::get<1>(local_x_direction_control_planes_);

						const double new_plane_d_x = std::get<1>(local_x_direction_control_planes_) - gap_x_plane;

						const Vec3 shift_x = {-new_plane_d_y, 0.0, 0.0};
						if (d_z > std::get<2>(local_z_direction_control_planes_))
						{
							// fakeCubeDiag[2]
							const Triangle z_triangle1 = control_cage_face_d_z_max_.first;
							const Triangle z_triangle2 = control_cage_face_d_z_max_.second;

							// find plane of "fake" cube
							const double gap_z_plane = std::get<2>(local_z_direction_control_planes_) -
													   std::get<1>(local_z_direction_control_planes_);

							const double new_plane_d_z = std::get<2>(local_z_direction_control_planes_) + gap_z_plane;

							const Vec3 shift_z = {0.0, 0.0, -new_plane_d_z};

							// find intersection edge between y and z triangle
							const std::vector<CMap2::Vertex> vertex_face_x;
							vertex_face_x.push_back(y_triangle1.vertices[2]);
							vertex_face_x.push_back(y_triangle1.vertices[0]);
							vertex_face_x.push_back(y_triangle2.vertices[1]);
							vertex_face_x.push_back(y_triangle1.vertices[1]);

							const std::vector<CMap2::Vertex> vertex_face_z;
							vertex_face_z.push_back(z_triangle1.vertices[2]);
							vertex_face_z.push_back(z_triangle1.vertices[0]);
							vertex_face_z.push_back(z_triangle2.vertices[1]);
							vertex_face_z.push_back(z_triangle1.vertices[1]);

							std::vector<CMap::Vertex> intersect_vertices;
							// TODO: find alternative to double loops
							for (std::size_t i = 0; i < 4; i++)
							{
								CMap2::Vertex target_x_vertex = vertex_face_x[i];
								for (std::size_t j = 0; j < 4; j++)
								{
									CMap2::Vertex target_z_vertex = vertex_face_z[j];

									if (target_x_vertex == target_z_vertex)
									{
										intersect_vertices.push_back(target_x_vertex);
									}
								}
							}

							// create new face composed of the intersecting edge and this edge shifted by y
							const Vec3 face_position0 =
								value<Vec3>(*control_cage_, control_cage_vertex_position_, intersect_vertices[0]);
							const Vec3 face_position1 =
								value<Vec3>(*control_cage_, control_cage_vertex_position_, intersect_vertices[1]);

							const Vec3 face_position2 = face_position1 + shift_x;
							const Vec3 face_position3 = face_position0 + shift_x;

							const Triangle local_triangle1, triangle2;
							local_triangle1.positions = {face_position1, face_position3, face_position0};
							local_triangle2.positions = {face_position1, face_position2, face_position3};

							std::vector<Triangle> virtual_cube_triangles =
								get_virtual_cube_triangles(local_triangle1, local_triangle2, shift_z);
						}
						else
						{
							// fakeCubeDiag[0]
							const Triangle z_triangle1 = control_cage_face_d_z_min_.first;
							const Triangle z_triangle2 = control_cage_face_d_z_min_.second;

							const double gap_z_plane = std::get<2>(local_z_direction_control_planes_) -
													   std::get<1>(local_z_direction_control_planes_);

							const double new_plane_d_z = std::get<1>(local_z_direction_control_planes_) - gap_z_plane;

							const Vec3 shift_z = {0.0, 0.0, -new_plane_d_z};

							// find intersection edge between y and z triangle
							const std::vector<CMap2::Vertex> vertex_face_x;
							vertex_face_x.push_back(x_triangle1.vertices[2]);
							vertex_face_x.push_back(x_triangle1.vertices[0]);
							vertex_face_x.push_back(x_triangle2.vertices[1]);
							vertex_face_x.push_back(x_triangle1.vertices[1]);

							const std::vector<CMap2::Vertex> vertex_face_z;
							vertex_face_z.push_back(z_triangle1.vertices[2]);
							vertex_face_z.push_back(z_triangle1.vertices[0]);
							vertex_face_z.push_back(z_triangle2.vertices[1]);
							vertex_face_z.push_back(z_triangle1.vertices[1]);

							std::vector<CMap::Vertex> intersect_vertices;
							// TODO: find alternative to double loops
							for (std::size_t i = 0; i < 4; i++)
							{
								CMap2::Vertex target_x_vertex = vertex_face_x[i];
								for (std::size_t j = 0; j < 4; j++)
								{
									CMap2::Vertex target_z_vertex = vertex_face_z[j];

									if (target_x_vertex == target_z_vertex)
									{
										intersect_vertices.push_back(target_x_vertex);
									}
								}
							}

							// create new face composed of the intersecting edge and this edge shifted by y
							const Vec3 face_position0 =
								value<Vec3>(*control_cage_, control_cage_vertex_position_, intersect_vertices[0]);
							const Vec3 face_position1 =
								value<Vec3>(*control_cage_, control_cage_vertex_position_, intersect_vertices[1]);

							const Vec3 face_position2 = face_position1 + shift_x;
							const Vec3 face_position3 = face_position0 + shift_x;

							const Triangle local_triangle1, triangle2;
							local_triangle1.positions = {face_position1, face_position3, face_position0};
							local_triangle2.positions = {face_position1, face_position2, face_position3};

							std::vector<Triangle> virtual_cube_triangles =
								get_virtual_cube_triangles(local_triangle1, local_triangle2, shift_z);
						}
					}
				}
			}
			else if (valid_z_dir)
			{
				if (d_y > std::get<2>(local_y_direction_control_planes_))
				{
					const Triangle y_triangle1 = control_cage_face_d_y_max_.first;
					const Triangle y_triangle2 = control_cage_face_d_y_max_.second;

					const double gap_y_plane =
						std::get<2>(local_y_direction_control_planes_) - std::get<1>(local_y_direction_control_planes_);

					const double new_plane_d_y = std::get<2>(local_y_direction_control_planes_) + gap_y_plane;

					const Vec3 shift_y = {0.0, -new_plane_d_y, 0.0};
					if (d_x > std::get<2>(local_x_direction_control_planes_))
					{
						// fakeCubeDiag[7]
						const Triangle x_triangle1 = control_cage_face_d_x_max_.first;
						const Triangle x_triangle2 = control_cage_face_d_x_max_.second;

						// find plane of "fake" cube
						const double gap_x_plane = std::get<2>(local_x_direction_control_planes_) -
												   std::get<1>(local_x_direction_control_planes_);

						const double new_plane_d_x = std::get<2>(local_x_direction_control_planes_) + gap_x_plane;

						const Vec3 shift_x = {-new_plane_d_x, 0.0, 0.0};

						// find intersection edge between y and x triangle
						const std::vector<CMap2::Vertex> vertex_face_y;
						vertex_face_y.push_back(y_triangle1.vertices[2]);
						vertex_face_y.push_back(y_triangle1.vertices[0]);
						vertex_face_y.push_back(y_triangle2.vertices[1]);
						vertex_face_y.push_back(y_triangle1.vertices[1]);

						const std::vector<CMap2::Vertex> vertex_face_x;
						vertex_face_x.push_back(x_triangle1.vertices[2]);
						vertex_face_x.push_back(x_triangle1.vertices[0]);
						vertex_face_x.push_back(x_triangle2.vertices[1]);
						vertex_face_x.push_back(x_triangle1.vertices[1]);

						std::vector<CMap::Vertex> intersect_vertices;
						// TODO: find alternative to double loops
						for (std::size_t i = 0; i < 4; i++)
						{
							CMap2::Vertex target_y_vertex = vertex_face_y[i];
							for (std::size_t j = 0; j < 4; j++)
							{
								CMap2::Vertex target_y_vertex = vertex_face_y[j];

								if (target_x_vertex == target_z_vertex)
								{
									intersect_vertices.push_back(target_x_vertex);
								}
							}
						}

						// create new face composed of the intersecting edge and this edge shifted by y
						const Vec3 face_position0 =
							value<Vec3>(*control_cage_, control_cage_vertex_position_, intersect_vertices[0]);
						const Vec3 face_position1 =
							value<Vec3>(*control_cage_, control_cage_vertex_position_, intersect_vertices[1]);

						const Vec3 face_position2 = face_position1 + shift_y;
						const Vec3 face_position3 = face_position0 + shift_y;

						const Triangle local_triangle1, triangle2;
						local_triangle1.positions = {face_position1, face_position3, face_position0};
						local_triangle2.positions = {face_position1, face_position2, face_position3};

						std::vector<Triangle> virtual_cube_triangles =
							get_virtual_cube_triangles(local_triangle1, local_triangle2, shift_x);
					}
					else
					{
						// fakeCubeDiag[5]
						const Triangle x_triangle1 = control_cage_face_d_x_min_.first;
						const Triangle x_triangle2 = control_cage_face_d_x_min_.second;

						const double gap_x_plane = std::get<2>(local_x_direction_control_planes_) -
													   std::get<1>(local_x_direction_control_planes_);

							const double new_plane_d_x = std::get<1>(local_x_direction_control_planes_) - gap_x_plane;

							const Vec3 shift_x = {-new_plane_d_z, 0.0, 0.0};

							// find intersection edge between y and z triangle
							const std::vector<CMap2::Vertex> vertex_face_y;
							vertex_face_y.push_back(y_triangle1.vertices[2]);
							vertex_face_y.push_back(y_triangle1.vertices[0]);
							vertex_face_y.push_back(y_triangle2.vertices[1]);
							vertex_face_y.push_back(y_triangle1.vertices[1]);

							const std::vector<CMap2::Vertex> vertex_face_x;
							vertex_face_x.push_back(x_triangle1.vertices[2]);
							vertex_face_x.push_back(x_triangle1.vertices[0]);
							vertex_face_x.push_back(x_triangle2.vertices[1]);
							vertex_face_x.push_back(x_triangle1.vertices[1]);

							std::vector<CMap::Vertex> intersect_vertices;
							// TODO: find alternative to double loops
							for (std::size_t i = 0; i < 4; i++)
							{
								CMap2::Vertex target_y_vertex = vertex_face_y[i];
								for (std::size_t j = 0; j < 4; j++)
								{
									CMap2::Vertex target_x_vertex = vertex_face_x[j];

									if (target_y_vertex == target_x_vertex)
									{
										intersect_vertices.push_back(target_x_vertex);
									}
								}
							}

							// create new face composed of the intersecting edge and this edge shifted by y
							const Vec3 face_position0 =
								value<Vec3>(*control_cage_, control_cage_vertex_position_, intersect_vertices[0]);
							const Vec3 face_position1 =
								value<Vec3>(*control_cage_, control_cage_vertex_position_, intersect_vertices[1]);

							const Vec3 face_position2 = face_position1 + shift_y;
							const Vec3 face_position3 = face_position0 + shift_y;

							const Triangle local_triangle1, triangle2;
							local_triangle1.positions = {face_position1, face_position3, face_position0};
							local_triangle2.positions = {face_position1, face_position2, face_position3};

							std::vector<Triangle> virtual_cube_triangles =
								get_virtual_cube_triangles(local_triangle1, local_triangle2, shift_x);

					}
				}
				else
				{
					const Triangle y_triangle1 = control_cage_face_d_y_min_.first;
						const Triangle y_triangle2 = control_cage_face_d_y_min_.second;

						const double gap_y_plane = std::get<2>(local_y_direction_control_planes_) -
												   std::get<1>(local_y_direction_control_planes_);

						const double new_plane_d_y = std::get<1>(local_y_direction_control_planes_) - gap_y_plane;

						const Vec3 shift_y = {0.0, -new_plane_d_y, 0.0};
					if (d_x > std::get<2>(local_x_direction_control_planes_))
					{
						// fakeCubeDiag[6]
						const Triangle x_triangle1 = control_cage_face_d_x_max_.first;
							const Triangle x_triangle2 = control_cage_face_d_x_max_.second;

							// find plane of "fake" cube
							const double gap_x_plane = std::get<2>(local_x_direction_control_planes_) -
													   std::get<1>(local_x_direction_control_planes_);

							const double new_plane_d_x = std::get<2>(local_x_direction_control_planes_) + gap_x_plane;

							const Vec3 shift_x = {-new_plane_d_x, 0.0, 0.0};

							// find intersection edge between y and z triangle
							const std::vector<CMap2::Vertex> vertex_face_y;
							vertex_face_y.push_back(y_triangle1.vertices[2]);
							vertex_face_y.push_back(y_triangle1.vertices[0]);
							vertex_face_y.push_back(y_triangle2.vertices[1]);
							vertex_face_y.push_back(y_triangle1.vertices[1]);

							const std::vector<CMap2::Vertex> vertex_face_x;
							vertex_face_x.push_back(z_triangle1.vertices[2]);
							vertex_face_x.push_back(z_triangle1.vertices[0]);
							vertex_face_x.push_back(z_triangle2.vertices[1]);
							vertex_face_x.push_back(z_triangle1.vertices[1]);

							std::vector<CMap::Vertex> intersect_vertices;
							// TODO: find alternative to double loops
							for (std::size_t i = 0; i < 4; i++)
							{
								CMap2::Vertex target_y_vertex = vertex_face_y[i];
								for (std::size_t j = 0; j < 4; j++)
								{
									CMap2::Vertex target_x_vertex = vertex_face_x[j];

									if (target_y_vertex == target_x_vertex)
									{
										intersect_vertices.push_back(target_y_vertex);
									}
								}
							}

							// create new face composed of the intersecting edge and this edge shifted by y
							const Vec3 face_position0 =
								value<Vec3>(*control_cage_, control_cage_vertex_position_, intersect_vertices[0]);
							const Vec3 face_position1 =
								value<Vec3>(*control_cage_, control_cage_vertex_position_, intersect_vertices[1]);

							const Vec3 face_position2 = face_position1 + shift_y;
							const Vec3 face_position3 = face_position0 + shift_y;

							const Triangle local_triangle1, triangle2;
							local_triangle1.positions = {face_position1, face_position3, face_position0};
							local_triangle2.positions = {face_position1, face_position2, face_position3};

							std::vector<Triangle> virtual_cube_triangles =
								get_virtual_cube_triangles(local_triangle1, local_triangle2, shift_x);
					}
					else
					{
						// fakeCubeDiag[4]
						const Triangle x_triangle1 = control_cage_face_d_x_min_.first;
							const Triangle x_triangle2 = control_cage_face_d_x_min_.second;

							const double gap_x_plane = std::get<2>(local_x_direction_control_planes_) -
													   std::get<1>(local_x_direction_control_planes_);

							const double new_plane_d_x = std::get<1>(local_x_direction_control_planes_) - gap_x_plane;

							const Vec3 shift_x = {-new_plane_d_x, 0.0, 0.0};

							// find intersection edge between y and z triangle
							const std::vector<CMap2::Vertex> vertex_face_y;
							vertex_face_y.push_back(y_triangle1.vertices[2]);
							vertex_face_y.push_back(y_triangle1.vertices[0]);
							vertex_face_y.push_back(y_triangle2.vertices[1]);
							vertex_face_y.push_back(y_triangle1.vertices[1]);

							const std::vector<CMap2::Vertex> vertex_face_x;
							vertex_face_x.push_back(x_triangle1.vertices[2]);
							vertex_face_x.push_back(x_triangle1.vertices[0]);
							vertex_face_x.push_back(x_triangle2.vertices[1]);
							vertex_face_x.push_back(x_triangle1.vertices[1]);

							std::vector<CMap::Vertex> intersect_vertices;
							// TODO: find alternative to double loops
							for (std::size_t i = 0; i < 4; i++)
							{
								CMap2::Vertex target_y_vertex = vertex_face_y[i];
								for (std::size_t j = 0; j < 4; j++)
								{
									CMap2::Vertex target_x_vertex = vertex_face_x[j];

									if (target_y_vertex == target_x_vertex)
									{
										intersect_vertices.push_back(target_x_vertex);
									}
								}
							}

							// create new face composed of the intersecting edge and this edge shifted by y
							const Vec3 face_position0 =
								value<Vec3>(*control_cage_, control_cage_vertex_position_, intersect_vertices[0]);
							const Vec3 face_position1 =
								value<Vec3>(*control_cage_, control_cage_vertex_position_, intersect_vertices[1]);

							const Vec3 face_position2 = face_position1 + shift_y;
							const Vec3 face_position3 = face_position0 + shift_y;

							const Triangle local_triangle1, triangle2;
							local_triangle1.positions = {face_position1, face_position3, face_position0};
							local_triangle2.positions = {face_position1, face_position2, face_position3};

							std::vector<Triangle> virtual_cube_triangles =
								get_virtual_cube_triangles(local_triangle1, local_triangle2, shift_x);
					}
				}
			}
			else
			{
				if (d_x > std::get<2>(local_x_direction_control_planes_))
				{
					if (d_y > std::get<2>(local_y_direction_control_planes_))
					{
						if (d_z > std::get<2>(local_z_direction_control_planes_))
						{
							// fakeCubeOut[7]
						}
						else
						{
							// fakeCubeOut[5]
						}
					}
					else
					{
						if (d_z > std::get<2>(local_z_direction_control_planes_))
						{
							// fakeCubeOut[3]
						}
						else
						{
							// fakeCubeOut[1]
						}
					}
				}
				else
				{
					if (d_y > std::get<2>(local_y_direction_control_planes_))
					{
						if (d_z > std::get<2>(local_z_direction_control_planes_))
						{
							// fakeCubeOut[6]
						}
						else
						{
							// fakeCubeOut[4]
						}
					}
					else
					{
						if (d_z > std::get<2>(local_z_direction_control_planes_))
						{
							// fakeCubeOut[2]
						}
						else
						{
							// fakeCubeOut[0]
						}
					}
				}
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

	std::pair<Triangle, Triangle> control_cage_face_d_x_min_; // pair of indices in cage_triangles_
	std::pair<Triangle, Triangle> control_cage_face_d_x_max_;

	std::pair<Triangle, Triangle> control_cage_face_d_y_min_;
	std::pair<Triangle, Triangle> control_cage_face_d_y_max_;

	std::pair<Triangle, Triangle> control_cage_face_d_z_min_;
	std::pair<Triangle, Triangle> control_cage_face_d_z_max_;

	// local frame direction, mininum and maximum positions of the planes from control cage (d from ax+by+cz+d = 0)
	std::tuple<Vec3, double, double> local_x_direction_control_planes_;
	std::tuple<Vec3, double, double> local_y_direction_control_planes_;
	std::tuple<Vec3, double, double> local_z_direction_control_planes_;

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
					control_cage_face_d_x_min_ = std::make_pair(triangle1, triangle2);
				}

				if (d > d_x_max)
				{
					d_x_max = d;
					control_cage_face_d_x_max_ = std::make_pair(triangle1, triangle2);
				}
			}
			else if (normal.dot(y_dir) != 0.0)
			{
				const double d = -(triangle_point.dot(y_dir));
				if (d < d_y_min)
				{
					d_y_min = d;
					control_cage_face_d_y_min_ = std::make_pair(triangle1, triangle2);
				}

				if (d > d_y_max)
				{
					d_y_max = d;
					control_cage_face_d_y_min_ = std::make_pair(triangle1, triangle2);
				}
			}
			else
			{
				const double d = -(triangle_point.dot(z_dir));
				if (d < d_z_min)
				{
					d_z_min = d;
					control_cage_face_d_z_min_ = std::make_pair(triangle1, triangle2);
				}

				if (d > d_z_max)
				{
					d_z_max = d;
					control_cage_face_d_z_min_ = std::make_pair(triangle1, triangle2);
				}
			}
		}

		local_x_direction_control_planes_ = std::make_tuple(x_dir, d_x_min, d_x_max);

		local_y_direction_control_planes_ = std::make_tuple(y_dir, d_y_min, d_y_max);

		local_z_direction_control_planes_ = std::make_tuple(z_dir, d_z_min, d_z_max);
	}

	std::vector<Triangle> get_virtual_cube_triangles(const Triangle& triangle1, const Triangle& triangle2,
													 const Vec3& shift_vector)
	{
		std::vector<Triangle> virtual_cube_triangles(6);
		virtual_cube_triangles[0] = triangle1;
		virtual_cube_triangles[1] = triangle2;

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
		// triangle3.normal = -1.0 * triangle1.normal;

		triangle4.positions = {new_positions.begin() + 3, new_positions.end()};
		// triangle4.normal = -1.0 * triangle2.normal;

		triangle5.positions = {triangle1.positions[0], triangle3.positions[2], triangle1.positions[2]};

		triangle6.positions = {triangle1.positions[0], triangle3.positions[0], triangle3.positions[2]};

		triangle7.positions = {triangle1.positions[2], triangle3.positions[1], triangle3.positions[2]};

		triangle8.positions = {triangle1.positions[2], triangle1.positions[1], triangle3.positions[1]};

		triangle9.positions = {triangle1.positions[0], triangle4.positions[1], triangle3.positions[0]};

		triangle10.positions = {triangle1.positions[0], triangle2.positions[1], triangle4.positions[1]};

		triangle11.positions = {triangle4.positions[1], triangle1.positions[1], triangle3.positions[1]};

		triangle12.positions = {triangle4.positions[1], triangle2.positions[1], triangle1.positions[1]};

		virtual_cube_triangles[3] = triangle3;
		virtual_cube_triangles[4] = triangle4;
		virtual_cube_triangles[5] = triangle5;
		virtual_cube_triangles[6] = triangle6;
		virtual_cube_triangles[7] = triangle7;
		virtual_cube_triangles[8] = triangle8;
		virtual_cube_triangles[9] = triangle9;
		virtual_cube_triangles[10] = triangle10;
		virtual_cube_triangles[11] = triangle11;
		virtual_cube_triangles[12] = triangle12;

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
