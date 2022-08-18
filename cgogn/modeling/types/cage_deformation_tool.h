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
#include <cgogn/geometry/types/vector_traits.h>

namespace cgogn
{

namespace modeling
{

template <typename MESH>
class CageDeformationTool
{
	template <typename T>
	using Attribute = typename mesh_traits<MESH>::template Attribute<T>;
	using Vertex = typename mesh_traits<MESH>::Vertex;
	using Face = typename mesh_traits<MESH>::Face;

public:
	CageDeformationTool()
	{
		influence_cage_ = nullptr;
		influence_set_ = nullptr;
		control_set_ = nullptr;
		local_def = false;
		b_simple_mix = false;
		m_hFactor = -1.0f;
	}

	~CageDeformationTool()
	{
	}

	void bind_object_mvc(MESH& object, const std::shared_ptr<MeshAttribute<Vec3>>& object_vertex_position, MESH& cage,
						 const std::shared_ptr<MeshAttribute<Vec3>>& cage_vertex_position)
	{

		std::shared_ptr<Attribute<uint32>> object_vertex_index = get_attribute<uint32, Vertex>(object, "weight_index");

		std::shared_ptr<Attribute<uint32>> cage_vertex_index = get_attribute<uint32, Vertex>(cage, "weight_index");

		std::shared_ptr<Attribute<bool>> cage_vertex_marked = get_attribute<bool, Vertex>(cage, "marked_vertex");

		uint32 nbv_object = nb_cells<Vertex>(object);
		uint32 nbv_cage = nb_cells<Vertex>(cage);

		coords_.resize(nbv_object, nbv_cage);

		parallel_foreach_cell(object, [&](MeshVertex v) -> bool {
			const Vec3& surface_point = value<Vec3>(object, object_vertex_position, v);
			uint32 surface_point_idx = value<uint32>(object, object_vertex_index, v);

			DartMarker dm(cage);
			float sumMVC = 0.0;
			for (Dart d = cage.begin(), end = cage.end(); d != end; d = cage.next(d))
			{
				MeshVertex cage_vertex = CMap2::Vertex(d);
				bool vc_marked = value<bool>(cage, cage_vertex_marked, cage_vertex);

				if (!dm.is_marked(d) && !vc_marked)
				{

					const Vec3& cage_point = value<Vec3>(cage, cage_vertex_position, cage_vertex);
					uint32 cage_point_idx = value<uint32>(cage, cage_vertex_index, cage_vertex);

					float mvc_value = compute_mvc(surface_point, d, cage, cage_point, cage_vertex_position.get());

					cd.coords_(surface_point_idx, cage_point_idx) = mvc_value;
					dm.mark(d);

					value<bool>(cage, cage_vertex_marked, cage_vertex) = true;

					sumMVC += mvc_value;
				}
			}

			float sum_lambda = 0.0;

			parallel_foreach_cell(cage, [&](MeshVertex vc) -> bool {
				uint32 cage_point_idx2 = value<uint32>(cage, cage_vertex_index, vc);

				cd.coords_(surface_point_idx, cage_point_idx2) =
					cd.coords_(surface_point_idx, cage_point_idx2) / sumMVC;

				value<bool>(cage, cage_vertex_marked, vc) = false;

				return true;
			});

			return true;
		});

		cd.cage_attribute_update_connection_ =
			boost::synapse::connect<typename MeshProvider<MESH>::template attribute_changed_t<Vec3>>(
				&cage, [&](MeshAttribute<Vec3>* attribute) {
					if (cd.cage_vertex_position_.get() == attribute)
					{

						std::shared_ptr<MeshAttribute<uint32>> object_vertex_index =
							cgogn::get_attribute<uint32, MeshVertex>(object, "weight_index");
						std::shared_ptr<MeshAttribute<uint32>> cage_vertex_index =
							cgogn::get_attribute<uint32, MeshVertex>(cage, "weight_index");

						parallel_foreach_cell(object, [&](MeshVertex v) -> bool {
							uint32 vidx = value<uint32>(object, object_vertex_index, v);

							Vec3 new_pos_ = {0.0, 0.0, 0.0};

							foreach_cell(cage, [&](MeshVertex cv) -> bool {
								const Vec3& cage_point = value<Vec3>(cage, cage_vertex_position, cv);
								uint32 cage_point_idx = value<uint32>(cage, cage_vertex_index, cv);

								new_pos_ += cd.coords_(vidx, cage_point_idx) * cage_point;

								return true;
							});

							value<Vec3>(object, object_vertex_position, v) = new_pos_;
							return true;
						});

						mesh_provider_->emit_attribute_changed(object, object_vertex_position.get());
					}
				});
	}



private:

	Cage3D influence_cage_;

	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> coords_;
	Eigen::Matrix<Vec2, Eigen::Dynamic, Eigen::Dynamic> n_coords_;

	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> ctrl_cage_coords_;

	Eigen::VectorXd attenuation_;
	Eigen::VectorXd mixing_factor_;

	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> global_matrix_;

	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> global_matrix_comp_;

	float smoothing_factor_;

	float m_hFactor;

	bool local_def;

	bool b_simple_mix;

	std::shared_ptr<boost::synapse::connection> cage_attribute_update_connection_;
	
};

} // namespace modeling

} // namespace cgogn

#endif // CGOGN_MODELING_CAGE_DEFORMATION_TOOL_H_
