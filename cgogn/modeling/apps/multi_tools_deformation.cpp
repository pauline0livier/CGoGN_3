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

#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/core/ui_modules/graph_provider.h>

#include <cgogn/rendering/ui_modules/surface_render.h>
#include <cgogn/rendering/ui_modules/graph_render.h>

#include <cgogn/modeling/ui_modules/multi_tools_deformation_module.h>

#include <cgogn/modeling/algos/deformation/creation_space_tool.h>

#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_DATA_PATH) "/meshes/"

using namespace cgogn::numerics;

using Mesh = cgogn::CMap2;
using Graph = cgogn::IncidenceGraph;

template <typename T>
using GraphAttribute = typename cgogn::mesh_traits<Graph>::Attribute<T>;

template <typename T>
using MeshAttribute = typename cgogn::mesh_traits<Mesh>::Attribute<T>;
using MeshVertex = typename cgogn::mesh_traits<Mesh>::Vertex;
using MeshFace = typename cgogn::mesh_traits<Mesh>::Face;

int main(int argc, char** argv)
{
	using Vec3 = cgogn::geometry::Vec3;
	using Scalar = cgogn::geometry::Scalar;

	std::string filename;
	 
	if (argc < 2)
		filename = std::string(DEFAULT_MESH_PATH) + std::string("obj/low-poly-fox-by-pixelmannen.obj"); 
	else
		filename = std::string(argv[1]);

	cgogn::thread_start();

	cgogn::ui::App app;
	app.set_window_title("Deformation Multi tools");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Mesh> mp(app);
	cgogn::ui::GraphProvider<Graph> gp(app);

	cgogn::ui::SurfaceRender<Mesh> sr(app);
	cgogn::ui::GraphRender<Graph> gr(app);

	cgogn::ui::MultiToolsDeformation<Mesh, Graph> sd2(app); 

	app.init_modules();

	cgogn::ui::View* v1 = app.current_view();

	v1->link_module(&mp);

	v1->link_module(&gp);

	v1->link_module(&gr);
	v1->link_module(&sr);

	v1->link_module(&sd2); 

	Mesh* m = mp.load_surface_from_file(filename);
	if (!m)
	{
		std::cout << "File could not be loaded" << std::endl;
		return 1;
	}

	auto vertex_position = cgogn::get_attribute<Vec3, cgogn::mesh_traits<Mesh>::Vertex>(*m, "position");

	auto vertex_normal = cgogn::add_attribute<Vec3, cgogn::mesh_traits<Mesh>::Vertex>(*m, "normal");

	auto object_position_indices = cgogn::add_attribute<uint32, cgogn::mesh_traits<Mesh>::Vertex>(*m, "vertex_index");
	cgogn::modeling::set_attribute_vertex_index(*m, object_position_indices.get()); 

	sd2.set_model(*m, vertex_position, vertex_normal);

	sr.set_vertex_position(*v1, *m, vertex_position);
	sr.set_vertex_normal(*v1, *m, vertex_normal);
	sr.set_render_vertices(*v1, *m, false);
	sr.set_render_edges(*v1, *m, false);

	return app.launch();
}
