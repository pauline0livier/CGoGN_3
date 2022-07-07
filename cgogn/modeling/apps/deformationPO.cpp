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
#include <cgogn/geometry/ui_modules/surface_differential_properties.h>
#include <cgogn/geometry/ui_modules/surface_selection.h>
#include <cgogn/modeling/ui_modules/surface_deformation.h>
#include <cgogn/modeling/ui_modules/cage_deformation.h>
#include <cgogn/rendering/ui_modules/surface_render.h>

#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_DATA_PATH) "/meshes/"

using namespace cgogn::numerics;

using Mesh = cgogn::CMap2;

template <typename T>
using Attribute = typename cgogn::mesh_traits<Mesh>::Attribute<T>;
using Vertex = typename cgogn::mesh_traits<Mesh>::Vertex;
using Face = typename cgogn::mesh_traits<Mesh>::Face;

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
	app.set_window_title("DeformationPO");
	app.set_window_size(1000, 800);

	cgogn::ui::MeshProvider<Mesh> mp(app);
	cgogn::ui::SurfaceRender<Mesh> sr(app);
	cgogn::ui::SurfaceDifferentialProperties<Mesh> sdp(app);
	cgogn::ui::SurfaceDeformation<Mesh> sd(app);
	cgogn::ui::SurfaceSelection<Mesh> ss(app);
	cgogn::ui::CageDeformation<Mesh> cd(app);

	app.init_modules();

	cgogn::ui::View* v1 = app.current_view();
	v1->link_module(&mp);
	v1->link_module(&sr);
	v1->link_module(&sd);
	v1->link_module(&ss);

	Mesh* m = mp.load_surface_from_file(filename);
	if (!m)
	{
		std::cout << "File could not be loaded" << std::endl;
		return 1;
	}

	std::shared_ptr<Attribute<Vec3>> vertex_position = cgogn::get_attribute<Vec3, Vertex>(*m, "position");
	std::shared_ptr<Attribute<Vec3>> vertex_normal = cgogn::add_attribute<Vec3, Vertex>(*m, "normal");

	sdp.compute_normal(*m, vertex_position.get(), vertex_normal.get());

	sr.set_vertex_position(*v1, *m, vertex_position);
	sr.set_vertex_normal(*v1, *m, vertex_normal);
	sr.set_render_vertices(*v1, *m, false);
	sr.set_render_edges(*v1, *m, false);

	ss.set_vertex_position(*m, vertex_position);

	Mesh* cage = cd.generate_cage(*m, vertex_position);
	std::shared_ptr<Attribute<Vec3>> cage_vertex_position = cgogn::get_attribute<Vec3, Vertex>(*cage, "position");

	std::shared_ptr<Attribute<Vec3>> cage_vertex_normal = cgogn::add_attribute<Vec3, Vertex>(*cage, "normal");

	sdp.compute_normal(*cage, cage_vertex_position.get(), cage_vertex_normal.get());



	sr.set_vertex_position(*v1, *cage, cage_vertex_position);
	sr.set_render_faces(*v1, *cage, false);

	cgogn::ui::MeshData<Mesh>& md = mp.mesh_data(*cage);
	cgogn::ui::CellsSet<Mesh, Vertex>& set = md.add_cells_set<Vertex>();

	ss.set_selected_mesh(*cage);
	ss.set_vertex_position(*cage, cage_vertex_position);
	ss.set_selected_cells_set(*cage, &set);

	sd.set_selected_mesh(*cage);
	sd.set_vertex_position(*cage, cage_vertex_position);
	sd.set_selected_handle_vertices_set(*cage, &set);

	return app.launch();
}
