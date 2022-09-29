cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(cgogn_rendering_examples
    LANGUAGES CXX
)

find_package(cgogn_core REQUIRED)
find_package(cgogn_ui REQUIRED)
find_package(cgogn_io REQUIRED)
find_package(cgogn_rendering REQUIRED)
find_package(cgogn_modeling REQUIRED)

add_executable(modeling modeling.cpp)
target_link_libraries(modeling
	cgogn::core
	cgogn::ui
	cgogn::io
	cgogn::rendering
	cgogn::modeling
	${CMAKE_DL_LIBS}
)

add_executable(tubular_mesh tubular_mesh.cpp)
target_link_libraries(tubular_mesh
	cgogn::core
	cgogn::ui
	cgogn::io
	cgogn::rendering
	cgogn::modeling
	simpleICP
	${CMAKE_DL_LIBS}
)

add_executable(deformation deformation.cpp)
target_link_libraries(deformation
	cgogn::core
	cgogn::ui
	cgogn::io
	cgogn::rendering
	cgogn::modeling
	${CMAKE_DL_LIBS}
)

add_executable(def_global_cage_PO def_global_cage_PO.cpp)
target_link_libraries(def_global_cage_PO
	cgogn::core
	cgogn::ui
	cgogn::io
	cgogn::rendering
	cgogn::modeling
	${CMAKE_DL_LIBS}
)

add_executable(def_local_tool_PO def_local_tool_PO.cpp)
target_link_libraries(def_local_tool_PO
	cgogn::core
	cgogn::ui
	cgogn::io
	cgogn::rendering
	cgogn::modeling
	${CMAKE_DL_LIBS}
)

add_executable(def_local_tool_PO_update def_local_tool_PO_update.cpp)
target_link_libraries(def_local_tool_PO_update
	cgogn::core
	cgogn::ui
	cgogn::io
	cgogn::rendering
	cgogn::modeling
	${CMAKE_DL_LIBS}
)

add_executable(skeleton_extractor skeleton_extractor.cpp)
target_link_libraries(skeleton_extractor
	cgogn::core
	cgogn::ui
	cgogn::io
	cgogn::rendering
	cgogn::modeling
	${CMAKE_DL_LIBS}
)

add_executable(convex_hull convex_hull.cpp)
target_link_libraries(convex_hull
	cgogn::core
	cgogn::ui
	cgogn::io
	cgogn::rendering
	cgogn::modeling
	${CMAKE_DL_LIBS}
)

if(APPLE)
	find_library(CORE_FOUNDATION CoreFoundation)
	find_library(CARBON Carbon)
	target_link_libraries(modeling
		${CORE_FOUNDATION}
		${CARBON}
	)
	target_link_libraries(tubular_mesh
		${CORE_FOUNDATION}
		${CARBON}
	)
	target_link_libraries(deformation
		${CORE_FOUNDATION}
		${CARBON}
	)
	target_link_libraries(def_global_cage_PO
		${CORE_FOUNDATION}
		${CARBON}
	)
	target_link_libraries(def_local_tool_PO
		${CORE_FOUNDATION}
		${CARBON}
	)
	target_link_libraries(def_local_tool_PO_update
		${CORE_FOUNDATION}
		${CARBON}
	)
	target_link_libraries(skeleton_extractor
		${CORE_FOUNDATION}
		${CARBON}
	)
	target_link_libraries(convex_hull
		${CORE_FOUNDATION}
		${CARBON}
	)
endif()

set_target_properties(def_global_cage_PO PROPERTIES FOLDER pauline)
set_target_properties(def_local_tool_PO PROPERTIES FOLDER pauline)