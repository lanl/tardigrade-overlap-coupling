# Project

# Build the support modules
set(PREVIOUS_MODULES "")
foreach(support_module ${LOCAL_SUPPORT_MODULES})
    add_library (${support_module} SHARED "${support_module}.cpp" "${support_module}.h")
    set_target_properties(${support_module} PROPERTIES CXX_STANDARD 11
                                                     PUBLIC_HEADER ${support_module}.h)
    target_link_libraries(${support_module} error_tools micromorphic_tools yaml-cpp Xdmf XdmfCore XdmfUtils xml2 ${PREVIOUS_MODULES})
    target_compile_options(${support_module} PUBLIC)
    set(PREVIOUS_MODULES ${PREVIOUS_MODULES} ${support_module})
endforeach(support_module)

## Build the main library
#add_library (${PROJECT_NAME} SHARED "${PROJECT_NAME}.cpp" "${PROJECT_NAME}.h")
#set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 11
#                                                 PUBLIC_HEADER ${PROJECT_NAME}.h)
#target_link_libraries(${PROJECT_NAME} error_tools solver_tools micromorphic_tools ${LOCAL_SUPPORT_MODULES})
#target_compile_options(${PROJECT_NAME} PUBLIC)
#
## Local builds of upstream projects require local include paths
#if(NOT cmake_build_type_lower STREQUAL "release")
#    target_include_directories(${PROJECT_NAME} PUBLIC
#                               "${vector_tools_SOURCE_DIR}/${CPP_SRC_PATH}"
#                               "${error_tools_SOURCE_DIR}/${CPP_SRC_PATH}"
#                               "${micromorphic_tools_SOURCE_DIR}/${CPP_SRC_PATH}")
#endif()
#
#
#foreach(package ${PROJECT_NAME} ${FILTER_NAME})
#    install(TARGETS ${package}
#            EXPORT ${package}_Targets
#            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
#            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
#            INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
#            PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
#endforeach(package)
