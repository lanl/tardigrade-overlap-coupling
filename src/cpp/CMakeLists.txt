# Project

# Build the micromorphic balance equations
add_library (${BALANCE_EQUATION_LIBRARY} SHARED "${BALANCE_EQUATION_LIBRARY_FILENAME}.cpp" "${BALANCE_EQUATION_LIBRARY_FILENAME}.h")
set_target_properties(${BALANCE_EQUATION_LIBRARY} PROPERTIES CXX_STANDARD 11
                                                 PUBLIC_HEADER ${BALANCE_EQUATION_LIBRARY_FILENAME}.h)
target_link_libraries(${BALANCE_EQUATION_LIBRARY} error_tools)
target_compile_options(${BALANCE_EQUATION_LIBRARY} PUBLIC)

# Local builds of upstream projects require local include paths
if(NOT cmake_build_type_lower STREQUAL "release")
    target_include_directories(${BALANCE_EQUATION_LIBRARY} PUBLIC
                               "${vector_tools_SOURCE_DIR}/${CPP_SRC_PATH}"
                               "${error_tools_SOURCE_DIR}/${CPP_SRC_PATH}")
endif()

# Build the micromorphic material model library
set(LOCAL_BUILD_INCLUDES "")
foreach(user_subroutine ${USER_SUBROUTINES})

    add_library("${user_subroutine}_interface" SHARED "${user_subroutine}_interface.h" "${user_subroutine}_interface.cpp")
    set_target_properties("${user_subroutine}_interface" PROPERTIES CXX_STANDARD 11
                                                         PUBLIC_HEADER "${user_subroutine}_interface.h")
#                                                         PREFIX ""
#                                                         SUFFIX ".o")
    target_compile_options("${user_subroutine}_interface" PUBLIC)
    target_link_libraries("${user_subroutine}_interface" -Wl,--whole-archive "${user_subroutine}" -Wl,--no-whole-archive)
    if(NOT cmake_build_type_lower STREQUAL "release")
        set(LOCAL_BUILD_INCLUDES ${LOCAL_BUILD_INCLUDES} "${${user_subroutine}_SOURCE_DIR}/${CPP_SRC_PATH}" )
        target_include_directories("${user_subroutine}_interface" PUBLIC "${${user_subroutine}_SOURCE_DIR}/${CPP_SRC_PATH}")
    endif()
endforeach(user_subroutine)

# Project
add_library (${MATERIAL_MODEL_LIBRARY} SHARED "${MATERIAL_MODEL_LIBRARY_FILENAME}.cpp" "${MATERIAL_MODEL_LIBRARY_FILENAME}.h")
set_target_properties(${MATERIAL_MODEL_LIBRARY} PROPERTIES CXX_STANDARD 11
                                                PUBLIC_HEADER "${MATERIAL_MODEL_LIBRARY_FILENAMEE}.h"
                                                SUFFIX ".so")
target_compile_options(${MATERIAL_MODEL_LIBRARY} PUBLIC)
target_link_libraries(${MATERIAL_MODEL_LIBRARY} error_tools ${USER_INTERFACES})

# Local builds of upstream projects require local include paths
if(NOT cmake_build_type_lower STREQUAL "release")
    target_include_directories(${MATERIAL_MODEL_LIBRARY} PUBLIC
                               "${vector_tools_SOURCE_DIR}/${CPP_SRC_PATH}"
                               "${error_tools_SOURCE_DIR}/${CPP_SRC_PATH}"
                               ${LOCAL_BUILD_INCLUDES})
endif()

foreach(package ${BALANCE_EQUATION_LIBRARY} ${MATERIAL_MODEL_LIBRARY})
    install(TARGETS ${package}
            EXPORT ${package}_Targets
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
            PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
endforeach(package)
