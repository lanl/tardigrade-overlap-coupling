#Set the tests for the support modules

#Copy over all of the supporting YAML files
file(GLOB YAML_TEST_FILES
    "./*.yaml"
)
file(COPY ${YAML_TEST_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(GLOB XDMF_TEST_FILES
    "./*.xdmf"
)
file(COPY ${XDMF_TEST_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(GLOB H5_TEST_FILES
    "./*.h5"
)
file(COPY ${XDMF_TEST_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(GLOB TXT_TEST_FILES
    "./*.txt"
)
file(COPY ${TXT_TEST_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY testFiles DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

foreach(support_module ${LOCAL_SUPPORT_MODULES}) 
    set(TEST_NAME "test_${support_module}")
    message("building ${TEST_NAME}")
    add_executable(${TEST_NAME} "${TEST_NAME}.cpp")
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
    set_target_properties(${TEST_NAME} PROPERTIES CXX_STANDARD 11)
#    target_compile_options(${TEST_NAME} PRIVATE "-lrt")
    target_link_libraries(${TEST_NAME} PRIVATE ${support_module} ${${TEST_NAME}_SUPPORT_LIBS})
    get_target_property(OUT "test_${support_module}" LINK_LIBRARIES)

    if(NOT cmake_build_type_lower STREQUAL "release")
        target_include_directories(${TEST_NAME} PUBLIC
                                   ${Boost_INCLUDE_DIRS}  # Required for MacOSX CMake builds. Not sure if it's a CMake or clang issue.
                                   ${vector_tools_SOURCE_DIR}/${CPP_SRC_PATH}
                                   ${error_tools_SOURCE_DIR}/${CPP_SRC_PATH}
                                   ${support_module_SOURCE_DIR}/${CPP_SRC_PATH})
    endif()
endforeach(support_module)
