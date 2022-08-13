# Project

# Build the support modules
message("LOCAL_SUPPORT_MODULES: ${LOCAL_SUPPORT_MODULES}")
foreach(support_module ${LOCAL_SUPPORT_MODULES})
    add_library (${support_module} SHARED "${support_module}.cpp" "${support_module}.h")
    set_target_properties(${support_module} PROPERTIES CXX_STANDARD 11
                                                     PUBLIC_HEADER ${support_module}.h)
    target_link_libraries(${support_module} error_tools micromorphic_tools yaml-cpp xml2)
    target_compile_options(${support_module} PUBLIC)
endforeach(support_module)
