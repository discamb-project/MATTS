MACRO(add_exec)
    add_executable(${ARGN})
    target_link_libraries(${ARGV0} optimized discamb debug discambd)

    SET_PROPERTY(TARGET ${ARGV0} PROPERTY CXX_STANDARD 17)
    SET_PROPERTY(TARGET ${ARGV0} PROPERTY FOLDER "MATTS_programs")
    if(MT_MSVC_RUNTIME_LIB AND MSVC)
        set_property(TARGET ${ARGV0} PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    endif(MT_MSVC_RUNTIME_LIB AND MSVC)
    
ENDMACRO(add_exec)

