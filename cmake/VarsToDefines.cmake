# MODULE:   VarsToDefines
#
# PROVIDES:
#   vars_to_defines(<target> <INTERFACE|PUBLIC|PRIVATE> [items...])
#

function(vars_to_defines target visibility)
    foreach(var ${ARGN})
        if(${var})
            target_compile_definitions(${target} ${visibility}
                    ${var})
        endif()
    endforeach()
endfunction()
