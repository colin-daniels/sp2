# MODULE:   ParseSources
#
# PROVIDES:
#   parse_sources(<output>
#       [[IF_<option>_ENABLED] file] ...)
#

function(parse_sources output)
    foreach(arg ${ARGN})
        if("${arg}" MATCHES "^IF_([^ \t]+)_ENABLED$")
            if(DEFINED SP2_ENABLE_${CMAKE_MATCH_1})
                # if the argument matches an option for SP2, enable or disable
                # the next argument accordingly
                if(SP2_ENABLE_${CMAKE_MATCH_1})
                    set(ignore_arg FALSE)
                else()
                    set(ignore_arg TRUE)
                endif()
            else()
                set(ignore_arg FALSE)
            endif()
        elseif(ignore_arg)
            set(ignore_arg FALSE)
        else()
            # everything is OK, add the argument to the output
            list(APPEND processed_sources "${arg}")
        endif()
    endforeach()

    set(${output} ${processed_sources} PARENT_SCOPE)
endfunction()