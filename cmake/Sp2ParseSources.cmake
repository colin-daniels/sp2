# MODULE:   ParseSources
#
# PROVIDES:
#   parse_sources(<output>
#       [[IF_<option>_ENABLED] file] ...)
#

function(parse_sources output)
    foreach(arg ${ARGN})
        if ("${arg}" MATCHES "^_+$")
            # spacer, do nothing
        elseif("${arg}" MATCHES "^IF_([^ \t]+)_ENABLED$")
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
        elseif("${arg}" MATCHES "^WITH_([^ \t]+)_FLAGS$")
            if(DEFINED SP2_${CMAKE_MATCH_1}_FLAGS)
                set(extra_flags "${extra_flags} ${SP2_${CMAKE_MATCH_1}_FLAGS}")
            endif()
        else()
            # arg is a file
            if(NOT ignore_arg)
                # everything is OK, add the argument to the output
                list(APPEND processed_sources "${arg}")
                if(DEFINED extra_flags)
                    set_property(SOURCE "${arg}" APPEND_STRING PROPERTY
                            COMPILE_FLAGS " ${extra_flags}")
                endif()
            endif()
            unset(extra_flags)
            set(ignore_arg FALSE)
        endif()
    endforeach()

    set(${output} ${processed_sources} PARENT_SCOPE)
endfunction()