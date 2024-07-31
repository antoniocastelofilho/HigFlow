
macro(SeparatePathFromInc var)
    string(REPLACE "-I" "" inc ${var})
    list(APPEND _INCS ${inc})
endmacro()