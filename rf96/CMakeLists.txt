enable_language(Fortran)
file ( GLOB rf96_fortran_functions *.f* )
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
Python_add_library(rf96 MODULE WITH_SOABI ${rf96_fortran_functions})

install(TARGETS rf96 DESTINATION ./pyrf96)
