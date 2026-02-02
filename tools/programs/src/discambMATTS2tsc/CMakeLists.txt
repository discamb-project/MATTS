add_exec(discambMATTS2tsc
               discambMATTS2tsc.cpp MATTS_Default.cpp MATTS_Default.h)

target_link_libraries(discambMATTS2tsc discamb OpenMP::OpenMP_CXX)
