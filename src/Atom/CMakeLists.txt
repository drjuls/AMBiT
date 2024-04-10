set(MODS_ATOM Atom.cpp
              Atom_Open.cpp
              Atom_Autoionization.cpp
              OutStreams.cpp
              CACHE INTERNAL "")

add_library(atom "${MODS_ATOM}")

target_link_libraries(atom PUBLIC externalfield)
