set (HEADERS
)
set (SOURCES
    cpptest.cpp
)

source_group ("Header Files" FILES ${HEADERS})
source_group ("Source Files" FILES ${SOURCES})
add_executable (cpptest ${HEADERS} ${SOURCES})

target_link_libraries(cpptest cholmod_wrapper)

include_directories(..)

