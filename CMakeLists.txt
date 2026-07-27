cmake_minimum_required(VERSION 3.21)

project(Stochia VERSION 0.2.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
option(STOCHIA_BUILD_APP "Build the Qt desktop application" ON)

add_library(stochia_core
    src/core/Distribution.cpp
    src/core/Expression.cpp
    src/core/ProbabilityModel.cpp
)
target_include_directories(stochia_core PUBLIC src)

if(STOCHIA_BUILD_APP)
    find_package(Qt6 6.4 QUIET COMPONENTS Widgets)
    if(Qt6_FOUND)
        set(CMAKE_AUTOMOC ON)
        add_executable(Stochia WIN32
            src/main.cpp
            src/ui/MainWindow.cpp
            src/ui/PlotWidget.cpp
            src/ui/VariableDialog.cpp
        )
        target_link_libraries(Stochia PRIVATE stochia_core Qt6::Widgets)
    else()
        message(WARNING "Qt 6.4+ Widgets was not found; building the tested mathematical core only.")
    endif()
endif()

if(MSVC)
    target_compile_options(stochia_core PRIVATE /W4 /permissive- /utf-8)
    if(TARGET Stochia)
        target_compile_options(Stochia PRIVATE /W4 /permissive- /utf-8)
    endif()
else()
    target_compile_options(stochia_core PRIVATE -Wall -Wextra -Wpedantic)
    if(TARGET Stochia)
        target_compile_options(Stochia PRIVATE -Wall -Wextra -Wpedantic)
    endif()
endif()

include(CTest)
if(BUILD_TESTING)
    add_executable(stochia_core_tests tests/core_tests.cpp)
    target_link_libraries(stochia_core_tests PRIVATE stochia_core)
    add_test(NAME stochia_core_tests COMMAND stochia_core_tests)
endif()

if(TARGET Stochia)
    install(TARGETS Stochia RUNTIME DESTINATION bin)
endif()
