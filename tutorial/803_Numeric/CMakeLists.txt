get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${PROJECT_NAME})

set(test_sources
    test_utils.cpp
    main.cpp
)

add_executable(${PROJECT_NAME}_bin ${test_sources})
target_link_libraries(${PROJECT_NAME}_bin igl::core tutorials)
# find_package(OpenMP)
# if(OpenMP_CXX_FOUND)
#     target_link_libraries(${PROJECT_NAME}_bin OpenMP::OpenMP_CXX)
# endif()
