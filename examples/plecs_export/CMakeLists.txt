add_executable(plecs_export_example main.cpp)
target_link_libraries(plecs_export_example PRIVATE MKF)
# MKF's target_include_directories already exposes src/ and src/processors/
# as PUBLIC, so plecs_export_example inherits the include paths it needs.
set_target_properties(plecs_export_example PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
)
