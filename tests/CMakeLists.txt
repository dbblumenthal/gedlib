add_custom_target(tests DEPENDS unit_tests pr2018 sspr2018)
add_subdirectory(unit_tests)
add_subdirectory(sspr2018)
add_subdirectory(pr2018)