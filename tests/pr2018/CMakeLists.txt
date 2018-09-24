# Set executable path.

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/tests/pr2018/bin)

# Add targets.

add_custom_target(pr2018 DEPENDS train_ring train_walks train_subgraph train_ml)

add_executable(train_ring src/train_ring.cpp)
target_link_libraries(train_ring gxlgedlib)
if(APPLE)
  add_custom_command(TARGET train_ring POST_BUILD COMMAND install_name_tool -change libnomad.so ${NOMAD_HOME}/lib/libnomad.so ${EXECUTABLE_OUTPUT_PATH}/train_ring)
  add_custom_command(TARGET train_ring POST_BUILD COMMAND install_name_tool -change libdoublefann.2.dylib ${FANN_HOME}/lib/libdoublefann.2.dylib ${EXECUTABLE_OUTPUT_PATH}/train_ring)
  add_custom_command(TARGET train_ring POST_BUILD COMMAND install_name_tool -change libsvm.so ${LIBSVM_HOME}/libsvm.so ${EXECUTABLE_OUTPUT_PATH}/train_ring)
endif()

add_executable(train_walks src/train_walks.cpp)
target_link_libraries(train_walks gxlgedlib)
if(APPLE)
  add_custom_command(TARGET train_walks POST_BUILD COMMAND install_name_tool -change libnomad.so ${NOMAD_HOME}/lib/libnomad.so ${EXECUTABLE_OUTPUT_PATH}/train_walks)
  add_custom_command(TARGET train_walks POST_BUILD COMMAND install_name_tool -change libdoublefann.2.dylib ${FANN_HOME}/lib/libdoublefann.2.dylib ${EXECUTABLE_OUTPUT_PATH}/train_walks)
  add_custom_command(TARGET train_walks POST_BUILD COMMAND install_name_tool -change libsvm.so ${LIBSVM_HOME}/libsvm.so ${EXECUTABLE_OUTPUT_PATH}/train_walks)
endif()

add_executable(train_subgraph src/train_subgraph.cpp)
target_link_libraries(train_subgraph gxlgedlib)
if(APPLE)
  add_custom_command(TARGET train_subgraph POST_BUILD COMMAND install_name_tool -change libnomad.so ${NOMAD_HOME}/lib/libnomad.so ${EXECUTABLE_OUTPUT_PATH}/train_subgraph)
  add_custom_command(TARGET train_subgraph POST_BUILD COMMAND install_name_tool -change libdoublefann.2.dylib ${FANN_HOME}/lib/libdoublefann.2.dylib ${EXECUTABLE_OUTPUT_PATH}/train_subgraph)
  add_custom_command(TARGET train_subgraph POST_BUILD COMMAND install_name_tool -change libsvm.so ${LIBSVM_HOME}/libsvm.so ${EXECUTABLE_OUTPUT_PATH}/train_subgraph)
endif()

add_executable(train_ml src/train_ml.cpp)
target_link_libraries(train_ml gxlgedlib)
if(APPLE)
  add_custom_command(TARGET train_ml POST_BUILD COMMAND install_name_tool -change libnomad.so ${NOMAD_HOME}/lib/libnomad.so ${EXECUTABLE_OUTPUT_PATH}/train_ml)
  add_custom_command(TARGET train_ml POST_BUILD COMMAND install_name_tool -change libdoublefann.2.dylib ${FANN_HOME}/lib/libdoublefann.2.dylib ${EXECUTABLE_OUTPUT_PATH}/train_ml)
  add_custom_command(TARGET train_ml POST_BUILD COMMAND install_name_tool -change libsvm.so ${LIBSVM_HOME}/libsvm.so ${EXECUTABLE_OUTPUT_PATH}/train_ml)
endif()