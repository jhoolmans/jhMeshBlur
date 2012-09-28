SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

cmake_minimum_required(VERSION 2.8)

PROJECT(jhMeshBlur)

# can optionally specify a version: FIND_PACKAGE( Maya 2012 REQUIRED )
FIND_PACKAGE( Maya REQUIRED )

SET( jhMeshBlur_SRCS
	src/jhMeshBlur.cpp
)

INCLUDE_DIRECTORIES( ${PUBLIC_INCLUDE_DIRS} ${jhMeshBlur_SOURCE_DIR}/src ${MAYA_INCLUDE_DIR} )

# create the plugin
ADD_LIBRARY( jhMeshBlur_plugin SHARED
  ${jhMeshBlur_SRCS}
)

TARGET_LINK_LIBRARIES( jhMeshBlur_plugin
   ${MAYA_Foundation_LIBRARY} ${MAYA_OpenMaya_LIBRARY} ${MAYA_OpenMayaAnim_LIBRARY}
)

SET_TARGET_PROPERTIES( jhMeshBlur_plugin PROPERTIES 
  OUTPUT_NAME jhMeshBlur
  CLEAN_DIRECT_OUTPUT 1
)

# this is a macro defined in FindMaya that sets up our plugin with standard Maya build settings
MAYA_SET_PLUGIN_PROPERTIES( jhMeshBlur_plugin )