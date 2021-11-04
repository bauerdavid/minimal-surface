cmake_minimum_required(VERSION 3.10.2)

project(MinArea_SITK)

# Find ITK.
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
set(SOURCES AreaEikonal.cpp AreaEikonal.h ChildFrm.cpp ChildFrm.h ChildView.cpp ChildView.h commontype.h FitPlane.cpp FitPlane.h ImageOp.cpp ImageOp.h MainFrm.cpp MainFrm.h MinArea.cpp MinArea.h resource.h stdafx.cpp stdafx.h targetver.h Transport.cpp Transport.h)
add_executable(MinArea ${SOURCES})

target_link_libraries(MinArea ${ITK_LIBRARIES})