CMAKE_MINIMUM_REQUIRED(VERSION 3.5.1)

PROJECT(PVblockMeshReader)

FIND_PACKAGE(ParaView REQUIRED)

IF(${ParaView_VERSION} VERSION_LESS 5.7)
    # ParaView-5.7.x and lower will prepend "lib" to the target name
    SET(TARGET_NAME PVblockMeshReader_SM)
ELSE()
    SET(TARGET_NAME libPVblockMeshReader_SM)
ENDIF()

IF(${ParaView_VERSION} VERSION_LESS 5.7)
    # ParaView-5.7.x and lower requires the paraview include file
    INCLUDE(${PARAVIEW_USE_FILE})
ENDIF()

LINK_DIRECTORIES(
    $ENV{FOAM_LIBBIN}
    $ENV{FOAM_EXT_LIBBIN}
    $ENV{PV_PLUGIN_PATH}
)

INCLUDE_DIRECTORIES(
    $ENV{WM_PROJECT_DIR}/src/OpenFOAM/lnInclude
    $ENV{WM_PROJECT_DIR}/src/OSspecific/$ENV{WM_OSTYPE}/lnInclude
    ${PROJECT_SOURCE_DIR}/../vtkPVblockMesh
    ${PROJECT_SOURCE_DIR}/vtk
)

ADD_DEFINITIONS(
    -std=c++17
    -DWM_$ENV{WM_PRECISION_OPTION}
    -DWM_LABEL_SIZE=$ENV{WM_LABEL_SIZE}
)

# Add the plugin
IF(${ParaView_VERSION} VERSION_LESS 5.7)
    ADD_PARAVIEW_PLUGIN(
        ${TARGET_NAME}
        "1.0"
        SERVER_MANAGER_XML PVblockMeshReader_SM.xml
        SERVER_MANAGER_SOURCES vtk/vtkPVblockMeshReader.cxx
    )
ELSE()
    # Paraview-5.7.x and higher builds the vtk module separately
    PARAVIEW_ADD_PLUGIN(
        ${TARGET_NAME}
        VERSION "1.0"
        SERVER_MANAGER_XML PVblockMeshReader_SM.xml
        MODULES PVblockMeshReader_VTK
        MODULE_FILES "vtk/vtk.module"
    )
ENDIF()

# Set the output library destination to the plugin directory
SET_TARGET_PROPERTIES(
    ${TARGET_NAME}
    PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "$ENV{PV_PLUGIN_PATH}"
)

# Build the plugin
IF(${ParaView_VERSION} VERSION_LESS 5.5)
    TARGET_LINK_LIBRARIES(
        ${TARGET_NAME}
        LINK_PUBLIC
        vtkPVblockMesh
        OpenFOAM
        LINK_PRIVATE
        pqCore
    )
ELSEIF(${ParaView_VERSION} VERSION_LESS 5.7)
    # Assume QT version 5 if it hasn't been defined
    IF(NOT PARAVIEW_QT_MAJOR_VERSION)
        set(PARAVIEW_QT_MAJOR_VERSION 5)
    ENDIF()
    # Paraview-5.5.x and higher needs QT libraries listed in the link command
    # in order for the headers to be available
    TARGET_LINK_LIBRARIES(
        ${TARGET_NAME}
        LINK_PUBLIC
        vtkPVblockMesh
        OpenFOAM
        Qt${PARAVIEW_QT_MAJOR_VERSION}::Core
        Qt${PARAVIEW_QT_MAJOR_VERSION}::Gui
        LINK_PRIVATE
        pqCore
    )
ELSE()
    # Assume QT version 5 if it hasn't been defined
    IF(NOT PARAVIEW_QT_MAJOR_VERSION)
        set(PARAVIEW_QT_MAJOR_VERSION 5)
    ENDIF()
    # Paraview-5.5.x and higher needs QT libraries listed in the link command
    # in order for the headers to be available
    TARGET_LINK_LIBRARIES(
        ${TARGET_NAME}
        LINK_PUBLIC
        vtkPVblockMesh
        OpenFOAM
        Qt${PARAVIEW_QT_MAJOR_VERSION}::Core
        Qt${PARAVIEW_QT_MAJOR_VERSION}::Gui
        LINK_PRIVATE
        ParaView::pqCore
    )
ENDIF()
