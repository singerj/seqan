# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /sandbox/jsinger/apps/comparePeptideFromNuc/CMakeLists.txt
#
# CMakeLists.txt file for comparePeptideFromNuc.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_sandbox_jsinger_apps_comparePeptideFromNuc)
message (STATUS "Configuring sandbox/jsinger/apps/comparePeptideFromNuc")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
#set (SEQAN_FIND_DEPENDENCIES NONE)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your application.
add_executable (comparePeptideFromNuc comparePeptideFromNuc.cpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (comparePeptideFromNuc ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (comparePeptideFromNuc)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install comparePeptideFromNuc in ${PREFIX}/bin directory
install (TARGETS comparePeptideFromNuc
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/comparePeptideFromNuc for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (comparePeptideFromNuc)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:comparePeptideFromNuc")
  set (CPACK_PACKAGE_NAME "comparePeptideFromNuc")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "comparePeptideFromNuc")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Your Name <your.email@example.net>")
  set (CPACK_PACKAGE_VENDOR "Your Name <your.email@example.net>")

  seqan_configure_cpack_app (comparePeptideFromNuc "comparePeptideFromNuc")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:comparePeptideFromNuc")

