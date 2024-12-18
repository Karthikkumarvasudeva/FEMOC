# Set the name of the target
SET(TARGET "p7")



# Declare source files
SET(TARGET_SRC
   ${TARGET}.cc
)

set(CMAKE_CXX_STANDARD 17)
find_package(deal.II REQUIRED)
deal_ii_initialize_cached_variables()

# Create executable
ADD_EXECUTABLE(${TARGET} ${TARGET_SRC})

# Set up deal.II target
DEAL_II_SETUP_TARGET(${TARGET})
