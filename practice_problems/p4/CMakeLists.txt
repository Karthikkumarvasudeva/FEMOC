# Set the name of the target
SET(TARGET "p4")



# Declare source files
SET(TARGET_SRC
   ${TARGET}.cc
)

# Create executable
ADD_EXECUTABLE(${TARGET} ${TARGET_SRC})

# Set up deal.II target
DEAL_II_SETUP_TARGET(${TARGET})
