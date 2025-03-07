# Set the name of the target
SET(TARGET "HM_Biharmonic")



# Declare source files
SET(TARGET_SRC
   ${TARGET}.cc
)

# Create executable
ADD_EXECUTABLE(${TARGET} ${TARGET_SRC})

# Set up deal.II target
DEAL_II_SETUP_TARGET(${TARGET})


#target_link_libraries(p11 deal.II::deal.II)
