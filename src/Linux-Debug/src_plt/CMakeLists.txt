include_directories(${CUDA_INCLUDE_DIRS})
include_directories(${THRUST_INCLUDE_DIR})
include_directories(${SRC}) #include for higher source files include

add_library(pltLib STATIC
	Advance_Positions_Plt.cu
	Advance_Positions_Plt.h
	Plt_Vol_Exc_Force.cu
	Plt_Vol_Exc_Force.h
	functor_plt_vol_exc.h
	Plt_Field_Plt_Force.cu
	Plt_Field_Plt_Force.h
	functor_plt_field_plt.h
	Plt_Field_Node_Force.cu
	Plt_Field_Node_Force.h
	functor_plt_field_node.h
	Plt_Arm_Plt_Force.cu
	Plt_Arm_Plt_Force.h
	functor_plt_arm_plt.h
	Plt_Arm_Node_Force.cu
	Plt_Arm_Node_Force.h
	functor_plt_arm_node.h
	Plt_Arm_Node_Force_Time.cu
	Plt_Arm_Node_Force_Time.h
	functor_plt_arm_node_time.h)

target_compile_features(pltLib PUBLIC cxx_std_11)