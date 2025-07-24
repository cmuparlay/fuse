set(COMMON_DEFS "")

function(add_example NAME SOURCE DEFS)
  add_executable(${NAME} ${SOURCE})
  target_compile_definitions(${NAME} PRIVATE ${COMMON_DEFS})
  target_compile_definitions(${NAME} PRIVATE ${DEFS})
  target_link_libraries(${NAME} PRIVATE fuse)
endfunction()

add_example(move move.cpp "")
add_example(move_mv_stm move_verlib.cpp "MV_STM;HashLock;NoHelp")
add_example(move_mv_tlf move_verlib.cpp "MV_TLF;HashLock;NoHelp")
add_example(move_mv_tlf_stm move_verlib.cpp "MV_TLF_STM;HashLock;NoHelp")
add_example(move_2plsf_stm move_verlib.cpp "TWOPLSF;SV_STM")
add_example(move_2plsf_tlf move_verlib.cpp "TWOPLSF;SV_TLF")
add_example(move_2plsf_tlf_stm move_verlib.cpp "TWOPLSF;SV_TLF_STM")

