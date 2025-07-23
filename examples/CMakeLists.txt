function(add_example NAME)
  add_executable(${NAME} ${NAME}.cpp)
  target_link_libraries(${NAME} PRIVATE fuse)
endfunction()

add_example(move)
