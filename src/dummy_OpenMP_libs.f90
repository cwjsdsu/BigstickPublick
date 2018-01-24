integer function omp_get_num_threads()
   omp_get_num_threads = 1
   return
end function omp_get_num_threads

integer function omp_get_thread_num()
   omp_get_thread_num = 0
   return
end function omp_get_thread_num

