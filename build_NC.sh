cmake -B build -S . -DENABLE_NETCDF=ON \
  -DNETCDF_ROOT=$EBROOTNETCDF \
  -DNETCDF_C_LIBRARY=$EBROOTNETCDF/lib/libnetcdf.so \
  -DHDF5_HL_LIBRARY=$EBROOTHDF5/lib/libhdf5_hl.so \
  -DHDF5_LIBRARY=$EBROOTHDF5/lib/libhdf5.so
cmake --build build -j 8