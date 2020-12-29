cmake_minimum_required(VERSION 3.8)

project(cg_method)
set(CMAKE_CXX_STANDARD 14)            # Enable c++14 standard

set(TAKAHASHI_SOURCE takahashi.c)
add_executable(takahashi ${TAKAHASHI_SOURCE})

set(EISENSTAT_SOURCE eisenstat.c)
add_executable(eisenstat ${EISENSTAT_SOURCE})

set(TRI_SOURCE tri.c)
add_executable(tri ${TRI_SOURCE})

set(CRSTRI_SOURCE crs_tri.c)
add_executable(crs_tri ${CRSTRI_SOURCE})

set(GPS_SOURCE gps_opt.c gps_opt.h test_gpsopt.c)
add_executable(gps ${GPS_SOURCE})

set(PARSE_GPS_SOURCE gps_opt.c mmio.c parse_and_gps.c)
add_executable(PG ${PARSE_GPS_SOURCE})

set(CACHE_SOURCE cachecache.c)
add_executable(cache ${CACHE_SOURCE})