cmake_minimum_required(VERSION 3.8)

project(cg_method)
set(CMAKE_CXX_STANDARD 14)            # Enable c++14 standard

set(TAKAHASHI_SOURCE takahashi.c print_util.c)
add_executable(takahashi ${TAKAHASHI_SOURCE})

set(EISENSTAT_SOURCE eisenstat.c print_util.c)
add_executable(eisenstat ${EISENSTAT_SOURCE})

set(TRI_SOURCE tri.c print_util.c)
add_executable(tri ${TRI_SOURCE})

set(CRSTRI_SOURCE crs_tri.c print_util.c)
add_executable(crs_tri ${CRSTRI_SOURCE})

set(GPS_SOURCE gps_opt.c gps_opt.h test_gpsopt.c print_util.c)
add_executable(gps ${GPS_SOURCE})

set(PARSE_GPS_SOURCE gps_opt.c mmio.c parse_and_gps.c cachecache_test.c print_util.c)
add_executable(PG ${PARSE_GPS_SOURCE})

set(CACHE_SOURCE cachecache.c cachecache_test.c print_util.c)
add_executable(cache ${CACHE_SOURCE})

set(GPS_CCE_SOURCE gps_cce.c gps_opt.c mmio.c cachecache.c print_util.c print_util.h)
add_executable(gps_cce ${GPS_CCE_SOURCE})