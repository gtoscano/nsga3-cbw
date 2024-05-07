cmake_minimum_required(VERSION 3.14)
project(nsga3)

include(FetchContent)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

SET(CMAKE_CXX_FLAGS "-O3")
#SET(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
set(THREADS_PREFER_PTHREAD_FLAG ON)

# -DBOOST_LOG_DYN_LINK
# definitions
add_definitions(-DBOOST_ALL_NO_LIB)
add_definitions(-DBOOST_LOG_DYN_LINK=1)

find_package(Boost 1.74.0 COMPONENTS thread log log_setup REQUIRED)
find_package(PkgConfig)

# Check libraries of loki: pkg-config --libs loki-cpp
pkg_check_modules(PROTOBUF REQUIRED IMPORTED_TARGET protobuf)
pkg_check_modules(SNAPPY REQUIRED IMPORTED_TARGET snappy)
pkg_check_modules(CURLPP REQUIRED IMPORTED_TARGET curlpp)
pkg_check_modules(LIBCURL REQUIRED IMPORTED_TARGET libcurl)

set(Parquet_DIR /usr/lib/x86_64-linux-gnu/cmake/Arrow)

find_package(Threads REQUIRED)
find_package(Arrow CONFIG REQUIRED)
find_package(Parquet CONFIG REQUIRED)
find_package(crossguid CONFIG REQUIRED)
find_package(aws-cpp-sdk-core CONFIG REQUIRED)
find_package(aws-cpp-sdk-s3 CONFIG REQUIRED)
find_package(Protobuf MODULE REQUIRED)

file(GLOB_RECURSE SRC_FILES src/*.cpp)
add_executable(nsga3 ${SRC_FILES})
target_include_directories(nsga3 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(nsga3 PUBLIC Boost::log Boost::log_setup hiredis redis++ fmt arrow_shared parquet_shared crossguid aws-cpp-sdk-core aws-cpp-sdk-s3 SimpleAmqpClient Threads::Threads PkgConfig::LIBCURL) #PkgConfig::LOKI-CPP 

install(TARGETS nsga3 RUNTIME DESTINATION /usr/local/bin)

