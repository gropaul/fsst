add_library(fsst_search INTERFACE)
target_include_directories(fsst_search INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(fsst_search INTERFACE cxx_std_17)

add_executable(prefilter_matcher_tests prefilter/scan/matcher/tests.cpp)
target_link_libraries(prefilter_matcher_tests fsst_search)
add_test(NAME prefilter_matcher COMMAND prefilter_matcher_tests)

add_executable(prefilter_resolver_tests prefilter/scan/resolver/tests.cpp)
target_link_libraries(prefilter_resolver_tests fsst_search)
add_test(NAME prefilter_resolver COMMAND prefilter_resolver_tests)

add_executable(prefilter_policy_tests prefilter/scan/policy/tests.cpp)
target_link_libraries(prefilter_policy_tests fsst_search)
add_test(NAME prefilter_policy COMMAND prefilter_policy_tests)

add_executable(prefilter_tests prefilter/tests.cpp)
target_link_libraries(prefilter_tests fsst_search fsst)
target_include_directories(prefilter_tests PRIVATE ${CMAKE_SOURCE_DIR})
add_test(NAME prefilter COMMAND prefilter_tests)

set(BENCH_DIR ${CMAKE_CURRENT_SOURCE_DIR}/prefilter/scan/bench)

add_executable(prefilter_corpus prefilter/scan/bench/corpus.cpp)
target_link_libraries(prefilter_corpus fsst)
target_include_directories(prefilter_corpus PRIVATE ${CMAKE_SOURCE_DIR})

add_executable(prefilter_matcher_sweep prefilter/scan/bench/matcher_fit.cpp)
target_link_libraries(prefilter_matcher_sweep fsst_search)
target_compile_definitions(prefilter_matcher_sweep PRIVATE BENCH_DIR="${BENCH_DIR}")

add_executable(prefilter_resolver_sweep prefilter/scan/bench/resolver_fit.cpp)
target_link_libraries(prefilter_resolver_sweep fsst_search)
target_compile_definitions(prefilter_resolver_sweep PRIVATE BENCH_DIR="${BENCH_DIR}")
