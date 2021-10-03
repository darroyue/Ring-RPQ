g++ -fpermissive -O9 -msse4.2 -std=c++11 -DNODEBUG -I ~/include -L ~/lib src/nfa/src/bitmasks.c src/nfa/src/options.c src/nfa/src/parser.c src/nfa/src/regular.c src/build-index.cpp -o build-index -lsdsl -ldivsufsort -ldivsufsort64

g++ -fpermissive -Ofast -msse4.2 -frename-registers -std=c++11 -DNODEBUG -I ~/include -L ~/lib src/nfa/src/bitmasks.c src/nfa/src/options.c src/nfa/src/parser.c src/nfa/src/regular.c src/query-index.cpp -o query-index -lsdsl -ldivsufsort -ldivsufsort64


