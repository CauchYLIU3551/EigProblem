# by R.Lie, Nov 01, 2002

#include ../../Make.global_options

source = $(wildcard *.cpp)
object = $(patsubst %.cpp, %.o, $(source))
#LDFLAGS += -L~/AFEPack/library/lib -lAFEPack
CXXFLAGS = -I/usr/local/dealii-8.1.0/base/include -I/usr/local/dealii-8.1.0/lac/include -I/usr/local/dealii-8.1.0/include -I/usr/local/dealii-8.1.0/include/deal.II  -g -O2 -D__SERIALIZATION__ -DMULTITHREAD -pthread -fPIC -std=c++11 -ftemplate-depth-256 -I/usr/local/AFEPack/library/include 
LIBS = -ldeal_II -L/usr/local/dealii-8.1.0/lib  -lm -ldl 



all : main

%.o : %.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS)

main : $(object)
	$(CXX) -o $@ $(object) $(LDFLAGS) $(LIBS)

clean :
	-rm -rf $(object)
	-rm -rf main
	-rm -f *.[nes]
	-rm -f *.dx

.PHONY : default clean
