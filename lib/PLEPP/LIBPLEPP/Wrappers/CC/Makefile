include ../../Makefile.in 
#

SHARED = $(FCOMPILER) -shared
STATIC = ar rcs 

FLAGS =


all:
	@echo
	#====================================================# 
	@echo "--> make [test01/test02/libcommdom/shared]"
	#====================================================#
	@echo 


libcommdom: 
	$(CPPCOMPILER) -c -g -fPIC $(INCLUDES)/build_octree.c $(INCLUDE_PLE) 
	$(CPPCOMPILER) -c -g -fPIC $(INCLUDES)/commdom.cpp $(INCLUDE_PLE) $(FLAGS)
	$(CPPCOMPILER) -c -g -fPIC commdom_wrapper.cpp $(INCLUDE_PLE) 


static: 
	$(CPPCOMPILER) -c $(INCLUDES)/build_octree.c $(INCLUDE_PLE) 
	$(CPPCOMPILER) -c $(INCLUDES)/commdom.cpp $(INCLUDE_PLE) $(FLAGS) 
	$(CCOMPILER) -c commdom_wrapper.cpp $(INCLUDE_PLE) $(FLAGS) 
	$(STATIC) libcommdom.a build_octree.o commdom.o commdom_wrapper.o $(ROOT_PLE)/lib/libple.a


clean:
	rm -rf *~ *.o *.x *.mod *.so *.a
 
