#cd src
#gcc -c scan.c
#cd cuda
#gcc -c -I /usr/local/cuda/include cudasextractor.c
#nvcc -c cudatypes.cu
#nvcc -c -I /usr/local/cuda-5.0/samples/common/inc cudainit.cu -arch sm_30
#nvcc -c -I /usr/local/cuda-5.0/samples/common/inc cudascan.cu -arch sm_30
#nvcc -c -I /usr/local/cuda-5.0/samples/common/inc cudadeblend.cu -arch sm_30
#nvcc -c -I /usr/local/cuda-5.0/samples/common/inc cudafilter.cu -arch sm_30
#ar cru libcuda.a cudatypes.o cudainit.o cudascan.o cudadeblend.o cudafilter.o

#CUDPP required version 2.2
#CUDA required version 6.5
#first ``configure'' following sextractor instructions on "./doc/SExtractor installation - MediaWiki.html", then execute rebuild script(current file)
#when rebuild, please replace "/home/zhao/gwac/cudpp/cudpp_build/lib/libcudpp.so" with the cudpp library path
#on your own machine
cd src
gcc -c *.c
cd cuda
make
cd ..
nvcc  -O3 -g --compiler-options -funroll-loops --compiler-options -fomit-frame-pointer --compiler-options -Wall --compiler-options -D_REENTRANT  -o sex analyse.o assoc.o astrom.o back.o bpro.o catout.o check.o clean.o extract.o fft.o field.o filter.o fitswcs.o flag.o graph.o growth.o globals.o image.o interpolate.o main.o makeit.o manobjlist.o misc.o neurro.o pattern.o pc.o photom.o plist.o prefs.o profit.o psf.o readimage.o refine.o retina.o scan.o som.o weight.o winpos.o xml.o ../src/fits/libfits.a ../src/wcs/libwcs_c.a ../src/levmar/liblevmar.a ../src/cuda/libcuda.a -L/usr/local/atlas/lib -lcudpp -llapack -lptcblas -lcblas -latlas -L/usr/local/lib -lfftw3_threads -lfftw3 -lpthread -lm
cd ../
sudo make install
