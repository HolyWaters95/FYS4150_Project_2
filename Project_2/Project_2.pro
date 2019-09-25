TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
#QMAKE_CXXFLAGS += -O3

SOURCES += \
        #main.cpp \
        p2_functions.cpp \
        quantum.cpp \
        quantum2particles.cpp

INCLUDEPATH += C:\Qt\armadillo-9.700.2\include
DEPENDPATH += C:\Qt\armadillo-9.700.2\include


LIBS += \
    -LC:\Qt\armadillo-9.700.2\examples\lib_win64 \
    -llapack_win64_MT \
    -lblas_win64_MT

HEADERS += \
    catch.hpp
