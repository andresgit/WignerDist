TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -lfftw3

INCLUDEPATH += /home/andres/Documents/LMU/CompPhys/CompPhys2017/TOOLS
INCLUDEPATH += /home/andres/Documents/LMU/CompPhys/CompPhys2017/LIBRARIES/alglib

SOURCES += main.cpp \
    waveFunction.cpp \
    wigner.cpp \
    basisSet.cpp \
    orthopol.cpp \
    basisboundary.cpp

HEADERS += \
    waveFunction.h \
    wigner.h \
    basisSet.h \
    orthopol.h \
    basisboundary.h


win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../CompPhys2017/LIBRARIES/alglib/release/ -lalglib
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../CompPhys2017/LIBRARIES/alglib/debug/ -lalglib
else:unix: LIBS += -L$$PWD/../CompPhys2017/LIBRARIES/alglib/ -lalglib

INCLUDEPATH += $$PWD/../CompPhys2017/LIBRARIES/alglib
DEPENDPATH += $$PWD/../CompPhys2017/LIBRARIES/alglib

win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../CompPhys2017/LIBRARIES/alglib/release/libalglib.a
else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../CompPhys2017/LIBRARIES/alglib/debug/libalglib.a
else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../CompPhys2017/LIBRARIES/alglib/release/alglib.lib
else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../CompPhys2017/LIBRARIES/alglib/debug/alglib.lib
else:unix: PRE_TARGETDEPS += $$PWD/../CompPhys2017/LIBRARIES/alglib/libalglib.a
