#
# Cross Platform Makefile
# Compatible with MSYS2/MINGW, Ubuntu 14.04.1 and Mac OS X
#
# You will need GLFW (http://www.glfw.org):
# Linux:
#   apt-get install libglfw-dev
# Mac OS X:
#   brew install glfw
# MSYS2:
#   pacman -S --noconfirm --needed mingw-w64-x86_64-toolchain mingw-w64-x86_64-glfw
#

CXX = g++

TARGET = opengl3_example
IMGUI_DIR = ../third_party/imgui
OBJS = main.o
OBJS += $(IMGUI_DIR)/imgui_impl_glfw_gl3.o $(IMGUI_DIR)/imgui.o $(IMGUI_DIR)/imgui_demo.o $(IMGUI_DIR)/imgui_draw.o
OBJS += ../third_party/gl3w/GL/gl3w.o

UNAME_S := $(shell uname -s)


ifeq ($(UNAME_S), Linux) #LINUX
	ECHO_MESSAGE = "Linux"
	LIBS = -lGL `pkg-config --static --libs glfw3`

	CXXFLAGS = -I../third_party/imgui -I../third_party/gl3w `pkg-config --cflags glfw3`
	CXXFLAGS += -Wall -Wformat
	CFLAGS = $(CXXFLAGS)
endif

ifeq ($(UNAME_S), Darwin) #APPLE
	ECHO_MESSAGE = "Mac OS X"
	LIBS = -framework OpenGL -framework Cocoa -framework IOKit -framework CoreVideo
	#LIBS += -L/usr/local/lib -lglfw3
	LIBS += -L/usr/local/lib -lglfw

	CXXFLAGS = -I../../ -I../libs/gl3w -I/usr/local/include
	CXXFLAGS += -Wall -Wformat
	CFLAGS = $(CXXFLAGS)
endif

ifeq ($(findstring MINGW,$(UNAME_S)),MINGW)
   ECHO_MESSAGE = "Windows"
   LIBS = -lglfw3 -lgdi32 -lopengl32 -limm32

   CXXFLAGS = -I../../ -I../libs/gl3w `pkg-config --cflags glfw3`
   CXXFLAGS += -Wall -Wformat
   CFLAGS = $(CXXFLAGS)
endif


.cpp.o:
	$(CXX) $(CXXFLAGS) -c -o $@ $<

all: $(TARGET)
	@echo Build complete for $(ECHO_MESSAGE)

$(TARGET): $(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(CXXFLAGS) $(LIBS)

clean:
	rm $(TARGET) $(OBJS)
