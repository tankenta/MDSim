#pragma once
#include <iostream>
#include <GL/gl3w.h>
#include <GLFW/glfw3.h>

class Window {
private:
    GLFWwindow *const window;
    GLfloat size[2];
    GLfloat scale;
    GLfloat location[2];

public:
    Window(int width=1280, int height=720, const char *title="ImGui OpenGL3 example")
    : window(glfwCreateWindow(width, height, title, NULL, NULL)),
        scale(100.0f), location{0, 0} {
        if (window == NULL) {
            std::cerr << "Failed to create GLFW window." << std::endl;
            exit(1);
        }

        glfwMakeContextCurrent(window);
        glfwSwapInterval(1);    // enable vsync
        glfwSetWindowSizeCallback(window, resize);
        glfwSetScrollCallback(window, wheel);
        glfwSetWindowUserPointer(window, this);
        gl3wInit();
        resize(window, width, height);
    }

    virtual ~Window() {
        glfwDestroyWindow(window);
    }

    int shouldClose() const {
        return glfwWindowShouldClose(window);
    }

    void swapBuffers() {
        glfwSwapBuffers(window);
        glfwPollEvents();

        if (glfwGetKey(window, GLFW_KEY_N) != GLFW_RELEASE) {
            scale -= 1.0f;
        } else if (glfwGetKey(window, GLFW_KEY_P) != GLFW_RELEASE) {
            scale += 1.0f;
        }

        if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_1) != GLFW_RELEASE) {
            double x, y;
            glfwGetCursorPos(window, &x, &y);

            location[0] = static_cast<GLfloat>(x) * 2.0f / size[0] - 1.0f;
            location[1] = 1.0f - static_cast<GLfloat>(y) * 2.0f / size[1];
        }
    }

    GLFWwindow *getWindow() const {
        return window;
    }

    const GLfloat *getSize() const {
        return size;
    }

    GLfloat getScale() const {
        return scale;
    }

    const GLfloat *getLocation() const {
        return location;
    }

    static void resize(GLFWwindow *const window, int width, int height) {
        glViewport(0, 0, width, height);

        Window *const instance = 
            static_cast<Window *>(glfwGetWindowUserPointer(window));
        if (instance != NULL) {
            instance->size[0] = static_cast<GLfloat>(width);
            instance->size[1] = static_cast<GLfloat>(height);
        }
    }

    static void wheel(GLFWwindow *const window, double x, double y) {
        Window *const instance = 
            static_cast<Window *>(glfwGetWindowUserPointer(window));
        if (instance != NULL) {
            instance->scale += static_cast<GLfloat>(y);
        }
    }
};
