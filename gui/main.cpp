// ImGui - standalone example application for GLFW + OpenGL 3, using programmable pipeline
// If you are new to ImGui, see examples/README.txt and documentation at the top of imgui.cpp.
// (GLFW is a cross-platform general purpose library for handling windows, inputs, OpenGL/Vulkan graphics context creation, etc.)
// (GL3W is a helper library to access OpenGL functions since there is no standard header to access modern OpenGL functions easily. Alternatives are GLEW, Glad, etc.)

#include "imgui.h"
#include "imgui_impl_glfw_gl3.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <cmath>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <GL/gl3w.h>    // This example is using gl3w to access OpenGL functions (because it is small). You may use glew/glad/glLoadGen/etc. whatever already works for you.
#include <GLFW/glfw3.h>

#include "Window.h"
#include "Shape.h"
#include "Uniform.h"
#include "Material.h"

#include <string>
#define EIGEN_NO_DEBUG
#define EIGEN_MPL2_ONLY
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "md_utils.hpp"

struct SLposPolar{
    float radius;
    float theta;
    float phi;
    glm::vec3 ortho() {
        return glm::vec3(
                radius * std::sin(theta) * std::cos(phi),
                radius * std::sin(theta) * std::sin(phi),
                radius * std::cos(theta)
        );
    }
};

GLboolean printShaderInfoLog(GLuint shader, const char *str) {
    GLint status;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &status);
    if (status == GL_FALSE) {
        std::cerr << "Compile Error in " << str << std::endl;
    }

    GLsizei buf_size;
    glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &buf_size);

    if (buf_size > 1) {
        std::vector<GLchar> info_log(buf_size);
        GLsizei length;
        glGetShaderInfoLog(shader, buf_size, &length, &info_log[0]);
        std::cerr << &info_log[0] << std::endl;
    }
    return static_cast<GLboolean>(status);
}

GLboolean printProgramInfoLog(GLuint program) {
    GLint status;
    glGetProgramiv(program, GL_LINK_STATUS, &status);
    if (status == GL_FALSE) {
        std::cerr << "Link Error." << std::endl;
    }

    GLsizei buf_size;
    glGetProgramiv(program, GL_INFO_LOG_LENGTH, &buf_size);

    if (buf_size > 1) {
        std::vector<GLchar> info_log(buf_size);
        GLsizei length;
        glGetProgramInfoLog(program, buf_size, &length, &info_log[0]);
        std::cerr << &info_log[0] << std::endl;
    }
    return static_cast<GLboolean>(status);
}

GLuint createProgram(const char *vsrc, const char *fsrc) {
    const GLuint program = glCreateProgram();

    if (vsrc != NULL) {
        const GLuint vobj = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(vobj, 1, &vsrc, NULL);
        glCompileShader(vobj);

        if (printShaderInfoLog(vobj, "vertex shader")) {
            glAttachShader(program, vobj);
        }
        glDeleteShader(vobj);
    }

    if (fsrc != NULL) {
        const GLuint fobj = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(fobj, 1, &fsrc, NULL);
        glCompileShader(fobj);

        if (printShaderInfoLog(fobj, "fragment shader")) {
            glAttachShader(program, fobj);
        }
        glDeleteShader(fobj);
    }

    glBindAttribLocation(program, 0, "position");
    glBindAttribLocation(program, 1, "normal");
    glBindFragDataLocation(program, 0, "fragment");
    glLinkProgram(program);

    if (printProgramInfoLog(program)) {
        return program;
    } else {
        glDeleteProgram(program);
        return 0;
    }
}

bool readShaderSource(const char *name, std::vector<GLchar> &buffer) {
    if (name == NULL) {
        return false;
    }

    std::ifstream file(name, std::ios::binary);
    if (file.fail()) {
        std::cerr << "Error: Failed to open source file: " << name << std::endl;
        return false;
    }

    file.seekg(0L, std::ios::end);
    GLsizei length = static_cast<GLsizei>(file.tellg());
    buffer.resize(length + 1);

    file.seekg(0L, std::ios::beg);
    file.read(buffer.data(), length);
    buffer[length] = '\0';

    if (file.fail()) {
        std::cerr << "Error: Failed to read source file: " << name << std::endl;
        file.close();
        return false;
    }

    file.close();
    return true;
}

GLuint loadProgram(const char *vert, const char *frag) {
    std::vector<GLchar> vsrc;
    const bool vstat = readShaderSource(vert, vsrc);
    std::vector<GLchar> fsrc;
    const bool fstat = readShaderSource(frag, fsrc);

    return vstat && fstat ? createProgram(vsrc.data(), fsrc.data()) : 0;
}

static void error_callback(int error, const char* description)
{
    fprintf(stderr, "Error %d: %s\n", error, description);
}

int main(int, char**)
{
    const double dt = 0.005;
    const double total_time = 10.0;
    const double temp_cont_time = 1.0;
    std::string bc_mode = "periodic";
    double number_density = 0.8;

    const double target_temp = 1.0;
    const int num_particles = 256;
    const double ptcl_mass = 1.0;
    const int RDF_hist_size = 200;
    Eigen::Vector3i cube_size = Eigen::Vector3i::Constant(4);

    MDSim md(
            dt, total_time, temp_cont_time, number_density, bc_mode,
            target_temp, num_particles, ptcl_mass, RDF_hist_size, cube_size);


    // Setup window
    glfwSetErrorCallback(error_callback);
    if (!glfwInit())
        return 1;
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#if __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif
    Window window;

    glFrontFace(GL_CCW);
    glCullFace(GL_BACK);
    glEnable(GL_CULL_FACE);

    glClearDepth(1.0);
    glDepthFunc(GL_LESS);
    glEnable(GL_DEPTH_TEST);

    const GLuint program = loadProgram("point.vert", "point.frag");
    const GLint modelview_loc = glGetUniformLocation(program, "modelview");
    const GLint projection_loc = glGetUniformLocation(program, "projection");
    const GLint normalmat_loc = glGetUniformLocation(program, "normalmat");
    const GLint Lpos_loc = glGetUniformLocation(program, "Lpos");
    const GLint Lamb_loc = glGetUniformLocation(program, "Lamb");
    const GLint Ldiff_loc = glGetUniformLocation(program, "Ldiff");
    const GLint Lspec_loc = glGetUniformLocation(program, "Lspec");

    const GLint material_loc = glGetUniformBlockIndex(program, "Material");
    glUniformBlockBinding(program, material_loc, 0);

    const int slices = 32, stacks = 16;
    std::vector<Object::Vertex> solidSphereVertex;
    for (int j = 0; j <= stacks; j++) {
        const float t = static_cast<float>(j) / static_cast<float>(stacks);
        const float y = std::cos(glm::pi<float>() * t);
        const float r = std::sin(glm::pi<float>() * t);
        for (int i = 0; i <= slices; i++) {
            const float s = static_cast<float>(i) / static_cast<float>(slices);
            const float z = r * std::cos(2.0f * glm::pi<float>() * s);
            const float x = r * std::sin(2.0f * glm::pi<float>() * s);
            const Object::Vertex v = { x, y, z, x, y, z };
            solidSphereVertex.emplace_back(v);
        }
    }
    std::vector<GLuint> solidSphereIndex;
    for (int j = 0; j < stacks; j++) {
        const int k = (slices + 1) * j;
        for (int i = 0; i < slices; i++) {
            const GLuint k0 = k + i;
            const GLuint k1 = k0 + 1;
            const GLuint k2 = k1 + slices;
            const GLuint k3 = k2 + 1;

            solidSphereIndex.emplace_back(k0);
            solidSphereIndex.emplace_back(k2);
            solidSphereIndex.emplace_back(k3);

            solidSphereIndex.emplace_back(k0);
            solidSphereIndex.emplace_back(k3);
            solidSphereIndex.emplace_back(k1);
        }
    }
    std::unique_ptr<const Shape> shape(
            new SolidShapeIndex(3,
                static_cast<GLsizei>(solidSphereVertex.size()), solidSphereVertex.data(),
                static_cast<GLsizei>(solidSphereIndex.size()), solidSphereIndex.data()));

    auto Lamb = glm::vec3(0.2f);
    auto Ldiff = glm::vec3(1.0f);
    auto Lspec = glm::vec3(1.0f);

    static constexpr Material color[] = {
        // Kamb              Kdiff              Kspec              Kshi
        { 0.9f, 0.1f, 0.1f,  0.9f, 0.1f, 0.1f,  0.3f, 0.3f, 0.3f,  30.0f },
        { 0.1f, 0.1f, 0.9f,  0.1f, 0.1f, 0.9f,  0.4f, 0.4f, 0.4f,  60.0f }
    };
    const Uniform<Material> material(color, 2);

    glfwSetTime(0.0);

    // Setup ImGui binding
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    ImGui_ImplGlfwGL3_Init(window.getWindow(), true);
    //io.NavFlags |= ImGuiNavFlags_EnableKeyboard;  // Enable Keyboard Controls
    //io.NavFlags |= ImGuiNavFlags_EnableGamepad;   // Enable Gamepad Controls

    // Setup style
    ImGui::StyleColorsDark();
    //ImGui::StyleColorsClassic();

    // Load Fonts
    // - If no fonts are loaded, dear imgui will use the default font. You can also load multiple fonts and use ImGui::PushFont()/PopFont() to select them. 
    // - AddFontFromFileTTF() will return the ImFont* so you can store it if you need to select the font among multiple. 
    // - If the file cannot be loaded, the function will return NULL. Please handle those errors in your application (e.g. use an assertion, or display an error and quit).
    // - The fonts will be rasterized at a given size (w/ oversampling) and stored into a texture when calling ImFontAtlas::Build()/GetTexDataAsXXXX(), which ImGui_ImplXXXX_NewFrame below will call.
    // - Read 'misc/fonts/README.txt' for more instructions and details.
    // - Remember that in C/C++ if you want to include a backslash \ in a string literal you need to write a double backslash \\ !
    //io.Fonts->AddFontDefault();
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/Roboto-Medium.ttf", 16.0f);
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/Cousine-Regular.ttf", 15.0f);
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/DroidSans.ttf", 16.0f);
    //io.Fonts->AddFontFromFileTTF("../../misc/fonts/ProggyTiny.ttf", 10.0f);
    //ImFont* font = io.Fonts->AddFontFromFileTTF("c:\\Windows\\Fonts\\ArialUni.ttf", 18.0f, NULL, io.Fonts->GetGlyphRangesJapanese());
    //IM_ASSERT(font != NULL);

    bool show_demo_window = true;
    bool show_another_window = false;
    ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);

    // Main loop
    while (!window.shouldClose())
    {
        // You can read the io.WantCaptureMouse, io.WantCaptureKeyboard flags to tell if dear imgui wants to use your inputs.
        // - When io.WantCaptureMouse is true, do not dispatch mouse input data to your main application.
        // - When io.WantCaptureKeyboard is true, do not dispatch keyboard input data to your main application.
        // Generally you may always pass all inputs to dear imgui, and hide them from your application based on those two flags.
        //glfwPollEvents();
        ImGui_ImplGlfwGL3_NewFrame();

        // 1. Show a simple window.
        // Tip: if we don't call ImGui::Begin()/ImGui::End() the widgets automatically appears in a window called "Debug".
        static float scale = 0.25f;
        static float cell_angle = 0.0f;
        static SLposPolar Lpos_polar = { 30.0f, 1.0f/2.0f*glm::pi<float>(), 1.0f/6.0f*glm::pi<float>() };
        {
            static float f = 0.0f;
            static int counter = 0;
            ImGui::Text("Hello, world!");                           // Display some text (you can use a format string too)
            ImGui::SliderFloat("float", &f, 0.0f, 1.0f);            // Edit 1 float using a slider from 0.0f to 1.0f    
            ImGui::ColorEdit3("clear color", (float*)&clear_color); // Edit 3 floats representing a color

            ImGui::Checkbox("Demo Window", &show_demo_window);      // Edit bools storing our windows open/close state
            ImGui::Checkbox("Another Window", &show_another_window);

            if (ImGui::Button("Button"))                            // Buttons return true when clicked (NB: most widgets return true when edited/activated)
                counter++;
            ImGui::SameLine();
            ImGui::Text("counter = %d", counter);

            ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);

            ImGui::DragFloat("particle radius", &scale, 0.01f, 0.0f, 100.0f, NULL, 3.0f);
            ImGui::SliderAngle("cell angle", &cell_angle);
            if (ImGui::TreeNode("Light position")) {
                ImGui::DragFloat("radius", &Lpos_polar.radius, 0.01f, 0.0f, 100.0f, NULL, 3.0f);
                ImGui::SliderAngle("angle1", &Lpos_polar.theta);
                ImGui::SliderAngle("angle2", &Lpos_polar.phi);
            }
        }

        // 2. Show another simple window. In most cases you will use an explicit Begin/End pair to name your windows.
        if (show_another_window)
        {
            ImGui::Begin("Another Window", &show_another_window);
            ImGui::Text("Hello from another window!");
            if (ImGui::Button("Close Me"))
                show_another_window = false;
            ImGui::End();
        }

        // 3. Show the ImGui demo window. Most of the sample code is in ImGui::ShowDemoWindow(). Read its code to learn more about Dear ImGui!
        if (show_demo_window)
        {
            ImGui::SetNextWindowPos(ImVec2(650, 20), ImGuiCond_FirstUseEver); // Normally user code doesn't need/want to call this because positions are saved in .ini file anyway. Here we just want to make the demo initial state a bit more friendly!
            ImGui::ShowDemoWindow(&show_demo_window);
        }

        // Rendering
        glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glUseProgram(program);

        const GLfloat *const size = window.getSize();
        //const GLfloat scale = window.getScale() * 2.0f;
        //const GLfloat w = size[0] / scale, h = size[1] / scale;
        const GLfloat fovy = window.getScale() * 0.01f;
        const GLfloat aspect = size[0] / size[1];
        const glm::mat4 projection = glm::perspective(fovy, aspect, 0.1f, 100.0f);

        const GLfloat *const location = window.getLocation();
        const glm::mat4 global_rot = glm::rotate(
                cell_angle, glm::vec3(0.0f, 1.0f, 0.0f));
        const glm::mat4 ptcl_scale = glm::scale(glm::mat4(1.0f), glm::vec3(scale));
        const glm::mat4 global_trans = glm::translate(
                glm::mat4(1.0f), glm::vec3(location[0], location[1], 0.0f));
        const glm::mat4 view = glm::lookAt(
                glm::vec3(15.0f, 20.0f, 25.0f),
                glm::vec3(0.0f, 0.0f, 0.0f),
                glm::vec3(0.0f, 1.0f, 0.0f));
        auto Lpos = glm::vec4(Lpos_polar.ortho(), 1.0f);

        glUniformMatrix4fv(projection_loc, 1, GL_FALSE, glm::value_ptr(projection));
        glUniform4fv(Lpos_loc, 1, glm::value_ptr(view * Lpos));
        glUniform3fv(Lamb_loc, 1, glm::value_ptr(Lamb));
        glUniform3fv(Ldiff_loc, 1, glm::value_ptr(Ldiff));
        glUniform3fv(Lspec_loc, 1, glm::value_ptr(Lspec));
        material.select(0);

        for (int i = 0; i < 256; i++) {
            const Eigen::Vector3d ptcl_pos = md.ptcls_pos.col(i);
            const glm::mat4 model = glm::translate(glm::mat4(1.0f),
                    glm::vec3(ptcl_pos.x(), ptcl_pos.y(), ptcl_pos.z())) * ptcl_scale;

            const glm::mat4 modelview = view * global_trans * global_rot * model;
            const glm::mat3 normalmat = glm::transpose(glm::inverse(glm::mat3(modelview)));

            glUniformMatrix4fv(modelview_loc, 1, GL_FALSE, glm::value_ptr(modelview));
            glUniformMatrix3fv(normalmat_loc, 1, GL_FALSE, glm::value_ptr(normalmat));
            shape->draw();
        }

        ImGui::Render();
        window.swapBuffers();
    }

    // Cleanup
    ImGui_ImplGlfwGL3_Shutdown();
    ImGui::DestroyContext();
    glfwTerminate();

    return 0;
}
