#pragma once
#include <GL/gl3w.h>

class Object {
private:
    GLuint vao;
    GLuint vbo;
    GLuint ibo;

    Object(const Object &o);
    Object &operator=(const Object &o);

public:
    struct Vertex {
        GLfloat position[3];
        GLfloat normal[3];
    };

    Object(GLint size, GLsizei vertex_cnt, const Vertex *vertex,
            GLsizei index_cnt=0, const GLuint *index=NULL) {
        glGenVertexArrays(1, &vao);
        glBindVertexArray(vao);

        glGenBuffers(1, &vbo);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(GL_ARRAY_BUFFER,
                vertex_cnt * sizeof(Vertex), vertex, GL_STATIC_DRAW);

        glVertexAttribPointer(0, size, GL_FLOAT, GL_FALSE, sizeof(Vertex), 0);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex),
                static_cast<char *>(0) + sizeof(vertex->position));
        glEnableVertexAttribArray(1);

        glGenBuffers(1, &ibo);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                index_cnt * sizeof(GLuint), index, GL_STATIC_DRAW);
    }

    virtual ~Object() {
        glDeleteBuffers(1, &vao);
        glDeleteBuffers(1, &vbo);
        glDeleteBuffers(1, &ibo);
    }

    void bind() const {
        glBindVertexArray(vao);
    }
};
