#pragma once
#include <memory>

#include "Object.h"

class Shape {
private:
    std::shared_ptr<const Object> object;

protected:
    const GLsizei vertex_cnt;

public:
    Shape(GLint size, GLsizei vertex_cnt, const Object::Vertex *vertex,
            GLsizei index_cnt=0, const GLuint *index=NULL)
    : object(new Object(size, vertex_cnt, vertex, index_cnt, index)),
            vertex_cnt(vertex_cnt) {
    }

    void draw() const {
        object->bind();
        execute();
    }

    virtual void execute() const {
        glDrawArrays(GL_LINE_LOOP, 0, vertex_cnt);
    }
};


class ShapeIndex : public Shape {
protected:
    const GLsizei index_cnt;

public:
    ShapeIndex(GLint size, GLsizei vertex_cnt, const Object::Vertex *vertex,
            GLsizei index_cnt, const GLuint *index)
    : Shape(size, vertex_cnt, vertex, index_cnt, index), index_cnt(index_cnt) {
    }

    virtual void execute() const {
        glDrawElements(GL_LINES, index_cnt, GL_UNSIGNED_INT, 0);
    }
};


class SolidShape : public Shape {
public:
    SolidShape(GLint size, GLsizei vertex_cnt, const Object::Vertex *vertex)
    : Shape(size, vertex_cnt, vertex) {
    }

    virtual void execute() const {
        glDrawArrays(GL_TRIANGLES, 0, vertex_cnt);
    }
};


class SolidShapeIndex : public ShapeIndex {
public:
    SolidShapeIndex(GLint size, GLsizei vertex_cnt, const Object::Vertex *vertex,
            GLsizei index_cnt, const GLuint *index)
    : ShapeIndex(size, vertex_cnt, vertex, index_cnt, index) {
    }

    virtual void execute() const {
        glDrawElements(GL_TRIANGLES, index_cnt, GL_UNSIGNED_INT, 0);
    }
};
