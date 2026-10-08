/*
Author: David Holmqvist <daae19@student.bth.se>
*/

#include "matrix.hpp"
#include "ppm.hpp"
#include <fstream>
#include <stdexcept>
#include <vector>

Matrix::Matrix(unsigned char* R, unsigned char* G, unsigned char* B, unsigned x_size, unsigned y_size, unsigned color_max)
    : R { R }
    , G { G }
    , B { B }
    , x_size { x_size }
    , y_size { y_size }
    , color_max { color_max }
{
}

Matrix::Matrix()
    : Matrix {
        nullptr,
        nullptr,
        nullptr,
        0,
        0,
        0,
    }
{
}

Matrix::Matrix(unsigned dimension)
    : R { new unsigned char[dimension * dimension] }
    , G { new unsigned char[dimension * dimension] }
    , B { new unsigned char[dimension * dimension] }
    , x_size { dimension }
    , y_size { dimension }
    , color_max { 0 }
{
}

Matrix::Matrix(unsigned x_dimension, unsigned y_dimension)
    : R { new unsigned char[x_dimension * y_dimension] }
    , G { new unsigned char[x_dimension * y_dimension] }
    , B { new unsigned char[x_dimension * y_dimension] }
    , x_size { x_dimension }
    , y_size { y_dimension }
    , color_max { 0 }
{
}

Matrix::Matrix(const Matrix& other)
    : R { new unsigned char[other.x_size * other.y_size] }
    , G { new unsigned char[other.x_size * other.y_size] }
    , B { new unsigned char[other.x_size * other.y_size] }
    , x_size { other.x_size }
    , y_size { other.y_size }
    , color_max { other.color_max }
{
    for (auto x { 0 }; x < x_size; x++) {
        for (auto y { 0 }; y < y_size; y++) {
            auto &r_val { r(x, y) }, &g_val { g(x, y) }, &b_val { b(x, y) };
            auto other_r_val { other.r(x, y) }, other_g_val { other.g(x, y) }, other_b_val { other.b(x, y) };

            r_val = other_r_val;
            g_val = other_g_val;
            b_val = other_b_val;
        }
    }
}

Matrix& Matrix::operator=(const Matrix& other) noexcept
{
    if (this == &other) {
        return *this;
    }

    this->~Matrix();

    R = new unsigned char[other.x_size * other.y_size];
    G = new unsigned char[other.x_size * other.y_size];
    B = new unsigned char[other.x_size * other.y_size];

    x_size = other.x_size;
    y_size = other.y_size;
    color_max = other.color_max;

    for (auto x { 0 }; x < x_size; x++) {
        for (auto y { 0 }; y < y_size; y++) {
            auto &r_val { r(x, y) }, &g_val { g(x, y) }, &b_val { b(x, y) };
            auto other_r_val { other.r(x, y) }, other_g_val { other.g(x, y) }, other_b_val { other.b(x, y) };

            r_val = other_r_val;
            g_val = other_g_val;
            b_val = other_b_val;
        }
    }

    return *this;
}

Matrix& Matrix::operator=(Matrix&& other) noexcept
{
    if (this != &other) 
    {
        delete[] R;
        delete[] G;
        delete[] B;

        R = other.R;
        G = other.G;
        B = other.B;

        x_size = other.x_size;
        y_size = other.y_size;
        color_max = other.color_max;

        other.R = nullptr;
        other.G = nullptr;
        other.B = nullptr;
        other.x_size = 0;
        other.y_size = 0;
    }

    return *this;
}

Matrix::~Matrix()
{
    if (R) {
        delete[] R;
        R = nullptr;
    }
    if (G) {
        delete[] G;
        G = nullptr;
    }
    if (B) {
        delete[] B;
        B = nullptr;
    }

    x_size = y_size = color_max = 0;
}

unsigned Matrix::get_x_size() const
{
    return x_size;
}

unsigned Matrix::get_y_size() const
{
    return y_size;
}

unsigned Matrix::get_color_max() const
{
    return color_max;
}

unsigned char* Matrix::get_R()
{
    return R;
}

unsigned char* Matrix::get_G()
{
    return G;
}

unsigned char* Matrix::get_B()
{
    return B;
}

unsigned char Matrix::r(unsigned x, unsigned y) const
{
    return R[y * x_size + x];
}

unsigned char Matrix::g(unsigned x, unsigned y) const
{
    return G[y * x_size + x];
}

unsigned char Matrix::b(unsigned x, unsigned y) const
{
    return B[y * x_size + x];
}

unsigned char& Matrix::r(unsigned x, unsigned y)
{
    return R[y * x_size + x];
}

unsigned char& Matrix::g(unsigned x, unsigned y)
{
    return G[y * x_size + x];
}

unsigned char& Matrix::b(unsigned x, unsigned y)
{
    return B[y * x_size + x];
}

void Matrix::transpose() {    
    const int width = x_size;
    const int height = y_size;
    const int size = width * height;

    std::vector<bool> visited(size, false);

    for (int start = 0; start < size; ++start)
    {
        if (visited[start])
            continue;

        int current = start;

        unsigned char held_r = R[current];
        unsigned char held_g = G[current];
        unsigned char held_b = B[current];

        while (true)
        {
            visited[current] = true;

            int x = current % width;
            int y = current / width;

            int next = x * height + y;

            if (next == start)
                break;

            std::swap(R[next], held_r);
            std::swap(G[next], held_g);
            std::swap(B[next], held_b);

            current = next;
        }

        R[start] = held_r;
        G[start] = held_g;
        B[start] = held_b;
    }

    std::swap(x_size, y_size);
}