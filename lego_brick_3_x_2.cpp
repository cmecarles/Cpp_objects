#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdint>

// 1. DATA STRUCTURE
struct Vertex {
    float px, py, pz;
    float nx, ny, nz;
    float u, v;
};

// 2. GEOMETRY HELPERS
static void addBox(std::vector<Vertex>& V, std::vector<uint32_t>& I,
                   float x0, float y0, float z0,
                   float x1, float y1, float z1)
{
    auto pushFace = [&](Vertex a, Vertex b, Vertex c, Vertex d) {
        uint32_t base = (uint32_t)V.size();
        V.push_back(a); V.push_back(b); V.push_back(c); V.push_back(d);
        I.push_back(base+0); I.push_back(base+1); I.push_back(base+2);
        I.push_back(base+0); I.push_back(base+2); I.push_back(base+3);
    };

    float X0=x0, X1=x1, Y0=y0, Y1=y1, Z0=z0, Z1=z1;

    pushFace({X0,Y0,Z1, 0,0,1, 0,0}, {X1,Y0,Z1, 0,0,1, 1,0}, {X1,Y1,Z1, 0,0,1, 1,1}, {X0,Y1,Z1, 0,0,1, 0,1}); // +Z
    pushFace({X0,Y1,Z0, 0,0,-1, 0,0}, {X1,Y1,Z0, 0,0,-1, 1,0}, {X1,Y0,Z0, 0,0,-1, 1,1}, {X0,Y0,Z0, 0,0,-1, 0,1}); // -Z
    pushFace({X1,Y0,Z0, 1,0,0, 0,0}, {X1,Y1,Z0, 1,0,0, 1,0}, {X1,Y1,Z1, 1,0,0, 1,1}, {X1,Y0,Z1, 1,0,0, 0,1}); // +X
    pushFace({X0,Y1,Z0, -1,0,0, 0,0}, {X0,Y0,Z0, -1,0,0, 1,0}, {X0,Y0,Z1, -1,0,0, 1,1}, {X0,Y1,Z1, -1,0,0, 0,1}); // -X
    pushFace({X1,Y1,Z0, 0,1,0, 0,0}, {X0,Y1,Z0, 0,1,0, 1,0}, {X0,Y1,Z1, 0,1,0, 1,1}, {X1,Y1,Z1, 0,1,0, 0,1}); // +Y
    pushFace({X0,Y0,Z0, 0,-1,0, 0,0}, {X1,Y0,Z0, 0,-1,0, 1,0}, {X1,Y0,Z1, 0,-1,0, 1,1}, {X0,Y0,Z1, 0,-1,0, 0,1}); // -Y
}

static void addDisk(std::vector<Vertex>& V, std::vector<uint32_t>& I,
                    float cx, float cy, float z, float radius, int segments, bool up) 
{
    uint32_t centerIdx = (uint32_t)V.size();
    float nz = up ? 1.0f : -1.0f;
    V.push_back({cx, cy, z, 0, 0, nz, 0.5f, 0.5f});

    for (int k = 0; k <= segments; k++) {
        float t = (float)k / (float)segments;
        float ang = t * 6.28318530718f;
        float x = cosf(ang), y = sinf(ang);
        V.push_back({cx + radius * x, cy + radius * y, z, 0, 0, nz, x * 0.5f + 0.5f, y * 0.5f + 0.5f});
        
        if (k > 0) {
            uint32_t curr = centerIdx + k + 1;
            uint32_t prev = centerIdx + k;
            if (up) { I.push_back(centerIdx); I.push_back(prev); I.push_back(curr); }
            else { I.push_back(centerIdx); I.push_back(curr); I.push_back(prev); }
        }
    }
}

static void addCylinder(std::vector<Vertex>& V, std::vector<uint32_t>& I,
                        float cx, float cy, float z0,
                        float radius, float height, int segments)
{
    addDisk(V, I, cx, cy, z0 + height, radius, segments, true);

    uint32_t base = (uint32_t)V.size();
    for (int k=0; k<=segments; k++) {
        float t = (float)k / (float)segments;
        float ang = t * 6.28318530718f;
        float x = cosf(ang), y = sinf(ang);
        V.push_back({cx + radius*x, cy + radius*y, z0, x, y, 0, t, 0});
        V.push_back({cx + radius*x, cy + radius*y, z0 + height, x, y, 0, t, 1});
    }

    for (int k=0; k<segments; k++) {
        uint32_t i0 = base + 2*k, i1 = base + 2*k + 1, i2 = base + 2*k + 2, i3 = base + 2*k + 3;
        I.push_back(i0); I.push_back(i2); I.push_back(i1);
        I.push_back(i1); I.push_back(i2); I.push_back(i3);
    }
}

static void addInnerTubes(std::vector<Vertex>& V, std::vector<uint32_t>& I, float H) {
    const float tubeR = 6.51f * 0.5f;
    const int seg = 32;
    float tubeX[2] = {-4.0f, 4.0f};
    
    for (int i = 0; i < 2; i++) {
        addCylinder(V, I, tubeX[i], 0.0f, 0.0f, tubeR, H - 1.2f, seg);
    }
}

// 3. BRICK BUILDER
static void buildLego3x2(std::vector<Vertex>& V, std::vector<uint32_t>& I)
{
    const float pitch = 8.0f;
    const int studsX = 3, studsY = 2;
    const float L = studsX * pitch, W = studsY * pitch;
    const float H = 9.6f, studH = 1.8f, studR = 2.4f;
    const int seg = 32;

    addBox(V, I, -L*0.5f, -W*0.5f, 0.0f, +L*0.5f, +W*0.5f, H);

    float shell = 1.2f; 
    addBox(V, I, +L*0.5f - shell, +W*0.5f - shell, 0.0f, 
                -L*0.5f + shell, -W*0.5f + shell, H - shell);

    for (int ix=0; ix<studsX; ix++) {
        float x = (ix - (studsX-1)*0.5f) * pitch;
        for (int iy=0; iy<studsY; iy++) {
            float y = (iy - (studsY-1)*0.5f) * pitch;
            addCylinder(V, I, x, y, H, studR, studH, seg);
        }
    }
    addInnerTubes(V, I, H);
}

// 4. MAIN ENTRY
int main() {
    std::vector<Vertex> vertices;
    std::vector<uint32_t> indices;

    buildLego3x2(vertices, indices);

    std::ofstream file("lego_brick_3_x_2.obj");
    if (!file) return 1;

    for (const auto& v : vertices) {
        file << "v " << v.px << " " << v.py << " " << v.pz << "\n";
        file << "vn " << v.nx << " " << v.ny << " " << v.nz << "\n";
    }
    for (size_t i = 0; i < indices.size(); i += 3) {
        file << "f " << indices[i]+1 << "//" << indices[i]+1 << " "
             << indices[i+1]+1 << "//" << indices[i+1]+1 << " "
             << indices[i+2]+1 << "//" << indices[i+2]+1 << "\n";
    }
    return 0;
}
