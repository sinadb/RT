#pragma once
#define INF 10000000.f
#include "ray.h"
#include "randomstuff.h"
#include "glm.hpp"
#include <algorithm>
#include "textures.h"


static Vec3 Point_Light = Vec3(0, 0, 10);
static Vec3 White = Vec3(1., 1., 1.);

struct Hit_Record {
    bool hit;
    ray r;
    bool collected_light;
    // frame_buffer test
    float z;
    Vec3 Normal;

    Hit_Record() {
        hit = false;
        r.colour() = Vec3(0, 0, 0);
        z = std::numeric_limits<float>::infinity();
        collected_light = false;
    }
    void reset() {
        r = ray();
        hit = false;
        r.colour() = Vec3(0, 0, 0);
        z = std::numeric_limits<float>::infinity();
        collected_light = false;
    }

};

class Hittable
{
public:
    virtual Hit_Record Hit(ray r) { Hit_Record place_holder; return place_holder; };
    virtual Vec3 Centre() { return Vec3(); }
    virtual Vec3 Colour(ray& r) { return Vec3(); }
    virtual float get_index() { return 0; }
    virtual Vec3 Normal(ray& r) { return Vec3(); }
    virtual float radius() { return 0; }
    virtual const Vec3& get_speed() { return Vec3(); }
    virtual Vec3 emit(ray& r) { return Vec3(0, 0, 0); }

};







enum class Material {
    Matte,
    Metal,
    Fuzzy,
    Refractable


};


class Fuzzy {

public:

    static ray init_ray(ray r, float solution, Hittable* object) {
        r.Bounces()++;

        // new origin
        r.origin() += r.direction() * solution;
        // Normal 
        Vec3 N = object->Normal(r);
        //avoid delf reflection
        r.origin() += N * 0.0001f;
        //new reflected ray
        r.direction() = r.direction() - N * (2 * (dot(N, r.direction())));
        // offset vector by the fuzziness factor
        Vec3 offset = Vec3(random_float(), random_float(), random_float());
        while (offset.length() > object->get_index()) {
            offset = Vec3(random_float(), random_float(), random_float());
        }
        r.direction() += offset;
        if (dot(r.direction(), N) > 0) {
            r.colour() = (r.colour()) * (object->Colour(r) * 1.f);
        }
        else {
            r.colour() = Vec3();
        }

        return r;

    }

public:

};

class Refractable {

public:

    static ray init_ray(ray r, float solution, Hittable* object) {
        r.Bounces()++;
        // new origin
        r.origin() += r.direction() * solution;
        // Normal 
        Vec3 N = object->Normal(r);

        //avoid delf reflection
       // r.origin() += N * 0.0001f;
        //avoid delf reflection

        //reverse direction for consistency throughout code . incident light is now facing away from object
        r.direction() = r.direction() * (-1.f);
        r.direction().normalise();

        float test = dot(N, r.direction());
        float rio;
        if (test > 0) {
            rio = 1.f / object->get_index();
            //printf("%f\n", test);
            //std::cout << "p";

        }
        else {

            rio = object->get_index();
            N = N * (-1.f);
            //std::cout << "n";
           // printf("%f\n", test);
        }


        // find out the cos of the angle between incident ray and Normal
        float COS_Theta = dot(r.direction(), N);
        // find the sin now 
        float SIN_Theta = sqrtf(1.f - powf(COS_Theta, 2));
        // create refracted ray but first test that we can refract 
        float Critical = SIN_Theta * rio;
        // only refracting if less than 1
        if (Critical < 1) {
            r.origin() += N * -0.0001f;
            Vec3 R_Para = N * (-sqrtf(1.f - powf(rio * SIN_Theta, 2)));
            Vec3 R_Perp = (r.direction() - N * dot(N, r.direction())) * (-rio);
            Vec3 R = R_Para + R_Perp;
            r.direction() = R;
            r.colour() = Vec3(1, 1, 1) * r.colour();
            //std::cout << "p";

            return r;
        }




        r.origin() += N * 0.0001f;

        //reset Normal if we are not refracting;
        //N = N * (-1.f);


        r.direction() = N * (2 * (dot(N, r.direction()))) - r.direction();
        //r.colour() = (r.colour()) * (object->Colour() * 0.8f);
        r.colour() = Vec3(1, 1, 1) * r.colour();

        return r;


    }
};

class Matte {

public:

    static ray init_ray(ray r, float solution, Hittable* object) {


     
        r.origin() += r.direction() * solution;
        // Normal 
        Vec3 N = object->Normal(r);
        
        
        r.origin() += N * 0.0001f;
        Vec3 Unit_Circle_Centre = r.origin() + N;
        Vec3 Random_p = Vec3(random_float(), random_float(), random_float());
        while (Random_p.length() > 1.f)
        {
            Random_p = Vec3(random_float(), random_float(), random_float());
        }
        Vec3 Random_p_on_U_Circle = Unit_Circle_Centre + Random_p.normalise();
        r.direction() = Random_p_on_U_Circle - r.origin();
        r.colour() =  object->Colour(r) * r.colour() * 0.5;
        r.colour() += r.colour() * object->emit(r);
        //std::cout << r.direction();
        return r;
    }

};

class Metal {


public:


    static ray init_ray(ray r, float solution, Hittable* object) {
        r.Bounces()++;
        // new origin
        r.origin() += r.direction() * solution;
        // Normal 
        Vec3 N = object->Normal(r);
        //avoid delf reflection
        r.origin() += N * 0.0001f;
        r.direction() = r.direction() * (-1.f);
       
        r.direction() = N * (2 * (dot(N, r.direction()))) - r.direction();
        r.colour() = (r.colour()) * (object->Colour(r) * 1.f);
#
        return r;


    }

};



//inline void create_Diffused_ray(ray& r, Vec3& Normal) {
//    Vec3 Unit_Circle_Centre = r.origin() + Normal;
//    Vec3 Random_p = Vec3(random_float(), random_float(), random_float());
//    while (Random_p.length() > 1.f)
//    {
//        Random_p = Vec3(random_float(), random_float(), random_float());
//    }
//    Vec3 Random_p_on_U_Circle = Unit_Circle_Centre + Random_p.normalise();
//    r.direction() = Random_p_on_U_Circle - r.origin();
//
//}

//inline 





class Circles : public Hittable {
private:
    float m_R;
    Vec3 m_Centre;
    Textures* m_texture;
    float m_Albedo;
    Material m_material;
    float m_index;
    Vec3 m_speed;
public:
    Circles(float R, const Vec3& Centre, Textures* tex, Material m, float index, Vec3 speed = Vec3()) :m_R(R), m_Centre(Centre), m_texture(tex), m_material(m), m_index(index), m_speed(speed){}
    Hit_Record Hit(ray r) {

        Hit_Record current;
        current.r = r;
        // check for intersection

        Vec3 oc = r.origin() - (m_Centre + m_speed * r.get_time());
        auto a = dot(r.direction(), r.direction());
        auto b = 2.0 * dot(oc, r.direction());
        auto c = dot(oc, oc) - m_R * m_R;
        auto discriminant = b * b - 4 * a * c;
        // calculate the intersection and make sure it is hitting from outside surface (n.p<0)


        if (discriminant > 0.0f) {


            float solution_2 = (-b + sqrt(discriminant)) / (2 * a);
            float solution_1 = (-b - sqrt(discriminant)) / (2 * a);
            if (solution_1 > 0.0f || solution_2 > 0.0f) {
                float final_solution = std::min(solution_1, solution_2);

                if (final_solution < 0) {
                    final_solution = std::max(solution_1, solution_2);
                }

                switch (m_material) {
                case Material::Matte:
                    current.r = Matte::init_ray(r, final_solution, this);
                    break;
                case Material::Metal:
                    current.r = Metal::init_ray(r, final_solution, this);
                    break;

                case Material::Fuzzy:
                    current.r = Fuzzy::init_ray(r, final_solution, this);
                    break;
                case Material::Refractable:
                    current.r = Refractable::init_ray(r, final_solution, this);
                    break;
                }

               
                current.z = (current.r.origin() - r.origin()).length();
                current.Normal = Normal(r);
                current.hit = true;
                //current.r.where_hit() = this;
                return current;
            }
            else {


                return current;

            }

        }
        else {

            return current;
        }
    }

    Vec3 Centre() { return m_Centre; }
    float Albedo() { return m_Albedo; }
    Vec3 Colour( ray& r) { return m_texture->col_val(r) ; }
    float get_index() { return m_index; }
    Vec3 Normal(ray& r) {
        Vec3 n = (r.origin() - m_Centre + m_speed * r.get_time()).normalise();
        return n;
    }
    float radius() { return m_R; }
    const Vec3& get_speed() { return m_speed; }


};
class Triangle : public Hittable {
private:
    // Anti clock-wise
    Vec3 m_vertices[3];
    Vec3 m_Normal;
    Vec3 m_colour;
    float d;
    Material m_material;
public:
    Triangle(const Vec3& v_0, const Vec3& v_1, const Vec3& v_2, Vec3 colour, Material material) {
        m_vertices[0] = v_0;
        m_vertices[1] = v_1;
        m_vertices[2] = v_2;
        m_colour = colour;
        m_Normal = Cross_Product(m_vertices[1] - m_vertices[0], m_vertices[2] - m_vertices[0]);
        d = -dot(m_Normal, m_vertices[0]);
        m_material = material;

    }
    Hit_Record Hit(ray r) {
        Hit_Record current;
        current.r = r;
        // normal.direction
        float test = dot(m_Normal, r.direction());

        if (test != 0.0) {
            float final_solution = -(d + dot(m_Normal, r.origin())) / test;

            Vec3 hit_point = r.direction() * final_solution;
            Vec3 R = hit_point - m_vertices[0];
            Vec3 Q_1 = m_vertices[1] - m_vertices[0];
            Vec3 Q_2 = m_vertices[2] - m_vertices[0];
            glm::vec2 RQ(dot(R, Q_1), dot(R, Q_2));
            float Q_1_mag = dot(Q_1, Q_1);
            float Q_2_mag = dot(Q_2, Q_2);
            float Q_1_dot_Q_2 = dot(Q_1, Q_2);
            float det = Q_1_mag * Q_2_mag - Q_1_dot_Q_2 * Q_1_dot_Q_2;

            det = 1 / det;

            glm::mat2x2 Matrix({ det * Q_2_mag , det * -Q_1_dot_Q_2 }, { det * -Q_1_dot_Q_2 , det * Q_1_mag });
            glm::vec2 w(Matrix * RQ);
            float w1 = w[0];
            float w2 = w[1];
            /*Vec3 hit_point = r.direction() * final_solution;
            Vec3 R = hit_point - m_vertices[0];
            Vec3 Q_1 = m_vertices[1] - m_vertices[0];
            Vec3 Q_2 = m_vertices[2] - m_vertices[0];
            float w1 = dot(Cross_Product(Q_1, (hit_point - m_vertices[0])), Normal);
            float w2 = dot(Cross_Product((m_vertices[2]-m_vertices[1]), (hit_point - m_vertices[1])), Normal);
            float w3 = dot(Cross_Product((m_vertices[0] - m_vertices[2]), (hit_point - m_vertices[2])), Normal);*/
            //std::cout << w1 << "\t" << w2 << "\n";
            if (w1 >= 0 && w2 >= 0 && w1 <= 1 && w2 <= 1 && 1 - w1 - w2 <= 1 && 1 - w1 - w2 >= 0) {
                switch (m_material) {
                case Material::Matte:
                    current.r = Matte::init_ray(r, final_solution, this);
                    break;
                case Material::Metal:
                    current.r = Metal::init_ray(r, final_solution, this);
                    break;

                case Material::Fuzzy:
                    current.r = Fuzzy::init_ray(r, final_solution, this);
                    break;
                case Material::Refractable:
                    current.r = Refractable::init_ray(r, final_solution, this);
                    break;
                }
                // current.r.oring - r.origin = z

                current.z = (current.r.origin() - r.origin()).length();
                //current.Normal = Normal;
                current.hit = true;
                return current;
            }
            return current;
        }

    }
    virtual Vec3 Centre() { return Vec3(); }
    virtual Vec3 Colour() { return m_colour; }
    virtual float Albedo() { return 1; }
    virtual float fuzzy_factor() { return 1; }
    virtual bool is_refractionable() { return true; }
    virtual float rio() { return 1; }
    Vec3 Normal(Vec3 Intersection) {
        Vec3 N = m_Normal;
        return N.normalise();
    }
    float radius() { return 0; }
};

struct rec {
    Vec3 min;
    Vec3 max;
    void update_box(rec b) {
        for (int i = 0; i != 3; i++) {
            min[i] = std::min(min[i], b.min[i]);
            max[i] = std::max(max[i], b.max[i]);
        }
    }
    bool point_inside(const Vec3& p) {
        if (p[0] < min[0] || p[0] > max[0]) {
            return false;
        }
        if (p[1] < min[1] || p[1] > max[1]) {
            return false;
        }
        if (p[2] < min[2] || p[2] > max[2]) {
            return false;
        }
        return true;
    }
    bool Hit(ray& r) {
        float t_min = -INF;
        float t_max = INF;
        Vec3 inv_Dir = r.inv_dir();
        if (r.direction()[0] != 0) {
            //float denum = 1 / r.direction()[0];
            float t_min_x = (min[0] - r.origin()[0]) * inv_Dir[0];
            float t_max_x = (max[0] - r.origin()[0]) * inv_Dir[0];

            if (t_min_x < 0 && t_max_x < 0) {
                return false;
            }
            if (t_min_x > t_max_x) {
                std::swap(t_min_x, t_max_x);
            }
            t_min = std::max(t_min_x, t_min);
            t_max = std::min(t_max_x, t_max);
            if (t_max < t_min || (t_max < 0 && t_min < 0)) {
                return false;
            }
           
        }
        else {
            if (r.origin()[0] < min[0] || r.origin()[0] > max[0]) {
                return false;
            }
        }
        if (r.direction()[1] != 0) {
            //float denum = 1 / r.direction()[1];
            float t_min_y = (min[1] - r.origin()[1]) * inv_Dir[1];
            float t_max_y = (max[1] - r.origin()[1]) * inv_Dir[1];
            if (t_min_y < 0 && t_max_y < 0) {
                return false;
            }
            if (t_min_y > t_max_y) {
                std::swap(t_min_y, t_max_y);
            }
            
            t_min = std::max(t_min_y, t_min);
            t_max = std::min(t_max_y, t_max);
            if (t_max < t_min || (t_max < 0 && t_min < 0)) {
                return false;
            }
        

        }
        else {
            if (r.origin()[1] < min[1] || r.origin()[1] > max[1]) {
                return false;
            }
        }
        if (r.direction()[2] != 0) {
            //float denum = 1 / r.direction()[2];
            float t_min_z = (min[2] - r.origin()[2]) * inv_Dir[2];
            float t_max_z = (max[2] - r.origin()[2]) * inv_Dir[2];
            if (t_min_z < 0 && t_max_z < 0) {
                return false;
            }
            if (t_min_z > t_max_z) {
                std::swap(t_min_z, t_max_z);
            }
            
            t_min = std::max(t_min_z, t_min);
            t_max = std::min(t_max_z, t_max);
            //std::cout << t_min_z << "\t" << t_max;
            if (t_max < t_min || (t_max < 0 && t_min < 0)) {
                return false;
            }
        }
        else {
            if (r.origin()[2] < min[2] || r.origin()[2] > max[2]) {
                return false;
            }
        }

        return t_max >= t_min;


    }
};

class xz_rect : public Hittable {
private:
    rec m_xz;
    Textures* m_emission;
    Textures* m_colour;
public:
    xz_rect() = default;
    xz_rect(float x0, float x1, float z0, float z1, float y, Textures* colour, Textures* light) {
        m_xz.min[0] = x0; m_xz.max[0] = x1;
        m_xz.min[1] = y - 0.001; m_xz.max[1] = y - 0.001;
        m_xz.min[2] = z0; m_xz.max[2] = z1;
        m_colour = colour;
        m_emission = light;

        // ensure the min and max are in the right order;

    }
    Hit_Record Hit(ray r) {
        Hit_Record current;
        current.r = r;
        float inv_y_dir = r.inv_dir()[1];
        const Vec3 dir = current.r.direction();
        const Vec3 ori = current.r.origin();
        float t = (m_xz.min[1] + 0.001 - ori[1]) * inv_y_dir;
        if (t < 0) {
            return current;
        }
        float hit_x = ori[0] + dir[0] * t;

        if (hit_x < m_xz.min[0] || hit_x > m_xz.max[0]) {
            return current;
        }
        float hit_z = ori[2] + dir[2] * t;
        if (hit_z < m_xz.min[2] || hit_z > m_xz.max[2]) {
            return current;
        }
        //std::cout << 2;
        current.r = Matte::init_ray(r, t, this);
        current.z = (current.r.origin() - r.origin()).length();
        current.hit = true;
       
        if (this->m_emission != nullptr) {
            current.collected_light = true;
        }
        if (this->m_colour != nullptr) {
            current.collected_light = false;
        }
        
        
       
        return current;





    }
    Vec3 Centre() { return Vec3(); };
    Vec3 Colour(ray& r) { 
        if (m_colour!=nullptr) {
            return m_colour->col_val(r);
        }
        return Vec3(1, 1, 1);
    }
    float get_index() { return 1; }
    Vec3 Normal(ray& r) {  
        if (dot(Vec3(0, 1, 0), r.direction()) < 0){
        return Vec3(0, 1, 0); }
        return Vec3(0, -1, 0);
    };
    float radius() { return 0; }
    const Vec3& get_speed() { return Vec3(); }
    Vec3 emit(ray& r) { 
        
        if (m_emission!=nullptr) {
            return m_emission->col_val(r);
        }
 
        return Vec3();
    
    }
    const rec& bounds() { return m_xz; }

};

class xy_rect : public Hittable {
private:
    rec m_xy;
    Textures* m_emission;
    Textures* m_colour;
public:
    xy_rect() = default;
    xy_rect(float x0, float x1, float y0, float y1, float z, Textures* colour, Textures* light) {
        m_xy.min[0] = x0; m_xy.max[0] = x1;
        m_xy.min[1] = y0; m_xy.max[1] = y1;
        m_xy.min[2] = z - 0.001; m_xy.max[2] = z + 0.001;
        m_emission = light;
        m_colour = colour;

        // ensure the min and max are in the right order;

    }
    Hit_Record Hit(ray r) {
        Hit_Record current;
        current.r = r;
        float inv_z_dir = r.inv_dir()[2];
        const Vec3 dir = current.r.direction();
        const Vec3 ori = current.r.origin();
        float t = (m_xy.min[2] + 0.001 - ori[2]) * inv_z_dir;
        if (t < 0) {
            return current;
        }
        float hit_x = ori[0] + dir[0] * t;

        if (hit_x < m_xy.min[0] || hit_x > m_xy.max[0]) {
            return current;
        }
        float hit_y = ori[1] + dir[1] * t;
        if (hit_y < m_xy.min[1] || hit_y > m_xy.max[1]) {
            return current;
        }
        //std::cout << 2;
        current.r = Matte::init_ray(r, t, this);
        current.z = (current.r.origin() - r.origin()).length();
        current.Normal = this->Normal(r);
        current.hit = true;
        if (this->m_emission != nullptr) {
            current.collected_light = true;
        }
        if (this->m_colour != nullptr) {
            current.collected_light = false;
        }
      
        
        return current;





    }
    Vec3 Centre() { return Vec3(); };
    Vec3 Colour(ray& r) {
        if (m_colour!= nullptr) {
            return m_colour->col_val(r);
        }
        return Vec3(1, 1, 1);
    }
    float get_index() { return 1; }
    Vec3 Normal(ray& r) {
        if (dot(Vec3(0, 0, 1), r.direction()) < 0) {
            return Vec3(0, 0, 1);
        }
        return Vec3(0, 0, -1);
    };
    float radius() { return 0; }
    const Vec3& get_speed() { return Vec3(); }
    Vec3 emit(ray& r) {

        if (m_emission!=nullptr) {
            return m_emission->col_val(r);
        }
        return Vec3();


    }
    const rec& bounds() { return m_xy; }

};

class yz_rect : public Hittable {
private:
    rec m_yz;
    Textures* m_emission;
    Textures* m_colour;
public:
    yz_rect() = default;
    yz_rect(float y0, float y1, float z0, float z1, float x, Textures* colour, Textures* light) {
        m_yz.min[0] = x-0.001; m_yz.max[0] = x+0.001;
        m_yz.min[1] = y0; m_yz.max[1] = y1;
        m_yz.min[2] = z0; m_yz.max[2] = z1;
        m_emission = light;
        m_colour = colour;

        // ensure the min and max are in the right order;

    }
    Hit_Record Hit(ray r) {
        Hit_Record current;
        current.r = r;
        float inv_x_dir = r.inv_dir()[0];
        const Vec3 dir = current.r.direction();
        const Vec3 ori = current.r.origin();
        float t = (m_yz.min[0] + 0.001 - ori[0]) * inv_x_dir;
        if (t < 0) {
            return current;
        }
        float hit_y = ori[1] + dir[1] * t;

        if (hit_y < m_yz.min[1] || hit_y > m_yz.max[1]) {
            return current;
        }
        float hit_z = ori[2] + dir[2] * t;
        if (hit_z < m_yz.min[2] || hit_z > m_yz.max[2]) {
            return current;
        }
        
        current.r = Matte::init_ray(r, t, this);
        current.z = (current.r.origin() - r.origin()).length();

        current.hit = true;
        if (this->m_emission != nullptr) {
            current.collected_light = true;
        }
        if(this->m_colour != nullptr) {
            current.collected_light = false;
        }

        return current;





    }
    Vec3 Centre() { return Vec3(); };
    Vec3 Colour(ray& r) {
        if (m_colour!=nullptr) {
            return m_colour->col_val(r);
        }
        return Vec3(1, 1, 1);
    }
    float get_index() { return 1; }
    Vec3 Normal(ray& r) {
        if (dot(Vec3(1, 0, 0), r.direction()) < 0) {
            return Vec3(1, 0, 0);
        }
        return Vec3(-1, 0, 0);
    };
    float radius() { return 0; }
    const Vec3& get_speed() { return Vec3(); }
    Vec3 emit(ray& r) {

        if (m_emission!=nullptr) {
            return m_emission->col_val(r);
        }
        return Vec3();


    }
    const rec& bounds() { return m_yz; }

};

class Box_primirive : public Hittable {
private:
    xy_rect s0, s1;
    yz_rect s2, s3;
    xz_rect s4, s5;
  
public:
    Box_primirive(const Vec3& min, const Vec3& max, Textures* colour, Textures* emission) {
        s0 = xy_rect(min[0], max[0], min[1], max[1], min[2], colour, emission);
        s0 = xy_rect(min[0], max[0], min[1], max[1], max[2], colour, emission);
        s2 = yz_rect(min[1], max[1], min[2], max[2], min[0], colour, emission);
        s3 = yz_rect(min[1], max[1], min[2], max[2], max[0], colour, emission);
        s4 = xz_rect(min[0], max[0], min[2], max[2], min[1], colour, emission);
        s4 = xz_rect(min[0], max[0], min[2], max[2], max[1], colour, emission);
    }
    Hit_Record Hit(ray r) {
        Hit_Record current;
        current.r = r;
        // go through all side and find the hit with smallest intersection
        Hit_Record temp = this->s0.Hit(r);
        if (temp.hit) {
            current = temp;
        }
        temp = this->s1.Hit(r);
        if (temp.hit && temp.z < current.z) {
            current = temp;
        }
        temp = this->s2.Hit(r);
        if (temp.hit && temp.z < current.z) {
            current = temp;
        }
        temp = this->s3.Hit(r);
        if (temp.hit && temp.z < current.z) {
            current = temp;
        }
        temp = this->s4.Hit(r);
        if (temp.hit && temp.z < current.z) {
            current = temp;
        }
        temp = this->s5.Hit(r);
        if (temp.hit && temp.z < current.z) {
            current = temp;
        }
        return current;

    }
};

class Box {
private:
    rec m_bounds;
    Hittable* m_obj;
public:
    Box() { m_bounds.min = Vec3(-INF, -INF, -INF); m_bounds.max = Vec3(INF, INF, INF); }
    Box(Circles* object) {

        m_bounds.min = object->Centre() - object->radius();
        m_bounds.max = (object->Centre() + object->get_speed()) + object->radius();
        m_obj = object;

        // sort my min and max 
        for (int i = 0; i != 3; i++) {
            if (m_bounds.min[i] > m_bounds.max[i]) {
                std::swap(m_bounds.min[i], m_bounds.max[i]);
            }
        }
    }
    Box(xy_rect* object) {
        m_obj = object;
        m_bounds = object->bounds();
    }
    Box(xz_rect* object) {
        m_obj = object;
        m_bounds = object->bounds();
    }
    Box(yz_rect* object) {
        m_obj = object;
        m_bounds = object->bounds();
    }


    Vec3 get_min() { return m_bounds.min; }
    Vec3 get_max() { return m_bounds.max; }
    Hittable* get_obj() { return m_obj; }
    rec get_rec() { return m_bounds; }

};

struct {
    bool operator()( Box& a,  Box& b) {
        if (a.get_min()[0] != b.get_min()[0]) {
            return (a.get_min()[0] < b.get_min()[0]);
        }
        
        else if (a.get_min()[1] != b.get_min()[1]){
            return (a.get_min()[1] < b.get_min()[1]);
        }
        else if (a.get_min()[2] != b.get_min()[2]) {
            return (a.get_min()[2] < b.get_min()[2]);
        }
        return false;
    }
}compare;

struct Node {
public:
    std::vector<Box> world;
    rec m_boundary;
    Node* left;
    Node* right;
    Node* parent;
    Node() {
        //std::cout << "yes";
        for (int i = 0; i != 3; i++) {
            m_boundary.min[i] = INF;
            m_boundary.max[i] = -INF;
        }
    }
    void update_boundary() {
        for (int j = 0; j != 3; j++) {
            for (int i = 0; i != world.size(); i++) {
                m_boundary.min[j] = std::min(m_boundary.min[j], world[i].get_rec().min[j]);
                m_boundary.max[j] = std::max(m_boundary.max[j], world[i].get_rec().max[j]);
            }
        }
    }


};

inline void sort_boxes(std::vector<Box>& data) {
    if (data.size() > 1) {
        std::sort(data.begin(), data.end(), compare);
    }
}

class BVH {
private:
    Node* head;
    std::vector<Hittable*> hitted_objects;
public:
    BVH() { head = nullptr; hitted_objects.reserve(10); }
    Node* populate_world(std::vector<Box>& objects) {
        if (head == nullptr) {
            head = new Node();
            head->world = objects;
            head->update_boundary();
            head->parent = nullptr;
            head->left = nullptr;
            head->right = nullptr;
            //left leaf
        }
        else {
            Node* current = new Node();
            current = head;
            while (current != nullptr) {
                if (compare(objects[0], current->world[current->world.size() / 2]) && current->left != nullptr) {
                    current = current->left;
                }
                else if (compare(objects[0], current->world[current->world.size() / 2]) && current->left == nullptr) {
                    Node* temp = new Node();
                    temp->parent = current;
                    current->left = temp;
                    temp->world = objects;
                    temp->update_boundary();
                    temp->left = temp->right = nullptr;
                    break;
                }
                if (!compare(objects[0], current->world[current->world.size() / 2]) && current->right != nullptr) {
                    current = current->right;

                }
                else if (!compare(objects[0], current->world[current->world.size() / 2]) && current->right == nullptr) {
                    Node* temp = new Node();
                    temp->parent = current;
                    current->right = temp;
                    temp->world = objects;
                    temp->update_boundary();
                    temp->left = temp->right = nullptr;
         
                    break;
                }

            }
        }
        return head;
    }
    Node* get_head() const { return head; }
    std::vector<Hittable*>& get_hittables() { return hitted_objects; }
  
    
     void reset_hittables() { hitted_objects.clear(); }
};
inline void split_vector(std::vector<std::vector<Box>>& result, std::vector<Box>& world) {
    int end = world.size();
    int mid = end / 2;
    if (result.size() == 0) {
        result.push_back(world);
       
    }
    std::vector<Box> left;
    std::vector<Box> right;
    for (int i = 0; i != mid; i++) {
        left.push_back(world[i]);
        

    }
    if (left.size() != 0) {
        result.push_back(left);
        
    }
    if (end != 1) {
        for (int i = mid; i != end; i++) {
            right.push_back(world[i]);


        }
        if (right.size() != 0) {
            result.push_back(right);
            
        }
    }
    
    if (mid > 1 && left.size() != 1) {
       
        split_vector(result, left);
    }
    if (end > 1 && right.size() != 1) {
       
        split_vector(result, right);
    }


};
inline Node* init_BVH(std::vector<Box>& data, BVH& bvh) {
    std::vector<std::vector<Box>> result;
    split_vector(result, data);
    Node* head = new Node();
    for (auto& i : result) {
        head = bvh.populate_world(i);
    }
    //std::cout << result.size();
    return head;

}

inline void Traverse_BT(ray& r, Node* head, BVH& bvh) {

    //Node* temp = new Node();
    
    bool hit_left = false;
    bool hit_right = false;
    bool hit = false;
   
    
        hit = head->m_boundary.Hit(r);
    
    
        if (hit) {
           
            if (head->left != nullptr) {
                if (head->left->m_boundary.Hit(r)) {
                    //std::cout << 1;
                    Traverse_BT(r, head->left, bvh);

                }
            }

            if (head->right != nullptr) {
                if (head->right->m_boundary.Hit(r)) {
                    //std::cout << 2;
                    Traverse_BT(r, head->right, bvh);
                }
            }


            else if (head->left == nullptr && head->right == nullptr) {
                bvh.get_hittables().push_back(head->world[0].get_obj());

            }
        }
    
}




