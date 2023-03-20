#include <SFML/Graphics.hpp>
#include <iostream>
#include <vector>
#include <chrono>
#include <thread>
#include <boost/multiprecision/cpp_dec_float.hpp>

/* UNITS USED (probably)
* distance/length - 1 fm (10^-15 m)
* velocity - 1 fm/fs (1 m/s)
* acceleration - 1 fm/fs^2 (1 m/s^2)
* mass - 1 e (~9.11 * 10^-31 kg) ? questionable
* energy - 1 eV (1.6 * 10^-19 J)
*/

namespace mp = boost::multiprecision;
using namespace sf;
using namespace std;
using boost::multiprecision::cpp_dec_float_100;

// SCIENTIFIC CONSTANTS
const double PI = 3.141592653589793238463;  // ain't it obvious?
const long double h = 6.582119514e-16;      // planck constant
const double _h = h / (2 * PI);             // reduced planck constant
const double k = 8.9875517873681764e9;      // coulomb's constant
const cpp_dec_float_100 dt = 1e-14;         // delta time
const long long c = 299792458;              // speed of light

// SPACE
const double HEIGHT = 2e-6, WIDTH = 2e-6;
// vector<vector<int>> space(WIDTH, vector<int>(HEIGHT));

// DISPLAY RELATED CONSTANTS
const double DBN = 20.f;            // distance between nodes
const double S_NONE = 5.f;          // size of an empty node
const double S_E = 7.f;             // size of a node containing an electron
const double S_S = 7.f;             // size of a node containing an immobile charge source
const double SHIFT = S_E - S_NONE;  // node shift to compensate nodes' size difference

class Source
{
public:
    vector<double> pos;
    double Q;
    Source(vector<double> _pos, double _Q)
    {
        pos = _pos;
        Q = _Q;
    }
};

class Particle
{
public:
	vector<cpp_dec_float_100> pos;
	vector<cpp_dec_float_100> v;
    cpp_dec_float_100 m;
    cpp_dec_float_100 q;
	double s;
    cpp_dec_float_100 e;
    cpp_dec_float_100 d;
    vector<cpp_dec_float_100> a;
    vector<cpp_dec_float_100> f;
    Color color;
	Particle(vector<cpp_dec_float_100> _pos, vector<cpp_dec_float_100> _v, 
        vector<cpp_dec_float_100> _a, cpp_dec_float_100 _m, cpp_dec_float_100 _q, 
        double _s, cpp_dec_float_100 _e, vector<cpp_dec_float_100> _f,
        Color _color)
	{
		pos = _pos;
		v = _v;
		m = _m;
		q = _q;
		s = _s;
		e = _e;
        d = 0;
        a = _a;
        f = _f;
        color = _color;
	}
    void update_forces(vector<Source> s)
    {
        for (int i = 0; i < s.size(); i++)
        {
            cpp_dec_float_100 dx = pos[0] - s[i].pos[0],
                dy = pos[1] - s[i].pos[1];
            cpp_dec_float_100 r = sqrt((dx * dx) + (dy * dy));
            cpp_dec_float_100 sf = k * q * s[i].Q / (r * r);
            cpp_dec_float_100 ux = dx / r,
                uy = dy / r;
            f[0] += sf * ux;
            f[1] += sf * uy;
        }
    }
    void move()
    {
        // this shit is called verlet's integration
        cpp_dec_float_100 x_new = pos[0] + v[0] * dt + 0.5 * a[0] * dt * dt,
            y_new = pos[1] + v[1] * dt + 0.5 * a[1] * dt * dt;
        cpp_dec_float_100 ax_new = f[0] / m,
            ay_new = f[1] / m;
        cpp_dec_float_100 vx_new = v[0] + 0.5 * (a[0] + ax_new) * dt,
            vy_new = v[1] + 0.5 * (a[1] + ay_new) * dt;
        if (vx_new > c)
        {
            vx_new = c;
        }
        else if (vx_new < -c)
        {
            vx_new = -c;
        }
        if (vy_new > c)
        {
            vy_new = c;
        }
        else if (vy_new < -c)
        {
            vy_new = -c;
        }
        if (x_new >= WIDTH)
        {
            x_new = WIDTH;
            vx_new = -abs(vx_new);
            ax_new = -abs(ax_new);
        }
        else if (x_new <= 0)
        {
            x_new = 0;
            vx_new = abs(vx_new);
            ax_new = abs(ax_new);
        }
        if (y_new >= HEIGHT)
        {
            y_new = HEIGHT;
            vy_new = -abs(vy_new);
            ay_new = -abs(ay_new);
        }
        else if (y_new <= 0)
        {
            y_new = 0;
            vy_new = abs(vy_new);
            ay_new = abs(ay_new);
        }
        pos[0] = x_new;
        pos[1] = y_new;
        v = { vx_new, vy_new };
        a = { ax_new, ay_new };
    }
    void log()
    {
        cout << "position: (" << pos[0] << ';' << pos[1] << ")\n";
        cout << "velocity: (" << v[0] << ';' << v[1] << ")\n";
        cout << "acceleration: (" << a[0] << ';' << a[1] << ")\n";
        cout << "resulting force: (" << f[0] << ';' << f[1] << ")\n";
    }
};

class Electron : public Particle {
public:
    Electron(vector<cpp_dec_float_100> _pos, vector<cpp_dec_float_100> _v,
        vector<cpp_dec_float_100> _a = { 0, 0 }, vector<cpp_dec_float_100> _f = { 0, 0 }) :
        Particle(_pos, _v, _a, 9.10938356E-31, -1.602176634E-19, 0.5, 0.0, _f,
            Color(173, 216, 230, 170)) {}
};

class Positron : public Particle {
public:
    Positron(vector<cpp_dec_float_100> _pos, vector<cpp_dec_float_100> _v,
        vector<cpp_dec_float_100> _a = { 0, 0 }, vector<cpp_dec_float_100> _f = { 0, 0 }) : 
        Particle(_pos, _v, _a, 9.10938356e-31, 1.6021766208e-19, 0.5, 0.0, _f,
            Color(0, 0, 139, 170)) {}
};

class ParticleSystem
{
public:
    vector<Particle> particles;
    ParticleSystem(vector<Particle> _particles)
    {
        particles = _particles;
    }
    /*void update_forces(vector<Source> s)
    {
        for (int i = 0; i < s.size(); i++)
        {
            cpp_dec_float_100 dx = pos[0] - s[i].pos[0],
                dy = pos[1] - s[i].pos[1];
            cpp_dec_float_100 r = sqrt((dx * dx) + (dy * dy));
            cpp_dec_float_100 sf = k * q * s[i].Q / (r * r);
            cpp_dec_float_100 ux = dx / r,
                uy = dy / r;
            f[0] += sf * ux;
            f[1] += sf * uy;
        }
    }*/
    void update_forces()
    {
        for (int i = 0; i < particles.size(); i++)
        {
            for (int j = 0; j < particles.size(); j++)
            {
                if (i != j)
                {
                    cpp_dec_float_100 dx = particles[i].pos[0] - particles[j].pos[0],
                        dy = particles[i].pos[1] - particles[j].pos[1];
                    cpp_dec_float_100 r = sqrt((dx * dx) + (dy * dy));
                    cpp_dec_float_100 sf = k * particles[i].q * particles[j].q / (r * r);
                    cpp_dec_float_100 ux = dx / r,
                        uy = dy / r;
                    particles[i].f[0] += sf * ux;
                    particles[i].f[1] += sf * uy;
                }
            }
        }
    }
    void process()
    {
        for (int i = 0; i < particles.size(); i++)
        {
            particles[i].move();
        }
    }
};

int main()
{
    cout << "STARTED" << endl;
    Electron x({ 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, 0 });
    Positron y({ WIDTH, HEIGHT }, { 0, 0 }, { 0, 0 }, { 0, 0 });
    Source s({ 1e-6, 0 }, -1e-23);
    Source s1({ 1e-6, 2e-6 }, 5e-21);
    sf::ContextSettings settings;
    settings.antialiasingLevel = 8;
    RenderWindow window(VideoMode(WIDTH * 205.f * 1e6, HEIGHT * 205.f * 1e6), 
        L"EPsim", Style::Default, settings);
    window.setFramerateLimit(60);
    window.setVerticalSyncEnabled(true);
    cpp_dec_float_100 mindist = 1;
    ParticleSystem ps({ x, y });
    while (window.isOpen())
    {
        Event event;
        while (window.pollEvent(event))
        {
            if (event.type == Event::Closed)
                window.close();
        }
        window.clear(Color::White);
        ps.update_forces();
        ps.process();
        Particle p1 = ps.particles[0];
        Particle p2 = ps.particles[1];
        cpp_dec_float_100 dist = sqrt((p1.pos[0] - p2.pos[0]) * (p1.pos[0] - p2.pos[0]) +
            (p1.pos[1] - p2.pos[1]) * (p1.pos[1] - p2.pos[1]));
        if (dist < mindist)
        {
            mindist = dist;
            cout << "current mindist: " << mindist << endl;
        }
        for (int i = 0; i < ps.particles.size(); i++)
        {
            Particle p = ps.particles[i];
            CircleShape shape(S_E);
            shape.setPosition(double(p.pos[0] * DBN * 1e7 - SHIFT),
                double(p.pos[1] * DBN * 1e7 - SHIFT));
            shape.setFillColor(p.color);
            window.draw(shape);
        }
        window.display();
        // this_thread::sleep_for(std::chrono::milliseconds(100));
    }
    return 0;
}