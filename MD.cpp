//该程序用于读取run.in文件，并将读取到的参数赋值给相应的变量。
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <array>
#include <format>

const int Ns = 100;             // output frequency
const double K_B = 8.617343e-5; // Boltzmann's constant in natural unit
const double TIME_UNIT_CONVERSION = 1.018051e+1; // from natural unit to fs

struct Config {
    //该结构体记录每一帧构型的信息
    int AtomNum{};
    double Box[6]{};   //此处的box有六个参数，分别是xyz方向的晶格常数，与其的一半，这样可以减小计算量
    double pe{};   // 该构型的势能
    std::vector<std::string> element;
    std::vector<double> x{}, y{}, z{}, vx{}, vy{}, vz{}, fx{}, fy{}, fz{}, mass{};
};

void readrun(int& numstep, double& timestep, double& temp, int& dump_thermo ,int& dump_xyz) {
    /*
    该函数使用来读取run.in 文件，并将读取到的参数赋值给相应的变量。
    */
    std::vector<std::string> tokens;
    std::string line;
    std::string token{};
    std::ifstream runin("run.in");
    if (!runin) {
        std::cerr << "Error: Cannot open run.in file" << std::endl;
        exit(1);
    }
    std::cout << "\nReading run.in file...\n" << std::endl;

    while (std::getline(runin, line)) {
        if (line.empty() || line[0] == '#') { continue; }
        std::istringstream iss(line);
        while (iss >> token) {
            tokens.push_back(token);
        }

        if (tokens[0] == "time_step") { timestep = std::stod(tokens[1])/TIME_UNIT_CONVERSION; }
        if (tokens[0] == "run") { numstep = std::stoi(tokens[1]); }
        if (tokens[0] == "dump_thermo") { dump_thermo = std::stoi(tokens[1]); }
        if (tokens[0] == "dump_xyz") { dump_xyz = std::stoi(tokens[1]); }
        if (tokens[0] == "temperature" || tokens[0] == "velocity") { temp = std::stod(tokens[1]); }
        //std::cout << tokens[0] << " = " << tokens[1] << std::endl;
        tokens.clear();
    }
    std::cout << "temperature = " << temp << std::endl;
    std::cout << "time_step = " << timestep << std::endl;
    std::cout << "run = " << numstep << std::endl;
    runin.close();
    std::cout << "\nrun.in file read successfully." << std::endl;
}

enum class State {
    FirstLine,
    SecondLine,
    DataLines
};

Config ReadXYZ(const std::string& filename) {
    std::ifstream xyzfile(filename);
    State state = State::FirstLine ;
    Config config{};
    std::string line{};
    int AtomIndex{ 0 };
    if (!xyzfile) {
        std::cerr << "Error : Cannot open" << filename << std::endl;
        exit(1);
    }
    std::cout << "\nReading " << filename << " file...\n" << std::endl;
    while (std::getline(xyzfile, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        std::istringstream iss{ line };
        std::string token{};
        std::vector<std::string> tokens{};
        while (iss >> token) {
            tokens.push_back(token);
        }
        switch (state) {
        case State::FirstLine:
            config.AtomNum = std::stoi(tokens[0]);
            config.x.resize(config.AtomNum);
            config.y.resize(config.AtomNum);
            config.z.resize(config.AtomNum);
            config.vx.resize(config.AtomNum);
            config.vy.resize(config.AtomNum);
            config.vz.resize(config.AtomNum);
            config.fx.resize(config.AtomNum);
            config.fy.resize(config.AtomNum);
            config.fz.resize(config.AtomNum);
            config.element.resize(config.AtomNum);
            config.mass.resize(config.AtomNum);
            state = State::SecondLine;
            std::cout << "AtomNum = " << config.AtomNum << std::endl;
            break;

        case State::SecondLine:
            config.Box[0] = std::stod(tokens[0]);
            config.Box[1] = std::stod(tokens[1]);
            config.Box[2] = std::stod(tokens[2]);
            config.Box[3] = 0.5 * config.Box[0];
            config.Box[4] = 0.5 * config.Box[1];
            config.Box[5] = 0.5 * config.Box[2];

            state = State::DataLines;
            std::cout << "Box = " << config.Box[0] << " " << config.Box[1] << " " << config.Box[2] << std::endl;
            break;

        case State::DataLines:
            if (AtomIndex < config.AtomNum && tokens.size() >= 5) {
                config.element[AtomIndex] = tokens[0];
                config.x[AtomIndex] = std::stod(tokens[1]);
                config.y[AtomIndex] = std::stod(tokens[2]);
                config.z[AtomIndex] = std::stod(tokens[3]);
                config.mass[AtomIndex] = std::stod(tokens[4]);
                AtomIndex++;
                break;
            }
        }
    }
    xyzfile.close();  
    return config;
}

double getEk(const Config& config) {
    double Ek{ 0 };
    for (int n{ 0 }; n < config.AtomNum; n++) {
        Ek += 0.5 * config.mass[n] * (config.vx[n] * config.vx[n] + config.vy[n] * config.vy[n] + config.vz[n] * config.vz[n]);
    }
    return Ek;
}

void scaleVelocity(const double T0, Config& config) {
    // 该函数使用来矫正初始化后的温度，让其符合初始设置的温度
    double Ek{ getEk(config) };
    double T{ 2 * Ek / (3 * config.AtomNum * K_B) };
    double scaleFactor{ std::sqrt(T0 / T) };
    for (int n{ 0 }; n < config.AtomNum; n++) {
        config.vx[n] *= scaleFactor;
        config.vy[n] *= scaleFactor;
        config.vz[n] *= scaleFactor;
    }
}
void Test_velocity(const Config& config) {
    double vxSum{}, vySum{}, vzSum{};
    
    std::cout << "vx init = " << config.vx[0] << std::endl;
    std::cout << "vy init = " << config.vy[0] << std::endl;
    std::cout << "vz init = " << config.vz[0] << std::endl;
}
void initializeVelocity(const double T0, Config& config) {
    //该函数使用来进行MD速度的初始化,先使其初始质心速度为零，而后将其速度初始化到指定的温度

    double centerOfMassVelocity[3]{};
    double totalMass{};
    for (int n = 0; n < config.AtomNum; n++) {
        totalMass += config.mass[n];
        config.vx[n] = -1 + (2 * rand()) / RAND_MAX;
        config.vy[n] = -1 + (2 * rand()) / RAND_MAX;
        config.vz[n] = -1 + (2 * rand()) / RAND_MAX;
        centerOfMassVelocity[0] += config.vx[n] * config.mass[n];
        centerOfMassVelocity[1] += config.vy[n] * config.mass[n];
        centerOfMassVelocity[2] += config.vz[n] * config.mass[n];
    }
    centerOfMassVelocity[0] /= totalMass;
    centerOfMassVelocity[1] /= totalMass;
    centerOfMassVelocity[2] /= totalMass;
    for (int n = 0; n < config.AtomNum; n++) {
        config.vx[n] -= centerOfMassVelocity[0];
        config.vy[n] -= centerOfMassVelocity[1];
        config.vz[n] -= centerOfMassVelocity[2];
    }
    scaleVelocity(T0, config);
}
void CheckPBC(double length, std::vector<double>& x) {
    // 该函数使用来检查其中的一边是否满足周期性边界条件
    for (double& i : x) {
        if (i < 0) i += length;
        if (i > length) i -= length;
    }
}
void ApplyPBC(Config& config) {
    //该函数使用来检查构型是否满足周期性边界条件
    CheckPBC(config.Box[0], config.x);
    CheckPBC(config.Box[0], config.y);
    CheckPBC(config.Box[0], config.z);
}

void CheckMic(const double& length, const double& halfLength, double& val) {
    if (val < -halfLength) { val += length; }
    if (val > halfLength) { val -= length; }
}
void ApplyMir(const double(&Box)[6], double& x, double& y, double& z) {
    // 使用最小镜像约定
    CheckMic(Box[0], Box[3], x);
    CheckMic(Box[1], Box[4], y);
    CheckMic(Box[2], Box[5], z);
}
void GetForce(Config& config) {
    // 该函数使用来读取构型，并通过LJ势来将其构型中的所有的原子的力计算出来
    const double epsilon = 1.032e-2;
    const double sigma = 3.405;
    const double cutoff = 15.0;
    const double cutoffSquare = cutoff * cutoff;
    const double sigma3 = sigma * sigma * sigma;
    const double sigma6 = sigma3 * sigma3;
    const double sigma12 = sigma6 * sigma6;
    const double e24s6 = 24 * epsilon * sigma6;
    const double e48s12 = 48 * epsilon * sigma12;
    const double e4s6 = 4 * epsilon * sigma6;
    const double e4s12 = 4 * epsilon * sigma12;

    // 先进行各种的初始化，防止出现奇怪的值
    config.pe = 0;
    for (int n = 0; n < config.AtomNum; n++) {
        config.fx[n] = config.fy[n] = config.fz[n] = 0;
    }

    for (int i{ 0 }; i < config.AtomNum - 1; i++) {
        for (int j{ i + 1 }; j < config.AtomNum; j++) {
            double dx{ config.x[j] - config.x[i] };
            double dy{ config.y[j] - config.y[i] };
            double dz{ config.z[j] - config.z[i] };
            ApplyMir(config.Box, dx, dy, dz);
            const double distanceSquare{ dx * dx + dy * dy + dz * dz };
            if (distanceSquare > cutoffSquare) continue;

            const double r2inv{ 1 / distanceSquare };
            const double r4inv{ r2inv * r2inv };
            const double r6inv{ r4inv * r2inv };
            const double r8inv{ r4inv * r4inv };
            const double r12inv{ r6inv * r6inv };
            const double r14inv{ r6inv * r8inv };
            const double f_ij{ e24s6 * r8inv - e48s12 * r14inv };

            config.pe += e4s12 * r12inv - e4s6 * r6inv;
            config.fx[i] += f_ij * dx;
            config.fy[i] += f_ij * dy;
            config.fz[i] += f_ij * dz;
            config.fx[j] -= f_ij * dx;
            config.fy[j] -= f_ij * dy;
            config.fz[j] -= f_ij * dz;
        }

    }

}

void MDintegrate(const bool flag, const double timeStep, Config& config) {
    //该函数使用来进行时间积分(该积分分为上下两部分，因为使用的积分方法中间需要进行一次求力)
    const double halfTimeStep{ 0.5 * timeStep };
    for (int n = 0; n < config.AtomNum; n++) {
        const double mass_inv = 1 / config.mass[n];
        const double ax = config.fx[n] * mass_inv;
        const double ay = config.fy[n] * mass_inv;
        const double az = config.fz[n] * mass_inv;
        config.vx[n] += ax * halfTimeStep;
        config.vy[n] += ay * halfTimeStep;
        config.vz[n] += az * halfTimeStep;
        if (flag) {
            config.x[n] += config.vx[n] * timeStep;
            config.y[n] += config.vy[n] * timeStep;
            config.z[n] += config.vz[n] * timeStep;
        }
    }
}


int main() {
    int numstep{};
    double timestep{};
    double temp{};
    int dump_thermo{};
    int dump_xyz{};
    readrun(numstep, timestep, temp,dump_thermo,dump_xyz);
    Config config = ReadXYZ("xyz.in");
    std::ofstream ofile("thermo.out");
    ofile << std::format("{:<10} {:^15} {:^15} {:^15}\n", "Step", "T", "Ek", "Ep");
    std::ofstream xyzfile("dump.xyz");
    std::ofstream vlog("v.log");
    std::ofstream pos_vlog("pos_v.log");

    const clock_t tStart = clock();
    initializeVelocity(temp, config);
    std::cout << "Start MD" << std::endl;
    for (int step{ 0 }; step <= numstep; ++step) {
        ApplyPBC(config);
        MDintegrate(true, timestep, config);
        GetForce(config);
        MDintegrate(false, timestep, config);
        if (dump_thermo >0) {
            if (step % dump_thermo == 0) {
                double Ek = getEk(config);
                double T = 2 * Ek / (3 * config.AtomNum * K_B);
                ofile << std::format("{:<10} {:>15.10f} {:>15.10f} {:>15.10f}\n", step, T, Ek, config.pe);  
            }
        }
        if (dump_xyz > 0) {
            if (step % dump_xyz == 0) {
                xyzfile << config.AtomNum << std::endl;
                xyzfile << std::format("{:<15.10f}  {:<15.10f}  {:<15.10f}\n",config.Box[0], config.Box[1], config.Box[2]);
                for (int n = 0; n < config.AtomNum; n++) {
                    xyzfile << std::format("{:<5}  {:>15.10f}   {:>15.10f}  {:>15.10f}  {:>15.10f}  {:>15.10f}  {:>15.10f}\n", config.element[n],config.x[n],config.y[n],config.z[n],config.vx[n],config.vy[n],config.vz[n]);
                    vlog << std::format("{:>15.10f}  {:>15.10f}  {:>15.10f}",config.vx[n],config.vy[n],config.vz[n]) << std::endl;
                    pos_vlog << std::format("{:>15.10f}  {:>15.10f}  {:>15.10f}  {:>15.10f}  {:>15.10f}  {:>15.10f}",config.x[n],config.y[n],config.z[n],config.vx[n],config.vy[n],config.vz[n]) << std::endl;
                }
            }
        }
        
    }
    ofile.close();
    xyzfile.close();
    vlog.close();
    clock_t tEnd = clock();
    float tElapsed = float(tEnd - tStart) / CLOCKS_PER_SEC;
    std::cout << "Time elapsed: " << tElapsed << " seconds." << std::endl;

    return 0;
}