#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>
#include <algorithm>
#include <limits>
#include "csvWriter.cpp"

class Particle {
public:
  std::vector<double> position;
  std::vector<double> velocity;
  std::vector<double> best_position;
  double best_fitness;

  Particle() {
    best_fitness = 1e9;
  }  

  Particle(int dimensions) {
    position.resize(dimensions);
    velocity.resize(dimensions);
    best_position.resize(dimensions);
    best_fitness = 1e9;
  }
};

class PSO {
public:
  int num_particles;
  int num_dimensions;
  double inertia_weight;
  double cognitive_weight;
  double social_weight;
  double a;
  double max_velocity;
  std::vector<Particle> particles;
  std::vector<double> global_best_position;
  double global_best_fitness;

  std::vector<double> realPValues;
  std::vector<double> realQValues;
  int functionCalls;
  int badLocalCounter;
  const double min_value = 0.1;
  const double max_value = 3;

  PSO(int n, int d, double w, double c1, double c2, double a)
    : num_particles(n), num_dimensions(d), inertia_weight(w),
      cognitive_weight(c1), social_weight(c2), a(a) {
    max_velocity = a*(max_value - min_value);
    global_best_position.resize(d);
    global_best_fitness = 1e9;
    functionCalls = 0;
    badLocalCounter = 0;
    loadData();


  }

  void loadData() {

    std::ifstream file("data.txt"); 
    float a, b; 

    while (file >> a >> b) {
        realPValues.emplace_back(a);
        realQValues.emplace_back(b);

    }

    file.close();
}


//   ROSENBROCK: CHANGE DIMENSION IN EXPERIMENT TO 2. SHOULD CHANGE ACCURACY AND NUMBER OF PARTICLES FOR BETTER RESULTS
//   double fitness(std::vector<double> x) {
//     double f = 0.0;
//     for (int i = 0; i < num_dimensions - 1; i++) {
//       double x_i = x[i];
//       double x_i1 = x[i + 1];
//       f += 100.0 * pow(x_i1 - x_i * x_i, 2) + pow(x_i - 1, 2);
//     }
//     functionCalls++;

//     return f;
//   }

  double calculateAbsoluteRelativeError(double prediction, double real_value){
    return abs(prediction - real_value ) / real_value;
}

  double fitness(std::vector<double>& x) {

        double r = x[0];
        double K = x[1];
        double s = x[2];
        double u = x[3];
        double v = x[4];
        double currentP = realPValues[0];
        double currentQ = realQValues[0];
        double deltaP = (r * (1 - realPValues[0] / K) - s*realQValues[0]) * realPValues[0];
        double deltaQ = (-u + v*realPValues[0])*realQValues[0];

        std::vector<double> errorsOfP;
        std::vector<double> errorsOfQ;

        for (int i = 1; i < 101; i++){
            currentP = currentP + deltaP;
            currentQ = currentQ + deltaQ;
            double errorP = fabs(currentP - realPValues[i] ) / realPValues[i];
            double errorQ = fabs(currentQ - realQValues[i] ) / realQValues[i];

            errorsOfP.emplace_back(errorP);
            errorsOfQ.emplace_back(errorQ);

            deltaP = (r * (1 - currentP / K) - s*currentQ) * currentP;
            deltaQ = (-u + v*currentP)*currentQ;
        }

        double maxErrorP = *max_element(errorsOfP.begin(), errorsOfP.end());
        double maxErrorQ = *max_element(errorsOfQ.begin(), errorsOfQ.end());

        errorsOfP.clear();
        errorsOfQ.clear();
        functionCalls++;
        double error = std::max(maxErrorP, maxErrorQ);

        return std::max(maxErrorP, maxErrorQ);

   }

  void initialize() {
    particles.resize(num_particles);
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<double> dist(0.1, 3.0);

    for (int i = 0; i < num_particles; i++) {
      particles[i] = Particle(num_dimensions); 

      for (int j = 0; j < num_dimensions; j++) {
        particles[i].position[j] = dist(generator);
        particles[i].velocity[j] = dist(generator) * max_velocity;
      }

      particles[i].best_position = particles[i].position;
      particles[i].best_fitness = fitness(particles[i].position);

      if (particles[i].best_fitness < global_best_fitness) {
        global_best_fitness = particles[i].best_fitness;
        global_best_position = particles[i].position;

      }
    }
}

  void update() {

    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    for (int i = 0; i < num_particles; i++) {

      Particle& p = particles[i];
      for (int j = 0; j < num_dimensions; j++) {
        double r1 = dist(generator);
        double r2 = dist(generator);
        p.velocity[j] = inertia_weight * p.velocity[j] +
                        cognitive_weight * r1 * (p.best_position[j] - p.position[j]) +
                        social_weight * r2 * (global_best_position[j] - p.position[j]);
        if (p.velocity[j] > max_velocity) {
          p.velocity[j] = max_velocity;
        } else if (p.velocity[j] < -max_velocity) {
          p.velocity[j] = -max_velocity;
        }

        p.position[j] += p.velocity[j];

        // Constrain values to a specific range of values
        if (p.position[j] < min_value) {
          p.position[j] = min_value;
          p.velocity[j] = - p.velocity[j];
        }
        
        if (p.position[j]> max_value) {
          p.position[j] = max_value;
          p.velocity[j] = - p.velocity[j];
        }

      }
    
      double f = fitness(p.position);
      if (f < p.best_fitness) {
        p.best_fitness = f;
        p.best_position = p.position;

        if (f < global_best_fitness) {
          global_best_fitness = f;
          global_best_position = p.best_position;
        }
      }
    }
  }

};



std::vector<double> psoExperiment(double accuracy, int maxFunctionCalls) {

      PSO pso(100, 5, 0.729, 2.05, 2.05, 0.1);
      pso.initialize();
      
      while (pso.global_best_fitness > accuracy && pso.functionCalls < maxFunctionCalls) {
        pso.update();
        pso.badLocalCounter++;
      }
      
      
      std::cout << "Global best position: (";
      for (int j = 0; j < pso.num_dimensions - 1; j++) {
        std::cout << pso.global_best_position[j] << ", ";
      }
      std::cout << pso.global_best_position[pso.num_dimensions - 1] << ")\n";
      std::cout << "Global best fitness: " << pso.global_best_fitness << "\n";
      std::cout << "Function calls: " << pso.functionCalls << "\n";
      
      std::vector<double> results;
      results.emplace_back(pso.global_best_fitness);
      for (auto &x: pso.global_best_position)
        results.emplace_back(x);
      results.emplace_back(pso.functionCalls);
      
      return results;
}

int main() {

  std::string csvFile = "PSO.csv";

  double accuracy = 0.035;
  int maxFunctionCalls = 1000000;
  
  for (int i=0; i<30; i++){
    std::vector<double> res = psoExperiment(accuracy, maxFunctionCalls);
    
    if(!fileExists(csvFile))
        writeCsvFile(csvFile, "fitness", "r", "K", "s", "u", "v", "functionCalls");
    double fitness = res[0];
    double r = res[1];
    double K = res[2];
    double s = res[3];
    double u = res[4];
    double v = res[5];
    int functionCalls = res[6];
    if (!writeCsvFile(csvFile, fitness,r, K, s, u, v, functionCalls)) {
        std::cerr << "Failed to write to file: " << csvFile << "\n";
        }
  }
  return 0;
}

