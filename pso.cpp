#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>
#include <algorithm>


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
  double max_velocity;
  std::vector<Particle> particles;
  std::vector<double> global_best_position;
  double global_best_fitness;

  std::vector<double> realPValues;
  std::vector<double> realQValues;
  int functionCalls;
  const double min_value = 0.1;
  const double max_value = 3.0;

  PSO(int n, int d, double w, double c1, double c2, double v_max)
    : num_particles(n), num_dimensions(d), inertia_weight(w),
      cognitive_weight(c1), social_weight(c2), max_velocity(v_max) {
    global_best_position.resize(d);
    global_best_fitness = 1e9;
    loadData();
    functionCalls = 0;
    for (auto& x: realPValues)
     std::cout<<x << " ";
  }

  void loadData() {

    std::ifstream file("data.txt");  // Open the input file
    float a, b;                       // Variables to store the data

    // Read data from the file
    while (file >> a >> b) {
        // Do something with the data
        realPValues.emplace_back(a);
        realQValues.emplace_back(b);
    }

    file.close();  // Close the file
}

  double calculateAbsoluteRelativeError(double prediction, double real_value){
    return abs(prediction - real_value ) / real_value;
}

//   double fitness(std::vector<double> x) {
//     double f = 0.0;
//     for (int i = 0; i < num_dimensions - 1; i++) {
//       double x_i = x[i];
//       double x_i1 = x[i + 1];
//       f += 100.0 * pow(x_i1 - x_i * x_i, 2) + pow(x_i - 1, 2);
//     }
//     return f;
//   }

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
            double errorP = calculateAbsoluteRelativeError(currentP, realPValues[i]);
            std::cout<<errorP<<std::endl;
            double errorQ = calculateAbsoluteRelativeError(currentQ, realQValues[i]);

            errorsOfP.emplace_back(errorP);
            errorsOfQ.emplace_back(errorQ);

            deltaP = (r * (1 - currentP / K) - s*currentQ) * currentP;
            deltaQ = (-u + v*currentP)*currentQ;
        }

        double maxErrorP = *max_element(errorsOfP.begin(), errorsOfP.end());
        double maxErrorQ = *max_element(errorsOfQ.begin(), errorsOfQ.end());

        // if (maxErrorP == 0) {
        //     for(auto& e: errorsOfP)
        //         std::cout << e << " ";
        // }
        // std::cout << "MaxErrorP: " << maxErrorP << std::endl;
        // std::cout << "MaxErrorQ: " << maxErrorQ << std::endl;

        errorsOfP.clear();
        errorsOfQ.clear();
        // cout<< maxErrorP << maxErrorQ << endl;
        functionCalls++;
        return std::max(maxErrorP, maxErrorQ);


  }

  void initialize() {
    particles.resize(num_particles);
    std::mt19937 generator(0);
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
    std::mt19937 generator(0);
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
        }
        if (p.position[j]> max_value) {
          p.position[j] = max_value;
        }
        // std::cout<< "Best fitness so far: "<< global_best_fitness << std::endl;

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




int main() {
  PSO pso(100, 5, 0.7298, 1.49618, 1.49618, 2.5);
  pso.initialize();
  for (int i = 0; i < 1000; i++) {
    pso.update();

  }

  
  std::cout << "Global best position: (";
  for (int j = 0; j < pso.num_dimensions - 1; j++) {
    std::cout << pso.global_best_position[j] << ", ";
  }
  std::cout << pso.global_best_position[pso.num_dimensions - 1] << ")\n";
  std::cout << "Global best fitness: " << pso.global_best_fitness << "\n";
  return 0;
}

