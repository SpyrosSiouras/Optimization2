#include <algorithm>
#include <functional>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>


using namespace std;


class Point {
    private:
        vector<double> coordinates;

    public:
        Point(size_t size) : coordinates(size) {};
        Point(initializer_list<double> coordinates): coordinates(coordinates.begin(), coordinates.end()) {};
        Point(vector<double> coordinates): coordinates(coordinates.begin(), coordinates.end()) {};

        size_t dimension() const { return coordinates.size(); };

        Point operator=(const Point &other) {
            coordinates.assign(other.coordinates.begin(), other.coordinates.end());
            return *this;
        };

        Point operator+(const Point &other) const {
            assert (dimension() == other.dimension());
            Point result(dimension());
            for (size_t i = 0; i < dimension(); i++)
                result.coordinates[i] = coordinates[i] + other.coordinates[i];
            return result;
        };

        Point operator-(const Point &other) const {
            assert (dimension() == other.dimension());
            Point result(dimension());
            for (size_t i = 0; i < dimension(); i++)
                result.coordinates[i] = coordinates[i] - other.coordinates[i];
            return result;
        };

        Point operator-() const {
            Point zero(dimension());
            return zero - *this;
        };

        Point operator*(const double &number) const {
            Point result(dimension());
            for (size_t i = 0; i < dimension(); i++)
                result.coordinates[i] = number * coordinates[i];
            return result;
        };

        Point operator/(const double &number) const {
            assert (number != 0);
            return *this * (1/number);
        };

        double operator[](int idx) const {
            return coordinates[idx];
        };

        void show() const {
            for (double i: coordinates)
                std::cout << i << ", ";
            cout<<endl;
        };

        static Point get_random_point(size_t dimension, double min, double max) {
            Point random_point(dimension);
            for ( size_t i=0; i<dimension ; i++ ) {
                double r = ((double(rand()) / double(RAND_MAX)) * (max - min)) + min;
                random_point.coordinates[i] = r;
            };
            return random_point;
        };
};

inline Point operator* (const double number, const Point &point){ return point * number; };


class Simplex {
    private:
        multimap<double, Point> vertices;

    public:
        Simplex() {};

        Simplex(vector<pair<double, Point>> vertices) {
            for ( auto pair: vertices ) {
                this -> vertices.emplace(pair);
            };
        };

        Simplex(function<double(Point)> func, vector<Point> vertices) {
            for ( Point point: vertices ) {
                this -> vertices.emplace(make_pair(func(point), point));
            };
        };

        size_t dimension() const { return vertices.size(); };

        size_t space_dimension() const {
            Point first_point = vertices.begin() -> second;
            return first_point.dimension();
        };

        Point centroid_except_worst() const {
            Point sum = Point(space_dimension());
            for (auto itr = vertices.begin(); itr != prev(vertices.end()); itr++)
                sum = sum + itr -> second;
            return sum / (dimension() - 1);
        };

        auto best_vertex_pos() const { return vertices.begin(); };

        auto worse_vertex_pos() const { return prev(vertices.end()); };

        auto second_worse_vertex_pos() const { return prev(vertices.end(), 2); };

    private:
        void erase_worse() { vertices.erase( worse_vertex_pos() ); };

    public:
        void insert(const double functional_value, const Point point) {
            erase_worse();
            vertices.emplace(functional_value, point);
        };

        void insert_best(const double functional_value, const Point &point) {
            erase_worse();
            vertices.emplace_hint(vertices.begin(), functional_value, point);
        };

        void insert_worst(const double functional_value, const Point &point) {
            erase_worse();
            vertices.emplace_hint(vertices.end(), functional_value, point);
        };

        void shrink_towards_best(const function<double(Point)> &func) {
            multimap<double, Point> shrunk_vertices;

            Point const best = best_vertex_pos() -> second;

            shrunk_vertices.emplace(best_vertex_pos() -> first, best);

            for( auto pair = next(vertices.begin()); pair != vertices.end(); ++pair ) {
                Point vertex = pair -> second;
                Point new_vertex = (best + vertex)/2;
                shrunk_vertices.emplace( func(new_vertex), new_vertex );
            };

            vertices = shrunk_vertices;
        };

        void print() {
            cout<<"simplex state:"<<endl;
            for (auto point: vertices){
                cout<<"    point ";
                point.second.show();
                cout<<"    with error = "<<point.first<<endl;
            };
        };
};


class NelderMead {
    private:
        function<double(Point)> objective_function;
        Simplex simplex;
        double min_bound, max_bound;

        const double rho_inc = -0.5;
        const double rho_exc = 0.5;
        const double rho_ref = 1;
        const double rho_exp = 2;

    public:
        Point best_known_minimizer() const { return simplex.best_vertex_pos() -> second; };
        double best_known_minimum() const { return simplex.best_vertex_pos() -> first; };

    private:
        Point worst_vertex() const { return simplex.worse_vertex_pos() -> second; };
        double worst_vertex_value() const { return simplex.worse_vertex_pos() -> first; };

        Point second_worst_vertex() const { return simplex.second_worse_vertex_pos() -> second; };
        double second_worst_vertex_value() const { return simplex.second_worse_vertex_pos() -> first; };

        Point centroid() const { return simplex.centroid_except_worst(); };

        bool is_within_bounds(Point point) const {
            for (size_t i=0; i < point.dimension(); i++)
                if (point[i] < min_bound || point[i] > max_bound)
                    return false;
            return true;
        };

        inline Point candidate_point_for(double rho) const {
            Point point = (1 + rho) * centroid() - rho * worst_vertex();
            if (is_within_bounds(point)) {
                return point;
            } else {
                return worst_vertex();
            }
        };

        inline Point reflection_point() const { return candidate_point_for(rho_ref); };
        inline Point expansion_point() const { return candidate_point_for(rho_exp); };
        inline Point external_contraction_point() const { return candidate_point_for(rho_exc); };
        inline Point internal_contraction_point() const { return candidate_point_for(rho_inc); };

    public:
        NelderMead(function<double(Point)> objective_function, vector<Point> initial_points, double min_bound, double max_bound) {
            this -> objective_function = objective_function;
            this -> simplex = Simplex(objective_function, initial_points);
            this -> min_bound = min_bound;
            this -> max_bound = max_bound;
        };

        Point iterate_once() {
            Point const reflection_point = this -> reflection_point();
            double const f_ref = objective_function(reflection_point);
            double const f_best = best_known_minimum();
            double const f_worst = worst_vertex_value();
            double const f_second_worst = second_worst_vertex_value();

            if (f_ref < f_best) {

                Point new_best_vertex = reflection_point;
                double f_new_best = f_ref;

                Point const expansion_point = this -> expansion_point();
                double const f_exp = objective_function(expansion_point);
                if (f_exp < f_ref) {
                    new_best_vertex = expansion_point;
                    f_new_best = f_exp;
                };

                simplex.insert_best(f_new_best, new_best_vertex);

            } else if (f_ref < f_second_worst) {

                simplex.insert(f_ref, reflection_point);

            } else if (f_ref < f_worst) {

                Point new_worst_vertex = reflection_point;
                double f_new_worst = f_ref;

                Point const external_contraction_point = this -> external_contraction_point();
                double const f_exc = objective_function(external_contraction_point);
                if (f_exc <= f_ref) {
                    new_worst_vertex = external_contraction_point;
                    f_new_worst = f_exc;
                };

                simplex.insert_worst(f_new_worst, new_worst_vertex);
            } else {

                Point const internal_contraction_point = this -> internal_contraction_point();
                double const f_inc = objective_function(internal_contraction_point);

                if (f_inc < f_worst) {
                    simplex.insert_worst(f_inc, internal_contraction_point);
                } else {
                    simplex.shrink_towards_best(objective_function);
                };

            };

            return best_known_minimizer();
        };

        Point iterate(int iterations = 1) {
            int badLocalMinCounter = 0;
            for (int i=0; i < iterations; i++) {
                
                iterate_once();
                
                if  (best_known_minimum() > 0.035) {
                    badLocalMinCounter ++;

                    if (badLocalMinCounter > 250) {

                        Point p = best_known_minimizer();
                        vector<double> v = {0.08, 0, 0, 0, -0.08};
                        Point p1 = p + (v);
                        double const f_inc = objective_function(p1);

                        simplex.insert(f_inc, p1);
                        badLocalMinCounter = 0;

                    }
                }

            }
            return best_known_minimizer();
        };

        void print() { simplex.print(); };
};


class PredatorPreyModel {
    public:
        double r, K, s, u, v;
        const vector<double> P_data, Q_data;
        double error = 0;
        double Pn = P_data[0], Qn = Q_data[0];

    public:
        PredatorPreyModel(vector<double> P_data, vector<double> Q_data) : P_data(P_data), Q_data(Q_data) {};

        double operator()(double r, double K, double s, double u, double v) {
            Pn = P_data[0];
            Qn = Q_data[0];
            error = 0;
            this -> r = r;
            this -> K = K;
            this -> s = s;
            this -> u = u;
            this -> v = v;

            // cout<<"initial state"<<": Pn="<<Pn<<", Pn_real="<<P_data[0]<<", Qn="<<Qn<<", Qn_real="<<Q_data[0]<<"\nerror="<<error<<endl;

            for (size_t i = 1; i < P_data.size(); i++){


                double next_Pn = Pn * (1 + r*(1 - Pn/K) - s*Qn);
                double next_Qn = Qn * (1 - u + v*Pn);

                if ( isnan(next_Pn) || isnan(next_Qn) ) {
                    error = numeric_limits<double>::infinity();
                    break;
                }

                Pn = next_Pn;
                Qn = next_Qn;

                double Pn_real = P_data[i];
                double Qn_real = Q_data[i];

                double abs_rel_error_of_Pn = abs(Pn-Pn_real)/abs(Pn_real);
                double abs_rel_error_of_Qn = abs(Qn-Qn_real)/abs(Qn_real);

                error = max(error, max(abs_rel_error_of_Pn, abs_rel_error_of_Qn));

                // cout<<"\n>>iteration "<<i<<": Pn="<<Pn<<", Pn_real="<<Pn_real<<", Qn="<<Qn<<", Qn_real="<<Qn_real<<"\nerror="<<error<<endl;
            };

            // cout<<"\ninput r="<<r<<", K="<<K<<", s="<<s<<", u="<<u<<", v="<<v<<"\n~~~output  ~~~>error="<<error<<endl;
            // cin.ignore();
            return error;
        };
};



double rosenbrock(double x, double y) {
    return pow(1.0-x,2) + 100.0 * pow((y-pow(x,2)),2);
}

void test_NelderMead_on_rosenbrock(int iterations = 200) {
    double x,y;

    x = 1.00000001;
    y = 1.00000001;
    cout << "\nRosenbrock" << "(" << x << ", " << y << ") = " << rosenbrock(x,y)<<endl;

    cout << "\nRosenbrock" << "(" << -2.06 << ", " << -24.8085 << ") = "<<rosenbrock(-2.06,-24.8085)<<endl<<endl<<endl<<endl<<endl;

    NelderMead algo(
        [](Point point){ return rosenbrock(point[0], point[1]); },
        vector<Point>{
            Point::get_random_point(2, -1000, 1000), //Point{-2.06,    -24.8085},
            Point::get_random_point(2, -1000, 1000), //Point{-59.6179, -91.7051},
            Point::get_random_point(2, -1000, 1000), //Point{-44.8103,  88.8241},
        },
        -10000,
        10000
    );

    cout<<"finding the rosenbrock minimum with random initial simplex vertices"<<endl<<endl;
    for (int i=1; i <= iterations; i++){
        cout<<"iteration "<<i<<"      ";
        algo.iterate_once().show();
    }

}





int main() {
    int const MAX_ITERATIONS = 500000; // change this to a suitable value

    // test_NelderMead_on_rosenbrock();

    // cout<< std::fixed << std::setprecision(64);
    std::srand( ( unsigned int )std::time( nullptr ) );rand();

    auto P_data = vector<double> {1.0000, 0.4982, 0.3432, 0.2815, 0.2946, 0.3597, 0.4772, 0.5915, 0.6239, 0.5527, 0.4510, 0.3666, 0.3231, 0.3439, 0.4091, 0.4859, 0.5604, 0.5442, 0.4893, 0.4178, 0.3612, 0.3487, 0.3787, 0.4364, 0.4970, 0.5243, 0.4984, 0.4436, 0.3930, 0.3679, 0.3779, 0.4075, 0.4718, 0.4941, 0.4963, 0.4645, 0.4261, 0.3870, 0.3863, 0.4100, 0.4408, 0.4724, 0.4839, 0.4706, 0.4479, 0.4126, 0.4001, 0.4073, 0.4224, 0.4524, 0.4655, 0.4787, 0.4521, 0.4287, 0.4116, 0.4074, 0.4170, 0.4380, 0.4742, 0.4700, 0.4561, 0.4431, 0.4231, 0.4097, 0.4219, 0.4316, 0.4501, 0.4586, 0.4488, 0.4456, 0.4310, 0.4171, 0.4194, 0.4316, 0.4377, 0.4550, 0.4546, 0.4520, 0.4416, 0.4265, 0.4132, 0.4305, 0.4355, 0.4402, 0.4481, 0.4487, 0.4357, 0.4326, 0.4260, 0.4287, 0.4319, 0.4355, 0.4472, 0.4478, 0.4403, 0.4425, 0.4313, 0.4306, 0.4372, 0.4370, 0.4323};
    auto Q_data = vector<double> {1.0000, 1.9028, 2.1218, 1.8335, 1.3446, 1.0102, 0.9077, 0.9540, 1.1918, 1.5487, 1.8263, 1.8779, 1.6391, 1.3386, 1.1287, 1.0844, 1.1673, 1.3785, 1.6247, 1.7691, 1.7025, 1.4942, 1.2976, 1.1845, 1.1856, 1.3216, 1.5170, 1.6517, 1.6875, 1.5795, 1.4075, 1.3022, 1.2783, 1.2999, 1.4244, 1.5325, 1.6206, 1.6495, 1.5065, 1.3565, 1.3149, 1.3135, 1.3938, 1.4803, 1.6075, 1.5909, 1.5559, 1.4472, 1.3731, 1.3559, 1.3722, 1.4389, 1.5131, 1.5668, 1.5371, 1.4761, 1.4009, 1.4188, 1.3895, 1.4020, 1.4766, 1.5437, 1.5390, 1.5075, 1.4558, 1.4116, 1.3650, 1.4037, 1.4377, 1.5174, 1.5321, 1.4820, 1.4604, 1.4557, 1.4227, 1.3902, 1.4438, 1.4898, 1.5189, 1.5172, 1.4708, 1.4435, 1.4465, 1.4037, 1.4387, 1.4475, 1.4847, 1.4945, 1.4768, 1.4363, 1.4517, 1.4269, 1.4228, 1.4567, 1.4872, 1.4961, 1.4745, 1.4590, 1.4632, 1.4333, 1.4181};


    for (int trial=0; trial < 1; trial++) {
        auto model = PredatorPreyModel(P_data, Q_data);

        cout<<"\n\n\niteration "<<trial<<endl<<endl;
        vector<Point> initial_points;
        for (int i=0;i<=5;i++){
            Point p = Point::get_random_point(5, 0.3, 2.5);
            initial_points.push_back(p);
            cout<<"    initial random point:   ";
            p.show();
        };

        NelderMead model_nelder_mead (
            [&model](Point point) {
                // cout<<"point (r,K,s,u,v) = ";point.show();
                double val = model(point[0], point[1], point[2], point[3], point[4]);
                // cout<<"  ==>  error = "<<val<<endl<<endl;
                return val;
            },
            initial_points,
            0.1,
            3
        );

        model_nelder_mead.iterate(MAX_ITERATIONS);

        cout<<"\n    after "<<MAX_ITERATIONS<<" iterations, the minimum found was "<<model_nelder_mead.best_known_minimum()<<"\n    at ";
        model_nelder_mead.best_known_minimizer().show();

        // ofstream minimum("./output/execution"+to_string(trial)+"/found_minimum.txt", ofstream::out);
        // minimum << fixed << setprecision(32);
        // minimum << model_nelder_mead.best_known_minimum();
        // minimum.close();

        // ofstream minimizer("./output/execution"+to_string(trial)+"/found_minimizer.txt", ofstream::out);
        // minimizer << fixed << setprecision(32);
        // for (int i=0; i<5; i++)
        //     minimizer << model_nelder_mead.best_known_minimizer()[i] << " ";
        // minimizer.close();

        // ofstream input("./output/execution"+to_string(trial)+"/input.txt", ofstream::out);
        // input << fixed << setprecision(32);
        // for (Point point: initial_points) {
        //     for (auto i=0; i<5; i++)
        //         input << point[i] << " ";
        //     input << endl;
        // };
        // input.close();

        // cout << "\ncontinue?" << endl;
        // cin.ignore();
    };

    return 0;
}