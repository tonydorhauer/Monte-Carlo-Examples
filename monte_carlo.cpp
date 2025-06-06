#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <ctime>
#include <vector>
#include <string>
#include <utility>



// --- Function Prototypes --- //

// Generates random numbers uniformly on some interval (defaults to [-1,1])
double generate_random_number(double a = -1.0, double b = 1.0);

// Finds the minimum between two numbers
double find_minimum(double a, double b);

// Approximates pi using the number of points inside a unit circle inscribed in a 2x2 box
double pi_approximation(long int N);

// Approximates e using MC
double e_approximation(long int N);

// Approximates the area under a curve using an MC algorithm
double MC_integral_calculator(double a, double b, long int N);

// Approximate pi using rejection sampling
double rejection_sampling(double a, double b, long int N);

// Applies importance sampling
std::vector<double> importance_sampling(long int N, int initializationSteps);

// Applies the Metropolis-Hastings algorithm
std::vector<double> markov_chain_MC(long int N, int initializationSteps);

// Metropolis algorithm 2D -- change sigmas for non-isotropic Gaussian
std::vector<std::pair<double, double>> metropolis2d(long int N, int initializationSteps, double sigma_x = 1.0, double sigma_y = 1.0);

// Saves the samples from markov_chain_MC() to a csv
void save_samples_to_csv(const std::vector<double>& samples, const std::string& filename);

// Overload: Save 2D samples to CSV: one x,y pair per line
void save_samples_to_csv(const std::vector<std::pair<double, double>>& samples, const std::string& filename);

// Energy function used for simulated annealing
double energy(double x, double y);

// Function for simulated annealing
std::pair<double, double> simulated_annealing(int steps, double T_initial, double T_final);

// Generates a matrix of values +/- 1 for the Ising model (Default is 10 x 10)
std::vector<std::vector<int>> generate_spin_matrix(int matrixSize);

// Calculates the energy for a matrix of spins of +/- 1 --- J is the coupling constant
double ising_energy(double J, const std::vector<std::vector<int>>& spin_matrix);

// Monte Carlo Ising Model (Default coupling constant is 1.0)
std::vector<std::vector<int>> MC_ising_model(long int N, double temperature, int matrixSize = 10, double J = 1.0);

// Compute magnetization (sum of spins)
double magnetization(const std::vector<std::vector<int>>& spin_matrix);

// Save matrix to CSV
void save_spin_matrix(const std::vector<std::vector<int>>& spins, const std::string& filename);

// Path integral Monte Carlo (PIMC) to find average position and ground state energy
double path_integral_MC(int numBeads, long int N, double beta);




// --- Operator Overloads --- //

// Overload the << operator for std::vector<double>
std::ostream& operator<<(std::ostream& os, const std::vector<double>& vec)
{
    os << "[";
    for (size_t i = 0; i < vec.size(); ++i)
    {
        os << vec[i];
        if (i != vec.size() - 1)
            os << ", ";
    }
    os << "]";
    return os;
}

// Overload for printing std::pair<T1, T2>
template <typename T1, typename T2>
std::ostream& operator<<(std::ostream& os, const std::pair<T1, T2>& p)
{
    // Format the pair as (first, second)
    os << "(" << p.first << ", " << p.second << ")";
    return os; // Return the stream to allow chaining (e.g., std::cout << p << "\n")
}

// Overload for printing std::vector<T>
template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec)
{
    os << "["; // Start of vector
    for (size_t i = 0; i < vec.size(); ++i)
    {
        os << vec[i]; // This works recursively for pairs because of the overload above

        // Print comma separator between elements, but not after the last one
        if (i != vec.size() - 1)
            os << ", ";
    }
    os << "]"; // End of vector
    return os; // Return the stream to allow chaining
}




// --- Main --- //

int main()
{
    std::ofstream outfile("magnetization_vs_temperature.csv");
    outfile << "Temperature,Magnetization\n";

    int J = 1.0;
    int matrixSize = 10;
    long int N = 1'000'000;

    for (double T = 1.0; T <= 50.0; T += 1.0)
    {
        auto spins = MC_ising_model(N, T, matrixSize, J);
        double M = magnetization(spins);  // optional: std::abs(magnetization(spins))
        outfile << T << "," << M << "\n";
    }

    outfile.close();

    return 0;
}





// --- Function Definitions --- //

double generate_random_number(double a, double b)
{
    static std::random_device rd;
    static std::mt19937 gen(rd()); // Using static to initialize RNG only once
    std::uniform_real_distribution<> dis(a, b);
    
    return dis(gen);
}

double find_minimum(double a, double b)
{
    if (a < b) return a;
    else return b;
}

double pi_approximation(long int N)
{
    int insideCircle = 0;

    for (int i = 0; i < N; i++)
    {
        double x = generate_random_number();
        double y = generate_random_number();

        if (x*x + y*y < 1)
            insideCircle++;
    }

    double pi_approx = 4.0 * static_cast<double>(insideCircle) / N;

    return pi_approx;
}

double e_approximation(long int N)
{
    double sum = 0.0;
    int count = 0;
    int total = 0;

    for (int i = 0; i < N; i++)
    {
        while (sum < 1)
        {
            double randomNumber = generate_random_number(0.0, 1.0);
            sum += randomNumber;
            count++;
        }

        total += count;  // Add the count to total
        count = 0;       // Reset count for the next trial
        sum = 0;         // Reset sum for the next trial
    }
    
    double e_approx = static_cast<double>(total) / N;

    return e_approx;
}

double MC_integral_calculator(double a, double b, long int N)
{
    double total = 0.0; // Sum of all f(x) values

    for (int i = 0; i < N; i++)
    {
        double x_values = generate_random_number(a, b);
        double f_of_x = std::sqrt(1 - (x_values*x_values)); // Make this whatever function you want to integrate
        total += f_of_x;
    }

    double average_f_of_x = static_cast<double>(total) / N;
    double integral_approx = average_f_of_x * (b-a);

    return integral_approx;
}

double rejection_sampling(double a, double b, long int N)
{
    int accepted = 0; // Counts how many points meet the conditional

    for (int i = 0; i < N; i++)
    {
        double x_values = generate_random_number(a, b);
        double y_values = generate_random_number(a, b);
        double f_of_x = std::sqrt(1 - (x_values*x_values));

        if (y_values < f_of_x)
            accepted++;
    }
    
    double pi_approx = 4.0 * static_cast<double>(accepted) / N;

    return pi_approx;
}

std::vector<double> importance_sampling(long int N, int initializationSteps)
{
    std::vector<double> samples;
    double sum = 0.0;

    for (int i = 0; i < N; i++)
    {
        double u = generate_random_number(0.0, 1.0);
        double x = -std::log(u);       // Inverse transform sampling
        double f_of_x = 1 / (1 + x*x); // Function you wish to integrate
        double p_of_x = std::exp(-x);  // Probability distribution

        sum += f_of_x / p_of_x;
        
        if (i >= initializationSteps)
            samples.push_back(sum / (i+1));
    }
    
    return samples;
}

std::vector<double> markov_chain_MC(long int N, int initializationSteps)
{
    std::vector<double> samples;     // Collect all samples to see the distrbution
    double x = 0.0;                  // Start at x = 0
    
    for (int i = 0; i < N; i++)
    {
        double delta = generate_random_number(-1.0, 1.0);
        double x_new = x + delta;
        double probability = std::exp(-(x*x) / 2); // Probability distribution
        double probability_new = std::exp(-(x_new*x_new) / 2);
        double acceptance_probability = find_minimum(1.0, probability_new / probability); // Probability that we accept x_new

        if (generate_random_number(0.0, 1.0) <= acceptance_probability)
            x = x_new;
        
        if (i >= initializationSteps) // Only appends after 'initializationSteps' trials
            samples.push_back(x);
    }
    
    return samples;
}

std::vector<std::pair<double, double>> metropolis2d(long int N, int initializationSteps, double sigma_x, double sigma_y)
{
    std::vector<std::pair<double, double>> samples;
    double x = 0.0;       // Initial position is (0,0)
    double y = 0.0;
    double epsilon = 0.2; // Vary this until you find best results

    for (int i = 0; i < N; i++)
    {
        double delta_x = generate_random_number(-epsilon, epsilon);
        double delta_y = generate_random_number(-epsilon, epsilon);
        double x_new = x + delta_x;
        double y_new = y + delta_y;
        double probability = std::exp(-0.5 * ((x*x) / (sigma_x*sigma_x) + (y*y) / (sigma_y*sigma_y)));
        double probability_new = std::exp(-0.5 * (((x_new*x_new) / (sigma_x*sigma_x)) + ((y_new*y_new) / (sigma_y*sigma_y))));
        double acceptance_probability = find_minimum(1.0, probability_new / probability);

        if (generate_random_number(0.0, 1.0) <= acceptance_probability)
        {
            x = x_new;
            y = y_new;   
        }

        if (i >= initializationSteps)
            samples.emplace_back(x, y);
    }
    
    return samples;
}

void save_samples_to_csv(const std::vector<double>& samples, const std::string& filename)
{
    std::ofstream file(filename);

    if (!file.is_open())
    {
        std::cerr << "Error opening file for writing: " << filename << "\n";
        return;
    }

    for (const auto& sample : samples)
        file << sample << "\n";

    file.close();
}

void save_samples_to_csv(const std::vector<std::pair<double, double>>& samples, const std::string& filename)
{
    std::ofstream file(filename);

    if (!file.is_open())
    {
        std::cerr << "Error opening file for writing: " << filename << "\n";
        return;
    }

    for (const auto& [x, y] : samples)
        file << x << "," << y << "\n";

    file.close();
}

double energy(double x, double y)
{
    return x*x + y*y + 4*std::sin(5*x) + 4*std::sin(5*y);
}

std::pair<double, double> simulated_annealing(int steps, double T_initial, double T_final)
{
    double x = 0.0;
    double y = 0.0;
    double epsilon = 0.1;  // Perturbation scale

    for (int i = 0; i < steps; ++i)
    {
        double T = T_initial * std::pow(T_final / T_initial, (double)i / steps);

        double x_new = x + generate_random_number(-epsilon, epsilon);
        double y_new = y + generate_random_number(-epsilon, epsilon);

        double E = energy(x, y);
        double E_new = energy(x_new, y_new);
        double delta_E = E_new - E;

        if (delta_E < 0 || generate_random_number(0.0, 1.0) <= std::exp(-delta_E / T))
        {
            x = x_new;
            y = y_new;
        }
    }

    return {x, y};  // Final best estimate
}

std::vector<std::vector<int>> generate_spin_matrix(int matrixSize)
{
    std::vector<std::vector<int>> spins(matrixSize, std::vector<int>(matrixSize));

    for (int i = 0; i < matrixSize; i++)
    {
        for (int j = 0; j < matrixSize; j++)
        {
            double rand = generate_random_number(0.0, 1.0);
            spins[i][j] = (rand < 0.5) ? 1 : -1;
        }
    }

    return spins;
}

double ising_energy(double J, const std::vector<std::vector<int>>& spin_matrix)
{
    int rows = spin_matrix.size();
    int cols = spin_matrix[0].size();
    double total_energy = 0.0;

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            int s = spin_matrix[i][j];
            int right = spin_matrix[i][(j + 1) % cols]; // Periodic boundary conditions
            int down  = spin_matrix[(i + 1) % rows][j];

            total_energy += -J * s * right;
            total_energy += -J * s * down;
        }
    }

    return total_energy;
}

std::vector<std::vector<int>> MC_ising_model(long int N, double temperature, int matrixSize, double J)
{
    std::vector<std::vector<int>> spins = generate_spin_matrix(matrixSize);

    for (long int i = 0; i < N; i++)
    {
        int row = std::rand() % matrixSize;
        int column = std::rand() % matrixSize;

        double E_before = ising_energy(J, spins);   // Calculate energy of old configuration
        spins[row][column] *= -1;                   // Flip spin
        double E_after = ising_energy(J, spins);    // Calculate energy of new configuration
        double delta_E = E_after - E_before;        // Calculate change in energy
        double acceptance_probability = find_minimum(1.0, std::exp(-delta_E / temperature));

        if (generate_random_number(0.0, 1.0) > acceptance_probability)
            spins[row][column] *= -1;               // Revert flip
    }

    return spins;
}

double magnetization(const std::vector<std::vector<int>>& spin_matrix)
{
    int total = 0;

    for (const auto& row : spin_matrix)
        for (int s : row)
            total += s;

    return static_cast<double>(total);
}

void save_spin_matrix(const std::vector<std::vector<int>>& spins, const std::string& filename)
{
    std::ofstream file(filename);
    for (const auto& row : spins)
    {
        for (size_t j = 0; j < row.size(); ++j)
        {
            file << row[j];
            if (j < row.size() - 1)
                file << ",";
        }
        file << "\n";
    }
}

double path_integral_MC(int numBeads, long int N, double beta)
{
    double total = 0.0; // Needed to find average position
    double beta_P = beta / static_cast<double>(numBeads);
    double mass = 1.0;  // Mass of particle
    double hbar = 1.0;  // Natural units
    double omega = 1.0; // Harmonic frequency (for harmonic oscillator)
    
    std::vector<double> position(numBeads, 0.0); // Holds the positions; initialized to zero

    for (long int i = 0; i < N; i++)
    {
        int bead = std::rand() % numBeads; // Picks a random bead
        int next = (bead + 1) % numBeads;  // Index for the next bead
        
        double x = position[bead];         // Grabs a random spot in the position vector
        double x_next = position[next];    // Grabs the next one
        
        double x_new = x + generate_random_number(-0.5, 0.5);  // Use symmetric proposal
        
        double V = 0.5 * mass * omega * omega * x * x; // Potential energy
        double V_new = 0.5 * mass * omega * omega * x_new * x_new;
        
        double kinetic = (mass / (2.0 * beta_P*beta_P * hbar*hbar)) * (x - x_next)*(x - x_next); // Kinetic energy
        double kinetic_new = (mass / (2.0 * beta_P*beta_P * hbar*hbar)) * (x_new - x_next)*(x_new - x_next);
        
        double S = beta_P * V + kinetic; // Action
        double S_new = beta_P * V_new + kinetic_new;

        double delta_S = S_new - S;
        double acceptance_probability = find_minimum(1.0, std::exp(-delta_S));

        if (generate_random_number(0.0, 1.0) <= acceptance_probability)
            position[bead] = x_new;

        total += position[bead];
    }

    double avg_position = total / static_cast<double>(N);
    return avg_position;
}