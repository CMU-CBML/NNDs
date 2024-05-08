#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <cstdlib> // For std::atoi

void generate_job_file(int numThreads, int numNeurons, int job_number, const std::string& identifier, const std::string& variableName, const std::string& value, std::ofstream &batchScript) {
    std::ofstream file;
    std::string filename = identifier + "_" + std::to_string(numNeurons) + "_" + std::to_string(job_number) + ".job";
    file.open(filename);

    file << "#!/bin/sh\n\n";
    file << "#SBATCH -p RM-shared\n";
    file << "#SBATCH --ntasks-per-node=" << numThreads << "\n";
    file << "#SBATCH --job-name=" << identifier << "_" << numNeurons << "_" << job_number << "\n";
    file << "#SBATCH --output=" << identifier << "_" << numNeurons << "_" << job_number << ".out\n";
    file << "#SBATCH --time=48:00:00\n";
    file << "mkdir -p ../io2D_" << identifier << "_" << numNeurons << "_" << job_number << "/outputs/\n";
    // file << "cd ../2D_solver_coarse_05072024/\n";
    file << "/ocean/projects/eng170006p/ussqww/petsc/arch-linux-c-opt/bin/mpiexec -np " << numThreads << " ./2DNG " << numNeurons << " 350000 ../io2D_" << identifier << "_" << numNeurons << "_" << job_number << "/\n\n";
    file << "#-- exit\n#\n";
    file.close();
    std::cout << "Generated file: " << filename << std::endl;

    // Write a batch submission command for the newly created job file to the batch script
    batchScript << "sbatch " << filename << std::endl;

    // Generate simulation_parameters.txt file
    std::map<std::string, std::string> params = {
        {"expandCK_invl", "100"},
        {"var_save_invl", "100"},
        {"numNeuron", "1"},
        {"gc_sz", "2"},
        {"aniso", "6"},
        {"gamma", "1"},
        {"seed_radius", "15"},
        {"kappa", "20"},
        {"dt", "5e-4"},
        {"Dc", "60"},
        {"kp75", "0"},
        {"k2", "0"},
        {"c_opt", "10"},
        {"alpha", "0.9"},
        {"M_phi", "60"},
        {"M_axon", "100"},
        {"M_neurite", "50"},
        {"s_coeff", "0.007"},
        {"delta", "0.20"},
        {"epsilonb", "0.04"},
        {"r", "5"},
        {"g", "0.1"},
        {"alphaT", "0.001"},
        {"betaT", "0.001"},
        {"Diff", "4"},
        {"source_coeff", "15"}
    };

    // Update the specified variable
    params[variableName] = value;

    std::ofstream paramFile("../io2D_" + identifier + "_" + std::to_string(numNeurons) + "_" + std::to_string(job_number) + "/simulation_parameters.txt");
    for (const auto& param : params) {
        paramFile << param.first << " = " << param.second << "\n";
    }
    paramFile.close();
    std::cout << "Generated simulation_parameters.txt in folder: ../io2D_" << identifier << "_" << numNeurons << "_" << job_number << std::endl;
}

int main(int argc, char* argv[]) {
    // if (argc != 7) {
    //     std::cerr << "Usage: " << argv[0] << " [#Threads] [#Neurons] [numCases] [variableName] [variableValue]" << std::endl;
    //     return 1;
    // }

    int numThreads = std::atoi(argv[1]);
    int numNeurons = std::atoi(argv[2]);
    int numCases = std::atoi(argv[3]);
    std::string variableName = argv[4];
    std::string variableValue = argv[5];
    std::string identifier = variableName + "_" + variableValue; // Composite identifier

    std::ofstream batchScript("submit_jobs.sh");
    batchScript << "#!/bin/bash\n\n";

    for (int i = 0; i < numCases; i++) {
        generate_job_file(numThreads, numNeurons, i, identifier, variableName, variableValue, batchScript);
    }

    batchScript.close();
    std::cout << "Generated submit_jobs.sh script for submitting all job files." << std::endl;

    return 0;
}
