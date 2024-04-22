#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib> // For std::atoi

void generate_job_file(int numThreads, int numNeurons, int job_number, std::ofstream &batchScript) {
    std::ofstream file;
    std::string filename = "Single" + std::to_string(job_number) + ".job";
    file.open(filename);

    file << "#!/bin/sh\n\n";
    file << "#SBATCH -p RM-shared\n";
    file << "#SBATCH --ntasks-per-node=" << numThreads << "\n";
    file << "#SBATCH --job-name=" << numNeurons << "_" << job_number << "\n";
    file << "#SBATCH --output=" << numNeurons << "_" << job_number << ".out\n";
    file << "#SBATCH --time=48:00:00\n";
    // file << "#SBATCH --mail-type=ALL\n";
    // file << "#SBATCH --mail-user=your_email@domain.com\n\n";
    file << "mkdir ../io2D_" << numNeurons << "_" << job_number << "/\n";
    file << "mkdir ../io2D_" << numNeurons << "_" << job_number  << "/outputs/\n";
    file << "cd ../2D_solver_coarse_smallGC_04182024/\n";
    file << "/ocean/projects/eng170006p/ussqww/petsc/arch-linux-c-opt/bin/mpiexec -np " << numThreads << "./2DNG " << numNeurons << " 350000 ../io2D_" << numNeurons << "_" << job_number << "/\n\n";
    file << "#-- exit\n#\n";

    file.close();
    std::cout << "Generated file: " << filename << std::endl;

    // Write a batch submission command for the newly created job file to the batch script
    batchScript << "sbatch " << filename << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " [end number]" << std::endl;
        return 1;
    }

    std::ofstream batchScript("submit_jobs.sh");  // File to hold batch submission commands
    batchScript << "#!/bin/bash\n\n"; // Start the batch script with a shebang

    int numThreads = std::atoi(argv[1]);  // Convert command-line input to integer
    int numNeurons = std::atoi(argv[2]);  // Convert command-line input to integer
    int end = std::atoi(argv[3]);  // Convert command-line input to integer

    for (int i = 0; i <= end; i++) {
        generate_job_file(numThreads, numNeurons, i, batchScript);
    }

    batchScript.close(); // Close the batch script file
    std::cout << "Generated submit_jobs.sh script for submitting all job files." << std::endl;

    return 0;
}