//
// Created by ozdalkiran-l on 1/14/26.
//

#include <omp.h>

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
extern "C" {
#include <lime.h>
}
#include "io.h"

namespace fs = std::filesystem;

// Saves a vector of doubles in ../data/filename/filename_plaquette.txt
void io::save_plaquette(const std::vector<double>& data, const std::string& filename,
                        const std::string& dirpath, int precision) {
    // Create a data folder if doesn't exists
    fs::path base_dir(dirpath);
    fs::path dir = base_dir / filename;

    try {
        if (!fs::exists(dir)) {
            fs::create_directories(dir);
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Couldn't create folder data : " << e.what() << std::endl;
        return;
    }
    fs::path filepath = dir / (filename + "_plaquette.txt");

    std::ofstream file(filepath, std::ios::out | std::ios::app);
    if (!file.is_open()) {
        std::cout << "Could not open file " << filepath << "\n";
        return;
    }
    file << std::fixed << std::setprecision(precision);
    for (const double& x : data) {
        file << x << "\n";
    }
    file.close();
    std::cout << "Plaquette written in " << filepath << "\n";
}

void io::save_event_nb(const std::vector<size_t>& event_nb, const std::string& filename,
                       const std::string& dirpath) {
    // Create a data folder if doesn't exists
    fs::path base_dir(dirpath);
    fs::path dir = base_dir / filename;

    try {
        if (!fs::exists(dir)) {
            fs::create_directories(dir);
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Couldn't create folder data : " << e.what() << std::endl;
        return;
    }
    fs::path filepath = dir / (filename + "_nbevent.txt");

    std::ofstream file(filepath, std::ios::out | std::ios::app);
    if (!file.is_open()) {
        std::cout << "Could not open file " << filepath << "\n";
        return;
    }
    for (const size_t& x : event_nb) {
        file << x << "\n";
    }
    file.close();
    std::cout << "Number of events written in " << filepath << "\n";
}

void io::save_event_nb(const std::vector<size_t>& event_nb, const std::vector<size_t>& lift_nb,
                       const std::vector<size_t> rev_nb, const std::vector<double> lambda,
                       const std::string& filename, const std::string& dirpath) {
    // Create a data folder if doesn't exists
    fs::path base_dir(dirpath);
    fs::path dir = base_dir / filename;

    try {
        if (!fs::exists(dir)) {
            fs::create_directories(dir);
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Couldn't create folder data : " << e.what() << std::endl;
        return;
    }
    fs::path filepath = dir / (filename + "_nbevent.txt");

    std::ofstream file(filepath, std::ios::out | std::ios::app);
    if (!file.is_open()) {
        std::cout << "Could not open file " << filepath << "\n";
        return;
    }
    for (size_t i = 0; i < event_nb.size(); i++) {
        file << event_nb[i] << " " << lift_nb[i] << " " << rev_nb[i] << " " << lambda[i] << "\n";
    }
    file.close();
    std::cout << "Number of events written in " << filepath << "\n";
}

// Saves the Mersenne Twister states into dirpath/filename/filename_seed/filename_seed.txt
void io::save_seed(std::mt19937_64& rng, const std::string& filename, const std::string& dirpath) {
    // Create a data folder if doesn't exists
    fs::path base_dir(dirpath);
    fs::path run_dir =
        base_dir / filename /
        (filename + "_seed");  // Utilise l'opérateur / pour gérer les slashs proprement

    try {
        // create_directories crée dirpath PUIS "data/run_name" si nécessaire
        if (!fs::exists(run_dir)) {
            fs::create_directories(run_dir);
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Error creating directory structure " << run_dir << " : " << e.what()
                  << std::endl;
        return;
    }
    fs::path filepath = run_dir / (filename + "_seed.txt");
    std::ofstream file(filepath);
    if (!file.is_open()) {
        std::cout << "Could not open file " << filepath << "\n";
        return;
    }
    file << rng;
    file.close();
    std::cout << "Seed saved in " << filepath << "\n";
};

// Saves the Mersenne Twister states into
// dirpath/filename/filename_seed/filename_seed_r[rank]_t[thread].txt
void io::save_seed(std::vector<std::mt19937_64>& rng, const std::string& filename,
                   const std::string& dirpath) {
    // Create a data folder if doesn't exists
    fs::path base_dir(dirpath);
    fs::path run_dir =
        base_dir / filename /
        (filename + "_seed");  // Utilise l'opérateur / pour gérer les slashs proprement
    try {
        // create_directories crée dirpath PUIS "data/run_name" si nécessaire_r
        if (!fs::exists(run_dir)) {
            fs::create_directories(run_dir);
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Error creating directory structure " << run_dir << " : " << e.what()
                  << std::endl;
        return;
    }
    for (size_t t = 0; t < rng.size(); ++t) {
        // Nom de fichier incluant le thread (t)
        // Exemple : ma_run_seed_t4.txt
        std::string seed_name = filename + "_seed_t" + std::to_string(t) + ".txt";
        fs::path filepath = run_dir / seed_name;

        std::ofstream file(filepath);
        if (!file.is_open()) {
            std::cerr << "Thread " << t << ": Could not open file " << filepath << "\n";
            continue;  // On essaie quand même les autres threads
        }

        // On sérialise l'état interne du générateur
        file << rng[t];
        file.close();
    }

    std::cout << "All threads seeds (" << rng.size() << ") saved in " << run_dir << "\n";
};

// Saves the local chain state of each core in
// dirpath/filename/filename_state/filename_state.txt
void io::save_state(const LocalChainState& state, const std::string& filename,
                    const std::string& dirpath) {
    // 1. Construction du chemin du dossier : dirpath/filename/filename_state/
    fs::path base_dir = fs::path(dirpath) / filename;
    fs::path state_dir = base_dir / (filename + "_state");

    // 2. Création du dossier (le rang 0 s'en occupe)
    if (!fs::exists(state_dir)) {
        fs::create_directories(state_dir);
    }

    // 3. Fichier spécifique au rang : filename_state[rank].txt
    std::string file_rank = filename + "_state.txt";
    fs::path full_path = state_dir / file_rank;

    std::ofstream ofs(full_path);
    if (!ofs.is_open()) {
        std::cerr << "Error: Could not open state file for writing: " << full_path << std::endl;
        return;
    }

    // 4. Écriture des données avec précision maximale
    ofs << std::setprecision(17);

    // Scalaires
    ofs << state.site << " " << state.mu << " " << state.epsilon << " "
        << state.theta_parcouru_refresh << " " << state.theta_refresh << " " << state.set_counter
        << " " << state.lift_counter << " " << state.rev_counter << " " << state.initialized
        << "\n";

    // Matrice SU3 R (3x3 nombres complexes)
    // On sauvegarde la partie réelle et imaginaire de chaque composante
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            auto val = state.R(i, j);
            ofs << val.real() << " " << val.imag() << " ";
        }
        ofs << "\n";
    }

    ofs.close();
    std::cout << "ECMC Chain saved in : " << full_path << std::endl;
}

void io::save_sweep_nb(int sweep_nb, const std::string& filename, const std::string& dirpath) {
    fs::path base_dir(dirpath);
    fs::path dir = base_dir / filename;

    try {
        if (!fs::exists(dir)) {
            fs::create_directories(dir);
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Couldn't create folder data : " << e.what() << std::endl;
        return;
    }
    fs::path filepath = dir / (filename + "_checkpoint.txt");

    std::ofstream file(filepath, std::ios::out);
    if (!file.is_open()) {
        std::cout << "Could not open file " << filepath << "\n";
        return;
    }

    file << sweep_nb;
    file.close();
}

void io::load_sweep_nb(int& sweep_nb, const std::string& filename, const std::string& dirpath) {
    fs::path full_path = fs::path(dirpath) / filename / (filename + "_checkpoint.txt");
    std::ifstream file(full_path);

    if (!file.is_open()) {
        throw std::runtime_error("Erreur : Impossible d'ouvrir le fichier " + full_path.string());
    }

    // Lecture de l'entier
    file >> sweep_nb;
    file.close();
}

// Loads the local chain state of each core
void io::load_state(LocalChainState& state, const std::string& filename,
                    const std::string& dirpath) {
    // 1. Construction du chemin du fichier
    fs::path state_path =
        fs::path(dirpath) / filename / (filename + "_state") / (filename + "_state.txt");

    // 2. Tentative d'ouverture
    std::ifstream ifs(state_path);
    if (!ifs.is_open()) {
        // Si le fichier n'existe pas, on considère que c'est un nouveau run
        state.initialized = false;
        return;
    }

    // 3. Lecture des scalaires
    ifs >> state.site >> state.mu >> state.epsilon >> state.theta_parcouru_refresh >>
        state.theta_refresh >> state.set_counter >> state.lift_counter >> state.rev_counter >>
        state.initialized;

    // 4. Lecture de la matrice SU3 R
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            double re, im;
            ifs >> re >> im;
            // On reconstruit le complexe et on l'assigne à la matrice
            state.R(i, j) = std::complex<double>(re, im);
        }
    }

    ifs.close();

    std::cout << "Successfully loaded ECMC chain state from: " << state_path << std::endl;
}

// Saves the run params into data/filename/filename_params.txt
void io::save_params(const RunParamsHbCB& rp, const std::string& filename,
                     const std::string& dirpath) {
    // Create a data/run_name folder if doesn't exists
    fs::path base_dir(dirpath);
    fs::path dir = base_dir / filename;

    try {
        if (!fs::exists(dir)) {
            fs::create_directories(dir);
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Couldn't create folder data : " << e.what() << std::endl;
        return;
    }
    fs::path filepath = dir / (filename + "_params.txt");

    std::ofstream file(filepath, std::ios::out | std::ios::app);
    if (!file.is_open()) {
        std::cout << "Could not open file " << filepath << "\n";
        return;
    }
    file << std::boolalpha;

    file << "\n#########################################################\n\n";

    file << "# Lattice params\n";
    file << "L=" << rp.L << "\n";
    file << "cold_start=" << rp.cold_start << "\n\n";

    file << "# Run params\n";
    file << "seed=" << rp.seed << "\n";
    file << "N_samples=" << rp.N_samples << "\n";

    file << "# Heatbath params\n";
    file << "beta = " << rp.hp.beta << "\n";
    file << "N_sweeps = " << rp.hp.N_sweeps << "\n";
    file << "N_hits = " << rp.hp.N_hits << "\n\n";

    file << "# Plaquette params\n";
    file << "N_plaquette= " << rp.N_plaquette << "\n";
    file << "T_min = " << rp.T_min << "\n";
    file << "T_max = " << rp.T_max << "\n\n";

    file << "#Save params\n";
    file << "save_each = " << rp.save_each << "\n\n";

    file.close();
    std::cout << "Parameters saved in " << filepath << "\n";
};

// Utilitary function to trim the spaces
std::string io::trim(const std::string& s) {
    size_t first = s.find_first_not_of(" \t");
    if (first == std::string::npos) return "";
    size_t last = s.find_last_not_of(" \t");
    return s.substr(first, (last - first + 1));
}

// Loads the params contained in filename into rp
void io::load_params(const std::string& filename, RunParamsECMC& rp) {
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Impossible d'ouvrir " + filename);

    std::map<std::string, std::string> config;
    std::string line;

    while (std::getline(file, line)) {
        // Ignore comments and empty lines
        if (line.empty() || line[0] == '#') continue;

        std::istringstream is_line(line);
        std::string key, value;
        if (std::getline(is_line, key, '=') && std::getline(is_line, value)) {
            config[trim(key)] = trim(value);
        }
    }

    // Lattice params
    if (config.count("T")) rp.T = std::stoi(config["T"]);
    if (config.count("L")) rp.L = std::stoi(config["L"]);
    if (config.count("cold_start")) rp.cold_start = (config["cold_start"] == "true");
    if (config.count("seed")) rp.seed = std::stoi(config["seed"]);
    if (config.count("N_samples")) rp.N_samples = std::stoi(config["N_samples"]);

    // ECMC params
    if (config.count("beta")) rp.ecmc_params.beta = std::stod(config["beta"]);
    // 1 sample each run bc frozen links increase artificially the correlation
    rp.ecmc_params.N_samples = 1;
    if (config.count("param_theta_sample"))
        rp.ecmc_params.param_theta_sample = std::stod(config["param_theta_sample"]);
    if (config.count("param_theta_refresh"))
        rp.ecmc_params.param_theta_refresh = std::stod(config["param_theta_refresh"]);
    if (config.count("poisson")) rp.ecmc_params.poisson = (config["poisson"] == "true");
    if (config.count("epsilon_set")) rp.ecmc_params.epsilon_set = std::stod(config["epsilon_set"]);
    // Plaquette
    if (config.count("N_plaquette")) rp.N_plaquette = std::stoi(config["N_plaquette"]);
    if (config.count("T_min")) rp.T_min = std::stoi(config["T_min"]);
    if (config.count("T_max")) rp.T_max = std::stoi(config["T_max"]);
    // Run name
    if (config.count("run_name")) rp.run_name = config["run_name"];
    if (config.count("run_dir")) rp.run_dir = config["run_dir"];
    if (config.count("save_each")) rp.save_each = std::stoi(config["save_each"]);
}

// Loads the params contained in filename into rp
void io::load_params(const std::string& filename, RunParamsHbCB& rp) {
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Can't open file " + filename);

    std::map<std::string, std::string> config;
    std::string line;

    while (std::getline(file, line)) {
        // Ignore comments and empty lines
        if (line.empty() || line[0] == '#') continue;

        std::istringstream is_line(line);
        std::string key, value;
        if (std::getline(is_line, key, '=') && std::getline(is_line, value)) {
            config[trim(key)] = trim(value);
        }
    }

    // Lattice params
    if (config.count("L")) rp.L = std::stoi(config["L"]);
    if (config.count("cold_start")) rp.cold_start = (config["cold_start"] == "true");
    if (config.count("seed")) rp.seed = std::stoi(config["seed"]);
    if (config.count("N_samples")) rp.N_samples = std::stoi(config["N_samples"]);

    // Hb params
    if (config.count("beta")) rp.hp.beta = std::stod(config["beta"]);
    // 1 sample bc frozen links increase artificially the correlation
    rp.hp.N_samples = 1;
    if (config.count("N_sweeps")) rp.hp.N_sweeps = std::stoi(config["N_sweeps"]);
    if (config.count("N_hits")) rp.hp.N_hits = std::stoi(config["N_hits"]);
    if (config.count("N_plaquette")) rp.N_plaquette = std::stoi(config["N_plaquette"]);
    // Plaquette
    if (config.count("N_plaquette")) rp.N_plaquette = std::stoi(config["N_plaquette"]);
    if (config.count("T_min")) rp.T_min = std::stoi(config["T_min"]);
    if (config.count("T_max")) rp.T_max = std::stoi(config["T_max"]);
    // Run name
    if (config.count("run_name")) rp.run_name = config["run_name"];
    if (config.count("run_dir")) rp.run_dir = config["run_dir"];
    if (config.count("save_each")) rp.save_each = std::stoi(config["save_each"]);
}

// Print parameters of the run
void print_parameters(const RunParamsHbCB& rp) {
    std::cout << "==========================================" << std::endl;
    std::cout << "Heatbath - Checkboard" << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << "Run name : " + rp.run_name << "\n";
    std::cout << "---Lattice params---\n";
    std::cout << "Total lattice size : " << rp.L << "^4\n";

    std::cout << "---Run params---\n";
    std::cout << "Initial seed : " << rp.seed << "\n";
    std::cout << "Number of samples : " << rp.N_samples << "\n";
    std::cout << "Save each : " << rp.save_each << " samples \n\n";

    std::cout << "---Heatbath params---\n";
    std::cout << "Beta : " << rp.hp.beta << "\n";
    std::cout << "Number of sweeps : " << rp.hp.N_sweeps << "\n";
    std::cout << "Number of hits : " << rp.hp.N_hits << "\n\n";

    std::cout << "---Measures params---\n";
    std::cout << "Number of <P> samples : " << rp.N_samples / rp.N_plaquette << "\n";
    std::cout << "Measure <P> each : " << rp.N_plaquette << " shifts\n";
    std::cout << "T_min = " << rp.T_min << ", T_max = " << rp.T_max << "\n";
    std::cout << "==========================================" << std::endl;
}

// Print time
void print_time(long elapsed) {
    std::cout << "==========================================" << std::endl;
    std::cout << "Elapsed time : " << elapsed << "s\n";
    std::cout << "==========================================" << std::endl;
}

// to_string for double with a fixed precision
std::string io::format_double(double val, int precision) {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(precision) << val;
    return ss.str();
}

// Reads the parameters of input file into RunParams struct
// Returns true if the necessary files are found
bool io::read_params(RunParamsHbCB& params, const std::string& input) {
    try {
        io::load_params(input, params);
    } catch (const std::exception& e) {
        std::cerr << "Error reading input : " << e.what() << std::endl;
    }

    // 4. Vérification de l'existence des fichiers pour la reprise (Resume)
    fs::path base_path = fs::path(params.run_dir) / params.run_name;
    fs::path config_file = base_path / params.run_name;  // Le fichier ILDG
    fs::path seed_dir = base_path / (params.run_name + "_seed");

    // 1. On vérifie d'abord si le fichier de configuration global existe
    bool local_existing = fs::exists(config_file);
    return local_existing;
}

void print_parameters(const RunParamsECMC& rp) {
    std::cout << "==========================================" << std::endl;
    std::cout << "ECMC - OBC " << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << "Run name : " + rp.run_name << "\n";
    std::cout << "---Lattice params---\n";
    std::cout << "Total lattice size : " << rp.T << "x" << rp.L << "^3\n";

    std::cout << "---Run params---\n";
    std::cout << "Initial seed : " << rp.seed << "\n";
    std::cout << "Number of samples : " << rp.N_samples << "\n";
    std::cout << "Save each : " << rp.save_each << " samples\n\n";

    std::cout << "---ECMC params---\n";
    std::cout << "Beta : " << rp.ecmc_params.beta << "\n";
    std::cout << "Theta sample : " << rp.ecmc_params.param_theta_sample << "\n";
    std::cout << "Theta refresh : " << rp.ecmc_params.param_theta_refresh<< "\n";
    std::cout << "Epsilon set : " << rp.ecmc_params.epsilon_set << "\n";
    std::cout << "Poisson law : " << (rp.ecmc_params.poisson ? "Yes" : "No") << "\n\n";

    std::cout << "---Measures params---\n";
    std::cout << "Number of <P> samples : " << rp.N_samples / rp.N_plaquette << "\n";
    std::cout << "Measure <P> each : " << rp.N_plaquette << " shifts\n";
    std::cout << "T_min = " << rp.T_min << ", T_max = " << rp.T_max << "\n";
    std::cout << "==========================================" << std::endl;
}

void io::save_params(const RunParamsECMC& rp, const std::string& filename,
                     const std::string& dirpath) {
    // Create a data folder if doesn't exists
    fs::path base_dir(dirpath);
    fs::path dir = base_dir / filename;

    try {
        if (!fs::exists(dir)) {
            fs::create_directories(dir);
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Couldn't create folder data : " << e.what() << std::endl;
        return;
    }
    fs::path filepath = dir / (filename + "_params.txt");

    std::ofstream file(filepath, std::ios::out | std::ios::app);
    if (!file.is_open()) {
        std::cout << "Could not open file " << filepath << "\n";
        return;
    }
    file << std::boolalpha;

    file << "\n#########################################################\n\n";

    file << "# Lattice params\n";
    file << "T = " << rp.T << "\n";
    file << "L = " << rp.L << "\n";
    file << "cold_start = " << rp.cold_start << "\n\n";

    file << "# Run params\n";
    file << "seed = " << rp.seed << "\n";
    file << "N_samples= " << rp.N_samples << "\n";

    file << "# ECMC params\n";
    file << "beta = " << rp.ecmc_params.beta << "\n";
    file << "theta_sample = " << rp.ecmc_params.param_theta_sample << "\n";
    file << "theta_refresh= " << rp.ecmc_params.param_theta_refresh << "\n";
    file << "epsilon_set = " << rp.ecmc_params.epsilon_set << "\n";
    file << "poisson = " << rp.ecmc_params.poisson << "\n\n";

    file << "# Plaquette params\n";
    file << "N_plaquette = " << rp.N_plaquette << "\n";
    file << "T_min = " << rp.T_min << "\n";
    file << "T_max = " << rp.T_max << "\n\n";

    file << "#Save params\n";
    file << "save_each= " << rp.save_each << "\n\n";

    file.close();
    std::cout << "Parameters saved in " << filepath << "\n";
};

// Reads the parameters of input file into RunParams struct
bool io::read_params(RunParamsECMC& params, const std::string& input) {
    try {
        io::load_params(input, params);
    } catch (const std::exception& e) {
        std::cerr << "Error reading input : " << e.what() << std::endl;
    }
    // 4. Vérification de l'existence des fichiers pour la reprise (Resume)
    fs::path base_path = fs::path(params.run_dir) / params.run_name;
    fs::path config_file = base_path / params.run_name;  // Le fichier ILDG
    fs::path seed_dir = base_path / (params.run_name + "_seed");

    bool local_existing = fs::exists(config_file);
    return local_existing;
}

// Adds the log of saved shift to params
void io::add_shift(int shift, const std::string& filename, const std::string& dirpath) {
    fs::path base_dir(dirpath);
    fs::path dir = base_dir / filename;

    try {
        if (!fs::exists(dir)) {
            fs::create_directories(dir);
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Couldn't create folder data : " << e.what() << std::endl;
        return;
    }
    fs::path filepath = dir / (filename + "_params.txt");

    std::ofstream file(filepath, std::ios::out | std::ios::app);
    if (!file.is_open()) {
        std::cout << "Could not open file " << filepath << "\n";
        return;
    }

    file << "Saved shift " << shift << "\n";
    file.close();
};

// Add finished to params
void io::add_finished(const std::string& filename, const std::string& dirpath) {
    fs::path base_dir(dirpath);
    fs::path dir = base_dir / filename;

    try {
        if (!fs::exists(dir)) {
            fs::create_directories(dir);
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Couldn't create folder data : " << e.what() << std::endl;
        return;
    }
    fs::path filepath = dir / (filename + "_params.txt");

    std::ofstream file(filepath, std::ios::out | std::ios::app);
    if (!file.is_open()) {
        std::cout << "Could not open file " << filepath << "\n";
        return;
    }

    file << "Saved final state !\n";
    file.close();
};
