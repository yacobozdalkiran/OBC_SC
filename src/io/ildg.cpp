#include "ildg.h"

#include <filesystem>
#include <iostream>

extern "C" {
#include <lime.h>
}

namespace fs = std::filesystem;

std::string generate_ildg_xml(int lx, int ly, int lz, int lt, int precision) {
    std::ostringstream oss;
    oss << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        << "<ildgFormat xmlns=\"http://www.lqcd.org/ildg\" "
        << "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
        << "xsi:schemaLocation=\"http://www.lqcd.org/ildg http://www.lqcd.org/ildg/ildg-format.xsd\">\n"
        << "  <version>1.1</version>\n"
        << "  <field>su3gauge</field>\n"
        << "  <precision>" << precision << "</precision>\n"
        << "  <lx>" << lx << "</lx>\n"
        << "  <ly>" << ly << "</ly>\n"
        << "  <lz>" << lz << "</lz>\n"
        << "  <lt>" << lt << "</lt>\n"
        << "</ildgFormat>";
    return oss.str();
}

void save_ildg(const GaugeField &field, const Geometry &geo, const std::string &run_name, const std::string& run_dir) {
    //Filepath for conf
    fs::path dir(run_dir);
    try {
        if (!fs::exists(dir)) {
            fs::create_directories(dir);
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Couldn't create folder : " << e.what() << std::endl;
        return;
    }
    fs::path final_path = dir/run_name/run_name;
    fs::path tmp_path = final_path;
    tmp_path += ".tmp";

    // 1. Prepare XML (double precision)
    std::string xml_header = generate_ildg_xml(geo.L, geo.L, geo.L, geo.T); // Lx, Ly, Lz, Lt

    // 2. Open C-LIME writer on temporary file
    FILE *fp = fopen(tmp_path.c_str(), "wb");
    if (!fp) throw std::runtime_error("Impossible de créer le fichier temporaire : " + tmp_path.string());
    LimeWriter *writer = limeCreateWriter(fp);

    // 3. Write the XML Record
    int MB_flag = 1, ME_flag = 0;
    n_uint64_t xml_len = xml_header.size();
    LimeRecordHeader *h_xml = limeCreateHeader(MB_flag, ME_flag, (char*)"ildg-format", xml_len);
    limeWriteRecordHeader(h_xml, writer);
    limeWriteRecordData((void*)xml_header.c_str(), &xml_len, writer);
    limeDestroyHeader(h_xml);

    // 4. Prepare and write the binary data
    // Physical volume V (no halos)
    n_uint64_t nbytes = (n_uint64_t) geo.V * 4 * 9 * sizeof(std::complex<double>);
    ME_flag = 1;
    LimeRecordHeader *h_bin = limeCreateHeader(0, ME_flag, (char*)"ildg-binary-data", nbytes);
    limeWriteRecordHeader(h_bin, writer);

    for (size_t site = 0; site < geo.V; site++) {
        for (int mu = 0; mu < 4; mu++) {
            // Copy of the matrix
            SU3 matrix = field.view_link_const(site, mu);

            // Endianess inversion
            double* raw_data = reinterpret_cast<double*>(matrix.data());
            for (int i = 0; i < 18; i++) {
                uint64_t *u = reinterpret_cast<uint64_t*>(&raw_data[i]);
                *u = __builtin_bswap64(*u);
            }

            // Matrix writing (18 doubles)
            n_uint64_t mat_size = 18 * sizeof(double);
            limeWriteRecordData(raw_data, &mat_size, writer);
        }
    }

    limeDestroyHeader(h_bin);
    limeDestroyWriter(writer);
    fclose(fp);

    // 5. Atomic rename
    try {
        fs::rename(tmp_path, final_path);
    } catch (const fs::filesystem_error& e) {
        throw std::runtime_error("Erreur lors du renommage du fichier ILDG : " + std::string(e.what()));
    }

    std::cout << "Configuration written in " << final_path << "\n";
}

void read_ildg(GaugeField &field, const Geometry &geo, const std::string &run_name, const std::string& run_dir) {
    fs::path dir(run_dir);
    fs::path filepath = dir/run_name/run_name;
    FILE *fp = fopen(filepath.c_str(), "rb");
    if (!fp) throw std::runtime_error("Impossible d'ouvrir le fichier : " + filepath.string());

    std::cout << "Reading configuration from " << filepath << "\n";
    LimeReader *reader = limeCreateReader(fp);
    bool data_found = false;

    // Parcourir les records du fichier LIME
    while (limeReaderNextRecord(reader) == LIME_SUCCESS) {
        char* type = limeReaderType(reader);

        // On cherche uniquement le record binaire
        if (std::string(type) == "ildg-binary-data") {
            data_found = true;

            n_uint64_t nbytes = limeReaderBytes(reader);
            n_uint64_t expected_bytes = (n_uint64_t) geo.V * 4 * 9 * sizeof(std::complex<double>);

            if (nbytes != expected_bytes) {
                throw std::runtime_error("Taille du fichier ILDG incohérente avec le volume du GaugeField !");
            }

            // Lecture site par site, mu par mu
            for (size_t site = 0; site < geo.V; ++site) {
                for (int mu = 0; mu < 4; ++mu) {

                    // On lit 18 doubles (9 complexes)
                    double buffer[18];
                    n_uint64_t to_read = 18 * sizeof(double);
                    limeReaderReadData(buffer, &to_read, reader);

                    // Conversion Big-Endian -> Little-Endian
                    for (int i = 0; i < 18; ++i) {
                        uint64_t *u = reinterpret_cast<uint64_t*>(&buffer[i]);
                        *u = __builtin_bswap64(*u);
                    }

                    // On "map" le buffer vers une matrice Eigen temporaire (Row-Major)
                    // Puis on l'affecte à notre view_link

                    // Méthode plus sûre pour copier les données vers votre vecteur 'links' via le Map
                    auto link = field.view_link(site, mu);
                    for(int i=0; i<3; ++i) {
                        for(int j=0; j<3; ++j) {
                            // ILDG stocke : Re(0,0), Im(0,0), Re(0,1), Im(0,1)...
                            double re = buffer[2 * (i * 3 + j)];
                            double im = buffer[2 * (i * 3 + j) + 1];
                            link(i, j) = std::complex<double>(re, im);
                        }
                    }
                }
            }
            break; // Données lues, on peut sortir
        }
    }

    if (!data_found) {
        throw std::runtime_error("Record 'ildg-binary-data' has not been found in the file.");
    }

    limeDestroyReader(reader);
    fclose(fp);
}
