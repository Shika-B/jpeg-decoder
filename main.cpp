#include <iostream>
#include <fstream>
#include <iterator>
#include <utility>
#include <vector>
#include <array>
#include <numeric>
#include <span>
#include <cassert>


#include "utils.h"

enum class DensityUnit {NoUnit, PixelPerInch, PixelPerCm};

struct RGB {
    unsigned char r;
    unsigned char g;
    unsigned char b;
};

struct JFIFVersion {
    unsigned char major;
    unsigned char minor;
};


struct JFIFData {
    JFIFVersion version;
    DensityUnit density_unit;
    unsigned short x_density;
    unsigned short y_density;
    unsigned char x_thumbnail;
    unsigned char y_thumbnail;
    std::vector<RGB> thumbnail_data;
};


struct QuantizationTable  {
    std::array<unsigned short, 64> data;
};

struct HuffmanCode {
    unsigned char length;
    unsigned short code;
    // The value the code is mapped to.
    unsigned char value;
};


struct HuffmanTable {
    // Codes sorted by code length
    std::vector<HuffmanCode> codes;
    static HuffmanTable from_size_data(std::array<unsigned char, 16> size_data, std::vector<unsigned char> data_area);

private:
    static std::vector<unsigned short> make_code_table(std::vector<unsigned char> sizes);
};

HuffmanTable HuffmanTable::from_size_data(std::array<unsigned char, 16> size_data, std::vector<unsigned char> data_area) {
    std::vector<unsigned char> codes_length {};
    for (int i = 0; i < size_data.size(); i++) {
        for (int j = 0; j < size_data[i]; j++) {
            codes_length.push_back(i+1);
        }
    }
    std::vector<unsigned short> code_table = HuffmanTable::make_code_table(codes_length);
    std::vector<HuffmanCode> huffman_codes {};
    for (int i = 0; i < data_area.size(); i++) {

        huffman_codes.push_back(HuffmanCode {codes_length[i], code_table[i], data_area[i]});
    }
    return HuffmanTable { huffman_codes };
}


std::vector<unsigned short> HuffmanTable::make_code_table(std::vector<unsigned char> sizes) {
    std::vector<unsigned short> code_table {};
    unsigned short code = 0;
    unsigned char current_size = sizes[0];
    for (unsigned char size: sizes) {

        while (size > current_size) {
            code <<= 1;
            current_size +=1;
        }
        code_table.push_back(code);
        if (current_size > 16 || code == 0xffff) {
            break;
        };
        code += 1;
    }
    return code_table;
}



struct JPEGEncoded {
    JFIFData metadata;
    std::array<QuantizationTable, 4> q_tables;
    unsigned char q_tables_nbr;
    std::array<HuffmanTable, 32> huffman_ac_tables;
    std::array<HuffmanTable, 32> huffman_dc_tables;
};

class JPEGParser {
public:
    JPEGEncoded parse();
    explicit JPEGParser(std::vector<unsigned char> raw_data) noexcept:
    raw_data(std::move(raw_data)), index(0){};


private:
    std::vector<unsigned char> raw_data;
    unsigned long long index;
    JFIFData parse_jfif_data();
    QuantizationTable parse_quantization_table();
    HuffmanTable parse_huffman_table();
};


JFIFData JPEGParser::parse_jfif_data() {
    index += 5; // Ignores the 5 (constant) identifier bytes

    JFIFVersion version {raw_data[index], raw_data[index+1]};
    index += 2;

    DensityUnit density_unit = (DensityUnit) raw_data[index];
    index += 1;

    unsigned short x_density = u8_to_u16(raw_data[index], raw_data[index+1]);
    unsigned short y_density = u8_to_u16(raw_data[index+2], raw_data[index+3]);
    index += 4;

    unsigned char x_thumbnail = raw_data[index];
    unsigned char y_thumbnail = raw_data[index+1];
    index += 2;

    std::vector<RGB> thumbnail_data;
    for (int i = 0; i < x_thumbnail*y_thumbnail; i++) {
        unsigned char r = raw_data[i];
        unsigned char g = raw_data[i+1];
        unsigned char b = raw_data[i+2];
        thumbnail_data.push_back(RGB {r, g, b});
    }
    return JFIFData {version, density_unit, x_density, y_density, x_thumbnail, y_thumbnail, thumbnail_data};
}

QuantizationTable JPEGParser::parse_quantization_table() {
    // Precision 0 if the Quantization table contains 8-bit integers, 1 if it contains 16-bit integers.
    // In the case of Baseline DCT encoding (the only one supported here), precision is always 0.
    [[maybe_unused]]
    unsigned char precision = (raw_data[index] & 0xf0) >> 4;

    // Index (ranging from 0 to 3) of the table
    [[maybe_unused]]
    unsigned char identifier = raw_data[index] & 0x0f;
    index += 1;

    std::array<unsigned short, 64> table_data {};
    std::copy(raw_data.begin() + index, raw_data.begin() + index + 64, table_data.begin());
    index += 64;

    QuantizationTable table { table_data };

    return table;
}

HuffmanTable JPEGParser::parse_huffman_table() {

    std::array<unsigned char, 16> size_data {};
    std::copy(raw_data.begin() + index, raw_data.begin() + index + 16, size_data.begin());
    index += 16;

    int codes_count = std::accumulate(size_data.begin(), size_data.end(), (int) 0);

    std::vector<unsigned char> data_area(raw_data.begin() + index, raw_data.begin() + index + codes_count);
    index += codes_count;

    return HuffmanTable::from_size_data(size_data, data_area);
}


JPEGEncoded JPEGParser::parse() {
    JFIFData jfif_data;

    std::array<QuantizationTable, 4> q_tables {};
    unsigned char q_tables_nbr = 0;

    std::array<HuffmanTable, 32> h_ac_tables {};
    std::array<HuffmanTable, 32> h_dc_tables {};
    while (index < raw_data.size()) {
        unsigned char b = raw_data[index];
        if (b == 0xff) {
            unsigned char marker = raw_data[index + 1];

            index += 2;

            // Those markers have no length, so they must be treated separately
            if (marker == 0xd8 or marker == 0xd9) {
                continue;
            }

            unsigned short length = ((unsigned short)raw_data[index] << 8) + raw_data[index + 1];
            size_t segment_end = index + length;
            index += 2;

            std::cout << "Parsing marker: " << std::hex << (int) marker << ", index: " << std::dec << index << std::endl;

            switch (marker) {
                case 0xe0:
                    jfif_data = parse_jfif_data();
                    break;
                case 0xdb: {
                    while (index < segment_end) {
                        QuantizationTable table = parse_quantization_table();
                        q_tables[q_tables_nbr] = table;
                        q_tables_nbr += 1;
                    };
                    break;
                }
                case 0xc4: {
                    while (index < segment_end) {
                        unsigned char table_class = (raw_data[index] & 0xf0) >> 4;
                        unsigned char table_dest_id = raw_data[index] & 0x0f;
                        index += 1;
                        HuffmanTable table = parse_huffman_table();

                        if (table_class == 0) {
                            h_dc_tables[table_dest_id] = table;
                        } else {
                            h_ac_tables[table_dest_id] = table;
                        }
                    };
                    break;
                }
                case 0xc0: {
                    index  += length - 2;
                    break;
                }
                case 0xda: {
                    // We are assuming there is a SINGLE Scan.
                    index += raw_data.size() - index;
                    break;
                }
                default:
                    std::cout << "Ignored unknown marker " << std::hex << (int) marker << " of length " << std::dec << length << std::endl;
                    index  += length - 2;
                    break;
            }
            // std::cout << std::dec << index << " " << segment_end << std::endl;
        } else {
            std::cout << std::hex << (int) b << std::endl;
            throw;
        }
    }
    return JPEGEncoded {jfif_data, q_tables, q_tables_nbr, h_ac_tables, h_dc_tables};
}



int main() {
    const char* input_file = R"(C:\Users\abdel\CLionProjects\jpeg_parser\sample1.jfif)";

    std::ifstream input(input_file, std::ios::binary);
    std::vector<unsigned char> buffer(std::istreambuf_iterator<char>(input), {});

    JPEGParser parser {buffer};
    JPEGEncoded jpeg_encoded = parser.parse();
    JFIFData metadata = jpeg_encoded.metadata;
    std::cout << "JFIF Version: " << (int) metadata.version.major << "." << (int) metadata.version.minor << std::endl;
    std::cout << "Thumbnail size: " << (int) metadata.x_thumbnail << "x" << (int) metadata.y_thumbnail << std::endl;
    std::cout << "XY density: " << metadata.x_density << "x" << metadata.y_density << std::endl;
    std::cout << "Quantization tables number: " << (int) jpeg_encoded.q_tables_nbr << std::endl;
    std::cout << "Finished" << std::endl;
    return 0;
}
